#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Why an exposure, or a whole visit, is an outlier.

koloa's mixture says which exposures and visits are outliers; this module
looks for why, in what was recorded with each exposure: the header keywords
that LBL copies into its .rdb files (S/N, airmass, exposure time, telluric
absorption, the shape of the image, the Fabry-Perot of the calibration, the
age of the wavelength solution...), the activity indicators (DTEMP, d2v,
dW, FWHM, CRX...) and the error bar of the velocity itself. For every
outlier it asks which of them are significantly off:

- a single exposure is compared with the good exposures of its own visit
  (the same night, the same star, the same sky), in units of the scatter of
  a good exposure about the others of its visit;
- a whole visit, or a visit of one exposure, is compared with the good
  visits nearest in time, in units of the scatter of a good visit about its
  neighbours, so that seasons, slow drifts and activity cycles are not
  taken for anomalies.

A number is significantly off when its deviation z is beyond the normal
threshold for the number of keys tested (a Bonferroni correction: with 30
keys at alpha = 0.01, |z| > 3.6) and at most 1 % of the good exposures (or
visits) deviate as much, a guard against keys with heavy tails. A text key
(DPRTYPE, the position of a rhomb) is off when its value is rare among the
good exposures (under 5 %). A permutation test then asks the same of all
the outliers of an instrument together: which keys they share (a lower
S/N, larger error bars...), Holm-corrected over the keys.

    report = explain(fit)                  # a koloa fit with the mixture
    print(report.text())
    report.units[0]['keys']                # the keys of the first outlier
    kplot.outlier_keys(report)             # the figure

The keys: keys='auto' takes the list of each instrument (SPIROU_KEYS,
NIRPS_KEYS or DACE_KEYS, recognised from its name or its columns), keeps
those the series has, and adds its activity indicators; a list (or the name
of a list) mixes names of those lists, of RVData.meta and of
RVData.indicators, and 'svrad', the error bar of the velocity.

Created on 2026-09-29

@author: artigau
"""
import datetime
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy import stats

from koloa.data import RVData

# =============================================================================
# Define variables
# =============================================================================
#: an exposure is an outlier above this probability, and good below GOOD
THRESHOLD = 0.5
GOOD = 0.1
#: a text value rarer than this among the good exposures is off
RARE = 0.05
#: a number is off only if at most this fraction of the good units deviate
#: as much
TAIL = 0.01
#: the association test clips every z at this, so that one extreme outlier
#: cannot make a trend alone
CLIP = 3.0


@dataclass
class Key:
    """
    A quantity recorded with each exposure

    :param name: the column (of RVData.meta or RVData.indicators), 'svrad'
                 for the error bar of the velocity, or the name of a
                 derived key
    :param label: what it is, in words (with its unit)
    :param bad: where a problem would lie: 'low', 'high' or 'either'
    :param value: for a derived key, a function of the series that returns
                  one value per exposure (a KeyError when the series lacks
                  what it needs)
    """
    name: str
    label: str
    bad: str = 'either'
    value: Optional[Callable[[RVData], np.ndarray]] = None

    def values(self, data: RVData
               ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        The value of every exposure (numbers or text) and its error

        :return: tuple, (values, errors or None), or (None, None) when the
                 series does not have the key
        """
        if self.value is not None:
            try:
                with np.errstate(invalid='ignore', divide='ignore'):
                    return np.asarray(self.value(data), dtype=float), None
            except (KeyError, TypeError, ValueError):
                return None, None
        if self.name == 'svrad':
            return data.err, None
        if self.name in data.indicators:
            val, err = data.indicators[self.name]
            return np.asarray(val, dtype=float), np.asarray(err, dtype=float)
        if self.name in data.meta:
            return np.asarray(data.meta[self.name]), None
        return None, None


def _meta(name: str) -> np.ndarray:
    """a column of the meta, as numbers (for the derived keys)"""
    return lambda data: np.asarray(data.meta[name], dtype=float)


def _snr_rate(snr: str, exptime: str):
    """S/N squared per second: the flux collected, whatever the exposure
    time (clouds, seeing, guiding lower it)"""
    return lambda data: (_meta(snr)(data) ** 2 / _meta(exptime)(data))


def _ratio(num: str, den: str):
    """one column over another"""
    return lambda data: _meta(num)(data) / _meta(den)(data)


def _wave_age(data: RVData) -> np.ndarray:
    """days between the exposure and its wavelength solution"""
    return np.abs(_meta('MJDMID')(data) - _meta('WAVETIME')(data))


#: the keys of every instrument: the error bar of the velocity
COMMON_KEYS = [Key('svrad', 'error bar of the velocity [m/s]', 'high')]
#: what the image and the reduction (APERO, LBL) record, SPIRou and NIRPS
_APERO_KEYS = [
    Key('TLPEH2O', 'telluric water absorption (TLPEH2O)', 'high'),
    Key('TLPEOTR', 'other telluric absorption (TLPEOTR)', 'high'),
    Key('TLPDVH2O', 'velocity of the water model (TLPDVH2O)'),
    Key('TLPDVOTR', 'velocity of the other absorbers (TLPDVOTR)'),
    Key('SHAPE_DX', 'shift of the image along x (SHAPE_DX)'),
    Key('SHAPE_DY', 'shift of the image along y (SHAPE_DY)'),
    Key('SHAPE_A', 'shape of the image (SHAPE_A)'),
    Key('SHAPE_B', 'shape of the image (SHAPE_B)'),
    Key('SHAPE_C', 'shape of the image (SHAPE_C)'),
    Key('SHAPE_D', 'shape of the image (SHAPE_D)'),
    Key('CCF_EW', 'width of the CCF (CCF_EW)'),
    Key('ITE_RV', 'LBL iterations (ITE_RV)', 'high'),
    Key('RESET_RV', 'LBL reset of the velocity (RESET_RV)', 'high'),
    Key('wave_age', 'age of the wavelength solution [d]', 'high', _wave_age),
    Key('DPRTYPE', 'type of observation (DPRTYPE)')]
#: SPIRou (APERO and LBL keywords)
SPIROU_KEYS = COMMON_KEYS + [
    Key('EXTSN035', 'S/N (order 35, EXTSN035)', 'low'),
    Key('snr_rate', 'S/N squared per second (EXTSN035^2 / EXPTIME)', 'low',
        _snr_rate('EXTSN035', 'EXPTIME')),
    Key('snr_goal', 'S/N over its goal (EXTSN035 / SNRGOAL)', 'low',
        _ratio('EXTSN035', 'SNRGOAL')),
    Key('EXPTIME', 'exposure time [s]'),
    Key('AIRMASS', 'airmass', 'high')] + _APERO_KEYS + [
    Key('SBCFPI_T', 'Fabry-Perot temperature (SBCFPI_T) [C]'),
    Key('SBCFPE_T', 'Fabry-Perot temperature (SBCFPE_T) [C]'),
    Key('SBCFPB_P', 'Fabry-Perot pressure (SBCFPB_P)'),
    Key('SBCDEN_P', 'density of the calibration (SBCDEN_P)'),
    Key('SBRHB1_P', 'position of rhomb 1 (SBRHB1_P)'),
    Key('SBRHB2_P', 'position of rhomb 2 (SBRHB2_P)')]
#: NIRPS (APERO and LBL keywords, and ESO's)
NIRPS_KEYS = COMMON_KEYS + [
    Key('EXTSN060', 'S/N (order 60, EXTSN060)', 'low'),
    Key('snr_rate', 'S/N squared per second (EXTSN060^2 / EXPTIME)', 'low',
        _snr_rate('EXTSN060', 'EXPTIME')),
    Key('EXPTIME', 'exposure time [s]'),
    Key('HIERARCH ESO TEL AIRM START', 'airmass', 'high'),
    Key('HIERARCH ESO TEL AMBI FWHM START', 'seeing [arcsec]', 'high')
    ] + _APERO_KEYS + [
    Key('HIERARCH ESO INS TEMP13 VAL', 'instrument temperature (INS TEMP13)'),
    Key('HIERARCH ESO INS TEMP14 VAL', 'instrument temperature (INS TEMP14)'),
    Key('HIERARCH ESO INS PRES108 VAL', 'instrument pressure (INS PRES108)')]
#: HARPS, HARPS-N, ESPRESSO, CORALIE... from DACE (koloa.dace.rvdata)
DACE_KEYS = COMMON_KEYS + [
    Key('spectro_sn50', 'S/N (order 50)', 'low'),
    Key('snr_rate', 'S/N squared per second (order 50)', 'low',
        _snr_rate('spectro_sn50', 'texp')),
    Key('texp', 'exposure time [s]'),
    Key('cal_drift_noise', 'noise of the drift of the reference', 'high'),
    Key('cal_drift_rv', 'drift of the simultaneous reference'),
    Key('cal_therror', 'error of the wavelength calibration', 'high'),
    Key('ccf_noise', 'photon noise of the CCF', 'high'),
    Key('ccf_asym', 'asymmetry of the CCF'),
    Key('spectro_ca', 'Ca II index'),
    Key('spectro_na', 'Na I index'),
    Key('dpr_type', 'type of observation (dpr_type)')]
LISTS = {'spirou': SPIROU_KEYS, 'nirps': NIRPS_KEYS, 'dace': DACE_KEYS}
NAMES = {'spirou': 'SPIRou', 'nirps': 'NIRPS', 'dace': 'DACE'}
#: the words for the activity indicators that LBL and DACE give
INDICATOR_LABELS = {
    'd2v': 'second derivative of the lines (d2v)',
    'd3v': 'third derivative of the lines (d3v)',
    'dW': 'change of the line width (dW)',
    'fwhm': 'FWHM of the lines',
    'contrast': 'contrast of the lines',
    'CRX': 'chromatic index (CRX)',
    'vrad_chromatic_slope': 'chromatic slope of the velocity',
    'bis': 'bisector span', 'rhk': "log R'HK", 'smw': 'S index',
    'halpha': 'H alpha index'}


# =============================================================================
# Define functions
# =============================================================================
def instrument_list(data: RVData, inst: str) -> Tuple[str, List[Key]]:
    """
    The list of keys of an instrument: from its name (SPIRou, NIRPS, the
    DACE eras HARPS03, ESPRESSO19...), or from the columns of its exposures

    :return: tuple, the name of the list ('spirou', 'nirps', 'dace' or
             'generic') and its keys
    """
    lower = inst.lower()
    for name in ('spirou', 'nirps'):
        if name in lower:
            return name, LISTS[name]
    if any(lower.startswith(era) for era in ('harps', 'espresso', 'coralie',
                                              'harpn', 'sophie')):
        return 'dace', LISTS['dace']
    cols = set(data.meta)
    if 'EXTSN035' in cols:
        return 'spirou', LISTS['spirou']
    if 'EXTSN060' in cols or any(col.startswith('HIERARCH ESO')
                                 for col in cols):
        return 'nirps', LISTS['nirps']
    if 'spectro_sn50' in cols:
        return 'dace', LISTS['dace']
    return 'generic', list(COMMON_KEYS)


def indicator_keys(data: RVData) -> List[Key]:
    """
    The activity indicators of a series as keys: every indicator but the
    velocities in parts of the spectrum (vrad_*, which hold the outlier as
    much as the velocity does), except the chromatic slope
    """
    out = []
    for name in data.indicators:
        if name.startswith(('vrad', 'svrad')) and name != \
                'vrad_chromatic_slope':
            continue
        label = INDICATOR_LABELS.get(name)
        if label is None and name.upper().startswith('DTEMP'):
            label = f'temperature change ({name}) [K]'
        out.append(Key(name, label or name))
    return out


def resolve_keys(data: RVData, inst: str,
                 keys: Union[str, Sequence[Union[str, Key]]] = 'auto'
                 ) -> Tuple[str, List[Key], List[str]]:
    """
    The keys to test for an instrument, those it has

    :param keys: 'auto', the name of a list ('spirou', 'nirps', 'dace'), or
                 a list of names (columns, indicators, keys of the lists,
                 'svrad') and Key objects

    :return: tuple, the list's name, the keys the instrument has, and the
             names of those it lacks
    """
    sel = data.select(data.inst == inst)
    if isinstance(keys, str):
        if keys == 'auto':
            lname, wanted = instrument_list(sel, inst)
            wanted = wanted + indicator_keys(sel)
        else:
            lname, wanted = keys, LISTS[keys.lower()] + indicator_keys(sel)
    else:
        lname, known = 'custom', {}
        _, own = instrument_list(sel, inst)
        for group in list(LISTS.values()) + [own]:
            known.update({key.name: key for key in group})
        wanted = []
        for key in keys:
            if isinstance(key, Key):
                wanted.append(key)
            elif key in known:
                wanted.append(known[key])
            else:
                wanted.append(Key(key, INDICATOR_LABELS.get(key, key)))
    have, lacking, seen = [], [], set()
    for key in wanted:
        if key.name in seen:
            continue
        seen.add(key.name)
        val, _ = key.values(sel)
        if val is None or not _usable(val):
            lacking.append(key.name)
        else:
            have.append(key)
    return lname, have, lacking


def _is_text(values: np.ndarray) -> bool:
    """a column of text"""
    return values.dtype.kind in 'USO'


def _usable(values: np.ndarray) -> bool:
    """a column with at least three values"""
    if _is_text(values):
        return np.sum(np.char.strip(values.astype(str)) != '') >= 3
    return np.sum(np.isfinite(values)) >= 3


def _scale(dev: np.ndarray) -> float:
    """the robust standard deviation of deviations (1.4826 MAD), or their
    standard deviation when most of them are equal"""
    dev = dev[np.isfinite(dev)]
    if len(dev) < 2:
        return np.nan
    mad = 1.4826 * np.median(np.abs(dev - np.median(dev)))
    return float(mad) if mad > 0 else float(np.std(dev))


def _zscore(delta: float, scale: float) -> float:
    """a deviation in units of the scale; infinite when the good units
    never vary and this one differs"""
    if not np.isfinite(delta):
        return np.nan
    if scale > 0:
        return float(delta / scale)
    return 0.0 if delta == 0 else float(np.sign(delta) * np.inf)


def _date(rjd: float, with_time: bool = False) -> str:
    """the calendar date (UTC, near enough) of a time in BJD - 2400000"""
    moment = datetime.datetime(1858, 11, 17) + datetime.timedelta(
        days=float(rjd) - 0.5)
    return moment.strftime('%Y-%m-%d %H:%M' if with_time else '%Y-%m-%d')


def _holm(pvals: Sequence[float]) -> np.ndarray:
    """Holm's step-down adjustment of p-values (for testing many keys)"""
    pvals = np.asarray(pvals, dtype=float)
    order = np.argsort(pvals)
    adj = np.empty(len(pvals))
    running = 0.0
    for rank, it in enumerate(order):
        running = max(running, (len(pvals) - rank) * pvals[it])
        adj[it] = min(running, 1.0)
    return adj


# =============================================================================
# The units: every outlier, with what it is compared with
# =============================================================================
def _units(data: RVData, prob: np.ndarray, inst: str, threshold: float,
           good: float) -> Tuple[List[Dict[str, Any]], np.ndarray]:
    """
    The outliers of an instrument: whole visits (every exposure an outlier)
    and single exposures; each with the exposures it is compared with

    :return: tuple, the units and the good visits (all of their exposures
             good)
    """
    units, good_visits = [], []
    idx = np.where(data.inst == inst)[0]
    for visit in np.unique(data.seq[idx]):
        members = idx[data.seq[idx] == visit]
        pv = prob[members]
        bad, fine = pv > threshold, pv < good
        if np.all(fine):
            good_visits.append(visit)
        if not np.any(bad):
            continue
        if np.all(bad):
            units.append(dict(kind='visit' if len(members) > 1 else
                              'exposure', members=members, visit=int(visit),
                              compare='neighbours'))
            continue
        for it in members[bad]:
            companions = members[fine]
            units.append(dict(kind='exposure', members=np.array([it]),
                              visit=int(visit),
                              compare='visit' if len(companions) else
                              'neighbours', companions=companions))
    return units, np.array(good_visits, dtype=int)


def _visit_values(data: RVData, values: np.ndarray, errors, idx: np.ndarray
                  ) -> Tuple[Dict[int, float], Dict[int, float]]:
    """the median value of each visit of the exposures idx, and the error
    of that median (the median error over the square root of the number)"""
    vals, errs = {}, {}
    for visit in np.unique(data.seq[idx]):
        members = idx[data.seq[idx] == visit]
        good = np.isfinite(values[members])
        if not np.any(good):
            continue
        vals[int(visit)] = float(np.median(values[members][good]))
        if errors is not None:
            errs[int(visit)] = float(np.median(errors[members][good])
                                     / np.sqrt(np.sum(good)))
    return vals, errs


def _numeric(data: RVData, key: Key, values: np.ndarray, errors,
             inst_idx: np.ndarray, prob: np.ndarray, units, good_visits,
             good: float, neighbours: int) -> Dict[str, Any]:
    """
    The deviations of a number: of every unit, and of the good exposures
    and visits (the null distributions)
    """
    fine = inst_idx[prob[inst_idx] < good]
    # within a visit: a good exposure against the other good exposures of
    #   its visit
    null_w, typ_err_w = [], None
    for visit in np.unique(data.seq[fine]):
        members = fine[data.seq[fine] == visit]
        members = members[np.isfinite(values[members])]
        if len(members) < 2:
            continue
        for it in members:
            others = members[members != it]
            null_w.append(values[it] - np.median(values[others]))
    null_w = np.array(null_w)
    scale_w = _scale(null_w) if len(null_w) >= 5 else np.nan
    if errors is not None and len(fine):
        typ_err_w = float(np.nanmedian(errors[fine]))
    # between visits: a good visit against the good visits nearest in time
    vals, verrs = _visit_values(data, values, errors, inst_idx)
    times = {int(v): float(np.median(data.time[inst_idx][
        data.seq[inst_idx] == v])) for v in vals}
    goods = [v for v in good_visits if v in vals]
    gtime = np.array([times[v] for v in goods])
    gval = np.array([vals[v] for v in goods])

    def reference(time, exclude=None):
        mask = np.ones(len(goods), dtype=bool)
        if exclude is not None:
            mask &= np.array(goods) != exclude
        if not np.any(mask):
            return np.nan
        near = np.argsort(np.abs(gtime[mask] - time))[:neighbours]
        return float(np.median(gval[mask][near]))
    null_b = np.array([vals[v] - reference(times[v], exclude=v)
                       for v in goods]) if len(goods) > 1 else np.array([])
    scale_b = _scale(null_b) if len(null_b) >= 5 else np.nan
    typ_err_b = (float(np.median([verrs[v] for v in goods]))
                 if errors is not None and goods else None)
    zs = []
    for unit in units:
        members = unit['members']
        if unit['compare'] == 'visit':
            comp = unit['companions']
            comp = comp[np.isfinite(values[comp])]
            value = float(values[members[0]])
            ref = float(np.median(values[comp])) if len(comp) else np.nan
            scale = scale_w
            if errors is not None and np.isfinite(errors[members[0]]):
                scale = np.sqrt(scale ** 2 + max(errors[members[0]] ** 2
                                                 - typ_err_w ** 2, 0.0))
        else:
            good_m = members[np.isfinite(values[members])]
            value = (float(np.median(values[good_m])) if len(good_m)
                     else np.nan)
            ref = reference(float(np.median(data.time[members])))
            scale = scale_b
            if errors is not None and len(good_m):
                err = float(np.median(errors[good_m]) / np.sqrt(len(good_m)))
                scale = np.sqrt(scale ** 2 + max(err ** 2 - typ_err_b ** 2,
                                                 0.0))
        zs.append(dict(value=value, reference=ref,
                       z=_zscore(value - ref, scale) if np.isfinite(scale)
                       else np.nan))
    return dict(units=zs,
                null={'visit': null_w / scale_w if scale_w > 0 else
                      np.zeros(len(null_w)),
                      'neighbours': null_b / scale_b if scale_b > 0 else
                      np.zeros(len(null_b))})


def _text(values: np.ndarray, inst_idx: np.ndarray, prob: np.ndarray,
          units, good: float) -> Dict[str, Any]:
    """The value of a text key for every unit, and how common it is among
    the good exposures"""
    values = np.char.strip(values.astype(str))
    fine = inst_idx[prob[inst_idx] < good]
    names, counts = np.unique(values[fine], return_counts=True)
    freq = dict(zip(names, counts / max(len(fine), 1)))
    out = []
    for unit in units:
        cats, num = np.unique(values[unit['members']], return_counts=True)
        cat = str(cats[np.argmax(num)])
        out.append(dict(value=cat, frequency=float(freq.get(cat, 0.0))))
    return dict(units=out, frequency=freq)


# =============================================================================
# The report
# =============================================================================
@dataclass
class OutlierReport:
    """
    What each outlier has that the good data do not

    units: every outlier (a whole visit, or a single exposure), in time
    order, with its keys sorted from the most to the least deviant; each key
    as a dict (key, label, value, reference, z, tail, significant,
    direction, bad; for text: value, frequency, significant).
    association: per instrument and key, what the outliers share (mean z,
    permutation p-value, Holm-adjusted p-value, significant).
    """
    name: str
    threshold: float
    alpha: float
    units: List[Dict[str, Any]] = field(default_factory=list)
    association: List[Dict[str, Any]] = field(default_factory=list)
    keys: Dict[str, List[str]] = field(default_factory=dict)
    lists: Dict[str, str] = field(default_factory=dict)
    missing: Dict[str, List[str]] = field(default_factory=dict)
    counts: Dict[str, Dict[str, int]] = field(default_factory=dict)
    zcrit: Dict[str, float] = field(default_factory=dict)

    def shared(self, inst: Optional[str] = None) -> List[Dict[str, Any]]:
        """the keys that the outliers of an instrument (or of all) share"""
        return [row for row in self.association if row['significant']
                and (inst is None or row['inst'] == inst)]

    def text(self, also: int = 3) -> str:
        """
        The report in words

        :param also: int, how many more keys to name per outlier, beyond
                     the significant ones (those with |z| > 2)
        """
        lines = [f'Outliers of {self.name} and why (p > {self.threshold}: '
                 f'an outlier; keys off beyond |z| > z_crit, Bonferroni at '
                 f'alpha = {self.alpha}, and rarer than {100 * TAIL:.0f} % '
                 f'of the good data)']
        for inst, count in self.counts.items():
            lines.append('')
            # a file without an instrument column: the instrument its keys
            #   are those of
            shown = inst if inst != 'inst' else NAMES.get(self.lists[inst],
                                                          'the series')
            lines.append(
                f'{shown}: {count["exposures"]} exposures in '
                f'{count["visits"]} visits; '
                f'{_plural(count["outliers"], "outlying exposure")} '
                f'({_plural(count["visit_units"], "whole visit")} '
                f'and {_plural(count["exposure_units"], "single exposure")}), '
                f'{count["borderline"]} borderline; '
                f'{_plural(len(self.keys.get(inst, [])), "key")} '
                f'({self.lists[inst]} list), |z| > '
                f'{self.zcrit.get(inst, np.nan):.1f} to be off')
            if self.missing.get(inst):
                lines.append(f'  not in the data: '
                             f'{", ".join(self.missing[inst])}')
            shared = self.shared(inst)
            tested = [row for row in self.association if row['inst'] == inst]
            if shared:
                lines.append('  What the outliers share (all of them against '
                             'the good data, Holm-corrected):')
                for row in shared:
                    lines.append(f'    {row["label"]}: {row["summary"]}')
            elif tested:
                lines.append('  Nothing that the outliers share, beyond '
                             'chance (Holm-corrected).')
            for unit in [u for u in self.units if u['inst'] == inst]:
                where = ('against the good exposures of its visit'
                         if unit['compare'] == 'visit' else
                         'against the good visits nearest in time')
                what = (f'visit of {unit["n"]} exposures'
                        if unit['kind'] == 'visit' else 'exposure')
                resid = ''
                if unit.get('resid') is not None and \
                        np.isfinite(unit['resid']):
                    resid = (f', residual {unit["resid"]:+.2f} m/s '
                             f'({unit["resid_sigma"]:+.1f} sigma)')
                lines.append(f'  {unit["date"]}, {what}, p = '
                             f'{unit["prob"]:.2f}{resid}, {where}:')
                sig = [key for key in unit['keys'] if key['significant']]
                rest = [key for key in unit['keys'] if not key['significant']
                        and key.get('z') is not None
                        and np.isfinite(key['z']) and abs(key['z']) > 2]
                if sig:
                    lines.append('    off: ' + '; '.join(
                        _describe(key) for key in sig))
                else:
                    lines.append('    nothing significantly off')
                if rest[:also]:
                    lines.append('    also unusual: ' + '; '.join(
                        _unusual(key, self.zcrit.get(inst, np.inf))
                        for key in rest[:also]))
        return '\n'.join(lines)

    def table(self) -> List[Dict[str, Any]]:
        """one row per outlier and key"""
        rows = []
        for unit in self.units:
            for key in unit['keys']:
                rows.append(dict(inst=unit['inst'], date=unit['date'],
                                 kind=unit['kind'], time=unit['time'],
                                 prob=unit['prob'], **{k: v for k, v in
                                                       key.items()}))
        return rows

    def to_dict(self) -> Dict[str, Any]:
        """everything, as numbers and text (for a JSON file)"""
        def clean(obj):
            if isinstance(obj, dict):
                return {str(k): clean(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple)):
                return [clean(v) for v in obj]
            if isinstance(obj, np.ndarray):
                return [clean(v) for v in obj.tolist()]
            if isinstance(obj, (np.floating, float)):
                return float(obj) if np.isfinite(obj) else None
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, np.bool_):
                return bool(obj)
            return obj
        return clean(dict(name=self.name, threshold=self.threshold,
                          alpha=self.alpha, units=self.units,
                          association=self.association, keys=self.keys,
                          lists=self.lists, missing=self.missing,
                          counts=self.counts, zcrit=self.zcrit))


def _unusual(key: Dict[str, Any], zcrit: float) -> str:
    """a key that deviates without being off, and why it is not"""
    text = f'{key["label"]} (z {key["z"]:+.1f}'
    if abs(key['z']) >= zcrit and np.isfinite(key.get('tail', np.nan)):
        text += (f', but {100 * key["tail"]:.0f} % of the good data '
                 f'deviate as much')
    return text + ')'


def _plural(count: int, word: str) -> str:
    """'1 visit', '2 visits'"""
    return f'{count} {word}' if count == 1 else f'{count} {word}s'


def _describe(key: Dict[str, Any]) -> str:
    """a key that is off, in words"""
    if key.get('text'):
        return (f'{key["label"]} = {key["value"]} ({100 * key["frequency"]:.0f}'
                f' % of the good exposures)')
    value, ref = _fmt_pair(key['value'], key['reference'])
    return f'{key["label"]} {value} against {ref} (z {key["z"]:+.1f})'


def _fmt_pair(value: float, ref: float) -> Tuple[str, str]:
    """two numbers with three significant digits, or as many more as it
    takes to tell them apart"""
    if value is None or ref is None or not (np.isfinite(value)
                                            and np.isfinite(ref)):
        return str(value), str(ref)
    for digits in range(3, 10):
        one, two = f'{value:.{digits}g}', f'{ref:.{digits}g}'
        if one != two:
            return one, two
    return one, two


def explain(fit: Any, keys: Union[str, Sequence[Union[str, Key]]] = 'auto',
            prob: Optional[np.ndarray] = None, threshold: float = THRESHOLD,
            good: float = GOOD, alpha: float = 0.01, neighbours: int = 10,
            npermutations: int = 4000, seed: int = 1) -> OutlierReport:
    """
    Why each outlier of a series is an outlier: which recorded quantities
    are significantly off for it, and which the outliers share

    :param fit: a FitResult of koloa with the mixture likelihood (its
                series, outlier probabilities and residuals are used), or an
                RVData (then prob is required)
    :param keys: 'auto', the name of a list ('spirou', 'nirps', 'dace'), or
                 a list of names and Key objects (see resolve_keys)
    :param prob: np.ndarray or None, the outlier probability of every
                 exposure (from the fit when None)
    :param threshold: float, an exposure is an outlier above it
    :param good: float, and good below it (the reference)
    :param alpha: float, the false-alarm rate of a key, for all the keys of
                  an outlier together (Bonferroni), and of the association
                  (Holm)
    :param neighbours: int, the good visits a visit is compared with
    :param npermutations: int, the random draws of the association test
    :param seed: int, the seed of those draws

    :return: OutlierReport
    """
    if isinstance(fit, RVData):
        data, resid, noise = fit, None, None
        if prob is None:
            raise ValueError('explain(data) needs prob, the outlier '
                             'probability of every exposure')
    else:
        data = fit.model.data
        prob = fit.outlier_prob if prob is None else prob
        resid = fit.residuals()
        diag, seq_var = fit.model.noise(fit.theta)
        noise = (diag, seq_var)
    prob = np.asarray(prob, dtype=float)
    if len(prob) != data.n:
        raise ValueError(f'{len(prob)} probabilities for {data.n} exposures')
    rng = np.random.default_rng(seed)
    report = OutlierReport(name=data.name, threshold=threshold, alpha=alpha)
    for inst in data.instruments:
        inst_idx = np.where(data.inst == inst)[0]
        lname, have, lacking = resolve_keys(data, inst, keys)
        report.lists[inst], report.missing[inst] = lname, lacking
        report.keys[inst] = [key.name for key in have]
        units, good_visits = _units(data, prob, inst, threshold, good)
        pin = prob[inst_idx]
        report.counts[inst] = dict(
            exposures=int(len(inst_idx)),
            visits=int(len(np.unique(data.seq[inst_idx]))),
            outliers=int(np.sum(pin > threshold)),
            borderline=int(np.sum((pin >= good) & (pin <= threshold))),
            visit_units=sum(u['kind'] == 'visit' for u in units),
            exposure_units=sum(u['kind'] == 'exposure' for u in units))
        numeric = [key for key in have
                   if not _is_text(key.values(data)[0])]
        mtests = max(len(numeric), 1)
        zcrit = float(stats.norm.isf(alpha / (2 * mtests)))
        report.zcrit[inst] = zcrit
        # every unit: its time, probability and residual
        for unit in units:
            members = unit['members']
            unit['inst'] = inst
            unit['n'] = int(len(members))
            unit['time'] = float(np.median(data.time[members]))
            unit['date'] = _date(unit['time'], unit['kind'] == 'exposure')
            unit['prob'] = float(np.mean(prob[members]))
            unit['index'] = [int(it) for it in members]
            unit['keys'] = []
            unit['resid'] = unit['resid_sigma'] = None
            if resid is not None:
                weight = 1.0 / noise[0][members]
                mean = float(np.sum(weight * resid[members]) / np.sum(weight))
                var = 1.0 / np.sum(weight)
                seq_var = noise[1]
                if np.ndim(seq_var):
                    seq_var = float(seq_var[unit['visit']])
                var += float(seq_var) if len(members) > 1 or \
                    unit['kind'] == 'visit' else 0.0
                unit['resid'] = mean
                unit['resid_sigma'] = mean / np.sqrt(var)
        nulls = {}
        for key in have:
            values, errors = key.values(data)
            if _is_text(values):
                res = _text(values, inst_idx, prob, units, good)
                for unit, row in zip(units, res['units']):
                    unit['keys'].append(dict(
                        key=key.name, label=key.label, text=True,
                        value=row['value'], frequency=row['frequency'],
                        significant=row['frequency'] < RARE, z=None,
                        bad=key.bad))
                nulls[key.name] = ('text', res['frequency'])
                continue
            res = _numeric(data, key, values, errors, inst_idx, prob, units,
                           good_visits, good, neighbours)
            nulls[key.name] = ('number', res['null'])
            for unit, row in zip(units, res['units']):
                null = res['null'][unit['compare']]
                zval = row['z']
                tail = (float(np.mean(np.abs(null) >= abs(zval)))
                        if len(null) and not np.isnan(zval) else np.nan)
                # an infinite z (a key that never varies among the good
                #   data, and changes here) is off too
                significant = bool(not np.isnan(zval) and abs(zval) >= zcrit
                                   and (not len(null) or tail <= TAIL))
                unit['keys'].append(dict(
                    key=key.name, label=key.label, text=False,
                    value=row['value'], reference=row['reference'], z=zval,
                    tail=tail, significant=significant,
                    direction=('high' if zval > 0 else 'low')
                    if np.isfinite(zval) and zval != 0 else '',
                    bad=key.bad))
        for unit in units:
            unit['keys'].sort(key=lambda row: (
                not row['significant'],
                -(abs(row['z']) if row.get('z') is not None
                  and np.isfinite(row['z']) else
                  (1 - row.get('frequency', 1.0)) if row.get('text') else 0)))
        report.units.extend(units)
        # what the outliers share: every key, all the units together
        rows = _association(data, have, units, nulls, prob, inst, inst_idx,
                            threshold, good, npermutations, rng)
        if rows:
            adj = _holm([row['p'] for row in rows])
            for row, padj in zip(rows, adj):
                row['p_holm'] = float(padj)
                row['significant'] = bool(padj < 0.05)
            rows.sort(key=lambda row: row['p'])
            report.association.extend(rows)
    report.units.sort(key=lambda unit: unit['time'])
    return report


def _association(data, have, units, nulls, prob, inst, inst_idx, threshold,
                 good, npermutations, rng) -> List[Dict[str, Any]]:
    """
    For every key, do the outliers of an instrument deviate together? A
    number: the mean z of the units, each clipped at +-3 so that one
    extreme outlier cannot make a trend alone, against the same mean for as
    many good units drawn at random (the within-visit and between-visit
    kinds in the same numbers). Text: is a value more common among the
    outliers than among the good exposures (the binomial tail of as many
    outliers having it, its frequency among the good exposures counted
    with half an exposure more; the most over-represented value, times the
    number of values)?
    """
    if len(units) < 2:
        return []
    rows = []
    kinds = [unit['compare'] for unit in units]
    for key in have:
        kind, null = nulls[key.name]
        if kind == 'text':
            freq = null
            fine = inst_idx[prob[inst_idx] < good]
            cats = [unit_key['value'] for unit in units
                    for unit_key in unit['keys'] if unit_key['key'] == key.name]
            best = (1.0, None, 0.0, 0.0, 0)
            for cat in set(cats):
                count = cats.count(cat)
                share = ((freq.get(cat, 0.0) * len(fine) + 0.5)
                         / (len(fine) + 1))
                pval = float(stats.binom.sf(count - 1, len(cats), share))
                if pval < best[0]:
                    best = (pval, cat, count / len(cats),
                            freq.get(cat, 0.0), count)
            if best[1] is None:
                continue
            pval = min(best[0] * max(len(freq), 1), 1.0)
            rows.append(dict(
                inst=inst, key=key.name, label=key.label, text=True, p=pval,
                value=best[1], frac_outliers=best[2], frac_good=best[3],
                summary=(f'{best[1]} for {best[4]} of the {len(cats)} '
                         f'outliers, {100 * best[3]:.0f} % of the good '
                         f'exposures (p = {pval:.2g})')))
            continue
        zs = np.array([unit_key['z'] for unit in units
                       for unit_key in unit['keys']
                       if unit_key['key'] == key.name])
        finite = np.isfinite(zs)
        if np.sum(finite) < 2:
            continue
        zs_c = np.clip(zs, -CLIP, CLIP)
        observed = float(np.mean(zs_c[finite]))
        used = [kinds[it] for it in np.where(finite)[0]]
        pools = {name: np.clip(null[name], -CLIP, CLIP) for name in null}
        if any(len(pools[name]) == 0 for name in set(used)):
            continue
        draws = np.zeros(npermutations)
        for name in set(used):
            count = used.count(name)
            pool = pools[name]
            pick = rng.integers(0, len(pool), size=(npermutations, count))
            draws += pool[pick].sum(axis=1)
        draws /= len(used)
        pval = float((1 + np.sum(np.abs(draws) >= abs(observed)))
                     / (1 + npermutations))
        direction = 'higher' if observed > 0 else 'lower'
        sign = 1 if observed > 0 else -1
        agree = int(np.sum(sign * zs[finite] > 2))
        rows.append(dict(
            inst=inst, key=key.name, label=key.label, text=False, p=pval,
            mean_z=observed, direction=direction, nunits=int(np.sum(finite)),
            nagree=agree,
            summary=(f'{direction} in the outliers ({agree} of '
                     f'{int(np.sum(finite))} beyond 2 sigma; mean z, clipped '
                     f'at {CLIP:.0f}, {observed:+.1f}; p = {pval:.2g})')))
    return rows


# =============================================================================
# End of code
# =============================================================================

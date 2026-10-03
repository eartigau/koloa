#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Everything koloa can say about a star, from one velocity file.

`detailed_analysis` starts from an LBL .rdb (or a csv, or an RVData) and:

1. reads the series; a file without an instrument column is named after
   the instrument its keywords are those of (SPIRou, NIRPS);
2. finds who the star is (CDS Sesame: its main name and aliases) and what
   is known of it: its planets in the NASA Exoplanet Archive, with every
   published solution (koloa.archive);
3. asks DACE for more velocities: the public ones of every other
   instrument (one instrument per era: HARPS03, HARPS15, ESPRESSO19...),
   merged with the file (koloa.data.merge); DACE answers from some networks
   only, and the analysis goes on without it; adds the published
   velocities given (literature=[files or series]) and those published
   with the known planets, found on VizieR from the papers of their
   solutions (koloa.literature; kept in outdir/literature);
4. fits the noise of every instrument without planets (a white and a visit
   jitter each), inflates the errors of each instrument to it, and runs the
   outlier-aware FIP with up to kmax signals (koloa.fip.inflate_to_fit,
   koloa.oafip);
5. fits every interval with FIP < 1 % together, with the known planets
   (tested at their periods whether the FIP finds them or not) and the
   periods given (each period free within 2 %; every exposure and visit may
   be an outlier; by MCMC with mcmc=True), and compares each with the known
   planets (every published
   K, the most recent first, and the published ephemeris carried to the
   fitted conjunction);
6. runs the FIP again, the errors inflated to the fit with the planets
   (the first pass had the planets in its noise);
7. looks at the activity indicators: the peaks of their outlier-aware
   periodograms, and whether a signal of the velocities sits at one of them
   (or at twice or half of it), or at the rotation period of the star the
   archive knows (or at its half or third);
   then fits the signals with a GP of the activity (the rotation period
   of the archive as its prior, a free period otherwise; chromatic with
   several instruments): the likelihood each signal adds over the GP, the
   periodograms whitened by the GP and their false-alarm probabilities,
   and the whole series with the GP (koloa.gpcheck);
8. puts every signal to the duck test;
9. explains every outlier (koloa.outliers: the header keywords, the
   indicators and the error bars that are off for it);
10. writes a text report, a JSON summary and the figures;
11. writes it all as a LaTeX report, compiled to <star>_report.pdf when
    pdflatex is there (koloa.latex).

    from koloa.detailed import detailed_analysis
    out = detailed_analysis('gl687_spirou.rdb', outdir='gl687')

    koloa gl687_spirou.rdb --detailed --outdir gl687

Created on 2026-09-29

@author: artigau
"""
import json
import os
import time
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np

from koloa import literature as klit
from koloa import kepler
from koloa import plotting as kplot
from koloa.archive import MATCH, compare, conjunction, known_planets, resolve
from koloa.dace import names as dace_names
from koloa.data import RVData, instrument_names, merge, robust_std
from koloa.diagnostics import duck_test
from koloa.fip import inflate_to_fit, oafip
from koloa.fit import RVModel, mcmc_orbits
from koloa.gpcheck import (figure_fold, figure_sequence, figure_whitened,
                           gp_signals)
from koloa.log import log, step
from koloa.outliers import NAMES, explain, instrument_list
from koloa.periodogram import find_peaks, frequency_grid, oap
from koloa.utils import blas_threads

# =============================================================================
# Define variables
# =============================================================================
#: a signal of the FIP is detected below this
THRESHOLD = 0.01


def _decisive(peak: Dict[str, Any]) -> float:
    """planet or no planet: the FIP of a peak's period OR any of its
    aliases (the FIP of the period alone for a FIP without a family)"""
    fam = peak.get('family_fip')
    return float(peak['fip'] if fam is None else fam)


def _found(res) -> List[float]:
    """the periods whose family FIP is below THRESHOLD, the most probable
    alias of each family only (an alias of a period already found is the
    same signal, not a second planet)"""
    from koloa.aliases import same_family
    out: List[float] = []
    # the width of an interval, 1/T (set by _fip)
    width = res.settings.get('width', 1.0 / 365.25)
    for peak in sorted(res.peaks, key=lambda pk: pk['fip']):
        if _decisive(peak) >= THRESHOLD:
            continue
        if any(same_family(per, peak['period'], width) for per in out):
            continue
        out.append(float(peak['period']))
    return out
#: the indicators looked at, in this order, when the series has them
INDICATORS = ('DTEMP3500', 'DTEMP', 'fwhm', 'd2v', 'dW', 'contrast', 'CRX',
              'bis', 'rhk', 'smw', 'halpha')
#: a signal sits at an activity period when within this fraction of it
ACTIVITY_MATCH = 0.03
#: a signal at the rotation period divided by one of these may be activity
ROTATION_HARMONICS = (1, 2, 3)
#: an indicator value this many robust sigma off is left out of its
#: periodogram
GROSS = 25.0


# =============================================================================
# Define functions
# =============================================================================
def _read(source: Union[str, RVData], name: Optional[str],
          label: Optional[str] = None) -> RVData:
    """the series of a file (or the series given), its instrument named:
    as given (label), else from the file's own columns (koloa.data.
    instrument_names: HARPS03 or HARPS15, ESPRESSO18 or 19, NIRPS, SPIRou,
    HARPN...); a CSV written by koloa.dace.fetch is read by
    koloa.dace.rvdata (one instrument per era)"""
    if isinstance(source, RVData):
        data = source
    else:
        with open(source) as handle:
            header = next((line for line in handle
                           if line.strip() and not line.startswith('#')), '')
        if 'instrument_name' in header and 'rv_err' in header:
            from koloa.dace import rvdata
            data = rvdata(source, name=name)
        else:
            data = RVData.from_csv(source, name=name)
    if name:
        data.name = name
    names = None
    if label and str(label).lower() != 'auto':
        names = np.full(data.n, str(label))
    elif data.instruments == ['inst']:
        names = instrument_names(data)
        if names is None:
            lname, _ = instrument_list(data, 'inst')
            if NAMES.get(lname):
                names = np.full(data.n, NAMES[lname])
    return data if names is None else relabel(data, names)


def read_files(files: Sequence[Any],
               instruments: Optional[Sequence[Optional[str]]] = None,
               name: Optional[str] = None) -> List[RVData]:
    """
    Several files (or series), each its instruments: as given, else from
    the file's own columns; an instrument already in an earlier file is
    named <inst>_<k> (k, the file's rank), each its own offset

    :return: list of RVData, one per file
    """
    labels = list(instruments or []) + [None] * len(files)
    parts: List[RVData] = []
    for ifile, (src, label) in enumerate(zip(files, labels)):
        part = _read(src, name, label)
        taken = {inst for prev in parts for inst in prev.instruments}
        clash = [inst for inst in part.instruments if inst in taken]
        if clash:
            new = np.array([f'{inst}_{ifile + 1}' if inst in clash
                            else inst for inst in part.inst])
            log(f'{_label(src)}: {", ".join(clash)} is in another file '
                f'too: named with _{ifile + 1} (instruments= to name them)',
                'warn')
            part = relabel(part, new)
        parts.append(part)
        log(f'{_label(src)}: {part.n} exposures, ' + ', '.join(
            f'{inst} {np.sum(part.inst == inst)}'
            for inst in part.instruments), 'value')
    return parts


def _label(source: Any) -> str:
    """a file by its name, a series by its own"""
    return (os.path.basename(source) if isinstance(source, str)
            else getattr(source, 'name', 'series'))


def relabel(data: RVData, names: np.ndarray) -> RVData:
    """the same series with other instruments (each its own zero point)"""
    vrad = data.rv + np.array([data.zero_point[str(inst)]
                               for inst in data.inst])
    return RVData(data.time, vrad, data.err, inst=np.asarray(names),
                  indicators=dict(data.indicators), name=data.name,
                  sequence_gap=data.sequence_gap, meta=dict(data.meta))


def _target(data: RVData, target: Optional[str],
            files: Sequence[Any] = ()) -> Optional[str]:
    """the name of the star: given; else its APERO name (the OBJECT column
    of the file, or its name, lbl_<OBJECT>_<TEMPLATE>.rdb, in APERO's
    database of names: koloa.apero_names) as SIMBAD knows it; else the
    OBJECT column or the series (a guess: SIMBAD may not know it, and the
    archives then add nothing)"""
    if target:
        return target
    from koloa import apero_names
    guess = data.name
    for col in ('OBJECT', 'object', 'target', 'TARGET'):
        if col in data.meta and data.meta[col].dtype.kind in 'USO':
            names, counts = np.unique(data.meta[col], return_counts=True)
            best = str(names[np.argmax(counts)]).strip()
            if best:
                guess = best
                break
    tries = [guess] + [name for path in files if isinstance(path, str)
                       for name in apero_names.name_candidates(path)]
    for name in tries:
        try:
            entry = apero_names.lookup(name)
        except (ImportError, OSError, ValueError) as err:
            log(f'APERO names: not read ({err})', 'warn')
            break
        if entry is not None:
            star = apero_names.simbad_target(entry)
            log(f'no SIMBAD name given: {name!r}, from the file, is '
                f'{entry["apero"]} in APERO\'s names, {star!r} for SIMBAD',
                'value')
            return star
    log(f'no SIMBAD name given: {guess!r}, from the file; give it '
        f'(target=, --target) for the archives to find the star', 'warn')
    return guess


def _carmenes(ident: Dict[str, Any], star: str, folder: str,
              data: Optional[RVData], refresh: bool = False):
    """
    The velocities of CARMENES DR1 that the series does not have (koloa.gather:
    corrected for the nightly zero points), and a note for the sources

    :return: tuple, RVData or None, and the note
    """
    from koloa.gather import carmenes_rv, carmenes_star
    try:
        cstar = carmenes_star(ident['ra'], ident['dec'])
        if cstar is None:
            return None, 'not in CARMENES DR1'
        more = carmenes_rv(cstar, star, os.path.join(folder, 'carmenes'),
                           refresh=refresh)
    except OSError as err:
        log(f'CARMENES DR1: {err}', 'warn')
        return None, f'unreachable ({err})'
    if more is None:
        return None, (f'{cstar["carmenes_id"]}: no velocity corrected for the '
                      f'nightly zero points')
    if data is not None:
        # set apart from the file's own CARMENES velocities, its spectra
        #   not counted twice
        more = distinct(more, data.instruments, data.time, 'DR1')
        if more is None:
            return None, (f'{cstar["carmenes_id"]}: nothing the file does '
                          f'not have')
    return more, (f'{cstar["carmenes_id"]}, velocities corrected for the '
                  f'nightly zero points (Ribas et al. 2023)')


def fetch_dace(names: Sequence[str], folder: str, exclude: Sequence[str] = (),
               times: Optional[np.ndarray] = None, refresh: bool = False
               ) -> Optional[RVData]:
    """
    The public velocities of a star on DACE, trying its names in turn

    :param names: list of str, the names to try
    :param folder: str, where the CSV is kept
    :param exclude: list of str, the instruments of the file: DACE's own
                    of the same name are set apart (distinct)
    :param times: np.ndarray or None, the times of the file's exposures: a
                  DACE exposure within a minute of one is the same spectrum,
                  and is left out
    :param refresh: bool, ask DACE even when a CSV is there

    :return: RVData or None (DACE unreachable, or nothing public)
    """
    from koloa.dace import fetch, rvdata
    for name in names:
        path = os.path.join(folder, f'{name}_dace.csv')
        try:
            fetch(name, path, refresh=refresh)
        except RuntimeError as err:
            text = str(err)
            log(f'DACE, {name}: {text.split(";")[0]}', 'warn')
            if 'could not be reached' in text:
                return None
            continue
        try:
            data = rvdata(path, name=name)
        except (ValueError, IndexError) as err:
            log(f'DACE, {name}: {err}', 'warn')
            continue
        return distinct(data, exclude, times, 'DACE')
    return None


def distinct(more: RVData, theirs: Sequence[str] = (),
             times: Optional[np.ndarray] = None, source: str = 'DACE'
             ) -> Optional[RVData]:
    """
    The velocities of an archive set apart from those of the file: an
    instrument the file has too (NIRPS in an LBL file and on DACE, or an era
    of it: HARPS and HARPS15) keeps its own name, <inst>_<source>, so that
    each pipeline has its offset and can be left out on its own; and the
    exposures that are the file's spectra (the same time, within a minute)
    are left out, not counted twice

    :param more: RVData, the archive's series
    :param theirs: list of str, the instruments of the file
    :param times: np.ndarray or None, the times of the file's exposures
    :param source: str, the archive (DACE, DR1...), for the names

    :return: RVData or None (nothing the file does not have)
    """
    mask = np.ones(more.n, dtype=bool)
    if times is not None and len(times):
        nearest = np.min(np.abs(more.time[:, None]
                                - np.asarray(times)[None, :]), axis=1)
        same = nearest < 1.0 / 1440
        if np.any(same):
            log(f'{source}: {int(np.sum(same))} exposures that the file '
                f'already has (the same spectra), left out', 'warn')
        mask &= ~same
    if not np.any(mask):
        log(f'{source}: nothing that the file does not have', 'warn')
        return None
    more = more.select(mask)
    for inst in list(more.instruments):
        if any(inst.upper().startswith(ex.upper())
               or ex.upper().startswith(inst.upper()) for ex in theirs):
            new = f'{inst}_{source}'
            more.inst = np.where(more.inst == inst, new, more.inst)
            more.zero_point[new] = more.zero_point.pop(inst)
            log(f'{source}: {inst} is in the file too, from another '
                f'pipeline: named {new}, with its own offset', 'value')
    return more


def _gp_spec(fip_gp: Any, known: Dict[str, Any]) -> Any:
    """
    The GP of the FIP: 'auto' is a local GP, with a rotation GP beside it
    when the archive knows the rotation period (its prior, +- 10 %); 'sho'
    an SHO at the rotation period (held) and at its half
    """
    if fip_gp in (None, False):
        return None
    if fip_gp == 'banded':
        return 'banded'
    prot = (known.get('star') or {}).get('rotation')
    if fip_gp == 'sho':
        if prot:
            # a rotation that can be trusted: its period held
            return [dict(kind='sho', period=float(prot))]
        log('an SHO GP needs the rotation period (rotation=, --rotation): '
            'the FIP by bands instead', 'warn')
        return 'banded'
    if fip_gp is not True and fip_gp != 'auto':
        return fip_gp
    if prot:
        return ['local', dict(kind='rotation', period=dict(
            mu=float(np.log(float(prot))), sd=0.1))]
    return 'local'


def _fip_gp_text(res) -> str:
    """the GP of a FIP, in words"""
    if res.settings.get('bands'):
        from koloa.bandfip import gp_text
        return 'GP by band (local, its shortest scale): ' + gp_text(
            res.settings['bands'])
    comps = res.settings.get('gp') or []
    if not comps:
        return 'no GP'
    from koloa.gpbasis import label
    kinds = ' + '.join(label(comp) for comp in comps)
    vals = '; '.join(f'{key.split("_", 1)[1]} {mid:.3g} ({low:.3g} to '
                     f'{high:.3g})' for key, (mid, low, high)
                     in res.gp_summary().items())
    return f'GP {kinds}: {vals}'


def _detection_map(data: RVData, fit: Any, nplanet: int, pmin: float,
                   ninj: int, workers: int, seed: int,
                   mstar: Optional[float],
                   prot: Optional[float] = None) -> Dict[str, Any]:
    """
    Which planets the series could have found: circular planets injected
    into it, the signals found taken out first, on a grid of periods (from
    pmin to the baseline) and semi-amplitudes (0.3 to 15 times the median
    error), each looked for by a blind outlier-aware search (koloa.
    completeness, test 'search', a false-alarm probability of 1 %, the
    threshold of each period bin calibrated on injections at K = 0, a
    slope fitted beside each sinusoid)

    The activity is left in the series, so a period bin that holds the
    rotation period (or its half) finds the activity's own peak: its
    recoveries are not a planet's, and the bin is marked unreliable.

    :return: dict, map (RecoveryMap), K50 and K90 per period bin [m/s],
             msini50 and msini90 [Earth masses] when the star's mass is
             known, unreliable (per bin: holds P_rot or P_rot/2)
    """
    from koloa.completeness import recovery_map
    series = data.select(np.ones(data.n, dtype=bool))
    if nplanet:
        series.rv = data.rv - sum(fit.model.planet_rv(fit.theta, ip)
                                  for ip in range(nplanet))
    med = float(np.median(data.err))
    period_edges = np.geomspace(pmin, max(data.baseline, 2 * pmin), 13)
    amp_edges = np.geomspace(0.3 * med, 15 * med, 11)
    unit = 'point' if data.nseq == data.n else 'sequence'
    rmap = recovery_map(series, period_edges, amp_edges, ninj=ninj,
                        test='search', outliers=True, unit=unit, fap=0.01,
                        nnull=max(20, 10 * ninj), seed=seed, workers=workers,
                        trend=1,
                        calibrate=True, quiet=True)
    k50, k90 = rmap.k_at(0.5), rmap.k_at(0.9)
    edges = rmap.period_edges
    unreliable = np.zeros(len(rmap.periods), dtype=bool)
    for per in ([prot, prot / 2] if prot else []):
        unreliable |= (edges[:-1] <= per) & (per < edges[1:])
    out = dict(map=rmap, K50=k50, K90=k90, msini50=None, msini90=None,
               unreliable=unreliable, prot=prot)
    if mstar:
        out['msini50'] = np.array([
            kepler.minimum_mass(amp, per, 0.0, mstar) if np.isfinite(amp)
            else np.nan for amp, per in zip(k50, rmap.periods)])
        out['msini90'] = np.array([
            kepler.minimum_mass(amp, per, 0.0, mstar) if np.isfinite(amp)
            else np.nan for amp, per in zip(k90, rmap.periods)])
    log('detection map: ' + rmap.summary().split('\n')[0], 'value')
    return out


def _fip_detection_map(data: RVData, fit: Any, nplanet: int, pmin: float,
                       gp: Any, workers: int, seed: int,
                       mstar: Optional[float],
                       prot: Optional[float] = None) -> Dict[str, Any]:
    """
    Which planets the series could have found, by the rule that decides
    (koloa.fipmap): planets injected into the series (the signals found
    taken out), each found when the FIP with the same GP, of its period or
    any of its aliases, is below 1 %; per period band, the K found 50 and
    90 % of the time, by adaptive rounds of injections

    :return: dict, from koloa.fipmap.fip_map, with msini50 and msini90
             [Earth masses] when the star's mass is known, prot, and
             unreliable (per band: holds P_rot or P_rot/2, where the GP of
             the rotation competes with a planet)
    """
    from koloa.fipmap import fip_map
    series = data.select(np.ones(data.n, dtype=bool))
    if nplanet:
        series.rv = data.rv - sum(fit.model.planet_rv(fit.theta, ip)
                                  for ip in range(nplanet))
    if isinstance(gp, list) and gp and isinstance(gp[0], dict) \
            and 'low' in gp[0]:
        # the FIP by bands: each planet with the GP of its band
        from koloa.bandfip import spec
        bands = gp

        def gp(period, bands=bands):
            for band in bands:
                if band['low'] <= period < band['high']:
                    return spec(band['gp'], data.baseline)
            edge = bands[0] if period >= bands[0]['high'] else bands[-1]
            return spec(edge['gp'], data.baseline)
    out = fip_map(series, pmin, gp=gp, workers=workers, seed=seed)
    _map_masses(out, mstar)
    edges = out['period_edges']
    unreliable = np.zeros(len(out['periods']), dtype=bool)
    for per in ([prot, prot / 2] if prot else []):
        unreliable |= (edges[:-1] <= per) & (per < edges[1:])
    out.update(unreliable=unreliable, prot=prot)
    ok = np.isfinite(out['K90'])
    if ok.any():
        log(f'detection map (FIP): K at 90 % from {np.min(out["K90"][ok]):.2f}'
            f' to {np.max(out["K90"][ok]):.2f} m/s' + (
                f'; not reached in {np.sum(~ok)} band(s)' if (~ok).any()
                else ''), 'value')
    return out


def _map_masses(dmap: Dict[str, Any], mstar: Optional[float]):
    """m sin i at the levels of a FIP map (NaN where not reached)"""
    dmap['msini50'] = dmap['msini90'] = None
    if mstar:
        for key in ('50', '90'):
            dmap[f'msini{key}'] = np.array([
                kepler.minimum_mass(amp, per, 0.0, mstar)
                if np.isfinite(amp) else np.nan
                for amp, per in zip(dmap[f'K{key}'], dmap['periods'])])


def _map_lines(dmap: Dict[str, Any]) -> List[str]:
    """the detection map in words: per period bin, the K (and m sin i)
    recovered 50 and 90 % of the time"""
    if dmap.get('kind') == 'fip':
        edges = dmap['period_edges']
        lines = [f'By the FIP with the GP of the analysis (a planet found when '
                 f'its period or any of its aliases has FIP < '
                 f'{dmap["threshold"]:g}): {len(dmap["injections"])} '
                 f'planets injected, adaptively, per period band']
        def level(ib, key):
            note = (dmap.get(f'notes{key}') or [None] * 99)[ib]
            if note:
                return f'{note} (not reached)'
            rng = dmap[f'K{key}_range'][ib]
            return (f'{dmap[f"K{key}"][ib]:5.2f} ({rng[0]:.2f} to '
                    f'{rng[1]:.2f})')
        for ib in range(len(dmap['periods'])):
            line = (f'  P {edges[ib]:7.2f} to {edges[ib + 1]:7.2f} d: '
                    f'K(50 %) = {level(ib, "50")}, K(90 %) = '
                    f'{level(ib, "90")} m/s')
            if dmap['msini50'] is not None and np.isfinite(
                    dmap['msini50'][ib]):
                line += (f'; m sin i {dmap["msini50"][ib]:5.2f} and '
                         f'{dmap["msini90"][ib]:5.2f} Me')
            if dmap['unreliable'][ib]:
                line += (f'  (holds the rotation, {dmap["prot"]:.0f} d, or '
                         f'its half: the activity is there)')
            lines.append(line)
        return lines
    rmap = dmap['map']
    lines = [rmap.summary().split('\n')[0]]
    if rmap.nnull:
        lines.append(f'  false alarms at K = 0: {rmap.false_alarm:.3f}')
    for ip, per in enumerate(rmap.periods):
        line = (f'  P = {per:8.2f} d: K(50 %) = {dmap["K50"][ip]:6.2f}, '
                f'K(90 %) = {dmap["K90"][ip]:6.2f} m/s')
        if dmap['msini50'] is not None:
            line += (f'; m sin i {dmap["msini50"][ip]:6.2f} and '
                     f'{dmap["msini90"][ip]:6.2f} Me')
        if dmap['unreliable'][ip]:
            line += (f'  (holds the rotation, {dmap["prot"]:.0f} d, or its '
                     f'half: the activity is found, not a planet)')
        lines.append(line)
    return lines


def _tois(toi: Any, ident: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """the TOIs asked for (True: all of the star's but the false positives;
    or their numbers), with their ephemerides from TESS"""
    if not toi:
        return []
    from koloa.archive import TOI_NOT_PLANETS, tois
    if not ident or not ident.get('tic'):
        log('TOI: SIMBAD gave no TIC number for the star: no TOI', 'warn')
        return []
    found = tois(ident['tic'])
    if toi is True or (isinstance(toi, str) and toi.lower() == 'all'):
        chosen = [item for item in found
                  if item['disposition'] not in TOI_NOT_PLANETS]
    else:
        names = [toi] if isinstance(toi, (str, float, int)) else list(toi)
        wanted = {f'{float(str(name).upper().replace("TOI", "").strip(" -")):.2f}'
                  for name in names}
        chosen = [item for item in found if item['toi'] in wanted]
        for name in sorted(wanted - {item['toi'] for item in chosen}):
            log(f'TOI-{name}: not a TOI of {ident["tic"]}', 'warn')
    chosen = [item for item in chosen if item['tc'] is not None]
    if not chosen:
        log(f'TOI: none of {ident["tic"]} to fit', 'warn')
    for item in chosen:
        log(f'TOI-{item["toi"]} ({item["disposition"]}): P = {item["P"]:.7f} '
            f'+- {item["P_err"] or 0:.1e} d, transit at {item["tc"]:.5f} '
            f'+- {item["tc_err"] or 0:.1e}: fitted with these priors',
            'value')
    return chosen


def _planet(period: float, toi: Optional[Dict[str, Any]] = None
            ) -> Dict[str, Any]:
    """an orbit to fit: its period free within 2 %, or a TOI's, with the
    ephemeris of TESS as priors (P and the time of a transit)"""
    if toi is None:
        return dict(period=period, period_range=(0.98 * period,
                                                 1.02 * period))
    return dict(period=toi['P'], period_err=toi['P_err'] or 1e-4 * toi['P'],
                tc=toi['tc'], tc_err=toi['tc_err'] or 0.01,
                period_range=(0.98 * toi['P'], 1.02 * toi['P']),
                toi=toi['toi'])


def _acceleration(fit, trend: int, sampled: bool) -> Optional[Dict[str, Any]]:
    """the acceleration of the star (dv/dt) a fit measured, and its change
    when fitted, with their errors (koloa.secular.acceleration: the
    percentiles of the chain, or the Laplace covariance of the maximum a
    posteriori, every other parameter marginalised)"""
    if trend < 1:
        return None
    from koloa.secular import acceleration
    try:
        acc = acceleration(fit)
    except (ValueError, np.linalg.LinAlgError) as err:
        log(f'the acceleration of the star: {err}', 'warn')
        return None
    acc.pop('draws', None)
    acc['errors'] = ('the posterior (MCMC), 16th to 84th percentiles'
                     if sampled else 'the Laplace covariance of the '
                     'maximum a posteriori')
    return acc


def _accel_lines(acc: Optional[Dict[str, Any]]) -> List[str]:
    """the acceleration in words, with its errors and significance"""
    if not acc:
        return []
    out = []
    for key, what, unit in (('accel', 'acceleration of the star, dv/dt',
                             'm/s/yr'),
                            ('jerk', 'its change, d2v/dt2', 'm/s/yr^2')):
        if key not in acc:
            continue
        val, low, high = acc[key]
        sig = abs(val) / max(0.5 * (low + high), 1e-30)
        err = (f'+- {low:.3g}' if abs(low - high) < 0.05 * max(low, high)
               else f'-{low:.3g} +{high:.3g}')
        out.append(f'{what} = {val:+.3g} {err} {unit} ({sig:.1f} sigma)')
    if out:
        out[0] += (f', at rjd {acc["tref"]:.1f}; errors from '
                   f'{acc["errors"]}')
    return out


def _fip(data, fit, kmax, nsweep, nburn, seed, label, gp=None,
         decided=None, trend=1):
    """the FIP of a series, the errors of each instrument inflated to a fit
    (and a GP of the activity inside it; 'banded': the FIP by period bands
    of koloa.bandfip, decided from the top down or as `decided` says)"""
    inflated, info = inflate_to_fit(fit)
    log(f'{label}: errors inflated by ' + ', '.join(
        f'{inst} {val:.2f}' for inst, val in info['inflation'].items())
        + ' m/s', 'value')
    with blas_threads(1):
        if gp == 'banded':
            from koloa.bandfip import as_result, banded_fip
            band = banded_fip(inflated, kmax=kmax, nsweep=nsweep, nburn=nburn,
                              nchains=2, seed=seed, decided=decided,
                              label=label, trend=trend)
            res = as_result(band)
            info['banded'] = band['decided']
        else:
            res = oafip(inflated, kmax=kmax, outliers='both', nsweep=nsweep,
                        nburn=nburn, nchains=2, seed=seed, progress=False,
                        gp=gp, label=label, trend=trend)
    res.settings['width'] = 1 / inflated.baseline
    if gp:
        log(f'{label}: {_fip_gp_text(res)}', 'value')
    found = [pk for pk in res.peaks if _decisive(pk) < THRESHOLD]
    log(f'{label}: P(k) = ' + ', '.join(f'{val:.2f}' for val in res.pk)
        + '; FIP (the period or any alias) < 1 %: '
        + (', '.join(f'{pk["period"]:.4f} d ({_decisive(pk):.1e}; alone '
                     f'{pk["fip"]:.1e})' for pk in found) or 'none'),
        'value')
    return res, info


def _indicators(data: RVData, pmin: float, pmax: Optional[float],
                names: Union[str, Sequence[str]] = 'auto',
                npeaks: int = 3) -> Dict[str, Any]:
    """the strongest peaks of the outlier-aware periodogram of each activity
    indicator, per instrument that has it, up to pmax (the baseline of the
    indicator itself when None: the published velocities have none)"""
    if names == 'auto':
        names = [name for name in INDICATORS if name in data.indicators]
        names += [name for name in data.indicators
                  if name.upper().startswith('DTEMP') and name not in names]
    out = {}
    for name in names:
        try:
            series = data.indicator(name)
        except KeyError:
            continue
        # an indicator that does not move beyond its errors says nothing
        #   (the DTEMP of a temperature far from the star's)
        scale = robust_std(series.rv)
        if scale < 0.1 * np.median(series.err):
            log(f'{name}: scatter {scale:.2g} for errors of '
                f'{np.median(series.err):.2g}, no information: left out',
                'warn')
            continue
        # a value hundreds of sigma off (a failed exposure) takes the
        #   outlier model itself away: out before the periodogram
        dev = np.abs(series.rv - np.median(series.rv))
        gross = dev > GROSS * max(scale, np.median(series.err))
        if np.any(gross):
            log(f'{name}: {int(np.sum(gross))} value(s) beyond {GROSS:.0f} '
                f'robust sigma left out of its periodogram', 'warn')
            series = series.select(~gross)
        if series.n < 20:
            continue
        freq = frequency_grid(series.time, pmin, pmax or series.baseline, 10)
        dlnl = np.asarray(oap(series, freq, unit='point',
                              outliers=True)['dlnl'])
        peaks = find_peaks(freq, dlnl, npeaks)
        out[name] = dict(n=int(series.n), freq=freq, dlnl=dlnl,
                         peaks=[dict(period=float(1 / freq[it]),
                                     dlnl=float(dlnl[it])) for it in peaks])
    return out


def _activity_match(period: float, indicators: Dict[str, Any],
                    min_dlnl: float = 10.0) -> List[str]:
    """the indicators with a strong peak at the period, twice it or half of
    it (a signal of the velocities that activity may make)"""
    out = []
    for name, res in indicators.items():
        for pk in res['peaks']:
            if pk['dlnl'] < min_dlnl:
                continue
            for factor, word in ((1.0, ''), (2.0, ' (twice)'),
                                 (0.5, ' (half)')):
                if abs(pk['period'] * factor / period - 1) < ACTIVITY_MATCH:
                    out.append(f'{name} at {pk["period"]:.2f} d{word}, '
                               f'Delta lnL {pk["dlnl"]:.0f}')
    return out


def _gp_summary(gpres: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """what the GP found, as numbers and text (no fit, no array)"""
    if gpres is None:
        return None
    return dict(label=gpres['label'], runtime=gpres['runtime'],
                gp_null=gpres['gp_null'], gp_all=gpres['gp_all'],
                orbits=gpres['orbits'],
                whitened={key: dict(peaks=res['peaks'], kept=res['kept'],
                                    n=res['n'], nsim=len(res['maxsim']))
                          for key, res in gpres['whitened'].items()})


def _gp_text(gpsum: Optional[Dict[str, Any]]) -> str:
    """what the GP found, in words"""
    if gpsum is None:
        return ''
    lines = ['', f'With a GP of the activity: {gpsum["label"]}',
             '  GP alone: ' + ', '.join(f'{key} {val:.3g}' for key, val
                                        in gpsum['gp_null'].items())]
    for orb in gpsum['orbits']:
        lines.append(f'  {orb["P"][0]:.4f} d: K = {orb["K"][0]:.2f} '
                     f'+{orb["K"][2]:.2f}/-{orb["K"][1]:.2f} m/s, adds '
                     f'Delta ln L {orb["dlnl"]:+.1f}; whitened by the GP, '
                     f'Delta chi2 {orb["whitened"]["dchi2"]:.1f} (FAP '
                     f'{orb["whitened"]["fap"]:.1e})')
    for key, res in gpsum['whitened'].items():
        lines.append(f'  whitened periodogram ({key}): ' + ', '.join(
            f'{pk["period"]:.3f} d (Delta chi2 {pk["dchi2"]:.1f}, FAP '
            f'{pk["fap"]:.1e})' for pk in res['peaks'][:3]))
    return '\n'.join(lines)


def _rotation_match(period: float, star: Dict[str, Any]) -> List[str]:
    """whether a signal sits at the rotation period of the star (from the
    archive), or at its half or third"""
    prot = star.get('rotation')
    if not prot:
        return []
    out = []
    for harm in ROTATION_HARMONICS:
        if abs(prot / harm / period - 1) < ACTIVITY_MATCH:
            word = 'P_rot' if harm == 1 else f'P_rot/{harm}'
            out.append(f'{word} = {prot / harm:.2f} d (rotation {prot} d, '
                       f'NASA Exoplanet Archive)')
    return out


def _ephemeris(orb: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """the fitted conjunction of a known planet against its published
    ephemeris carried to it (the errors of the published time and period
    added, their covariance ignored)"""
    conj = conjunction(orb['known']) if 'known' in orb else None
    if conj is None or 'tc' not in orb:
        return None
    tc, tc_err = orb['tc'][0], 0.5 * (orb['tc'][1] + orb['tc'][2])
    ncycle = int(np.round((tc - conj['tc']) / conj['P']))
    pred = conj['tc'] + ncycle * conj['P']
    pred_err = float(np.hypot(conj['tc_err'] or 0.0,
                              ncycle * (conj['P_err'] or 0.0)))
    sig = float(np.hypot(tc_err, pred_err))
    diff = tc - pred
    return dict(tc=tc, tc_err=tc_err, ncycle=ncycle, predicted=pred,
                predicted_err=pred_err, diff=diff,
                nsigma=float(diff / sig) if sig > 0 else np.nan,
                dphase=float(diff / conj['P']), source=conj['source'])


def detailed_analysis(source: Union[str, RVData, Sequence[Any],
                                    None] = None,
                      outdir: str = 'koloa_detailed',
                      name: Optional[str] = None, target: Optional[str] = None,
                      archive: bool = True, dace: bool = False,
                      dace_folder: Optional[str] = None,
                      carmenes: bool = False,
                      literature: Optional[Sequence[Union[str, RVData]]]
                      = None, vizier: bool = False,
                      periods: Optional[Sequence[float]] = None,
                      gp: bool = True, gp_workers: int = 4,
                      gp_nsim: int = 20000,
                      refresh: bool = False, kmax: int = 4,
                      nsweep: int = 1500, nburn: Optional[int] = None,
                      pmin: float = 1.1, pmax: Optional[float] = None,
                      mcmc: bool = False, nsteps: int = 4000,
                      duck: bool = True, keys: Union[str, Sequence] = 'auto',
                      style: str = 'paper', seed: int = 1,
                      latex: bool = True, tess: bool = True,
                      site: Optional[str] = None,
                      fip_gp: Any = 'banded',
                      nightly: Optional[bool] = None,
                      detection_map: Any = False,
                      map_ninj: int = 10,
                      exclude: Optional[Sequence[str]] = None,
                      rotation: Optional[float] = None,
                      trend: bool = True, curvature: bool = False,
                      instruments: Optional[Sequence[Optional[str]]] = None,
                      toi: Any = None
                      ) -> Dict[str, Any]:
    """
    Everything koloa can say about a star (see the module's docstring)

    :param source: str, RVData, a list of them, or None: an LBL .rdb, a
                   csv, or a series, or several (SPIRou and NIRPS, HARPS
                   again through LBL, HARPS-N...: each file its
                   instruments, named from its own columns, or as
                   instruments= says); None for the archives only (DACE,
                   CARMENES DR1, VizieR), from the star's name (target)
    :param toi: None, True or a list of str: the TESS Objects of Interest
                of the star (the toi table of the NASA Exoplanet Archive)
                fitted with the ephemerides of TESS, P and the time of a
                transit as gaussian priors (K free, the phase held by the
                transit): True for all of them but the false positives,
                or their numbers ('175.01'); each is fitted whether the
                FIP finds it or not, its fold at the phase of the transit
    :param instruments: list of str or None, the instrument of each file,
                        in order ('auto', or None, for the name the file
                        gives)
    :param outdir: str, where the report, the summary and the figures go
    :param name: str or None, the name of the series
    :param target: str or None, the SIMBAD name of the star, for SIMBAD, the
                   archive, DACE, CARMENES and TESS (from the OBJECT column
                   of the file when None; needed without a file)
    :param archive: bool, ask the NASA Exoplanet Archive
    :param dace: bool, add the velocities DACE has of the star (off by
                 default: the velocities of the files only)
    :param dace_folder: str or None, where the DACE and CARMENES files are
                        kept (outdir when None)
    :param carmenes: bool, add the velocities of CARMENES DR1 (Ribas et al.
                     2023, the GTO of 2016 to 2020, about 360 M dwarfs of
                     the north), corrected for the nightly zero points
                     (koloa.gather; off by default)
    :param literature: list or None, published velocities to add: files (a
                       VizieR .dat, time velocity error instrument, without
                       a header; or a csv or .rdb with named columns) or
                       series (see koloa.literature.read)
    :param vizier: bool, add the velocities published with the known
                   planets, from VizieR (the tables of the papers of their
                   solutions in the archive; off by default)
    :param periods: list of float or None, more periods to test [days] (a
                    candidate of a paper that the archive does not list);
                    the known planets are tested at their periods anyway
    :param gp: bool, fit the signals with a GP of the activity (its prior
               the rotation period of the archive), and show the series
               with it (koloa.gpcheck)
    :param gp_workers: int, the GP fits at a time
    :param gp_nsim: int, the simulations of each whitened periodogram
    :param refresh: bool, ask DACE, the archive and VizieR again
    :param kmax: int, the largest number of signals of the FIP
    :param nsweep: int, recorded sweeps per chain of the FIP
    :param nburn: int or None, burn-in sweeps (a fifth of nsweep, at least
                  300, when None)
    :param pmin: float, the shortest period [days]
    :param pmax: float or None, the longest (the baseline when None)
    :param mcmc: bool, sample the orbits by MCMC (the maximum a posteriori
                 and its Laplace errors otherwise)
    :param nsteps: int, steps per walker of the MCMC
    :param duck: bool, run the duck test on every signal
    :param keys: the keys of koloa.outliers.explain
    :param style: str, paper (PDF) or web (SVG; the LaTeX report then has
                  no figure)
    :param seed: int, the seed
    :param latex: bool, write the report as LaTeX (<star>_report.tex) and
                  compile it to PDF (with pdflatex, when there is one)
    :param tess: bool, fetch the TESS light curves of the star (koloa.tess)
                 and look for a photometric peak at each signal, its half,
                 third or double (a check of the duck test, and a figure)
    :param site: str or None, the observatory of the plans that lift an
                 alias (a key of koloa.aliases.SITES; from the instruments
                 when None)
    :param fip_gp: a GP of the activity inside both FIPs (koloa.gpbasis):
                   'banded' (the default: the FIP by period bands, a local
                   GP that cannot reach the periods of each band and is only
                   as flexible as the data ask for, koloa.bandfip), 'auto'
                   (one local GP, and a rotation GP beside it when the
                   archive knows the rotation period, its prior +- 10 %),
                   'sho' (an SHO at the rotation period and one at its
                   half, celerite's RotationTerm, over every period: no
                   bands; the default when a rotation is given), None or
                   False (none), or a spec of oafip(gp=)
    :param detection_map: 'fip', 'search' or False (the default): map
                          which planets the series could have found, with
                          the signals found taken out. 'fip' looks for each
                          injected planet with the rule that decides (the
                          FIP with the GP, the period or any alias; adaptive
                          rounds per period band, koloa.fipmap; hours on a
                          long series); 'search' is the quicker blind
                          periodogram search of koloa.completeness (no GP,
                          an alias counts as missed)
    :param map_ninj: int, the injections per cell of the map (and
                     gp_workers the processes)
    :param rotation: float or None, a rotation period that can be trusted
                     [days] (a published one): the star's rotation for
                     every check, and the GP of the FIP an SHO at it
                     (fip_gp='sho', instead of the bands)
    :param trend: bool, fit a trend in time with the planets, in the
                  likelihood (and in the FIP), one for every instrument:
                  the acceleration of the star, dv/dt, reported in m/s/yr
                  with its errors (the default)
    :param curvature: bool, fit its change too, d2v/dt2 [m/s/yr^2] (a
                      second-order trend; implies the trend)
    :param exclude: list of str or None, instruments left out of the
                    analysis once the series is assembled (the file, DACE,
                    CARMENES, VizieR), by name (any case): NIRPS, HARPS03...
    :param nightly: bool, analyse the nightly means of the series once it is
                    assembled (RVData.nightly, koloa's default: an outlier
                    is then a night); False keeps the exposures; None is
                    koloa.data.NIGHTLY

    Planet or no planet is decided on the FIP of each period OR any of its
    aliases; which alias it is, is reported apart (the velocities folded at
    each, their share of the probability, and the nights to observe to
    lift an ambiguous one).

    :return: dict, every result (and outdir/<star>_report.txt,
             <star>_summary.json, <star>_report.tex and .pdf, and the
             figures)
    """
    start = time.time()
    os.makedirs(outdir, exist_ok=True)
    kplot.set_style(style)
    out: Dict[str, Any] = dict(figures=[])
    # the polynomial in time fitted with the planets: the acceleration of the
    #   star, and its change
    degree = 2 if curvature else (1 if trend else 0)
    step('the series')
    # 1. the series: a file, or none (the archives only, from the name)
    if source is None:
        if not target:
            raise ValueError('detailed_analysis: give a file of velocities, '
                             'or the SIMBAD name of the star (target=), or '
                             'both')
        if not (dace or carmenes or vizier or literature):
            raise ValueError('detailed_analysis: no file, and no archive '
                             'asked for: the velocities of an archive are '
                             'used only when asked (dace=, carmenes=, '
                             'vizier=; --dace, --carmenes, --vizier)')
        data, star = None, target
        log(f'koloa, detailed: {star}, no file: the velocities of the '
            f'archives')
    else:
        files = (list(source) if isinstance(source, (list, tuple))
                 else [source])
        parts = read_files(files, instruments, name)
        data = parts[0] if len(parts) == 1 else merge(parts,
                                                      name=parts[0].name)
        star = _target(data, target, files)
        log(f'koloa, detailed: {star} ({data.name}), {data.n} exposures in '
            f'{data.nseq} visits over {data.baseline:.0f} d, '
            f'{", ".join(data.instruments)}')
    step('SIMBAD and the archive')
    # 2. who the star is, and what is known of it
    ident, known = None, dict(host=None, planets=[])
    if archive or dace or carmenes or data is None:
        try:
            ident = resolve(star)
            log(f'SIMBAD: {ident["main"]} ('
                + ', '.join(val for val in (ident['gj'], ident['hd'],
                                            ident['hip'], ident['tic'])
                            if val) + ')', 'value')
        except (OSError, ValueError) as err:
            log(f'SIMBAD could not resolve {star}: {err}', 'warn')
    if archive and ident is not None:
        log(f'asking the NASA Exoplanet Archive for the planets of {star}')
        try:
            known = known_planets(
                name=star, path=os.path.join(outdir, 'archive.json'),
                refresh=refresh)
            for pl in known['planets']:
                last = (pl['solutions'] or [pl])[-1]
                log(f'  known: {pl["name"]}, P = {pl["P"]:.4f} d, K = '
                    f'{last.get("K")} m/s ({last.get("reference")})', 'value')
        except OSError as err:
            log(f'NASA Exoplanet Archive: {err}', 'warn')
    tois = _tois(toi, ident)
    out['tois'] = tois
    if rotation:
        known['star'] = dict(known.get('star') or {},
                             rotation=float(rotation),
                             rotation_source='given')
        if fip_gp == 'banded':
            fip_gp = 'sho'
        log(f'rotation period given: {float(rotation):g} d'
            + (', the GP of the FIP an SHO at it and at its half (no bands)'
               if fip_gp == 'sho' else ''), 'value')
    out['ident'], out['known'] = ident, known
    step('more velocities: DACE, CARMENES, VizieR')
    # 3. more velocities: DACE, CARMENES, the ones given, the ones on VizieR
    sources = [] if data is None else [
        dict(kind='file', label=_label(src), n=int(part.n),
             instruments=_counts(part), note='')
        for src, part in zip(files, parts)]
    if not dace:
        sources.append(dict(kind='DACE', label='DACE', n=0, instruments={},
                            note='not asked'))
    elif ident is None:
        sources.append(dict(kind='DACE', label='DACE', n=0, instruments={},
                            note='not asked: SIMBAD did not resolve the star'))
    else:
        names = dace_names(ident, star)
        more = fetch_dace(names, dace_folder or outdir,
                          exclude=data.instruments if data else (),
                          times=data.time if data else None,
                          refresh=refresh)
        if more is not None:
            data = (more if data is None
                    else merge([data, more], name=data.name))
            log(f'with DACE: {data.n} exposures in {data.nseq} visits, '
                + ', '.join(f'{inst} {np.sum(data.inst == inst)}'
                            for inst in data.instruments), 'value')
        sources.append(dict(
            kind='DACE', label='DACE', n=int(more.n) if more else 0,
            instruments=_counts(more) if more else {},
            note=(('what the file does not have; an instrument of the file '
                   'too is named <inst>_DACE (another pipeline, its own '
                   'offset)') if more else
                  'nothing added (unreachable, or no velocity the file does '
                  'not have) under ' + ', '.join(names))))
    if not carmenes:
        sources.append(dict(kind='CARMENES', label='CARMENES DR1', n=0,
                            instruments={}, note='not asked'))
    elif ident is None or ident.get('ra') is None:
        sources.append(dict(kind='CARMENES', label='CARMENES DR1', n=0,
                            instruments={}, note='not asked: SIMBAD did not '
                            'resolve the star'))
    else:
        more, note = _carmenes(ident, star, dace_folder or outdir, data,
                               refresh)
        if more is not None:
            data = (more if data is None
                    else merge([data, more], name=data.name))
            log(f'with CARMENES DR1: {data.n} exposures, ' + ', '.join(
                f'{inst} {np.sum(data.inst == inst)}'
                for inst in data.instruments), 'value')
        sources.append(dict(kind='CARMENES', label='CARMENES DR1',
                            n=int(more.n) if more else 0,
                            instruments=_counts(more) if more else {},
                            note=note))
    others, labels, kinds = [], [], []
    for item in literature or []:
        try:
            series = klit.read(item)
        except (OSError, ValueError) as err:
            log(f'literature, {item}: {err}', 'warn')
            continue
        others.append(series)
        labels.append(os.path.basename(item) if isinstance(item, str)
                      else series.name)
        kinds.append(('given', '', ''))
    if vizier and ident is not None:
        # the surveys and the star's papers on VizieR (koloa.published),
        #   then the papers of its known planets
        from koloa import published as kpub
        try:
            pnotes = kpub.fetch(ident, os.path.join(dace_folder or outdir,
                                                    'published'),
                                refresh=refresh)
        except (OSError, ValueError) as err:
            pnotes = []
            log(f'published velocities: {err}', 'warn')
        for note in pnotes:
            if note.get('file'):
                others.append(RVData.from_csv(os.path.join(
                    dace_folder or outdir, 'published', note['file']),
                    inst='inst', name=note['reference']))
                labels.append(note['reference'])
                kinds.append(('VizieR', note['catalogue'], note['note']))
            elif note['kind'] == 'survey':
                sources.append(dict(kind='VizieR', label=note['reference'],
                                    n=0, instruments={},
                                    note=f'{note["catalogue"]}: '
                                         f'{note["note"]}'))
    if vizier and known.get('planets') and ident is not None:
        pubs, notes = klit.fetch(known, [ident['main']] + ident['aliases'],
                                 os.path.join(outdir, 'literature'),
                                 refresh=refresh)
        found = iter(pubs)
        for note in notes:
            if note.get('file'):
                others.append(next(found))
                labels.append(note['reference'])
                kinds.append(('VizieR', note['catalogue'], note['note']))
            else:
                sources.append(dict(kind='VizieR', label=note['reference'],
                                    n=0, instruments={},
                                    note=f'{note["bibcode"]}: {note["note"]}'))
    if others and data is None:
        # no file and nothing on DACE or CARMENES: the first published
        #   series is the base
        data = others.pop(0)
        sources.append(dict(kind=kinds[0][0], label=labels.pop(0),
                            n=int(data.n), instruments=_counts(data),
                            note=kinds.pop(0)[2]))
        data.name = name or star
    if data is None:
        raise ValueError(f'no velocities of {star}: no file, and nothing on '
                         f'DACE, CARMENES DR1 or VizieR')
    if source is None:
        data.name = name or star
    if others:
        data, added = klit.add(data, others, labels)
        for series, info, (kind, cat, note) in zip(others, added, kinds):
            words = [cat] if cat else []
            if info['same']:
                words.append(f'{info["same"]} left out: the same spectra as '
                             f'the file')
            if info['renamed']:
                words.append('renamed ' + ', '.join(
                    f'{old} to {new}' for old, new in info['renamed'].items()))
            sources.append(dict(kind=kind, label=info['label'], n=info['n'],
                                instruments=_counts(series), note='; '.join(
                                    words) or note))
        log(f'with the literature: {data.n} exposures in {data.nseq} visits '
            f'over {data.baseline:.0f} d, ' + ', '.join(
                f'{inst} {np.sum(data.inst == inst)}'
                for inst in data.instruments), 'value')
    if exclude:
        wanted = {str(name).strip().upper() for name in exclude}
        drop = [inst for inst in data.instruments if inst.upper() in wanted]
        missing = wanted - {inst.upper() for inst in data.instruments}
        if missing:
            log(f'not in the series, nothing to leave out: '
                f'{", ".join(sorted(missing))}', 'warn')
        if drop:
            keep = ~np.isin(data.inst, drop)
            if not np.any(keep):
                raise ValueError('every instrument is left out (exclude=)')
            counts = {inst: int(np.sum(data.inst == inst)) for inst in drop}
            data = data.select(keep)
            sources.append(dict(kind='left out', label=', '.join(drop),
                                n=-sum(counts.values()), instruments=counts,
                                note='left out of the analysis, as asked'))
            log('left out, as asked: ' + ', '.join(
                f'{inst} ({num})' for inst, num in counts.items())
                + f'; {data.n} exposures remain', 'value')
    out['sources'] = sources
    nexp = data.n
    if nightly is None:
        from koloa import data as kdata
        nightly = kdata.NIGHTLY
    if nightly:
        means = data.nightly()
        if means is not data:
            log(f'the nightly means: {data.n} exposures in {means.n} nights '
                f'(nightly=False keeps the exposures)', 'value')
            data = means
    out['data'] = data
    multi = len(data.instruments) > 1
    # a jitter per visit, when visits hold several points
    seq_jitter = 'instrument' if multi and data.nseq < data.n else None
    nburn = nburn if nburn is not None else max(300, nsweep // 5)
    step('the FIP, first pass')
    # 4. the noise without planets, and the FIP
    with blas_threads(1):
        noise = RVModel(data, [], likelihood='mixture', unit='both',
                        trend=degree, seq_jitter=seq_jitter).fit(nstart=2,
                                                                quiet=True)
    gpspec = _gp_spec(fip_gp, known)
    fip1, info1 = _fip(data, noise, kmax, nsweep, nburn, seed,
                       'FIP, noise without planets', gp=gpspec, trend=degree)
    # planet or no planet: the period or any of its aliases; one period
    #   per family of aliases
    found = sorted(_found(fip1))
    step('the signals, fitted together')
    # 5. the signals, the known planets and the periods given, fitted
    #   together (a known planet the FIP does not find is tested all the
    #   same), against the known planets
    # the TOIs first: their ephemeris wins over a period found near theirs
    tests = [(item['P'], f'TOI {item["toi"]}') for item in tois]
    tests += [(per, 'FIP') for per in found]
    tests += [(float(pl['P']), f'known ({pl["name"]})')
              for pl in known.get('planets', [])
              if pl.get('P') and pmin <= pl['P'] < data.baseline]
    tests += [(float(per), 'given') for per in periods or []]
    chosen: List[Any] = []
    for per, origin in tests:
        if all(abs(per / other - 1) > 0.02 for other, _ in chosen):
            chosen.append((per, origin))
    chosen.sort()
    mstar = (known.get('star') or {}).get('mass')
    fit, orbits, fip2, info2 = noise, [], fip1, info1
    if chosen:
        bytoi = {f'TOI {item["toi"]}': item for item in tois}
        planets = [_planet(per, bytoi.get(origin)) for per, origin in chosen]
        with blas_threads(1):
            if mcmc:
                fit = mcmc_orbits(data, planets, likelihood='mixture',
                                  unit='both', trend=degree,
                                  seq_jitter=seq_jitter,
                                  nsteps=nsteps, nburn=nsteps // 3, nstart=2,
                                  seed=seed, quiet=True)
            else:
                fit = RVModel(data, planets, likelihood='mixture',
                              unit='both', trend=degree,
                              seq_jitter=seq_jitter).fit(nstart=2,
                                                         quiet=True)
        for (per, origin), orb in zip(chosen, fit.orbits(mstar=mstar)):
            orbits.append(compare({key: [float(val) for val in orb[key]]
                                   for key in ('P', 'K', 'e', 'msini', 'tc')},
                                  known, MATCH))
            orbits[-1]['origin'] = origin
            orbits[-1]['tested'] = float(per)
            eph = _ephemeris(orbits[-1])
            if eph is not None:
                orbits[-1]['ephemeris'] = eph
                log(f'{orbits[-1]["known"]["name"]}: conjunction '
                    f'{eph["diff"]:+.3f} d ({eph["nsigma"]:+.1f} sigma, '
                    f'{eph["dphase"]:+.3f} in phase) from the published '
                    f'ephemeris carried {eph["ncycle"]} cycles', 'value')
        # 6. the FIP again, with the noise of the fit with the planets
        step('the FIP, second pass')
        fip2, info2 = _fip(data, fit, kmax, nsweep, nburn, seed + 1,
                           'FIP, noise with the planets', gp=gpspec,
                           decided=info1.get('banded'), trend=degree)
    else:
        log('no interval with FIP < 1 %, no known planet and no period '
            'given: nothing to fit', 'warn')
    for orb in orbits:
        # the FIP of the best interval that holds each signal
        orb['fip'] = fip2.fip_containing(orb['P'][0], 1 / data.baseline)
        orb['family_fip'] = fip2.family_containing(orb['P'][0],
                                                   1 / data.baseline)
    out.update(fit=fit, orbits=orbits, fip_first=fip1, fip_second=fip2)
    # the acceleration of the star: the trend fitted with the planets
    out['acceleration'] = accel = _acceleration(fit, degree, mcmc and orbits)
    for line in _accel_lines(accel):
        log(line, 'value')
    step('the activity indicators')
    # 7. the activity indicators, and the rotation of the star
    indic = _indicators(data, pmin, pmax)
    out['indicators'] = indic
    for orb in orbits:
        orb['activity'] = _activity_match(orb['P'][0], indic)
        orb['rotation'] = _rotation_match(orb['P'][0], known.get('star') or {})
    # 7b. the signals against a GP of the activity
    gpres = None
    if gp:
        step('the GP of the activity')
        try:
            gpres = gp_signals(
                data, [orb['P'][0] for orb in orbits],
                star=known.get('star'), reference=next(
                    iter(sources[0]['instruments']), None),
                seq_jitter=seq_jitter, mstar=mstar, workers=gp_workers,
                nsim=gp_nsim, pmin=pmin, pmax=pmax, seed=seed,
                ranges=[(0.98 * orb['tested'], 1.02 * orb['tested'])
                        for orb in orbits])
            for orb, gorb in zip(orbits, gpres['orbits']):
                orb['gp'] = gorb
        except Exception as err:  # the GP is a help, not a stop
            log(f'GP: {err}', 'warn')
    out['gp'] = gpres
    # 8. the duck test, with the TESS light curves of the star (fetched
    #    once, for every signal)
    lcs, phot = None, None
    # TESS needs the network, as the archive does
    if tess and archive:
        step('TESS')
        from koloa import tess as ktess
        try:
            lcs = ktess.light_curves((ident or {}).get('tic') or star)
            phot = ktess.periodicity(lcs, [orb['P'][0] for orb in orbits])
        except Exception as err:  # no network, or an unknown star
            log(f'TESS: {err}', 'warn')
    out['tess'] = phot
    ducks = {}
    if duck:
        step('the duck tests')
        for orb in orbits:
            try:
                report = duck_test(data, orb['P'][0], fipres=fip2, gp=False,
                                   unit='both', quiet=True,
                                   tess=lcs if lcs is not None else False,
                                   aliases=True, target=star,
                                   archive=archive, site=site)
                ducks[f'{orb["P"][0]:.4f}'] = report
            except Exception as err:  # the test is a help, not a stop
                log(f'duck test at {orb["P"][0]:.4f} d: {err}', 'warn')
    out['duck'] = ducks
    # 8c. the detection map: which planets the series could have found
    dmap = None
    if detection_map:
        step('the detection map')
        prot = (known.get('star') or {}).get('rotation')
        try:
            if detection_map == 'search':
                dmap = _detection_map(data, fit, len(orbits), pmin, map_ninj,
                                      gp_workers, seed, mstar, prot)
            else:
                dmap = _fip_detection_map(
                    data, fit, len(orbits), pmin,
                    fip2.settings['bands'] if gpspec == 'banded' else gpspec,
                    gp_workers, seed, mstar, prot)
        except Exception as err:  # the map is a help, not a stop
            log(f'detection map: {err}', 'warn')
    out['detection_map'] = dmap
    step('why each outlier is one')
    # 9. why each outlier is one
    why = explain(fit, keys=keys)
    out['outliers'] = why
    step('the figures')
    # 10. the figures, the report and the summary
    safe = data.name.replace(' ', '_')
    figs, named = out['figures'], {}

    def keep(fig, role):
        """a figure saved as <star>_<role>, and known by its role"""
        named[role] = kplot.savefig(fig, os.path.join(outdir,
                                                      f'{safe}_{role}'))
        figs.append(named[role])

    marks = [pl['P'] for pl in known.get('planets', []) if pl.get('P')]
    keep(kplot.timeseries(
        data, fit.outlier_prob,
        title=f'{star}: every instrument, each exposure coloured by its '
              f'outlier probability'), 'rv')
    if orbits or accel:
        keep(kplot.model_series(
            fit, accel, title=f'{star}: the velocities and the best model '
                              f'(the Keplerians and the trend, fitted '
                              f'together)'), 'model')
    if gpres is not None:
        # the whole series with the GP, the periodograms whitened by it,
        #   and every signal folded without it
        best = gpres['fits'].get('all', gpres['fits']['null'])
        keep(figure_sequence(best, f'{star}: the whole series, the GP of the '
                                   f'activity (one sigma shaded) and the GP '
                                   f'with the signals'), 'gp_sequence')
        prot = (known.get('star') or {}).get('rotation')
        wmarks = [(orb['P'][0], f'{orb["P"][0]:.2f} d', 'outlier')
                  for orb in orbits]
        if prot:
            wmarks += [(prot, 'P_rot', 'koloa'), (prot / 2, 'P_rot/2', 'koloa')]
        keep(figure_whitened(gpres['whitened'], wmarks,
                             f'{star}: periodograms whitened by the GP'),
             'gp_whitened')
        for ip, orb in enumerate(orbits):
            published = None
            if orb.get('known', {}).get('K') is not None:
                published = dict(K=orb['known']['K'])
                if orb.get('ephemeris'):
                    published['dphase'] = -orb['ephemeris']['dphase']
            keep(figure_fold(best, ip, f'{star}, {orb["P"][0]:.3f} d, the GP '
                                       f'and the other signals removed',
                             published), f'gp_fold_{ip}')
    for label, res in (('first', fip1), ('second', fip2)):
        keep(kplot.fip_family(
            res, marks=marks, threshold=THRESHOLD,
            title=f'{star}: FIP ({label} pass; known planets dashed)'),
            f'fip_{label}')
    for ip, orb in enumerate(orbits):
        keep(kplot.phase(
            fit, planet=ip, level='point',
            title=f'{star}, {orb["P"][0]:.3f} d: K = {orb["K"][0]:.2f} m/s'
                  + (f' ({orb["origin"]}; phase 0, the transit of TESS)'
                     if str(orb.get('origin', '')).startswith('TOI') else '')),
            f'phase_{ip}')
        report = ducks.get(f'{orb["P"][0]:.4f}')
        if report is None:
            continue
        # the duck test, in figures: which visits hold the peak up, and
        #   whether the signal keeps its amplitude and phase
        try:
            jack = report.details.get('jackknife')
            jfreq = report.details.get('jackknife_freq')
            if jack is not None and jfreq is not None:
                keep(kplot.jackknife(jfreq, jack, orb['P'][0],
                                     report.details.get('jackknife_data',
                                                        data),
                                     unit='sequence',
                                     title=f'{star}: leave one visit out'),
                     f'duck_{ip}_jackknife')
            for split, coh in report.details.get('coherence', {}).items():
                if len(coh['chunks']) >= 2:
                    keep(kplot.coherence(coh, title=f'{star}, '
                                         f'{orb["P"][0]:.3f} d: coherent? '
                                         f'(by {split})'),
                         f'duck_{ip}_coherence_{split}')
            # the velocities folded at each alias, and the plan
            from koloa import aliases as kal
            if report.details.get('alias_solutions'):
                keep(kal.figure(report.details['alias_solutions'],
                                title=f'{star}, {orb["P"][0]:.3f} d: the '
                                      f'period and its aliases'),
                     f'duck_{ip}_aliases')
            pln = report.details.get('alias_plan')
            if pln and pln.get('nights'):
                keep(kal.plan_figure(pln, report.details['alias_solutions'],
                                     title=f'{star}: when to observe to '
                                           f'lift the alias'),
                     f'duck_{ip}_plan')
        except Exception as err:  # a figure is a help, not a stop
            log(f'duck test figures at {orb["P"][0]:.4f} d: {err}', 'warn')
    if dmap is not None and dmap.get('kind') == 'fip':
        from koloa.fipmap import figure as fipmap_figure
        marks = [dict(period=orb['P'][0], K=orb['K'][0], kind='signal')
                 for orb in orbits]
        marks += [dict(period=pl['P'], K=pl['K'], kind='known')
                  for pl in known.get('planets', []) if pl.get('P')
                  and pl.get('K')]
        keep(fipmap_figure(dmap, marks, dmap.get('prot'),
                           title=f'{star}: which planets the series could '
                                 f'have found (the FIP with the GP)'),
             'detection_map')
    elif dmap is not None:
        fig = kplot.recovery_map(dmap['map'], title=f'{star}: which planets '
                                 f'the series could have found (blind '
                                 f'search)')
        axm = fig.axes[0]
        for orb in orbits:
            axm.plot(orb['P'][0], orb['K'][0], marker='*', ms=11,
                     color=kplot.C['outlier'], mec='white', ls='none')
        for pl in known.get('planets', []):
            if pl.get('P') and pl.get('K'):
                axm.plot(pl['P'], pl['K'], marker='o', ms=7, mfc='none',
                         mec=kplot.C['text'], ls='none')
        # the rotation and its half, where the activity is found
        if dmap.get('prot'):
            for per in (dmap['prot'], dmap['prot'] / 2):
                axm.axvline(per, color=kplot.C['muted'], lw=1.0,
                            ls=(0, (1, 2)))
        keep(fig, 'detection_map')
    if phot is not None and phot.get('sectors'):
        keep(ktess.figure(phot, [orb['P'][0] for orb in orbits],
                          title=f'{star}: TESS, the periodogram of each '
                                f'sector (the signals and their harmonics '
                                f'in red)'), 'tess')
    if indic:
        keep(_indicator_figure(indic, orbits, star,
                               (known.get('star') or {}).get('rotation')),
             'indicators')
    if why.units:
        keep(kplot.outlier_keys(why), 'outliers')
    gpsum = _gp_summary(gpres)
    text = _report(star, data, ident, known, fip1, fip2, orbits, indic,
                   ducks, why, time.time() - start, sources) + _gp_text(gpsum)
    if tois:
        text += ('\n\nTESS Objects of Interest, fitted with the ephemerides '
                 'of TESS (P and the time of a transit as gaussian priors, K '
                 'free; phase 0 of their folds is the transit)')
        for item in tois:
            orb = next((orb for orb in orbits if orb.get('origin')
                        == f'TOI {item["toi"]}'), None)
            text += (f'\n  TOI-{item["toi"]} ({item["disposition"]}): P = '
                     f'{item["P"]:.7f} d, transit at {item["tc"]:.5f}'
                     + (f'; K = {orb["K"][0]:.2f} -{orb["K"][1]:.2f} '
                        f'+{orb["K"][2]:.2f} m/s' if orb else
                        '; not fitted'))
    if accel:
        text += ('\n\nThe trend fitted with the planets (one for every '
                 'instrument)\n' + '\n'.join(f'  {line}' for line in
                                             _accel_lines(accel)))
    if phot is not None:
        text += ('\n\nTESS photometry (TIC ' + str(phot.get('tic')) + ')\n  '
                 + phot['summary'] + '\n'
                 + ''.join(f'  {chk["period"]:.4f} d: [{chk["status"]}] '
                           f'{chk["summary"]}\n' for chk in phot['checks'])
                 + '\n'.join('  ' + line for line in ktess.table(phot)))
    if dmap is not None:
        text += ('\n\nDetection map (which planets the series could have '
                 'found, the signals found taken out)\n'
                 + '\n'.join(_map_lines(dmap)))
    if data.n < nexp:
        text += (f'\n\nThe analysis ran on the nightly means: {nexp} '
                 f'exposures in {data.n} nights (nightly=False keeps the '
                 f'exposures).')
    with open(os.path.join(outdir, f'{safe}_report.txt'), 'w') as handle:
        handle.write(text + '\n')
    out['report'] = text
    # 11. the report, in LaTeX and PDF
    paths = {}
    if latex:
        step('the report, LaTeX and PDF')
        from koloa.latex import compile_pdf, detailed_report
        tex = os.path.join(outdir, f'{safe}_report.tex')
        rep = dict(
            star=star, source=(', '.join(os.path.abspath(src) for src in files
                                     if isinstance(src, str)) or None)
            if source is not None else None,
            data=data, ident=ident, known=known, sources=sources,
            fip_first=fip1, fip_second=fip2,
            inflation_first=info1.get('inflation', {}),
            inflation_second=info2.get('inflation', {}), orbits=orbits,
            indicators=indic, duck=ducks, outliers=why, figures=named,
            gp=gpsum, threshold=THRESHOLD, runtime=time.time() - start,
            tess=phot, detection_map=dmap, acceleration=accel, tois=tois,
            settings=dict(archive=archive, dace=dace, carmenes=carmenes,
                          toi=', '.join(f'TOI-{item["toi"]}' for item in tois)
                          or 'none',
                          trend=degree,
                          rotation=(f'{float(rotation):g} d, given' if rotation
                                    else 'the archive\'s'),
                          exclude=', '.join(exclude or []) or 'none',
                          vizier=vizier,
                          literature=len(literature or []),
                          periods=', '.join(f'{per}' for per in periods or [])
                          or '--', gp=gp, gp_nsim=gp_nsim, kmax=kmax,
                          nsweep=nsweep, nburn=nburn, pmin=pmin,
                          pmax=pmax or 'the baseline', mcmc=mcmc,
                          nsteps=nsteps if mcmc else '--', seed=seed,
                          style=style, tess=tess,
                          fip_gp=_fip_gp_text(fip2),
                          nightly=(f'{nexp} exposures in {data.n} nights'
                                   if data.n < nexp else 'no')),
            files=sorted({name for name in os.listdir(outdir)
                          if not name.startswith('.')} | {
                f'{safe}_summary.json', f'{safe}_report.tex',
                f'{safe}_report.pdf'}))
        try:
            paths['tex'] = detailed_report(rep, tex)
            paths['pdf'] = compile_pdf(tex)
        except Exception as err:  # the report is a help, not a stop
            log(f'LaTeX report: {err}', 'warn')
    out['latex'] = paths
    summary = _summary(star, data, ident, known, fip1, fip2, orbits, indic,
                       ducks, why, figs)
    summary.update(sources=sources, gp=gpsum, report_tex=paths.get('tex'),
                   report_pdf=paths.get('pdf'), acceleration=accel)
    if phot is not None:
        summary['tess'] = ktess.light(phot)
    if dmap is not None and dmap.get('kind') == 'fip':
        summary['detection_map'] = {key: val for key, val in dmap.items()}
    elif dmap is not None:
        rmap = dmap['map']
        summary['detection_map'] = dict(
            period_edges=rmap.period_edges, amp_edges=rmap.amp_edges,
            rate=rmap.rate, threshold=rmap.threshold,
            false_alarm=rmap.false_alarm, K50=dmap['K50'], K90=dmap['K90'],
            msini50=dmap['msini50'], msini90=dmap['msini90'],
            unreliable=dmap['unreliable'])
    with open(os.path.join(outdir, f'{safe}_summary.json'), 'w') as handle:
        json.dump(summary, handle, indent=1, default=_json)
    log(f'done in {(time.time() - start) / 60:.1f} min: the report'
        + (' (PDF)' if paths.get('pdf') else '') + f', the summary and '
        f'{len(figs)} figures in {outdir}')
    return out


def _json(obj):
    """numbers and arrays for json.dump"""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return str(obj)


def _counts(data: Optional[RVData]) -> Dict[str, int]:
    """the exposures of each instrument of a series"""
    if data is None:
        return {}
    return {str(inst): int(np.sum(data.inst == inst))
            for inst in data.instruments}


def _indicator_figure(indic: Dict[str, Any], orbits, star: str,
                      prot: Optional[float] = None):
    """the outlier-aware periodogram of each indicator, the signals of the
    velocities marked (dashed), and the rotation period of the star and
    its half (dotted)"""
    fig, axes = kplot.plt.subplots(len(indic), 1, sharex=True,
                                   figsize=(7.0, 1.3 * len(indic) + 0.8))
    axes = np.atleast_1d(axes)
    for ax, (name, res) in zip(axes, indic.items()):
        ax.plot(1 / res['freq'], res['dlnl'], color=kplot.C['koloa'], lw=0.8)
        for orb in orbits:
            ax.axvline(orb['P'][0], color=kplot.C['muted'], ls='--', lw=0.8)
        for harm in ((1, 2) if prot else ()):
            ax.axvline(prot / harm, color=kplot.C['gaussian'], ls=':',
                       lw=1.0)
        ax.set_ylabel(r'$\Delta\ln L$', fontsize=8)
        ax.text(0.01, 0.85, name, transform=ax.transAxes, fontsize=8,
                color=kplot.C['text'])
    axes[-1].set_xscale('log')
    axes[-1].set_xlabel('period [d]')
    axes[0].set_title(f'{star}: activity indicators (the signals of the '
                      f'velocities dashed' + (f'; the rotation, {prot} d, and '
                                              f'its half dotted' if prot
                                              else '') + ')',
                      loc='left', fontsize=9)
    fig.tight_layout()
    return fig


def _report(star, data, ident, known, fip1, fip2, orbits, indic, ducks, why,
            runtime, sources=()) -> str:
    """the report, in words"""
    lines = [f'koloa, detailed analysis of {star}', '=' * 70, '']
    lines.append(f'{data.n} exposures in {data.nseq} visits over '
                 f'{data.baseline:.0f} d: ' + ', '.join(
                     f'{inst} {np.sum(data.inst == inst)}'
                     for inst in data.instruments))
    for src in sources:
        lines.append(f'  {src["kind"]}, {src["label"]}: {src["n"]} velocities'
                     + (f' ({src["note"]})' if src.get('note') else ''))
    if ident:
        lines.append(f'SIMBAD: {ident["main"]}; ' + ', '.join(
            val for val in (ident['gj'], ident['hd'], ident['hip'],
                            ident['tic'], ident['gaia_dr3']) if val))
    params = {key: val for key, val in (known.get('star') or {}).items()
              if val is not None}
    if params:
        lines.append('The star (NASA Exoplanet Archive): ' + ', '.join(
            f'{key} {val}' for key, val in params.items()))
    lines += ['', 'Known planets (NASA Exoplanet Archive)']
    if not known.get('planets'):
        lines.append('  none')
    for pl in known.get('planets', []):
        sols = '; '.join(f'K {sol["K"]} +- {sol["K_err"]} ({sol["reference"]})'
                         for sol in reversed(pl.get('solutions') or [])
                         if sol.get('K_err'))
        lines.append(f'  {pl["name"]}: P = {pl["P"]} d; {sols}')
    for label, res in (('FIP, noise without planets', fip1),
                       ('FIP, noise with the planets', fip2)):
        lines += ['', label, '  P(k) = ' + ', '.join(f'{val:.2f}'
                                                     for val in res.pk),
                  '  ' + _fip_gp_text(res)]
        for pk in res.peaks:
            flag = 'DETECTED' if _decisive(pk) < THRESHOLD else ''
            fam = pk.get('family_fip')
            lines.append(f'  {pk["period"]:12.4f} d  FIP {pk["fip"]:.2e}'
                         + (f'  FIP (the period or any alias) {fam:.2e}'
                            if fam is not None else '') + f'  {flag}')
    lines += ['', 'Signals (planet or no planet: the FIP of the period or '
                  'any of its aliases)']
    if not orbits:
        lines.append('  none with FIP < 1 %')
    for orb in orbits:
        lines.append(f'  {orb["P"][0]:.4f} d: K = {orb["K"][0]:.2f} '
                     f'+{orb["K"][2]:.2f}/-{orb["K"][1]:.2f} m/s, e = '
                     f'{orb["e"][0]:.2f}, m sin i = {orb["msini"][0]:.2f} Me')
        if 'known' in orb:
            lines.append(f'    {orb["known"]["name"]}: ' + '; '.join(
                f'{sol["K"]} +- {sol["K_err"]:.2f} m/s ({sol["reference"]}, '
                f'{sol["z"]:+.1f} sigma)' for sol in orb['comparisons']))
        else:
            lines.append('    not in the archive')
        eph = orb.get('ephemeris')
        if eph:
            lines.append(f'    conjunction {eph["diff"]:+.3f} d '
                         f'({eph["nsigma"]:+.1f} sigma, {eph["dphase"]:+.3f} '
                         f'in phase) from the published ephemeris carried '
                         f'{eph["ncycle"]} cycles')
        if orb.get('activity'):
            lines.append('    ACTIVITY at this period: '
                         + '; '.join(orb['activity']))
        if orb.get('rotation'):
            lines.append('    ROTATION at this period: '
                         + '; '.join(orb['rotation']))
        if orb.get('family_fip') is not None:
            lines.append(f'    FIP (the period or any alias) '
                         f'{orb["family_fip"]:.2e}; the period alone '
                         f'{orb["fip"]:.2e}')
        report = ducks.get(f'{orb["P"][0]:.4f}')
        if report is not None:
            lines.append(f'    duck test: {report.verdict}')
            from koloa.aliases import describe, plan_text
            if report.details.get('alias_solutions'):
                lines += ['    the period and its aliases:'] + [
                    '      ' + line for line in
                    describe(report.details['alias_solutions'])]
            if report.details.get('alias_plan'):
                lines += ['    ' + line for line in
                          plan_text(report.details['alias_plan'])]
    lines += ['', 'Activity indicators (strongest peaks of their outlier-'
                  'aware periodograms)']
    for name, res in indic.items():
        lines.append(f'  {name}: ' + ', '.join(
            f'{pk["period"]:.2f} d (Delta lnL {pk["dlnl"]:.0f})'
            for pk in res['peaks']))
    lines += ['', why.text(), '', f'Run in {runtime / 60:.1f} min.']
    return '\n'.join(lines)


def _summary(star, data, ident, known, fip1, fip2, orbits, indic, ducks,
             why, figs) -> Dict[str, Any]:
    """everything, as numbers and text"""
    def fipsum(res):
        return dict(pk=[float(val) for val in res.pk],
                    gp=res.gp_summary(), gp_components=res.settings.get('gp'),
                    peaks=[dict(period=float(pk['period']),
                                fip=float(pk['fip'])) for pk in res.peaks])
    return dict(
        star=star, n=int(data.n), nvisits=int(data.nseq),
        instruments={inst: int(np.sum(data.inst == inst))
                     for inst in data.instruments},
        baseline=float(data.baseline), ident=ident, known=known,
        fip_first=fipsum(fip1), fip_second=fipsum(fip2), orbits=orbits,
        indicators={name: res['peaks'] for name, res in indic.items()},
        duck={key: report.verdict for key, report in ducks.items()},
        aliases={key: dict(
            family_fip=report.details.get('family_fip'),
            odds=report.details.get('alias_odds'),
            solutions=_light_aliases(report.details.get('alias_solutions')),
            plan=report.details.get('alias_plan'))
            for key, report in ducks.items()},
        outliers=why.to_dict(), figures=figs)


def _light_aliases(sols):
    """the alias solutions without their fits"""
    from koloa.aliases import light
    return light(sols) if sols else None


# =============================================================================
# End of code
# =============================================================================

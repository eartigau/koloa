#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The TESS light curves of a star, and whether its brightness varies at the
period of a velocity signal.

Spots that rotate with the star modulate its brightness and its velocities
at the rotation period and its harmonics (P/2, P/3); a planet does not
change the brightness, unless it transits. A photometric peak at the period
of a velocity signal, at its half or third, or at its double, therefore
speaks against a planet.

- light_curves(name): the light curve of every sector at MAST (by the TIC
  number, from Sesame), one per sector: SPOC's 2-minute light curve when
  there is one, else TESS-SPOC's or QLP's from the full frames. The FITS
  files are kept in a folder and read back.
- periodicity(lcs, periods): each sector binned to 30 minutes (flares and
  bad points clipped first), its GLS periodogram, its strongest periods and
  their amplitudes, and at each velocity period and its harmonics the power
  there against the sector's highest peak.
  A match is a local maximum of the periodogram within one peak width
  (1/T) of the tested period that holds half the sector's highest power;
  the check flags it when it is the strongest peak of two sectors or more
  (or of the only one), unless it sits at the orbit of TESS (13.7 d) or its
  half, where the light curves keep systematics.
- figure(res): the periodograms of the sectors with the periods marked, and
  the light curve of the sector that best shows the strongest period.

What TESS cannot tell: a sector lasts 27 days and SPOC's PDCSAP flux
flattens slower changes, so periods longer than about half a sector (13 d)
are not tested; putting sectors side by side does not help, since each is
normalised on its own. The rotation of most M dwarfs (tens of days) is out
of reach; its harmonics may not be.

Created on 2026-09-30

@author: artigau
"""
import json
import os
import time
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np

from koloa.log import log
from koloa.periodogram import gls
from koloa.paths import cache

# =============================================================================
# Define variables
# =============================================================================
MAST = 'https://mast.stsci.edu/api/v0/invoke'
#: which sectors' full frames hold a position (MAST's TESScut)
TESSCUT = 'https://mast.stsci.edu/tesscut/api/v0.1/sector'
DOWNLOAD = 'https://mast.stsci.edu/api/v0.1/Download/file?uri='
#: where the light curves are kept (they are public, and shared by runs)
CACHE = cache('tess')
#: the light curves of a sector, the best first: (provenance, exposure [s])
PREFERENCE = [('SPOC', 120), ('TESS-SPOC', None), ('QLP', None)]
#: the flux column of each provenance (QLP's KSPSAP is flattened by a
#:  spline, which removes rotation: its raw SAP flux is used)
FLUX = {'SPOC': 'PDCSAP_FLUX', 'TESS-SPOC': 'PDCSAP_FLUX', 'QLP': 'SAP_FLUX'}
#: the bins of the periodograms [days]
BIN = 30.0 / 1440
#: the shortest period searched [days]
PMIN = 0.1
#: the longest period searched, as a fraction of a sector's length
PMAX_FRACTION = 0.5
#: a peak that holds at least this fraction of the sector's highest power
#:  counts, as long as that highest power is itself above POWER_MIN
PEAK_FRACTION = 0.5
POWER_MIN = 0.1
#: the orbit of TESS [days]: scattered light and momentum dumps leave
#:  systematics at it and its half, where a match is only a note
TESS_ORBIT = 13.7
#: the multiples of a velocity period that rotation would also put power at
HARMONICS = ((1, 'P'), (0.5, 'P/2'), (1 / 3, 'P/3'), (2, '2P'))


# =============================================================================
# Define functions
# =============================================================================
def _mast(service: str, params: Dict[str, Any],
          timeout: float = 60.0) -> List[Dict[str, Any]]:
    """the rows of a MAST query"""
    req = dict(service=service, format='json', params=params)
    body = urllib.parse.urlencode(dict(request=json.dumps(req))).encode()
    with urllib.request.urlopen(urllib.request.Request(MAST, data=body),
                                timeout=timeout) as resp:
        return json.loads(resp.read()).get('data', [])


def tic_number(name: Union[str, int]) -> int:
    """
    The TIC number of a star: given as such ('TIC 274626553', 274626553),
    or from its name through Sesame

    :param name: str or int

    :return: int
    """
    text = str(name).strip()
    if text.upper().startswith('TIC'):
        text = text[3:].strip(' _-')
    if text.isdigit():
        return int(text)
    from koloa.archive import resolve
    tic = resolve(str(name)).get('tic')
    if not tic:
        raise ValueError(f'SIMBAD knows no TIC number of {name}')
    return int(tic.split()[-1])


def sectors_at(ra: Optional[float], dec: Optional[float],
               timeout: float = 30.0) -> Optional[List[int]]:
    """
    The sectors of TESS whose full frames hold a position (MAST's TESScut):
    none means TESS has not looked there (its sectors leave gaps: GJ 1214,
    at an ecliptic latitude of 28 deg, is in one)

    :param ra: float, right ascension [deg, J2000]
    :param dec: float, declination [deg, J2000]

    :return: list of int (empty: not observed), or None when it could not
             be asked (no position, no network)
    """
    if ra is None or dec is None:
        return None
    url = TESSCUT + '?' + urllib.parse.urlencode(dict(
        ra=f'{ra:.5f}', dec=f'{dec:.5f}', radius='0.01'))
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            rows = json.loads(resp.read()).get('results') or []
    except Exception as err:  # not asked: unknown
        log(f'TESS: the sectors at this position could not be asked ({err})',
            'warn')
        return None
    return sorted({int(row['sector']) for row in rows})


def planned_sectors(ra: Optional[float], dec: Optional[float]
                    ) -> Optional[List[int]]:
    """
    The sectors whose pointing holds a position, those to come too (the
    mission's pointing table, tess-point, when it is installed: pip install
    tess-point): GJ 1214, which no sector held before, is in sectors 118
    and 131

    :return: list of int, or None without tess-point (or a position)
    """
    if ra is None or dec is None:
        return None
    try:
        from tess_stars2px import tess_stars2px_function_entry as entry
    except ImportError:
        return None
    try:
        return sorted({int(sec) for sec in entry(0, ra, dec)[3] if sec > 0})
    except Exception:  # a help, not a need
        return None


def _choose(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """one observation per sector, the best provenance"""
    best = {}
    for row in rows:
        for rank, (prov, texp) in enumerate(PREFERENCE):
            if row.get('provenance_name') != prov:
                continue
            if texp is not None and row.get('t_exptime') != texp:
                continue
            sector = int(row['sequence_number'])
            if sector not in best or rank < best[sector][0]:
                best[sector] = (rank, row)
    return [best[sector][1] for sector in sorted(best)]


def _read(path: str, prov: str) -> Optional[Dict[str, np.ndarray]]:
    """the time [BJD - 2400000], flux and error [ppt] of a light curve"""
    from astropy.io import fits
    with fits.open(path) as hdul:
        tbl = hdul[1].data
        cols = tbl.columns.names
        fcol = FLUX[prov] if FLUX[prov] in cols else 'SAP_FLUX'
        ecol = next((col for col in (fcol + '_ERR', 'KSPSAP_FLUX_ERR',
                                     'SAP_FLUX_ERR') if col in cols), None)
        time_ = np.asarray(tbl['TIME'], float)
        flux = np.asarray(tbl[fcol], float)
        err = (np.asarray(tbl[ecol], float) if ecol
               else np.full(len(flux), np.nan))
        qual = (np.asarray(tbl['QUALITY']) if 'QUALITY' in cols
                else np.zeros(len(flux), int))
    good = (qual == 0) & np.isfinite(time_) & np.isfinite(flux)
    if good.sum() < 100:
        return None
    med = np.nanmedian(flux[good])
    err = np.where(np.isfinite(err), err, np.nanstd(flux[good]))
    # TESS time is BJD - 2457000
    return dict(time=time_[good] + 57000.0,
                flux=1e3 * (flux[good] / med - 1), err=1e3 * err[good] / med,
                column=fcol)


def light_curves(name: Union[str, int], folder: Optional[str] = None,
                 max_sectors: int = 30, refresh: bool = False,
                 timeout: float = 120.0) -> Dict[str, Any]:
    """
    The TESS light curves of a star, one per sector (see the module)

    :param name: str or int, the star's name or TIC number
    :param folder: str or None, where the FITS files are kept (CACHE/<TIC>
                   when None)
    :param max_sectors: int, the most recent sectors kept at most
    :param refresh: bool, download the files again
    :param timeout: float [s], of each download

    :return: dict, tic, sectors (list of dict: sector, provenance, exposure,
             file, column, time, flux, err) and, when MAST has no light
             curve, observed: the sectors whose full frames hold the star
             (empty: TESS has not looked at it yet; None: not asked), and
             planned: the sectors its pointings hold it in (tess-point)
    """
    tic = tic_number(name)
    folder = folder or os.path.join(CACHE, str(tic))
    log(f'TESS: asking MAST for the light curves of TIC {tic}')
    rows = _mast('Mast.Caom.Filtered', dict(
        columns='obsid,provenance_name,sequence_number,t_exptime',
        filters=[dict(paramName='target_name', values=[str(tic)]),
                 dict(paramName='dataproduct_type', values=['timeseries'])]))
    chosen = _choose(rows)[-max_sectors:]
    out = dict(tic=tic, sectors=[])
    if not chosen:
        # not observed at all, or in the full frames only?
        try:
            from koloa.archive import resolve
            ident = resolve(name if not str(name).strip().isdigit()
                            else f'TIC {name}')
            out['observed'] = sectors_at(ident.get('ra'), ident.get('dec'))
            if out['observed'] == []:
                out['planned'] = planned_sectors(ident.get('ra'),
                                                 ident.get('dec'))
        except Exception:  # its position not known: not asked
            out['observed'] = None
        if out['observed'] == []:
            log(f'TESS has not observed TIC {tic} yet: no sector at MAST '
                f'holds its position'
                + (f' (its pointings hold it in sectors '
                   f'{", ".join(map(str, out["planned"]))})'
                   if out.get('planned') else ''), 'warn')
        elif out['observed']:
            log(f'TESS: no light curve of TIC {tic} at MAST (in the full '
                f'frames of sectors '
                f'{", ".join(map(str, out["observed"]))} only)', 'warn')
        else:
            log(f'TESS: no light curve of TIC {tic} at MAST', 'warn')
        return out
    os.makedirs(folder, exist_ok=True)
    files = _mast('Mast.Caom.Products', dict(
        obsid=','.join(str(row['obsid']) for row in chosen)))
    for row in chosen:
        prov = row['provenance_name']
        names = [fl for fl in files if fl.get('obsID') == str(row['obsid'])
                 or fl.get('obsID') == row['obsid']]
        lcf = next((fl for fl in names
                    if fl['productFilename'].endswith(('_lc.fits',
                                                       '_llc.fits'))
                    and 'fast-lc' not in fl['productFilename']), None)
        if lcf is None:
            continue
        path = os.path.join(folder, lcf['productFilename'])
        if refresh or not os.path.exists(path):
            url = DOWNLOAD + urllib.parse.quote(lcf['dataURI'], safe=':/')
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                content = resp.read()
            with open(path, 'wb') as handle:
                handle.write(content)
        try:
            lc = _read(path, prov)
        except Exception as err:  # a broken file: the next sector
            log(f'TESS: sector {row["sequence_number"]} unreadable ({err})',
                'warn')
            continue
        if lc is None:
            continue
        lc.update(sector=int(row['sequence_number']), provenance=prov,
                  exposure=float(row['t_exptime']), file=path)
        out['sectors'].append(lc)
    log(f'TESS: {len(out["sectors"])} sectors of TIC {tic} ('
        + ', '.join(f'{lc["sector"]} {lc["provenance"]}'
                    for lc in out['sectors']) + ')', 'value')
    return out


def _clean_bin(time_: np.ndarray, flux: np.ndarray, err: np.ndarray,
               binsize: float = BIN):
    """flares and bad points clipped (against a 2-hour running median),
    then bins of binsize: time, flux, error"""
    from scipy.ndimage import median_filter
    order = np.argsort(time_)
    time_, flux, err = time_[order], flux[order], err[order]
    step = np.median(np.diff(time_)) if len(time_) > 1 else BIN
    width = max(3, int(round((2.0 / 24) / max(step, 1e-6))) | 1)
    resid = flux - median_filter(flux, size=width, mode='nearest')
    sig = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    keep = (resid < 4 * sig) & (resid > -6 * sig) if sig > 0 else \
        np.ones(len(flux), bool)
    time_, flux, err = time_[keep], flux[keep], err[keep]
    idx = np.floor((time_ - time_[0]) / binsize).astype(int)
    nbin = np.bincount(idx)
    full = nbin > 0
    btime = np.bincount(idx, time_)[full] / nbin[full]
    bflux = np.bincount(idx, flux)[full] / nbin[full]
    # the error of a bin: the scatter of its points, or their errors
    berr = np.sqrt(np.bincount(idx, err ** 2)[full]) / nbin[full]
    scatter = sig / np.sqrt(nbin[full]) if sig > 0 else berr
    return btime, bflux, np.maximum(berr, scatter)


def _amplitude(time_, flux, err, freq):
    """the semi-amplitude of the best sinusoid at a frequency [ppt]"""
    arg = 2 * np.pi * freq * (time_ - time_.mean())
    design = np.array([np.cos(arg), np.sin(arg), np.ones(len(time_))]).T
    wts = 1 / err ** 2
    coef = np.linalg.lstsq(design * np.sqrt(wts)[:, None],
                           flux * np.sqrt(wts), rcond=None)[0]
    return float(np.hypot(coef[0], coef[1]))


def periodicity(lcs: Dict[str, Any], periods: Sequence[float] = (),
                pmin: float = PMIN, npeaks: int = 3,
                oversample: int = 5) -> Dict[str, Any]:
    """
    The periodicity of each sector, and the velocity periods against it

    :param lcs: dict, from light_curves()
    :param periods: list of float, the periods of the velocity signals [d]
    :param pmin: float, the shortest period searched [d]
    :param npeaks: int, the strongest peaks kept per sector
    :param oversample: int, frequencies per 1/T

    :return: dict, tic, sectors (per sector: sector, provenance, span, rms
             [ppt], pmax, freq, power, peaks (period, power, amplitude)),
             checks (per velocity period: the matches, per multiple, of the
             sectors whose peaks fall there, a status and a summary), and a
             summary
    """
    out = dict(tic=lcs.get('tic'), sectors=[], checks=[])
    for lc in lcs.get('sectors', []):
        btime, bflux, berr = _clean_bin(lc['time'], lc['flux'], lc['err'])
        if len(btime) < 50:
            continue
        span = float(btime[-1] - btime[0])
        pmax = PMAX_FRACTION * span
        fmin, fmax = 1 / pmax, 1 / pmin
        freq = np.arange(fmin, fmax, 1 / (oversample * span))
        power = gls(btime, bflux, berr, freq)
        # the local maxima, strongest first
        isl = np.r_[False, (power[1:-1] > power[:-2])
                    & (power[1:-1] >= power[2:]), False]
        order = np.argsort(power[isl])[::-1][:npeaks]
        # also kept: every local maximum, for the matches
        peaks = [dict(period=float(1 / freq[isl][ii]),
                      power=float(power[isl][ii]),
                      amplitude=_amplitude(btime, bflux, berr,
                                           freq[isl][ii]))
                 for ii in order]
        out['sectors'].append(dict(
            sector=lc['sector'], provenance=lc['provenance'],
            column=lc.get('column'), span=span, pmax=pmax,
            rms=float(np.std(bflux)), npoints=int(len(btime)),
            time=btime, flux=bflux, err=berr, freq=freq, power=power,
            peaks=peaks, maxima=np.where(isl)[0]))
    secs = out['sectors']
    for period in periods:
        matches, tested = [], set()
        for mult, label in HARMONICS:
            ptest = period * mult
            for sec in secs:
                if not pmin <= ptest <= sec['pmax']:
                    continue
                tested.add(label)
                ftest = 1 / ptest
                # a peak (a local maximum, not the flank of another)
                #   within one peak width, 1/T, of the tested period
                maxi = sec['maxima']
                near = maxi[np.abs(sec['freq'][maxi] - ftest)
                            < 0.5 / sec['span']]
                top = float(np.max(sec['power']))
                if top < POWER_MIN or not len(near):
                    continue
                inear = near[np.argmax(sec['power'][near])]
                pnear = float(sec['power'][inear])
                if pnear >= PEAK_FRACTION * top:
                    pmatch = 1 / sec['freq'][inear]
                    width = pmatch ** 2 * 0.5 / sec['span']
                    orbital = any(abs(pmatch - TESS_ORBIT / nn) < width
                                  for nn in (1, 2))
                    matches.append(dict(orbital=orbital,
                        multiple=label, sector=sec['sector'],
                        period=float(1 / sec['freq'][inear]), power=pnear,
                        top=top, strongest=bool(pnear >= top * 0.999),
                        amplitude=_amplitude(sec['time'], sec['flux'],
                                             sec['err'],
                                             sec['freq'][inear])))
        nsec = len(secs)
        strong = [mt for mt in matches if mt['strongest']]
        if not secs:
            status, summary = 'info', 'no TESS light curve'
        elif not tested:
            status = 'info'
            summary = (f'P = {period:.2f} d and its harmonics are beyond '
                       f'what a TESS sector tests (up to '
                       f'{max(sc["pmax"] for sc in secs):.1f} d)')
        elif not matches:
            status = 'pass'
            summary = (f'no photometric peak at {", ".join(sorted(tested))}'
                       f' in {nsec} sector{"s" if nsec > 1 else ""}')
        else:
            by = {}
            for mt in matches:
                by.setdefault(mt['multiple'], []).append(mt)
            parts = []
            for label, mts in by.items():
                amps = np.median([mt['amplitude'] for mt in mts])
                parts.append(f'{label} ({np.median([mt["period"] for mt in mts]):.2f} d) '
                             f'in {len(set(mt["sector"] for mt in mts))} of '
                             f'{nsec} sectors, {amps:.2f} ppt')
            # the strongest peak of several sectors, or of the only one,
            #   away from the orbit of TESS
            real = [mt for mt in strong if not mt['orbital']]
            nstrong = len(set(mt['sector'] for mt in real))
            status = 'flag' if (nstrong >= 2 or (nstrong and nsec == 1)) \
                else 'info'
            summary = 'photometric peak at ' + '; '.join(parts)
            if any(mt['orbital'] for mt in matches):
                summary += (f' (at the orbit of TESS, {TESS_ORBIT} d, or its'
                            f' half, where its light curves keep '
                            f'systematics: a note, not a flag)')
        out['checks'].append(dict(period=float(period), status=status,
                                  summary=summary, matches=matches,
                                  tested=sorted(tested)))
    tops = [sec['peaks'][0] for sec in secs if sec['peaks']]
    if tops:
        pers = np.array([pk['period'] for pk in tops])
        out['summary'] = (f'{len(secs)} sectors; strongest photometric '
                          f'period per sector: median {np.median(pers):.2f} d'
                          f' (range {pers.min():.2f} to {pers.max():.2f} d), '
                          f'amplitude median '
                          f'{np.median([pk["amplitude"] for pk in tops]):.2f}'
                          f' ppt; periods tested up to '
                          f'{max(sc["pmax"] for sc in secs):.1f} d')
    else:
        out['summary'] = 'no TESS light curve'
    return out


def light(res: Dict[str, Any]) -> Dict[str, Any]:
    """
    A result of periodicity() without its arrays (for a JSON file)

    :param res: dict, from periodicity()

    :return: dict
    """
    arrays = ('time', 'flux', 'err', 'freq', 'power', 'maxima')
    out = {key: val for key, val in res.items() if key != 'sectors'}
    out['sectors'] = [{key: val for key, val in sec.items()
                       if key not in arrays} for sec in res.get('sectors', [])]
    return out


def table(res: Dict[str, Any]) -> List[str]:
    """
    One line per sector: its source, its scatter and its strongest periods

    :param res: dict, from periodicity()

    :return: list of str
    """
    return [f'{sec["sector"]:>4d}  {sec["provenance"]:<9s} rms '
            f'{sec["rms"]:5.2f} ppt  '
            + '  '.join(f'{pk["period"]:6.2f} d ({pk["power"]:.2f}, '
                        f'{pk["amplitude"]:.2f} ppt)' for pk in sec['peaks'])
            for sec in res.get('sectors', [])]


def check(name: Union[str, int], periods: Sequence[float],
          folder: Optional[str] = None, **kwargs) -> Dict[str, Any]:
    """
    light_curves() then periodicity(), in one call

    :param name: str or int, the star's name or TIC number
    :param periods: list of float, the velocity periods [d]
    :param folder: str or None, where the FITS files are kept
    :param kwargs: passed to periodicity()

    :return: dict, from periodicity()
    """
    start = time.time()
    lcs = light_curves(name, folder=folder)
    res = periodicity(lcs, periods, **kwargs)
    log(f'TESS: {res["summary"]} ({time.time() - start:.0f} s)')
    return res


def figure(res: Dict[str, Any], periods: Sequence[float] = (),
           title: Optional[str] = None):
    """
    The periodograms of the sectors (each thin, their median thick) with the
    velocity periods and their harmonics marked, and the light curve of the
    sector whose strongest peak is highest, with its sinusoid

    :param res: dict, from periodicity()
    :param periods: list of float, the velocity periods to mark [d]
    :param title: str or None

    :return: matplotlib figure
    """
    from koloa import plotting as kplot
    plt, C = kplot.plt, kplot.C
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.6))
    secs = res['sectors']
    if not secs:
        ax1.text(0.5, 0.5, 'no TESS light curve', ha='center',
                 transform=ax1.transAxes)
        return fig
    grid = np.geomspace(PMIN, max(sc['pmax'] for sc in secs), 2000)
    stack = []
    for sec in secs:
        per = 1 / sec['freq'][::-1]
        ax1.plot(per, sec['power'][::-1], color=C['muted'], lw=0.6,
                 alpha=0.5)
        stack.append(np.interp(grid, per, sec['power'][::-1], left=np.nan,
                               right=np.nan))
    with np.errstate(all='ignore'):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            median = np.nanmedian(np.array(stack), axis=0)
    ax1.plot(grid, median, color=C['koloa'], lw=1.6,
             label=f'median of {len(secs)} sectors')
    for nn in (1, 2):
        ax1.axvspan(TESS_ORBIT / nn * 0.95, TESS_ORBIT / nn * 1.05,
                    color=C['muted'], alpha=0.12, lw=0,
                    label='the orbit of TESS and its half' if nn == 1
                    else None)
    for ip, period in enumerate(periods):
        for mult, label in HARMONICS:
            ptest = period * mult
            if PMIN <= ptest <= grid[-1]:
                ax1.axvline(ptest, color=C['outlier'], lw=1.0 if mult == 1
                            else 0.7, ls='-' if mult == 1 else (0, (3, 2)))
                # the labels of each period on their own line
                ax1.text(ptest, 1.0 - 0.07 * ip, f' {label}'
                         if len(periods) == 1 else f' {label}$_{ip + 1}$',
                         color=C['outlier'], fontsize=7, va='top',
                         transform=ax1.get_xaxis_transform())
    ax1.set_xscale('log')
    kplot.plain_log_ticks(ax1, 'x')
    ax1.set_xlabel('period [d]')
    ax1.set_ylabel('GLS power')
    ax1.set_xlim(PMIN, grid[-1])
    ax1.set_ylim(bottom=0)
    ax1.legend(loc='upper left')
    # the sector whose strongest peak is the highest
    best = max(secs, key=lambda sc: sc['peaks'][0]['power']
               if sc['peaks'] else 0)
    ax2.errorbar(best['time'], best['flux'], best['err'], fmt='o', ms=2,
                 color=C['koloa'], elinewidth=0.5, alpha=0.7)
    if best['peaks']:
        pk = best['peaks'][0]
        arg = 2 * np.pi * (best['time'] - best['time'].mean()) / pk['period']
        design = np.array([np.cos(arg), np.sin(arg),
                           np.ones(len(arg))]).T
        coef = np.linalg.lstsq(design, best['flux'], rcond=None)[0]
        tfine = np.linspace(best['time'].min(), best['time'].max(), 2000)
        afine = 2 * np.pi * (tfine - best['time'].mean()) / pk['period']
        ax2.plot(tfine, coef[0] * np.cos(afine) + coef[1] * np.sin(afine)
                 + coef[2], color=C['text'], lw=1.2,
                 label=f'sector {best["sector"]}: {pk["period"]:.2f} d, '
                       f'{pk["amplitude"]:.2f} ppt')
        ax2.legend(loc='upper left')
    ax2.set_xlabel('time [BJD - 2400000]')
    ax2.set_ylabel(f'flux [ppt] ({best["provenance"]}, 30 min)')
    if title:
        ax1.set_title(title, loc='left')
    fig.tight_layout()
    return fig


# =============================================================================
# End of code
# =============================================================================

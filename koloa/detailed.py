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
from koloa import plotting as kplot
from koloa.archive import MATCH, compare, conjunction, known_planets, resolve
from koloa.data import RVData, merge, robust_std
from koloa.diagnostics import duck_test
from koloa.fip import inflate_to_fit, oafip
from koloa.fit import RVModel, mcmc_orbits
from koloa.gpcheck import (figure_fold, figure_sequence, figure_whitened,
                           gp_signals)
from koloa.log import log
from koloa.outliers import NAMES, explain, instrument_list
from koloa.periodogram import find_peaks, frequency_grid, oap
from koloa.utils import blas_threads

# =============================================================================
# Define variables
# =============================================================================
#: a signal of the FIP is detected below this
THRESHOLD = 0.01
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
def _read(source: Union[str, RVData], name: Optional[str]) -> RVData:
    """the series of a file (or the series given), its instrument named; a
    CSV written by koloa.dace.fetch is read by koloa.dace.rvdata (one
    instrument per era)"""
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
    if data.instruments == ['inst']:
        lname, _ = instrument_list(data, 'inst')
        label = NAMES.get(lname)
        if label:
            data.inst = np.array([label] * data.n)
            data.zero_point = {label: data.zero_point.get('inst', 0.0)}
    return data


def _target(data: RVData, target: Optional[str]) -> Optional[str]:
    """the name of the star: given, from the OBJECT column, or the series"""
    if target:
        return target
    for col in ('OBJECT', 'object', 'target', 'TARGET'):
        if col in data.meta and data.meta[col].dtype.kind in 'USO':
            names, counts = np.unique(data.meta[col], return_counts=True)
            best = str(names[np.argmax(counts)]).strip()
            if best:
                return best
    return data.name


def dace_names(ident: Dict[str, Any], target: str) -> List[str]:
    """the names DACE may know a star by: its catalogue names without
    spaces (HD69830, GJ687, HIP86162...), then as SIMBAD writes them"""
    out = []
    for key in ('hd', 'gj', 'hip', 'main'):
        if ident.get(key):
            out.append(ident[key].replace('NAME ', '').replace(' ', ''))
    out.append(target.replace(' ', ''))
    for alias in ident.get('aliases', []):
        if alias.split()[0] in ('GJ', 'Gl', 'HD', 'HIP', 'TOI', 'K2',
                                'Kepler', 'LHS', 'Wolf', 'Ross'):
            out.append(alias.replace(' ', ''))
    return list(dict.fromkeys(name for name in out if name))


def fetch_dace(names: Sequence[str], folder: str, exclude: Sequence[str] = (),
               times: Optional[np.ndarray] = None, refresh: bool = False
               ) -> Optional[RVData]:
    """
    The public velocities of a star on DACE, trying its names in turn

    :param names: list of str, the names to try
    :param folder: str, where the CSV is kept
    :param exclude: list of str, instruments to leave out (those of the file:
                    DACE's pipeline velocities of the same spectra)
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
        keep = [inst for inst in data.instruments
                if not any(inst.upper().startswith(ex.upper())
                           for ex in exclude)]
        mask = np.isin(data.inst, keep)
        if times is not None and len(times):
            # the same spectra as the file's: the same time, within a minute
            nearest = np.min(np.abs(data.time[:, None] - np.asarray(times)[None, :]),
                             axis=1)
            same = nearest < 1.0 / 1440
            if np.any(same & mask):
                log(f'DACE, {name}: {int(np.sum(same & mask))} exposures that '
                    f'the file already has, left out', 'warn')
            mask &= ~same
        if not np.any(mask):
            log(f'DACE, {name}: nothing that the file does not have', 'warn')
            return None
        return data.select(mask)
    return None


def _fip(data, fit, kmax, nsweep, nburn, seed, label):
    """the FIP of a series, the errors of each instrument inflated to a fit"""
    inflated, info = inflate_to_fit(fit)
    log(f'{label}: errors inflated by ' + ', '.join(
        f'{inst} {val:.2f}' for inst, val in info['inflation'].items())
        + ' m/s', 'value')
    with blas_threads(1):
        res = oafip(inflated, kmax=kmax, outliers='both', nsweep=nsweep,
                    nburn=nburn, nchains=2, seed=seed, progress=False)
    found = [pk for pk in res.peaks if pk['fip'] < THRESHOLD]
    log(f'{label}: P(k) = ' + ', '.join(f'{val:.2f}' for val in res.pk)
        + '; FIP < 1 %: ' + (', '.join(f'{pk["period"]:.4f} d '
                                       f'({pk["fip"]:.1e})' for pk in found)
                             or 'none'), 'value')
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


def detailed_analysis(source: Union[str, RVData], outdir: str = 'koloa_detailed',
                      name: Optional[str] = None, target: Optional[str] = None,
                      archive: bool = True, dace: bool = True,
                      dace_folder: Optional[str] = None,
                      literature: Optional[Sequence[Union[str, RVData]]]
                      = None, vizier: bool = True,
                      periods: Optional[Sequence[float]] = None,
                      gp: bool = True, gp_workers: int = 4,
                      gp_nsim: int = 20000,
                      refresh: bool = False, kmax: int = 4,
                      nsweep: int = 1500, nburn: Optional[int] = None,
                      pmin: float = 1.1, pmax: Optional[float] = None,
                      mcmc: bool = False, nsteps: int = 4000,
                      duck: bool = True, keys: Union[str, Sequence] = 'auto',
                      style: str = 'paper', seed: int = 1,
                      latex: bool = True, tess: bool = True
                      ) -> Dict[str, Any]:
    """
    Everything koloa can say about a star (see the module's docstring)

    :param source: str or RVData, an LBL .rdb, a csv, or a series
    :param outdir: str, where the report, the summary and the figures go
    :param name: str or None, the name of the series
    :param target: str or None, the name of the star for SIMBAD, the
                   archive and DACE (from the OBJECT column when None)
    :param archive: bool, ask the NASA Exoplanet Archive
    :param dace: bool, ask DACE for more velocities
    :param dace_folder: str or None, where the DACE CSV is kept (outdir when
                        None)
    :param literature: list or None, published velocities to add: files (a
                       VizieR .dat, time velocity error instrument, without
                       a header; or a csv or .rdb with named columns) or
                       series (see koloa.literature.read)
    :param vizier: bool, add the velocities published with the known
                   planets, from VizieR (the tables of the papers of their
                   solutions in the archive)
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

    :return: dict, every result (and outdir/<star>_report.txt,
             <star>_summary.json, <star>_report.tex and .pdf, and the
             figures)
    """
    start = time.time()
    os.makedirs(outdir, exist_ok=True)
    kplot.set_style(style)
    out: Dict[str, Any] = dict(figures=[])
    # 1. the series
    data = _read(source, name)
    star = _target(data, target)
    log(f'koloa, detailed: {star} ({data.name}), {data.n} exposures in '
        f'{data.nseq} visits over {data.baseline:.0f} d, '
        f'{", ".join(data.instruments)}')
    # 2. who the star is, and what is known of it
    ident, known = None, dict(host=None, planets=[])
    if archive or dace:
        try:
            ident = resolve(star)
            log(f'SIMBAD: {ident["main"]} ('
                + ', '.join(val for val in (ident['gj'], ident['hd'],
                                            ident['hip'], ident['tic'])
                            if val) + ')', 'value')
        except (OSError, ValueError) as err:
            log(f'SIMBAD could not resolve {star}: {err}', 'warn')
    if archive and ident is not None:
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
    out['ident'], out['known'] = ident, known
    # 3. more velocities: DACE, the ones given, the ones on VizieR
    sources = [dict(kind='file', label=(os.path.basename(source)
                                        if isinstance(source, str)
                                        else data.name),
                    n=int(data.n), instruments=_counts(data), note='')]
    if not dace:
        sources.append(dict(kind='DACE', label='DACE', n=0, instruments={},
                            note='not asked'))
    elif ident is None:
        sources.append(dict(kind='DACE', label='DACE', n=0, instruments={},
                            note='not asked: SIMBAD did not resolve the star'))
    else:
        names = dace_names(ident, star)
        more = fetch_dace(names, dace_folder or outdir,
                          exclude=data.instruments, times=data.time,
                          refresh=refresh)
        if more is not None:
            data = merge([data, more], name=data.name)
            log(f'with DACE: {data.n} exposures in {data.nseq} visits, '
                + ', '.join(f'{inst} {np.sum(data.inst == inst)}'
                            for inst in data.instruments), 'value')
        sources.append(dict(
            kind='DACE', label='DACE', n=int(more.n) if more else 0,
            instruments=_counts(more) if more else {},
            note=('public velocities of the other instruments' if more else
                  'nothing added (unreachable, or no public velocity the '
                  'file does not have) under ' + ', '.join(names))))
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
    out['data'], out['sources'] = data, sources
    multi = len(data.instruments) > 1
    seq_jitter = 'instrument' if multi else None
    nburn = nburn if nburn is not None else max(300, nsweep // 5)
    # 4. the noise without planets, and the FIP
    with blas_threads(1):
        noise = RVModel(data, [], likelihood='mixture', unit='both', trend=1,
                        seq_jitter=seq_jitter).fit(nstart=2, quiet=True)
    fip1, info1 = _fip(data, noise, kmax, nsweep, nburn, seed,
                       'FIP, noise without planets')
    found = sorted(pk['period'] for pk in fip1.peaks
                   if pk['fip'] < THRESHOLD)
    # 5. the signals, the known planets and the periods given, fitted
    #   together (a known planet the FIP does not find is tested all the
    #   same), against the known planets
    tests = [(per, 'FIP') for per in found]
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
        planets = [dict(period=per, period_range=(0.98 * per, 1.02 * per))
                   for per, _ in chosen]
        with blas_threads(1):
            if mcmc:
                fit = mcmc_orbits(data, planets, likelihood='mixture',
                                  unit='both', trend=1, seq_jitter=seq_jitter,
                                  nsteps=nsteps, nburn=nsteps // 3, nstart=2,
                                  seed=seed, quiet=True)
            else:
                fit = RVModel(data, planets, likelihood='mixture',
                              unit='both', trend=1,
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
        fip2, info2 = _fip(data, fit, kmax, nsweep, nburn, seed + 1,
                           'FIP, noise with the planets')
    else:
        log('no interval with FIP < 1 %, no known planet and no period '
            'given: nothing to fit', 'warn')
    for orb in orbits:
        # the FIP of the best interval that holds each signal
        orb['fip'] = fip2.fip_containing(orb['P'][0], 1 / data.baseline)
    out.update(fit=fit, orbits=orbits, fip_first=fip1, fip_second=fip2)
    # 7. the activity indicators, and the rotation of the star
    indic = _indicators(data, pmin, pmax)
    out['indicators'] = indic
    for orb in orbits:
        orb['activity'] = _activity_match(orb['P'][0], indic)
        orb['rotation'] = _rotation_match(orb['P'][0], known.get('star') or {})
    # 7b. the signals against a GP of the activity
    gpres = None
    if gp:
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
    if tess:
        from koloa import tess as ktess
        try:
            lcs = ktess.light_curves((ident or {}).get('tic') or star)
            phot = ktess.periodicity(lcs, [orb['P'][0] for orb in orbits])
        except Exception as err:  # no network, or an unknown star
            log(f'TESS: {err}', 'warn')
    out['tess'] = phot
    ducks = {}
    if duck:
        for orb in orbits:
            try:
                report = duck_test(data, orb['P'][0], fipres=fip2, gp=False,
                                   unit='both', quiet=True,
                                   tess=lcs if lcs is not None else False)
                ducks[f'{orb["P"][0]:.4f}'] = report
            except Exception as err:  # the test is a help, not a stop
                log(f'duck test at {orb["P"][0]:.4f} d: {err}', 'warn')
    out['duck'] = ducks
    # 9. why each outlier is one
    why = explain(fit, keys=keys)
    out['outliers'] = why
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
        keep(kplot.periodograms(
            res.freq, fips=dict(koloa=res), mark=marks,
            title=f'{star}: FIP ({label} pass; known planets dashed)'),
            f'fip_{label}')
    for ip, orb in enumerate(orbits):
        keep(kplot.phase(
            fit, planet=ip, level='point',
            title=f'{star}, {orb["P"][0]:.3f} d: K = {orb["K"][0]:.2f} m/s'),
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
        except Exception as err:  # a figure is a help, not a stop
            log(f'duck test figures at {orb["P"][0]:.4f} d: {err}', 'warn')
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
    if phot is not None:
        text += ('\n\nTESS photometry (TIC ' + str(phot.get('tic')) + ')\n  '
                 + phot['summary'] + '\n'
                 + ''.join(f'  {chk["period"]:.4f} d: [{chk["status"]}] '
                           f'{chk["summary"]}\n' for chk in phot['checks'])
                 + '\n'.join('  ' + line for line in ktess.table(phot)))
    with open(os.path.join(outdir, f'{safe}_report.txt'), 'w') as handle:
        handle.write(text + '\n')
    out['report'] = text
    # 11. the report, in LaTeX and PDF
    paths = {}
    if latex:
        from koloa.latex import compile_pdf, detailed_report
        tex = os.path.join(outdir, f'{safe}_report.tex')
        rep = dict(
            star=star, source=source if isinstance(source, str) else None,
            data=data, ident=ident, known=known, sources=sources,
            fip_first=fip1, fip_second=fip2,
            inflation_first=info1.get('inflation', {}),
            inflation_second=info2.get('inflation', {}), orbits=orbits,
            indicators=indic, duck=ducks, outliers=why, figures=named,
            gp=gpsum, threshold=THRESHOLD, runtime=time.time() - start,
            tess=phot,
            settings=dict(archive=archive, dace=dace, vizier=vizier,
                          literature=len(literature or []),
                          periods=', '.join(f'{per}' for per in periods or [])
                          or '--', gp=gp, gp_nsim=gp_nsim, kmax=kmax,
                          nsweep=nsweep, nburn=nburn, pmin=pmin,
                          pmax=pmax or 'the baseline', mcmc=mcmc,
                          nsteps=nsteps if mcmc else '--', seed=seed,
                          style=style, tess=tess),
            files=sorted(set(os.listdir(outdir)) | {
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
                   report_pdf=paths.get('pdf'))
    if phot is not None:
        summary['tess'] = ktess.light(phot)
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
                                                     for val in res.pk)]
        for pk in res.peaks:
            flag = 'DETECTED' if pk['fip'] < THRESHOLD else ''
            lines.append(f'  {pk["period"]:12.4f} d  FIP {pk["fip"]:.2e}  '
                         f'{flag}')
    lines += ['', 'Signals']
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
        report = ducks.get(f'{orb["P"][0]:.4f}')
        if report is not None:
            lines.append(f'    duck test: {report.verdict}')
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
        outliers=why.to_dict(), figures=figs)


# =============================================================================
# End of code
# =============================================================================

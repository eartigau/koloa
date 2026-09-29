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
   only, and the analysis goes on without it;
4. fits the noise of every instrument without planets (a white and a visit
   jitter each), inflates the errors of each instrument to it, and runs the
   outlier-aware FIP with up to kmax signals (koloa.fip.inflate_to_fit,
   koloa.oafip);
5. fits every interval with FIP < 1 % together (each period free within its
   interval; every exposure and visit may be an outlier; by MCMC with
   mcmc=True), and compares each with the known planets (every published
   K, the most recent first);
6. runs the FIP again, the errors inflated to the fit with the planets
   (the first pass had the planets in its noise);
7. looks at the activity indicators: the peaks of their outlier-aware
   periodograms, and whether a signal of the velocities sits at one of them
   (or at twice or half of it);
8. puts every signal to the duck test;
9. explains every outlier (koloa.outliers: the header keywords, the
   indicators and the error bars that are off for it);
10. writes a text report, a JSON summary and the figures.

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

from koloa import plotting as kplot
from koloa.archive import MATCH, compare, known_planets, resolve
from koloa.data import RVData, merge
from koloa.diagnostics import duck_test
from koloa.fip import inflate_to_fit, oafip
from koloa.fit import RVModel, mcmc_orbits
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


def _indicators(data: RVData, pmin: float, pmax: float,
                names: Union[str, Sequence[str]] = 'auto',
                npeaks: int = 3) -> Dict[str, Any]:
    """the strongest peaks of the outlier-aware periodogram of each activity
    indicator, per instrument that has it"""
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
        if series.n < 20:
            continue
        freq = frequency_grid(series.time, pmin, pmax, 10)
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


def detailed_analysis(source: Union[str, RVData], outdir: str = 'koloa_detailed',
                      name: Optional[str] = None, target: Optional[str] = None,
                      archive: bool = True, dace: bool = True,
                      dace_folder: Optional[str] = None,
                      refresh: bool = False, kmax: int = 4,
                      nsweep: int = 1500, nburn: Optional[int] = None,
                      pmin: float = 1.1, pmax: Optional[float] = None,
                      mcmc: bool = False, nsteps: int = 4000,
                      duck: bool = True, keys: Union[str, Sequence] = 'auto',
                      style: str = 'paper', seed: int = 1) -> Dict[str, Any]:
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
    :param refresh: bool, ask DACE and the archive again
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
    :param style: str, paper (PDF) or web (SVG)
    :param seed: int, the seed

    :return: dict, every result (and outdir/<star>_report.txt,
             <star>_summary.json and the figures)
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
    # 3. more velocities from DACE
    if dace and ident is not None:
        more = fetch_dace(dace_names(ident, star), dace_folder or outdir,
                          exclude=data.instruments, times=data.time,
                          refresh=refresh)
        if more is not None:
            data = merge([data, more], name=data.name)
            log(f'with DACE: {data.n} exposures in {data.nseq} visits, '
                + ', '.join(f'{inst} {np.sum(data.inst == inst)}'
                            for inst in data.instruments), 'value')
    out['data'] = data
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
    # 5. the signals, fitted together, against the known planets
    fit, orbits, fip2, info2 = noise, [], fip1, info1
    if found:
        planets = [dict(period=per, period_range=(0.98 * per, 1.02 * per))
                   for per in found]
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
        mstar = (known.get('star') or {}).get('mass')
        for orb in fit.orbits(mstar=mstar):
            orbits.append(compare({key: [float(val) for val in orb[key]]
                                   for key in ('P', 'K', 'e', 'msini')},
                                  known, MATCH))
        # 6. the FIP again, with the noise of the fit with the planets
        fip2, info2 = _fip(data, fit, kmax, nsweep, nburn, seed + 1,
                           'FIP, noise with the planets')
    else:
        log('no interval with FIP < 1 %: nothing to fit', 'warn')
    out.update(fit=fit, orbits=orbits, fip_first=fip1, fip_second=fip2)
    # 7. the activity indicators
    top = pmax or data.baseline
    indic = _indicators(data, pmin, top)
    out['indicators'] = indic
    for orb in orbits:
        orb['activity'] = _activity_match(orb['P'][0], indic)
    # 8. the duck test
    ducks = {}
    if duck:
        for orb in orbits:
            try:
                report = duck_test(data, orb['P'][0], fipres=fip2, gp=False,
                                   unit='both', quiet=True)
                ducks[f'{orb["P"][0]:.4f}'] = report
            except Exception as err:  # the test is a help, not a stop
                log(f'duck test at {orb["P"][0]:.4f} d: {err}', 'warn')
    out['duck'] = ducks
    # 9. why each outlier is one
    why = explain(fit, keys=keys)
    out['outliers'] = why
    # 10. the figures, the report and the summary
    safe = data.name.replace(' ', '_')
    figs = out['figures']
    marks = [pl['P'] for pl in known.get('planets', []) if pl.get('P')]
    figs.append(kplot.savefig(kplot.timeseries(
        data, fit.outlier_prob, fit=fit if orbits else None,
        title=f'{star}: every instrument, each exposure coloured by its '
              f'outlier probability'), os.path.join(outdir, f'{safe}_rv')))
    for label, res in (('first', fip1), ('second', fip2)):
        figs.append(kplot.savefig(kplot.periodograms(
            res.freq, fips=dict(koloa=res), mark=marks,
            title=f'{star}: FIP ({label} pass; known planets dashed)'),
            os.path.join(outdir, f'{safe}_fip_{label}')))
    for ip, orb in enumerate(orbits):
        figs.append(kplot.savefig(kplot.phase(
            fit, planet=ip, level='point',
            title=f'{star}, {orb["P"][0]:.3f} d: K = {orb["K"][0]:.2f} m/s'),
            os.path.join(outdir, f'{safe}_phase_{ip}')))
    if indic:
        figs.append(kplot.savefig(_indicator_figure(indic, orbits, star),
                                  os.path.join(outdir, f'{safe}_indicators')))
    if why.units:
        figs.append(kplot.savefig(kplot.outlier_keys(why),
                                  os.path.join(outdir, f'{safe}_outliers')))
    text = _report(star, data, ident, known, fip1, fip2, orbits, indic,
                   ducks, why, time.time() - start)
    with open(os.path.join(outdir, f'{safe}_report.txt'), 'w') as handle:
        handle.write(text + '\n')
    summary = _summary(star, data, ident, known, fip1, fip2, orbits, indic,
                       ducks, why, figs)
    with open(os.path.join(outdir, f'{safe}_summary.json'), 'w') as handle:
        json.dump(summary, handle, indent=1, default=_json)
    out['report'] = text
    log(f'done in {(time.time() - start) / 60:.1f} min: the report, the '
        f'summary and {len(figs)} figures in {outdir}')
    return out


def _json(obj):
    """numbers and arrays for json.dump"""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return str(obj)


def _indicator_figure(indic: Dict[str, Any], orbits, star: str):
    """the outlier-aware periodogram of each indicator, the signals of the
    velocities marked"""
    fig, axes = kplot.plt.subplots(len(indic), 1, sharex=True,
                                   figsize=(7.0, 1.3 * len(indic) + 0.8))
    axes = np.atleast_1d(axes)
    for ax, (name, res) in zip(axes, indic.items()):
        ax.plot(1 / res['freq'], res['dlnl'], color=kplot.C['koloa'], lw=0.8)
        for orb in orbits:
            ax.axvline(orb['P'][0], color=kplot.C['muted'], ls='--', lw=0.8)
        ax.set_ylabel(r'$\Delta\ln L$', fontsize=8)
        ax.text(0.01, 0.85, name, transform=ax.transAxes, fontsize=8,
                color=kplot.C['text'])
    axes[-1].set_xscale('log')
    axes[-1].set_xlabel('period [d]')
    axes[0].set_title(f'{star}: activity indicators (the signals of the '
                      f'velocities dashed)', loc='left', fontsize=9)
    fig.tight_layout()
    return fig


def _report(star, data, ident, known, fip1, fip2, orbits, indic, ducks, why,
            runtime) -> str:
    """the report, in words"""
    lines = [f'koloa, detailed analysis of {star}', '=' * 70, '']
    lines.append(f'{data.n} exposures in {data.nseq} visits over '
                 f'{data.baseline:.0f} d: ' + ', '.join(
                     f'{inst} {np.sum(data.inst == inst)}'
                     for inst in data.instruments))
    if ident:
        lines.append(f'SIMBAD: {ident["main"]}; ' + ', '.join(
            val for val in (ident['gj'], ident['hd'], ident['hip'],
                            ident['tic'], ident['gaia_dr3']) if val))
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
        if orb.get('activity'):
            lines.append('    ACTIVITY at this period: '
                         + '; '.join(orb['activity']))
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

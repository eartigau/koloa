#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Which planets a series could have found, by the rule that decides: a
detection map made with the FIP.

Circular planets are injected into the series (the signals found taken out
first) and looked for exactly as koloa decides on a signal: the
outlier-aware FIP, with the same GP of the activity inside it, a planet
counting as found when the FIP of its period OR any of its aliases is below
the threshold (1 %). A blind periodogram search (koloa.completeness) is
harsher: it counts a planet found on its daily alias as missed, has no GP
for the activity, and sets its threshold for the whole grid; on G 203-42 it
found 8 of 12 planets of K = 5 m/s near 7 d that this rule finds 11 times.

A FIP costs tens of seconds, so the map is adaptive rather than a grid: in
each period band, rounds of injections, each at the K where a logistic curve
of the detection probability against ln K, fitted to every injection of the
band so far, puts 50 % and 90 %. The result per band is the K detected 50
and 90 % of the time, with the 16-84 % range of a bootstrap of the
injections.

Created on 2026-09-30

@author: artigau
"""
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from scipy.optimize import minimize

from koloa.data import RVData
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: the FIP below which a planet (its period or any of its aliases) is found
THRESHOLD = 0.01
#: the logistic slope prior: ln K from 50 to 73 % of detection, a gaussian
#:  on its log (keeps the fit finite when the injections are all found or
#:  all missed)
SLOPE_PRIOR = (np.log(0.3), 1.0)


# =============================================================================
# Define functions
# =============================================================================
def _one(job: Dict[str, Any]) -> Dict[str, Any]:
    """one injection: the planet added, the FIP run, found or not"""
    from koloa.completeness import inject
    from koloa.fip import oafip
    from koloa.utils import blas_threads
    data = inject(job['series'], job['period'], job['amp'], job['phase'])
    with blas_threads(1):
        res = oafip(data, kmax=job['kmax'], nsweep=job['nsweep'],
                    nburn=job['nburn'], nchains=1, seed=job['seed'],
                    progress=False, gp=job['gp'])
    width = 1.0 / data.baseline
    fam = res.family_containing(job['period'], width)
    fam = res.fip_containing(job['period'], width) if fam is None else fam
    out = {key: val for key, val in job.items() if key not in ('series', 'gp')}
    out.update(fip=float(fam), found=bool(fam < job['threshold']))
    return out


def _logistic_fit(lnk: np.ndarray, found: np.ndarray):
    """the centre (ln K at 50 %) and the slope of a logistic curve, by
    maximum a posteriori"""
    lnk, found = np.asarray(lnk, float), np.asarray(found, float)

    def neg(pars):
        centre, lslope = pars
        prob = 1.0 / (1.0 + np.exp(-(lnk - centre) / np.exp(lslope)))
        prob = np.clip(prob, 1e-9, 1 - 1e-9)
        like = np.sum(found * np.log(prob) + (1 - found) * np.log(1 - prob))
        prior = -0.5 * ((lslope - SLOPE_PRIOR[0]) / SLOPE_PRIOR[1]) ** 2
        return -(like + prior)
    start = np.array([np.median(lnk), SLOPE_PRIOR[0]])
    res = minimize(neg, start, method='Nelder-Mead',
                   options=dict(xatol=1e-4, fatol=1e-6, maxiter=4000))
    return float(res.x[0]), float(np.exp(res.x[1]))


def _levels(centre: float, slope: float):
    """K at 50 and 90 % of a logistic curve"""
    return float(np.exp(centre)), float(np.exp(centre + slope * np.log(9.0)))


def _next(injections: List[Dict[str, Any]], kmin: float, kmax: float):
    """the K of the next round: at 50 and 90 % of the fit, or beyond the
    injections when they are all found or all missed"""
    amps = np.array([inj['amp'] for inj in injections])
    found = np.array([inj['found'] for inj in injections])
    if found.all():
        low = 0.6 * amps.min()
        return [max(low, kmin)] * 2
    if not found.any():
        high = 1.7 * amps.max()
        return [min(high, kmax)] * 2
    k50, k90 = _levels(*_logistic_fit(np.log(amps), found))
    return [float(np.clip(k50, kmin, kmax)), float(np.clip(k90, kmin, kmax))]


def fip_map(series: RVData, pmin: float, pmax: Optional[float] = None,
            gp: Any = None, nband: int = 6, nround: int = 4,
            per_round: int = 6, kmax: int = 2, nsweep: int = 600,
            nburn: int = 300, workers: int = 4, threshold: float = THRESHOLD,
            seed: int = 1, nboot: int = 200) -> Dict[str, Any]:
    """
    The detection map of a series by the FIP (see the module)

    :param series: RVData, the series with the signals found taken out (the
                   nightly means when the analysis runs on them)
    :param pmin: float, the shortest period [days]
    :param pmax: float or None, the longest (the baseline when None)
    :param gp: the GP of the FIP (oafip(gp=)), as in the analysis
    :param nband: int, the period bands (log-spaced)
    :param nround: int, the rounds of injections per band
    :param per_round: int, the injections per band and round (half at the
                      50 %, half at the 90 % level of the current fit)
    :param kmax: int, the signals of each FIP
    :param nsweep: int, the recorded sweeps of each FIP (one chain)
    :param nburn: int, its burn-in
    :param workers: int, the FIPs at a time
    :param threshold: float, the FIP of a detection
    :param seed: int, the seed
    :param nboot: int, the bootstrap resamples of the levels' ranges

    :return: dict, kind ('fip'), period_edges, periods (band centres), K50,
             K90 (per band, m/s), K50_range and K90_range (16-84 %),
             injections (period, amp, fip, found, band, round), threshold,
             and the settings
    """
    rng = np.random.default_rng(seed)
    pmax = pmax or max(series.baseline, 2 * pmin)
    edges = np.geomspace(pmin, pmax, nband + 1)
    # the K range: from a fraction of the error of a visit to many of them
    med = float(np.median(series.err))
    ampmin, ampmax = 0.1 * med, 30 * med
    # the first round: around the K a sinusoid needs to stand out
    rms = float(np.std(series.rv))
    sig_k = rms * np.sqrt(2.0 / series.n)
    start = [3 * sig_k, 5 * sig_k, 8 * sig_k]
    injections: List[Dict[str, Any]] = []
    for iround in range(nround):
        jobs = []
        for band in range(nband):
            mine = [inj for inj in injections if inj['band'] == band]
            amps = (start if iround == 0 else _next(mine, ampmin, ampmax))
            for k in range(per_round):
                jobs.append(dict(
                    series=series, gp=gp, kmax=kmax, nsweep=nsweep,
                    nburn=nburn, threshold=threshold, band=band,
                    round=iround, amp=float(amps[k % len(amps)]),
                    period=float(np.exp(rng.uniform(np.log(edges[band]),
                                                    np.log(edges[band + 1])))),
                    phase=float(rng.uniform(0, 2 * np.pi)),
                    seed=int(rng.integers(1, 2 ** 31))))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            injections += list(pool.map(_one, jobs))
        nfound = sum(inj['found'] for inj in injections)
        log(f'detection map (FIP): round {iround + 1}/{nround}, '
            f'{len(injections)} injections, {nfound} found', 'value')
    out = summarize(injections, edges, nboot=nboot, seed=seed)
    out.update(threshold=threshold,
               settings=dict(nband=nband, nround=nround, per_round=per_round,
                             kmax=kmax, nsweep=nsweep, nburn=nburn,
                             gp=gp is not None))
    return out


def summarize(injections: List[Dict[str, Any]], edges: Sequence[float],
              nboot: int = 200, seed: int = 1) -> Dict[str, Any]:
    """
    The levels of each band from its injections: K found 50 and 90 % of the
    time and their bootstrap ranges. A level beyond the K injected is not
    extrapolated: it is NaN, with a note ('> K' when even the largest
    planet injected was not found often enough, '< K' when the smallest
    already was)

    :return: dict, kind, period_edges, periods, K50, K90, K50_range,
             K90_range, notes50, notes90, injections
    """
    rng = np.random.default_rng(seed)
    edges = np.asarray(edges, dtype=float)
    out = dict(kind='fip', period_edges=edges,
               periods=np.sqrt(edges[1:] * edges[:-1]), K50=[], K90=[],
               K50_range=[], K90_range=[], notes50=[], notes90=[],
               injections=injections)
    for band in range(len(edges) - 1):
        mine = [inj for inj in injections if inj['band'] == band]
        amps = np.array([inj['amp'] for inj in mine])
        lnk, found = np.log(amps), np.array([inj['found'] for inj in mine])
        lev = _levels(*_logistic_fit(lnk, found))
        boots = []
        for _ in range(nboot):
            pick = rng.integers(0, len(mine), len(mine))
            boots.append(_levels(*_logistic_fit(lnk[pick], found[pick])))
        boots = np.array(boots)
        for level, idx, key in ((lev[0], 0, '50'), (lev[1], 1, '90')):
            note = None
            if not found.any() or level > 1.2 * amps.max():
                note, level = f'> {amps.max():.1f}', np.nan
            elif level < amps.min() / 1.2:
                note, level = f'< {amps.min():.1f}', np.nan
            out[f'K{key}'].append(level)
            out[f'notes{key}'].append(note)
            out[f'K{key}_range'].append(
                (np.nan, np.nan) if note else
                tuple(np.percentile(boots[:, idx], [16, 84])))
    for key in ('K50', 'K90'):
        out[key] = np.array(out[key])
    return out


def figure(dmap: Dict[str, Any], marks: Sequence[Dict[str, Any]] = (),
           prot: Optional[float] = None, title: Optional[str] = None):
    """
    The injections (found filled, missed open), the K detected 50 % (dashed)
    and 90 % (solid) of the time per period band with their ranges, the
    signals found (stars) and the known planets (circles)

    :param dmap: dict, from fip_map()
    :param marks: list of dict, period, K and kind ('signal' or 'known')
    :param prot: float or None, the rotation period (dotted, with its half)

    :return: matplotlib figure
    """
    from koloa import plotting as kplot
    plt, C = kplot.plt, kplot.C
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    inj = dmap['injections']
    for found, style in ((True, dict(marker='o', mfc=C['koloa'], mec=C['koloa'],
                                     label='injected, found')),
                         (False, dict(marker='o', mfc='none', mec=C['muted'],
                                      label='injected, missed'))):
        sel = [i for i in inj if i['found'] == found]
        if sel:
            ax.plot([i['period'] for i in sel], [i['amp'] for i in sel],
                    ls='none', ms=4, alpha=0.8, **style)
    per = np.asarray(dmap['periods'])
    for key, rng_key, ls, lab in (('K50', 'K50_range', (0, (4, 2)), '50 % found'),
                                  ('K90', 'K90_range', '-', '90 % found')):
        vals = np.asarray(dmap[key], dtype=float)
        ax.plot(per, vals, color=C['text'], lw=1.6, ls=ls, label=lab)
        lo = np.array([r[0] for r in dmap[rng_key]], dtype=float)
        hi = np.array([r[1] for r in dmap[rng_key]], dtype=float)
        ax.fill_between(per, lo, hi, color=C['text'], alpha=0.08, lw=0)
    # the bands where 90 % is never reached: written, not extrapolated
    for ib, note in enumerate(dmap.get('notes90') or []):
        if note and note.startswith('>'):
            ax.text(per[ib], 0.97, f'not reached\n(K {note} m/s)', ha='center',
                    va='top', fontsize=6.5, color=C['outlier'],
                    transform=ax.get_xaxis_transform(),
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.85,
                              pad=1.5))
    for mark in marks:
        star = mark.get('kind') == 'signal'
        ax.plot(mark['period'], mark['K'], marker='*' if star else 'o',
                ms=12 if star else 8, mfc=C['outlier'] if star else 'none',
                mec='white' if star else C['text'], ls='none')
    if prot:
        for p in (prot, prot / 2):
            ax.axvline(p, color=C['muted'], lw=1.0, ls=(0, (1, 2)))
    for edge in dmap['period_edges']:
        ax.axvline(edge, color=C['muted'], lw=0.4, alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('period [d]')
    ax.set_ylabel('K [m s$^{-1}$]')
    ax.legend(loc='upper left', fontsize=7)
    if title:
        ax.set_title(title, loc='left')
    fig.tight_layout()
    return fig


# =============================================================================
# End of code
# =============================================================================

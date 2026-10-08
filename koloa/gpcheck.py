#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The signals of a series against a GP of the stellar activity

`gp_signals` fits the velocities with a GP of the activity (a rotation
kernel with the rotation period of the star as its prior when the archive
knows it, a SHO with a free period otherwise; chromatic with several
instruments), jointly with the signals: the GP alone, the GP with every
signal, and the GP with every signal but one, so that each signal has the
likelihood it adds over the activity. koloa's FIP has no GP, so the
significance with the GP comes from a periodogram whitened by the
covariance of the GP-only model (generalised least squares), its
false-alarm probability from simulations of that noise (the highest peak
of each); a second one, on what the GP and the signals leave, says what is
left.

The figures: the whole series, season by season, with the GP and the GP
with the signals; the whitened periodograms; each signal folded, the GP and
the other signals removed.

Created on 2026-09-30

@author: artigau
"""
import time
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from scipy.linalg import solve_triangular

from koloa import gp as kgp
from koloa import kepler
from koloa import plotting as kplot
from koloa.data import RVData
from koloa.diagnostics import seasons
from koloa.fit import Prior, RVModel
from koloa.log import log
from koloa.periodogram import find_peaks, frequency_grid
from koloa.utils import blas_threads

# =============================================================================
# Define variables
# =============================================================================
#: the smallest relative error of the rotation period in the prior (the
#: period a star shows drifts with the latitude of its spots)
MIN_ROT_ERR = 0.05
#: the range of a free GP period [days]
FREE_RANGE = (5.0, 200.0)
#: a gap that starts a new panel of the whole series [days]
SEASON_GAP = 60.0
#: a panel holds at least this many exposures (sparse seasons share one),
#: and there are at most MAX_PANELS
MIN_PANEL = 8
MAX_PANELS = 12
#: in a panel longer than this, the GP is drawn around the visits [days]
LONG_PANEL = 400.0
#: the most frequencies of a whitened periodogram (the oversampling is
#: lowered to stay below)
MAX_FREQ = 20000
#: the colours of the instruments, in turn
INST_KEYS = ('koloa', 'gaussian', 'soft', 'hard', 'white', 'visit')


# =============================================================================
# The GP and the fits
# =============================================================================
def gp_model(data: RVData, star: Optional[Dict[str, Any]] = None,
             reference: Optional[str] = None):
    """
    The GP of the activity: a rotation kernel with the rotation period of
    the star as a gaussian prior (its error at least MIN_ROT_ERR), or a SHO
    with a free period; chromatic with several instruments

    :param data: RVData, the series
    :param star: dict or None, the star (koloa.archive: rotation,
                 rotation_err)
    :param reference: str or None, the instrument whose amplitude the GP
                      has (the others relative to it; the first when None)

    :return: tuple, the GP of RVModel and its description
    """
    star = star or {}
    prot = star.get('rotation')
    if prot:
        err = max((star.get('rotation_err') or 0.0) / prot, MIN_ROT_ERR)
        gp = dict(kernel='rotation',
                  prior={'log_period': (float(np.log(prot)), float(err))},
                  init=dict(log_period=float(np.log(prot))))
        label = (f'rotation kernel (two SHOs, at P and P/2), P = {prot} d '
                 f'+- {100 * err:.0f} % (the rotation period of the NASA '
                 f'Exoplanet Archive) as a prior')
    else:
        top = min(FREE_RANGE[1], 0.5 * data.baseline)
        gp = dict(kernel='sho', prior={'log_period': Prior(
            'uniform', np.log(FREE_RANGE[0]), np.log(top))},
            init=dict(log_period=float(np.log(min(30.0, top)))))
        label = (f'SHO with a free period ({FREE_RANGE[0]:.0f} to {top:.0f} d;'
                 f' no rotation period known)')
    if len(data.instruments) > 1:
        ref = reference if reference in data.instruments \
            else data.instruments[0]
        gp.update(scale='instrument', reference=ref)
        label += (f'; chromatic: each instrument its own amplitude, relative '
                  f'to {ref}')
    return gp, label


def _fit(job):
    """one fit, in a worker"""
    start = time.time()
    with blas_threads(1):
        fit = RVModel(job['data'], job['planets'], gp=job['gp'],
                      **job['kwargs']).fit(nstart=job.get('nstart', 2),
                                           quiet=True)
    log(f'GP, {job["key"]}: ln post {fit.logpost:.2f} '
        f'({(time.time() - start) / 60:.1f} min)')
    return job['key'], fit


def loglike(fit) -> float:
    """the log likelihood of a fit (its log posterior without the prior)"""
    return float(fit.logpost - fit.model.log_prior(fit.theta))


def gp_values(fit) -> Dict[str, float]:
    """the GP parameters of a fit, in natural units"""
    out = {}
    for name, idx in fit.model.index.items():
        if name.startswith('gp_log_'):
            out[name[7:]] = float(np.exp(fit.theta[idx]))
        elif name.startswith('gp_logit_'):
            out[name[9:]] = float(1 / (1 + np.exp(-fit.theta[idx])))
    return out


# =============================================================================
# The periodogram whitened by the noise of a fit
# =============================================================================
def covariance(fit, keep: np.ndarray) -> np.ndarray:
    """the covariance of the noise of a fit (white and visit jitters, and
    its GP), on the points kept"""
    model, theta, data = fit.model, fit.theta, fit.model.data
    diag, seq_var = model.noise(theta)
    seq = np.asarray(data.seq)
    var = seq_var[seq] if np.ndim(seq_var) else np.full(data.n, seq_var)
    cov = np.diag(diag) + (seq[:, None] == seq[None, :]) * np.sqrt(
        var[:, None] * var[None, :])
    if model.gp is not None:
        kern = kgp.kernel_matrix(model.gp['kernel'],
                                 data.time[:, None] - data.time[None, :],
                                 model.gp_pars(theta))
        scale = model.gp_scale(theta)
        if scale is not None:
            kern = kern * scale[:, None] * scale[None, :]
        cov = cov + kern
    return cov[np.ix_(keep, keep)]


def whitened_periodogram(fit, freq: np.ndarray, nsim: int = 20000,
                         subtract: bool = False, seed: int = 1,
                         batch: int = 250) -> Dict[str, Any]:
    """
    Delta chi^2 of a sinusoid at each frequency, the noise that of a fit
    (generalised least squares; an offset per instrument and a linear trend
    fitted with it), and the highest peak of simulations of that noise

    :param fit: FitResult, with a GP (or without: white and visit jitters)
    :param freq: np.ndarray, the frequencies [1/day]
    :param nsim: int, the simulations
    :param subtract: bool, take the signals of the fit out first (what the
                     GP and the signals leave)
    :param seed: int, the seed
    :param batch: int, the simulations at a time

    :return: dict, freq, power, maxsim, kept (the points that are not
             outliers) and n
    """
    model, data = fit.model, fit.model.data
    keep = fit.outlier_prob < 0.5
    value = data.rv.copy()
    if subtract:
        for ip in range(len(model.planets)):
            value = value - model.planet_rv(fit.theta, ip)
    time_, value = data.time[keep], value[keep]
    inst = np.asarray(data.inst)[keep]
    lchol = np.linalg.cholesky(covariance(fit, keep))

    def white(arr):
        return solve_triangular(lchol, arr, lower=True)

    cols = [(inst == name).astype(float) for name in np.unique(inst)]
    cols.append((time_ - np.mean(time_)) / 365.25)
    base, _ = np.linalg.qr(white(np.array(cols).T))

    def project(arr):
        return arr - base @ (base.T @ arr)

    yw = project(white(value))
    bases = np.zeros((len(freq), 2, len(time_)))
    power = np.zeros(len(freq))
    for ifr, fr in enumerate(freq):
        arg = 2 * np.pi * fr * time_
        qq, _ = np.linalg.qr(project(white(np.array([np.cos(arg),
                                                     np.sin(arg)]).T)))
        bases[ifr] = qq.T
        power[ifr] = np.sum((qq.T @ yw) ** 2)
    rng = np.random.default_rng(seed)
    flat = bases.reshape(-1, len(time_))
    maxsim = []
    for start in range(0, nsim, batch):
        draws = project(rng.standard_normal((len(time_),
                                             min(batch, nsim - start))))
        proj = (flat @ draws) ** 2
        maxsim.append(np.max(proj.reshape(len(freq), 2, -1).sum(axis=1),
                             axis=0))
    return dict(freq=freq, power=power, maxsim=np.concatenate(maxsim),
                kept=int(np.sum(keep)), n=int(data.n))


def fap(res: Dict[str, Any], value: float) -> float:
    """the false-alarm probability of a peak: the fraction of simulations
    whose highest peak is higher (1/nsim when none is: a limit)"""
    return max(int(np.sum(res['maxsim'] >= value)), 1) / len(res['maxsim'])


def peaks(res: Dict[str, Any], npeaks: int = 5) -> List[Dict[str, float]]:
    """the highest peaks of a whitened periodogram, with their FAP"""
    out = []
    for ip in find_peaks(res['freq'], res['power'], npeaks):
        out.append(dict(period=float(1 / res['freq'][ip]),
                        dchi2=float(res['power'][ip]),
                        fap=fap(res, res['power'][ip])))
    return out


def at_period(res: Dict[str, Any], period: float,
              width: float = 0.01) -> Dict[str, float]:
    """the highest whitened power within a fraction of a period, and its
    FAP (as the highest peak of the whole range)"""
    per = 1 / res['freq']
    near = np.abs(per / period - 1) < width
    if not np.any(near):
        return dict(period=period, dchi2=np.nan, fap=np.nan)
    it = np.flatnonzero(near)[np.argmax(res['power'][near])]
    return dict(period=float(per[it]), dchi2=float(res['power'][it]),
                fap=fap(res, res['power'][it]))


# =============================================================================
# Everything
# =============================================================================
def gp_signals(data: RVData, periods: Sequence[float],
               star: Optional[Dict[str, Any]] = None,
               reference: Optional[str] = None,
               seq_jitter: Optional[str] = None, mstar: Optional[float] = None,
               workers: int = 4, nsim: int = 20000, pmin: float = 1.1,
               pmax: Optional[float] = None, seed: int = 1,
               ranges: Optional[Sequence[Sequence[float]]] = None
               ) -> Dict[str, Any]:
    """
    The signals of a series against a GP of the activity (see the module's
    docstring)

    :param data: RVData, the series
    :param periods: list of float, the periods of the signals [days]
    :param star: dict or None, the star (koloa.archive): its rotation period
                 sets the prior of the GP
    :param reference: str or None, the reference instrument of a chromatic
                      GP
    :param seq_jitter: str or None, as in RVModel
    :param mstar: float or None, the stellar mass [solar masses]
    :param workers: int, the fits at a time
    :param nsim: int, the simulations of each whitened periodogram
    :param pmin: float, the shortest period of the periodograms [days]
    :param pmax: float or None, the longest (the baseline when None)
    :param seed: int, the seed
    :param ranges: list or None, the range of each period [days] (2 % of it
                   when None): a period tested (a known planet's) stays
                   where it was tested

    :return: dict, gp (the RVModel GP), label, fits (null, all, without_<i>),
             orbits (per signal: the orbit with the GP and the ln L it adds),
             gp_null and gp_all (the GP parameters), whitened (null and
             residual periodograms), and runtime [s]
    """
    start = time.time()
    gp, label = gp_model(data, star, reference)
    kwargs = dict(likelihood='mixture', unit='both', trend=1,
                  seq_jitter=seq_jitter)
    ranges = ranges or [(0.98 * per, 1.02 * per) for per in periods]
    planets = [dict(period=float(np.clip(per, low, high)),
                    period_range=(float(low), float(high)))
               for per, (low, high) in zip(periods, ranges)]
    jobs = [dict(key='all', planets=planets)] if planets else []
    jobs.append(dict(key='null', planets=[]))
    if len(planets) > 1:
        jobs += [dict(key=f'without_{ip}',
                      planets=planets[:ip] + planets[ip + 1:])
                 for ip in range(len(planets))]
    for job in jobs:
        job.update(data=data, gp=gp, kwargs=kwargs)
    nproc = max(1, min(workers, len(jobs)))
    log(f'GP: {label}; {len(jobs)} fits on {nproc} processes')
    with ProcessPoolExecutor(nproc) as pool:
        fits = dict(pool.map(_fit, jobs))
    orbits = []
    if planets:
        full = fits['all']
        for ip, orb in enumerate(full.orbits(mstar=mstar)):
            other = fits['null'] if len(planets) == 1 else \
                fits[f'without_{ip}']
            entry = {key: [float(val) for val in orb[key]]
                     for key in ('P', 'K', 'e', 'msini', 'tc')}
            entry['dlnl'] = loglike(full) - loglike(other)
            orbits.append(entry)
    top = pmax or data.baseline
    nfreq = (1 / pmin - 1 / top) * data.baseline
    oversample = int(np.clip(MAX_FREQ / max(nfreq, 1.0), 2, 10))
    freq = frequency_grid(data.time, pmin, top, oversample)
    log(f'GP: whitened periodograms, {len(freq)} frequencies, {nsim} '
        f'simulations each')
    whitened = dict(null=whitened_periodogram(fits['null'], freq, nsim,
                                              seed=seed))
    if planets:
        whitened['residual'] = whitened_periodogram(
            fits['all'], freq, nsim, subtract=True, seed=seed + 1)
    for key, res in whitened.items():
        res['peaks'] = peaks(res)
    for orb in orbits:
        orb['whitened'] = at_period(whitened['null'], orb['P'][0])
    out = dict(gp=gp, label=label, fits=fits, orbits=orbits,
               gp_null=gp_values(fits['null']),
               gp_all=gp_values(fits['all']) if planets else None,
               whitened=whitened, runtime=time.time() - start)
    log(f'GP: done in {out["runtime"] / 60:.1f} min; ' + '; '.join(
        f'{orb["P"][0]:.4f} d: K = {orb["K"][0]:.2f} m/s, ln L '
        f'{orb["dlnl"]:+.1f}' for orb in orbits), 'value')
    return out


# =============================================================================
# Figures
# =============================================================================
def _colours(insts: Sequence[str]) -> Dict[str, str]:
    """a colour per instrument"""
    return {inst: kplot.C[INST_KEYS[it % len(INST_KEYS)]]
            for it, inst in enumerate(insts)}


def _planets_at(fit, times: np.ndarray) -> np.ndarray:
    """every orbit of a fit, at any time"""
    model = fit.model
    total = np.zeros(len(times))
    for ip in range(len(model.planets)):
        period, tperi, ecc, omega, amp = model.orbit(fit.theta, ip)
        total += kepler.rv_keplerian(times, period, tperi, ecc, omega, amp)
    return total


def _panels(time: np.ndarray, min_points: int = MIN_PANEL,
            max_panels: int = MAX_PANELS) -> np.ndarray:
    """the panel of each point: a season of min_points or more has its own,
    the seasons of fewer points in between share one (min_points of them
    at least; fewer than max_panels in all)"""
    season = seasons(time, SEASON_GAP)
    groups = [np.flatnonzero(season == sid) for sid in range(season.max() + 1)]
    while True:
        merged: List[np.ndarray] = []
        pooled: List[bool] = []
        pool: List[np.ndarray] = []
        for group in groups:
            if len(group) >= min_points:
                if pool:
                    # a few sparse points before a full season: with the
                    #   sparse panel before them, if there is one
                    if sum(len(item) for item in pool) < min_points \
                            and pooled and pooled[-1]:
                        merged[-1] = np.concatenate([merged[-1]] + pool)
                    else:
                        merged.append(np.concatenate(pool))
                        pooled.append(True)
                    pool = []
                merged.append(group)
                pooled.append(False)
                continue
            pool.append(group)
            if sum(len(item) for item in pool) >= min_points:
                merged.append(np.concatenate(pool))
                pooled.append(True)
                pool = []
        if pool:
            # the last few sparse points: with the sparse panel before, if
            #   it is one
            if pooled and pooled[-1]:
                merged[-1] = np.concatenate([merged[-1]] + pool)
            else:
                merged.append(np.concatenate(pool))
        if len(merged) <= max_panels:
            break
        min_points *= 2
    out = np.zeros(len(time), dtype=int)
    for ipanel, group in enumerate(merged):
        out[group] = ipanel
    return out


def _segments(times: np.ndarray) -> List[np.ndarray]:
    """where to draw the GP in a panel: all of it when it is short, around
    each visit when it spans years (the GP would fill it)"""
    t0, t1 = times.min(), times.max()
    if t1 - t0 <= LONG_PANEL:
        pad = max(2.0, 0.03 * (t1 - t0))
        return [np.linspace(t0 - pad, t1 + pad, 800)]
    half = max(2.0, 0.005 * (t1 - t0))
    out: List[List[float]] = []
    for tt in np.unique(np.round(np.sort(times), 1)):
        if out and tt - half <= out[-1][1]:
            out[-1][1] = tt + half
        else:
            out.append([tt - half, tt + half])
    return [np.linspace(low, high, 60) for low, high in out]


def figure_sequence(fit, title: str):
    """
    The whole series, season by season: every exposure (the offsets and the
    trend removed; outliers as crosses, those beyond the axes counted), the
    GP and its one sigma, and the GP with the signals; the seasons of few
    exposures share a panel, and in a panel that spans years the GP is
    drawn around the visits only

    :param fit: FitResult, with a GP (and signals)
    :param title: str, the title

    :return: matplotlib figure
    """
    from koloa.outliers import _date
    model, theta, data = fit.model, fit.theta, fit.model.data
    clean = data.rv - model.systematics(theta)
    diag, _ = model.noise(theta)
    err = np.sqrt(diag)
    prob = fit.outlier_prob
    season = _panels(data.time)
    nseason = int(season.max()) + 1
    ncol = 1 if nseason <= 3 else 2
    nrow = int(np.ceil(nseason / ncol))
    fig, axes = kplot.plt.subplots(nrow, ncol, squeeze=False,
                                   figsize=(11, (3.4 if nrow <= 2 else 2.4)
                                            * nrow + 0.8))
    colours = _colours(data.instruments)
    chromatic = bool(getattr(model, 'gp_scale_insts', []))
    good_all = prob < 0.5
    for iseason in range(nseason):
        ax = axes.flat[iseason]
        sel = season == iseason
        t0, t1 = data.time[sel].min(), data.time[sel].max()
        lows, highs = [], []
        for inst in data.instruments:
            isel = sel & (data.inst == inst)
            if not np.any(isel):
                continue
            col = colours[inst]
            for iseg, grid in enumerate(_segments(data.time[isel])):
                mean, sig = fit.gp_prediction(
                    time=grid, inst=inst if chromatic else None)
                first = iseason == 0 and iseg == 0
                ax.fill_between(grid, mean - sig, mean + sig, color=col,
                                alpha=0.18, lw=0)
                ax.plot(grid, mean, color=col, lw=1.3,
                        label=f'GP ({inst})' if first else None)
                if model.planets:
                    # faint, under the exposures: it must not hide them
                    ax.plot(grid, mean + _planets_at(fit, grid),
                            color=kplot.C['text'], lw=0.6, alpha=0.4,
                            zorder=2, label='GP + signals' if first
                            and inst == data.instruments[0] else None)
                lows.append(np.min(mean - sig))
                highs.append(np.max(mean + sig))
            good = isel & good_all
            # each exposure seen through its neighbours
            ax.errorbar(data.time[good], clean[good], err[good], fmt='none',
                        ecolor=col, elinewidth=0.7, alpha=kplot.ERR_ALPHA,
                        zorder=4)
            ax.plot(data.time[good], clean[good], ls='none', marker='o',
                    ms=3.6, mfc=col, mec='none', alpha=kplot.POINT_ALPHA,
                    zorder=5, label=inst if iseason == 0 else None)
            if np.any(good):
                lows.append(np.min(clean[good] - err[good]))
                highs.append(np.max(clean[good] + err[good]))
        low, high = min(lows), max(highs)
        margin = 0.08 * (high - low)
        low, high = low - margin, high + margin
        bad = sel & ~good_all
        inside = bad & (clean > low) & (clean < high)
        ax.plot(data.time[inside], clean[inside], 'x',
                color=kplot.C['outlier'], ms=4,
                label='outliers' if iseason == 0 else None)
        off = int(np.sum(bad & ~inside))
        if off:
            ax.text(0.99, 0.03, f'{off} outlier(s) off the axes',
                    transform=ax.transAxes, ha='right', fontsize=7,
                    color=kplot.C['outlier'])
        ax.set_ylim(low, high)
        pad = max(2.0, 0.03 * (t1 - t0))
        ax.set_xlim(t0 - pad, t1 + pad)
        ax.set_title(f'{_date(t0)} to {_date(t1)} ({int(np.sum(sel))} '
                     f'exposures)' + (', the GP around the visits'
                                      if t1 - t0 > LONG_PANEL else ''),
                     fontsize=8, loc='left')
        ax.set_ylabel('RV [m/s]', fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes.flat[nseason:]:
        ax.axis('off')
    for ax in axes[-1]:
        ax.set_xlabel('BJD - 2400000', fontsize=8)
    kplot.solid_legend(axes.flat[0].legend(fontsize=6.5, ncol=4,
                                           loc='upper left'))
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    return fig


def figure_whitened(whitened: Dict[str, Any], marks: Sequence, title: str):
    """
    The whitened periodograms, their false-alarm levels, and the periods
    that matter marked

    :param whitened: dict, from gp_signals
    :param marks: list of (period, label, colour key)
    :param title: str
    """
    labels = dict(null='the GP alone as the noise',
                  residual='what the GP and the signals leave')
    fig, axes = kplot.plt.subplots(len(whitened), 1, sharex=True,
                                   figsize=(10, 2.8 * len(whitened) + 0.6))
    axes = np.atleast_1d(axes)
    for ax, (key, res) in zip(axes, whitened.items()):
        ax.plot(1 / res['freq'], res['power'], color=kplot.C['text'], lw=0.6)
        for lev, style in ((0.01, ':'), (0.001, '--')):
            if len(res['maxsim']) * lev >= 5:
                ax.axhline(np.quantile(res['maxsim'], 1 - lev), ls=style,
                           color=kplot.C['threshold'], lw=0.8,
                           label=f'FAP {lev:g}')
        for per, label, ckey in marks:
            ax.axvline(per, color=kplot.C[ckey], ls='--', lw=0.8)
            if ax is axes[0]:
                ax.text(per, 1.0, f' {label}', color=kplot.C[ckey],
                        fontsize=7, rotation=90, va='top', ha='right',
                        transform=ax.get_xaxis_transform())
        ax.set_xscale('log')
        ax.set_ylabel(r'$\Delta\chi^2$')
        ax.set_title(labels.get(key, key), fontsize=9, loc='left')
        ax.legend(fontsize=7, loc='upper right')
    axes[-1].set_xlabel('period [d]')
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return fig


def _visit_means(time_, value, err, seq):
    """the weighted mean of each visit"""
    out = []
    for sid in np.unique(seq):
        sel = seq == sid
        wgt = 1 / err[sel] ** 2
        out.append((np.sum(time_[sel] * wgt) / np.sum(wgt),
                    np.sum(value[sel] * wgt) / np.sum(wgt),
                    1 / np.sqrt(np.sum(wgt))))
    return np.array(out).T if out else np.zeros((3, 0))


def figure_fold(fit, ip: int, title: str,
                published: Optional[Dict[str, Any]] = None):
    """
    A signal folded, the GP and the other signals removed: the visits of
    each instrument, the orbit, and the published one when there is one
    (its K, at its conjunction carried to the fit)

    :param fit: FitResult, with the GP and the signals
    :param ip: int, the signal
    :param title: str
    :param published: dict or None, K and dphase (the published conjunction
                      minus the fitted one, in phase)
    """
    model, theta, data = fit.model, fit.theta, fit.model.data
    period, tperi, ecc, omega, amp = model.orbit(theta, ip)
    tconj = kepler.tp_to_tc(tperi, period, ecc, omega)
    resid = fit.residuals(ip)
    diag, _ = model.noise(theta)
    keep = fit.outlier_prob < 0.5
    colours = _colours(data.instruments)
    fig, ax = kplot.plt.subplots(figsize=(7.5, 3.8))
    phi = np.linspace(0, 1, 500)
    for inst in data.instruments:
        sel = keep & (data.inst == inst)
        if not np.any(sel):
            continue
        tt, vv, ee = _visit_means(data.time[sel], resid[sel],
                                  np.sqrt(diag[sel]),
                                  np.asarray(data.seq)[sel])
        phase = ((tt - tconj) / period) % 1
        # each visit seen through its neighbours
        ax.errorbar(phase, vv, ee, fmt='none', ecolor=colours[inst],
                    elinewidth=0.6, alpha=kplot.ERR_ALPHA)
        ax.plot(phase, vv, ls='none', marker='o', ms=3.4, mfc=colours[inst],
                mec='none', alpha=kplot.POINT_ALPHA,
                label=f'{inst} ({len(tt)} visits)')
    curve = kepler.rv_keplerian(tconj + phi * period, period, tperi, ecc,
                                omega, amp)
    ax.plot(phi, curve, color=kplot.C['text'], lw=1.5,
            label=f'fit: K = {amp:.2f} m/s')
    if published:
        ax.plot(phi, -published['K'] * np.sin(2 * np.pi * (
            phi - published.get('dphase', 0.0))), color=kplot.C['outlier'],
            lw=1.2, ls='--', label=f'published: K = {published["K"]} m/s'
            + (' at its ephemeris' if 'dphase' in published else ''))
    ax.set_xlabel(f'phase (P = {period:.4f} d, 0 = conjunction)')
    ax.set_ylabel('RV [m/s]')
    ax.set_title(title, fontsize=9, loc='left')
    # where the fold of an orbit leaves room: phase 0 is its conjunction,
    #   the velocity falling
    kplot.solid_legend(ax.legend(fontsize=7, loc='upper left'))
    fig.tight_layout()
    return fig


# =============================================================================
# End of code
# =============================================================================

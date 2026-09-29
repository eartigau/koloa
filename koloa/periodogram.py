#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Periodograms that know about outliers.

- `gls`: the generalised Lomb-Scargle periodogram (Zechmeister & Kurster
  2009), the reference everything else is compared to.
- `profile_periodogram`: at each frequency, an offset and a sinusoid fitted
  by maximum likelihood with a free jitter, and optionally with the
  outlier mixture of koloa.noise. The power is the gain in log likelihood
  over the offset alone, so the gaussian and the outlier-aware versions sit
  on the same axis and can be read against each other. With outliers on,
  this is the outlier-aware periodogram (OAP).
- `jackknife`: the periodogram with each sequence (or point) left out in
  turn. The envelope shows which peaks depend on a handful of units, and
  the influence of each unit on a peak says which ones.
- `window` and `aliases`: what the sampling alone puts at a period, and
  the periods a peak can be confused with.

Created on 2026-09-27

@author: artigau
"""
from typing import Dict, List, Optional, Sequence

import numpy as np

from koloa.data import RVData, robust_std

# =============================================================================
# Define variables
# =============================================================================
#: the sidereal day, the year and the synodic month [days]
SIDEREAL_DAY = 0.99726957
YEAR = 365.25
SYNODIC_MONTH = 29.530589

#: sampling frequencies whose aliases are checked [1/day]
ALIAS_FREQUENCIES = {'1 day': 1.0 / SIDEREAL_DAY, '1 year': 1.0 / YEAR,
                     '1 month': 1.0 / SYNODIC_MONTH}


# =============================================================================
# Frequency grids
# =============================================================================
def frequency_grid(time: np.ndarray, pmin: float = 1.1,
                   pmax: Optional[float] = None, oversample: int = 10
                   ) -> np.ndarray:
    """
    A grid regular in frequency, `oversample` points per peak width (1/T)

    :param time: np.ndarray, the time [days]
    :param pmin: float, the shortest period [days]
    :param pmax: float or None, the longest period (twice the baseline when
                 None) [days]
    :param oversample: int, points per 1/T

    :return: np.ndarray, the frequencies [1/day], increasing
    """
    baseline = float(np.ptp(time))
    if pmax is None:
        pmax = 2 * baseline
    step = 1.0 / (oversample * baseline)
    fmin, fmax = 1.0 / pmax, 1.0 / pmin
    return np.arange(fmin, fmax + 0.5 * step, step)


def _trig(time: np.ndarray, freq: np.ndarray, tref: float = 0.0):
    """
    cos and sin of 2 pi f (t - tref), (n x G)

    :param time: np.ndarray, the time [days]
    :param freq: np.ndarray, the frequencies [1/day]
    :param tref: float, the reference time [days]

    :return: tuple, the cos and sin matrices
    """
    phase = 2 * np.pi * np.outer(time - tref, freq)
    return np.cos(phase), np.sin(phase)


# =============================================================================
# The generalised Lomb-Scargle periodogram
# =============================================================================
def _gls_sums(weight: np.ndarray, value: np.ndarray, cos: np.ndarray,
              sin: np.ndarray) -> Dict[str, np.ndarray]:
    """
    The weighted sums the GLS power is made of (additive over points)

    :return: dict, the sums
    """
    return dict(w=np.sum(weight), y=weight @ value, yy=weight @ value ** 2,
                c=weight @ cos, s=weight @ sin, yc=(weight * value) @ cos,
                ys=(weight * value) @ sin, cc=weight @ cos ** 2,
                ss=weight @ sin ** 2, cs=weight @ (cos * sin))


def _gls_power(sums: Dict[str, np.ndarray]) -> np.ndarray:
    """
    The GLS power from its sums (Zechmeister & Kurster 2009, eq. 5)

    :param sums: dict, the output of _gls_sums (or a difference of them)

    :return: np.ndarray, the power (0 to 1)
    """
    wsum = sums['w']
    ymean, cmean, smean = sums['y'] / wsum, sums['c'] / wsum, sums['s'] / wsum
    yy = sums['yy'] / wsum - ymean ** 2
    yc = sums['yc'] / wsum - ymean * cmean
    ys = sums['ys'] / wsum - ymean * smean
    cc = sums['cc'] / wsum - cmean ** 2
    ss = sums['ss'] / wsum - smean ** 2
    cs = sums['cs'] / wsum - cmean * smean
    det = cc * ss - cs ** 2
    with np.errstate(invalid='ignore', divide='ignore'):
        power = (ss * yc ** 2 + cc * ys ** 2 - 2 * cs * yc * ys) / (yy * det)
    return np.clip(np.nan_to_num(power), 0, 1)


def gls(time: np.ndarray, value: np.ndarray, err: np.ndarray,
        freq: np.ndarray) -> np.ndarray:
    """
    The generalised Lomb-Scargle periodogram (floating mean, weighted)

    :param time: np.ndarray, the time [days]
    :param value: np.ndarray, the values
    :param err: np.ndarray, the error bars
    :param freq: np.ndarray, the frequencies [1/day]

    :return: np.ndarray, the power (0 to 1)
    """
    cos, sin = _trig(time, freq, float(np.mean(time)))
    return _gls_power(_gls_sums(1.0 / err ** 2, value, cos, sin))


# =============================================================================
# The profile likelihood periodogram, gaussian or outlier-aware
# =============================================================================
def profile_periodogram(time: np.ndarray, value: np.ndarray, err: np.ndarray,
                        freq: np.ndarray, outliers: bool = True,
                        width: Optional[float] = None,
                        frac_prior: Sequence[float] = (1.0, 20.0),
                        niter: int = 30, chunk: int = 4000, trend: int = 0
                        ) -> Dict[str, np.ndarray]:
    """
    Maximum likelihood of offset + sinusoid at each frequency, with a free
    jitter and (optionally) the outlier mixture

    At every frequency the model value = c + A cos + B sin is fitted by EM:
    the E step gives each point its probability of being an outlier, the M
    step solves the weighted normal equations, the outlier fraction and the
    jitter (a Fisher scoring step). The outlier width W is common to all
    frequencies. The null model (offset only) is fitted the same way, and
    the power is the gain in log likelihood, Delta ln L.

    :param time: np.ndarray, the time [days]
    :param value: np.ndarray, the values
    :param err: np.ndarray, the error bars
    :param freq: np.ndarray, the frequencies [1/day]
    :param outliers: bool, the mixture (True) or a plain gaussian (False)
    :param width: float or None, the outlier width W (8 robust sigma of the
                  data when None)
    :param frac_prior: tuple, the beta prior (a, b) of the outlier fraction,
                       used as a MAP penalty
    :param niter: int, the number of EM iterations
    :param chunk: int, how many frequencies at a time (memory)
    :param trend: int, the degree of a polynomial fitted in both the null
                  and the sinusoid model (a drift or a companion's pull
                  otherwise makes every long period look significant)

    :return: dict, dlnl (G), lnl (G), lnl0, frac (G), jitter (G), amp (G),
             and outprob0 (n, the outlier probability of each point under
             the null model)
    """
    value = np.asarray(value, dtype=float)
    err = np.asarray(err, dtype=float)
    scale = max(robust_std(value), float(np.median(err)))
    if width is None:
        width = 8.0 * scale
    tref = float(np.mean(time))
    span = max(float(np.max(time) - np.min(time)), 1e-9)
    poly = [((time - tref) / span) ** deg for deg in range(1, trend + 1)]
    # the null model, the same code with no sinusoid
    null = _profile_em(np.stack([np.ones(len(time))] + poly,
                                axis=1)[:, :, None], value, err, outliers,
                       width, frac_prior, niter)
    out = dict(dlnl=np.zeros(len(freq)), lnl=np.zeros(len(freq)),
               frac=np.zeros(len(freq)), jitter=np.zeros(len(freq)),
               amp=np.zeros(len(freq)), lnl0=float(null['lnl'][0]),
               outprob0=null['resp'][:, 0], width=width)
    for start in range(0, len(freq), chunk):
        sub = freq[start:start + chunk]
        cos, sin = _trig(time, sub, tref)
        columns = [np.ones_like(cos)] + [np.repeat(col[:, None], len(sub),
                                                   axis=1) for col in poly]
        design = np.stack(columns + [cos, sin], axis=1)
        fit = _profile_em(design, value, err, outliers, width, frac_prior,
                          niter)
        sl = slice(start, start + len(sub))
        out['lnl'][sl] = fit['lnl']
        out['frac'][sl] = fit['frac']
        out['jitter'][sl] = fit['jitter']
        out['amp'][sl] = np.hypot(fit['beta'][1 + trend],
                                  fit['beta'][2 + trend])
    out['dlnl'] = out['lnl'] - out['lnl0']
    return out


def _profile_em(design: np.ndarray, value: np.ndarray, err: np.ndarray,
                outliers: bool, width: float, frac_prior: Sequence[float],
                niter: int) -> Dict[str, np.ndarray]:
    """
    EM fit of a linear model with jitter and outliers, many models at once

    :param design: np.ndarray, (n x p x G) one design matrix per model
    :param value: np.ndarray, (n) the values
    :param err: np.ndarray, (n) the error bars
    :param outliers: bool, whether the outlier component is on
    :param width: float, the outlier width
    :param frac_prior: tuple, the beta prior of the fraction
    :param niter: int, the number of iterations

    :return: dict, beta (p x G), lnl (G), frac (G), jitter (G), resp (n x G)
    """
    npts, npar, nmod = design.shape
    err2 = err[:, None] ** 2
    # start: weighted least squares with the robust excess as a jitter
    excess = max(robust_std(value) ** 2 - np.median(err) ** 2, 0.0)
    jit2 = np.full(nmod, excess)
    frac = np.full(nmod, 0.02 if outliers else 0.0)
    resp = np.zeros((npts, nmod))
    abeta, bprior = frac_prior
    beta = np.zeros((npar, nmod))
    for it in range(niter):
        vin = err2 + jit2[None, :]
        vout = vin + width ** 2
        weight = (1 - resp) / vin + resp / vout
        # the normal equations, one small system per model
        amat = np.einsum('ipg,iqg,ig->gpq', design, design, weight)
        bvec = np.einsum('ipg,i,ig->gp', design, value, weight)
        amat += 1e-12 * np.eye(npar)[None, :, :]
        beta = np.linalg.solve(amat, bvec[:, :, None])[:, :, 0].T
        resid = value[:, None] - np.einsum('ipg,pg->ig', design, beta)
        r2 = resid ** 2
        lin = -0.5 * (r2 / vin + np.log(vin))
        if outliers:
            lout = -0.5 * (r2 / vout + np.log(vout))
            lgood = np.log1p(-frac)[None, :] + lin
            lbad = np.log(frac)[None, :] + lout
            total = np.logaddexp(lgood, lbad)
            resp = np.exp(lbad - total)
            # the MAP of the fraction under its beta prior
            frac = ((np.sum(resp, axis=0) + abeta - 1)
                    / (npts + abeta + bprior - 2))
            frac = np.clip(frac, 1e-6, 0.5)
        # Fisher scoring for the jitter, the outliers weighted out of it
        good = 1 - resp
        score = (np.sum(good * (r2 - vin) / vin ** 2, axis=0)
                 + np.sum(resp * (r2 - vout) / vout ** 2, axis=0))
        info = (np.sum(good / vin ** 2, axis=0)
                + np.sum(resp / vout ** 2, axis=0))
        jit2 = np.clip(jit2 + score / info, 0.0, None)
    vin = err2 + jit2[None, :]
    lin = -0.5 * (r2 / vin + np.log(vin) + np.log(2 * np.pi))
    if outliers:
        vout = vin + width ** 2
        lout = -0.5 * (r2 / vout + np.log(vout) + np.log(2 * np.pi))
        lnl = np.sum(np.logaddexp(np.log1p(-frac)[None, :] + lin,
                                  np.log(frac)[None, :] + lout), axis=0)
    else:
        lnl = np.sum(lin, axis=0)
    return dict(beta=beta, lnl=lnl, frac=frac, jitter=np.sqrt(jit2),
                resp=resp)


def oap(data: RVData, freq: np.ndarray, unit: str = 'sequence',
        outliers: bool = True, **kwargs) -> Dict[str, np.ndarray]:
    """
    The outlier-aware periodogram of a series

    With unit = 'sequence' each visit is first reduced to its weighted mean
    (a visit lasts minutes, a period days), so an outlier is a visit, and
    the exposures of one night do not count as independent evidence.

    :param data: RVData, the series
    :param freq: np.ndarray, the frequencies [1/day]
    :param unit: str, sequence or point
    :param outliers: bool, the mixture (True) or the gaussian (False)
    :param kwargs: passed to profile_periodogram

    :return: dict, see profile_periodogram
    """
    src = data.binned() if unit == 'sequence' else data
    out = profile_periodogram(src.time, src.rv, src.err, freq,
                              outliers=outliers, **kwargs)
    out['unit'] = unit
    return out


# =============================================================================
# Jackknife: which units hold a peak up
# =============================================================================
def jackknife(data: RVData, freq: np.ndarray, unit: str = 'sequence'
              ) -> Dict[str, np.ndarray]:
    """
    The GLS with each unit left out in turn

    The GLS power is a ratio of weighted sums over points, so leaving a unit
    out is a subtraction, and the whole jackknife costs about one
    periodogram.

    :param data: RVData, the series
    :param freq: np.ndarray, the frequencies [1/day]
    :param unit: str, sequence or point

    :return: dict, power (G, all units), low and high (G, the envelope),
             loo (nunit x G, each leave-one-out power), unit (n, the unit
             of each point)
    """
    units = data.seq if unit == 'sequence' else np.arange(data.n)
    nunit = int(np.max(units)) + 1
    cos, sin = _trig(data.time, freq, float(np.mean(data.time)))
    weight = 1.0 / data.err ** 2
    total = _gls_sums(weight, data.rv, cos, sin)
    power = _gls_power(total)
    loo = np.zeros((nunit, len(freq)))
    for uu in range(nunit):
        mask = units == uu
        part = _gls_sums(weight[mask], data.rv[mask], cos[mask], sin[mask])
        loo[uu] = _gls_power({key: total[key] - part[key] for key in total})
    return dict(power=power, low=np.min(loo, axis=0),
                high=np.max(loo, axis=0), loo=loo, unit=units)


def influence(jack: Dict[str, np.ndarray], freq: np.ndarray, period: float
              ) -> np.ndarray:
    """
    How much each unit holds up the peak at a period

    :param jack: dict, the output of jackknife
    :param freq: np.ndarray, the frequencies [1/day]
    :param period: float, the period of the peak [days]

    :return: np.ndarray, (nunit) power with the unit minus power without it
             (positive: the unit supports the peak)
    """
    idx = int(np.argmin(np.abs(freq - 1.0 / period)))
    return jack['power'][idx] - jack['loo'][:, idx]


# =============================================================================
# Window, aliases, peaks
# =============================================================================
def window(time: np.ndarray, freq: np.ndarray) -> np.ndarray:
    """
    The spectral window of the sampling, |sum exp(2 pi i f t)|^2 / n^2

    :param time: np.ndarray, the time [days]
    :param freq: np.ndarray, the frequencies [1/day]

    :return: np.ndarray, the window power (1 at zero frequency)
    """
    cos, sin = _trig(time, np.atleast_1d(freq), float(np.mean(time)))
    return (np.sum(cos, axis=0) ** 2 + np.sum(sin, axis=0) ** 2) / len(time) ** 2


def aliases(freq0: float, fmin: float = 0.0, fmax: float = np.inf,
            names: Optional[Sequence[str]] = None) -> List[Dict[str, float]]:
    """
    The frequencies a peak can be the alias of, or be aliased to

    For a sampling frequency fs, a signal at f also shows at |f - fs| and
    f + fs (Dawson & Fabrycky 2010).

    :param freq0: float, the frequency of the peak [1/day]
    :param fmin: float, the lowest frequency to keep [1/day]
    :param fmax: float, the highest frequency to keep [1/day]
    :param names: list of str or None, which sampling frequencies (all when
                  None)

    :return: list of dict, name, freq and period of each alias
    """
    out = []
    for name, fsamp in ALIAS_FREQUENCIES.items():
        if names is not None and name not in names:
            continue
        for sign in (-1, 1):
            falias = abs(freq0 + sign * fsamp)
            if fmin <= falias <= fmax and falias > 0:
                out.append(dict(name=name, freq=falias, period=1.0 / falias))
    return out


def find_peaks(freq: np.ndarray, score: np.ndarray, npeaks: int = 5,
               separation: Optional[float] = None) -> List[int]:
    """
    The highest peaks of a periodogram, at least one peak width apart

    :param freq: np.ndarray, the frequencies [1/day]
    :param score: np.ndarray, higher is better
    :param npeaks: int, how many
    :param separation: float or None, the smallest separation [1/day]
                       (two grid steps when None)

    :return: list of int, the indices of the peaks, best first
    """
    if separation is None:
        separation = 2 * float(np.median(np.diff(freq)))
    score = np.array(score, dtype=float)
    out = []
    for _ in range(npeaks):
        if not np.any(np.isfinite(score)):
            break
        idx = int(np.nanargmax(score))
        if not np.isfinite(score[idx]):
            break
        out.append(idx)
        score[np.abs(freq - freq[idx]) < separation] = -np.inf
    return out


# =============================================================================
# End of code
# =============================================================================

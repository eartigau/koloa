#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Velocity series with known planets, known activity and known outliers.

The outliers are the ones a real survey produces, not gaussian points
drawn a little wide:

- 'spike': one exposure off by several sigma (a cosmic ray on a key line,
  a bad pixel in the reference, a wrong order merge);
- 'visit': every exposure of a visit moved together (a drift of the
  calibration, moonlight on the fibre, a thin cloud with a coloured
  transmission): the outlier is the sequence, not the point;
- 'flare': a positive excursion that decays over the next visits, the
  shape of a stellar flare seen through the line profiles;
- 'heavy': heavy-tailed noise everywhere (a Student-t of few degrees of
  freedom), which a mixture of two gaussians only approximates.

Clear outliers (six or eight median errors) are easy for any method. Real
series also hold borderline ones, a visit or an exposure two to four median
errors off, which a 3-sigma clip may or may not catch and which a mixture
can only half-flag. An amplitude given as a range (low, high) draws the
size uniformly in it: BORDERLINE holds the recipes the demos use.

The observing times can be drawn (seasons, visits at a similar hour, a few
exposures per visit) or taken from a real series, whose error bars are then
kept as well: that is the most honest injection.

Created on 2026-09-27

@author: artigau
"""
from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np

from koloa import gp as kgp
from koloa import kepler
from koloa.data import RVData

# =============================================================================
# Define variables
# =============================================================================
#: clear outliers: whole visits and single exposures well away from the bulk
#: (6 and 8 median errors, times 1 + E with E exponential of mean 0.5)
CLEAR = [dict(kind='visit', frac=0.05, amplitude=6.0, label='clear'),
         dict(kind='spike', frac=0.02, amplitude=8.0, label='clear')]
#: borderline outliers: visits off by 1.5 to 3 median errors (2 to 3.5
#: sigma of a visit mean with three exposures and a 2 m/s visit jitter) and
#: exposures off by 2.5 to 4 median errors (2 to 3.5 sigma of one exposure)
BORDERLINE = [dict(kind='visit', frac=0.05, amplitude=(1.5, 3.0),
                   label='borderline'),
              dict(kind='spike', frac=0.02, amplitude=(2.5, 4.0),
                   label='borderline')]
#: what the demos inject: both
REALISTIC = CLEAR + BORDERLINE


# =============================================================================
# Define functions
# =============================================================================
def observing_times(nvisits: int = 60, per_visit: int = 3,
                    baseline: float = 1000.0, season: float = 240.0,
                    cadence: float = 0.0104, hour_scatter: float = 0.08,
                    seed: int = 1) -> np.ndarray:
    """
    Times of a realistic campaign: seasons, one visit a night, a few
    exposures back to back

    Visits happen at about the same hour every night, which is what makes
    the 1-day aliases of real data, and only while the target is up.

    :param nvisits: int, the number of visits
    :param per_visit: int, the exposures per visit
    :param baseline: float, the span of the campaign [days]
    :param season: float, the length of an observing season [days]
    :param cadence: float, the time between exposures of a visit [days]
    :param hour_scatter: float, the scatter of the visit time in the night
                         [days]
    :param seed: int, the seed

    :return: np.ndarray, the times [days], sorted
    """
    rng = np.random.default_rng(seed)
    nights = np.arange(int(baseline))
    visible = (nights % 365.25) < season
    chosen = np.sort(rng.choice(nights[visible], nvisits, replace=False))
    start = chosen + 0.1 + hour_scatter * rng.normal(size=nvisits)
    times = (start[:, None] + cadence * np.arange(per_visit)[None, :])
    return np.sort(times.ravel())


def keplerian_signal(time: np.ndarray, planets: Sequence[Dict[str, float]]
                     ) -> np.ndarray:
    """
    The velocity of a set of orbits

    :param time: np.ndarray, the time [days]
    :param planets: list of dict, each with P, K, and optionally e, omega
                    [rad] and tp [days] (random phase when absent)

    :return: np.ndarray, the velocity [m/s]
    """
    out = np.zeros(len(time))
    for planet in planets:
        out += kepler.rv_keplerian(time, planet['P'],
                                   planet.get('tp', float(np.min(time))),
                                   planet.get('e', 0.0),
                                   planet.get('omega', 0.0), planet['K'])
    return out


def activity_signal(time: np.ndarray, kernel: str = 'sho',
                    sigma: float = 3.0, period: float = 20.0,
                    quality: float = 2.0, seed: int = 1,
                    pars: Optional[np.ndarray] = None) -> np.ndarray:
    """
    One realisation of a gaussian process, as stellar activity

    :param time: np.ndarray, the time [days]
    :param kernel: str, the kernel
    :param sigma: float, its amplitude [m/s]
    :param period: float, its period [days]
    :param quality: float, its quality factor (SHO)
    :param seed: int, the seed
    :param pars: np.ndarray or None, the kernel parameters directly

    :return: np.ndarray, the velocity [m/s]
    """
    rng = np.random.default_rng(seed)
    if pars is None:
        pars = np.log([sigma, period, quality])
    cov = kgp.kernel_matrix(kernel, time[:, None] - time[None, :], pars)
    cov[np.diag_indices_from(cov)] += 1e-8 * sigma ** 2
    return np.linalg.cholesky(cov) @ rng.normal(size=len(time))


def _sizes(amplitude: Union[float, Tuple[float, float]], nout: int,
           rng: np.random.Generator) -> np.ndarray:
    """The sizes of nout outliers, in median errors: amplitude (1 + E)
    with E exponential of mean 0.5, or uniform in a (low, high) range"""
    if np.ndim(amplitude) == 0:
        return float(amplitude) * (1 + rng.exponential(0.5, nout))
    low, high = amplitude
    return rng.uniform(low, high, nout)


def add_outliers(time: np.ndarray, seq: np.ndarray, err: np.ndarray,
                 kind: str = 'visit', frac: float = 0.05,
                 amplitude: Union[float, Tuple[float, float]] = 8.0,
                 rng: Optional[np.random.Generator] = None,
                 avoid: Optional[np.ndarray] = None,
                 label: str = '') -> Dict[str, np.ndarray]:
    """
    Outliers of a given kind

    :param time: np.ndarray, the time [days]
    :param seq: np.ndarray, the visit of each point
    :param err: np.ndarray, the error bars
    :param kind: str, spike, visit, flare or heavy
    :param frac: float, the fraction of points (spike) or visits (visit,
                 flare) affected
    :param amplitude: float or (low, high), the typical size in units of
                      the median error, drawn as amplitude (1 + E) with E
                      exponential of mean 0.5, or uniformly in the range
                      (for heavy: the degrees of freedom)
    :param rng: np.random.Generator or None
    :param avoid: np.ndarray or None, points already made outliers: no
                  spike lands on them and no bad visit contains them
    :param label: str, a name for these outliers (clear, borderline)

    :return: dict, offset (n, what is added), mask (n, True where a point
             was made an outlier), label
    """
    rng = np.random.default_rng() if rng is None else rng
    npts = len(time)
    nseq = int(np.max(seq)) + 1
    scale = float(np.median(err))
    offset = np.zeros(npts)
    mask = np.zeros(npts, dtype=bool)
    busy = avoid is not None and bool(np.any(avoid))
    if kind == 'spike':
        nout = max(int(round(frac * npts)), 1)
        # without anything to avoid, the draws are those of earlier versions
        pool = np.where(~avoid)[0] if busy else npts
        idx = rng.choice(pool, nout, replace=False)
        sign = rng.choice([-1, 1], nout)
        offset[idx] = sign * scale * _sizes(amplitude, nout, rng)
        mask[idx] = True
    elif kind == 'visit':
        nout = max(int(round(frac * nseq)), 1)
        pool = (np.array([ss for ss in range(nseq)
                          if not np.any(avoid[seq == ss])]) if busy else nseq)
        bad = rng.choice(pool, nout, replace=False)
        sign = rng.choice([-1, 1], nout)
        size = sign * scale * _sizes(amplitude, nout, rng)
        for bb, val in zip(bad, size):
            sel = seq == bb
            # the exposures of a bad visit move together, with a little
            #   scatter of their own
            offset[sel] = val * (1 + 0.1 * rng.normal(size=np.sum(sel)))
            mask[sel] = True
    elif kind == 'flare':
        nout = max(int(round(frac * nseq)), 1)
        bad = rng.choice(nseq, nout, replace=False)
        for bb in bad:
            tflare = float(np.min(time[seq == bb]))
            decay = rng.uniform(0.5, 3.0)
            after = time >= tflare
            bump = scale * _sizes(amplitude, 1, rng)[0] * np.exp(
                -(time - tflare) / decay)
            offset[after] += bump[after]
            mask |= after & (bump > 2 * scale)
    elif kind == 'heavy':
        dof = amplitude
        extra = rng.standard_t(dof, npts) * err - rng.normal(size=npts) * err
        offset += extra
        mask = np.abs(extra) > 3 * err
    else:
        raise ValueError(f'Unknown kind of outlier: {kind}')
    return dict(offset=offset, mask=mask, label=label)


def simulate(planets: Sequence[Dict[str, float]] = (),
             template: Optional[RVData] = None,
             outliers: Sequence[Dict[str, Any]] = (),
             activity: Optional[Dict[str, Any]] = None,
             jitter: float = 0.0, visit_jitter: float = 0.0,
             err: float = 2.0, seed: int = 1,
             name: str = 'simulation', **timekw) -> Dict[str, Any]:
    """
    A velocity series with everything known

    :param planets: list of dict, the orbits (see keplerian_signal)
    :param template: RVData or None, a real series whose times, visits and
                     error bars are kept
    :param outliers: list of dict, each passed to add_outliers
    :param activity: dict or None, passed to activity_signal
    :param jitter: float, a white jitter added to the noise [m/s]
    :param visit_jitter: float, a noise shared by the exposures of a visit
                         (one draw per visit) [m/s]
    :param err: float, the error bars when there is no template [m/s]
    :param seed: int, the seed
    :param name: str, the name of the series
    :param timekw: passed to observing_times

    :return: dict, data (RVData), truth (the planets), signal, activity,
             outlier_offset, outlier_mask, and outlier_label (per point:
             the label of its outlier recipe, '' for a good point)
    """
    rng = np.random.default_rng(seed)
    if template is not None:
        time, errs = template.time.copy(), template.err.copy()
        seq = template.seq.copy()
        inst = template.inst.copy()
    else:
        time = observing_times(seed=seed, **timekw)
        errs = err * (0.8 + 0.4 * rng.random(len(time)))
        seq = None
        inst = None
    if seq is None:
        seq = np.concatenate([[0], np.cumsum(np.diff(time) > 0.3)])
    signal = keplerian_signal(time, planets)
    act = np.zeros(len(time))
    if activity is not None:
        act = activity_signal(time, seed=seed + 7, **activity)
    noise = rng.normal(size=len(time)) * np.sqrt(errs ** 2 + jitter ** 2)
    if visit_jitter > 0:
        nvisit = int(np.max(seq)) + 1
        noise += visit_jitter * rng.normal(size=nvisit)[seq]
    offset = np.zeros(len(time))
    mask = np.zeros(len(time), dtype=bool)
    labels = np.full(len(time), '', dtype=object)
    for it, spec in enumerate(outliers):
        spec = dict(spec)
        spec.setdefault('label', f'outlier{it}')
        out = add_outliers(time, seq, errs, rng=rng, avoid=mask, **spec)
        offset += out['offset']
        labels[out['mask'] & ~mask] = out['label']
        mask |= out['mask']
    rv = signal + act + noise + offset
    data = RVData(time=time, rv=rv, err=errs, inst=inst, seq=seq, name=name,
                  zero_point={name_: 0.0 for name_ in
                              (np.unique(inst) if inst is not None
                               else ['inst'])})
    return dict(data=data, truth=list(planets), signal=signal, activity=act,
                outlier_offset=offset, outlier_mask=mask,
                outlier_label=labels.astype(str))


# =============================================================================
# End of code
# =============================================================================

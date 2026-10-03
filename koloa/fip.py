#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The false inclusion probability, outlier-aware (OAFIP).

Hara et al. (2022) define the true inclusion probability of a frequency
interval I as the posterior probability that at least one signal has its
frequency in I, summed over the number of signals k:

    TIP_I = sum_k Pr{exists i <= k : f_i in I | y, k} Pr{k | y},
    FIP_I = 1 - TIP_I,

and the FIP periodogram as the FIP of sliding intervals of width 1/T.

koloa computes it for circular orbits with a sampler built so that every
step is an exact conditional draw:

1. Up to kmax signal slots, each on or off, each with a frequency on a
   grid. The on/off state and the frequency of one slot are drawn jointly,
   with every linear parameter (offsets, trend, decorrelation terms, the
   amplitudes of all signals) integrated out analytically (koloa.linear),
   so the conditional is computed over the WHOLE grid at once and a slot
   can jump from one alias to another in a single step.
2. The prior on the two amplitudes of a slot is a gaussian of variance
   tau^2, with tau itself on a log-uniform ladder: a scale mixture of
   Rayleigh distributions, close to the log-uniform prior on K that is
   common in radial velocity work, and a gaussian mixture on the linear
   parameters as Hara et al. (2022) recommend. The ladder is summed in the
   same step.
3. The noise is a white jitter, a jitter per observing sequence, and the
   outlier mixture of koloa.noise (and, with gp=, a GP of the activity as a
   finite basis whose weights are integrated out with the other linear
   parameters and whose hyperparameters are sampled with the jitters,
   koloa.gpbasis): each unit (point, or whole sequence) is
   an outlier or not, a latent indicator drawn from its exact conditional.
   Given the indicators the noise is gaussian again, which is what keeps
   step 1 exact. The outlier fraction and width are sampled too.

The FIP itself is Rao-Blackwellised: every slot update knows the exact
probability of every interval given the rest of the state, and it is that
probability, not the visited state, that is averaged. FIPs of 1e-8 are
therefore measured, not rounded to zero by a finite chain.

The number of signals has a uniform prior over 0..kmax by default.
Outlier probabilities are averaged the same way, and one minus them is the
reliability of each point that the plots use for colour.

Created on 2026-09-27

@author: artigau
"""
import multiprocessing as _mp
import sys as _sys
import threading as _threading
import time as _time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor
from concurrent.futures import wait as _wait
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.special import gammaln

from koloa import gpbasis
from koloa.data import RVData, night_index, robust_std
from koloa.linear import BaseProjection, LinearModel, sinusoid_gain
from koloa.log import log
from koloa.noise import BlockCov, DenseCov, block_gauss_loglike
from koloa.noise import block_indicator
from koloa.noise import _point_mixture_blocks
from koloa.periodogram import aliases, find_peaks, frequency_grid, window
from koloa.utils import blas_threads
from koloa.weights import hard_clip, soft_clip

# =============================================================================
# Define variables
# =============================================================================
#: a trial largest number of bytes for the cached trigonometric matrices
TRIG_CACHE_BYTES = 6.0e8
#: the most points for which a GP enters the FIP by its kernel, in the
#: covariance of the noise (koloa.noise.DenseCov), and not by a basis
KERNEL_MAX = 2000


# =============================================================================
# The grid of periods
# =============================================================================
class PeriodGrid:
    """
    The frequencies a signal may have, their prior, and their sinusoids

    The grid is regular in frequency with `oversample` points per 1/T; the
    prior is log-uniform in period (so proportional to 1/f per grid point);
    an interval of the FIP periodogram is the 2h + 1 grid points around
    each point, h = oversample // 2, which is a width of about 1/T.
    """

    def __init__(self, time: np.ndarray, pmin: float = 1.1,
                 pmax: Optional[float] = None, oversample: int = 10,
                 tref: Optional[float] = None,
                 freq: Optional[np.ndarray] = None):
        """
        :param time: np.ndarray, the time of the data [days]
        :param pmin: float, the shortest period [days]
        :param pmax: float or None, the longest period [days] (2T if None)
        :param oversample: int, grid points per 1/T
        :param tref: float or None, the reference time of the phases
        :param freq: np.ndarray or None, a grid to use as it is
        """
        self.time = np.asarray(time, dtype=float)
        self.baseline = float(np.ptp(self.time))
        if freq is None:
            freq = frequency_grid(self.time, pmin, pmax, oversample)
        self.freq = np.asarray(freq, dtype=float)
        self.period = 1.0 / self.freq
        self.oversample = oversample
        self.halfbin = max(int(oversample) // 2, 0)
        self.tref = float(np.mean(self.time)) if tref is None else tref
        logp = -np.log(self.freq)
        self.logprior = logp - np.logaddexp.reduce(logp)
        self._trig = None
        self._fam = None

    @property
    def size(self) -> int:
        """The number of grid points"""
        return len(self.freq)

    def trig(self):
        """
        The sinusoids of the grid and their products, cached

        :return: tuple, (n x 2G) [cos | sin] and (n x 3G)
                 [cos^2 | sin^2 | cos sin]
        """
        if self._trig is None:
            phase = 2 * np.pi * np.outer(self.time - self.tref, self.freq)
            cos, sin = np.cos(phase), np.sin(phase)
            trig_cs = np.ascontiguousarray(np.hstack([cos, sin]))
            trig_sq = np.ascontiguousarray(np.hstack([cos * cos, sin * sin,
                                                      cos * sin]))
            self._trig = (trig_cs, trig_sq)
        return self._trig

    def columns(self, idx: int) -> np.ndarray:
        """
        The two columns of one grid frequency

        :param idx: int, the grid index

        :return: np.ndarray, (n x 2) cos and sin
        """
        trig_cs = self.trig()[0]
        return trig_cs[:, [idx, idx + self.size]]

    def family(self) -> Dict[str, np.ndarray]:
        """
        The alias family of every grid frequency, as segments of the grid

        The family of f is f and its aliases |f - fs| and f + fs for every
        sampling frequency fs of ALIAS_FREQUENCIES (1 day, 1 year, 1 month):
        up to seven intervals of the grid, merged where they overlap. What
        lies outside them is a set of gaps whose bounds depend on the grid
        only, so they are computed once.

        :return: dict, centre (G x M grid indices), inside (G x M, the alias
                 is on the grid), and the gaps: start, end (G x M+1) and ok
                 (the gap exists)
        """
        if self._fam is None:
            from koloa.periodogram import ALIAS_FREQUENCIES
            size, hbin = self.size, self.halfbin
            members = [self.freq]
            for fsamp in ALIAS_FREQUENCIES.values():
                members += [np.abs(self.freq - fsamp), self.freq + fsamp]
            members = np.array(members).T
            idx = np.clip(np.searchsorted(self.freq, members), 1, size - 1)
            idx = np.where(np.abs(self.freq[idx - 1] - members)
                           < np.abs(self.freq[idx] - members), idx - 1, idx)
            half = 0.5 * np.median(np.diff(self.freq)) if size > 1 else 0.0
            inside = ((members >= self.freq[0] - half)
                      & (members <= self.freq[-1] + half))
            start = np.where(inside, np.clip(idx - hbin, 0, size - 1), size)
            end = np.where(inside, np.clip(idx + hbin, 0, size - 1), -1)
            order = np.argsort(start, axis=1, kind='stable')
            start = np.take_along_axis(start, order, axis=1)
            end = np.take_along_axis(end, order, axis=1)
            nvalid = inside.sum(axis=1)
            run = np.maximum.accumulate(end, axis=1)
            nrow = len(self.freq)
            gstart = np.hstack([np.zeros((nrow, 1), int), run + 1])
            gend = np.hstack([start - 1, np.full((nrow, 1), size - 1)])
            # the gaps up to the one after the last interval
            jcol = np.arange(members.shape[1] + 1)[None, :]
            gend = np.where(jcol == nvalid[:, None], size - 1, gend)
            self._fam = dict(centre=idx, inside=inside,
                             start=np.clip(gstart, 0, size),
                             end=np.clip(gend, -1, size - 1),
                             ok=jcol <= nvalid[:, None])
        return self._fam

    def bin_family(self, pon: np.ndarray, poff: float, total: float,
                   covered: Sequence[int] = ()) -> np.ndarray:
        """
        The FIP of the alias family of every grid frequency, given one
        slot's conditional: the probability that the slot is off or outside
        every interval of the family (the sum of the gaps between them), and
        zero if another active slot is already in one of them

        :param pon: np.ndarray, (G) the slot's unnormalised probability at
                    each grid point
        :param poff: float, its unnormalised probability of being off
        :param total: float, the normalisation
        :param covered: list of int, the grid indices of the other active
                        slots

        :return: np.ndarray, (G) the FIP of each family
        """
        fam = self.family()
        size = len(pon)
        cum = np.concatenate([[0.0], np.cumsum(pon)])
        rcum = np.concatenate([np.cumsum(pon[::-1])[::-1], [0.0]])
        gstart, gend = fam['start'], fam['end']
        # a gap that runs to the end is summed from the right (precision)
        seg = np.where(gend >= size - 1, rcum[gstart],
                       cum[gend + 1] - cum[gstart])
        use = fam['ok'] & (gend >= gstart)
        out = (poff + np.sum(np.where(use, seg, 0.0), axis=1)) / total
        for gidx in covered:
            hit = np.any(fam['inside'] & (np.abs(fam['centre'] - gidx)
                                          <= self.halfbin), axis=1)
            out[hit] = 0.0
        return out

    def bin_fip(self, pon: np.ndarray, poff: float, total: float,
                covered: Sequence[int] = ()) -> np.ndarray:
        """
        The FIP of every sliding interval, given one slot's conditional

        The probability that the interval around grid point k is empty is
        the probability that the slot is off or outside it, computed as the
        sum of what lies left of it and right of it (never as one minus a
        number close to one), and zero if another active slot is already
        inside it.

        :param pon: np.ndarray, (G) the slot's unnormalised probability at
                    each grid point
        :param poff: float, its unnormalised probability of being off
        :param total: float, the normalisation
        :param covered: list of int, the grid indices of the other active
                        slots

        :return: np.ndarray, (G) the FIP of each interval
        """
        size, hbin = len(pon), self.halfbin
        cum = np.concatenate([[0.0], np.cumsum(pon)])
        rcum = np.concatenate([np.cumsum(pon[::-1])[::-1], [0.0]])
        kidx = np.arange(size)
        left = cum[np.clip(kidx - hbin, 0, size)]
        right = rcum[np.clip(kidx + hbin + 1, 0, size)]
        fip = (poff + left + right) / total
        for gidx in covered:
            fip[max(gidx - hbin, 0):gidx + hbin + 1] = 0.0
        return fip


# =============================================================================
# The design matrix that is not a signal
# =============================================================================
def base_design(data: RVData, trend: int = 1,
                regressors: Optional[Dict[str, np.ndarray]] = None
                ) -> Dict[str, Any]:
    """
    Offsets, a polynomial trend and decorrelation terms, as columns

    Every one of these has a wide gaussian prior. They are common to every
    model with and without signals, so the width cancels in the FIP.

    :param data: RVData, the series
    :param trend: int, the degree of the polynomial in time (0: none)
    :param regressors: dict or None, name: values of the decorrelation
                       terms (each is standardised)

    :return: dict, design (n x p), prior_var (p), names (p)
    """
    cols, names = [], []
    for it, inst in enumerate(data.instruments):
        cols.append((data.inst == inst).astype(float))
        names.append(f'offset {inst}')
    tnorm = (data.time - data.tref) / max(data.baseline, 1e-9)
    for deg in range(1, trend + 1):
        cols.append(tnorm ** deg)
        names.append(f'trend t^{deg}')
    for name, values in (regressors or {}).items():
        values = np.asarray(values, dtype=float)
        std = np.std(values)
        cols.append((values - np.mean(values)) / (std if std > 0 else 1.0))
        names.append(f'decorrelation {name}')
    scale = max(robust_std(data.rv), float(np.median(data.err)))
    prior_var = np.full(len(cols), (100.0 * scale) ** 2)
    return dict(design=np.array(cols).T, prior_var=prior_var, names=names)


def prior_on_k(kmax: int, prior_k: Union[str, Sequence[float]] = 'uniform'
               ) -> np.ndarray:
    """
    The log prior of each number of signals, 0..kmax

    :param kmax: int, the largest number of signals
    :param prior_k: 'uniform', or a list of kmax + 1 probabilities

    :return: np.ndarray, (kmax + 1) log probabilities
    """
    if isinstance(prior_k, str):
        if prior_k != 'uniform':
            raise ValueError(f'Unknown prior on k: {prior_k}')
        prob = np.ones(kmax + 1)
    else:
        prob = np.asarray(prior_k, dtype=float)
        if len(prob) != kmax + 1:
            raise ValueError('prior_k needs kmax + 1 values')
    return np.log(prob / np.sum(prob))


def _log_binom(nval: int, kval: int) -> float:
    """log of the binomial coefficient"""
    return float(gammaln(nval + 1) - gammaln(kval + 1)
                 - gammaln(nval - kval + 1))


# =============================================================================
# The result
# =============================================================================
@dataclass
class FIPResult:
    """
    What the FIP computation returns
    """
    freq: np.ndarray
    fip: np.ndarray
    density: np.ndarray
    method: str
    pk: np.ndarray = None
    outlier_prob: np.ndarray = None
    unit: str = 'none'
    chains: Dict[str, np.ndarray] = field(default_factory=dict)
    signals: List[Dict[str, Any]] = field(default_factory=list)
    peaks: List[Dict[str, Any]] = field(default_factory=list)
    settings: Dict[str, Any] = field(default_factory=dict)
    chain_fip: np.ndarray = None
    rhat: Dict[str, float] = field(default_factory=dict)
    runtime: float = 0.0
    #: the FIP of the alias family of each frequency: no signal at it NOR
    #:  at any of its aliases (1 day, 1 year, 1 month); whether there is a
    #:  planet, whichever alias it is (None for a FIP without a sampler)
    family: np.ndarray = None
    #: every active slot of every recorded sweep: (grid index, K, phase),
    #:  the phase that of A cos + B sin about the grid's reference time,
    #:  atan2(B, A)
    slots: np.ndarray = None

    @property
    def period(self) -> np.ndarray:
        """The periods of the grid [days]"""
        return 1.0 / self.freq

    @property
    def log10fip(self) -> np.ndarray:
        """log10 of the FIP, floored at 1e-300"""
        return np.log10(np.clip(self.fip, 1e-300, 1.0))

    @property
    def reliability(self) -> Optional[np.ndarray]:
        """One minus the outlier probability of each point"""
        if self.outlier_prob is None:
            return None
        return 1.0 - self.outlier_prob

    def fip_at(self, period: float) -> float:
        """
        The FIP of the interval around a period

        :param period: float, the period [days]

        :return: float, the FIP
        """
        return float(self.fip[np.argmin(np.abs(self.freq - 1.0 / period))])

    def fip_containing(self, period: float, width: float) -> float:
        """
        The FIP of the best interval that holds a period: the lowest FIP of
        the intervals centred within half a width of it (a fitted period
        sits anywhere in its peak, not on the grid point of the best
        interval, and the grid point nearest it can belong to a weaker one)

        :param period: float, the period [days]
        :param width: float, the width of the intervals, 1/T [1/day]

        :return: float, the FIP
        """
        near = np.abs(self.freq - 1.0 / period) <= 0.5 * width
        if not np.any(near):
            return self.fip_at(period)
        return float(np.min(self.fip[near]))

    def family_containing(self, period: float, width: float
                          ) -> Optional[float]:
        """
        The FIP of the period OR any of its aliases (1 day, 1 year, 1
        month): the probability that there is no signal at any of them. It
        says whether there is a planet, whichever alias it is; which alias
        is alias_odds()

        :param period: float, the period [days]
        :param width: float, the width of the intervals, 1/T [1/day]

        :return: float, or None when the FIP has no family
        """
        if self.family is None:
            return None
        near = np.abs(self.freq - 1.0 / period) <= 0.5 * width
        if not np.any(near):
            near = np.argmin(np.abs(self.freq - 1.0 / period))
        return float(np.min(self.family[near]))

    def alias_odds(self, period: float, width: float
                   ) -> List[Dict[str, Any]]:
        """
        The period and each of its aliases on the grid, with the share of
        the probability that each holds: its true inclusion probability
        over the sum of those of the family (they are almost exclusive: a
        signal is at one alias or another). An alias within 1.5 widths of
        the period (the yearly one, with a baseline of about a year) is
        part of its peak, not an alternative, and is left out

        :param period: float, the period [days]
        :param width: float, the width of the intervals, 1/T [1/day]

        :return: list of dict, name ('P' for the period itself, else the
                 sampling of the alias, e.g. '1 day'), period, fip, tip,
                 share; the most probable first
        """
        from koloa.periodogram import aliases
        fmin, fmax = float(self.freq[0]), float(self.freq[-1])
        members = [dict(name='P', freq=1.0 / period)]
        members += [dict(name=al['name'], freq=al['freq'])
                    for al in aliases(1.0 / period, fmin, fmax)
                    if abs(al['freq'] - 1.0 / period) > 1.5 * width]
        out = []
        for mem in members:
            fipv = self.fip_containing(1.0 / mem['freq'], width)
            out.append(dict(name=mem['name'], period=1.0 / mem['freq'],
                            fip=fipv, tip=max(1.0 - fipv, 0.0)))
        total = sum(mem['tip'] for mem in out)
        for mem in out:
            mem['share'] = mem['tip'] / total if total > 0 else np.nan
        out.sort(key=lambda mem: -mem['tip'])
        return out

    def gp_summary(self) -> Dict[str, Tuple[float, float, float]]:
        """
        The hyperparameters of the GP, as (median, 16th, 84th percentiles)
        of the chains (empty without a GP)

        :return: dict, 'gp<i>_<name>': tuple
        """
        out = {}
        for key, val in self.chains.items():
            if key.startswith('gp') and len(val):
                pct = np.percentile(val, [50, 16, 84])
                out[key] = tuple(float(item) for item in pct)
        return out

    def best(self) -> Dict[str, Any]:
        """
        The most significant peak

        :return: dict, see `peaks`
        """
        return self.peaks[0] if self.peaks else {}

    def summary(self) -> str:
        """
        A few lines that say what was found

        :return: str, the summary
        """
        lines = [f'FIP ({self.method})']
        if self.pk is not None:
            pks = ', '.join(f'P(k={it})={val:.3g}'
                            for it, val in enumerate(self.pk))
            lines.append(f'  number of signals: {pks}')
        for peak in self.peaks[:5]:
            fam = peak.get('family_fip')
            lines.append(f'  P = {peak["period"]:.4f} d  FIP = '
                         f'{peak["fip"]:.2e}  window = '
                         f'{peak["window"]:.3f}'
                         + (f'  FIP of P or an alias = {fam:.2e}'
                            if fam is not None else ''))
        for key, (mid, low, high) in self.gp_summary().items():
            lines.append(f'  {key}: {mid:.3g} ({low:.3g} to {high:.3g})')
        if self.outlier_prob is not None:
            lines.append(f'  expected outlier points ({self.unit} model): '
                         f'{self.settings.get("expected_outliers", 0):.1f}')
        return '\n'.join(lines)


# =============================================================================
# One chain of the sampler
# =============================================================================
class _Chain:
    """
    The state of one chain and the moves that update it
    """

    def __init__(self, data: RVData, grid: PeriodGrid, cfg: Dict[str, Any],
                 seed: int):
        self.cfg = cfg
        self.rng = np.random.default_rng(seed)
        self.y, self.err = data.rv, data.err
        self.npts = data.n
        self.block, self.nblock = data.seq, data.nseq
        self.indicator = block_indicator(self.block, self.nblock)
        base = base_design(data, cfg['trend'], cfg['regressors'])
        self.fixed_base, self.fixed_var = base['design'], base['prior_var']
        # the GP of the activity: more base columns, their prior from its
        #   hyperparameters (koloa.gpbasis)
        self.time = data.time
        self.inst = data.inst
        self.gp_comp = cfg.get('gp') or []
        self.gp_val = gpbasis.initial(self.gp_comp)
        # the GP as basis columns, or (kernel mode) its kernel in the noise
        self.gp_mode = cfg.get('gp_mode', 'basis')
        self.kmat = None
        self._set_base()
        self.grid = grid
        self.trig_cs, self.trig_sq = grid.trig()
        self.gsize = grid.size
        self.tau2 = np.asarray(cfg['tau']) ** 2
        self.logprior_tau = np.full(len(self.tau2), -np.log(len(self.tau2)))
        self.kmax = cfg['kmax']
        self.logprior_k = prior_on_k(self.kmax, cfg['prior_k'])
        self.unit = cfg['unit']
        # which kinds of outliers: single points, whole sequences, or both
        self.use_pt = self.unit in ('point', 'both')
        self.use_seq = self.unit in ('sequence', 'both')
        # the state
        self.active = np.zeros(self.kmax, dtype=bool)
        self.gidx = np.zeros(self.kmax, dtype=int)
        self.tidx = np.zeros(self.kmax, dtype=int)
        init = cfg['init']
        self.jit = init['jitter']
        self.sjit = init['seq_jitter']
        self.wpt, self.wseq = init['width'], init['width']
        self.frac_pt, self.frac_seq = init['frac'], init['frac']
        self.q_pt = init['q_pt'].copy()
        self.q_seq = init['q_seq'].copy()
        self.steps = dict(jit=0.3, sjit=0.3, wpt=0.3, wseq=0.3)
        self.gp_keys = gpbasis.names(self.gp_comp)
        self.steps.update({key: 0.2 for key in self.gp_keys})
        self.accept = {key: [0, 0] for key in self.steps}
        # the accumulators
        size = grid.size
        self.fip_acc = np.zeros(size)
        self.dens_acc = np.zeros(size)
        self.fam_acc = np.zeros(size)
        self.nacc = 0
        self.qprob_acc = np.zeros(self.npts)
        self.qseq_acc = np.zeros(self.nblock)
        self.nq = 0
        self.records = dict(k=[], jitter=[], seq_jitter=[], width=[],
                            width_seq=[], frac=[], frac_seq=[], noutlier=[],
                            logz=[])
        self.records.update({key: [] for key in self.gp_keys})
        self.slot_records = []

    def _set_base(self):
        """the base columns: offsets, trend and regressors, then the GP (or,
        in the kernel mode, the GP's covariance at the points)"""
        if not self.gp_comp:
            self.base, self.base_var = self.fixed_base, self.fixed_var
            return
        if self.gp_mode == 'kernel':
            self.base, self.base_var = self.fixed_base, self.fixed_var
            self.kmat = gpbasis.matrix(self.gp_comp, self.gp_val, self.time,
                                       self.inst)
            return
        cols, var = gpbasis.columns(self.gp_comp, self.gp_val, self.time,
                                    self.inst)
        self.base = np.hstack([self.fixed_base, cols])
        self.base_var = np.concatenate([self.fixed_var, var])

    def update_gp(self, cov: BlockCov, model: LinearModel, logz: float):
        """
        One Metropolis step on each hyperparameter of the GP (in its
        logarithm), the linear parameters (the GP's weights among them)
        integrated out

        :return: tuple, the noise covariance (new in the kernel mode), the
                 model and its log evidence after the steps
        """
        for ic, comp in enumerate(self.gp_comp):
            for name, prior in comp['priors'].items():
                if prior['kind'] == 'fixed':
                    continue
                key = f'gp{ic}_{name}'
                current = self.gp_val[ic][name]
                prop = float(np.exp(np.log(current) + self.steps[key]
                                    * self.rng.normal()))
                self.accept[key][1] += 1
                lp_new = gpbasis.log_prior(prior, prop)
                if not np.isfinite(lp_new):
                    continue
                lp_old = gpbasis.log_prior(prior, current)
                old = (self.base, self.base_var, self.kmat)
                self.gp_val[ic][name] = prop
                self._set_base()
                # the kernel mode: the GP is in the noise, which changes
                cov_new = cov if self.kmat is None else self.make_cov()
                model_new = self.full_model(cov_new)
                if (np.log(self.rng.random())
                        < model_new.logz + lp_new - logz - lp_old):
                    cov, model, logz = cov_new, model_new, model_new.logz
                    self.accept[key][0] += 1
                else:
                    self.gp_val[ic][name] = current
                    self.base, self.base_var, self.kmat = old
        return cov, model, logz

    @property
    def nbad(self) -> int:
        """The number of points that are outliers or in an outlier sequence"""
        bad = self.q_pt > 0
        if self.use_seq:
            bad = bad | (self.q_seq[self.block] > 0)
        return int(np.sum(bad))

    # -------------------------------------------------------------------------
    def make_cov(self, jit=None, sjit=None, wpt=None, wseq=None
                 ) -> BlockCov:
        """The noise covariance of the current (or a proposed) state"""
        jit = self.jit if jit is None else jit
        sjit = self.sjit if sjit is None else sjit
        wpt = self.wpt if wpt is None else wpt
        wseq = self.wseq if wseq is None else wseq
        diag = self.err ** 2 + jit ** 2
        blockval = np.full(self.nblock, sjit ** 2)
        if self.use_pt:
            diag = diag + self.q_pt * wpt ** 2
        if self.use_seq:
            blockval = blockval + self.q_seq * wseq ** 2
        cov = BlockCov(diag, self.block, blockval, self.indicator)
        return cov if self.kmat is None else DenseCov(cov, self.kmat)

    def design(self, exclude: Optional[int] = None):
        """The design matrix and prior variances of the active slots"""
        cols, var = [self.base], [self.base_var]
        for jj in np.where(self.active)[0]:
            if jj == exclude:
                continue
            cols.append(self.grid.columns(self.gidx[jj]))
            var.append(np.full(2, self.tau2[self.tidx[jj]]))
        return np.hstack(cols), np.concatenate(var)

    def _slot_products(self, cov: BlockCov, gidx: int):
        """V^-1 [cos sin] of one frequency, and its products with the grid"""
        wcol = cov.solve(self.grid.columns(gidx))
        return wcol, wcol.T @ self.trig_cs

    def full_model(self, cov: BlockCov) -> LinearModel:
        """The linear model of the current state"""
        design, var = self.design()
        return LinearModel(design, self.y, cov, var)

    # -------------------------------------------------------------------------
    def sweep(self, record: bool):
        """One update of every part of the state"""
        cov = self.make_cov()
        size = self.gsize
        wy = cov.solve(self.y)
        cc, ss, cs = cov.grid_quad(self.trig_sq, self.trig_cs, size)
        ry = wy @ self.trig_cs
        ry_c, ry_s = ry[:size], ry[size:]
        # what the base (and the GP) takes from every sinusoid of the grid,
        #   once for all the slots
        proj = BaseProjection(self.base, cov.solve(self.base), self.base_var,
                              self.y, self.trig_cs, size)
        # what each active slot contributes to the products
        wslot, rslot = {}, {}
        for jj in np.where(self.active)[0]:
            wslot[jj], rslot[jj] = self._slot_products(cov, self.gidx[jj])
        # ---------------------------------------------------------------------
        # 1. every slot: on/off and frequency (and tau), linear parameters
        #    integrated out
        # ---------------------------------------------------------------------
        ntau = len(self.tau2)
        for jj in range(self.kmax):
            others = [ii for ii in np.where(self.active)[0] if ii != jj]
            if others:
                gain = proj.gain(
                    cc, ss, cs, ry_c, ry_s, self.tau2,
                    others=np.hstack([self.grid.columns(self.gidx[ii])
                                      for ii in others]),
                    wothers=np.hstack([wslot[ii] for ii in others]),
                    rothers=np.vstack([rslot[ii] for ii in others]),
                    var_others=np.concatenate([
                        np.full(2, self.tau2[self.tidx[ii]])
                        for ii in others]),
                    value=self.y)
            else:
                gain = proj.gain(cc, ss, cs, ry_c, ry_s, self.tau2)
            kother = len(others)
            log_on = (self.logprior_k[kother + 1]
                      - _log_binom(self.kmax, kother + 1))
            log_off = (self.logprior_k[kother]
                       - _log_binom(self.kmax, kother))
            gain += log_on + self.logprior_tau[0]
            gain += self.grid.logprior[:, None]
            lmax = max(float(np.max(gain)), log_off)
            wgt = np.exp(gain - lmax)
            woff = float(np.exp(log_off - lmax))
            pon = np.sum(wgt, axis=1)
            total = woff + float(np.sum(pon))
            if record:
                covered = [self.gidx[ii] for ii in others]
                self.fip_acc += self.grid.bin_fip(pon, woff, total, covered)
                self.fam_acc += self.grid.bin_family(pon, woff, total,
                                                     covered)
                self.dens_acc += pon / total
                self.nacc += 1
            # the draw: off, or a frequency, then a rung of the tau ladder
            uval = self.rng.random() * total
            if uval < woff:
                self.active[jj] = False
                wslot.pop(jj, None)
                rslot.pop(jj, None)
            else:
                cum = np.cumsum(pon)
                gnew = min(int(np.searchsorted(cum, uval - woff)), size - 1)
                tcum = np.cumsum(wgt[gnew])
                tnew = min(int(np.searchsorted(tcum, self.rng.random()
                                               * tcum[-1])), ntau - 1)
                self.active[jj] = True
                self.gidx[jj], self.tidx[jj] = gnew, tnew
                wslot[jj], rslot[jj] = self._slot_products(cov, gnew)
        # ---------------------------------------------------------------------
        # 2. the noise parameters, linear parameters still integrated out
        # ---------------------------------------------------------------------
        model = self.full_model(cov)
        logz = model.logz
        nout = dict(wpt=int(np.sum(self.q_pt)), wseq=int(np.sum(self.q_seq)))
        for name in ('jit', 'sjit', 'wpt', 'wseq'):
            if not self.cfg['free'][name]:
                continue
            low, high = self.cfg['bounds'][name]
            if name in nout and nout[name] == 0:
                # no outlier of this kind: its width is not in the
                #   likelihood at all, and is drawn from its prior
                setattr(self, name, float(np.exp(self.rng.uniform(
                    np.log(low), np.log(high)))))
                continue
            current = getattr(self, name)
            prop = float(np.exp(np.log(current)
                                + self.steps[name] * self.rng.normal()))
            self.accept[name][1] += 1
            if not low <= prop <= high:
                continue
            cov_new = self.make_cov(**{name: prop})
            model_new = self.full_model(cov_new)
            if np.log(self.rng.random()) < model_new.logz - logz:
                setattr(self, name, prop)
                cov, model, logz = cov_new, model_new, model_new.logz
                self.accept[name][0] += 1
        if self.gp_comp:
            cov, model, logz = self.update_gp(cov, model, logz)
        # ---------------------------------------------------------------------
        # 3. the outlier indicators, linear parameters integrated out, then
        #    a draw of the linear parameters for the amplitude records
        # ---------------------------------------------------------------------
        if self.use_pt or self.use_seq:
            if self.npts <= self.cfg.get('dense_max', 3000):
                self.update_outliers_collapsed(cov, model, record)
                cov = self.make_cov()
                model = self.full_model(cov)
            else:
                beta = model.draw(self.rng)
                self.update_outliers(self.y - model.design @ beta, record)
        beta = model.draw(self.rng)
        # ---------------------------------------------------------------------
        # what is kept
        # ---------------------------------------------------------------------
        if record:
            self.records['k'].append(int(np.sum(self.active)))
            self.records['jitter'].append(self.jit)
            self.records['seq_jitter'].append(self.sjit)
            self.records['width'].append(self.wpt if self.use_pt
                                         else self.wseq)
            self.records['width_seq'].append(self.wseq)
            self.records['frac'].append(self.frac_pt if self.use_pt
                                        else self.frac_seq)
            self.records['frac_seq'].append(self.frac_seq)
            self.records['noutlier'].append(self.nbad)
            self.records['logz'].append(logz)
            for ic, comp in enumerate(self.gp_comp):
                for name, prior in comp['priors'].items():
                    if prior['kind'] != 'fixed':
                        self.records[f'gp{ic}_{name}'].append(
                            self.gp_val[ic][name])
            nbase = self.base.shape[1]
            col = nbase
            for jj in np.where(self.active)[0]:
                amp_a, amp_b = beta[col], beta[col + 1]
                self.slot_records.append((int(self.gidx[jj]),
                                          float(np.hypot(amp_a, amp_b)),
                                          float(np.arctan2(amp_b, amp_a))))
                col += 2

    def update_outliers_collapsed(self, cov: BlockCov, model: LinearModel,
                                  record: bool):
        """
        Draw every outlier indicator from its exact conditional with the
        linear parameters integrated out, one at a time

        Flipping an indicator adds +-W^2 u u^T to the covariance of the data
        once the linear parameters are integrated out,
            Sigma_y = V + X Sigma X^T,
        with u the indicator vector of the visit (or of the point). The
        change of the evidence then follows from the determinant lemma and
        Sherman-Morrison:
            Delta ln Z = -1/2 ln(1 + d s) + 1/2 d t^2 / (1 + d s),
            s = u^T Sigma_y^-1 u,  t = u^T Sigma_y^-1 y,
        and Sigma_y^-1 is carried along by a rank-one update when an
        indicator flips. Without the linear parameters in the conditioning,
        an indicator no longer has to wait for the orbit to move before it
        can change, which is what made the indicators mix slowly when they
        were drawn given a draw of the orbit.

        Visits are updated first, then points, and the probability that
        each point is bad (an outlier itself or in an outlier visit) is
        Rao-Blackwellised on the point's own indicator.
        """
        # Sigma_y^-1 = V^-1 - V^-1 X K^-1 X^T V^-1, dense
        vinv = cov.dense_inv()
        wdes = model.wdesign
        sinv = vinv - wdes @ model.kinv() @ wdes.T
        aval = sinv @ self.y
        apri, bpri = self.cfg['frac_prior']
        pseq = np.zeros(self.nblock)
        members = self._members()
        if self.use_seq:
            lodds = np.log(self.frac_seq) - np.log1p(-self.frac_seq)
            w2 = self.wseq ** 2
            for bb, idx in enumerate(members):
                svec = np.sum(sinv[:, idx], axis=1)
                sval = float(np.sum(svec[idx]))
                tval = float(np.sum(aval[idx]))
                # the change of ln Z of turning this visit into an outlier,
                #   from its current state
                delta = w2 if self.q_seq[bb] == 0 else -w2
                gain = (-0.5 * np.log1p(delta * sval)
                        + 0.5 * delta * tval ** 2 / (1 + delta * sval))
                logit_on = lodds + (gain if self.q_seq[bb] == 0 else -gain)
                pseq[bb] = 1.0 / (1.0 + np.exp(-logit_on))
                new = float(self.rng.random() < pseq[bb])
                if new != self.q_seq[bb]:
                    denom = 1 + delta * sval
                    sinv -= delta * np.outer(svec, svec) / denom
                    aval -= delta * svec * tval / denom
                    self.q_seq[bb] = new
        ppt = np.zeros(self.npts)
        if self.use_pt:
            lodds = np.log(self.frac_pt) - np.log1p(-self.frac_pt)
            w2 = self.wpt ** 2
            for ii in range(self.npts):
                sval = float(sinv[ii, ii])
                tval = float(aval[ii])
                delta = w2 if self.q_pt[ii] == 0 else -w2
                gain = (-0.5 * np.log1p(delta * sval)
                        + 0.5 * delta * tval ** 2 / (1 + delta * sval))
                logit_on = lodds + (gain if self.q_pt[ii] == 0 else -gain)
                ppt[ii] = 1.0 / (1.0 + np.exp(-logit_on))
                new = float(self.rng.random() < ppt[ii])
                if new != self.q_pt[ii]:
                    svec = sinv[:, ii].copy()
                    denom = 1 + delta * sval
                    sinv -= delta * np.outer(svec, svec) / denom
                    aval -= delta * svec * tval / denom
                    self.q_pt[ii] = new
        if record:
            if self.use_seq and self.use_pt:
                inseq = self.q_seq[self.block] > 0
                self.qprob_acc += np.where(inseq, 1.0, ppt)
            elif self.use_seq:
                self.qprob_acc += pseq[self.block]
            else:
                self.qprob_acc += ppt
            self.qseq_acc += pseq
            self.nq += 1
        if self.use_pt:
            nout = float(np.sum(self.q_pt))
            self.frac_pt = float(np.clip(self.rng.beta(
                apri + nout, bpri + self.npts - nout), 1e-9, 1 - 1e-9))
        if self.use_seq:
            nout = float(np.sum(self.q_seq))
            self.frac_seq = float(np.clip(self.rng.beta(
                apri + nout, bpri + self.nblock - nout), 1e-9, 1 - 1e-9))

    def _members(self) -> List[np.ndarray]:
        """The point indices of each visit (cached)"""
        if not hasattr(self, '_member_cache'):
            order = np.argsort(self.block, kind='stable')
            bounds = np.cumsum(np.bincount(self.block,
                                           minlength=self.nblock))
            self._member_cache = np.split(order, bounds[:-1])
        return self._member_cache

    def update_outliers(self, resid: np.ndarray, record: bool):
        """
        Draw the outlier indicators given a draw of the linear parameters
        (the path for series too long for a dense inverse), then the
        outlier fractions

        Sequence indicators first, given the point ones (a Gaussian block
        likelihood each); then point indicators given the sequence ones,
        a whole block's configuration at a time when a sequence term
        couples its points.
        """
        diag = self.err ** 2 + self.jit ** 2
        apri, bpri = self.cfg['frac_prior']
        pseq = np.zeros(self.nblock)
        if self.use_seq:
            dfull = diag + (self.q_pt * self.wpt ** 2 if self.use_pt else 0.0)
            a0 = np.full(self.nblock, self.sjit ** 2)
            good = block_gauss_loglike(resid, dfull, self.block, a0,
                                       self.nblock)
            bad = block_gauss_loglike(resid, dfull, self.block,
                                      a0 + self.wseq ** 2, self.nblock)
            lgood = np.log1p(-self.frac_seq) + good
            lbad = np.log(self.frac_seq) + bad
            pseq = np.exp(lbad - np.logaddexp(lgood, lbad))
            self.q_seq = (self.rng.random(self.nblock) < pseq).astype(float)
        ppt = np.zeros(self.npts)
        if self.use_pt:
            ablock = np.full(self.nblock, self.sjit ** 2)
            if self.use_seq:
                ablock = ablock + self.q_seq * self.wseq ** 2
            _, logp, _, drawn = _point_mixture_blocks(
                resid, diag, self.block, self.nblock, ablock,
                np.log(self.frac_pt), np.log1p(-self.frac_pt), self.wpt,
                rng=self.rng)
            ppt = np.exp(logp)
            self.q_pt = drawn
        if record:
            # a point is bad when its sequence is, or else with the
            #   probability of its own indicator given the sequences
            inseq = (self.q_seq[self.block] > 0) if self.use_seq else \
                np.zeros(self.npts, dtype=bool)
            self.qprob_acc += np.where(inseq, 1.0, ppt)
            self.qseq_acc += pseq
            self.nq += 1
        if self.use_pt:
            nout = float(np.sum(self.q_pt))
            self.frac_pt = float(np.clip(self.rng.beta(
                apri + nout, bpri + self.npts - nout), 1e-9, 1 - 1e-9))
        if self.use_seq:
            nout = float(np.sum(self.q_seq))
            self.frac_seq = float(np.clip(self.rng.beta(
                apri + nout, bpri + self.nblock - nout), 1e-9, 1 - 1e-9))

    def adapt(self):
        """Tune the proposal steps towards an acceptance of about 0.4"""
        for name, (acc, tot) in self.accept.items():
            if tot < 20:
                continue
            rate = acc / tot
            if rate > 0.55:
                self.steps[name] *= 1.4
            elif rate < 0.25:
                self.steps[name] /= 1.4
            self.accept[name] = [0, 0]

    def output(self) -> Dict[str, Any]:
        """Everything the chain measured"""
        return dict(fip_acc=self.fip_acc, dens_acc=self.dens_acc,
                    fam_acc=self.fam_acc,
                    nacc=self.nacc, qprob_acc=self.qprob_acc,
                    qseq_acc=self.qseq_acc, nq=self.nq,
                    records={key: np.array(val)
                             for key, val in self.records.items()},
                    slot_records=np.array(self.slot_records).reshape(-1, 3))


#: a function told the progress of every FIP (label, done, total, seconds),
#: when a program runs FIPs in its own process (koloa's GUI)
PROGRESS_HOOK = None
#: a function that says (True) when the FIP running should stop, asked
#: every second while its chains run in processes (koloa's GUI): they are
#: ended at once and FIPCancelled raised
CANCEL_HOOK = None


class FIPCancelled(RuntimeError):
    """A FIP stopped from outside (CANCEL_HOOK)"""
#: the sweeps done by each chain of the FIP that runs, shared with the
#: processes of the chains (set by _share_counter in each of them)
_COUNTER = None


def _share_counter(counter):
    """in the process of a chain: where it counts its sweeps"""
    global _COUNTER
    _COUNTER = counter


class _Progress:
    """
    The sweeps of the chains of a FIP as they go: a bar (tqdm) on a
    terminal, a log line at every tenth otherwise (a pipe); for koloa's GUI
    (KOLOA_PROGRESS=gui), also a 'progress: label | done | total | seconds'
    line every two seconds, which it draws as a bar
    """

    def __init__(self, label: str, total: int):
        import os
        from koloa import log as klog
        self.label, self.total = label, total
        self.start, self.shown = _time.time(), 0
        self.gui = os.environ.get('KOLOA_PROGRESS', '') == 'gui'
        self.last = 0.0
        self.bar = None
        if klog.VERBOSE and _sys.stderr.isatty():
            try:
                from tqdm import tqdm
                self.bar = tqdm(total=total, desc=label, unit='sweep',
                                dynamic_ncols=True, leave=True)
            except ImportError:
                self.bar = None

    def update(self, done: int):
        done = min(int(done), self.total)
        now = _time.time()
        if PROGRESS_HOOK is not None:
            PROGRESS_HOOK(self.label, done, self.total, now - self.start)
        if self.gui and (now - self.last >= 2.0 or done >= self.total):
            self.last = now
            log(f'progress: {self.label} | {done} | {self.total} | '
                f'{now - self.start:.0f}')
        if self.bar is not None:
            self.bar.update(done - self.bar.n)
            return
        tenth = 10 * done // self.total
        if tenth > self.shown and done < self.total:
            self.shown = tenth
            elapsed = _time.time() - self.start
            left = elapsed * (self.total - done) / max(done, 1)
            log(f'{self.label}: {10 * tenth} % of the sweeps, '
                f'{_clock(elapsed)} so far, about {_clock(left)} left')

    def close(self, done: int):
        self.update(done)
        if self.bar is not None:
            self.bar.close()


def _clock(sec: float) -> str:
    """a duration in words: 45 s, 12 min, 2 h 05"""
    sec = int(round(sec))
    if sec < 90:
        return f'{sec} s'
    if sec < 5400:
        return f'{round(sec / 60)} min'
    return f'{sec // 3600} h {(sec % 3600) // 60:02d}'


def _run_chain(args) -> Dict[str, Any]:
    """
    Run one chain (a separate process when there are several)

    :param args: tuple, (data, grid settings, config, seed, index)

    :return: dict, the chain output
    """
    data, gridargs, cfg, seed, index = args
    with blas_threads(1):
        return _run_chain_body(data, gridargs, cfg, seed, index)


def _run_chain_body(data, gridargs, cfg, seed, index) -> Dict[str, Any]:
    """The loop of _run_chain, with the threads already limited"""
    grid = PeriodGrid(**gridargs)
    chain = _Chain(data, grid, cfg, seed)
    total = cfg['nburn'] + cfg['nsweep']
    for it in range(total):
        chain.sweep(record=it >= cfg['nburn'])
        if it < cfg['nburn'] and (it + 1) % 50 == 0:
            chain.adapt()
        if _COUNTER is not None and ((it + 1) % 5 == 0 or it + 1 == total):
            _COUNTER[index] = it + 1
    return chain.output()


# =============================================================================
# The public functions
# =============================================================================
def _default_setup(data: RVData, unit: str, seq_jitter: Optional[bool],
                   jitter: bool, tau_range: Optional[Sequence[float]],
                   ntau: int, width_range: Optional[Sequence[float]]
                   ) -> Dict[str, Any]:
    """
    Priors and starting point that follow from the scale of the data

    :return: dict, tau ladder, bounds, free flags and starting values
    """
    med_err = float(np.median(data.err))
    rstd = robust_std(data.rv)
    scale = max(rstd, med_err)
    if tau_range is None:
        tau_range = (med_err / np.sqrt(data.n), 5.0 * scale)
    tau = np.exp(np.linspace(np.log(tau_range[0]), np.log(tau_range[1]),
                             ntau))
    if width_range is None:
        width_range = (3.0 * scale, 300.0 * scale)
    has_seq = data.nseq < data.n
    if seq_jitter is None:
        seq_jitter = has_seq
    bounds = dict(jit=(1e-3 * med_err, 10 * scale),
                  sjit=(1e-3 * med_err, 10 * scale),
                  wpt=tuple(width_range), wseq=tuple(width_range))
    use_pt, use_seq = unit in ('point', 'both'), unit in ('sequence', 'both')
    free = dict(jit=bool(jitter), sjit=bool(seq_jitter and has_seq),
                wpt=use_pt, wseq=use_seq)
    # point outliers in a correlated sequence are summed configuration by
    #   configuration: a longer sequence is taken whole (koloa.noise)
    if use_pt and (free['sjit'] or use_seq):
        _long_sequences(data, use_seq)
    # the start: the robust excess as jitter, the obvious outliers flagged
    excess = np.sqrt(max(rstd ** 2 - med_err ** 2, (0.3 * med_err) ** 2))
    tnorm = (data.time - data.tref) / max(data.baseline, 1e-9)
    coeffs = np.polyfit(tnorm, data.rv, 1, w=1 / data.err)
    resid = data.rv - np.polyval(coeffs, tnorm)
    zval = np.abs(resid - np.median(resid)) / max(robust_std(resid), 1e-9)
    seqz = np.bincount(data.seq, weights=zval, minlength=data.nseq)
    seqz /= np.maximum(np.bincount(data.seq, minlength=data.nseq), 1)
    q_seq = (seqz > 5).astype(float) if use_seq else np.zeros(data.nseq)
    q_pt = np.zeros(data.n)
    if use_pt:
        q_pt = (zval > 5).astype(float)
        if use_seq:
            q_pt[q_seq[data.seq] > 0] = 0.0
    nflag = np.sum(q_pt) + np.sum(q_seq)
    init = dict(jitter=float(np.clip(excess if jitter else 0.0, 0,
                                     bounds['jit'][1])),
                seq_jitter=float(0.3 * med_err if free['sjit'] else 0.0),
                width=float(np.clip(8 * scale, *width_range)),
                frac=float(np.clip(nflag / max(data.nseq, 1), 0.01, 0.2)),
                q_pt=q_pt, q_seq=q_seq)
    if not jitter:
        init['jitter'] = 0.0
    return dict(tau=tau, bounds=bounds, free=free, init=init)


#: the series already told that their long sequences are taken whole
_LONG_TOLD = set()


def _long_sequences(data: RVData, use_seq: bool):
    """say so, once per series, when sequences are too long for their
    exposures to be outliers one by one"""
    from koloa.noise import MAX_POINT_BLOCK
    sizes = np.bincount(data.seq)
    nlong = int(np.sum(sizes > MAX_POINT_BLOCK))
    key = (data.name, data.n, data.nseq, nlong, bool(use_seq))
    if nlong and key not in _LONG_TOLD:
        _LONG_TOLD.add(key)
        log(f'{nlong} sequence{"s" if nlong > 1 else ""} of more than '
            f'{MAX_POINT_BLOCK} exposures (the longest has {np.max(sizes)}) '
            + ('taken whole: accepted or rejected as a sequence, their '
               'exposures never outliers one by one' if use_seq else
               'taken whole: their exposures are never outliers (use '
               'outliers="both" to accept or reject them as sequences)'),
            'warn')


def oafip(data: RVData, kmax: int = 3, outliers: Optional[str] = 'both',
          seq_jitter: Optional[bool] = None, jitter: bool = True,
          trend: int = 1, regressors: Optional[Dict[str, np.ndarray]] = None,
          pmin: float = 1.1, pmax: Optional[float] = None,
          oversample: int = 10, nsweep: int = 1500, nburn: int = 300,
          nchains: int = 2, ntau: int = 12,
          tau_range: Optional[Sequence[float]] = None,
          width_range: Optional[Sequence[float]] = None,
          frac_prior: Sequence[float] = (1.0, 20.0),
          prior_k: Union[str, Sequence[float]] = 'uniform', seed: int = 1,
          progress: bool = True, label: Optional[str] = None,
          npeaks: int = 5,
          freq: Optional[np.ndarray] = None,
          gp: Union[None, str, Dict[str, Any], Sequence[Any]] = None,
          nightly: Optional[bool] = None, gp_mode: str = 'auto'
          ) -> FIPResult:
    """
    The outlier-aware FIP periodogram of a series

    :param data: RVData, the series
    :param kmax: int, the largest number of signals
    :param outliers: str or None, 'both' (an exposure or a whole visit can
                     be an outlier), 'sequence' (only whole visits),
                     'point' (only single exposures), or None (gaussian)
    :param seq_jitter: bool or None, a jitter shared by the exposures of a
                       visit (on when the visits have several exposures)
    :param jitter: bool, a white jitter
    :param trend: int, the degree of the polynomial trend in time
    :param regressors: dict or None, decorrelation terms, name: values
    :param pmin: float, the shortest period [days]
    :param pmax: float or None, the longest period [days] (2T if None)
    :param oversample: int, grid points per 1/T
    :param nsweep: int, the number of recorded sweeps per chain
    :param nburn: int, the number of sweeps thrown away first
    :param nchains: int, the number of chains (run in parallel)
    :param ntau: int, the number of rungs of the amplitude prior ladder
    :param tau_range: tuple or None, the lowest and highest tau [m/s]
    :param width_range: tuple or None, the range of the outlier width [m/s]
    :param frac_prior: tuple, the beta prior (a, b) of the outlier fraction
    :param prior_k: 'uniform' or kmax + 1 prior probabilities of k
    :param seed: int, the seed of the first chain
    :param progress: bool, log the progress
    :param label: str or None, what the progress of the sweeps is shown as
                  (a bar on a terminal, a line at every tenth otherwise);
                  'OAFIP' when progress is True, none when neither
    :param npeaks: int, how many peaks to report
    :param freq: np.ndarray or None, a frequency grid to use as it is
    :param gp: None, 'local', 'rotation' (with its period), a dict or a list
               of them: a GP of the activity, as a finite basis whose
               weights are integrated out and whose hyperparameters are
               sampled with the jitters (koloa.gpbasis.setup: e.g.
               dict(kind='rotation', period=dict(mu=np.log(116), sd=0.1)),
               or dict(kind='local', length=(10, 100)), or dict(kind='sho',
               period=2.704) for a rotation known well)
    :param gp_mode: str, how the GP enters: 'kernel' (its exact kernel in
                    the covariance of the noise, dense: whatever its scales,
                    for up to a few thousand points), 'basis' (the finite
                    basis, for longer series; not for an SHO), or 'auto'
                    (the kernel up to KERNEL_MAX points)
    :param nightly: bool, run on the nightly means (RVData.nightly, koloa's
                    default): a night is then one point, and an outlier a
                    whole night; the outlier probability of each night is
                    given back to each of its exposures (outlier_prob, of
                    the length of data). False keeps the exposures, where an
                    exposure or a whole visit can be an outlier; None is
                    koloa.data.NIGHTLY (True)

    :return: FIPResult, the FIP periodogram and everything around it
    """
    start = _time.time()
    night = None
    if nightly is None:
        from koloa import data as kdata
        nightly = kdata.NIGHTLY
    if nightly:
        means = data.nightly()
        if means is not data:
            night, data = night_index(data), means
    unit = outliers if outliers in ('point', 'sequence', 'both') else 'none'
    if unit == 'both' and data.nseq == data.n:
        # no visit holds two points (nightly means): points and visits are
        #   the same
        unit = 'point'
    setup = _default_setup(data, unit, seq_jitter, jitter, tau_range, ntau,
                           width_range)
    scale = max(robust_std(data.rv), float(np.median(data.err)))
    cfg = dict(kmax=kmax, unit=unit, trend=trend, regressors=regressors,
               prior_k=prior_k, frac_prior=tuple(frac_prior), nsweep=nsweep,
               nburn=nburn, progress=progress,
               gp=gpbasis.setup(gp, data.time, scale, data.inst), **setup)
    if gp_mode == 'auto':
        gp_mode = ('kernel' if data.n <= KERNEL_MAX
                   or gpbasis.needs_kernel(cfg['gp']) else 'basis')
    if gp_mode == 'basis' and gpbasis.needs_kernel(cfg['gp']):
        raise ValueError('an SHO GP has no finite basis: gp_mode="kernel"')
    cfg['gp_mode'] = gp_mode
    gridargs = dict(time=data.time, pmin=pmin, pmax=pmax,
                    oversample=oversample, tref=data.tref, freq=freq)
    grid = PeriodGrid(**gridargs)
    method = ('gaussian' if unit == 'none' else f'mixture ({unit} outliers)')
    if night is not None:
        method += ', nightly means'
    if cfg['gp']:
        method += ', GP ' + ' + '.join(gpbasis.label(comp)
                                       for comp in cfg['gp'])
    if progress:
        log(f'OAFIP {method}: {data.n} points, {data.nseq} sequences, '
            f'{grid.size} frequencies, kmax = {kmax}, {nchains} chain(s) of '
            f'{nburn} + {nsweep} sweeps')
    jobs = [(data, gridargs, cfg, seed + 1000 * it, it)
            for it in range(nchains)]
    label = label or ('OAFIP' if progress else None)
    outs = _run_chains(jobs, label, nchains * (nburn + nsweep))
    result = _combine(outs, data, grid, cfg, method, npeaks)
    if night is not None:
        # the probability of each night, given back to its exposures
        result.settings['nightly'] = True
        result.settings['night_index'] = night
        if result.outlier_prob is not None:
            result.settings['nightly_outlier_prob'] = result.outlier_prob
            result.outlier_prob = result.outlier_prob[night]
    result.runtime = _time.time() - start
    if progress:
        log(f'OAFIP done in {result.runtime:.0f} s')
        for peak in result.peaks[:3]:
            log(f'  peak at {peak["period"]:.4f} d, FIP = {peak["fip"]:.2e}',
                'value')
    return result


def _run_chains(jobs: List[tuple], label: Optional[str], total: int
                ) -> List[Dict[str, Any]]:
    """
    The chains of a FIP, in processes of their own when there are several,
    their sweeps shown as they go when there is a label

    :return: list of dict, the output of each chain
    """
    global _COUNTER
    progress = _Progress(label, total) if label else None
    if len(jobs) == 1:
        _COUNTER = [0]
        stop = _threading.Event()

        def watch():
            while not stop.wait(1.0):
                progress.update(_COUNTER[0])
        if progress is not None:
            _threading.Thread(target=watch, daemon=True).start()
        try:
            outs = [_run_chain(jobs[0])]
        finally:
            stop.set()
            done, _COUNTER = _COUNTER[0], None
        if progress is not None:
            progress.close(done)
        return outs
    ctx = _mp.get_context('spawn')
    counter = ctx.Array('i', len(jobs), lock=False)
    with ProcessPoolExecutor(max_workers=len(jobs), mp_context=ctx,
                             initializer=_share_counter,
                             initargs=(counter,)) as pool:
        futures = [pool.submit(_run_chain, job) for job in jobs]
        pending = set(futures)
        while pending:
            _, pending = _wait(pending, timeout=1.0,
                               return_when=FIRST_COMPLETED)
            if progress is not None:
                progress.update(sum(counter))
            if CANCEL_HOOK is not None and CANCEL_HOOK():
                # the chains ended now, not when they are done
                for proc in list((getattr(pool, '_processes', None)
                                  or {}).values()):
                    proc.terminate()
                pool.shutdown(wait=False, cancel_futures=True)
                raise FIPCancelled('the FIP was stopped')
        outs = [fut.result() for fut in futures]
    if progress is not None:
        progress.close(sum(counter))
    return outs


def _combine(outs: List[Dict[str, Any]], data: RVData, grid: PeriodGrid,
             cfg: Dict[str, Any], method: str, npeaks: int) -> FIPResult:
    """
    Merge the chains into one result, with the convergence numbers

    :return: FIPResult
    """
    nacc = sum(out['nacc'] for out in outs)
    fip = sum(out['fip_acc'] for out in outs) / nacc
    density = sum(out['dens_acc'] for out in outs) / nacc
    family = sum(out['fam_acc'] for out in outs) / nacc
    chain_fip = np.array([out['fip_acc'] / out['nacc'] for out in outs])
    records = {key: np.concatenate([out['records'][key] for out in outs])
               for key in outs[0]['records']}
    kmax = cfg['kmax']
    pk = np.bincount(records['k'], minlength=kmax + 1) / len(records['k'])
    # the outlier probability of each point
    outlier_prob, expected = None, 0.0
    seq_prob = None
    if cfg['unit'] in ('point', 'sequence', 'both'):
        nq = sum(out['nq'] for out in outs)
        outlier_prob = sum(out['qprob_acc'] for out in outs) / nq
        seq_prob = sum(out['qseq_acc'] for out in outs) / nq
        expected = float(np.sum(outlier_prob))
    # convergence: split R-hat of k and of the jitter across chains
    rhat = {}
    if len(outs) > 1:
        for key in ('k', 'jitter', 'noutlier'):
            rhat[key] = _rhat([out['records'][key] for out in outs])
    settings = dict(cfg)
    settings.pop('init', None)
    settings['expected_outliers'] = expected
    settings['sequence_outlier_prob'] = seq_prob
    settings['nchains'] = len(outs)
    result = FIPResult(freq=grid.freq, fip=fip, density=density,
                       method=method, pk=pk, outlier_prob=outlier_prob,
                       unit=cfg['unit'], chains=records, settings=settings,
                       chain_fip=chain_fip, rhat=rhat, family=family)
    slots = np.vstack([out['slot_records'] for out in outs])
    result.slots = slots
    result.peaks = describe_peaks(data, grid, fip, npeaks, slots,
                                  chain_fip=chain_fip)
    for peak in result.peaks:
        peak['family_fip'] = result.family_containing(peak['period'],
                                                      1 / data.baseline)
    return result


def _rhat(chains: List[np.ndarray]) -> float:
    """
    The Gelman-Rubin statistic of a scalar across chains

    :param chains: list of np.ndarray, one trace per chain

    :return: float, R-hat (1 when converged)
    """
    length = min(len(chain) for chain in chains)
    if length < 4:
        return np.nan
    arr = np.array([np.asarray(chain[:length], dtype=float)
                    for chain in chains])
    within = np.mean(np.var(arr, axis=1, ddof=1))
    between = length * np.var(np.mean(arr, axis=1), ddof=1)
    if within <= 0:
        return 1.0 if between <= 0 else np.inf
    var = (length - 1) / length * within + between / length
    return float(np.sqrt(var / within))


def describe_peaks(data: RVData, grid: PeriodGrid, fip: np.ndarray,
                   npeaks: int = 5, slots: Optional[np.ndarray] = None,
                   chain_fip: Optional[np.ndarray] = None
                   ) -> List[Dict[str, Any]]:
    """
    The most significant intervals, with what is needed to believe them

    Next to each FIP: the window power at that period (a peak where the
    sampling has power is suspect), the FIP of its 1-day, 1-month and
    1-year aliases, and, when the sampler recorded them, the amplitude of
    the signals that fell in the interval.

    :param data: RVData, the series
    :param grid: PeriodGrid, the grid
    :param fip: np.ndarray, the FIP of each interval
    :param npeaks: int, how many peaks
    :param slots: np.ndarray or None, (grid index, K, phase) of every
                  recorded active slot
    :param chain_fip: np.ndarray or None, (nchain x G) the FIP per chain

    :return: list of dict, one per peak
    """
    score = -np.log10(np.clip(fip, 1e-300, 1))
    sep = (2 * grid.halfbin + 1) * float(np.median(np.diff(grid.freq)))
    idxs = find_peaks(grid.freq, score, npeaks, separation=sep)
    peaks = []
    for idx in idxs:
        freq0 = grid.freq[idx]
        peak = dict(index=int(idx), freq=float(freq0),
                    period=float(1.0 / freq0), fip=float(fip[idx]),
                    log10fip=float(np.log10(max(fip[idx], 1e-300))),
                    window=float(window(data.time, freq0)[0]))
        if chain_fip is not None and len(chain_fip) > 1:
            peak['fip_chains'] = [float(val) for val in chain_fip[:, idx]]
        alist = []
        for alias in aliases(freq0, grid.freq[0], grid.freq[-1]):
            aidx = int(np.argmin(np.abs(grid.freq - alias['freq'])))
            alias['fip'] = float(fip[aidx])
            alist.append(alias)
        peak['aliases'] = alist
        if slots is not None and len(slots):
            inside = np.abs(slots[:, 0] - idx) <= grid.halfbin
            if np.sum(inside) > 5:
                amp = slots[inside, 1]
                per = grid.period[slots[inside, 0].astype(int)]
                peak['amplitude'] = tuple(float(val) for val in
                                          np.percentile(amp, [50, 16, 84]))
                peak['period_post'] = tuple(float(val) for val in
                                            np.percentile(per, [50, 16, 84]))
        peaks.append(peak)
    return peaks


# =============================================================================
# The single-signal FIP with fixed noise (Hara-comparable, and fast)
# =============================================================================
def fip_single(data: RVData, err: Optional[np.ndarray] = None,
               jitter: Union[str, float] = 'robust', tau: Any = 'ladder',
               trend: int = 1, pmin: float = 1.1,
               pmax: Optional[float] = None, oversample: int = 10,
               prior_signal: float = 0.5, ntau: int = 12,
               freq: Optional[np.ndarray] = None, npeaks: int = 5,
               label: str = 'gaussian, one signal') -> FIPResult:
    """
    The FIP of one circular signal, with the noise held fixed: no sampling

    This is the calculation of LBL's report (lbl.science.report,
    fip_periodogram) in koloa's machinery: the jitter is the dispersion the
    error bars do not explain, the linear parameters are integrated out,
    and the only sum left, over the grid and the tau ladder, is done
    exactly. It is fast and comparable with a single-planet FIP of Hara et
    al. (2022). It is what the outlier-aware result is compared to.

    :param data: RVData, the series
    :param err: np.ndarray or None, error bars to use instead of the data's
                (for instance inflated by a soft clip)
    :param jitter: 'robust' (the excess of the robust rms over the median
                   error), or a value [m/s]
    :param tau: 'ladder' (the same log-uniform ladder as oafip), or a fixed
                prior width of the amplitudes [m/s]; LBL used ten times the
                rms of the data
    :param trend: int, the degree of the polynomial trend
    :param pmin: float, the shortest period [days]
    :param pmax: float or None, the longest period [days]
    :param oversample: int, grid points per 1/T
    :param prior_signal: float, the prior probability that a signal exists
    :param ntau: int, the rungs of the ladder
    :param freq: np.ndarray or None, a grid to use as it is
    :param npeaks: int, how many peaks to report
    :param label: str, the name of the method in the result

    :return: FIPResult
    """
    start = _time.time()
    err = data.err if err is None else np.asarray(err, dtype=float)
    if jitter == 'robust':
        excess = robust_std(data.rv) ** 2 - np.median(err) ** 2
        jitter = float(np.sqrt(excess)) if excess > 0 else 0.0
    grid = PeriodGrid(data.time, pmin, pmax, oversample, data.tref, freq)
    trig_cs, trig_sq = grid.trig()
    size = grid.size
    cov = BlockCov(err ** 2 + jitter ** 2)
    base = base_design(data, trend)
    model = LinearModel(base['design'], data.rv, cov, base['prior_var'])
    if isinstance(tau, str):
        setup = _default_setup(data, 'none', False, True, None, ntau, None)
        tau2 = setup['tau'] ** 2
    else:
        tau2 = np.atleast_1d(float(tau) ** 2)
    prod = np.vstack([model.wdesign.T, model.wvalue[None, :]]) @ trig_cs
    nbase = model.npar
    gain = sinusoid_gain(model, prod[:nbase, :size], prod[:nbase, size:],
                         prod[nbase, :size], prod[nbase, size:],
                         *cov.grid_quad(trig_sq, trig_cs, size), tau2)
    logw = (np.log(prior_signal) + grid.logprior[:, None]
            - np.log(len(tau2)) + gain)
    log_off = np.log(1 - prior_signal)
    lmax = max(float(np.max(logw)), log_off)
    pon = np.sum(np.exp(logw - lmax), axis=1)
    woff = float(np.exp(log_off - lmax))
    total = woff + float(np.sum(pon))
    fip = grid.bin_fip(pon, woff, total)
    result = FIPResult(freq=grid.freq, fip=fip, density=pon / total,
                       method=label,
                       pk=np.array([woff / total, 1 - woff / total]),
                       settings=dict(jitter=jitter, trend=trend,
                                     tau=np.sqrt(tau2)))
    result.peaks = describe_peaks(data, grid, fip, npeaks)
    result.runtime = _time.time() - start
    return result


def fip_comparison(data: RVData, clip: float = 3.0, **kwargs
                   ) -> Dict[str, FIPResult]:
    """
    The single-signal FIP without cleaning, with a soft clip, with a hard
    clip: what each way of handling the outliers buys

    The point of reporting them side by side (and next to the OAFIP) is
    that the reader sees what the cleaning did to the answer, rather than
    only the cleaned answer.

    :param data: RVData, the series
    :param clip: float, the clip threshold [sigma]
    :param kwargs: passed to fip_single

    :return: dict, gaussian, soft and hard results
    """
    out = dict()
    out['gaussian'] = fip_single(data, label='gaussian, no clip', **kwargs)
    soft = soft_clip(data.time, data.rv, data.err, clip=clip, quiet=True)
    res = fip_single(data, err=soft['err'],
                     label=f'gaussian, {clip:g}-sigma soft clip', **kwargs)
    res.settings.update(removed=soft['removed'], ndown=soft['ndown'],
                        over_budget=soft['over_budget'])
    out['soft'] = res
    keep = hard_clip(data.time, data.rv, data.err, clip=clip)
    grid_freq = out['gaussian'].freq
    sub = data.select(keep)
    res = fip_single(sub, label=f'gaussian, {clip:g}-sigma hard clip',
                     freq=grid_freq, **{key: val for key, val in
                                        kwargs.items() if key != 'freq'})
    res.settings.update(nremoved=int(np.sum(~keep)))
    out['hard'] = res
    return out

def inflate_to_fit(fit: Any):
    """
    The series with its errors inflated instrument by instrument, so that
    the FIP sees the noise that a fit with a jitter per instrument found

    The FIP sampler has one white jitter and one visit jitter for every
    instrument, whereas instruments differ: the optical ones see more of a
    star's activity than the near-infrared ones, and an old instrument is
    noisier than a new one. The smallest visit jitter of the fit, s_ref, is
    left to the sampler; the errors of every other instrument grow so that
    the variance of its visits matches the fit: by n (s^2 - s_ref^2) + j^2
    for an instrument with its own visit jitter s, n exposures per visit
    and white jitter j, and by max(j^2 - s_ref^2, 0) for one without.
    Without any visit jitter (nightly means, one exposure per visit), the
    smallest white jitter j_ref is the one left to the sampler, and the
    errors grow by max(j^2 - j_ref^2, 0): one instrument is not inflated at
    all, and its jitter is sampled with the signals (a fit without planets
    would otherwise put a planet's variance in the errors).

    :param fit: FitResult, a fit of the series with a jitter per instrument
                (RVModel(..., seq_jitter='instrument')), with or without
                planets

    :return: tuple, the RVData with the inflated errors, and a dict with the
             white and visit jitters of the fit, s_ref and the error added
             in quadrature to each instrument [m/s]
    """
    model, theta = fit.model, fit.theta
    data = model.data

    def value(pname):
        """a jitter of the fit, or zero if the model has none"""
        return (float(np.exp(theta[model.index[pname]]))
                if pname in model.index else 0.0)
    insts = [str(inst) for inst in data.instruments]
    jit = {inst: value(f'log_jit_{inst}') for inst in insts}
    sjit = {inst: value(f'log_sjit_{inst}') for inst in insts
            if f'log_sjit_{inst}' in model.index}
    sref = min(sjit.values()) if sjit else 0.0
    jref = min(jit.values()) if (jit and not sjit) else 0.0
    added = {}
    for inst in insts:
        if not sjit:
            added[inst] = max(jit[inst] ** 2 - jref ** 2, 0.0)
        elif inst in sjit:
            # the median number of exposures per visit of this instrument
            seqs = data.seq[data.inst == inst]
            nexp = float(np.median(np.bincount(seqs)[np.unique(seqs)]))
            added[inst] = nexp * (sjit[inst] ** 2 - sref ** 2) + jit[inst] ** 2
        else:
            added[inst] = max(jit[inst] ** 2 - sref ** 2, 0.0)
    extra = np.array([added[str(inst)] for inst in data.inst])
    # a copy of the series: the velocities as given (RVData takes the zero
    #   point of each instrument out again), the errors inflated
    zero = np.array([data.zero_point[str(inst)] for inst in data.inst])
    out = RVData(data.time.copy(), data.rv + zero,
                 np.sqrt(data.err ** 2 + extra), inst=data.inst.copy(),
                 seq=data.seq.copy(),
                 indicators={key: (val[0].copy(), val[1].copy())
                             for key, val in data.indicators.items()},
                 name=data.name, zero_point=dict(data.zero_point),
                 sequence_gap=data.sequence_gap,
                 meta={key: val.copy() for key, val in data.meta.items()})
    info = dict(jitters=jit, visit_jitters=sjit, visit_jitter_ref=sref,
                inflation={inst: float(np.sqrt(val))
                           for inst, val in added.items()})
    return out, info


# =============================================================================
# End of code
# =============================================================================

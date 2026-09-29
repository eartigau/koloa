#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Keplerians, a gaussian process and outliers, fitted together.

One model, `RVModel`, holds every part of what a velocity series is made
of: an offset and a jitter per instrument, a trend, decorrelation terms,
any number of orbits (circular or eccentric), optionally a GP for the
activity, a jitter shared by the exposures of a sequence, and a noise
model that says what an outlier is:

- 'gaussian': no outliers (the reference, and what a GP-only fit uses);
- 'student': a Student-t on each point, heavy tails with no latent state;
- 'mixture': a unit (point or sequence) is an outlier with probability f,
  and then its variance grows by W^2 (Box & Tiao 1968; Hogg et al. 2010).

Without a GP the units are independent and the mixture is summed exactly
in the likelihood, so any sampler (emcee here) works on it. With a GP the
units are coupled: the maximum is found by alternating the parameters and
the indicators, and the posterior is sampled by Metropolis within Gibbs,
the indicators drawn from their exact conditional (koloa.gp.loo_blocks).

A fit is joint by construction. `sequential_fit` exists only to show why it
should not be used: fitting the GP first and an orbit to what is left lets
the GP take the signal, and the function says so loudly.

Created on 2026-09-27

@author: artigau
"""
import time as _time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, gammaln, logit

from koloa import gp as kgp
from koloa import kepler
from koloa.data import RVData, robust_std
from koloa.log import log, loud_warning
from koloa.noise import block_gauss_loglike, mixture_loglike
from koloa.noise import mixture_loglike_both
from koloa.periodogram import find_peaks, frequency_grid, gls
from koloa.utils import blas_threads

# =============================================================================
# Define variables
# =============================================================================
#: the largest eccentricity an orbit may take
EMAX = 0.95


# =============================================================================
# Priors
# =============================================================================
@dataclass
class Prior:
    """
    The prior of one sampled parameter

    kind is 'uniform' (a, b are the bounds), 'normal' (a is the mean, b the
    sigma), or 'beta_logit' (the parameter is logit(f) and f has a beta
    prior of shape a, b; with the jacobian f (1 - f) the density of the
    logit is f^a (1 - f)^b).
    """
    kind: str
    a: float
    b: float

    def logpdf(self, value: float) -> float:
        """The log density, up to a constant"""
        if self.kind == 'uniform':
            return 0.0 if self.a <= value <= self.b else -np.inf
        if self.kind == 'normal':
            return -0.5 * ((value - self.a) / self.b) ** 2
        if self.kind == 'beta_logit':
            # p(u) = Beta(f) df/du = f^a (1 - f)^b up to a constant
            return float(-self.a * np.logaddexp(0, -value)
                         - self.b * np.logaddexp(0, value))
        raise ValueError(f'Unknown prior: {self.kind}')

    def curvature(self) -> float:
        """
        The inverse variance of the prior (a uniform one counts as a
        gaussian of the same variance), which keeps a Laplace covariance
        finite along a direction the data do not constrain
        """
        if self.kind == 'uniform':
            return 12.0 / (self.b - self.a) ** 2
        if self.kind == 'normal':
            return 0.0
        return 1.0 / (np.pi ** 2 / 3 + 1.0)

    def bounds(self) -> Tuple[Optional[float], Optional[float]]:
        """Hard bounds, for the optimiser"""
        if self.kind == 'uniform':
            return self.a, self.b
        return None, None

    def draw(self, rng: np.random.Generator) -> float:
        """One draw (used to start walkers far from a maximum)"""
        if self.kind == 'uniform':
            return float(rng.uniform(self.a, self.b))
        if self.kind == 'normal':
            return float(rng.normal(self.a, self.b))
        return float(logit(rng.beta(self.a, self.b)))


# =============================================================================
# The model
# =============================================================================
class RVModel:
    """
    Everything a velocity series is made of, with its parameters and priors
    """

    def __init__(self, data: RVData, planets: Sequence[Dict[str, Any]] = (),
                 gp: Optional[Dict[str, Any]] = None,
                 likelihood: str = 'mixture', unit: str = 'both',
                 jitter: bool = True,
                 seq_jitter: Optional[Union[bool, str]] = None,
                 trend: int = 1,
                 regressors: Optional[Dict[str, np.ndarray]] = None,
                 dof: float = 4.0, frac_prior: Sequence[float] = (1.0, 20.0),
                 kmax: Optional[float] = None, emax: float = EMAX,
                 priors: Optional[Dict[str, Prior]] = None,
                 perspective: Optional[Sequence[Optional[float]]] = None,
                 secular: Optional[Sequence[Optional[float]]] = None):
        """
        :param data: RVData, the series
        :param planets: list of dict, one per orbit: period (the start
                        [days]), eccentric (bool, default False),
                        period_range (tuple [days], default +-10 %); for a
                        transiting planet, tc [days], tc_err and
                        period_err turn the period and the phase into
                        gaussian priors from the transits
        :param gp: dict or None, kernel (sho, rotation, qp, se, matern32),
                   and optionally prior ({name: Prior or (mean, sigma) in
                   the log}) and init ({name: value}); scale='instrument'
                   makes it chromatic: the kernel's amplitude is that of
                   the instrument reference (the first by default), and
                   each other instrument has its own relative to it
                   (gp_log_scale_<inst>)
        :param likelihood: str, gaussian, student or mixture
        :param unit: str, what an outlier is: point, sequence, or both
                     (a lone exposure or a whole visit)
        :param jitter: bool, a white jitter per instrument
        :param seq_jitter: bool, 'instrument' or None, a jitter per sequence
                           (on when the sequences have several exposures);
                           'instrument' gives each instrument whose visits
                           hold several exposures its own (sjit_<inst>),
                           for series whose excess noise differs, the
                           others keeping their white jitter
        :param trend: int, the degree of the polynomial trend: its first
                      two derivatives are the acceleration of the star and
                      its change, which koloa.secular.acceleration reads
                      in m/s/yr and m/s/yr^2 with their errors
        :param regressors: dict or None, decorrelation terms, name: values
        :param dof: float, the degrees of freedom of the Student-t
        :param frac_prior: tuple, the beta prior of the outlier fraction
        :param kmax: float or None, the largest semi-amplitude [m/s]
        :param emax: float, the largest eccentricity
        :param priors: dict or None, priors that replace the defaults
        :param perspective: tuple or None, the perspective acceleration of
                            the star, mu^2 d, as (value, error) [m/s/yr]
                            (koloa.secular gives it from the astrometry): a
                            parameter secacc with that gaussian prior, times
                            (t - tref); with an error of None it is free
                            (uniform), a drift in physical units. A trend
                            fitted with it takes the acceleration of the
                            star itself. Leave it None for velocities whose
                            barycentric correction already removed it
                            (APERO's, through barycorrpy).
        :param secular: tuple or None, the name of perspective in the first
                        version (used when perspective is None)
        """
        if likelihood not in ('gaussian', 'student', 'mixture'):
            raise ValueError(f'Unknown likelihood: {likelihood}')
        if unit not in ('point', 'sequence', 'both'):
            raise ValueError(f'Unknown outlier unit: {unit}')
        if unit == 'both' and data.nseq == data.n:
            # no visit holds two exposures: points and visits are the same
            unit = 'point'
        self.data = data
        self.likelihood, self.unit = likelihood, unit
        if gp is None:
            self.gp = None
        elif isinstance(gp, dict):
            self.gp = dict(gp)
        else:
            # a GP per group of instruments, each its own parameters
            self.gp = dict(kernel='multi', groups=[dict(gg) for gg in gp])
        self.dof = dof
        self.emax = emax
        self.trend = trend
        self.regressors = {key: np.asarray(val, dtype=float)
                           for key, val in (regressors or {}).items()}
        self.planets = [dict(pl) for pl in planets]
        has_seq = data.nseq < data.n
        if isinstance(seq_jitter, str) and seq_jitter != 'instrument':
            raise ValueError(f'Unknown seq_jitter: {seq_jitter}')
        self.use_seq_jitter = has_seq if seq_jitter is None else (
            bool(seq_jitter) and has_seq)
        # the instrument of each sequence, and the instruments with their
        #   own sequence jitter (those whose visits hold several exposures)
        self.seq_inst = np.zeros(data.nseq, dtype=int)
        self.seq_inst[data.seq] = data.inst_index
        nper = np.bincount(data.seq, minlength=data.nseq)
        self.seq_jitter_insts = (
            [data.instruments[it] for it in np.unique(self.seq_inst[nper > 1])]
            if seq_jitter == 'instrument' and self.use_seq_jitter else [])
        if likelihood == 'student' and self.use_seq_jitter:
            raise ValueError('The Student-t likelihood treats points one by '
                             'one: set seq_jitter=False')
        scale = max(robust_std(data.rv), float(np.median(data.err)))
        self.scale = scale
        self.kmax = kmax if kmax is not None else 10 * scale + 10
        self.tref = data.tref
        self.inst_index = data.inst_index
        self.ninst = len(data.instruments)
        self.tnorm = (data.time - self.tref) / max(data.baseline, 1e-9)
        self.units = (np.arange(data.nseq) if unit == 'sequence'
                      else np.arange(data.n))
        if unit == 'both' and data.nseq < data.n and likelihood == 'mixture':
            # a visit of more than 10 exposures is taken whole: accepted or
            #   rejected as a visit (koloa.noise.MAX_POINT_BLOCK)
            from koloa.fip import _long_sequences
            _long_sequences(data, True)
        # ---------------------------------------------------------------------
        # the parameters, in order, with their priors and starting values
        # ---------------------------------------------------------------------
        self.names, self.priors, self.init = [], [], []
        wide = 100 * scale
        for it, inst in enumerate(data.instruments):
            mask = data.inst == inst
            self._add(f'offset_{inst}', Prior('uniform', -wide, wide),
                      float(np.median(data.rv[mask])))
        for deg in range(1, trend + 1):
            self._add(f'trend_{deg}', Prior('uniform', -wide, wide), 0.0)
        if perspective is None:
            perspective = secular
        self.secular = None if perspective is None else tuple(perspective)
        if self.secular is not None:
            sec_value, sec_error = self.secular
            if sec_error is None or not np.isfinite(sec_error) or \
                    sec_error <= 0:
                # free: as wide as a drift of 100 scales over the baseline
                span = wide * 365.25 / max(data.baseline, 1.0)
                self._add('secacc', Prior('uniform', -span, span),
                          float(sec_value or 0.0))
            else:
                self._add('secacc', Prior('normal', float(sec_value),
                                          float(sec_error)),
                          float(sec_value))
        for key in self.regressors:
            self._add(f'coef_{key}', Prior('uniform', -wide, wide), 0.0)
        for ip, planet in enumerate(self.planets):
            period = float(planet['period'])
            prange = planet.get('period_range', (period * 0.9, period * 1.1))
            planet['period_range'] = prange
            if 'tc' in planet:
                # a transiting planet: period and conjunction from the
                #   transits, as gaussian priors, and the amplitude free
                perr = float(planet.get('period_err', 1e-4 * period))
                self._add(f'logP_{ip}', Prior('normal', np.log(period),
                                              perr / period),
                          float(np.log(period)))
                tcerr = float(planet.get('tc_err', 0.01))
                self._add(f'tc_{ip}', Prior('normal', float(planet['tc']),
                                            tcerr), float(planet['tc']))
                self._add(f'k_{ip}', Prior('uniform', 0.0, self.kmax), 1.0)
                if planet.get('eccentric', False):
                    sqe = np.sqrt(emax)
                    self._add(f'xe_{ip}', Prior('uniform', -sqe, sqe), 0.1)
                    self._add(f'ye_{ip}', Prior('uniform', -sqe, sqe), 0.1)
                continue
            self._add(f'logP_{ip}', Prior('uniform', np.log(prange[0]),
                                          np.log(prange[1])),
                      float(np.log(period)))
            sqk = np.sqrt(self.kmax)
            self._add(f'xk_{ip}', Prior('uniform', -sqk, sqk), 1.0)
            self._add(f'yk_{ip}', Prior('uniform', -sqk, sqk), 0.0)
            if planet.get('eccentric', False):
                sqe = np.sqrt(emax)
                self._add(f'xe_{ip}', Prior('uniform', -sqe, sqe), 0.1)
                self._add(f'ye_{ip}', Prior('uniform', -sqe, sqe), 0.1)
        med_err = float(np.median(data.err))
        excess = np.sqrt(max(robust_std(data.rv) ** 2 - med_err ** 2,
                             (0.3 * med_err) ** 2))
        if jitter:
            for inst in data.instruments:
                self._add(f'log_jit_{inst}',
                          Prior('uniform', np.log(1e-3 * med_err),
                                np.log(10 * scale)), float(np.log(excess)))
        if self.use_seq_jitter:
            for name in ([f'log_sjit_{inst}' for inst in self.seq_jitter_insts]
                         or ['log_sjit']):
                self._add(name, Prior('uniform', np.log(1e-3 * med_err),
                                      np.log(10 * scale)),
                          float(np.log(0.3 * med_err)))
        if likelihood == 'mixture':
            self._add('logit_f', Prior('beta_logit', frac_prior[0],
                                       frac_prior[1]),
                      float(logit(0.03)))
            self._add('log_W', Prior('uniform', np.log(3 * scale),
                                     np.log(300 * scale)),
                      float(np.log(8 * scale)))
            if unit == 'both':
                # logit_f and log_W are the point outliers, these the
                #   sequence ones
                self._add('logit_fs', Prior('beta_logit', frac_prior[0],
                                            frac_prior[1]),
                          float(logit(0.03)))
                self._add('log_Ws', Prior('uniform', np.log(3 * scale),
                                          np.log(300 * scale)),
                          float(np.log(8 * scale)))
        if self.gp is not None:
            self._setup_gp(scale)
        for name, prior in (priors or {}).items():
            self.priors[self.names.index(name)] = prior
        self.init = np.array(self.init, dtype=float)
        self.index = {name: it for it, name in enumerate(self.names)}
        # the latent outlier indicators, used only with a GP and a mixture:
        #   qunit for the units (the visits when unit is 'both'), qpoint
        #   for the exposures of unit='both'
        self.qunit = np.zeros(data.nseq if unit in ('sequence', 'both')
                              else data.n)
        self.qpoint = np.zeros(data.n)
        self._smart_start()

    # -------------------------------------------------------------------------
    def _add(self, name: str, prior: Prior, init: float):
        """Register one parameter"""
        self.names.append(name)
        self.priors.append(prior)
        self.init.append(init)

    @staticmethod
    def _gp_defaults(scale: float, baseline: float) -> Dict[str, Any]:
        """The priors and starts of the kernel parameters, for a series of
        this scale [m/s] and baseline [days]"""
        return {
            'log_sigma': (Prior('uniform', np.log(1e-2 * scale),
                                np.log(20 * scale)), np.log(scale)),
            'log_period': (Prior('uniform', np.log(1.0), np.log(baseline)),
                           np.log(min(20.0, baseline / 4))),
            'log_quality': (Prior('uniform', np.log(0.3), np.log(1e3)),
                            np.log(2.0)),
            'log_q0': (Prior('uniform', np.log(0.05), np.log(50.0)),
                       np.log(1.0)),
            'log_dq': (Prior('uniform', np.log(0.01), np.log(100.0)),
                       np.log(1.0)),
            'logit_mix': (Prior('normal', 0.0, 2.0), 0.0),
            'log_decay': (Prior('uniform', np.log(1.0), np.log(3 * baseline)),
                          np.log(40.0)),
            'log_smooth': (Prior('uniform', np.log(0.1), np.log(3.0)),
                           np.log(0.5)),
            'log_length': (Prior('uniform', np.log(0.3), np.log(baseline)),
                           np.log(10.0)),
        }

    def _add_gp_params(self, kernel: str, prefix: str, defaults: Dict,
                       user_prior: Dict, user_init: Dict) -> List[str]:
        """The parameters of one kernel, named prefix + name"""
        names = []
        for pname in kgp.KERNEL_PARAMS[kernel]:
            prior, init = defaults[pname]
            if pname in user_prior:
                val = user_prior[pname]
                prior = val if isinstance(val, Prior) else Prior(
                    'normal', float(val[0]), float(val[1]))
                if prior.kind == 'normal':
                    init = prior.a
            if pname in user_init:
                init = float(user_init[pname])
            self._add(prefix + pname, prior, float(init))
            names.append(prefix + pname)
        return names

    def _setup_gp_groups(self):
        """One GP per group of instruments, each with its own kernel and
        parameters (gp_<name>_<parameter>): the covariance is block
        diagonal, the instruments of two groups uncorrelated"""
        self.gp_groups, taken = [], set()
        inst = self.data.inst.astype(str)
        # the velocities about each instrument's median, for the scale
        centred = self.data.rv.astype(float).copy()
        for name in self.data.instruments:
            centred[inst == name] -= np.median(centred[inst == name])
        for ig, group in enumerate(self.gp['groups']):
            insts = [str(ii) for ii in group.get('instruments', [])]
            unknown = [ii for ii in insts if ii not in self.data.instruments]
            if not insts or unknown:
                raise ValueError(f'GP group {ig}: no instruments, or unknown '
                                 f'ones: {unknown}')
            if taken & set(insts):
                raise ValueError(f'GP group {ig}: an instrument is already '
                                 f'in another group')
            taken |= set(insts)
            name = str(group.get('name', '_'.join(insts))).replace(' ', '')
            kernel = group.get('kernel', 'matern32')
            group['kernel'] = kernel
            idx = np.where(np.isin(inst, insts))[0]
            gscale = max(robust_std(centred[idx]),
                         float(np.median(self.data.err[idx])))
            baseline = float(np.ptp(self.data.time[idx])) or \
                self.data.baseline
            params = self._add_gp_params(
                kernel, f'gp_{name}_', self._gp_defaults(gscale, baseline),
                group.get('prior', {}) or {}, group.get('init', {}) or {})
            self.gp_groups.append(dict(name=name, kernel=kernel, idx=idx,
                                       params=params, instruments=insts))

    def _setup_gp(self, scale: float):
        """The parameters of the GP kernel, their priors and starts"""
        self.gp_scale_insts = []
        if self.gp.get('kernel') == 'multi':
            self._setup_gp_groups()
            return
        kernel = self.gp.get('kernel', 'sho')
        self.gp['kernel'] = kernel
        user_prior = self.gp.get('prior', {}) or {}
        user_init = self.gp.get('init', {}) or {}
        self._add_gp_params(kernel, 'gp_',
                            self._gp_defaults(scale, self.data.baseline),
                            user_prior, user_init)
        # a chromatic GP: the kernel's amplitude is that of a reference
        #   instrument, and each other instrument has its own relative to it
        #   (activity is weaker in the near infrared than in the optical)
        self.gp_scale_insts = []
        if self.gp.get('scale') == 'instrument' and \
                len(self.data.instruments) > 1:
            ref = self.gp.get('reference', self.data.instruments[0])
            if ref not in self.data.instruments:
                raise ValueError(f'Unknown GP reference instrument: {ref}')
            self.gp['reference'] = ref
            self.gp_scale_insts = [inst for inst in self.data.instruments
                                   if inst != ref]
            for inst in self.gp_scale_insts:
                pname = f'log_scale_{inst}'
                prior = Prior('uniform', np.log(0.01), np.log(100.0))
                if pname in user_prior:
                    val = user_prior[pname]
                    prior = val if isinstance(val, Prior) else Prior(
                        'normal', float(val[0]), float(val[1]))
                self._add(f'gp_{pname}', prior,
                          float(user_init.get(pname, 0.0)))
        elif self.gp.get('scale') not in (None, 'instrument'):
            raise ValueError(f'Unknown GP scale: {self.gp["scale"]}')

    def _linear_fill(self, theta: np.ndarray,
                     mask: Optional[np.ndarray] = None
                     ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Offsets, trend, regressors, and the amplitude and phase of each
        orbit, from a linear fit at the periods of theta

        :param theta: np.ndarray, the parameters (the periods are kept)
        :param mask: np.ndarray or None, the points to fit (all when None)

        :return: tuple, the filled parameters, and the residual of every
                 point
        """
        theta = np.array(theta, dtype=float)
        design, cols = self._linear_design(theta)
        err = np.sqrt(self.data.err ** 2 + self.scale ** 2)
        target = self.data.rv.copy()
        if 'secacc' in self.index:
            target -= theta[self.index['secacc']] * \
                (self.data.time - self.tref) / 365.25
        sel = np.ones(self.data.n, dtype=bool) if mask is None else mask
        coeffs = np.linalg.lstsq(design[sel] / err[sel, None],
                                 target[sel] / err[sel], rcond=None)[0]
        for name, val in zip(cols, coeffs):
            if name.startswith('cos_') or name.startswith('sin_'):
                continue
            theta[self.index[name]] = val
        for ip in range(len(self.planets)):
            aval = coeffs[cols.index(f'cos_{ip}')]
            bval = coeffs[cols.index(f'sin_{ip}')]
            amp = max(np.hypot(aval, bval), 1e-3)
            if f'k_{ip}' in self.index:
                theta[self.index[f'k_{ip}']] = amp
                continue
            # A cos(x) + B sin(x) = K cos(x + lambda0) with lambda0 = atan2(-B, A)
            lam = np.arctan2(-bval, aval)
            amp = min(amp, 0.99 * self.kmax)
            theta[self.index[f'xk_{ip}']] = np.sqrt(amp) * np.cos(lam)
            theta[self.index[f'yk_{ip}']] = np.sqrt(amp) * np.sin(lam)
        return theta, target - design @ coeffs

    def _smart_start(self):
        """Offsets, amplitudes and phases from a linear fit at the starts"""
        self.init, _ = self._linear_fill(self.init)

    def _trimmed_starts(self, start: np.ndarray,
                        keeps: Sequence[float] = (0.7, 0.5),
                        niter: int = 6) -> List[np.ndarray]:
        """
        Starts that no outlier can drag, however many there are

        The robust start of polyband (Artigau 2026): the linear part is
        fitted on the least deviant fraction of the units only, re-selected
        a few times; the jitter comes from the dispersion of the units kept
        (corrected for the truncation), and the outlier fraction and width
        from the units left out. A start from no outliers stays there when
        outliers are many (40 % and more): the likelihood then compares
        these starts with the others.

        :param start: np.ndarray, the parameters to start from (periods)
        :param keeps: list of float, the fractions of units kept
        :param niter: int, the re-selections

        :return: list of np.ndarray, one start per fraction
        """
        from scipy.stats import norm
        if self.likelihood != 'mixture' or self.gp is not None:
            return []
        data = self.data
        # judge the units that are outliers: whole visits, or exposures
        if self.unit == 'sequence':
            members = self._unit_members()
        else:
            members = [np.array([it]) for it in range(data.n)]
        nunit = len(members)
        med_err = float(np.median(data.err))
        starts = []
        for keep in keeps:
            nkeep = int(np.ceil(keep * nunit))
            if nkeep < len(self._linear_design(start)[1]) + 3:
                continue
            theta, resid = self._linear_fill(start)
            mask = np.ones(data.n, dtype=bool)
            for _ in range(niter):
                score = np.array([abs(np.mean(resid[mm])) for mm in members])
                new = np.zeros(data.n, dtype=bool)
                for unit in np.argsort(score, kind='stable')[:nkeep]:
                    new[members[unit]] = True
                if np.array_equal(new, mask):
                    break
                mask = new
                theta, resid = self._linear_fill(theta, mask)
            # the kept points are a gaussian truncated at the keep quantile
            quant = norm.ppf(0.5 + 0.5 * keep)
            shrink = np.sqrt(1 - 2 * quant * norm.pdf(quant) / keep)
            sigma = float(np.std(resid[mask])) / max(shrink, 1e-3)
            jit = np.sqrt(max(sigma ** 2 - med_err ** 2, (0.1 * med_err) ** 2))
            for inst in data.instruments:
                name = f'log_jit_{inst}'
                if name in self.index:
                    low, high = self.priors[self.index[name]].bounds()
                    theta[self.index[name]] = float(np.clip(
                        np.log(jit), low, high))
            # what was left out, as a fraction and a width
            score = np.array([abs(np.mean(resid[mm])) for mm in members])
            size = np.array([len(mm) for mm in members])
            bad = score > 4 * sigma / np.sqrt(size)
            frac = float(np.clip(np.mean(bad), 0.02, 0.6))
            width = float(np.std(resid[np.concatenate(
                [members[it] for it in np.where(bad)[0]])])) if np.any(bad) \
                else 8 * sigma
            for fname, wname in (('logit_f', 'log_W'), ('logit_fs', 'log_Ws')):
                if fname in self.index:
                    theta[self.index[fname]] = float(logit(frac))
                    low, high = self.priors[self.index[wname]].bounds()
                    theta[self.index[wname]] = float(np.clip(
                        np.log(max(width, 1e-3)), low + 1e-6, high - 1e-6))
            if np.isfinite(self.log_posterior(theta)):
                starts.append(theta)
        return starts

    def _linear_design(self, theta: np.ndarray) -> Tuple[np.ndarray, List[str]]:
        """Offsets, trend, regressors and a sinusoid per planet, as columns"""
        cols, names = [], []
        for it, inst in enumerate(self.data.instruments):
            cols.append((self.inst_index == it).astype(float))
            names.append(f'offset_{inst}')
        for deg in range(1, self.trend + 1):
            cols.append(self.tnorm ** deg)
            names.append(f'trend_{deg}')
        for key, val in self.regressors.items():
            cols.append(val - np.mean(val))
            names.append(f'coef_{key}')
        for ip in range(len(self.planets)):
            period = np.exp(theta[self.index[f'logP_{ip}']])
            phase = 2 * np.pi * (self.data.time - self.tref) / period
            cols += [np.cos(phase), np.sin(phase)]
            names += [f'cos_{ip}', f'sin_{ip}']
        return np.array(cols).T, names

    # -------------------------------------------------------------------------
    @property
    def ndim(self) -> int:
        """The number of parameters"""
        return len(self.names)

    def value(self, theta: np.ndarray, name: str, default: float = 0.0
              ) -> float:
        """One parameter by name"""
        idx = self.index.get(name)
        return default if idx is None else float(theta[idx])

    def orbit(self, theta: np.ndarray, ip: int
              ) -> Tuple[float, float, float, float, float]:
        """
        One orbit as (P, tp, e, omega, K)

        A circular orbit has e = 0 exactly: it has no eccentricity
        parameter, so no optimiser can move it.
        """
        logp = theta[self.index[f'logP_{ip}']]
        if f'tc_{ip}' in self.index:
            period = float(np.exp(logp))
            tconj = float(theta[self.index[f'tc_{ip}']])
            amp = float(theta[self.index[f'k_{ip}']])
            ecc, omega = 0.0, 0.0
            if f'xe_{ip}' in self.index:
                xe = theta[self.index[f'xe_{ip}']]
                ye = theta[self.index[f'ye_{ip}']]
                ecc = float(xe ** 2 + ye ** 2)
                omega = float(np.arctan2(ye, xe)) if ecc > 0 else 0.0
            tperi = kepler.tc_to_tp(tconj, period, ecc, omega)
            return period, tperi, ecc, omega, amp
        xk = theta[self.index[f'xk_{ip}']]
        yk = theta[self.index[f'yk_{ip}']]
        if f'xe_{ip}' in self.index:
            xe = theta[self.index[f'xe_{ip}']]
            ye = theta[self.index[f'ye_{ip}']]
        else:
            xe = ye = 0.0
        return kepler.basis_to_orbit(logp, xk, yk, xe, ye, self.tref)

    def planet_rv(self, theta: np.ndarray, ip: int,
                  time: Optional[np.ndarray] = None) -> np.ndarray:
        """The velocity of one orbit"""
        time = self.data.time if time is None else np.asarray(time)
        period, tperi, ecc, omega, amp = self.orbit(theta, ip)
        return kepler.rv_keplerian(time, period, tperi, ecc, omega, amp)

    def systematics(self, theta: np.ndarray) -> np.ndarray:
        """Offsets, trend and decorrelation at the data"""
        out = np.zeros(self.data.n)
        for it, inst in enumerate(self.data.instruments):
            out[self.inst_index == it] += theta[self.index[f'offset_{inst}']]
        for deg in range(1, self.trend + 1):
            out += theta[self.index[f'trend_{deg}']] * self.tnorm ** deg
        if 'secacc' in self.index:
            # m/s/yr times years from the reference time
            out += theta[self.index['secacc']] * \
                (self.data.time - self.tref) / 365.25
        for key, val in self.regressors.items():
            out += theta[self.index[f'coef_{key}']] * (val - np.mean(val))
        return out

    def mean_model(self, theta: np.ndarray) -> np.ndarray:
        """Everything but the GP and the noise, at the data"""
        out = self.systematics(theta)
        for ip in range(len(self.planets)):
            out += self.planet_rv(theta, ip)
        return out

    def noise(self, theta: np.ndarray) -> Tuple[np.ndarray, float]:
        """The noise diagonal and the sequence variance (a float, or one
        per sequence with a sequence jitter per instrument)"""
        diag = self.data.err ** 2
        for it, inst in enumerate(self.data.instruments):
            name = f'log_jit_{inst}'
            if name in self.index:
                diag = diag + np.where(self.inst_index == it,
                                       np.exp(2 * theta[self.index[name]]),
                                       0.0)
        if getattr(self, 'seq_jitter_insts', []):
            # one variance per sequence, from its instrument's jitter
            seq_var = np.zeros(self.data.nseq)
            for inst in self.seq_jitter_insts:
                it = self.data.instruments.index(inst)
                seq_var[self.seq_inst == it] = np.exp(
                    2 * theta[self.index[f'log_sjit_{inst}']])
            return diag, seq_var
        seq_var = (float(np.exp(2 * theta[self.index['log_sjit']]))
                   if 'log_sjit' in self.index else 0.0)
        return diag, seq_var

    def outlier_params(self, theta: np.ndarray, kind: str = 'point'
                       ) -> Tuple[float, float]:
        """
        The outlier fraction and width

        :param kind: str, 'point' (or the only kind), or 'sequence' (the
                     second kind of unit='both')
        """
        sfx = 's' if (kind == 'sequence' and self.unit == 'both') else ''
        frac = float(np.clip(expit(theta[self.index['logit_f' + sfx]]), 1e-12,
                             1 - 1e-12))
        return frac, float(np.exp(theta[self.index['log_W' + sfx]]))

    def gp_scale(self, theta: np.ndarray) -> Optional[np.ndarray]:
        """The relative amplitude of the GP at each point, or None when the
        GP is the same for every instrument"""
        insts = getattr(self, 'gp_scale_insts', [])
        if self.gp is None or not insts or self.gp.get('kernel') == 'multi':
            return None
        scale = np.ones(self.data.n)
        for inst in insts:
            it = self.data.instruments.index(inst)
            scale[self.inst_index == it] = np.exp(
                theta[self.index[f'gp_log_scale_{inst}']])
        return scale

    def dense_gp(self, theta: np.ndarray, diag: np.ndarray,
                 block: np.ndarray, blockval: np.ndarray) -> 'kgp.DenseGP':
        """The GP with the noise, factorised: one kernel (chromatic or
        not), or one per group of instruments"""
        if self.gp.get('kernel') == 'multi':
            kern = kgp.MultiKernel([
                (group['kernel'], [theta[self.index[name]]
                                   for name in group['params']],
                 group['idx']) for group in self.gp_groups])
            return kgp.DenseGP(self.data.time, kern, None, diag, block,
                               blockval)
        return kgp.DenseGP(self.data.time, self.gp['kernel'],
                           self.gp_pars(theta), diag, block, blockval,
                           scale=self.gp_scale(theta))

    def gp_pars(self, theta: np.ndarray) -> np.ndarray:
        """The parameters of the GP kernel"""
        return np.array([theta[self.index[f'gp_{name}']] for name in
                         kgp.KERNEL_PARAMS[self.gp['kernel']]])

    # -------------------------------------------------------------------------
    def log_prior(self, theta: np.ndarray) -> float:
        """The log prior, including the discs of K and e"""
        total = 0.0
        for prior, val in zip(self.priors, theta):
            total += prior.logpdf(val)
            if not np.isfinite(total):
                return -np.inf
        for ip, planet in enumerate(self.planets):
            if f'xk_{ip}' in self.index:
                xk = theta[self.index[f'xk_{ip}']]
                yk = theta[self.index[f'yk_{ip}']]
                if xk ** 2 + yk ** 2 > self.kmax:
                    return -np.inf
            if f'xe_{ip}' in self.index:
                xe = theta[self.index[f'xe_{ip}']]
                ye = theta[self.index[f'ye_{ip}']]
                if xe ** 2 + ye ** 2 > self.emax:
                    return -np.inf
        return total

    def gp_noise_blocks(self, theta: np.ndarray, qunit: np.ndarray,
                        qpoint: Optional[np.ndarray] = None
                        ) -> Tuple[np.ndarray, np.ndarray]:
        """
        The noise of a GP fit, with the outliers inflated

        :param qunit: np.ndarray, the indicators of the units (the visits
                      when unit is 'both')
        :param qpoint: np.ndarray or None, the indicators of the exposures
                       when unit is 'both' (the model's own when None)
        """
        diag, seq_var = self.noise(theta)
        blockval = np.full(self.data.nseq, seq_var)
        if self.likelihood == 'mixture':
            if self.unit == 'both':
                _, wseq = self.outlier_params(theta, 'sequence')
                _, wpt = self.outlier_params(theta, 'point')
                qpoint = self.qpoint if qpoint is None else qpoint
                blockval = blockval + qunit * wseq ** 2
                diag = diag + qpoint * wpt ** 2
            else:
                _, width = self.outlier_params(theta)
                if self.unit == 'sequence':
                    blockval = blockval + qunit * width ** 2
                else:
                    diag = diag + qunit * width ** 2
        return diag, blockval

    def log_likelihood(self, theta: np.ndarray,
                       qunit: Optional[np.ndarray] = None,
                       qpoint: Optional[np.ndarray] = None) -> float:
        """
        The log likelihood

        Without a GP, the outlier indicators are summed out exactly. With a
        GP and a mixture, it is the likelihood given the indicators
        (qunit, or the model's current ones), plus their prior.
        """
        resid = self.data.rv - self.mean_model(theta)
        diag, seq_var = self.noise(theta)
        block, nseq = self.data.seq, self.data.nseq
        if self.gp is None:
            if self.likelihood == 'gaussian':
                return float(np.sum(block_gauss_loglike(
                    resid, diag, block, np.full(nseq, seq_var), nseq)))
            if self.likelihood == 'student':
                nu = self.dof
                scale2 = diag
                return float(np.sum(
                    gammaln((nu + 1) / 2) - gammaln(nu / 2)
                    - 0.5 * np.log(nu * np.pi * scale2)
                    - (nu + 1) / 2 * np.log1p(resid ** 2 / (nu * scale2))))
            frac, width = self.outlier_params(theta)
            if self.unit == 'both':
                fseq, wseq = self.outlier_params(theta, 'sequence')
                per_block, _, _ = mixture_loglike_both(
                    resid, diag, block, nseq, seq_var, frac, width, fseq,
                    wseq)
                return float(np.sum(per_block))
            per_block, _, _ = mixture_loglike(resid, diag, block, nseq,
                                              seq_var, frac, width, self.unit)
            return float(np.sum(per_block))
        if self.likelihood == 'student':
            raise ValueError('A GP with a Student-t likelihood is not '
                             'supported; use the mixture')
        qunit = self.qunit if qunit is None else qunit
        qpoint = self.qpoint if qpoint is None else qpoint
        diag, blockval = self.gp_noise_blocks(theta, qunit, qpoint)
        try:
            dgp = self.dense_gp(theta, diag, block, blockval)
        except np.linalg.LinAlgError:
            return -np.inf
        total = dgp.loglike(resid)
        if self.likelihood == 'mixture':
            kinds = ([(qunit, 'sequence'), (qpoint, 'point')]
                     if self.unit == 'both' else [(qunit, 'point')])
            for arr, kind in kinds:
                frac, _ = self.outlier_params(theta, kind)
                nout = float(np.sum(arr))
                total += nout * np.log(frac) + (len(arr) - nout) * np.log1p(
                    -frac)
        return total

    def log_posterior(self, theta: np.ndarray) -> float:
        """The log posterior (up to a constant)"""
        lprior = self.log_prior(theta)
        if not np.isfinite(lprior):
            return -np.inf
        lnl = self.log_likelihood(theta)
        if not np.isfinite(lnl):
            return -np.inf
        return lprior + lnl

    # -------------------------------------------------------------------------
    def outlier_probability(self, theta: np.ndarray) -> np.ndarray:
        """
        The probability that each point belongs to an outlier unit

        Without a GP, exact given the parameters. With a GP, the exact
        conditional of each unit given the others at their current state.

        :return: np.ndarray, (n) the probability per point
        """
        if self.likelihood != 'mixture':
            return np.zeros(self.data.n)
        resid = self.data.rv - self.mean_model(theta)
        diag, seq_var = self.noise(theta)
        frac, width = self.outlier_params(theta)
        if self.gp is None and self.unit == 'both':
            fseq, wseq = self.outlier_params(theta, 'sequence')
            _, logp, _ = mixture_loglike_both(resid, diag, self.data.seq,
                                              self.data.nseq, seq_var, frac,
                                              width, fseq, wseq)
            return np.exp(logp)
        if self.gp is None:
            _, logp, _ = mixture_loglike(resid, diag, self.data.seq,
                                         self.data.nseq, seq_var, frac,
                                         width, self.unit)
            prob = np.exp(logp)
        else:
            prob = self._gp_outlier_conditionals(theta, resid,
                                                 update=False)[0]
        return prob[self.data.seq] if self.unit == 'sequence' else prob

    def _unit_members(self) -> List[np.ndarray]:
        """The point indices of each unit"""
        if self.unit == 'sequence':
            order = np.argsort(self.data.seq, kind='stable')
            bounds = np.cumsum(np.bincount(self.data.seq,
                                           minlength=self.data.nseq))
            return np.split(order, bounds[:-1])
        return [np.array([it]) for it in range(self.data.n)]

    def _visit_members(self) -> List[np.ndarray]:
        """The point indices of each visit"""
        order = np.argsort(self.data.seq, kind='stable')
        bounds = np.cumsum(np.bincount(self.data.seq,
                                       minlength=self.data.nseq))
        return np.split(order, bounds[:-1])

    def _gp_kinds(self) -> List[Tuple[str, List[np.ndarray], str]]:
        """The latent indicators of a GP fit: (kind, members, attribute)"""
        points = [np.array([it]) for it in range(self.data.n)]
        if self.unit == 'sequence':
            return [('sequence', self._visit_members(), 'qunit')]
        if self.unit == 'point':
            return [('point', points, 'qunit')]
        return [('sequence', self._visit_members(), 'qunit'),
                ('point', points, 'qpoint')]

    def _gp_outlier_conditionals(self, theta: np.ndarray, resid: np.ndarray,
                                 update: bool,
                                 rng: Optional[np.random.Generator] = None,
                                 mode: str = 'draw', by_kind: bool = False):
        """
        Go through the units: exact conditional of each indicator given all
        the others, under the GP (with unit='both', the visits first, then
        the exposures)

        :param update: bool, change the indicators as they are visited
                       ('draw': a Gibbs draw, 'max': the most probable)
        :param by_kind: bool, return the probabilities of each kind apart
        :return: tuple, the conditional probability of each unit (with
                 unit='both', of each exposure being in an outlier, visit
                 or its own; a dict per kind with by_kind), and the number
                 of flips
        """
        diag, blockval = self.gp_noise_blocks(theta, self.qunit)
        dgp = self.dense_gp(theta, diag, self.data.seq, blockval)
        cinv = dgp.inverse()
        alpha = cinv @ resid
        flips = 0
        probs = {}
        for kind, members, attr in self._gp_kinds():
            frac, width = self.outlier_params(theta, kind)
            lf, l1f = np.log(frac), np.log1p(-frac)
            w2 = width ** 2
            state = getattr(self, attr)
            prob = np.zeros(len(members))
            for uu, idx in enumerate(members):
                sub = np.linalg.inv(cinv[np.ix_(idx, idx)])
                mean = resid[idx] - sub @ alpha[idx]
                # the predictive with this unit's own outlier term taken out
                vec = np.ones(len(idx))
                base = sub - state[uu] * w2 * np.outer(vec, vec)
                logl = []
                for qval in (0.0, 1.0):
                    cov = base + qval * w2 * np.outer(vec, vec)
                    sign, logdet = np.linalg.slogdet(cov)
                    diff = resid[idx] - mean
                    logl.append(-0.5 * (diff @ np.linalg.solve(cov, diff)
                                        + logdet))
                lgood, lbad = l1f + logl[0], lf + logl[1]
                prob[uu] = float(np.exp(lbad - np.logaddexp(lgood, lbad)))
                if not update:
                    continue
                if mode == 'draw':
                    new = float(rng.random() < prob[uu])
                else:
                    new = float(prob[uu] > 0.5)
                if new != state[uu]:
                    # C changes by +-W^2 on this unit: carry the inverse along
                    uvec = np.zeros(self.data.n)
                    uvec[idx] = 1.0
                    cinv = kgp.sherman_morrison(cinv, uvec,
                                                w2 * (new - state[uu]))
                    alpha = cinv @ resid
                    state[uu] = new
                    flips += 1
            probs[kind] = prob
        if by_kind:
            return probs, flips
        if self.unit == 'both':
            per_visit = probs['sequence'][self.data.seq]
            return 1 - (1 - per_visit) * (1 - probs['point']), flips
        return probs[self._gp_kinds()[0][0]], flips

    # -------------------------------------------------------------------------
    def _bounds(self, around: Optional[np.ndarray] = None):
        """
        Box bounds for the optimiser

        With around, every free period is also kept inside the periodogram
        peak it starts in (1/T in frequency either side): the line searches
        of Powell and the first simplex of Nelder-Mead would otherwise leap
        across a +-10 % period range to whichever alias the other
        parameters happen to favour at that moment.

        :param around: np.ndarray or None, the start

        :return: list of (low, high)
        """
        bounds = [prior.bounds() for prior in self.priors]
        if around is None:
            return bounds
        for ip, planet in enumerate(self.planets):
            if 'tc' in planet:
                continue
            idx = self.index[f'logP_{ip}']
            low, high = bounds[idx]
            half = float(np.exp(around[idx])) / max(self.data.baseline, 1e-9)
            bounds[idx] = (max(low, around[idx] - half),
                           min(high, around[idx] + half))
        return bounds

    def _optimise(self, start: np.ndarray, maxiter: int = 20000,
                  local: bool = True) -> Tuple[np.ndarray, float]:
        """
        Powell then Nelder-Mead on minus the log posterior

        :param start: np.ndarray, the start
        :param maxiter: int, the iterations of each method
        :param local: bool, keep every period inside the peak it starts in
                      (see _bounds)

        :return: tuple, the maximum and its log posterior
        """
        def _neg(theta):
            val = self.log_posterior(theta)
            return 1e25 if not np.isfinite(val) else -val

        best = np.array(start, dtype=float)
        bounds = self._bounds(around=best if local else None)
        out = minimize(_neg, best, method='Powell', bounds=bounds,
                       options=dict(maxiter=maxiter, xtol=1e-6, ftol=1e-10))
        if out.fun < _neg(best):
            best = out.x
        # Nelder-Mead takes bounds only where every one is finite
        nm_bounds = None
        if all(low is not None and high is not None for low, high in bounds):
            nm_bounds = bounds
        elif local:
            # the free periods alone: a first simplex of 5 % in log P is a
            #   jump of about 12 % in period
            big = 1e10
            nm_bounds = [(low if low is not None else -big,
                          high if high is not None else big)
                         for low, high in bounds]
        out = minimize(_neg, best, method='Nelder-Mead', bounds=nm_bounds,
                       options=dict(maxiter=maxiter, xatol=1e-7, fatol=1e-10,
                                    adaptive=True))
        if out.fun < _neg(best):
            best = out.x
        return best, -_neg(best)

    def fit(self, nstart: int = 4, seed: int = 1, quiet: bool = False
            ) -> 'FitResult':
        """
        The maximum a posteriori, with a Laplace covariance

        A few starts (the phase of each eccentric orbit's periastron turned
        around), then the best is kept. With a GP and a mixture, the
        parameters and the indicators are alternated until the indicators
        stop changing.

        :param nstart: int, the number of starts
        :param seed: int, for the random parts of the starts
        :param quiet: bool, no log lines

        :return: FitResult
        """
        tstart = _time.time()
        rng = np.random.default_rng(seed)
        with blas_threads(1):
            eccentric = [ip for ip in range(len(self.planets))
                         if f'xe_{ip}' in self.index]
            starts = self._period_starts()
            # starts that ignore the most deviant units (many outliers)
            starts += self._trimmed_starts(starts[0])
            if self.gp is not None and self.likelihood == 'mixture':
                best, logpost = self._fit_gp_mixture(starts[:max(nstart, 1)])
            else:
                # round 1: every period candidate, quickly, each inside its
                #   own peak
                results = [self._optimise(st, maxiter=2500) for st in starts]
                order = np.argsort([-res[1] for res in results])
                # round 2: from the two best candidates at distinct periods
                #   (an eccentric orbit can lose round 1 before its
                #   eccentricity is found), the periastron turned around,
                #   and everything polished
                keep = []
                for idx in order:
                    cand = results[idx][0]
                    if all(not self._same_periods(cand, other)
                           for other in keep):
                        keep.append(cand)
                    if len(keep) == 2:
                        break
                trials = []
                for cand in keep:
                    trials.append(cand)
                    for it in range(1, nstart if eccentric else 1):
                        trial = cand.copy()
                        ang = 2 * np.pi * it / nstart
                        for ip in eccentric:
                            trial[self.index[f'xe_{ip}']] = 0.3 * np.cos(ang)
                            trial[self.index[f'ye_{ip}']] = 0.3 * np.sin(ang)
                        trials.append(trial)
                results = [self._optimise(tr) for tr in trials]
                best, logpost = max(results, key=lambda res: res[1])
            cov = self.laplace_covariance(best)
        result = FitResult(self, best, logpost, cov=cov)
        result.runtime = _time.time() - tstart
        if not quiet:
            log(f'Fit ({self.describe()}): log posterior = {logpost:.2f}, '
                f'{result.runtime:.1f} s')
        return result

    def _same_periods(self, theta1: np.ndarray, theta2: np.ndarray) -> bool:
        """Whether two points hold every orbit in the same periodogram peak
        (within 1/T in frequency)"""
        for ip in range(len(self.planets)):
            idx = self.index[f'logP_{ip}']
            dfreq = abs(np.exp(-theta1[idx]) - np.exp(-theta2[idx]))
            if dfreq * self.data.baseline > 1.0:
                return False
        return True

    def _period_starts(self, npeak: int = 3) -> List[np.ndarray]:
        """
        Starting points at the best peaks inside each orbit's period range

        A period range of +-10 % over a long baseline holds several
        aliases (the yearly ones are 1/T_year apart in frequency), and a
        local optimiser keeps whichever it starts in. The outlier-aware
        periodogram of the data picks the candidates; each orbit is moved
        to each of its candidates in turn, the others staying put.

        :param npeak: int, candidates per orbit

        :return: list of np.ndarray, the starting points
        """
        from koloa.periodogram import profile_periodogram
        starts = [self.init.copy()]
        data = self.data
        for ip, planet in enumerate(self.planets):
            if 'tc' in planet:
                continue
            plo, phi = planet['period_range']
            step = 1 / (10 * data.baseline)
            freq = np.arange(1 / phi, 1 / plo + step, step)
            if len(freq) < 20:
                continue
            src = data.binned() if data.nseq < data.n else data
            prof = profile_periodogram(src.time, src.rv, src.err, freq,
                                       outliers=self.likelihood != 'gaussian',
                                       niter=15)
            for idx in find_peaks(freq, prof['dlnl'], npeak,
                                  separation=1 / data.baseline):
                trial = self.init.copy()
                trial[self.index[f'logP_{ip}']] = -np.log(freq[idx])
                self.init = trial
                self._smart_start()
                starts.append(self.init.copy())
        self.init = starts[0]
        return starts

    def _fit_gp_mixture(self, starts: List[np.ndarray]
                        ) -> Tuple[np.ndarray, float]:
        """
        Alternate the parameters and the indicators (hard EM), from two
        seeds of the indicators, keeping the better joint posterior

        With a GP, bad visits can hide: fitted without outliers, the GP and
        the visit jitter grow until nothing looks bad, and the alternation
        never leaves that state. So the indicators are also seeded from the
        leave-one-unit-out residuals of the gaussian fit: each visit judged
        against what the GP predicts from all the OTHER points, flagged
        when it stands out from the other visits by more than 3 robust
        sigma.
        """
        best_all, logp_all, q_all = None, -np.inf, None
        nunit, npts = len(self.qunit), self.data.n
        for start in starts:
            # the gaussian GP first (no indicator on)
            self.qunit, self.qpoint = np.zeros(nunit), np.zeros(npts)
            theta0, _ = self._optimise(start, maxiter=4000)
            # no outliers; the units that stand out by 3 robust sigma; and
            #   a trimmed start, the most deviant 30 % of the visits and
            #   10 % of the exposures (inflated jitters can hide the
            #   others)
            seeds = [(np.zeros(nunit), np.zeros(npts)),
                     self._loo_seed(theta0),
                     self._loo_seed(theta0, top=dict(sequence=0.3,
                                                     point=0.1))]
            unique = []
            for seed in seeds:
                if not any(all(np.array_equal(aa, bb) for aa, bb in
                               zip(seed, other)) for other in unique):
                    unique.append(seed)
            seeds = unique
            for seed_unit, seed_point in seeds:
                self.qunit = seed_unit.copy()
                self.qpoint = seed_point.copy()
                theta = theta0
                for _ in range(12):
                    theta, logpost = self._optimise(theta, maxiter=4000)
                    resid = self.data.rv - self.mean_model(theta)
                    _, flips = self._gp_outlier_conditionals(
                        theta, resid, update=True, mode='max')
                    if flips == 0:
                        break
                logpost = self.log_posterior(theta)
                if logpost > logp_all:
                    best_all, logp_all = theta, logpost
                    q_all = (self.qunit.copy(), self.qpoint.copy())
        self.qunit, self.qpoint = q_all
        return best_all, logp_all

    def _loo_seed(self, theta: np.ndarray, nsig: float = 3.0,
                  top: Optional[Dict[str, float]] = None
                  ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Units whose leave-one-unit-out residual stands out from the others

        :param theta: np.ndarray, the parameters (of a gaussian GP fit)
        :param nsig: float, the threshold [robust sigma of the z of all
                     units]
        :param top: dict or None, instead of a threshold, the fraction of
                    each kind to flag, the most deviant first (a trimmed
                    start: the alternation then releases the units that do
                    not need to be outliers)

        :return: tuple, the seed indicators of the units (qunit) and of the
                 exposures (qpoint, used with unit='both')
        """
        resid = self.data.rv - self.mean_model(theta)
        diag, blockval = self.gp_noise_blocks(theta, np.zeros(len(
            self.qunit)), np.zeros(self.data.n))
        dgp = self.dense_gp(theta, diag, self.data.seq, blockval)
        cinv = dgp.inverse()
        seeds = dict(qunit=np.zeros(len(self.qunit)),
                     qpoint=np.zeros(self.data.n))
        for _, members, attr in self._gp_kinds():
            preds = kgp.loo_blocks(cinv, resid, members)
            zval = np.zeros(len(members))
            for uu, (idx, (mean, cov)) in enumerate(zip(members, preds)):
                # the mean offset of the unit, in units of its uncertainty
                ones = np.ones(len(idx))
                zval[uu] = (ones @ (resid[idx] - mean)) / np.sqrt(
                    ones @ cov @ ones)
            scale = 1.4826 * np.median(np.abs(zval - np.median(zval)))
            zrob = np.abs(zval - np.median(zval)) / max(scale, 1e-9)
            if top is None:
                seeds[attr] = (zrob > nsig).astype(float)
            else:
                kind = 'point' if attr == 'qpoint' or self.unit == 'point' \
                    else 'sequence'
                nflag = int(round(top.get(kind, 0.0) * len(members)))
                flag = np.zeros(len(members))
                flag[np.argsort(-zrob)[:nflag]] = 1.0
                seeds[attr] = flag
        return seeds['qunit'], seeds['qpoint']

    def _initial_outliers(self, theta: np.ndarray) -> np.ndarray:
        """Units more than 5 robust sigma off the mean model"""
        resid = self.data.rv - self.mean_model(theta)
        zval = np.abs(resid - np.median(resid)) / max(robust_std(resid), 1e-9)
        if self.unit == 'sequence':
            zseq = (np.bincount(self.data.seq, weights=zval,
                                minlength=self.data.nseq)
                    / np.bincount(self.data.seq, minlength=self.data.nseq))
            return (zseq > 5).astype(float)
        return (zval > 5).astype(float)

    def laplace_covariance(self, theta: np.ndarray) -> Optional[np.ndarray]:
        """
        The covariance of the parameters from the curvature at the maximum

        :param theta: np.ndarray, the maximum

        :return: np.ndarray or None, the covariance
        """
        ndim = self.ndim
        f0 = self.log_posterior(theta)
        # a step per parameter that changes the log posterior by about 0.1:
        #   a fixed fraction would be far too big for log P over a long
        #   baseline and far too small for an offset
        step = np.maximum(1e-3 * np.abs(theta), 1e-3)
        for ii in range(ndim):
            for _ in range(30):
                trial = theta.copy()
                trial[ii] += step[ii]
                fval = self.log_posterior(trial)
                if not np.isfinite(fval):
                    trial[ii] -= 2 * step[ii]
                    fval = self.log_posterior(trial)
                drop = f0 - fval if np.isfinite(fval) else np.inf
                if drop > 1.0:
                    step[ii] /= 3
                elif drop < 0.01 and step[ii] < 1e3:
                    step[ii] *= 3
                else:
                    break
        hess = np.zeros((ndim, ndim))
        for ii in range(ndim):
            for jj in range(ii, ndim):
                tpp = theta.copy()
                tpm = theta.copy()
                tmp = theta.copy()
                tmm = theta.copy()
                tpp[ii] += step[ii]
                tpp[jj] += step[jj]
                tpm[ii] += step[ii]
                tpm[jj] -= step[jj]
                tmp[ii] -= step[ii]
                tmp[jj] += step[jj]
                tmm[ii] -= step[ii]
                tmm[jj] -= step[jj]
                vals = [self.log_posterior(tt) for tt in (tpp, tpm, tmp, tmm)]
                if not np.all(np.isfinite(vals)):
                    vals = [f0 if not np.isfinite(val) else val
                            for val in vals]
                hess[ii, jj] = (vals[0] - vals[1] - vals[2] + vals[3]) / (
                    4 * step[ii] * step[jj])
                hess[jj, ii] = hess[ii, jj]
        # a uniform prior has no curvature, so a direction the data do not
        #   constrain (a jitter below the error bars, say) is given the
        #   width of its prior rather than an infinite one
        curv = -hess + np.diag([pr.curvature() for pr in self.priors])
        try:
            evals, evecs = np.linalg.eigh(curv)
        except np.linalg.LinAlgError:
            return None
        # a direction without curvature (or a saddle the optimiser stopped
        #   near) gets at most the variance of the widest prior
        prior_curv = [pr.curvature() for pr in self.priors]
        floor = min([val for val in prior_curv if val > 0] or [1e-12])
        evals = np.where(evals > floor, evals, floor)
        return (evecs / evals) @ evecs.T

    def describe(self) -> str:
        """A few words on what the model is"""
        parts = [f'{len(self.planets)} orbit(s)']
        if self.gp is not None and self.gp.get('kernel') == 'multi':
            parts.append('GP ' + ', '.join(
                f'{group["kernel"]} ({group["name"]})'
                for group in self.gp_groups))
        elif self.gp is not None:
            parts.append(f'GP {self.gp["kernel"]}')
        like = self.likelihood
        if like == 'mixture':
            like += f' ({self.unit} outliers)'
        parts.append(like)
        return ', '.join(parts)

    # -------------------------------------------------------------------------
    def sample(self, start: Optional['FitResult'] = None, nsteps: int = 4000,
               nburn: int = 1500, nwalkers: Optional[int] = None,
               seed: int = 1, thin: int = 10, quiet: bool = False,
               moves: str = 'de', converge: bool = False,
               max_steps: Optional[int] = None,
               indicators: str = 'sample') -> 'FitResult':
        """
        The posterior: emcee without a GP (or with a gaussian one),
        Metropolis within Gibbs with a GP and outliers

        :param start: FitResult or None, the maximum to start from (fitted
                      when None)
        :param nsteps: int, the number of steps
        :param nburn: int, the number of steps thrown away (at least five
                      autocorrelation times when converge is set)
        :param nwalkers: int or None, walkers of emcee (4 x ndim, at least 32)
        :param seed: int, the seed
        :param thin: int, keep one step in thin
        :param quiet: bool, no log lines
        :param moves: str, the proposals of emcee: 'de' (differential
                      evolution, 80 %, and its snooker variant, 20 %) or
                      'stretch' (the affine-invariant stretch move)
        :param converge: bool, keep running until the chain holds 50
                         autocorrelation times of every parameter after the
                         burn-in (or max_steps): emcee is extended in place;
                         Metropolis within Gibbs is run again at the length
                         a pilot run asks for
        :param max_steps: int or None, the most steps when converging
                          (20 x nsteps when None)
        :param indicators: str, with a GP and a mixture: 'sample' (the
                           indicators drawn by Gibbs, Metropolis within
                           Gibbs, exact but slow to mix when the GP and the
                           orbit are correlated) or 'map' (held at their
                           most probable state, found by fit, and the other
                           parameters sampled by emcee: conditional on that
                           state, fast, and exact when no indicator is in
                           doubt)

        :return: FitResult, with the chain
        """
        tstart = _time.time()
        if start is None:
            start = self.fit(quiet=quiet)
        rng = np.random.default_rng(seed)
        with blas_threads(1):
            if self.gp is not None and self.likelihood == 'mixture' and \
                    indicators == 'sample':
                chain, lnp, qprob, info = self._mwg(start, nsteps, nburn, rng,
                                                    thin)
                if converge:
                    # a pilot run measures the autocorrelation time; the
                    #   chain is then run again, long enough for 60 of it
                    #   after a burn-in of five
                    tau = np.asarray(info['tau'], dtype=float)
                    taumax = (float(np.nanmax(tau)) if np.any(np.isfinite(tau))
                              else np.nan)
                    if not np.isfinite(taumax):
                        # a pilot too short to measure it: take its own
                        #   length, the cautious guess
                        taumax = float(nsteps)
                    burn2 = max(nburn, int(np.ceil(5 * taumax)))
                    need = burn2 + int(np.ceil(60 * taumax))
                    if need > nsteps:
                        nsteps = min(need, max_steps or 20 * nsteps)
                        nburn = min(burn2, nsteps // 2)
                        chain, lnp, qprob, info = self._mwg(
                            start, nsteps, nburn, rng, thin)
            else:
                chain, lnp, info = self._emcee(
                    start, nsteps, nburn, nwalkers, rng, thin, moves,
                    (max_steps or 20 * nsteps) if converge else nsteps)
                qprob = None
            if qprob is None and self.likelihood == 'mixture':
                pick = chain[rng.choice(len(chain), min(400, len(chain)),
                                        replace=False)]
                qprob = np.mean([self.outlier_probability(th)
                                 for th in pick], axis=0)
        best = chain[int(np.argmax(lnp))]
        if lnp.max() < start.logpost:
            best = start.theta
        result = FitResult(self, best, float(max(lnp.max(), start.logpost)),
                           cov=np.cov(chain.T), chain=chain, lnp=lnp)
        if qprob is not None:
            result.outlier_prob = qprob
        result.runtime = _time.time() - tstart
        result.extra['mcmc'] = self._diagnostics(
            info, info.get('nkept', nsteps - nburn), quiet)
        if not quiet:
            diag = result.extra['mcmc']
            log(f'Posterior ({self.describe()}): {len(chain)} samples, '
                f'{result.runtime:.0f} s')
            log(f'  {diag["sampler"]}: acceptance {diag["acceptance"]:.2f}, '
                f'longest autocorrelation time {diag["tau_max"]:.0f} steps '
                f'({diag["tau_worst"]}), {diag["n_eff_min"]:.0f} effective '
                f'samples', 'value')
        return result

    def _diagnostics(self, info: Dict[str, Any], nkept: int,
                     quiet: bool) -> Dict[str, Any]:
        """
        Acceptance, autocorrelation times and effective sample sizes

        The chain is long enough when every parameter has run for 50
        autocorrelation times after the burn-in (the rule of thumb of
        emcee); otherwise a warning says which parameter mixes slowest.
        """
        tau = np.asarray(info['tau'], dtype=float)
        nchain = info.get('nwalkers', 1)
        worst = int(np.nanargmax(tau)) if np.any(np.isfinite(tau)) else 0
        n_eff = nchain * nkept / np.maximum(tau, 1.0)
        diag = dict(sampler=info['sampler'], nwalkers=nchain, nkept=nkept,
                    nsteps=int(info.get('nsteps', 0)),
                    nburn=int(info.get('nburn', 0)),
                    acceptance=float(info['acceptance']),
                    tau={name: float(val)
                         for name, val in zip(self.names, tau)},
                    tau_max=float(tau[worst]), tau_worst=self.names[worst],
                    n_eff={name: float(val)
                           for name, val in zip(self.names, n_eff)},
                    n_eff_min=float(np.min(n_eff)),
                    converged=bool(nkept >= 50 * tau[worst]))
        if not diag['converged'] and not quiet:
            log(f'The chain ran {nkept} steps after burn-in, fewer than 50 '
                f'autocorrelation times of {diag["tau_worst"]} '
                f'({diag["tau_max"]:.0f} steps): the posterior may be '
                f'under-sampled; increase nsteps', 'warn')
        return diag

    def _emcee(self, start: 'FitResult', nsteps: int, nburn: int,
               nwalkers: Optional[int], rng: np.random.Generator, thin: int,
               moves: str = 'de', max_steps: Optional[int] = None):
        """
        The ensemble sampler of emcee, started around the MAP

        With max_steps above nsteps, the chain is extended until it holds
        50 autocorrelation times of every parameter after a burn-in of at
        least five.
        """
        try:
            import emcee
        except ImportError:
            loud_warning('emcee is not installed: the posterior is sampled '
                         'with koloa\'s own Metropolis, which is slower to '
                         'mix (pip install emcee)')
            chain, lnp, info = self._metropolis(start, nsteps * 8, nburn * 8,
                                                rng, thin)
            return chain, lnp, info
        ndim = self.ndim
        nwalkers = nwalkers or max(4 * ndim, 32)
        pos = self._ball(start, nwalkers, rng)
        if moves == 'de':
            proposal = [(emcee.moves.DEMove(), 0.8),
                        (emcee.moves.DESnookerMove(), 0.2)]
        elif moves == 'stretch':
            proposal = None
        else:
            raise ValueError(f'Unknown moves: {moves}')
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_posterior,
                                        moves=proposal)
        sampler.random_state = np.random.RandomState(
            int(rng.integers(2 ** 31))).get_state()
        state = sampler.run_mcmc(pos, nsteps, progress=False)
        total, burn = nsteps, nburn
        max_steps = max(max_steps or nsteps, nsteps)
        while True:
            tau = integrated_time(sampler.get_chain(discard=burn))
            taumax = float(np.nanmax(tau))
            burn = min(max(nburn, int(np.ceil(5 * taumax))), total // 2)
            need = burn + int(np.ceil(50 * taumax))
            if total >= need or total >= max_steps:
                break
            more = min(max(int(1.2 * need) - total, nsteps // 4),
                       max_steps - total)
            state = sampler.run_mcmc(state, more, progress=False)
            total += more
        tau = integrated_time(sampler.get_chain(discard=burn))
        chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
        lnp = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
        self.autocorr = tau
        info = dict(sampler=f'emcee ({moves} moves)', nwalkers=nwalkers,
                    tau=tau, nkept=total - burn, nsteps=total, nburn=burn,
                    acceptance=float(np.mean(sampler.acceptance_fraction)))
        return chain, lnp, info

    def _ball(self, start: 'FitResult', nwalkers: int,
              rng: np.random.Generator) -> np.ndarray:
        """
        Walkers around the maximum, all inside the prior

        Each parameter is spread by a tenth of its Laplace sigma, clipped to
        between 1e-6 and 1e-2 of the width of its prior, so a curvature the
        finite differences got wrong can neither scatter the walkers out of
        the prior nor stack them on top of each other; the ensemble then
        finds the shape of the posterior during the burn-in.
        """
        ndim = self.ndim
        width = np.empty(ndim)
        for it, prior in enumerate(self.priors):
            low, high = prior.bounds()
            if low is not None and high is not None:
                width[it] = high - low
            elif prior.kind == 'normal':
                width[it] = 10 * prior.b
            else:
                width[it] = 10.0
        sig = np.full(ndim, np.nan)
        if start.cov is not None:
            var = np.diag(start.cov)
            good = np.isfinite(var) & (var > 0)
            sig[good] = 0.1 * np.sqrt(var[good])
        sig = np.where(np.isfinite(sig), sig, 1e-3 * width)
        sig = np.clip(sig, 1e-6 * width, 1e-2 * width)
        pos, tries, scale = [], 0, 1.0
        while len(pos) < nwalkers:
            tries += 1
            trial = start.theta + scale * sig * rng.normal(size=ndim)
            if np.isfinite(self.log_posterior(trial)):
                pos.append(trial)
            elif tries % (10 * nwalkers) == 0:
                scale /= 2
        return np.array(pos)

    def _metropolis(self, start: 'FitResult', nsteps: int, nburn: int,
                    rng: np.random.Generator, thin: int, gibbs: bool = False):
        """Adaptive random-walk Metropolis (optionally with Gibbs outliers)"""
        ndim = self.ndim
        theta = start.theta.copy()
        cov = start.cov if start.cov is not None else np.eye(ndim) * 1e-4
        chol = np.linalg.cholesky(cov + 1e-12 * np.eye(ndim))
        scale = 2.38 / np.sqrt(ndim)
        lnp = self.log_posterior(theta)
        chain, lnps = [], []
        accepted, history, kept_accept = 0, [], 0
        nq = np.zeros(self.data.n if self.unit == 'both'
                      else len(self.qunit))
        for it in range(nsteps):
            prop = theta + scale * chol @ rng.normal(size=ndim)
            lprop = self.log_posterior(prop)
            if np.log(rng.random()) < lprop - lnp:
                theta, lnp = prop, lprop
                accepted += 1
                kept_accept += int(it >= nburn)
            if gibbs:
                resid = self.data.rv - self.mean_model(theta)
                prob, _ = self._gp_outlier_conditionals(theta, resid, True,
                                                        rng, 'draw')
                lnp = self.log_posterior(theta)
                if it >= nburn:
                    nq += prob
            if it < nburn:
                history.append(theta.copy())
                if (it + 1) % 200 == 0:
                    rate = accepted / 200
                    scale *= np.exp(rate - 0.234)
                    accepted = 0
                    if len(history) > 4 * ndim:
                        emp = np.cov(np.array(history[len(history) // 2:]).T)
                        try:
                            chol = np.linalg.cholesky(
                                emp + 1e-10 * np.eye(ndim))
                        except np.linalg.LinAlgError:
                            pass
            elif (it - nburn) % thin == 0:
                chain.append(theta.copy())
                lnps.append(lnp)
        chain, lnps = np.array(chain), np.array(lnps)
        nkept = max(nsteps - nburn, 1)
        # the chain is thinned: its autocorrelation time is in kept samples
        info = dict(sampler='metropolis within gibbs' if gibbs
                    else 'metropolis', nwalkers=1,
                    tau=thin * integrated_time(chain[:, None, :]),
                    acceptance=kept_accept / nkept, nkept=nkept,
                    nsteps=nsteps, nburn=nburn)
        if gibbs:
            return chain, lnps, nq / nkept, info
        return chain, lnps, info

    def _mwg(self, start: 'FitResult', nsteps: int, nburn: int,
             rng: np.random.Generator, thin: int):
        """Metropolis within Gibbs: parameters, then every indicator"""
        chain, lnp, qunit, info = self._metropolis(start, nsteps, nburn,
                                                   rng, thin, gibbs=True)
        qpoint = qunit[self.data.seq] if self.unit == 'sequence' else qunit
        return chain, lnp, qpoint, info


# =============================================================================
# MCMC diagnostics
# =============================================================================
def autocorrelation(series: np.ndarray) -> np.ndarray:
    """
    The normalised autocorrelation function of one series, by FFT

    :param series: np.ndarray, (n)

    :return: np.ndarray, (n), one at lag zero
    """
    series = np.asarray(series, dtype=float)
    nn = len(series)
    size = 1 << int(np.ceil(np.log2(2 * nn)))
    dev = series - np.mean(series)
    spec = np.fft.rfft(dev, n=size)
    acf = np.fft.irfft(spec * np.conj(spec), n=size)[:nn]
    if acf[0] <= 0:
        return np.ones(nn)
    return acf / acf[0]


def integrated_time(chain: np.ndarray, window: float = 5.0) -> np.ndarray:
    """
    The integrated autocorrelation time of every parameter

    The autocorrelation function is averaged over the walkers of an
    ensemble, and summed up to the smallest lag M with M >= window * tau(M)
    (the automatic window of Sokal 1997, as in emcee).

    :param chain: np.ndarray, (nstep x nwalker x ndim), or (nstep x ndim)
                  for a single chain
    :param window: float, the window constant

    :return: np.ndarray, (ndim), in steps
    """
    chain = np.asarray(chain, dtype=float)
    if chain.ndim == 2:
        chain = chain[:, None, :]
    nstep, nwalker, ndim = chain.shape
    tau = np.full(ndim, np.nan)
    if nstep < 4:
        return tau
    for ip in range(ndim):
        acf = np.mean([autocorrelation(chain[:, iw, ip])
                       for iw in range(nwalker)], axis=0)
        taus = 2.0 * np.cumsum(acf) - 1.0
        lags = np.arange(len(taus))
        ok = lags >= window * taus
        tau[ip] = taus[np.argmax(ok)] if np.any(ok) else taus[-1]
    return tau


# =============================================================================
# The result
# =============================================================================
@dataclass
class FitResult:
    """
    A fitted RVModel: the best parameters, their covariance, the chain
    """
    model: RVModel
    theta: np.ndarray
    logpost: float
    cov: Optional[np.ndarray] = None
    chain: Optional[np.ndarray] = None
    lnp: Optional[np.ndarray] = None
    outlier_prob: Optional[np.ndarray] = None
    runtime: float = 0.0
    extra: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.outlier_prob is None:
            self.outlier_prob = self.model.outlier_probability(self.theta)

    # -------------------------------------------------------------------------
    @property
    def reliability(self) -> np.ndarray:
        """One minus the outlier probability of each point"""
        return 1.0 - self.outlier_prob

    def samples(self, nsample: int = 300, seed: int = 2) -> np.ndarray:
        """
        Parameter draws: from the chain, or from the Laplace covariance

        :param nsample: int, how many
        :param seed: int, the seed

        :return: np.ndarray, (nsample x ndim)
        """
        rng = np.random.default_rng(seed)
        if self.chain is not None and len(self.chain):
            pick = rng.choice(len(self.chain), min(nsample, len(self.chain)),
                              replace=False)
            return self.chain[pick]
        if self.cov is None:
            return self.theta[None, :]
        draws = rng.multivariate_normal(self.theta, self.cov, size=nsample)
        keep = [dr for dr in draws
                if np.isfinite(self.model.log_prior(dr))]
        # a parameter at a bound of its prior (a jitter at zero, say) loses
        #   half of the draws, and one whose Laplace spread is wider than its
        #   prior nearly all: when too few are left for percentiles, the
        #   draws are brought inside the bounds of the uniform priors instead
        if len(keep) < max(30, 0.02 * nsample):
            low, high = np.array([[-np.inf if lo is None else lo,
                                   np.inf if hi is None else hi]
                                  for lo, hi in (pr.bounds() for pr in
                                                 self.model.priors)]).T
            keep = [dr for dr in np.clip(draws, low, high)
                    if np.isfinite(self.model.log_prior(dr))]
        return np.array(keep) if keep else self.theta[None, :]

    @property
    def diagnostics(self) -> Optional[Dict[str, Any]]:
        """The MCMC diagnostics: sampler, acceptance, autocorrelation times
        (tau, per parameter, in steps), effective sample sizes (n_eff), and
        whether the chain ran for 50 autocorrelation times (converged);
        None for a maximum without a chain"""
        return self.extra.get('mcmc')

    def _natural_name(self, name: str):
        """A noise or linear parameter in natural units: (name, function),
        or None for the parameters of an orbit"""
        unit = self.model.unit
        if name == 'logit_f':
            return ('f_visit' if unit == 'sequence' else 'f_point'), expit
        if name == 'log_W':
            return ('W_visit' if unit == 'sequence' else 'W_point'), np.exp
        if name == 'logit_fs':
            return 'f_visit', expit
        if name == 'log_Ws':
            return 'W_visit', np.exp
        if name == 'log_sjit':
            return 'sjit', np.exp
        if name.startswith('log_sjit_'):
            return 'sjit_' + name[len('log_sjit_'):], np.exp
        if name.startswith('log_jit_'):
            return 'jit_' + name[len('log_jit_'):], np.exp
        if name.startswith('gp_log_'):
            return 'gp_' + name[len('gp_log_'):], np.exp
        if name.startswith('gp_') and '_log_' in name:
            # a GP of a group of instruments: gp_<group>_log_<parameter>
            head, tail = name.split('_log_', 1)
            return f'{head}_{tail}', np.exp
        if name.startswith(('offset_', 'trend_', 'coef_')) or \
                name == 'secacc':
            return name, None
        return None

    def posterior(self, mstar: Optional[float] = None,
                  nsample: Optional[int] = None, seed: int = 2
                  ) -> Dict[str, np.ndarray]:
        """
        The posterior in physical units, one array of draws per quantity

        Orbit ip gives P_ip [days], K_ip [m/s], e_ip, omega_ip [deg, 0 to
        360, the argument of periastron of the star], tp_ip and tc_ip (the
        times of periastron and conjunction, in the time of the data) and,
        when the stellar mass is given, msini_ip [Earth masses]. The noise
        is in natural units: jit_<inst> and sjit (the white and visit
        jitters) [m/s], f_point, W_point, f_visit and W_visit (the outlier
        fractions and widths [m/s]) and gp_<name>; the offsets, trends and
        decorrelation coefficients are as fitted. The arrays are aligned:
        entry j of every key comes from the same draw.

        :param mstar: float or None, the stellar mass [solar masses]
        :param nsample: int or None, at most this many draws (all of the
                        chain when None)
        :param seed: int, the seed of the subsampling

        :return: dict, name: np.ndarray of draws
        """
        if self.chain is not None and len(self.chain):
            draws = self.chain
            if nsample is not None and nsample < len(draws):
                pick = np.random.default_rng(seed).choice(
                    len(draws), nsample, replace=False)
                draws = draws[pick]
        else:
            log('No chain: the posterior is drawn from the Laplace '
                'approximation at the maximum (run RVModel.sample for the '
                'real one)', 'warn')
            draws = self.samples(nsample or 2000, seed=seed)
        out: Dict[str, np.ndarray] = {}
        for ip in range(len(self.model.planets)):
            vals = np.array([self.model.orbit(theta, ip) for theta in draws])
            period, tperi, ecc, omega, amp = vals.T
            out[f'P_{ip}'] = period
            out[f'K_{ip}'] = amp
            out[f'e_{ip}'] = ecc
            out[f'omega_{ip}'] = np.degrees(omega) % 360
            out[f'tp_{ip}'] = tperi
            out[f'tc_{ip}'] = np.array([kepler.tp_to_tc(*row) for row in
                                        zip(tperi, period, ecc, omega)])
            if mstar:
                out[f'msini_{ip}'] = np.array(
                    [kepler.minimum_mass(kk, pp, ee, mstar)
                     for kk, pp, ee in zip(amp, period, ecc)])
        for it, name in enumerate(self.model.names):
            natural = self._natural_name(name)
            if natural is None:
                continue
            key, func = natural
            out[key] = draws[:, it] if func is None else func(draws[:, it])
        return out

    def param(self, name: str) -> Tuple[float, float, float]:
        """
        One parameter: value, lower and upper one sigma

        :param name: str, the parameter

        :return: tuple, (value, minus, plus)
        """
        idx = self.model.index[name]
        if self.chain is not None:
            low, mid, high = np.percentile(self.chain[:, idx], [16, 50, 84])
            return float(mid), float(mid - low), float(high - mid)
        sig = np.sqrt(self.cov[idx, idx]) if self.cov is not None else np.nan
        return float(self.theta[idx]), float(sig), float(sig)

    def orbits(self, mstar: Optional[float] = None
               ) -> List[Dict[str, Tuple[float, float, float]]]:
        """
        Every orbit in physical units, with uncertainties

        :param mstar: float or None, the stellar mass, for m sin i
                      [solar masses]

        :return: list of dict, P, K, e, omega (deg), tp, tc and msini
                 [Earth masses], each as (value, minus, plus): the median
                 and the 16th and 84th percentiles of the chain when there
                 is one, else the maximum and the Laplace spread
        """
        draws = self.samples(2000)
        has_chain = self.chain is not None and len(self.chain) > 0
        out = []
        for ip in range(len(self.model.planets)):
            vals = {key: [] for key in ('P', 'K', 'e', 'omega', 'tp', 'tc',
                                        'msini')}
            for theta in list(draws) + [self.theta]:
                period, tperi, ecc, omega, amp = self.model.orbit(theta, ip)
                vals['P'].append(period)
                vals['K'].append(amp)
                vals['e'].append(ecc)
                vals['omega'].append(np.degrees(omega) % 360)
                vals['tp'].append(tperi)
                vals['tc'].append(kepler.tp_to_tc(tperi, period, ecc, omega))
                vals['msini'].append(kepler.minimum_mass(amp, period, ecc,
                                                         mstar)
                                     if mstar else np.nan)
            summary = {}
            for key, arr in vals.items():
                arr = np.array(arr)
                if key == 'omega':
                    # an angle: unwrapped around its circular mean
                    rad = np.radians(arr)
                    centre = np.degrees(np.arctan2(np.mean(np.sin(rad)),
                                                   np.mean(np.cos(rad))))
                    arr = centre + (arr - centre + 180) % 360 - 180
                best = arr[-1]
                if len(arr) > 3:
                    low, mid, high = np.percentile(arr[:-1], [16, 50, 84])
                    if has_chain:
                        best = mid
                    summary[key] = (float(best), float(best - low),
                                    float(high - best))
                else:
                    summary[key] = (float(best), np.nan, np.nan)
                if key == 'omega':
                    summary[key] = (summary[key][0] % 360,) + \
                        summary[key][1:]
            out.append(summary)
        return out

    def bic(self) -> float:
        """The Bayesian information criterion at the maximum"""
        lnl = self.model.log_likelihood(self.theta)
        return float(self.model.ndim * np.log(self.model.data.n) - 2 * lnl)

    def gp_prediction(self, time: Optional[np.ndarray] = None,
                      theta: Optional[np.ndarray] = None,
                      inst: Optional[str] = None
                      ) -> Tuple[np.ndarray, np.ndarray]:
        """
        The GP conditioned on the data minus the mean model

        The outlier units are inflated in the noise, so they do not pull
        the prediction towards themselves.

        :param time: np.ndarray or None, where [days] (the data if None)
        :param theta: np.ndarray or None, the parameters (the best if None)
        :param inst: str or None, with a chromatic GP, the instrument whose
                     amplitude the prediction at new times takes (the
                     reference when None)

        :return: tuple, the mean and one sigma of the GP
        """
        model = self.model
        at_data = time is None
        if model.gp is None:
            time = model.data.time if time is None else time
            return np.zeros(len(time)), np.zeros(len(time))
        theta = self.theta if theta is None else theta
        time = model.data.time if time is None else np.asarray(time)
        resid = model.data.rv - model.mean_model(theta)
        qpoint = np.zeros(model.data.n)
        if model.likelihood == 'mixture':
            prob = self.outlier_prob
            if model.unit == 'both':
                # a flagged exposure inflated on its own, whatever its kind
                qunit = np.zeros(model.data.nseq)
                qpoint = (prob > 0.5).astype(float)
            elif model.unit == 'sequence':
                unitq = np.array([prob[model.data.seq == ss].max()
                                  for ss in range(model.data.nseq)])
                qunit = (unitq > 0.5).astype(float)
            else:
                qunit = (prob > 0.5).astype(float)
        else:
            qunit = np.zeros(len(model.qunit))
        diag, blockval = model.gp_noise_blocks(theta, qunit, qpoint)
        dgp = model.dense_gp(theta, diag, model.data.seq, blockval)
        if at_data:
            # each point with its own amplitude, or its own group's GP
            return dgp.predict_data(resid)
        if model.gp.get('kernel') == 'multi':
            # the GP of the group of the instrument asked
            groups = [ig for ig, group in enumerate(model.gp_groups)
                      if inst in group['instruments']]
            if inst is None or not groups:
                raise ValueError('A GP per group of instruments: give the '
                                 'instrument (inst) of a group')
            return dgp.predict(resid, time, group=groups[0])
        # the amplitude of the instrument asked (the reference by default)
        if inst is not None and inst in getattr(model, 'gp_scale_insts',
                                                []):
            scale_new = float(np.exp(theta[model.index[
                f'gp_log_scale_{inst}']]))
        else:
            scale_new = None
        return dgp.predict(resid, time, scale_new=scale_new)

    def residuals(self, keep_planet: Optional[int] = None) -> np.ndarray:
        """
        The data minus everything, or minus everything but one orbit

        :param keep_planet: int or None, the orbit to leave in

        :return: np.ndarray, (n)
        """
        model = self.model
        out = model.data.rv - model.systematics(self.theta)
        for ip in range(len(model.planets)):
            if ip != keep_planet:
                out = out - model.planet_rv(self.theta, ip)
        if model.gp is not None:
            out = out - self.gp_prediction()[0]
        return out

    def summary(self, mstar: Optional[float] = None) -> str:
        """
        A few lines on the fit

        :param mstar: float or None, the stellar mass [solar masses]

        :return: str, the summary
        """
        lines = [f'{self.model.describe()}: log posterior '
                 f'{self.logpost:.2f}, BIC {self.bic():.1f}']
        for ip, orb in enumerate(self.orbits(mstar)):
            text = (f'  orbit {ip}: P = {orb["P"][0]:.4f} '
                    f'+{orb["P"][2]:.4f}/-{orb["P"][1]:.4f} d, '
                    f'K = {orb["K"][0]:.2f} +{orb["K"][2]:.2f}/'
                    f'-{orb["K"][1]:.2f} m/s, e = {orb["e"][0]:.3f}')
            if mstar:
                text += f', m sin i = {orb["msini"][0]:.2f} Me'
            lines.append(text)
        if self.model.trend >= 1:
            from koloa.secular import acceleration
            acc = acceleration(self)
            val, low, high = acc['accel']
            lines.append(f'  acceleration of the star: {val:.3f} (-{low:.3f} '
                         f'+{high:.3f}) m/s/yr at tref')
            if 'jerk' in acc:
                val, low, high = acc['jerk']
                lines.append(f'  its change: {val:.3f} (-{low:.3f} '
                             f'+{high:.3f}) m/s/yr^2')
        if 'secacc' in self.model.index:
            val, low, high = self.param('secacc')
            prior = self.model.priors[self.model.index['secacc']]
            told = (f'prior {prior.a:.4f} +- {prior.b:.2g}'
                    if prior.kind == 'normal' else 'free')
            lines.append(f'  perspective acceleration: {val:.4f} (-{low:.4f} '
                         f'+{high:.4f}) m/s/yr ({told})')
        for name in self.model.names:
            if name.startswith(('log_jit', 'log_sjit', 'gp_', 'logit_f',
                                'log_W')):
                val, low, high = self.param(name)
                shown = expit(val) if name.startswith('logit_f') else (
                    np.exp(val) if name.startswith('log') or
                    name.startswith('gp_log') or
                    (name.startswith('gp_') and '_log_' in name) else val)
                lines.append(f'  {name}: {val:.3f} (-{low:.3f} +{high:.3f})'
                             f'  -> {shown:.3f}')
        if self.model.likelihood == 'mixture':
            lines.append(f'  expected outlier points: '
                         f'{np.sum(self.outlier_prob):.1f}')
        diag = self.diagnostics
        if diag is not None:
            lines.append(f'  {diag["sampler"]}: acceptance '
                         f'{diag["acceptance"]:.2f}, longest autocorrelation '
                         f'time {diag["tau_max"]:.0f} steps '
                         f'({diag["tau_worst"]}), {diag["n_eff_min"]:.0f} '
                         f'effective samples'
                         + ('' if diag['converged'] else
                            ' [fewer than 50 autocorrelation times]'))
        return '\n'.join(lines)


# =============================================================================
# Convenience functions
# =============================================================================
def fit_planets(data: RVData, periods: Sequence[float],
                eccentric: bool = False, **kwargs) -> FitResult:
    """
    Orbits at the given periods, fitted with the model's defaults

    :param data: RVData, the series
    :param periods: list of float, the starting periods [days]
    :param eccentric: bool, eccentric orbits
    :param kwargs: passed to RVModel

    :return: FitResult, the maximum a posteriori
    """
    planets = [dict(period=per, eccentric=eccentric) for per in periods]
    return RVModel(data, planets, **kwargs).fit()


def mcmc_orbits(data: RVData, planets: Sequence[Any],
                likelihood: str = 'mixture', eccentric: bool = True,
                nsteps: int = 10000, nburn: int = 2000,
                nwalkers: Optional[int] = None, thin: int = 10,
                max_steps: Optional[int] = None, moves: str = 'de',
                nstart: int = 4, seed: int = 1, quiet: bool = False,
                **kwargs) -> FitResult:
    """
    Keplerian orbits by MCMC: the maximum a posteriori, then the posterior

    The one call from a series to orbital elements with their
    uncertainties. On the result, posterior(mstar) gives every draw in
    physical units (P, K, e, omega, tp, tc, m sin i and the noise),
    orbits(mstar) their summaries, diagnostics the autocorrelation times
    and effective sample sizes, and koloa.plotting.corner draws one or
    several of them on the same axes.

    Without a GP the sampler is emcee on the likelihood with the outlier
    indicators summed out, with differential-evolution moves, and it runs
    until the chain holds 50 autocorrelation times of every parameter
    after a burn-in of at least five (or max_steps). With a GP and
    outliers it is Metropolis within Gibbs, with exact draws of the
    indicators, run again at the length its pilot run's autocorrelation
    times ask for.

    :param data: RVData, the series
    :param planets: list, one per orbit: a period [days] to start from, or
                    a dict as in RVModel (period, period_range, eccentric,
                    or tc, tc_err and period_err for a transiting planet)
    :param likelihood: str, 'mixture' (koloa: outlier-aware), 'gaussian'
                       (the same noise model without outliers) or 'student'
    :param eccentric: bool, eccentric orbits (for the planets that do not
                      say otherwise)
    :param nsteps: int, the steps of every walker (at least)
    :param nburn: int, the steps thrown away (at least)
    :param nwalkers: int or None, the walkers of emcee (4 x ndim, at least
                     32)
    :param thin: int, keep one step in thin
    :param max_steps: int or None, the most steps (20 x nsteps when None)
    :param moves: str, 'de' or 'stretch' (see RVModel.sample)
    :param nstart: int, the starts of the maximum a posteriori
    :param seed: int, the seed
    :param quiet: bool, no log lines
    :param kwargs: passed to RVModel (unit, gp, seq_jitter, trend,
                   frac_prior, priors, ...)

    :return: FitResult, with the chain
    """
    orbits = []
    for planet in planets:
        entry = (dict(period=float(planet)) if np.isscalar(planet)
                 else dict(planet))
        entry.setdefault('eccentric', eccentric)
        orbits.append(entry)
    model = RVModel(data, orbits, likelihood=likelihood, **kwargs)
    start = model.fit(nstart=nstart, seed=seed, quiet=quiet)
    return model.sample(start=start, nsteps=nsteps, nburn=nburn,
                        nwalkers=nwalkers, seed=seed, thin=thin, quiet=quiet,
                        moves=moves, converge=True, max_steps=max_steps)


def sequential_fit(data: RVData, period: float, kernel: str = 'sho',
                   gp_prior: Optional[Dict[str, Any]] = None,
                   **kwargs) -> Dict[str, FitResult]:
    """
    A GP fitted alone, then an orbit fitted to what the GP left: DON'T

    This is kept so that the damage can be shown. A GP flexible enough to
    describe activity is flexible enough to describe a planet; fitted
    first, it takes whatever it can, and the orbit fitted afterwards is
    measured on what is left. On a real NIRPS series this takes K from
    about 6 m/s to a fraction of 1 m/s. Fit the two together (RVModel with
    planets and a GP), with a prior on the GP period from an activity
    indicator.

    :param data: RVData, the series
    :param period: float, the period of the orbit [days]
    :param kernel: str, the GP kernel
    :param gp_prior: dict or None, a prior on the GP parameters
    :param kwargs: passed to RVModel

    :return: dict, gp (the GP-only fit) and orbit (the orbit on the residual)
    """
    loud_warning('SEQUENTIAL FIT: the GP is fitted first and the orbit to '
                 'its residual. This is NOT a joint fit: the GP absorbs '
                 'whatever part of the signal it can, and the amplitude of '
                 'the orbit is biased low (it can collapse). Use a '
                 'joint RVModel(planets=..., gp=...) instead.')
    gp_fit = RVModel(data, [], gp=dict(kernel=kernel, prior=gp_prior),
                     **kwargs).fit(quiet=True)
    resid = data.rv - gp_fit.gp_prediction()[0]
    orbit_fit = RVModel(data.with_values(resid), [dict(period=period)],
                        **kwargs).fit(quiet=True)
    return dict(gp=gp_fit, orbit=orbit_fit)


def period_prior_from_indicator(data: RVData, name: str, kernel: str = 'sho',
                                likelihood: str = 'mixture',
                                period_guess: Optional[float] = None,
                                nsteps: int = 3000, nburn: int = 1000,
                                quiet: bool = False,
                                period_range: Optional[Sequence[float]]
                                = None, trend: int = 1
                                ) -> Tuple[Dict[str, Tuple[float, float]],
                                           FitResult]:
    """
    The period of the activity, measured on an indicator, as a prior

    A GP is fitted to an indicator (DTEMP, dW, FWHM...) with the same
    outlier model as the velocities, its posterior sampled, and the period
    handed back as a gaussian prior in log(P), ready for
    RVModel(gp=dict(kernel=..., prior=...)).

    :param data: RVData, the series (with the indicator)
    :param name: str, the indicator
    :param kernel: str, the kernel
    :param likelihood: str, gaussian or mixture
    :param period_guess: float or None, where the period starts [days]
                         (the highest GLS peak between 2 d and T/3 if None)
    :param nsteps: int, the MCMC steps
    :param nburn: int, the burn-in
    :param quiet: bool, no log lines
    :param period_range: tuple or None, the range of the period [days]
                         (1 day to the baseline when None): an indicator
                         that drifts over the years (a focus change moves
                         the FWHM) can otherwise take the GP to the drift
    :param trend: int, the degree of the polynomial fitted with the GP

    :return: tuple, the prior dict and the fit of the indicator
    """
    series = data.indicator(name)
    low, high = period_range or (2.0, series.baseline / 3)
    if period_guess is None:
        freq = frequency_grid(series.time, low, high, 5)
        power = gls(series.time, series.rv, series.err, freq)
        period_guess = float(1 / freq[find_peaks(freq, power, 1)[0]])
    gpdef = dict(kernel=kernel, init=dict(log_period=np.log(period_guess)))
    if period_range is not None:
        gpdef['prior'] = {'log_period': Prior('uniform', np.log(low),
                                              np.log(high))}
    model = RVModel(series, [], gp=gpdef, likelihood=likelihood,
                    trend=trend)
    post = model.sample(nsteps=nsteps, nburn=nburn, quiet=quiet)
    values = post.chain[:, model.index['gp_log_period']]
    prior = {'log_period': (float(np.mean(values)), float(np.std(values)))}
    if not quiet:
        mid = np.exp(np.mean(values))
        log(f'Activity period from {name}: {mid:.2f} d '
            f'(log sigma {np.std(values):.3f})', 'value')
    return prior, post


# =============================================================================
# End of code
# =============================================================================

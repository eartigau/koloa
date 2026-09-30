#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A GP of the stellar activity inside the FIP, as a finite basis.

The FIP integrates out every linear parameter analytically (koloa.linear),
with the noise kept block-diagonal (a jitter per exposure and per visit, and
the outliers). A GP enters that machinery without changing it when it is
written as a set of basis functions whose weights have a gaussian prior:

- 'local': gaussian bumps of width L every L/3; with the prior variance
  sigma^2 dt / (sqrt(pi/2) L) on each weight, their sum is a GP with the
  squared-exponential covariance sigma^2 exp(-tau^2 / 2 L^2) (to 1e-3 of
  sigma^2). A slow, non-oscillating activity: locally coherent, nothing
  periodic in it; a planet keeps one phase over the whole series and stays
  in the signal slots;
- 'rotation': the same bumps times cos and sin of 2 pi t / P_rot, and of
  4 pi t / P_rot (the first harmonic), each pair with its own amplitude:
  sigma_1^2 exp(-tau^2 / 2 L^2) cos(2 pi tau / P_rot) + sigma_2^2 exp(...)
  cos(4 pi tau / P_rot), a quasi-periodic GP of rotating, evolving spots.

The hyperparameters (sigma, L, P_rot) are sampled in the chain with the
jitters, the linear parameters integrated out, so they are fitted jointly
with the signals: fitted without them, the GP bends toward a planet and
competes with it. Their priors: log-uniform within bounds, or for P_rot a
gaussian on its logarithm (a rotation period from photometry or the
indicators, with its uncertainty).

The cost is a few dozen to a few hundred more columns in the linear model.

Created on 2026-09-30

@author: artigau
"""
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np

# =============================================================================
# Define variables
# =============================================================================
#: the spacing of the bumps, as a fraction of their width
STEP_FRAC = 1.0 / 3
#: the kinds of component
KINDS = ('local', 'rotation')


# =============================================================================
# Define functions
# =============================================================================
def _bounds(value: Any, default: Sequence[float]) -> Dict[str, Any]:
    """a prior: (low, high) log-uniform, a number (fixed), or dict(mu=,
    sd=) a gaussian on the log with bounds low, high"""
    if value is None:
        return dict(kind='loguniform', low=float(default[0]),
                    high=float(default[1]))
    if isinstance(value, dict):
        mu, sd = float(value['mu']), float(value['sd'])
        low = float(value.get('low', np.exp(mu - 5 * sd)))
        high = float(value.get('high', np.exp(mu + 5 * sd)))
        return dict(kind='lognormal', mu=mu, sd=sd, low=low, high=high)
    if np.ndim(value) == 0:
        return dict(kind='fixed', low=float(value), high=float(value))
    return dict(kind='loguniform', low=float(value[0]), high=float(value[1]))


def setup(spec: Union[str, Dict[str, Any], Sequence[Any], None],
          time: np.ndarray, scale: float) -> List[Dict[str, Any]]:
    """
    The components of a GP, with the priors of their hyperparameters

    :param spec: 'local', 'rotation', a dict (kind, and the priors: sigma,
                 length, period, sigma2 for the harmonic of a rotation), or
                 a list of them. A prior is (low, high) (log-uniform), a
                 number (fixed) or dict(mu=ln P, sd=) (gaussian on the log).
                 A rotation needs its period: a dict with period
    :param time: np.ndarray, the times of the series [days]
    :param scale: float, the scale of the velocities [m/s] (the default
                  amplitudes run from 0.05 to 10 times it)

    :return: list of dict, kind and priors (name: prior)
    """
    if spec is None:
        return []
    specs = [spec] if isinstance(spec, (str, dict)) else list(spec)
    base = float(np.ptp(time))
    out = []
    for item in specs:
        item = dict(kind=item) if isinstance(item, str) else dict(item)
        kind = item.get('kind', 'local')
        if kind not in KINDS:
            raise ValueError(f'a GP component is {KINDS}, not {kind}')
        amp = (0.05 * scale, 10 * scale)
        comp = dict(kind=kind, priors={})
        if kind == 'local':
            comp['priors']['length'] = _bounds(
                item.get('length'), (5.0, max(15.0, base / 3)))
            comp['priors']['sigma'] = _bounds(item.get('sigma'), amp)
        else:
            if item.get('period') is None:
                raise ValueError('a rotation GP needs its period: a range '
                                 '(low, high), or dict(mu=ln P, sd=)')
            per = _bounds(item['period'], (1.0, 1.0))
            comp['priors']['period'] = per
            comp['priors']['length'] = _bounds(
                item.get('length'), (0.5 * per['low'],
                                     max(5 * per['high'], 1.0)))
            comp['priors']['sigma'] = _bounds(item.get('sigma'), amp)
            comp['priors']['sigma2'] = _bounds(item.get('sigma2'), amp)
        out.append(comp)
    return out


def initial(components: List[Dict[str, Any]]) -> List[Dict[str, float]]:
    """the starting values: the geometric middle of each prior (its mode
    for a gaussian)"""
    out = []
    for comp in components:
        vals = {}
        for name, prior in comp['priors'].items():
            if prior['kind'] == 'lognormal':
                vals[name] = float(np.exp(prior['mu']))
            else:
                vals[name] = float(np.sqrt(prior['low'] * prior['high']))
        out.append(vals)
    return out


def log_prior(prior: Dict[str, Any], value: float) -> float:
    """the log prior of a hyperparameter, in its logarithm (-inf outside)"""
    if not prior['low'] * (1 - 1e-12) <= value <= prior['high'] * (1 + 1e-12):
        return -np.inf
    if prior['kind'] == 'lognormal':
        return float(-0.5 * ((np.log(value) - prior['mu']) / prior['sd'])
                     ** 2)
    return 0.0


def bumps(time: np.ndarray, length: float,
          step_frac: float = STEP_FRAC) -> np.ndarray:
    """
    Gaussian bumps exp(-(t - c)^2 / L^2) every step_frac L, those that
    touch the data

    :return: np.ndarray, (n x m)
    """
    step = step_frac * length
    cen = np.arange(time.min() - 2 * length, time.max() + 2 * length + step,
                    step)
    phi = np.exp(-((time[:, None] - cen[None, :]) / length) ** 2)
    return phi[:, phi.max(axis=0) > 1e-3]


def columns(components: List[Dict[str, Any]],
            values: List[Dict[str, float]], time: np.ndarray):
    """
    The basis of a GP and the prior variance of each weight

    :param components: list of dict, from setup()
    :param values: list of dict, the hyperparameters of each component
    :param time: np.ndarray, the times [days]

    :return: tuple, (n x m) columns, (m) prior variances
    """
    cols, var = [], []
    for comp, val in zip(components, values):
        phi = bumps(time, val['length'])
        unit = STEP_FRAC / np.sqrt(np.pi / 2)
        if comp['kind'] == 'local':
            cols.append(phi)
            var.append(np.full(phi.shape[1], val['sigma'] ** 2 * unit))
            continue
        for harm, name in ((1, 'sigma'), (2, 'sigma2')):
            arg = 2 * np.pi * harm * time / val['period']
            cols += [phi * np.cos(arg)[:, None], phi * np.sin(arg)[:, None]]
            var += [np.full(2 * phi.shape[1], val[name] ** 2 * unit)]
    if not cols:
        return np.zeros((len(time), 0)), np.zeros(0)
    return np.hstack(cols), np.concatenate(var)


def kernel(components: List[Dict[str, Any]],
           values: List[Dict[str, float]], tau: np.ndarray) -> np.ndarray:
    """
    The covariance the basis stands for, at time lags tau (for the tests
    and the figures)

    :return: np.ndarray
    """
    out = np.zeros_like(np.asarray(tau, dtype=float))
    for comp, val in zip(components, values):
        env = np.exp(-0.5 * (tau / val['length']) ** 2)
        if comp['kind'] == 'local':
            out = out + val['sigma'] ** 2 * env
        else:
            out = out + env * (val['sigma'] ** 2 * np.cos(
                2 * np.pi * tau / val['period']) + val['sigma2'] ** 2
                * np.cos(4 * np.pi * tau / val['period']))
    return out


def names(components: List[Dict[str, Any]]) -> List[str]:
    """the names of the sampled hyperparameters, 'gp<i>_<name>'"""
    return [f'gp{ic}_{name}' for ic, comp in enumerate(components)
            for name, prior in comp['priors'].items()
            if prior['kind'] != 'fixed']


# =============================================================================
# End of code
# =============================================================================

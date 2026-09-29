#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Soft down-weighting of outliers: the classical rule, and its accounting.

The rule is the one pca2d-preclean applies to its spectra (clip_weights)
and LBL's report applies to its velocities: the weight of a point is
untouched inside `clip` sigma and falls as (clip / |z|)^2 outside it. That
exponent is what inflating the variance of an outlying point until it IS a
clip-sigma point gives, so a 6 sigma point keeps a quarter of its weight.
Nothing is hard-rejected: a flare is real data, and a hard cut on the
residual sculpts the variability being measured.

Two habits make it safe to use:
- the scale of z is robust (1.4826 x MAD of the normalised residual), and
  the factor is recomputed from the ORIGINAL error bars at every pass, never
  compounded, so it is iterative re-weighting and not a ratchet;
- the clip is handed over as inflated error bars, error / sqrt(factor), so
  every fit downstream takes it without knowing.

What this module adds is the accounting that the rule does not do on its
own: how much weight was removed (and a refusal above a budget), and the
hard clip, kept only so that the difference can be shown.

In koloa this is the quick, classical path. The mixture model of
koloa.noise is the one that pays for its outliers in the evidence.

Created on 2026-09-27

@author: artigau
"""
from typing import Dict, Optional, Tuple

import numpy as np

from koloa.log import loud_warning

# =============================================================================
# Define variables
# =============================================================================
#: the largest fraction of the total weight a clip may remove before koloa
#:   refuses to hide it
WEIGHT_BUDGET = 0.05


# =============================================================================
# Define functions
# =============================================================================
def soft_clip_factor(resid: np.ndarray, err: np.ndarray, clip: float = 3.0
                     ) -> np.ndarray:
    """
    The weight factor of each point: 1 inside the clip, (clip / z)^2 outside

    :param resid: np.ndarray, the residual of a model
    :param err: np.ndarray, the ORIGINAL error bars
    :param clip: float, where the weight starts to fall [robust sigma]

    :return: np.ndarray, a factor between 0 and 1 per point
    """
    factor = np.ones(len(resid))
    good = np.isfinite(resid) & np.isfinite(err) & (err > 0)
    if np.sum(good) < 5:
        return factor
    znorm = resid[good] / err[good]
    scale = 1.4826 * np.median(np.abs(znorm - np.median(znorm)))
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    zval = np.abs(znorm - np.median(znorm)) / scale
    factor[good] = np.where(zval > clip, (clip / np.maximum(zval, 1e-12)) ** 2,
                            1.0)
    return factor


def _trend_residual(time: np.ndarray, value: np.ndarray, err: np.ndarray,
                    degree: int = 1) -> np.ndarray:
    """
    The residual of a weighted polynomial in time

    :param time: np.ndarray, the time [days]
    :param value: np.ndarray, the value
    :param err: np.ndarray, the (current) error bars
    :param degree: int, the degree of the polynomial

    :return: np.ndarray, the residual
    """
    tnorm = (time - np.mean(time)) / max(np.ptp(time), 1e-9)
    coeffs = np.polyfit(tnorm, value, degree, w=1.0 / err)
    return value - np.polyval(coeffs, tnorm)


def soft_clip(time: np.ndarray, value: np.ndarray, err: np.ndarray,
              clip: float = 3.0, niter: int = 5, degree: int = 1,
              model: Optional[np.ndarray] = None,
              budget: Optional[float] = WEIGHT_BUDGET,
              quiet: bool = False) -> Dict[str, np.ndarray]:
    """
    Iterative soft clip, returned as inflated error bars with its accounting

    The residual is taken against a polynomial in time refitted at each pass
    with the weights of the pass before (so a handful of outliers cannot
    drag it onto themselves), or against a model given by the caller. The
    factor always comes from the original error bars.

    :param time: np.ndarray, the time [days]
    :param value: np.ndarray, the value
    :param err: np.ndarray, the error bars
    :param clip: float, where the weight starts to fall [sigma]
    :param niter: int, how many passes
    :param degree: int, the degree of the polynomial in time
    :param model: np.ndarray or None, a model to take the residual against
                  instead of the polynomial
    :param budget: float or None, the largest fraction of the weight that
                   may be removed before a loud warning (None: never warn)
    :param quiet: bool, no warning even above the budget

    :return: dict, err (inflated), factor, removed (fraction of the weight
             removed), ndown (points down-weighted), over_budget (bool)
    """
    current = np.array(err, dtype=float)
    factor = np.ones(len(value))
    for _ in range(max(niter, 1)):
        if model is None:
            resid = _trend_residual(time, value, current, degree)
        else:
            resid = value - model
        factor = soft_clip_factor(resid, err, clip=clip)
        current = err / np.sqrt(np.clip(factor, 1.0e-12, None))
    weight = 1.0 / err ** 2
    removed = float(np.sum(weight * (1 - factor)) / np.sum(weight))
    over = budget is not None and removed > budget
    if over and not quiet:
        loud_warning(f'The {clip:.1f}-sigma soft clip removes '
                     f'{100 * removed:.1f}% of the weight, above the budget '
                     f'of {100 * budget:.0f}%. A clip this strong decides '
                     f'the answer; report the result with and without it, '
                     f'or use the mixture model, which pays for its '
                     f'outliers.')
    return dict(err=current, factor=factor, removed=removed,
                ndown=int(np.sum(factor < 1)), over_budget=bool(over))


def hard_clip(time: np.ndarray, value: np.ndarray, err: np.ndarray,
              clip: float = 3.0, niter: int = 10, degree: int = 1
              ) -> np.ndarray:
    """
    The classical iterative sigma clip (points removed), for comparison only

    koloa never uses this for an answer: a removed point leaves the
    likelihood, so the evidence of whatever remains is overstated. It is
    here so that a demo or a paper can show what it does.

    :param time: np.ndarray, the time [days]
    :param value: np.ndarray, the value
    :param err: np.ndarray, the error bars
    :param clip: float, the threshold [robust sigma]
    :param niter: int, the largest number of passes
    :param degree: int, the degree of the polynomial in time

    :return: np.ndarray, boolean mask of the points kept
    """
    keep = np.ones(len(value), dtype=bool)
    for _ in range(niter):
        tnorm = (time - np.mean(time)) / max(np.ptp(time), 1e-9)
        coeffs = np.polyfit(tnorm[keep], value[keep], degree,
                            w=1.0 / err[keep])
        znorm = (value - np.polyval(coeffs, tnorm)) / err
        scale = 1.4826 * np.median(np.abs(znorm[keep]
                                          - np.median(znorm[keep])))
        new = np.abs(znorm - np.median(znorm[keep])) / scale < clip
        if np.array_equal(new, keep):
            break
        keep = new
    return keep


def student_weights(resid: np.ndarray, err: np.ndarray, dof: float = 4.0
                    ) -> np.ndarray:
    """
    The weights a Student-t likelihood gives each point at its mode

    (dof + 1) / (dof + z^2): the iteratively re-weighted least squares view
    of a heavy-tailed likelihood, close to the soft clip in spirit but
    derived from a density rather than chosen.

    :param resid: np.ndarray, the residual
    :param err: np.ndarray, the error bars (with any jitter)
    :param dof: float, the degrees of freedom

    :return: np.ndarray, the weight factor per point
    """
    znorm = resid / err
    return (dof + 1.0) / (dof + znorm ** 2)


def weight_summary(err: np.ndarray, factor: np.ndarray) -> Tuple[float, int]:
    """
    How much weight a factor removes, and from how many points

    :param err: np.ndarray, the original error bars
    :param factor: np.ndarray, the weight factor

    :return: tuple, the fraction of the weight removed, and the number of
             points with a factor below 1
    """
    weight = 1.0 / err ** 2
    removed = float(np.sum(weight * (1 - factor)) / np.sum(weight))
    return removed, int(np.sum(factor < 1))


# =============================================================================
# End of code
# =============================================================================

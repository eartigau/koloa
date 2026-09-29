#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The evidence of a linear model, and of that model plus one more sinusoid.

Offsets, trends, decorrelation coefficients and the two amplitudes of every
circular orbit enter the velocity linearly. With gaussian priors on them
and a gaussian noise (which the outlier indicators of koloa.noise keep
gaussian) they integrate out exactly (Hara et al. 2022, their appendix C):

    Z = N(y | 0, V + X Sigma X^T),
    log Z = -1/2 (y^T V^-1 y - b^T K^-1 b) - 1/2 log|V| - 1/2 log|Sigma|
            - 1/2 log|K| - n/2 log 2 pi,
    K = X^T V^-1 X + Sigma^-1,   b = X^T V^-1 y.

Adding two columns C = [cos, sin] with a prior variance tau^2 each changes
this by the Schur complement of the new block:

    Delta log Z = 1/2 u^T S^-1 u - 1/2 log|S| - log tau^2,
    S = C^T V^-1 C + I / tau^2 - B^T K^-1 B,   B = X^T V^-1 C,
    u = C^T V^-1 y - B^T K^-1 b,

which is a two by two matrix per trial period, so a whole grid of periods
and a whole set of tau values cost one matrix product.

Created on 2026-09-27

@author: artigau
"""
from typing import Optional, Tuple

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

from koloa.noise import LOG2PI, BlockCov


# =============================================================================
# Define classes
# =============================================================================
class LinearModel:
    """
    y = X beta + noise, beta ~ N(0, diag(prior_var)), noise ~ N(0, V)
    """

    def __init__(self, design: np.ndarray, value: np.ndarray, cov: BlockCov,
                 prior_var: np.ndarray, wdesign: Optional[np.ndarray] = None,
                 wvalue: Optional[np.ndarray] = None,
                 logdet_v: Optional[float] = None):
        """
        :param design: np.ndarray, X (n x p)
        :param value: np.ndarray, y (n)
        :param cov: BlockCov, the noise covariance V
        :param prior_var: np.ndarray, the prior variance of each column (p)
        :param wdesign: np.ndarray or None, V^-1 X if already known
        :param wvalue: np.ndarray or None, V^-1 y if already known
        :param logdet_v: float or None, log|V| if already known
        """
        self.design = design
        self.value = value
        self.prior_var = np.asarray(prior_var, dtype=float)
        self.wdesign = cov.solve(design) if wdesign is None else wdesign
        self.wvalue = cov.solve(value) if wvalue is None else wvalue
        self.logdet_v = cov.logdet() if logdet_v is None else logdet_v
        self.npar = design.shape[1]
        kmat = design.T @ self.wdesign
        kmat[np.diag_indices_from(kmat)] += 1.0 / self.prior_var
        self.kmat = kmat
        self.bvec = design.T @ self.wvalue
        self.chol = cho_factor(kmat, lower=True)
        self.mean = cho_solve(self.chol, self.bvec)
        self.logdet_k = 2.0 * float(np.sum(np.log(np.diag(self.chol[0]))))
        self.yvy = float(value @ self.wvalue)
        self.logz = float(-0.5 * (self.yvy - self.bvec @ self.mean)
                          - 0.5 * self.logdet_v
                          - 0.5 * np.sum(np.log(self.prior_var))
                          - 0.5 * self.logdet_k
                          - 0.5 * len(value) * LOG2PI)

    def kinv(self) -> np.ndarray:
        """
        K^-1, the posterior covariance of the linear parameters

        :return: np.ndarray, (p x p)
        """
        return cho_solve(self.chol, np.eye(self.npar))

    def draw(self, rng: np.random.Generator) -> np.ndarray:
        """
        One draw of the linear parameters from their posterior

        :param rng: np.random.Generator, the random numbers

        :return: np.ndarray, (p) the parameters
        """
        lower = np.tril(self.chol[0])
        # K = L L^T, so beta = mean + L^-T xi has covariance K^-1
        return self.mean + solve_triangular(lower.T, rng.normal(
            size=self.npar), lower=False)


# =============================================================================
# Define functions
# =============================================================================
def sinusoid_gain(model: LinearModel, bcos: np.ndarray, bsin: np.ndarray,
                  dcos: np.ndarray, dsin: np.ndarray, cc: np.ndarray,
                  ss: np.ndarray, cs: np.ndarray, tau2: np.ndarray
                  ) -> np.ndarray:
    """
    Delta log Z of adding a sinusoid, for every period and every tau

    :param model: LinearModel, the model without the sinusoid
    :param bcos: np.ndarray, (p x G) X^T V^-1 cos
    :param bsin: np.ndarray, (p x G) X^T V^-1 sin
    :param dcos: np.ndarray, (G) y^T V^-1 cos
    :param dsin: np.ndarray, (G) y^T V^-1 sin
    :param cc: np.ndarray, (G) cos^T V^-1 cos
    :param ss: np.ndarray, (G) sin^T V^-1 sin
    :param cs: np.ndarray, (G) cos^T V^-1 sin
    :param tau2: np.ndarray, (M) the prior variances of the amplitudes

    :return: np.ndarray, (G x M) the gain in log evidence
    """
    if model.npar > 0:
        mbc = cho_solve(model.chol, bcos)
        mbs = cho_solve(model.chol, bsin)
        s11 = cc - np.sum(bcos * mbc, axis=0)
        s22 = ss - np.sum(bsin * mbs, axis=0)
        s12 = cs - np.sum(bcos * mbs, axis=0)
        u1 = dcos - model.mean @ bcos
        u2 = dsin - model.mean @ bsin
    else:
        s11, s22, s12, u1, u2 = cc, ss, cs, dcos, dsin
    lam = 1.0 / np.asarray(tau2)[None, :]
    a11 = s11[:, None] + lam
    a22 = s22[:, None] + lam
    a12 = s12[:, None]
    det = a11 * a22 - a12 ** 2
    quad = (u1[:, None] ** 2 * a22 - 2 * u1[:, None] * u2[:, None] * a12
            + u2[:, None] ** 2 * a11) / det
    return 0.5 * quad - 0.5 * np.log(det) - np.log(np.asarray(tau2))[None, :]


def weighted_lstsq(design: np.ndarray, value: np.ndarray, err: np.ndarray
                   ) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Weighted least squares, with the covariance of the parameters

    :param design: np.ndarray, X (n x p)
    :param value: np.ndarray, y (n)
    :param err: np.ndarray, the error bars (n)

    :return: tuple, the parameters, their covariance, the chi2
    """
    wdesign = design / err[:, None]
    wvalue = value / err
    amat = wdesign.T @ wdesign
    coeffs = np.linalg.solve(amat, wdesign.T @ wvalue)
    cov = np.linalg.inv(amat)
    chi2 = float(np.sum((wvalue - wdesign @ coeffs) ** 2))
    return coeffs, cov, chi2


# =============================================================================
# End of code
# =============================================================================

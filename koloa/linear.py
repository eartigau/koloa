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
    return _gain(s11, s22, s12, u1, u2, tau2)


def _gain(s11: np.ndarray, s22: np.ndarray, s12: np.ndarray, u1: np.ndarray,
          u2: np.ndarray, tau2: np.ndarray) -> np.ndarray:
    """Delta log Z of a sinusoid from its Schur complement (s) and its
    projection on the residuals (u), for every period and every tau"""
    lam = 1.0 / np.asarray(tau2)[None, :]
    a11 = s11[:, None] + lam
    a22 = s22[:, None] + lam
    a12 = s12[:, None]
    det = a11 * a22 - a12 ** 2
    quad = (u1[:, None] ** 2 * a22 - 2 * u1[:, None] * u2[:, None] * a12
            + u2[:, None] ** 2 * a11) / det
    return 0.5 * quad - 0.5 * np.log(det) - np.log(np.asarray(tau2))[None, :]


class BaseProjection:
    """
    What the base columns (offsets, trends, regressors, the GP) take from
    every sinusoid of a grid, once for all the slots of a sweep

    With K = X^T V^-1 X + Sigma^-1 = L L^T, the base's part of B^T K^-1 B for
    a column c is |W c|^2, W = L^-1 X^T V^-1 (p x n), and its part of
    b^T K^-1 B is (W y)^T (W c). When the base has more columns than there
    are points (a GP of many bumps), a thin QR, W = Q R, leaves |W c| = |R c|
    with R (n x n): the grid costs min(p, n) n G, whatever the number of
    bumps, instead of p^2 G per slot.

    The other active slots of a model then enter as a small block (two
    columns each) after the base: their Schur complement S_O = K_OO -
    M^T M, M = W O, and the whitened products L_S^-1 (O^T V^-1 c -
    M^T W c).
    """

    def __init__(self, design: np.ndarray, wdesign: np.ndarray,
                 prior_var: np.ndarray, value: np.ndarray,
                 trig_cs: np.ndarray, size: int):
        """
        :param design: np.ndarray, X (n x p) the base columns
        :param wdesign: np.ndarray, V^-1 X
        :param prior_var: np.ndarray, (p) their prior variances
        :param value: np.ndarray, y (n)
        :param trig_cs: np.ndarray, (n x 2G) [cos | sin] of the grid
        :param size: int, G
        """
        npts, npar = design.shape
        if npar == 0:
            self.wmat = np.zeros((0, npts))
        else:
            kmat = design.T @ wdesign
            kmat[np.diag_indices_from(kmat)] += 1.0 / np.asarray(prior_var)
            chol = cho_factor(kmat, lower=True)[0]
            self.wmat = solve_triangular(chol, wdesign.T, lower=True)
            if npar > npts:
                self.wmat = np.linalg.qr(self.wmat, mode='r')
        zed = self.wmat @ trig_cs
        self.zc, self.zs = zed[:, :size], zed[:, size:]
        self.zy = self.wmat @ value
        self.qcc = np.einsum('ij,ij->j', self.zc, self.zc)
        self.qss = np.einsum('ij,ij->j', self.zs, self.zs)
        self.qcs = np.einsum('ij,ij->j', self.zc, self.zs)
        self.uc = self.zy @ self.zc
        self.us = self.zy @ self.zs

    def gain(self, cc: np.ndarray, ss: np.ndarray, cs: np.ndarray,
             dcos: np.ndarray, dsin: np.ndarray, tau2: np.ndarray,
             others: Optional[np.ndarray] = None,
             wothers: Optional[np.ndarray] = None,
             rothers: Optional[np.ndarray] = None,
             var_others: Optional[np.ndarray] = None,
             value: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Delta log Z of adding a sinusoid to the base and the other slots,
        for every period and every tau (sinusoid_gain, the same numbers)

        :param cc: np.ndarray, (G) cos^T V^-1 cos
        :param ss: np.ndarray, (G) sin^T V^-1 sin
        :param cs: np.ndarray, (G) cos^T V^-1 sin
        :param dcos: np.ndarray, (G) y^T V^-1 cos
        :param dsin: np.ndarray, (G) y^T V^-1 sin
        :param tau2: np.ndarray, (M) the prior variances of the amplitudes
        :param others: np.ndarray or None, O (n x q) the other slots
        :param wothers: np.ndarray or None, V^-1 O
        :param rothers: np.ndarray or None, (q x 2G) O^T V^-1 [cos | sin]
        :param var_others: np.ndarray or None, (q) their prior variances
        :param value: np.ndarray or None, y (with the others)

        :return: np.ndarray, (G x M) the gain in log evidence
        """
        s11, s22, s12 = cc - self.qcc, ss - self.qss, cs - self.qcs
        u1, u2 = dcos - self.uc, dsin - self.us
        if others is not None and others.shape[1] > 0:
            size = len(cc)
            mmat = self.wmat @ others
            kmat = others.T @ wothers - mmat.T @ mmat
            kmat[np.diag_indices_from(kmat)] += 1.0 / np.asarray(var_others)
            chol = cho_factor(kmat, lower=True)[0]
            wc = solve_triangular(chol, rothers[:, :size] - mmat.T @ self.zc,
                                  lower=True)
            ws = solve_triangular(chol, rothers[:, size:] - mmat.T @ self.zs,
                                  lower=True)
            wy = solve_triangular(chol, wothers.T @ value - mmat.T @ self.zy,
                                  lower=True)
            s11 = s11 - np.einsum('ij,ij->j', wc, wc)
            s22 = s22 - np.einsum('ij,ij->j', ws, ws)
            s12 = s12 - np.einsum('ij,ij->j', wc, ws)
            u1 = u1 - wy @ wc
            u2 = u2 - wy @ ws
        return _gain(s11, s22, s12, u1, u2, tau2)


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

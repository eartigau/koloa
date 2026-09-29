#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gaussian processes for stellar activity, dense and dependency-free.

A few hundred velocities make a covariance matrix that numpy factorises in
a millisecond, so there is no need for celerite here; a thousand points
cost about 30 ms. The kernels are the ones activity work uses:

- 'sho': a stochastically driven, damped harmonic oscillator (celerite's
  SHOTerm), period P, quality factor Q, amplitude sigma;
- 'rotation': two SHOs at P and P/2 (celerite2's RotationTerm), the usual
  model of a spotted rotating star;
- 'qp': the quasi-periodic kernel (Haywood et al. 2014), a periodic term
  times a squared exponential decay;
- 'se' and 'matern32': smooth, no period at all, as a null model.

What makes the GP outlier-aware is `loo_blocks`: the predictive of each
unit (point or sequence) given every other point. An outlier indicator
under a GP is not independent of the others, but its exact conditional
only needs this predictive (Sundararajan & Keerthi 2001), and a rank-one
update carries the inverse along when an indicator flips.

Created on 2026-09-27

@author: artigau
"""
from typing import List, Optional, Tuple

import numpy as np
from scipy.linalg import cho_factor, cho_solve

# =============================================================================
# Define variables
# =============================================================================
#: the parameters of each kernel, in the log of the natural units
KERNEL_PARAMS = {
    'sho': ['log_sigma', 'log_period', 'log_quality'],
    'rotation': ['log_sigma', 'log_period', 'log_q0', 'log_dq', 'logit_mix'],
    'qp': ['log_sigma', 'log_period', 'log_decay', 'log_smooth'],
    'se': ['log_sigma', 'log_length'],
    'matern32': ['log_sigma', 'log_length'],
}


# =============================================================================
# Kernels
# =============================================================================
def sho_kernel(tau: np.ndarray, sigma: float, period: float,
               quality: float) -> np.ndarray:
    """
    The covariance of a stochastically driven damped harmonic oscillator

    celerite's SHOTerm with S0 omega0 Q = sigma^2, written so that it is
    continuous through the critical damping Q = 1/2.

    :param tau: np.ndarray, the time lags [days]
    :param sigma: float, the standard deviation of the process
    :param period: float, the undamped period 2 pi / omega0 [days]
    :param quality: float, the quality factor Q

    :return: np.ndarray, the covariance
    """
    tau = np.abs(tau)
    omega0 = 2 * np.pi / period
    arg = omega0 * tau
    damp = arg / (2 * quality)
    if quality >= 0.5:
        eta = np.sqrt(max(1 - 1 / (4 * quality ** 2), 0.0))
        # sin(eta x) / eta, written as x sinc so that eta -> 0 is safe
        return sigma ** 2 * np.exp(-damp) * (
            np.cos(eta * arg) + damp * np.sinc(eta * arg / np.pi))
    eta = np.sqrt(1 / (4 * quality ** 2) - 1)
    earg = eta * arg
    # exp(-d) cosh(e) and exp(-d) sinh(e) / e without overflow (e < d)
    ecosh = 0.5 * (np.exp(earg - damp) + np.exp(-earg - damp))
    with np.errstate(invalid='ignore', divide='ignore'):
        esinhc = np.where(earg > 1e-8,
                          0.5 * (np.exp(earg - damp) - np.exp(-earg - damp))
                          / np.where(earg > 1e-8, earg, 1.0),
                          np.exp(-damp))
    return sigma ** 2 * (ecosh + damp * esinhc)


def kernel_matrix(kernel: str, tau: np.ndarray, pars: np.ndarray
                  ) -> np.ndarray:
    """
    The covariance of a kernel at the given lags

    :param kernel: str, sho, rotation, qp, se or matern32
    :param tau: np.ndarray, the time lags [days]
    :param pars: np.ndarray, the parameters, as in KERNEL_PARAMS

    :return: np.ndarray, the covariance
    """
    if kernel == 'sho':
        sigma, period, quality = np.exp(pars[:3])
        return sho_kernel(tau, sigma, period, quality)
    if kernel == 'rotation':
        # celerite2's RotationTerm: SHOs at P and P/2 sharing the variance
        sigma, period, q0, dq = np.exp(pars[:4])
        mix = 1 / (1 + np.exp(-pars[4]))
        q1, q2 = 0.5 + q0 + dq, 0.5 + q0
        amp1 = sigma / np.sqrt(1 + mix)
        amp2 = amp1 * np.sqrt(mix)
        # P is the period of the damped oscillation, as in celerite2: the
        #   undamped period of each oscillator is P sqrt(1 - 1 / 4Q^2)
        und1 = period * np.sqrt(1 - 1 / (4 * q1 ** 2))
        und2 = 0.5 * period * np.sqrt(1 - 1 / (4 * q2 ** 2))
        return (sho_kernel(tau, amp1, und1, q1)
                + sho_kernel(tau, amp2, und2, q2))
    if kernel == 'qp':
        sigma, period, decay, smooth = np.exp(pars[:4])
        return sigma ** 2 * np.exp(-0.5 * (tau / decay) ** 2
                                   - np.sin(np.pi * tau / period) ** 2
                                   / (2 * smooth ** 2))
    if kernel == 'se':
        sigma, length = np.exp(pars[:2])
        return sigma ** 2 * np.exp(-0.5 * (tau / length) ** 2)
    if kernel == 'matern32':
        sigma, length = np.exp(pars[:2])
        arg = np.sqrt(3) * np.abs(tau) / length
        return sigma ** 2 * (1 + arg) * np.exp(-arg)
    raise ValueError(f'Unknown kernel: {kernel}')


def kernel_period(kernel: str, pars: np.ndarray) -> float:
    """
    The period of a periodic kernel (nan for the others)

    :param kernel: str, the kernel
    :param pars: np.ndarray, its parameters

    :return: float, the period [days]
    """
    if kernel in ('sho', 'rotation', 'qp'):
        return float(np.exp(pars[1]))
    return np.nan


# =============================================================================
# The dense machinery
# =============================================================================
class MultiKernel:
    """
    Independent GPs on groups of points (one per instrument, or per set of
    instruments), each with its own kernel and parameters: the covariance
    is block diagonal, the points of two groups uncorrelated
    """

    def __init__(self, parts: List[Tuple[str, np.ndarray, np.ndarray]]):
        """
        :param parts: list of (kernel, pars, idx): the kernel, its
                      parameters (as in KERNEL_PARAMS) and the indices of
                      its points
        """
        self.parts = [(kernel, np.asarray(pars, dtype=float),
                       np.asarray(idx, dtype=int))
                      for kernel, pars, idx in parts]

    def matrix(self, time: np.ndarray) -> np.ndarray:
        """The covariance of the GPs at the points (n x n)"""
        out = np.zeros((len(time), len(time)))
        for kernel, pars, idx in self.parts:
            sub = time[idx]
            out[np.ix_(idx, idx)] = kernel_matrix(
                kernel, sub[:, None] - sub[None, :], pars)
        return out

    def cross(self, newtime: np.ndarray, time: np.ndarray,
              group: int) -> np.ndarray:
        """The covariance of one group's GP at new times with the points
        (m x n, zero for the points of the other groups)"""
        kernel, pars, idx = self.parts[group]
        out = np.zeros((len(newtime), len(time)))
        out[:, idx] = kernel_matrix(
            kernel, np.asarray(newtime)[:, None] - time[idx][None, :], pars)
        return out

    def prior(self, group: int) -> float:
        """The variance of one group's GP"""
        kernel, pars, _ = self.parts[group]
        return float(kernel_matrix(kernel, np.zeros(1), pars)[0])


class DenseGP:
    """
    A GP plus a noise covariance, factorised once

    The noise is the same diagonal plus blocks as koloa.noise.BlockCov,
    made dense here.
    """

    def __init__(self, time: np.ndarray, kernel: str, pars: np.ndarray,
                 diag: np.ndarray, block: np.ndarray, blockval: np.ndarray,
                 scale: Optional[np.ndarray] = None):
        """
        :param time: np.ndarray, the time [days]
        :param kernel: str, the kernel, or a MultiKernel (a GP per group
                       of points, its own parameters: pars is then None)
        :param pars: np.ndarray, its parameters
        :param diag: np.ndarray, the noise diagonal (n)
        :param block: np.ndarray, the sequence of each point (n)
        :param blockval: np.ndarray, the constant of each block (nb)
        :param scale: np.ndarray or None, the amplitude of the GP at each
                      point relative to the kernel's (one per instrument for
                      a chromatic GP): the covariance is
                      scale_i scale_j k(t_i - t_j)
        """
        self.time = time
        self.kernel, self.pars = kernel, pars
        self.scale = None if scale is None else np.asarray(scale, dtype=float)
        if isinstance(kernel, MultiKernel):
            self.kgp = kernel.matrix(time)
        else:
            self.kgp = kernel_matrix(kernel, time[:, None] - time[None, :],
                                     pars)
        if self.scale is not None:
            self.kgp = self.kgp * np.outer(self.scale, self.scale)
        cov = self.kgp.copy()
        cov[np.diag_indices_from(cov)] += diag
        if blockval is not None and np.any(blockval > 0):
            cov += blockval[block][:, None] * (block[:, None] == block[None, :])
        self.cov = cov
        self.chol = cho_factor(cov, lower=True, check_finite=False)

    def logdet(self) -> float:
        """log|C|"""
        return 2.0 * float(np.sum(np.log(np.diag(self.chol[0]))))

    def loglike(self, resid: np.ndarray) -> float:
        """
        log N(resid | 0, C)

        :param resid: np.ndarray, the residual of the mean model

        :return: float, the log likelihood
        """
        alpha = cho_solve(self.chol, resid, check_finite=False)
        return float(-0.5 * (resid @ alpha + self.logdet()
                             + len(resid) * np.log(2 * np.pi)))

    def predict_data(self, resid: np.ndarray
                     ) -> Tuple[np.ndarray, np.ndarray]:
        """
        The GP conditioned on the residual, at the points themselves (each
        with its own amplitude or group)

        :param resid: np.ndarray, the residual of the mean model

        :return: tuple, the mean and its one sigma
        """
        alpha = cho_solve(self.chol, resid, check_finite=False)
        mean = self.kgp @ alpha
        vsol = cho_solve(self.chol, self.kgp, check_finite=False)
        var = np.diag(self.kgp) - np.sum(self.kgp * vsol.T, axis=1)
        return mean, np.sqrt(np.clip(var, 0, None))

    def predict(self, resid: np.ndarray, newtime: np.ndarray,
                return_std: bool = True, scale_new=None,
                group: Optional[int] = None
                ) -> Tuple[np.ndarray, np.ndarray]:
        """
        The GP conditioned on the residual, at new times

        :param resid: np.ndarray, the residual of the mean model
        :param newtime: np.ndarray, where to predict [days]
        :param return_std: bool, also return the one sigma
        :param scale_new: float, np.ndarray or None, the relative amplitude
                          of the GP at the new times (1 when None: the
                          kernel's own)
        :param group: int or None, with a MultiKernel, the group whose GP
                      is predicted

        :return: tuple, the mean and its one sigma
        """
        if isinstance(self.kernel, MultiKernel):
            if group is None:
                raise ValueError('A GP per group of points: say which '
                                 'group to predict')
            cross = self.kernel.cross(newtime, self.time, group)
            alpha = cho_solve(self.chol, resid, check_finite=False)
            mean = cross @ alpha
            if not return_std:
                return mean, None
            vsol = cho_solve(self.chol, cross.T, check_finite=False)
            var = self.kernel.prior(group) - np.sum(cross.T * vsol, axis=0)
            return mean, np.sqrt(np.clip(var, 0, None))
        cross = kernel_matrix(self.kernel,
                              np.asarray(newtime)[:, None] - self.time[None, :],
                              self.pars)
        if self.scale is not None:
            cross = cross * self.scale[None, :]
        snew = 1.0 if scale_new is None else np.asarray(scale_new,
                                                          dtype=float)
        cross = cross * (snew[:, None] if np.ndim(snew) else snew)
        alpha = cho_solve(self.chol, resid, check_finite=False)
        mean = cross @ alpha
        if not return_std:
            return mean, None
        vsol = cho_solve(self.chol, cross.T, check_finite=False)
        prior = kernel_matrix(self.kernel, np.zeros(1), self.pars)[0] * \
            snew ** 2
        var = prior - np.sum(cross.T * vsol, axis=0)
        return mean, np.sqrt(np.clip(var, 0, None))

    def inverse(self) -> np.ndarray:
        """C^-1, dense"""
        return cho_solve(self.chol, np.eye(len(self.time)), check_finite=False)


def loo_blocks(cinv: np.ndarray, resid: np.ndarray, units: List[np.ndarray]
               ) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    The predictive of each unit given all the other points

    For a unit B, the conditional of y_B given y_-B has covariance
    [(C^-1)_BB]^-1 and mean r_B - [(C^-1)_BB]^-1 (C^-1 r)_B: no refit, only
    the blocks of the inverse (Sundararajan & Keerthi 2001, for one point).

    :param cinv: np.ndarray, C^-1 (n x n)
    :param resid: np.ndarray, the residual (n)
    :param units: list of np.ndarray, the indices of each unit

    :return: list of tuple, (mean, covariance) of each unit
    """
    alpha = cinv @ resid
    out = []
    for idx in units:
        sub = np.linalg.inv(cinv[np.ix_(idx, idx)])
        out.append((resid[idx] - sub @ alpha[idx], sub))
    return out


def sherman_morrison(cinv: np.ndarray, vec: np.ndarray, scale: float
                     ) -> np.ndarray:
    """
    (C + scale v v^T)^-1 from C^-1

    :param cinv: np.ndarray, C^-1
    :param vec: np.ndarray, v
    :param scale: float, the scale of the rank-one term

    :return: np.ndarray, the updated inverse
    """
    cv = cinv @ vec
    return cinv - scale * np.outer(cv, cv) / (1 + scale * vec @ cv)


# =============================================================================
# End of code
# =============================================================================

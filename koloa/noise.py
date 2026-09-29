#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The noise of a velocity series, and what an outlier does to it.

koloa's noise covariance is a diagonal plus one constant block per
observing sequence:

    V = diag(d) + sum_b a_b 1_b 1_b^T

The diagonal carries the photon noise, a jitter, and the inflation of a
point that is an outlier on its own. The blocks carry what the exposures of
one visit share: a sequence jitter, and the inflation of a whole sequence
that is an outlier (a bad night moves all its exposures together). Every
inverse and determinant follows from the Sherman-Morrison formula block by
block, so nothing of the size of the data squared is ever built.

The outlier model itself is a two-component mixture (Box & Tiao 1968): a
unit (point or sequence) is good with probability 1 - f, or an outlier with
probability f, in which case its variance grows by W^2. Given which units
are outliers, the noise stays gaussian, and this is what lets every linear
parameter be integrated out exactly.

Created on 2026-09-27

@author: artigau
"""
from typing import Optional, Tuple

import numpy as np
from scipy import sparse

# =============================================================================
# Define variables
# =============================================================================
LOG2PI = float(np.log(2 * np.pi))


# =============================================================================
# Define classes
# =============================================================================
class BlockCov:
    """
    V = diag(d) + sum_b a_b 1_b 1_b^T, with the algebra that goes with it

    For one block with diagonal D and constant a, Sherman-Morrison gives
        V^-1 = D^-1 - w D^-1 1 1^T D^-1,   w = a / (1 + a S),
        log|V| = log|D| + log(1 + a S),    S = sum 1/d.
    """

    def __init__(self, diag: np.ndarray, block: Optional[np.ndarray] = None,
                 blockval: Optional[np.ndarray] = None,
                 indicator: Optional[sparse.csr_matrix] = None):
        """
        :param diag: np.ndarray, the diagonal (n)
        :param block: np.ndarray or None, the block of each point (n)
        :param blockval: np.ndarray or None, the constant of each block (nb)
        :param indicator: sparse matrix or None, the (nb x n) membership,
                          passed in when it is already built
        """
        self.diag = np.asarray(diag, dtype=float)
        self.inv = 1.0 / self.diag
        self.block = block
        self.has_blocks = (block is not None and blockval is not None
                           and np.any(np.asarray(blockval) > 0))
        if self.has_blocks:
            self.blockval = np.asarray(blockval, dtype=float)
            nb = len(self.blockval)
            if indicator is None:
                indicator = block_indicator(block, nb)
            self.indicator = indicator
            self.ssum = np.bincount(block, weights=self.inv, minlength=nb)
            self.wblock = self.blockval / (1.0 + self.blockval * self.ssum)

    # -------------------------------------------------------------------------
    def logdet(self) -> float:
        """
        log|V|

        :return: float, the log determinant
        """
        value = float(np.sum(np.log(self.diag)))
        if self.has_blocks:
            value += float(np.sum(np.log1p(self.blockval * self.ssum)))
        return value

    def solve(self, x: np.ndarray) -> np.ndarray:
        """
        V^-1 x, for a vector or a (n x m) matrix

        :param x: np.ndarray, the right-hand side

        :return: np.ndarray, V^-1 x
        """
        x = np.asarray(x, dtype=float)
        if x.ndim == 1:
            out = x * self.inv
            if self.has_blocks:
                bsum = np.bincount(self.block, weights=out,
                                   minlength=len(self.blockval))
                out = out - self.inv * (self.wblock * bsum)[self.block]
            return out
        out = x * self.inv[:, None]
        if self.has_blocks:
            bsum = self.indicator @ out
            out = out - self.inv[:, None] * (self.wblock[:, None]
                                             * bsum)[self.block]
        return out

    def quad_diag(self, cols: np.ndarray, other: Optional[np.ndarray] = None
                  ) -> np.ndarray:
        """
        The diagonal terms c_g^T V^-1 c'_g for every column g of a matrix

        This is what the evidence of a new pair of sinusoids needs, for every
        period of a grid at once: sum_i c_i c'_i / d_i, minus the part the
        blocks take, sum_b w_b (sum_b c/d)(sum_b c'/d).

        :param cols: np.ndarray, (n x G) columns c
        :param other: np.ndarray or None, (n x G) columns c' (c when None)

        :return: np.ndarray, (G) the quadratic forms
        """
        if other is None:
            other = cols
        out = self.inv @ (cols * other)
        if self.has_blocks:
            bsum_c = self.indicator @ (cols * self.inv[:, None])
            if other is cols:
                bsum_o = bsum_c
            else:
                bsum_o = self.indicator @ (other * self.inv[:, None])
            out = out - np.sum(self.wblock[:, None] * bsum_c * bsum_o, axis=0)
        return out

    def scaled_indicator(self) -> sparse.csr_matrix:
        """
        The block membership with each entry scaled by 1/d_i

        (E D^-1) @ C gives the block sums of C / d without ever building
        C / d, which is the size of the whole grid.

        :return: sparse matrix, (nb x n)
        """
        scaled = self.indicator.copy()
        scaled.data = self.inv[scaled.indices].copy()
        return scaled

    def grid_quad(self, trig_sq: np.ndarray, trig_cs: np.ndarray, size: int
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        cos^T V^-1 cos, sin^T V^-1 sin and cos^T V^-1 sin over a grid

        The products cos^2, sin^2 and cos sin of the grid never change, only
        the noise does, so they are computed once (trig_sq) and reused here.

        :param trig_sq: np.ndarray, (n x 3G) [cos^2 | sin^2 | cos sin]
        :param trig_cs: np.ndarray, (n x 2G) [cos | sin]
        :param size: int, G

        :return: tuple, the three (G) quadratic forms
        """
        diag = self.inv @ trig_sq
        cc, ss, csum = diag[:size], diag[size:2 * size], diag[2 * size:]
        if self.has_blocks:
            bsum = self.scaled_indicator() @ trig_cs
            bc, bs = bsum[:, :size], bsum[:, size:]
            wcol = self.wblock[:, None]
            cc = cc - np.sum(wcol * bc * bc, axis=0)
            ss = ss - np.sum(wcol * bs * bs, axis=0)
            csum = csum - np.sum(wcol * bc * bs, axis=0)
        return cc, ss, csum

    def quad_diag_squares(self, c2: np.ndarray, s2: np.ndarray,
                          cs: np.ndarray, cos: np.ndarray, sin: np.ndarray
                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        cc, ss and cs of a grid of sinusoids, from separate matrices

        :param c2: np.ndarray, (n x G) cos^2
        :param s2: np.ndarray, (n x G) sin^2
        :param cs: np.ndarray, (n x G) cos sin
        :param cos: np.ndarray, (n x G) cos
        :param sin: np.ndarray, (n x G) sin

        :return: tuple, the three (G) quadratic forms
        """
        size = cos.shape[1]
        return self.grid_quad(np.hstack([c2, s2, cs]), np.hstack([cos, sin]),
                              size)

    def loglike(self, resid: np.ndarray) -> float:
        """
        log N(resid | 0, V)

        :param resid: np.ndarray, the residual

        :return: float, the log likelihood
        """
        quad = float(resid @ self.solve(resid))
        return -0.5 * (quad + self.logdet() + len(resid) * LOG2PI)


# =============================================================================
# Define functions
# =============================================================================
def block_indicator(block: np.ndarray, nblock: int) -> sparse.csr_matrix:
    """
    The (nblock x n) membership matrix of the blocks

    :param block: np.ndarray, the block of each point
    :param nblock: int, the number of blocks

    :return: sparse matrix, one where point i is in block b
    """
    npts = len(block)
    return sparse.csr_matrix((np.ones(npts), (block, np.arange(npts))),
                             shape=(nblock, npts))


def block_gauss_loglike(resid: np.ndarray, diag: np.ndarray,
                        block: np.ndarray, blockval: np.ndarray,
                        nblock: int) -> np.ndarray:
    """
    log N(r_b | 0, D_b + a_b 1 1^T), for every block at once

    :param resid: np.ndarray, the residual (n)
    :param diag: np.ndarray, the diagonal (n)
    :param block: np.ndarray, the block of each point (n)
    :param blockval: np.ndarray, the constant of each block (nb), or of
                     each point's block (n) when of length n
    :param nblock: int, the number of blocks

    :return: np.ndarray, (nb) the log likelihood of each block
    """
    inv = 1.0 / diag
    ssum = np.bincount(block, weights=inv, minlength=nblock)
    rsum = np.bincount(block, weights=resid * inv, minlength=nblock)
    r2sum = np.bincount(block, weights=resid ** 2 * inv, minlength=nblock)
    ldsum = np.bincount(block, weights=np.log(diag), minlength=nblock)
    npts = np.bincount(block, minlength=nblock)
    aval = np.asarray(blockval, dtype=float)
    quad = r2sum - aval * rsum ** 2 / (1.0 + aval * ssum)
    logdet = ldsum + np.log1p(aval * ssum)
    return -0.5 * (quad + logdet + npts * LOG2PI)


def mixture_unit_loglike(resid: np.ndarray, err: np.ndarray, block: np.ndarray,
                         nblock: int, jitter: float, seq_jitter: float,
                         frac: float, width: float, unit: str = 'point'
                         ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    The log likelihood of each unit under the good and the outlier model

    For unit = 'sequence' the unit is a block, and an outlier sequence has
    its whole block constant raised by W^2 (the exposures move together).
    For unit = 'point' the unit is a point, an outlier point has its own
    variance raised by W^2, and the sequence jitter (if any) couples the
    points of a block: the configurations of each block are then summed
    exactly, which blocks of a few exposures make cheap.

    :param resid: np.ndarray, the residual of the model (n)
    :param err: np.ndarray, the error bars (n)
    :param block: np.ndarray, the sequence of each point (n)
    :param nblock: int, the number of sequences
    :param jitter: float or np.ndarray, the white jitter (per point allowed)
    :param seq_jitter: float, the sequence jitter
    :param frac: float, the outlier fraction f
    :param width: float, the outlier width W
    :param unit: str, point or sequence

    :return: tuple, 1. the log likelihood of each block (mixture summed),
             2. the log probability that each unit is an outlier,
             3. the unit of each point (the block, or the point itself)
    """
    diag = err ** 2 + np.asarray(jitter) ** 2
    return mixture_loglike(resid, diag, block, nblock, seq_jitter ** 2, frac,
                           width, unit)


def mixture_loglike(resid: np.ndarray, diag: np.ndarray, block: np.ndarray,
                    nblock: int, seq_var: float, frac: float, width: float,
                    unit: str = 'point'
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    mixture_unit_loglike, with the diagonal of the noise given directly

    :param resid: np.ndarray, the residual of the model (n)
    :param diag: np.ndarray, the noise variance of each point (n)
    :param block: np.ndarray, the sequence of each point (n)
    :param nblock: int, the number of sequences
    :param seq_var: float or np.ndarray, the sequence variance (jitter
                    squared), or one per sequence (nblock)
    :param frac: float, the outlier fraction f
    :param width: float, the outlier width W
    :param unit: str, point or sequence

    :return: see mixture_unit_loglike
    """
    a0 = seq_var
    lf, l1f = np.log(frac), np.log1p(-frac)
    if unit == 'sequence':
        good = block_gauss_loglike(resid, diag, block,
                                   np.full(nblock, a0), nblock)
        bad = block_gauss_loglike(resid, diag, block,
                                  np.full(nblock, a0 + width ** 2), nblock)
        total = np.logaddexp(l1f + good, lf + bad)
        return total, lf + bad - total, block
    # point outliers without a sequence jitter: every point on its own
    if np.all(np.asarray(a0) <= 0):
        var_in, var_out = diag, diag + width ** 2
        good = -0.5 * (resid ** 2 / var_in + np.log(var_in) + LOG2PI)
        bad = -0.5 * (resid ** 2 / var_out + np.log(var_out) + LOG2PI)
        total = np.logaddexp(l1f + good, lf + bad)
        per_block = np.bincount(block, weights=total, minlength=nblock)
        return per_block, lf + bad - total, np.arange(len(resid))
    # point outliers with a sequence jitter: sum over the configurations
    return _point_mixture_blocks(resid, diag, block, nblock, a0, lf, l1f,
                                 width)


#: the grouping of points into blocks, by block size, for the arrays in use
_BLOCK_CACHE = {}
#: the longest sequence whose exposures can be outliers one by one: its 2^n
#: configurations are summed. A longer sequence is taken whole: all of its
#: exposures good (and, with sequence outliers, the sequence good or bad as
#: a whole), never an outlier exposure by exposure inside it
MAX_POINT_BLOCK = 10


def _block_groups(block: np.ndarray, nblock: int):
    """
    The blocks grouped by size, with the indices of their points and every
    outlier configuration of that size (cached: a fit calls this thousands
    of times with the same sequences)

    :return: dict, size: (block indices, (nblocks x size) point indices,
             (2^size x size) configurations, or (1 x size) for a sequence
             longer than MAX_POINT_BLOCK: only its all-good configuration)
    """
    key = (id(block), len(block), nblock)
    if key in _BLOCK_CACHE and _BLOCK_CACHE[key][0] is block:
        return _BLOCK_CACHE[key][1]
    sizes = np.bincount(block, minlength=nblock)
    order = np.argsort(block, kind='stable')
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    groups = {}
    for size in np.unique(sizes):
        if size == 0:
            continue
        blocks = np.where(sizes == size)[0]
        idx = order[starts[blocks][:, None] + np.arange(size)[None, :]]
        if size <= MAX_POINT_BLOCK:
            conf = ((np.arange(2 ** size)[:, None]
                     >> np.arange(size)[None, :]) & 1).astype(float)
        else:
            conf = np.zeros((1, size))
        groups[int(size)] = (blocks, idx, conf)
    if len(_BLOCK_CACHE) > 64:
        _BLOCK_CACHE.clear()
    _BLOCK_CACHE[key] = (block, groups)
    return groups


def _logsumexp(arr: np.ndarray, axis: int) -> np.ndarray:
    """log(sum(exp(arr))) along an axis, faster than logaddexp.reduce"""
    amax = np.max(arr, axis=axis, keepdims=True)
    amax = np.where(np.isfinite(amax), amax, 0.0)
    out = np.log(np.sum(np.exp(arr - amax), axis=axis, keepdims=True)) + amax
    return np.squeeze(out, axis=axis)


def _point_mixture_blocks(resid: np.ndarray, diag: np.ndarray,
                          block: np.ndarray, nblock: int, a0,
                          lf: float, l1f: float, width: float,
                          rng: Optional[np.random.Generator] = None):
    """
    Point outliers inside correlated blocks, all configurations summed

    With a random generator, one configuration per block is also drawn
    from its exact conditional, which is what a Gibbs sampler needs (the
    points of a block are not independent once a sequence term couples
    them). Blocks whose constant is zero are independent points and are
    done point by point.

    :param a0: float or np.ndarray, the constant of each block (nb)

    :return: see mixture_unit_loglike, plus the drawn configuration (n)
             when rng is given
    """
    npts = len(resid)
    a0 = np.broadcast_to(np.asarray(a0, dtype=float), (nblock,))
    per_block = np.zeros(nblock)
    logp_out = np.zeros(npts)
    drawn = np.zeros(npts)
    # the points of uncoupled blocks, one by one
    free = a0[block] <= 0
    if np.any(free):
        var_in, var_out = diag[free], diag[free] + width ** 2
        res = resid[free]
        good = l1f - 0.5 * (res ** 2 / var_in + np.log(var_in) + LOG2PI)
        bad = lf - 0.5 * (res ** 2 / var_out + np.log(var_out) + LOG2PI)
        total = np.logaddexp(good, bad)
        logp_out[free] = bad - total
        per_block += np.bincount(block[free], weights=total,
                                 minlength=nblock)
        if rng is not None:
            drawn[free] = (rng.random(np.sum(free))
                           < np.exp(bad - total)).astype(float)
    coupled = a0 > 0
    for size, (gblocks, gidx, conf) in _block_groups(block, nblock).items():
        sel_blocks = coupled[gblocks]
        if not np.any(sel_blocks):
            continue
        blocks = gblocks[sel_blocks]
        idx = gidx[sel_blocks]
        res = resid[idx]
        dia = diag[idx]
        ablk = a0[blocks][:, None]
        # (nblocks, nconf, size)
        dconf = dia[:, None, :] + conf[None, :, :] * width ** 2
        inv = 1.0 / dconf
        ssum = np.sum(inv, axis=2)
        rsum = np.sum(res[:, None, :] * inv, axis=2)
        r2sum = np.sum(res[:, None, :] ** 2 * inv, axis=2)
        quad = r2sum - ablk * rsum ** 2 / (1.0 + ablk * ssum)
        logdet = np.sum(np.log(dconf), axis=2) + np.log1p(ablk * ssum)
        nout = np.sum(conf, axis=1)
        logprior = nout * lf + (size - nout) * l1f
        if size > MAX_POINT_BLOCK:
            # a long sequence, taken whole: its one configuration is certain
            logprior = np.zeros(1)
        logl = -0.5 * (quad + logdet + size * LOG2PI) + logprior[None, :]
        total = _logsumexp(logl, axis=1)
        per_block[blocks] = total
        # the probability that each point is an outlier: the configurations
        #   where it is, over all of them
        prob = np.exp(logl - total[:, None])
        pout = prob @ conf
        logp_out[idx] = np.log(np.clip(pout, 1e-300, None))
        if rng is not None:
            cum = np.cumsum(prob, axis=1)
            pick = np.sum(cum < rng.random(len(blocks))[:, None]
                          * cum[:, -1:], axis=1)
            pick = np.clip(pick, 0, len(conf) - 1)
            drawn[idx] = conf[pick]
    if rng is not None:
        return per_block, logp_out, np.arange(npts), drawn
    return per_block, logp_out, np.arange(npts)


def mixture_loglike_both(resid: np.ndarray, diag: np.ndarray,
                         block: np.ndarray, nblock: int, seq_var: float,
                         frac_pt: float, width_pt: float, frac_seq: float,
                         width_seq: float):
    """
    Point outliers and sequence outliers at once, everything summed

    A point is an outlier with probability f_p (its variance grows by
    W_p^2), and independently a whole sequence is an outlier with
    probability f_s (its block constant grows by W_s^2). A lone spike and a
    bad night are then each explained by the right kind of outlier: three
    exposures moved together cost one sequence outlier, not three points.

    :param resid: np.ndarray, the residual (n)
    :param diag: np.ndarray, the noise variance of each point (n)
    :param block: np.ndarray, the sequence of each point (n)
    :param nblock: int, the number of sequences
    :param seq_var: float or np.ndarray, the sequence variance (jitter
                    squared), or one per sequence (nblock)
    :param frac_pt: float, the point outlier fraction
    :param width_pt: float, the point outlier width
    :param frac_seq: float, the sequence outlier fraction
    :param width_seq: float, the sequence outlier width

    :return: tuple, 1. the log likelihood of each block, 2. the log
             probability that each POINT is bad (an outlier itself or in
             an outlier sequence), 3. the log probability that each
             sequence is an outlier
    """
    lfp, l1fp = np.log(frac_pt), np.log1p(-frac_pt)
    lfs, l1fs = np.log(frac_seq), np.log1p(-frac_seq)
    good, lpt_good, _ = _point_mixture_blocks(
        resid, diag, block, nblock, np.full(nblock, seq_var), lfp, l1fp,
        width_pt)
    bad, lpt_bad, _ = _point_mixture_blocks(
        resid, diag, block, nblock, np.full(nblock, seq_var + width_seq ** 2),
        lfp, l1fp, width_pt)
    total = np.logaddexp(l1fs + good, lfs + bad)
    logp_seq = lfs + bad - total
    pseq = np.exp(logp_seq)[block]
    pbad = pseq + (1 - pseq) * np.exp(lpt_good)
    return total, np.log(np.clip(pbad, 1e-300, 1)), logp_seq


# =============================================================================
# End of code
# =============================================================================

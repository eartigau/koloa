#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The FIP machinery against brute force.

Every shortcut of koloa.fip (Sherman-Morrison blocks, Schur complements,
collapsed Gibbs steps, Rao-Blackwellised averages) is checked here against
the slow and obvious computation: dense matrices, and posteriors summed
over every configuration of signals and outliers.

Created on 2026-09-27

@author: artigau
"""
import itertools

import numpy as np
import pytest

from koloa import log as klog
from koloa.data import RVData
from koloa.fip import PeriodGrid, _default_setup, base_design, fip_single
from koloa.fip import oafip
from koloa.linear import LinearModel, sinusoid_gain
from koloa.noise import BlockCov, block_gauss_loglike, mixture_unit_loglike

klog.VERBOSE = False


# =============================================================================
# Helpers
# =============================================================================
def dense_logz(design, value, cov_dense, prior_var):
    """log N(y | 0, V + X Sigma X^T), the slow way"""
    full = cov_dense + design @ np.diag(prior_var) @ design.T
    sign, logdet = np.linalg.slogdet(full)
    quad = value @ np.linalg.solve(full, value)
    return -0.5 * (quad + logdet + len(value) * np.log(2 * np.pi))


def dense_cov(diag, block, blockval):
    """The dense V = diag(d) + sum_b a_b 1 1^T"""
    cov = np.diag(diag).astype(float)
    for bb, aval in enumerate(blockval):
        idx = np.where(block == bb)[0]
        cov[np.ix_(idx, idx)] += aval
    return cov


def toy_data(npts=40, nseq=None, seed=3, period=7.3, amp=4.0, noise=2.0):
    """A small series with one sinusoid"""
    rng = np.random.default_rng(seed)
    time = np.sort(rng.uniform(0, 120, npts))
    if nseq is not None:
        # visits of a few exposures
        starts = np.sort(rng.uniform(0, 120, nseq))
        time = np.concatenate([st + 0.01 * np.arange(npts // nseq)
                               for st in starts])
    err = noise * (0.8 + 0.4 * rng.random(len(time)))
    rv = amp * np.sin(2 * np.pi * time / period) + err * rng.normal(
        size=len(time))
    return RVData(time=time, rv=rv, err=err)


# =============================================================================
# The algebra
# =============================================================================
def test_blockcov_matches_dense():
    rng = np.random.default_rng(1)
    npts, nblock = 12, 4
    block = np.repeat(np.arange(nblock), 3)
    diag = rng.uniform(1, 3, npts)
    blockval = rng.uniform(0, 5, nblock)
    cov = BlockCov(diag, block, blockval)
    dense = dense_cov(diag, block, blockval)
    vec = rng.normal(size=npts)
    mat = rng.normal(size=(npts, 3))
    assert np.allclose(cov.solve(vec), np.linalg.solve(dense, vec))
    assert np.allclose(cov.solve(mat), np.linalg.solve(dense, mat))
    assert np.isclose(cov.logdet(), np.linalg.slogdet(dense)[1])
    # the block likelihoods sum to the dense one
    per_block = block_gauss_loglike(vec, diag, block, blockval, nblock)
    full = -0.5 * (vec @ np.linalg.solve(dense, vec)
                   + np.linalg.slogdet(dense)[1] + npts * np.log(2 * np.pi))
    assert np.isclose(np.sum(per_block), full)


def test_grid_quad_matches_dense():
    rng = np.random.default_rng(2)
    npts, nblock, size = 12, 4, 7
    block = np.repeat(np.arange(nblock), 3)
    diag = rng.uniform(1, 3, npts)
    blockval = rng.uniform(0, 5, nblock)
    cov = BlockCov(diag, block, blockval)
    dense = dense_cov(diag, block, blockval)
    cos, sin = rng.normal(size=(npts, size)), rng.normal(size=(npts, size))
    cc, ss, cs = cov.quad_diag_squares(cos * cos, sin * sin, cos * sin,
                                       cos, sin)
    inv = np.linalg.inv(dense)
    assert np.allclose(cc, np.einsum('ig,ij,jg->g', cos, inv, cos))
    assert np.allclose(ss, np.einsum('ig,ij,jg->g', sin, inv, sin))
    assert np.allclose(cs, np.einsum('ig,ij,jg->g', cos, inv, sin))


def test_linear_evidence_matches_dense():
    rng = np.random.default_rng(4)
    npts, nblock = 15, 5
    block = np.repeat(np.arange(nblock), 3)
    diag = rng.uniform(1, 3, npts)
    blockval = rng.uniform(0, 2, nblock)
    cov = BlockCov(diag, block, blockval)
    design = rng.normal(size=(npts, 3))
    prior_var = np.array([100.0, 4.0, 9.0])
    value = rng.normal(size=npts) * 3
    model = LinearModel(design, value, cov, prior_var)
    expected = dense_logz(design, value, dense_cov(diag, block, blockval),
                          prior_var)
    assert np.isclose(model.logz, expected)


def test_sinusoid_gain_matches_difference_of_evidences():
    rng = np.random.default_rng(5)
    npts, nblock = 18, 6
    block = np.repeat(np.arange(nblock), 3)
    diag = rng.uniform(1, 3, npts)
    blockval = rng.uniform(0, 2, nblock)
    cov = BlockCov(diag, block, blockval)
    time = np.sort(rng.uniform(0, 50, npts))
    value = rng.normal(size=npts) * 3
    base = np.stack([np.ones(npts), time / 50], axis=1)
    base_var = np.array([1e4, 1e4])
    model = LinearModel(base, value, cov, base_var)
    freq = np.array([0.05, 0.13, 0.31])
    cos = np.cos(2 * np.pi * np.outer(time, freq))
    sin = np.sin(2 * np.pi * np.outer(time, freq))
    tau2 = np.array([1.0, 10.0, 100.0])
    wdesign = model.wdesign
    cc, ss, cs = cov.quad_diag_squares(cos * cos, sin * sin, cos * sin,
                                       cos, sin)
    gain = sinusoid_gain(model, wdesign.T @ cos, wdesign.T @ sin,
                         model.wvalue @ cos, model.wvalue @ sin,
                         cc, ss, cs, tau2)
    for gg in range(len(freq)):
        for mm in range(len(tau2)):
            design = np.hstack([base, cos[:, [gg]], sin[:, [gg]]])
            var = np.concatenate([base_var, [tau2[mm]] * 2])
            full = LinearModel(design, value, cov, var)
            assert np.isclose(gain[gg, mm], full.logz - model.logz)


# =============================================================================
# The samplers against enumeration
# =============================================================================
def test_gibbs_with_one_slot_equals_single_signal_fip():
    data = toy_data()
    single = fip_single(data, jitter=0.0, pmin=2.0, oversample=5)
    gibbs = oafip(data, kmax=1, outliers=None, jitter=False, pmin=2.0,
                  oversample=5, nsweep=30, nburn=5, nchains=1,
                  progress=False, seq_jitter=False)
    # with one slot and the noise fixed, every Rao-Blackwellised term IS the
    #   single-signal computation
    assert np.allclose(gibbs.fip, single.fip, rtol=1e-8, atol=1e-12)


def _brute_force_tip(data, grid, tau, kmax, cov, trend=1):
    """TIP of every interval, every configuration of kmax slots summed"""
    base = base_design(data, trend)
    size, ntau = grid.size, len(tau)
    logk = np.log(np.full(kmax + 1, 1.0 / (kmax + 1)))
    from scipy.special import comb
    configs, logws = [], []
    for kval in range(kmax + 1):
        for active in itertools.combinations(range(kmax), kval):
            for gsel in itertools.product(range(size), repeat=kval):
                for tsel in itertools.product(range(ntau), repeat=kval):
                    cols = [base['design']]
                    var = [base['prior_var']]
                    for gg, tt in zip(gsel, tsel):
                        cols.append(grid.columns(gg))
                        var.append(np.full(2, tau[tt] ** 2))
                    model = LinearModel(np.hstack(cols), data.rv, cov,
                                        np.concatenate(var))
                    logw = (logk[kval] - np.log(comb(kmax, kval))
                            + sum(grid.logprior[gg] for gg in gsel)
                            - kval * np.log(ntau) + model.logz)
                    configs.append(gsel)
                    logws.append(logw)
    logws = np.array(logws)
    prob = np.exp(logws - np.logaddexp.reduce(logws))
    tip = np.zeros(size)
    for kk in range(size):
        inside = [any(abs(gg - kk) <= grid.halfbin for gg in gsel)
                  for gsel in configs]
        tip[kk] = np.sum(prob[np.array(inside)])
    return tip


def test_gibbs_two_slots_against_enumeration():
    data = toy_data(npts=30, seed=11, amp=3.0, noise=2.5)
    freq = np.linspace(0.02, 0.4, 9)
    grid = PeriodGrid(data.time, oversample=3, tref=data.tref, freq=freq)
    setup = _default_setup(data, 'none', False, False, None, 3, None)
    cov = BlockCov(data.err ** 2)
    tip = _brute_force_tip(data, grid, setup['tau'], 2, cov)
    gibbs = oafip(data, kmax=2, outliers=None, jitter=False, freq=freq,
                  oversample=3, ntau=3, nsweep=4000, nburn=200, nchains=1,
                  progress=False, seq_jitter=False)
    assert np.max(np.abs((1 - gibbs.fip) - tip)) < 0.02


def test_outlier_probabilities_against_enumeration():
    # eight points, one of them far off, point outliers with a fixed
    #   fraction and width, one possible signal on a small grid
    rng = np.random.default_rng(7)
    npts = 8
    time = np.sort(rng.uniform(0, 30, npts))
    err = np.full(npts, 1.0)
    rv = 2.0 * np.sin(2 * np.pi * time / 9.0) + rng.normal(size=npts)
    rv[3] += 15.0
    data = RVData(time=time, rv=rv, err=err)
    frac, width = 0.1, 8.0
    freq = np.linspace(0.05, 0.3, 6)
    grid = PeriodGrid(data.time, oversample=3, tref=data.tref, freq=freq)
    setup = _default_setup(data, 'point', False, False, None, 2, None)
    tau = setup['tau']
    base = base_design(data, 1)
    # enumeration over q (2^8) and over the signal (off, or a frequency
    #   and a rung)
    logws, qs, gsel_all = [], [], []
    for qvec in itertools.product([0, 1], repeat=npts):
        qvec = np.array(qvec, dtype=float)
        cov = BlockCov(err ** 2 + qvec * width ** 2)
        logq = np.sum(qvec) * np.log(frac) + (npts - np.sum(qvec)) * np.log(
            1 - frac)
        options = [(None, None)] + [(gg, tt) for gg in range(grid.size)
                                    for tt in range(len(tau))]
        for gg, tt in options:
            cols, var = [base['design']], [base['prior_var']]
            logp = np.log(0.5)
            if gg is not None:
                cols.append(grid.columns(gg))
                var.append(np.full(2, tau[tt] ** 2))
                logp += grid.logprior[gg] - np.log(len(tau))
            model = LinearModel(np.hstack(cols), rv, cov, np.concatenate(var))
            logws.append(logq + logp + model.logz)
            qs.append(qvec)
            gsel_all.append(gg)
    logws = np.array(logws)
    prob = np.exp(logws - np.logaddexp.reduce(logws))
    qprob = np.sum(prob[:, None] * np.array(qs), axis=0)
    tip = np.zeros(grid.size)
    for kk in range(grid.size):
        inside = np.array([gg is not None and abs(gg - kk) <= grid.halfbin
                           for gg in gsel_all])
        tip[kk] = np.sum(prob[inside])
    big = 1e7
    gibbs = oafip(data, kmax=1, outliers='point', jitter=False,
                  seq_jitter=False, freq=freq, oversample=3, ntau=2,
                  width_range=(width, width),
                  frac_prior=(frac * big, (1 - frac) * big), nsweep=6000,
                  nburn=300, nchains=1, progress=False)
    assert np.max(np.abs(gibbs.outlier_prob - qprob)) < 0.02
    assert np.max(np.abs((1 - gibbs.fip) - tip)) < 0.02
    # the planted outlier is found
    assert qprob[3] > 0.8


def test_sequence_mixture_counts_a_visit_once():
    # the same offset on the three exposures of a visit is one outlier
    #   sequence, not three outlier points
    rng = np.random.default_rng(8)
    nseq, per = 20, 3
    block = np.repeat(np.arange(nseq), per)
    err = np.ones(nseq * per)
    resid = rng.normal(size=nseq * per)
    resid[block == 4] += 15.0
    _, logp, unit = mixture_unit_loglike(resid, err, block, nseq, 0.0, 0.0,
                                         0.05, 15.0, 'sequence')
    assert len(logp) == nseq
    assert np.exp(logp[4]) > 0.99
    assert np.all(np.exp(np.delete(logp, 4)) < 0.1)


# =============================================================================
# End of code
# =============================================================================


# =============================================================================
# Both kinds of outliers at once
# =============================================================================
def test_mixture_both_matches_enumeration():
    from koloa.noise import mixture_loglike_both
    rng = np.random.default_rng(12)
    nseq, per = 5, 3
    block = np.repeat(np.arange(nseq), per)
    npts = len(block)
    diag = rng.uniform(1, 2, npts)
    resid = rng.normal(size=npts) * 2
    resid[4] += 12.0
    resid[block == 3] += 9.0
    fpt, wpt, fsq, wsq, svar = 0.08, 10.0, 0.1, 12.0, 0.5
    per_block, logp_pt, logp_seq = mixture_loglike_both(
        resid, diag, block, nseq, svar, fpt, wpt, fsq, wsq)
    for bb in range(nseq):
        idx = np.where(block == bb)[0]
        terms, bad_pt = [], np.zeros(per)
        seq_bad = 0.0
        for qseq in (0, 1):
            for conf in itertools.product([0, 1], repeat=per):
                conf = np.array(conf, dtype=float)
                cov = np.diag(diag[idx] + conf * wpt ** 2) + (
                    svar + qseq * wsq ** 2) * np.ones((per, per))
                sign, logdet = np.linalg.slogdet(cov)
                quad = resid[idx] @ np.linalg.solve(cov, resid[idx])
                logl = -0.5 * (quad + logdet + per * np.log(2 * np.pi))
                logl += (np.sum(conf) * np.log(fpt) + (per - np.sum(conf))
                         * np.log(1 - fpt))
                logl += np.log(fsq) if qseq else np.log(1 - fsq)
                terms.append((logl, qseq, conf))
        total = np.logaddexp.reduce([tt[0] for tt in terms])
        assert np.isclose(per_block[bb], total)
        for logl, qseq, conf in terms:
            prob = np.exp(logl - total)
            seq_bad += prob * qseq
            bad_pt += prob * np.maximum(conf, qseq)
        assert np.isclose(np.exp(logp_seq[bb]), seq_bad)
        assert np.allclose(np.exp(logp_pt[idx]), bad_pt)


def test_sampler_with_both_kinds_against_enumeration():
    # two visits of three exposures and two lone points; one visit moved
    #   together, one lone exposure spiked
    rng = np.random.default_rng(21)
    time = np.array([0.0, 0.01, 0.02, 7.0, 7.01, 7.02, 13.0, 21.0])
    npts = len(time)
    err = np.full(npts, 1.0)
    rv = 2.0 * np.sin(2 * np.pi * time / 11.0) + rng.normal(size=npts)
    rv[3:6] += 10.0
    rv[6] -= 12.0
    data = RVData(time=time, rv=rv, err=err)
    assert data.nseq == 4
    # the sampler has one prior shape for both fractions: the same here
    frac, width = 0.1, 10.0
    freq = np.linspace(0.04, 0.2, 5)
    grid = PeriodGrid(data.time, oversample=3, tref=data.tref, freq=freq)
    tau = _default_setup(data, 'both', False, False, None, 2, None)['tau']
    base = base_design(data, 1)
    options = [(None, None)] + [(g, t) for g in range(grid.size)
                                for t in range(len(tau))]
    logws, bads, gsel = [], [], []
    for qseq in itertools.product([0, 1], repeat=data.nseq):
        qseq = np.array(qseq, dtype=float)
        for qpt in itertools.product([0, 1], repeat=npts):
            qpt = np.array(qpt, dtype=float)
            cov = BlockCov(err ** 2 + qpt * width ** 2, data.seq,
                           qseq * width ** 2)
            nq = np.sum(qpt) + np.sum(qseq)
            logq = nq * np.log(frac) + (npts + data.nseq - nq) * np.log(
                1 - frac)
            bad = np.maximum(qpt, qseq[data.seq])
            for gg, tt in options:
                cols, var = [base['design']], [base['prior_var']]
                logp = np.log(0.5)
                if gg is not None:
                    cols.append(grid.columns(gg))
                    var.append(np.full(2, tau[tt] ** 2))
                    logp += grid.logprior[gg] - np.log(len(tau))
                model = LinearModel(np.hstack(cols), rv, cov,
                                    np.concatenate(var))
                logws.append(logq + logp + model.logz)
                bads.append(bad)
                gsel.append(gg)
    logws = np.array(logws)
    prob = np.exp(logws - np.logaddexp.reduce(logws))
    pbad = np.sum(prob[:, None] * np.array(bads), axis=0)
    inside = np.array([[gg is not None and abs(gg - kk) <= grid.halfbin
                        for gg in gsel] for kk in range(grid.size)])
    tip = inside.astype(float) @ prob
    big = 1e7
    gibbs = oafip(data, kmax=1, outliers='both', jitter=False,
                  seq_jitter=False, freq=freq, oversample=3, ntau=2,
                  width_range=(width, width),
                  frac_prior=(frac * big, (1 - frac) * big), nsweep=8000,
                  nburn=300, nchains=1, progress=False)
    assert np.max(np.abs(gibbs.outlier_prob - pbad)) < 0.02
    assert np.max(np.abs((1 - gibbs.fip) - tip)) < 0.02
    # the moved visit and the spike stand above the clean points (with eight
    #   points and a free sinusoid, the exact probabilities are 0.4 to 0.45)
    assert np.min(pbad[3:7]) > np.max(pbad[[0, 1, 2]])


def test_a_long_sequence_is_taken_whole():
    """a sequence of more than 10 exposures cannot hold point outliers one
    by one (2^n configurations): with point and sequence outliers, its
    likelihood is that of the sequence alone, good or bad as a whole"""
    from koloa.noise import mixture_loglike, mixture_loglike_both
    rng = np.random.default_rng(4)
    npts = 12
    resid = rng.normal(0, 1.5, npts)
    resid[5] += 9.0
    diag = np.full(npts, 1.0)
    block = np.zeros(npts, dtype=int)
    both, _, _ = mixture_loglike_both(resid, diag, block, 1, 0.4, 0.05, 20.0,
                                      0.1, 15.0)
    whole, _, _ = mixture_loglike(resid, diag, block, 1, 0.4, 0.1, 15.0,
                                  'sequence')
    assert np.allclose(both, whole)


def test_long_visits_do_not_stop_a_fit_or_a_fip():
    """visits of three exposures and one of fifteen: the fit and the FIP
    run, the long visit made bad is flagged as a whole, and a single bad
    exposure in a short visit is still flagged alone"""
    from koloa.fit import RVModel
    rng = np.random.default_rng(2)
    time, seq = [], []
    for night, nexp in enumerate([3] * 30 + [15] + [3] * 10):
        time.extend(10.0 * night + 0.01 * np.arange(nexp))
    time = np.array(time)
    rv = 3.0 * np.sin(2 * np.pi * time / 7.3) + rng.normal(0, 1.0, len(time))
    data = RVData(time, rv, np.full(len(time), 1.0))
    long_visit = np.where(np.bincount(data.seq)[data.seq] == 15)[0]
    data.rv[long_visit] += 12.0
    single = np.where(data.seq == 5)[0][1]
    data.rv[single] += 12.0
    fit = RVModel(data, [dict(period=7.3, period_range=(7.0, 7.6))],
                  likelihood='mixture', unit='both').fit(nstart=1, quiet=True)
    assert np.all(fit.outlier_prob[long_visit] > 0.9)
    assert fit.outlier_prob[single] > 0.9
    res = oafip(data, kmax=1, outliers='both', nsweep=30, nburn=10,
                nchains=1, progress=False)
    assert res.pk is not None

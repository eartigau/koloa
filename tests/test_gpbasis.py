#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The GP of the activity inside the FIP, as a finite basis

Created on 2026-09-30

@author: artigau
"""
import numpy as np
import pytest

from koloa import gpbasis
from koloa.fip import oafip
from koloa.simulate import simulate


def test_the_basis_is_the_kernel():
    time = np.sort(np.random.default_rng(1).uniform(0, 300, 150))
    for spec in ('local', dict(kind='rotation', period=(20.0, 30.0)),
                 ['local', dict(kind='rotation', period=(20.0, 30.0))]):
        comp = gpbasis.setup(spec, time, 5.0)
        val = gpbasis.initial(comp)
        cols, var = gpbasis.columns(comp, val, time)
        approx = (cols * var) @ cols.T
        exact = gpbasis.kernel(comp, val, time[:, None] - time[None, :])
        assert np.max(np.abs(approx - exact)) < 1e-3 * np.max(exact)


def test_priors():
    with pytest.raises(ValueError):
        gpbasis.setup('rotation', np.arange(10.0), 1.0)
    comp = gpbasis.setup(dict(kind='rotation', period=dict(
        mu=np.log(30.0), sd=0.1), length=40.0), np.arange(100.0), 1.0)
    assert comp[0]['priors']['length']['kind'] == 'fixed'
    assert 'gp0_length' not in gpbasis.names(comp)
    prior = comp[0]['priors']['period']
    assert gpbasis.log_prior(prior, 30.0) == 0.0
    assert gpbasis.log_prior(prior, 30.0 * np.exp(0.1)) == pytest.approx(-0.5)
    assert gpbasis.log_prior(prior, 1e4) == -np.inf


def test_fip_with_a_gp_samples_its_hyperparameters():
    sim = simulate(planets=[dict(P=5.3, K=5.0, e=0.0, tp=0.0)],
                   activity=dict(kernel='sho', sigma=4.0, period=60.0,
                                 quality=2.0), seed=4, err=1.5)
    res = oafip(sim['data'], kmax=2, nsweep=150, nburn=100, nchains=1,
                progress=False, gp='local')
    summ = res.gp_summary()
    assert set(summ) == {'gp0_length', 'gp0_sigma'}
    assert all(low <= mid <= high for mid, low, high in summ.values())
    assert 'GP local' in res.method
    assert res.fip_containing(5.3, 1 / sim['data'].baseline) < 0.5


def test_banded_fip_puts_the_bands_together():
    """two bands from the top down: a GP kept or not, every frequency of the
    periodogram from its own band"""
    from koloa.bandfip import banded_fip, fip_at
    sim = simulate(planets=[dict(P=5.3, K=6.0, e=0.0, tp=0.0)], seed=4,
                   err=1.5)
    res = banded_fip(sim['data'], nband=2, kmax=1, nsweep=120, nburn=80,
                     nchains=1, quiet=True)
    assert len(res['bands']) == 2 and len(res['freq']) == len(res['fip'])
    assert np.all(np.diff(res['freq']) > 0)
    assert fip_at(res, 5.3) < 0.5


def test_one_gp_per_instrument():
    """two instruments: two independent GPs, each zero on the other's
    points, with their own hyperparameters"""
    time = np.sort(np.random.default_rng(2).uniform(0, 200, 80))
    inst = np.where(np.arange(80) % 2 == 0, 'SPIRou', 'HARPS')
    comp = gpbasis.setup('local', time, 3.0, inst)
    assert [c['instrument'] for c in comp] == ['HARPS', 'SPIRou']
    assert len(gpbasis.names(comp)) == 4
    val = gpbasis.initial(comp)
    cols, var = gpbasis.columns(comp, val, time, inst)
    cov = (cols * var) @ cols.T
    # no covariance between the instruments, the SE kernel within each
    assert np.allclose(cov[np.ix_(inst == 'SPIRou', inst == 'HARPS')], 0)
    own = inst == 'SPIRou'
    exact = gpbasis.kernel([comp[1]], [val[1]],
                           time[own][:, None] - time[own][None, :])
    assert np.max(np.abs(cov[np.ix_(own, own)] - exact)) < 1e-3 * exact.max()


@pytest.mark.parametrize('npar', [6, 80])
def test_the_base_projection_is_the_full_model(npar):
    """the base projected once per sweep, the other slots as a small block
    after it: the same gain as the whole model solved for every slot, with
    fewer base columns than points and with more"""
    from koloa.linear import BaseProjection, LinearModel, sinusoid_gain
    from koloa.noise import BlockCov
    rng = np.random.default_rng(3)
    npts, size = 40, 300
    time = np.sort(rng.uniform(0, 200, npts))
    value = rng.normal(0, 3, npts)
    block = np.repeat(np.arange(20), 2)
    cov = BlockCov(rng.uniform(1, 2, npts), block, np.full(20, 0.5))
    design = rng.normal(size=(npts, npar))
    var = rng.uniform(0.5, 5, npar)
    freq = np.linspace(0.01, 0.5, size)
    arg = 2 * np.pi * time[:, None] * freq[None, :]
    trig_cs = np.hstack([np.cos(arg), np.sin(arg)])
    trig_sq = np.hstack([np.cos(arg) ** 2, np.sin(arg) ** 2,
                         np.cos(arg) * np.sin(arg)])
    others = trig_cs[:, [10, 10 + size, 200, 200 + size]]
    var_o = np.array([4.0, 4.0, 9.0, 9.0])
    tau2 = np.array([1.0, 10.0, 100.0])
    cc, ss, cs = cov.grid_quad(trig_sq, trig_cs, size)
    ry = cov.solve(value) @ trig_cs
    proj = BaseProjection(design, cov.solve(design), var, value, trig_cs,
                          size)
    for cols, cvar in ((design, var), (np.hstack([design, others]),
                                       np.concatenate([var, var_o]))):
        model = LinearModel(cols, value, cov, cvar)
        prod = model.wdesign.T @ trig_cs
        old = sinusoid_gain(model, prod[:, :size], prod[:, size:],
                            ry[:size], ry[size:], cc, ss, cs, tau2)
        if cols is design:
            new = proj.gain(cc, ss, cs, ry[:size], ry[size:], tau2)
        else:
            wothers = cov.solve(others)
            new = proj.gain(cc, ss, cs, ry[:size], ry[size:], tau2,
                            others=others, wothers=wothers,
                            rothers=wothers.T @ trig_cs, var_others=var_o,
                            value=value)
        assert np.allclose(new, old, rtol=1e-8, atol=1e-8)

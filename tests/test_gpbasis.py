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

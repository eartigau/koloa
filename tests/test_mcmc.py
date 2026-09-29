#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The MCMC of the orbits: autocorrelation times, the start of the walkers,
the maximum it starts from, the posterior in physical units, and the corner
plot that draws it.

Created on 2026-09-27

@author: artigau
"""
import os

import numpy as np
import pytest

from koloa import kepler
from koloa import log as klog
from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fit import FitResult, RVModel, integrated_time, mcmc_orbits
from koloa.simulate import simulate

klog.VERBOSE = False
#: a real NIRPS sampling: the times and error bars of 181 exposures in
#: 58 visits (its velocities are all zero)
TEMPLATE = os.path.join(os.path.dirname(__file__), '..', 'data',
                        'nirps_template.csv')


# =============================================================================
# Autocorrelation
# =============================================================================
def test_integrated_time_of_an_ar1_ensemble():
    # an AR(1) chain of coefficient phi has tau = (1 + phi) / (1 - phi)
    rng = np.random.default_rng(3)
    phi = 0.9
    chain = np.zeros((6000, 16, 1))
    for step in range(1, len(chain)):
        chain[step] = phi * chain[step - 1] + rng.normal(size=(16, 1))
    tau = integrated_time(chain)[0]
    assert abs(tau / ((1 + phi) / (1 - phi)) - 1) < 0.15


def test_integrated_time_matches_emcee():
    emcee = pytest.importorskip('emcee')
    rng = np.random.default_rng(5)
    chain = np.cumsum(rng.normal(size=(3000, 8, 2)), axis=0) * 0.01 + \
        rng.normal(size=(3000, 8, 2))
    ours = integrated_time(chain)
    theirs = emcee.autocorr.integrated_time(chain, quiet=True)
    assert np.allclose(ours, theirs, rtol=1e-8)


# =============================================================================
# The start
# =============================================================================
@pytest.fixture(scope='module')
def eccentric_series():
    # the case that exposed the jump: an eccentric planet on a real NIRPS
    #   sampling, whose +-10 % period range holds two aliases
    tpl = RVData.from_csv(TEMPLATE, name='NIRPS template')
    rng = np.random.default_rng(1001)
    period = 11.2
    tperi = float(tpl.time[0] + rng.uniform(0, period))
    sim = simulate(planets=[dict(P=period, K=10.0, e=0.3,
                                 omega=np.radians(60), tp=tperi)],
                   template=tpl, visit_jitter=2.0, seed=1,
                   outliers=[dict(kind='visit', frac=0.07, amplitude=6.0),
                             dict(kind='spike', frac=0.03, amplitude=8.0)])
    # the series without its outliers (the same noise)
    data = sim['data'].with_values(sim['data'].rv - sim['outlier_offset'])
    return data, tperi


def test_fit_stays_in_the_peak_it_starts_in(eccentric_series):
    # Powell's line search over the whole period range used to leap from
    #   the true period to an alias 38 lower in log posterior
    data, _ = eccentric_series
    model = RVModel(data, [dict(period=11.2, eccentric=True)],
                    likelihood='gaussian')
    fit = model.fit(nstart=4, quiet=True)
    orbit = fit.orbits()[0]
    assert abs(orbit['P'][0] - 11.2) < 0.02
    assert abs(orbit['e'][0] - 0.3) < 0.15
    # and the maximum is higher than anything in the alias peaks
    for prange in ((10.7, 11.1), (11.3, 11.7)):
        alias = RVModel(data, [dict(period=float(np.mean(prange)),
                                    eccentric=True, period_range=prange)],
                        likelihood='gaussian').fit(nstart=4, quiet=True)
        assert fit.logpost > alias.logpost


def test_walkers_are_distinct_even_with_a_broken_covariance(
        eccentric_series):
    data, _ = eccentric_series
    model = RVModel(data, [dict(period=11.2, eccentric=True)],
                    likelihood='gaussian')
    theta = model.init.copy()
    broken = np.diag(np.full(model.ndim, 1e12))
    start = FitResult(model, theta, model.log_posterior(theta), cov=broken)
    pos = model._ball(start, 40, np.random.default_rng(1))
    assert pos.shape == (40, model.ndim)
    assert all(np.isfinite(model.log_posterior(pp)) for pp in pos)
    # emcee refuses an ensemble whose spread is (nearly) singular
    spread = (pos - pos.mean(axis=0)) / pos.std(axis=0)
    assert np.linalg.cond(spread) < 1e8


# =============================================================================
# The posterior
# =============================================================================
def test_posterior_in_physical_units(eccentric_series):
    pytest.importorskip('emcee')
    data, tperi = eccentric_series
    res = mcmc_orbits(data, [11.2], likelihood='gaussian', nsteps=800,
                      nburn=200, max_steps=800, quiet=True)
    post = res.posterior(mstar=0.6)
    nn = len(res.chain)
    for key in ('P_0', 'K_0', 'e_0', 'omega_0', 'tp_0', 'tc_0', 'msini_0',
                'offset_inst', 'jit_inst', 'sjit'):
        assert len(post[key]) == nn
    assert np.all((post['omega_0'] >= 0) & (post['omega_0'] < 360))
    assert np.all((post['e_0'] >= 0) & (post['e_0'] < 1))
    assert np.all(post['jit_inst'] > 0)
    # orbits() summarises the same draws with their median
    orbit = res.orbits(mstar=0.6)[0]
    assert abs(orbit['K'][0] - np.median(post['K_0'])) < \
        0.05 * np.std(post['K_0'])
    # one draw, by hand
    theta = res.chain[0]
    period, tp, ecc, omega, amp = res.model.orbit(theta, 0)
    assert post['P_0'][0] == period and post['K_0'][0] == amp
    assert abs(post['tc_0'][0] - kepler.tp_to_tc(tp, period, ecc,
                                                 omega)) < 1e-9
    assert abs(post['msini_0'][0] - kepler.minimum_mass(amp, period, ecc,
                                                        0.6)) < 1e-9
    diag = res.diagnostics
    for key in ('sampler', 'acceptance', 'tau', 'tau_max', 'n_eff',
                'n_eff_min', 'converged'):
        assert key in diag
    assert set(diag['tau']) == set(res.model.names)


def test_circular_posterior_has_zero_eccentricity(eccentric_series):
    pytest.importorskip('emcee')
    data, _ = eccentric_series
    res = mcmc_orbits(data, [11.2], likelihood='gaussian', eccentric=False,
                      nsteps=400, nburn=100, max_steps=400, quiet=True)
    assert np.all(res.posterior()['e_0'] == 0.0)


def test_corner_draws_several_posteriors():
    rng = np.random.default_rng(2)
    posts = [dict(samples=dict(a=rng.normal(0, 1, 4000),
                               b=rng.normal(1, 2, 4000),
                               c=rng.normal(-1, 0.5, 4000)),
                  label=label, color=color, dashed=dashed)
             for label, color, dashed in (('one', 'gaussian', False),
                                          ('two', 'koloa', False),
                                          ('three', 'neutral', True))]
    fig = kplot.corner(posts, ['a', 'b', 'c'], truths=dict(a=0, b=1, c=-1))
    axes = np.array(fig.axes).reshape(3, 3)
    assert not axes[0, 1].get_visible() and axes[1, 0].get_visible()
    # every posterior draws two contours in each lower panel
    assert len(axes[1, 0].collections) >= 3
    kplot.plt.close(fig)


def test_corner_draws_a_posterior_in_the_panels_of_its_names():
    rng = np.random.default_rng(3)
    full = dict(samples=dict(a=rng.normal(0, 1, 4000),
                             b=rng.normal(1, 2, 4000)), label='full')
    part = dict(samples=dict(b=rng.normal(5, 1, 4000)), label='b only',
                color='gaussian', dashed=True)
    fig = kplot.corner([full, part], ['a', 'b'])
    axes = np.array(fig.axes).reshape(2, 2)
    # the diagonal of b holds both densities, that of a only the full one
    assert len(axes[1, 1].lines) == 2 and len(axes[0, 0].lines) == 1
    # the range of b covers both
    assert axes[1, 1].get_xlim()[1] > 7.0
    kplot.plt.close(fig)


# =============================================================================
# End of code
# =============================================================================

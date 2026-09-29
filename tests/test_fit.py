#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Orbits, GPs and outlier-aware fits.

Created on 2026-09-27

@author: artigau
"""
import itertools

import numpy as np
import pytest

from koloa import gp as kgp
from koloa import kepler
from koloa import log as klog
from koloa.data import RVData
from koloa.fit import RVModel
from koloa.simulate import simulate

klog.VERBOSE = False


# =============================================================================
# Orbits
# =============================================================================
def test_kepler_solver_is_accurate():
    rng = np.random.default_rng(1)
    manom = rng.uniform(-20, 20, 2000)
    for ecc in (0.0, 0.1, 0.5, 0.9, 0.99):
        eanom = kepler.eccentric_anomaly(manom, ecc)
        assert np.max(np.abs(eanom - ecc * np.sin(eanom) - manom)) < 1e-9


def test_keplerian_matches_radvel():
    radvel = pytest.importorskip('radvel')
    time = np.linspace(0, 300, 500)
    for ecc, omega in ((0.0, 0.0), (0.3, 1.1), (0.8, -2.0)):
        orbel = [23.4, 12.0, ecc, omega, 7.5]
        ours = kepler.rv_keplerian(time, 23.4, 12.0, ecc, omega, 7.5)
        theirs = radvel.kepler.rv_drive(time, np.array(orbel))
        assert np.max(np.abs(ours - theirs)) < 1e-6


def test_conjunction_matches_radvel():
    radvel = pytest.importorskip('radvel')
    for ecc, omega in ((0.1, 0.4), (0.5, 2.5), (0.0, 0.0)):
        ours = kepler.tc_to_tp(100.0, 12.3, ecc, omega)
        theirs = radvel.orbit.timetrans_to_timeperi(100.0, 12.3, ecc, omega)
        assert np.isclose(ours, theirs)
        assert np.isclose(kepler.tp_to_tc(ours, 12.3, ecc, omega), 100.0)


def test_basis_round_trip():
    orbit = (17.2, 55.0, 0.37, 1.9, 4.2)
    basis = kepler.orbit_to_basis(*orbit, tref=100.0)
    back = kepler.basis_to_orbit(*basis, tref=100.0)
    period, tperi, ecc, omega, amp = back
    assert np.isclose(period, orbit[0]) and np.isclose(amp, orbit[4])
    assert np.isclose(ecc, orbit[2]) and np.isclose(omega, orbit[3])
    # the periastron is the same modulo the period
    assert np.isclose(((tperi - orbit[1]) / orbit[0] + 0.5) % 1, 0.5)


def test_minimum_mass_of_the_earth():
    # the Earth moves the Sun by 8.95 cm/s over a year
    msini = kepler.minimum_mass(0.0895, 365.25, 0.0167, 1.0)
    assert abs(msini - 1.0) < 0.01


# =============================================================================
# Gaussian processes
# =============================================================================
def test_sho_matches_celerite2():
    celerite2 = pytest.importorskip('celerite2')
    tau = np.linspace(0, 40, 300)
    for quality in (0.3, 0.7, 3.0, 50.0):
        term = celerite2.terms.SHOTerm(sigma=2.5, rho=9.0, Q=quality)
        assert np.allclose(term.get_value(tau),
                           kgp.sho_kernel(tau, 2.5, 9.0, quality), atol=1e-10)


def test_gp_outlier_conditional_is_exact():
    # the leave-one-unit-out conditional of each indicator, against the
    #   ratio of two dense likelihoods
    rng = np.random.default_rng(3)
    npts = 7
    time = np.sort(rng.uniform(0, 30, npts))
    data = RVData(time=time, rv=rng.normal(size=npts) * 3,
                  err=np.full(npts, 1.0), seq=np.arange(npts))
    model = RVModel(data, [], gp=dict(kernel='sho'), likelihood='mixture',
                    unit='point', trend=0, seq_jitter=False)
    theta = model.init.copy()
    model.qunit = (rng.random(npts) < 0.3).astype(float)
    resid = data.rv - model.mean_model(theta)
    prob, _ = model._gp_outlier_conditionals(theta, resid, update=False)
    frac, width = model.outlier_params(theta)
    for uu in range(npts):
        logl = []
        for qval in (0.0, 1.0):
            qtest = model.qunit.copy()
            qtest[uu] = qval
            diag, blockval = model.gp_noise_blocks(theta, qtest)
            dense = kgp.DenseGP(time, 'sho', model.gp_pars(theta), diag,
                                data.seq, blockval)
            logl.append(dense.loglike(resid))
        lgood, lbad = np.log1p(-frac) + logl[0], np.log(frac) + logl[1]
        assert np.isclose(prob[uu], np.exp(lbad - np.logaddexp(lgood, lbad)))


def test_gp_both_conditionals_are_exact():
    # visits and exposures as outliers under a GP: each conditional against
    #   the ratio of two dense likelihoods, the other indicators held
    rng = np.random.default_rng(4)
    nvisit, per = 4, 2
    time = np.sort(np.repeat(rng.uniform(0, 40, nvisit), per)
                   + np.tile(np.arange(per) * 0.02, nvisit))
    seq = np.repeat(np.arange(nvisit), per)
    npts = len(time)
    data = RVData(time=time, rv=rng.normal(size=npts) * 3,
                  err=np.full(npts, 1.0), seq=seq)
    model = RVModel(data, [], gp=dict(kernel='sho'), likelihood='mixture',
                    unit='both', trend=0)
    assert model.unit == 'both'
    theta = model.init.copy()
    model.qunit = (rng.random(nvisit) < 0.4).astype(float)
    model.qpoint = (rng.random(npts) < 0.3).astype(float)
    resid = data.rv - model.mean_model(theta)
    probs, _ = model._gp_outlier_conditionals(theta, resid, update=False,
                                              by_kind=True)
    for kind, count, attr in (('sequence', nvisit, 'qunit'),
                              ('point', npts, 'qpoint')):
        frac, _ = model.outlier_params(theta, kind)
        for uu in range(count):
            logl = []
            for qval in (0.0, 1.0):
                qunit, qpoint = model.qunit.copy(), model.qpoint.copy()
                (qunit if attr == 'qunit' else qpoint)[uu] = qval
                diag, blockval = model.gp_noise_blocks(theta, qunit, qpoint)
                dense = kgp.DenseGP(time, 'sho', model.gp_pars(theta), diag,
                                    seq, blockval)
                logl.append(dense.loglike(resid))
            lgood = np.log1p(-frac) + logl[0]
            lbad = np.log(frac) + logl[1]
            assert np.isclose(probs[kind][uu],
                              np.exp(lbad - np.logaddexp(lgood, lbad)))


# =============================================================================
# Fits
# =============================================================================
@pytest.fixture(scope='module')
def planet_with_bad_visits():
    return simulate(planets=[dict(P=12.3, K=5.0, e=0.3, omega=1.0,
                                  tp=3.0)],
                    outliers=[dict(kind='visit', frac=0.08, amplitude=8.0)],
                    err=2.0, jitter=1.0, seed=5, nvisits=50, per_visit=3,
                    baseline=600)


def test_circular_fit_has_zero_eccentricity(planet_with_bad_visits):
    data = planet_with_bad_visits['data']
    for likelihood in ('gaussian', 'mixture'):
        fit = RVModel(data, [dict(period=12.3)],
                      likelihood=likelihood).fit(quiet=True)
        assert fit.orbits()[0]['e'][0] == 0.0
        assert all(val == 0.0 for val in
                   [fit.model.orbit(th, 0)[2] for th in fit.samples(50)])


def test_mixture_recovers_the_orbit_and_the_bad_visits(
        planet_with_bad_visits):
    sim = planet_with_bad_visits
    data = sim['data']
    fit = RVModel(data, [dict(period=12.25, eccentric=True)],
                  likelihood='mixture', unit='sequence').fit(quiet=True)
    orbit = fit.orbits()[0]
    assert abs(orbit['P'][0] - 12.3) < 0.05
    assert abs(orbit['K'][0] - 5.0) < 1.0
    assert abs(orbit['e'][0] - 0.3) < 0.15
    bad = sim['outlier_mask']
    assert np.mean(fit.outlier_prob[bad] > 0.5) > 0.8
    assert np.mean(fit.outlier_prob[~bad] < 0.5) > 0.95


def test_gaussian_fit_is_hurt_by_the_bad_visits(planet_with_bad_visits):
    data = planet_with_bad_visits['data']
    gauss = RVModel(data, [dict(period=12.25)],
                    likelihood='gaussian').fit(quiet=True)
    mix = RVModel(data, [dict(period=12.25)],
                  likelihood='mixture').fit(quiet=True)
    # the relative error on K is larger when the outliers are in the noise
    gk, mk = gauss.orbits()[0]['K'], mix.orbits()[0]['K']
    assert gk[1] > 1.5 * mk[1]


def two_instruments(sjit_a: float, sjit_b: float, seed: int = 3) -> RVData:
    """visits of three exposures (A) and of two (B), each instrument with its
    own visit jitter, over 1 m/s of photon noise and no signal"""
    rng = np.random.default_rng(seed)
    time, rv, inst = [], [], []
    for name, nexp, sjit, nvis in (('A', 3, sjit_a, 60), ('B', 2, sjit_b, 60)):
        for night in np.sort(rng.uniform(0, 800, nvis)):
            shift = rng.normal(0, sjit)
            for iexp in range(nexp):
                time.append(night + 0.01 * iexp)
                rv.append(shift + rng.normal(0, 1.0))
                inst.append(name)
    order = np.argsort(time)
    return RVData(np.array(time)[order], np.array(rv)[order],
                  np.ones(len(time)), inst=np.array(inst)[order])


def test_sequence_jitter_per_instrument_equals_one_when_equal():
    data = two_instruments(3.0, 3.0)
    one = RVModel(data, [], likelihood='mixture', unit='both', trend=0)
    per = RVModel(data, [], likelihood='mixture', unit='both', trend=0,
                  seq_jitter='instrument')
    assert 'log_sjit' in one.index and 'log_sjit' not in per.index
    assert per.seq_jitter_insts == ['A', 'B']
    theta_one = np.array(one.init, dtype=float)
    theta_per = np.array(per.init, dtype=float)
    for name in one.names:
        if name == 'log_sjit':
            for inst in ('A', 'B'):
                theta_per[per.index[f'log_sjit_{inst}']] = np.log(3.0)
            theta_one[one.index[name]] = np.log(3.0)
        else:
            theta_per[per.index[name]] = theta_one[one.index[name]]
    assert np.isclose(one.log_likelihood(theta_one),
                      per.log_likelihood(theta_per))
    assert np.allclose(one.outlier_probability(theta_one),
                       per.outlier_probability(theta_per))


def test_sequence_jitter_per_instrument_is_recovered():
    data = two_instruments(2.0, 7.0)
    fit = RVModel(data, [], likelihood='mixture', unit='both', trend=0,
                  seq_jitter='instrument').fit(quiet=True)
    sjit = {inst: float(np.exp(fit.theta[fit.model.index[f'log_sjit_{inst}']]))
            for inst in ('A', 'B')}
    assert abs(sjit['A'] - 2.0) < 0.8
    assert abs(sjit['B'] - 7.0) < 2.0
    post = fit.posterior(nsample=200)
    assert 'sjit_A' in post and 'sjit_B' in post


def test_inflate_to_fit_gives_each_instrument_the_noise_of_the_fit():
    from koloa.fip import inflate_to_fit
    data = two_instruments(2.0, 7.0)
    fit = RVModel(data, [], likelihood='mixture', unit='both', trend=0,
                  seq_jitter='instrument').fit(quiet=True)
    inflated, info = inflate_to_fit(fit)
    sjit, jit = info['visit_jitters'], info['jitters']
    # the smaller visit jitter (A) is the reference, left to the FIP
    assert info['visit_jitter_ref'] == min(sjit.values()) == sjit['A']
    # B, two exposures per visit: 2 (s_B^2 - s_A^2) + j_B^2 added
    expect_b = np.sqrt(2 * (sjit['B'] ** 2 - sjit['A'] ** 2) + jit['B'] ** 2)
    assert np.isclose(info['inflation']['B'], expect_b)
    sel = inflated.inst == 'B'
    assert np.allclose(inflated.err[sel] ** 2,
                       data.err[sel] ** 2 + expect_b ** 2)
    # the velocities, instruments and visits are those of the series
    assert np.allclose(inflated.rv, data.rv)
    assert np.array_equal(inflated.seq, data.seq)


def test_chromatic_gp_equals_the_gp_when_the_scales_are_one():
    data = two_instruments(2.0, 2.0, seed=4)
    gp = dict(kernel='matern32')
    plain = RVModel(data, [], likelihood='mixture', unit='both', trend=0,
                    gp=dict(gp))
    chrom = RVModel(data, [], likelihood='mixture', unit='both', trend=0,
                    gp=dict(gp, scale='instrument', reference='A'))
    assert chrom.gp_scale_insts == ['B']
    assert 'gp_log_scale_B' in chrom.index
    theta_c = np.array(chrom.init, dtype=float)
    theta_p = np.array([theta_c[chrom.index[name]] for name in plain.names])
    theta_c[chrom.index['gp_log_scale_B']] = 0.0
    assert np.isclose(plain.log_likelihood(theta_p),
                      chrom.log_likelihood(theta_c))
    # a scale of 2 multiplies the covariance of B by 4 and its cross terms
    #   with A by 2
    theta_c[chrom.index['gp_log_scale_B']] = np.log(2.0)
    scale = chrom.gp_scale(theta_c)
    assert np.allclose(np.unique(scale), [1.0, 2.0])
    diag, blockval = chrom.gp_noise_blocks(theta_c, chrom.qunit, chrom.qpoint)
    dgp = kgp.DenseGP(data.time, 'matern32', chrom.gp_pars(theta_c), diag,
                      data.seq, blockval, scale=scale)
    base = kgp.kernel_matrix('matern32', data.time[:, None]
                             - data.time[None, :], chrom.gp_pars(theta_c))
    assert np.allclose(dgp.kgp, base * np.outer(scale, scale))
    # the prediction at the data, each point with its own amplitude
    resid = data.rv - data.rv.mean()
    mean, _ = dgp.predict(resid, data.time, scale_new=scale)
    assert np.allclose(mean, dgp.kgp @ dgp.inverse() @ resid)


def test_gp_groups_reduce_to_the_gp_and_are_independent():
    data = two_instruments(2.0, 4.0, seed=6)
    plain = RVModel(data, [], likelihood='gaussian', unit='both', trend=0,
                    gp=dict(kernel='matern32'))
    one = RVModel(data, [], likelihood='gaussian', unit='both', trend=0,
                  gp=[dict(kernel='matern32', instruments=['A', 'B'],
                           name='all')])
    # one group over everything is the GP itself
    theta_p = np.array(plain.init, dtype=float)
    theta_o = np.array(one.init, dtype=float)
    for name in plain.names:
        other = name.replace('gp_', 'gp_all_') if name.startswith('gp_') \
            else name
        theta_o[one.index[other]] = theta_p[plain.index[name]]
    assert np.isclose(plain.log_likelihood(theta_p),
                      one.log_likelihood(theta_o))
    # two groups: a block-diagonal covariance, each block its own kernel
    two = RVModel(data, [], likelihood='gaussian', unit='both', trend=0,
                  gp=[dict(kernel='matern32', instruments=['A']),
                      dict(kernel='matern32', instruments=['B'])])
    assert [group['name'] for group in two.gp_groups] == ['A', 'B']
    theta = np.array(two.init, dtype=float)
    theta[two.index['gp_A_log_sigma']] = np.log(3.0)
    theta[two.index['gp_B_log_sigma']] = np.log(5.0)
    diag, blockval = two.gp_noise_blocks(theta, two.qunit, two.qpoint)
    dgp = two.dense_gp(theta, diag, data.seq, blockval)
    is_a = data.inst.astype(str) == 'A'
    assert np.allclose(dgp.kgp[np.ix_(is_a, ~is_a)], 0.0)
    assert np.isclose(dgp.kgp[is_a, is_a][0], 9.0)
    assert np.isclose(dgp.kgp[~is_a, ~is_a][0], 25.0)
    # the prediction of a group at its own points is that at the data
    resid = data.rv - data.rv.mean()
    at_data, _ = dgp.predict_data(resid)
    mean_b, _ = dgp.predict(resid, data.time[~is_a], group=1)
    assert np.allclose(mean_b, at_data[~is_a])


# =============================================================================
# End of code
# =============================================================================


def test_laplace_errors_survive_parameters_at_their_bounds():
    """a maximum a posteriori with parameters at a bound of their priors (a
    jitter at zero, say) still gives finite errors: the Laplace draws that
    fall outside are brought inside the bounds when too few survive"""
    data = two_instruments(2.0, 2.0)
    data = data.with_values(data.rv + 3.0 * np.sin(2 * np.pi * data.time
                                                   / 17.0))
    fit = RVModel(data, [dict(period=17.0)], likelihood='gaussian', trend=0,
                  seq_jitter='instrument').fit(nstart=1, quiet=True)
    # every jitter at the lower bound of its prior, with a covariance far
    #   wider than the prior: nearly every draw would be rejected
    for name in fit.model.names:
        if name.startswith(('log_jit', 'log_sjit')):
            it = fit.model.index[name]
            fit.theta[it] = fit.model.priors[it].a
            fit.cov[it, :] = fit.cov[:, it] = 0.0
            fit.cov[it, it] = 1e4
    draws = fit.samples(500)
    assert len(draws) >= 100
    assert all(np.isfinite(fit.model.log_prior(dr)) for dr in draws)
    kval = fit.orbits()[0]['K']
    assert np.all(np.isfinite(kval)) and kval[1] > 0 and kval[2] > 0

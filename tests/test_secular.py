#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The acceleration of a star (a trend read in m/s/yr, and what it says of a
companion) and the perspective acceleration: values, errors, and their
place in a fit.

Created on 2026-09-27

@author: artigau
"""
import os

import numpy as np

from koloa import log as klog
from koloa.data import RVData
from koloa.fit import RVModel
from koloa.secular import (ACCEL_FACTOR_MIN, ACCEL_MASS, JULIAN_YEAR, MAS,
                           PARSEC, acceleration, astrometry_covariance,
                           companion_acceleration, companion_min_mass,
                           perspective_from_astrometry, remove_secular,
                           secular_acceleration, secular_drift,
                           secular_from_astrometry)
from koloa.simulate import REALISTIC, simulate

klog.VERBOSE = False
#: a real NIRPS sampling: the times and error bars of 181 exposures in
#: 58 visits (its velocities are all zero)
TEMPLATE = os.path.join(os.path.dirname(__file__), '..', 'data',
                        'nirps_template.csv')
#: Barnard's star (GJ 699), Gaia DR3 4472832130942575872
BARNARD = dict(pmra=-801.5509783684709, pmra_error=0.031820867,
               pmdec=10362.394206546573, pmdec_error=0.036070455,
               parallax=546.975939730948, parallax_error=0.040116355,
               parallax_pmra_corr=-0.05781335,
               parallax_pmdec_corr=0.5593565, pmra_pmdec_corr=0.03058501)


def test_barnards_star():
    # 4.5 m/s/yr (Kurster et al. 2003, with the Hipparcos astrometry)
    value, error = secular_from_astrometry(BARNARD)
    assert abs(value - 4.50) < 0.05
    assert 0 < error < 1e-3


def test_constant_from_first_principles():
    # v_t^2 / d, with v_t = mu d, in SI, for 1 arcsec/yr at 1 pc
    mu = 1000 * MAS / JULIAN_YEAR
    accel = mu ** 2 * PARSEC * JULIAN_YEAR
    assert abs(secular_acceleration(1000.0, 0.0, 1000.0)[0] - accel) < 1e-12


def test_error_propagation_matches_draws():
    astro = dict(BARNARD, parallax=20.0, parallax_error=0.4,
                 pmra=300.0, pmra_error=0.5, pmdec=-200.0, pmdec_error=0.5)
    cov = astrometry_covariance(astro)
    _, first = secular_acceleration(300.0, -200.0, 20.0, cov)
    _, drawn = secular_acceleration(300.0, -200.0, 20.0, cov, ndraw=200000)
    assert abs(first / drawn - 1) < 0.03


def test_fit_measures_and_uses_the_acceleration():
    tpl = RVData.from_csv(TEMPLATE, name='NIRPS template')
    accel = 4.54
    sim = simulate(template=tpl, visit_jitter=1.0, seed=4)
    data = sim['data']
    data = data.with_values(data.rv + secular_drift(data.time, accel,
                                                    data.tref))
    # the acceleration measured from the velocities alone
    free = RVModel(data, [], trend=0, secular=(0.0, None),
                   likelihood='gaussian').fit(nstart=1, quiet=True)
    val, low, high = free.param('secacc')
    assert abs(val - accel) < 3 * max(low, high)
    # with the prior, a trend fitted at the same time finds nothing left
    both = RVModel(data, [], trend=1, secular=(accel, 0.001),
                   likelihood='gaussian').fit(nstart=1, quiet=True)
    trend = both.theta[both.model.index['trend_1']]
    slope = trend / data.baseline * 365.25          # m/s/yr
    assert abs(slope) < 0.5
    # and taking the drift out by hand leaves a flat series
    flat = remove_secular(data, accel)
    coef = np.polyfit(flat.time - flat.tref, flat.rv, 1)[0] * 365.25
    assert abs(coef) < 1.0


def test_old_names_are_the_perspective_functions():
    assert secular_from_astrometry(BARNARD) == \
        perspective_from_astrometry(BARNARD)


def test_torres_constant_and_bound():
    # au^2 / (yr G Msun) = 5.34e-6 Msun (Torres 1999; Liu et al. 2002, Eq. 1)
    assert abs(ACCEL_MASS / 5.34e-6 - 1) < 1e-3
    assert abs(ACCEL_FACTOR_MIN - 2.598) < 1e-3
    # HR 7672: 0.79 arcsec at 17.7 pc, 24 m/s/yr (Liu et al. 2002)
    mmin = companion_min_mass(24.0, 17.7, 0.79)
    assert abs(mmin - 0.0651) < 5e-4
    # the bound is reached at depth rho / sqrt(2), and never exceeded
    assert abs(companion_acceleration(mmin, 17.7, 0.79) - 24.0) < 1e-9
    rho = 17.7 * 0.79
    for depth in np.linspace(-5 * rho, 5 * rho, 101):
        assert companion_acceleration(mmin, 17.7, 0.79, depth) <= 24.0 + 1e-9


def test_acceleration_with_outliers():
    # a companion drift and its change, under bad visits and spikes: koloa's
    #   likelihood finds both within their errors, in m/s/yr and m/s/yr^2
    tpl = RVData.from_csv(TEMPLATE, name='NIRPS template')
    accel, jerk = 3.0, 1.5
    sim = simulate(template=tpl, visit_jitter=1.0, outliers=REALISTIC,
                   seed=6)
    data = sim['data']
    dt = (data.time - data.tref) / 365.25
    data = data.with_values(data.rv + accel * dt + 0.5 * jerk * dt ** 2)
    fit = RVModel(data, [], trend=2, likelihood='mixture',
                  unit='both').fit(nstart=2, quiet=True)
    acc = acceleration(fit)
    val, low, high = acc['accel']
    assert abs(val - accel) < 3 * max(low, high)
    val, low, high = acc['jerk']
    assert abs(val - jerk) < 3 * max(low, high)
    assert acc['draws'] is None
    # a perspective term fitted beside it is reported apart, and in the sum
    both = RVModel(data, [], trend=2, perspective=(0.5, 0.01),
                   likelihood='mixture', unit='both').fit(nstart=2,
                                                          quiet=True)
    acc = acceleration(both)
    assert abs(acc['perspective'][0] - 0.5) < 0.05
    assert abs(acc['total'][0] - accel) < 3 * max(acc['total'][1:])


# =============================================================================
# End of code
# =============================================================================

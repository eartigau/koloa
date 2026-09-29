#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The acceleration of a star, with its error.

"Secular acceleration" has two senses. The acceleration of the star itself,
from a companion too far out to complete an orbit, is measured: it is the
trend of a fit, which acceleration() reads in m/s/yr (and in m/s/yr^2 for a
trend of degree two) with its errors, kept honest under outliers by koloa's
likelihood; companion_min_mass() gives the smallest companion that makes it
at a given separation (Torres 1999). The perspective acceleration, mu^2 d,
is geometry, known from Gaia with its error; APERO (SPIRou, NIRPS) already
removes it in its barycentric correction, so the last part adds that of
Barnard's star to a simulated series, as velocities that still hold it
would carry it, and fits it beside the acceleration of the star. The query
of SIMBAD and the Gaia archive needs the network.

Created on 2026-09-27

@author: artigau
"""
# koloa's series, model, log (timestamped lines, 'value' for a number),
#   acceleration tools (GM_SUN and GM_JUP in m^3/s^2) and simulation
from pathlib import Path

from koloa.data import RVData
from koloa.fit import RVModel
from koloa.log import log
from koloa.secular import (GM_JUP, GM_SUN, acceleration, companion_min_mass,
                           gaia_astrometry, perspective_acceleration,
                           perspective_drift, perspective_from_astrometry)
from koloa.simulate import REALISTIC, simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]
# a real NIRPS sampling: its times, visits and error bars for the
#   simulation
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')

# 1. a companion far out accelerates the star by 2 m/s/yr, on a series with
#   a visit jitter and bad visits and spikes
# (no planet; REALISTIC: whole visits and single exposures made outliers,
#   clear ones and borderline ones)
accel = 2.0
sim = simulate(template=tpl, visit_jitter=2.0, outliers=REALISTIC, seed=3)
data = sim['data']
# the drift of a constant acceleration, accel (t - tref) / 365.25 [m/s],
#   added to the velocities
data = data.with_values(data.rv + perspective_drift(data.time, accel,
                                                    data.tref))
for likelihood in ('gaussian', 'mixture'):
    # no orbit and a straight trend (trend=1), without or with outliers
    #   (unit='both': a lone exposure or a whole visit); nstart=2: the
    #   number of starts (4 by default)
    fit = RVModel(data, [], trend=1, likelihood=likelihood,
                  unit='both').fit(nstart=2, quiet=True)
    # the slope of the trend in m/s/yr, as (value, minus, plus): from the
    #   curvature at the maximum here, from the chain for a sampled fit
    val, low, high = acceleration(fit)['accel']
    log(f'{likelihood:8s} acceleration {val:5.2f} +- {low:.2f} m/s/yr '
        f'(true {accel})', 'value')

# what it asks of a companion 2 arcsec away at the distance of Barnard's
#   star (1.828 pc, from its Gaia DR3 parallax): the bound holds for any
#   orbit
# (val: the acceleration of the last fit, the mixture [m/s/yr]; then the
#   distance [pc] and the separation [arcsec]; the mass in solar masses)
mmin = companion_min_mass(val, 1.828, 2.0)
# GM_SUN / GM_JUP turns solar masses into Jupiter masses
log(f'a companion 2 arcsec away has at least {mmin:.5f} Msun, '
    f'{mmin * GM_SUN / GM_JUP:.2f} MJup', 'value')

# 2. the perspective acceleration: Barnard's star, the classic case, from
#   its Gaia DR3 proper motion and parallax (mas/yr, mas/yr, mas)
# (with no covariance given, the error it returns is zero)
acc, _ = perspective_acceleration(pmra=-801.551, pmdec=10362.394,
                                  parallax=546.976)
log(f"Barnard's star: {acc:.3f} m/s per year of perspective", 'value')
# the same by name: SIMBAD finds the Gaia DR3 source, the Gaia archive
#   gives its astrometry with errors and correlations (this needs the
#   network)
astro = gaia_astrometry('GJ 699')
# mu^2 d and its error, from the covariance of the astrometry [m/s/yr]
acc, err = perspective_from_astrometry(astro)
log(f'GJ 699 is Gaia DR3 {astro["source_id"]}: perspective acceleration '
    f'{acc:.4f} +- {err:.1e} m/s/yr', 'value')
# velocities that still hold it, as a barycentric correction that ignores
#   the motion of the star leaves them (APERO's do not): the series of the
#   first part, with its own 2 m/s/yr, plus the drift of Barnard's
#   perspective
held = data.with_values(data.rv + perspective_drift(data.time, acc,
                                                    data.tref))
# perspective=(value, error): a parameter with that gaussian prior
#   [m/s/yr], so that the trend takes only the star's own acceleration
fit = RVModel(held, [], trend=1, unit='both',
              perspective=(acc, err)).fit(nstart=2, quiet=True)
# perspective, accel and total (their sum), each (value, minus, plus)
#   [m/s/yr]
res = acceleration(fit)
log(f'fit with the Gaia prior: perspective {res["perspective"][0]:.3f}, '
    f'the star itself {res["accel"][0]:+.2f} +- {res["accel"][1]:.2f} '
    f'(true {accel}), together {res["total"][0]:+.2f} m/s/yr', 'value')

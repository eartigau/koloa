#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Keplerian orbits.

The Kepler solver, the velocity of an orbit in radvel's convention (omega is
the argument of periastron of the star), the conversions between the time
of periastron and of conjunction, and m sin i from the mass function
without assuming m << M.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's Kepler tools and log (timestamped lines, 'value' for
#   a number)
import numpy as np

from koloa import kepler
from koloa.log import log

# Kepler's equation E - e sin E = M, solved to 1e-12 even at e = 0.95
# M, the mean anomaly [rad]: one turn in 1000 steps
mean_anomaly = np.linspace(0, 2 * np.pi, 1001)
for ecc in (0.1, 0.5, 0.95):
    # E, the eccentric anomaly [rad], and how far E - e sin E misses M
    eanom = kepler.eccentric_anomaly(mean_anomaly, ecc)
    worst = np.max(np.abs(eanom - ecc * np.sin(eanom) - mean_anomaly))
    log(f'e = {ecc:4.2f}: largest residual of the solution {worst:.1e} rad',
        'value')

# an eccentric orbit: the velocity is asymmetric about zero
# P [days], e, omega [rad] (the argument of periastron of the star, as in
#   radvel) and K [m/s]
period, ecc, omega, amp = 11.2, 0.3, np.radians(60), 6.0
# the time of inferior conjunction (a transit, if one is seen) [days], and
#   the time of periastron that goes with it
tconj = 60000.0
tperi = kepler.tc_to_tp(tconj, period, ecc, omega)
# one period of the velocity, v = K [cos(nu + omega) + e cos(omega)], nu
#   the true anomaly
time = np.linspace(tconj, tconj + period, 2001)
rv = kepler.rv_keplerian(time, period, tperi, ecc, omega, amp)
log(f'P = {period} d, e = {ecc}, K = {amp} m/s: velocity from '
    f'{rv.min():.2f} to {rv.max():.2f} m/s', 'value')
# tp_to_tc goes back: the round trip should land on tconj
log(f'periastron {tconj - tperi:.3f} d before the conjunction; back: '
    f'{kepler.tp_to_tc(tperi, period, ecc, omega) - tconj:+.1e} d', 'value')

# koloa's sampling basis and back
# (log P, sqrt(K) cos l, sqrt(K) sin l, sqrt(e) cos w, sqrt(e) sin w), with
#   l the mean longitude at tref [days] and w = omega; basis_to_orbit
#   returns (P, tp, e, omega, K)
basis = kepler.orbit_to_basis(period, tperi, ecc, omega, amp, tref=60100.0)
back = kepler.basis_to_orbit(*basis, tref=60100.0)
log(f'basis (log P, sqrt(K) cos l, sqrt(K) sin l, sqrt(e) cos w, '
    f'sqrt(e) sin w) = ({", ".join(f"{val:.3f}" for val in basis)})')
log(f'and back: P = {back[0]:.4f} d, e = {back[2]:.4f}, K = {back[4]:.4f} '
    f'm/s', 'value')

# m sin i: the Earth, Jupiter, and 51 Peg b (Rosenthal et al. 2021, who
#   give 0.464 Jupiter masses)
# each: K [m/s], P [days], e, and the mass of the star [solar masses]
for name, kval, per, eval_, mstar in (
        ('Earth', 0.0894, 365.256, 0.0167, 1.0),
        ('Jupiter', 12.47, 4332.59, 0.0489, 1.0),
        ('51 Peg b', 55.73, 4.2307969, 0.0042, 1.069)):
    # the mass function solved without assuming m << M [Earth masses]
    msini = kepler.minimum_mass(kval, per, eval_, mstar)
    # the ratio of the masses turns Earth masses into Jupiter masses
    log(f'{name:8s} K = {kval:>6g} m/s, P = {per:>9g} d: m sin i = '
        f'{msini:6.1f} Earth masses, '
        f'{msini * kepler.MEARTH_MSUN / kepler.MJUP_MSUN:.3f} Jupiter masses',
        'value')

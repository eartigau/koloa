#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Keplerian orbits: the solver, the velocity, and the conversions around them.

The conventions are radvel's (Fulton et al. 2018), so that an orbit fitted
here reads the same there: omega is the argument of periastron of the
star's orbit, and

    v(t) = K [cos(nu + omega) + e cos(omega)].

The sampling basis of koloa is (log P, sqrt(K) cos(lambda),
sqrt(K) sin(lambda), sqrt(e) cos(omega), sqrt(e) sin(omega)), with lambda
the mean longitude at the reference epoch: uniform priors on the last four
are uniform priors on K, lambda, e and omega, and neither K = 0 nor e = 0
is a boundary the sampler has to sit against.

Created on 2026-09-27

@author: artigau
"""
from typing import Tuple

import numpy as np

# =============================================================================
# Define variables
# =============================================================================
#: the tolerance of the Kepler solver [rad]
KEPLER_TOL = 1.0e-12

#: the gravitational constant times a solar mass [m^3 s^-2]
GM_SUN = 1.32712440018e20
#: the mass of Jupiter and of the Earth over the Sun's: the ratios of the
#: IAU 2015 nominal GMs (Jupiter the planet, not its system with its
#: satellites, 1/1047.35)
MJUP_MSUN = 9.545942e-4
MEARTH_MSUN = 3.003489e-6
DAY = 86400.0


# =============================================================================
# Define functions
# =============================================================================
def eccentric_anomaly(mean_anomaly: np.ndarray, ecc) -> np.ndarray:
    """
    Solve Kepler's equation E - e sin E = M, vectorised

    Newton's method from Danby's starting point E0 = M + 0.85 e sign(sin M),
    which converges for every eccentricity below one in a handful of steps.

    :param mean_anomaly: np.ndarray, the mean anomaly [rad]
    :param ecc: float or np.ndarray, the eccentricity (0 <= e < 1)

    :return: np.ndarray, the eccentric anomaly [rad]
    """
    mean_anomaly = np.asarray(mean_anomaly, dtype=float)
    ecc = np.broadcast_to(np.asarray(ecc, dtype=float), mean_anomaly.shape)
    if np.all(ecc == 0):
        return mean_anomaly.copy()
    manom = np.mod(mean_anomaly, 2 * np.pi)
    eanom = manom + 0.85 * ecc * np.sign(np.sin(manom))
    for _ in range(100):
        sin_e, cos_e = np.sin(eanom), np.cos(eanom)
        func = eanom - ecc * sin_e - manom
        # Halley's correction on top of Newton's for the high eccentricities
        dfunc = 1 - ecc * cos_e
        d2func = ecc * sin_e
        step = func / dfunc
        step = func / (dfunc - 0.5 * step * d2func)
        eanom = eanom - step
        if np.max(np.abs(step)) < KEPLER_TOL:
            break
    # back to the same turn as the input
    return eanom + (mean_anomaly - manom)


def true_anomaly(time: np.ndarray, period: float, tperi: float, ecc: float
                 ) -> np.ndarray:
    """
    The true anomaly at each time

    :param time: np.ndarray, the time [days]
    :param period: float, the period [days]
    :param tperi: float, the time of periastron [days]
    :param ecc: float, the eccentricity

    :return: np.ndarray, the true anomaly [rad]
    """
    manom = 2 * np.pi * (np.asarray(time, dtype=float) - tperi) / period
    eanom = eccentric_anomaly(manom, ecc)
    return 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(eanom / 2),
                          np.sqrt(1 - ecc) * np.cos(eanom / 2))


def rv_keplerian(time: np.ndarray, period: float, tperi: float, ecc: float,
                 omega: float, amp: float) -> np.ndarray:
    """
    The radial velocity of one orbit, radvel's convention

    :param time: np.ndarray, the time [days]
    :param period: float, the period [days]
    :param tperi: float, the time of periastron [days]
    :param ecc: float, the eccentricity
    :param omega: float, the argument of periastron of the star [rad]
    :param amp: float, the semi-amplitude K [m/s]

    :return: np.ndarray, the velocity [m/s]
    """
    if ecc == 0:
        # a circular orbit: nu is the mean anomaly
        phase = 2 * np.pi * (np.asarray(time, dtype=float) - tperi) / period
        return amp * np.cos(phase + omega)
    nu = true_anomaly(time, period, tperi, ecc)
    return amp * (np.cos(nu + omega) + ecc * np.cos(omega))


def tc_to_tp(tconj: float, period: float, ecc: float, omega: float) -> float:
    """
    Time of periastron from the time of inferior conjunction

    (radvel.orbit.timetrans_to_timeperi)

    :param tconj: float, the time of conjunction (transit) [days]
    :param period: float, the period [days]
    :param ecc: float, the eccentricity
    :param omega: float, the argument of periastron [rad]

    :return: float, the time of periastron [days]
    """
    fanom = np.pi / 2 - omega
    eanom = 2 * np.arctan(np.tan(fanom / 2) * np.sqrt((1 - ecc) / (1 + ecc)))
    return float(tconj - period / (2 * np.pi) * (eanom - ecc * np.sin(eanom)))


def tp_to_tc(tperi: float, period: float, ecc: float, omega: float) -> float:
    """
    Time of inferior conjunction from the time of periastron

    :param tperi: float, the time of periastron [days]
    :param period: float, the period [days]
    :param ecc: float, the eccentricity
    :param omega: float, the argument of periastron [rad]

    :return: float, the time of conjunction [days]
    """
    fanom = np.pi / 2 - omega
    eanom = 2 * np.arctan(np.tan(fanom / 2) * np.sqrt((1 - ecc) / (1 + ecc)))
    return float(tperi + period / (2 * np.pi) * (eanom - ecc * np.sin(eanom)))


def basis_to_orbit(logp: float, xk: float, yk: float, xe: float, ye: float,
                   tref: float) -> Tuple[float, float, float, float, float]:
    """
    koloa's sampling basis to (P, tp, e, omega, K)

    :param logp: float, log of the period [log days]
    :param xk: float, sqrt(K) cos(lambda)
    :param yk: float, sqrt(K) sin(lambda)
    :param xe: float, sqrt(e) cos(omega)
    :param ye: float, sqrt(e) sin(omega)
    :param tref: float, the reference epoch of lambda [days]

    :return: tuple, period, time of periastron, e, omega, K
    """
    period = float(np.exp(logp))
    amp = xk ** 2 + yk ** 2
    lam = float(np.arctan2(yk, xk))
    ecc = xe ** 2 + ye ** 2
    omega = float(np.arctan2(ye, xe)) if ecc > 0 else 0.0
    # the mean anomaly at the reference epoch is lambda - omega
    manom = lam - omega
    tperi = tref - manom * period / (2 * np.pi)
    return period, float(tperi), float(ecc), omega, float(amp)


def orbit_to_basis(period: float, tperi: float, ecc: float, omega: float,
                   amp: float, tref: float) -> Tuple[float, ...]:
    """
    (P, tp, e, omega, K) to koloa's sampling basis

    :return: tuple, log P, sqrt(K) cos(lambda), sqrt(K) sin(lambda),
             sqrt(e) cos(omega), sqrt(e) sin(omega)
    """
    manom = 2 * np.pi * (tref - tperi) / period
    lam = manom + omega
    sqk, sqe = np.sqrt(max(amp, 0.0)), np.sqrt(max(ecc, 0.0))
    return (float(np.log(period)), float(sqk * np.cos(lam)),
            float(sqk * np.sin(lam)), float(sqe * np.cos(omega)),
            float(sqe * np.sin(omega)))


def minimum_mass(amp: float, period: float, ecc: float, mstar: float
                 ) -> float:
    """
    m sin i of the companion, from the mass function, exactly

    The mass function f = (m sin i)^3 / (M + m)^2 = P K^3 (1 - e^2)^1.5 /
    (2 pi G) is solved for m sin i without assuming m << M (a fixed point
    iteration, a few steps).

    :param amp: float, K [m/s]
    :param period: float, P [days]
    :param ecc: float, e
    :param mstar: float, the stellar mass [solar masses]

    :return: float, m sin i [Earth masses]
    """
    massfn = (period * DAY * amp ** 3 * (1 - ecc ** 2) ** 1.5
              / (2 * np.pi * GM_SUN))
    msini = (massfn * mstar ** 2) ** (1.0 / 3.0)
    for _ in range(50):
        msini = (massfn * (mstar + msini) ** 2) ** (1.0 / 3.0)
    return float(msini / MEARTH_MSUN)


# =============================================================================
# End of code
# =============================================================================

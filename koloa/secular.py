#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The acceleration of a star, with its error.

"Secular acceleration" is used in two senses, and this module has both.

1. The acceleration of the star itself, dv_r/dt, from a companion too far
   out to complete an orbit during the series (or an instrument). It is
   measured, not known: the trend of RVModel is a polynomial in time, and
   acceleration() reads it in physical units, dv_r/dt in m/s/yr at the
   reference time and, with a trend of degree two, its change in
   m/s/yr^2, with their errors from the fit (outlier-aware with koloa's
   likelihood). A companion at a projected separation rho that makes an
   acceleration a has a mass of at least

       M >= (sqrt(27) / 2) |a| rho^2 / G

   whatever its orbit (Torres 1999; Liu et al. 2002: 5.34e-6 Msun times
   (d/pc rho/arcsec)^2 |a|/(m/s/yr) times a factor of at least sqrt(27)/2,
   reached when the companion lies rho / sqrt(2) behind or in front of the
   star); companion_min_mass() gives it, for draws of a too.

2. The perspective acceleration: a star that crosses the sky at a transverse
   velocity v_t = mu d sees its radial velocity grow even at constant space
   velocity, because the line of sight turns towards its motion:

       dv_r / dt = v_t^2 / d = mu^2 d

   (Kurster et al. 2003 for Barnard's star, 4.5 m/s/yr; also called the
   secular acceleration). It is geometry, known in advance from the proper
   motion and the parallax, so it enters a fit as a parameter with a
   gaussian prior from the astrometry (RVModel(perspective=(value,
   error))) rather than as a free trend: a trend fitted at the same time
   takes the acceleration of the star itself.

   With mu in mas/yr and the parallax varpi in mas,

       a [m/s/yr] = C mu^2 / varpi,   C = mas^2 pc 1000 / yr = 2.2983e-5,

   and its error follows from the covariance of (pmra, pmdec, parallax),
   correlations included (Gaia gives them). The first-order propagation is
   exact enough whenever the parallax is measured to better than a few per
   cent; perspective_acceleration also draws the astrometry to check it.
   The second-order terms (the change of the acceleration itself, of order
   3 a v_r / d) are below a cm/s over a decade for every star but a few,
   and are neglected.

   Check first that the velocities still hold it. A barycentric correction
   that is given the proper motion and the parallax of the star removes it
   already: barycorrpy does (Kanodia & Wright 2018, Eq. 28 of Wright &
   Eastman 2014), and APERO, which reduces SPIRou and NIRPS, passes it the
   Gaia astrometry of every target, so APERO velocities carry no
   perspective acceleration and adding this term to their fit would count
   it twice (the SPIRou series of Barnard's star drifts by 1.2 m/s/yr, not
   4.5). The term is for velocities whose barycentric correction ignores
   the motion of the star.

   gaia_astrometry fetches the astrometry of a star by name: SIMBAD gives
   its Gaia DR3 identifier, the Gaia archive its astrometry.

The first version of this module called the perspective functions
secular_acceleration, secular_from_astrometry, secular_drift and
remove_secular; those names still work.

Created on 2026-09-27

@author: artigau
"""
import csv
import io
import urllib.parse
import urllib.request
from typing import Any, Dict, Optional, Tuple

import numpy as np

from koloa.data import RVData

# =============================================================================
# Define variables
# =============================================================================
#: one milliarcsecond [rad]
MAS = np.pi / (180.0 * 3600.0 * 1000.0)
#: one parsec [m] (IAU 2015 B2)
PARSEC = 3.0856775814913673e16
#: one Julian year [s], the unit of the Gaia proper motions
JULIAN_YEAR = 365.25 * 86400.0
#: a [m/s/yr] = SECULAR_CONSTANT mu^2 / parallax, mu in mas/yr, parallax
#: in mas (the perspective acceleration)
SECULAR_CONSTANT = MAS ** 2 * PARSEC * 1000.0 / JULIAN_YEAR
#: one astronomical unit [m] (IAU 2012 B2)
AU = 1.495978707e11
#: the nominal solar and jovian mass parameters G M [m^3/s^2] (IAU 2015 B3)
GM_SUN = 1.3271244e20
GM_JUP = 1.2668653e17
#: the mass [Msun] that pulls a star by 1 m/s/yr from 1 au along the line of
#: sight: au^2 / (yr G Msun) = 5.34e-6 (Torres 1999)
ACCEL_MASS = AU ** 2 / (JULIAN_YEAR * GM_SUN)
#: the smallest value of the geometric factor of Torres (1999), reached for
#: a companion rho / sqrt(2) in front of or behind the star
ACCEL_FACTOR_MIN = np.sqrt(27.0) / 2.0
#: the services
SIMBAD_TAP = 'https://simbad.cds.unistra.fr/simbad/sim-tap/sync'
GAIA_TAP = 'https://gea.esac.esa.int/tap-server/tap/sync'
GAIA_COLUMNS = ['source_id', 'ra', 'dec', 'parallax', 'parallax_error',
                'pmra', 'pmra_error', 'pmdec', 'pmdec_error',
                'parallax_pmra_corr', 'parallax_pmdec_corr',
                'pmra_pmdec_corr', 'radial_velocity',
                'radial_velocity_error', 'ruwe']


# =============================================================================
# Define functions
# =============================================================================
def acceleration(res) -> Dict[str, Any]:
    """
    The acceleration of the star that a fit measured, with its errors

    The trend of RVModel is a polynomial in (t - tref) / T, T the baseline;
    at the reference time its first two derivatives are

        dv/dt = c1 / T,    d2v/dt2 = 2 c2 / T^2

    (T in years): the acceleration of the star [m/s/yr] and its change
    [m/s/yr^2]. A perspective acceleration fitted beside it
    (RVModel(perspective=...)) is given apart, and in the sum that the
    velocities measure (total).

    :param res: FitResult, sampled (errors from the chain) or a maximum a
                posteriori (errors from its Laplace covariance)

    :return: dict: accel (dv/dt) and, when fitted, jerk (d2v/dt2),
             perspective and total, each (value, minus, plus); draws, the
             acceleration of every draw of the chain (None for a maximum);
             tref, the reference time [days]
    :raises ValueError: when the fit has no trend
    """
    model = res.model
    if model.trend < 1:
        raise ValueError('The fit has no trend: fit with trend=1 or more')
    years = max(model.data.baseline, 1e-9) / 365.25
    npar = len(model.names)
    # the linear maps from the parameters to what is reported
    maps = {}
    maps['accel'] = np.zeros(npar)
    maps['accel'][model.index['trend_1']] = 1.0 / years
    if model.trend >= 2:
        maps['jerk'] = np.zeros(npar)
        maps['jerk'][model.index['trend_2']] = 2.0 / years ** 2
    if 'secacc' in model.index:
        maps['perspective'] = np.zeros(npar)
        maps['perspective'][model.index['secacc']] = 1.0
        maps['total'] = maps['accel'] + maps['perspective']
    out: Dict[str, Any] = dict(tref=float(model.tref), draws=None)
    sampled = res.chain is not None and len(res.chain) > 0
    for key, vec in maps.items():
        if sampled:
            vals = res.chain @ vec
            low, mid, high = np.percentile(vals, [16, 50, 84])
            out[key] = (float(mid), float(mid - low), float(high - mid))
            if key == 'accel':
                out['draws'] = vals
        else:
            sig = (float(np.sqrt(vec @ res.cov @ vec))
                   if res.cov is not None else float('nan'))
            out[key] = (float(res.theta @ vec), sig, sig)
    return out


def companion_min_mass(accel, distance: float, separation: float):
    """
    The smallest mass of a companion that makes a star accelerate by accel
    along the line of sight from a projected separation (Torres 1999)

        M_min = (sqrt(27) / 2) |accel| (d rho)^2 au^2 / (yr G)

    whatever its orbit; the bound is reached when the companion lies
    rho / sqrt(2) in front of or behind the star.

    :param accel: float or np.ndarray, the acceleration [m/s/yr] (draws of
                  it give draws of the bound)
    :param distance: float, the distance of the star [pc]
    :param separation: float, the angular separation of the companion
                       [arcsec]

    :return: float or np.ndarray, the minimum mass [Msun] (times
             GM_SUN / GM_JUP in jovian masses)
    """
    rho = float(distance) * float(separation)
    return ACCEL_FACTOR_MIN * ACCEL_MASS * rho ** 2 * np.abs(accel)


def companion_acceleration(mass: float, distance: float, separation: float,
                           depth: Optional[float] = None) -> float:
    """
    The acceleration along the line of sight that a companion gives a star

        a = G M z / (rho^2 + z^2)^(3/2)

    :param mass: float, the mass of the companion [Msun]
    :param distance: float, the distance of the star [pc]
    :param separation: float, the angular separation [arcsec]
    :param depth: float or None, how far the companion lies behind the star
                  along the line of sight [au] (negative in front); the
                  largest acceleration, at rho / sqrt(2), when None

    :return: float, the acceleration of the star towards the companion
             along the line of sight [m/s/yr]
    """
    rho = float(distance) * float(separation)
    depth = rho / np.sqrt(2.0) if depth is None else float(depth)
    return float(mass * depth / (rho ** 2 + depth ** 2) ** 1.5
                 / ACCEL_MASS)


def astrometry_covariance(astro: Dict[str, float]) -> np.ndarray:
    """
    The covariance of (pmra, pmdec, parallax) from their errors and
    correlations (Gaia's names; a missing correlation counts as zero)

    :param astro: dict, pmra_error, pmdec_error, parallax_error and
                  optionally pmra_pmdec_corr, parallax_pmra_corr,
                  parallax_pmdec_corr

    :return: np.ndarray, (3 x 3), in (mas/yr, mas/yr, mas) squared
    """
    sig = np.array([astro['pmra_error'], astro['pmdec_error'],
                    astro['parallax_error']], dtype=float)
    corr = np.eye(3)
    for (ii, jj), key in (((0, 1), 'pmra_pmdec_corr'),
                          ((0, 2), 'parallax_pmra_corr'),
                          ((1, 2), 'parallax_pmdec_corr')):
        val = astro.get(key)
        if val is not None and np.isfinite(float(val)):
            corr[ii, jj] = corr[jj, ii] = float(val)
    return corr * np.outer(sig, sig)


def perspective_acceleration(pmra: float, pmdec: float, parallax: float,
                             cov: Optional[np.ndarray] = None,
                             ndraw: int = 0, seed: int = 1
                             ) -> Tuple[float, float]:
    """
    The perspective acceleration, mu^2 d, and its one-sigma error

    :param pmra: float, the proper motion in right ascension, times cos(dec)
                 [mas/yr]
    :param pmdec: float, the proper motion in declination [mas/yr]
    :param parallax: float, the parallax [mas]
    :param cov: np.ndarray or None, the covariance of (pmra, pmdec,
                parallax) (see astrometry_covariance); no error when None
    :param ndraw: int, when positive, the error from this many draws of the
                  astrometry instead of the first-order propagation
    :param seed: int, the seed of the draws

    :return: tuple, the acceleration and its error [m/s/yr]
    """
    if parallax <= 0:
        raise ValueError(f'The parallax must be positive: {parallax} mas')
    mu2 = pmra ** 2 + pmdec ** 2
    value = SECULAR_CONSTANT * mu2 / parallax
    if cov is None:
        return float(value), 0.0
    cov = np.asarray(cov, dtype=float)
    if ndraw > 0:
        rng = np.random.default_rng(seed)
        draws = rng.multivariate_normal([pmra, pmdec, parallax], cov, ndraw)
        draws = draws[draws[:, 2] > 0]
        acc = SECULAR_CONSTANT * (draws[:, 0] ** 2 + draws[:, 1] ** 2) / \
            draws[:, 2]
        return float(value), float(np.std(acc))
    grad = SECULAR_CONSTANT / parallax * np.array(
        [2 * pmra, 2 * pmdec, -mu2 / parallax])
    return float(value), float(np.sqrt(grad @ cov @ grad))


def perspective_from_astrometry(astro: Dict[str, float], ndraw: int = 0
                                ) -> Tuple[float, float]:
    """
    The perspective acceleration of a star from its astrometry (Gaia's
    names)

    :param astro: dict, pmra, pmdec, parallax, their errors and
                  correlations (as gaia_astrometry returns them)
    :param ndraw: int, see perspective_acceleration

    :return: tuple, the acceleration and its error [m/s/yr]
    """
    return perspective_acceleration(
        float(astro['pmra']), float(astro['pmdec']), float(astro['parallax']),
        astrometry_covariance(astro), ndraw=ndraw)


def perspective_drift(time: np.ndarray, acceleration: float,
                      tref: float) -> np.ndarray:
    """
    The velocity a constant acceleration adds, relative to tref (for the
    perspective acceleration, or any other)

    :param time: np.ndarray, the times [days]
    :param acceleration: float, the acceleration [m/s/yr]
    :param tref: float, the time at which the drift is zero [days]

    :return: np.ndarray, [m/s]
    """
    return acceleration * (np.asarray(time, dtype=float) - tref) / 365.25


def remove_perspective(data: RVData, acceleration: float,
                       tref: Optional[float] = None) -> RVData:
    """
    A copy of a series with the drift of a known acceleration taken out
    (the perspective acceleration, for velocities that still hold it)

    The error of the acceleration is not added to the error bars: it moves
    every point together along a slope, which no diagonal error can say.
    Fit it instead (RVModel(perspective=(value, error))).

    :param data: RVData, the series
    :param acceleration: float, the acceleration [m/s/yr]
    :param tref: float or None, where the drift is zero (the reference time
                 of the series when None)

    :return: RVData, the corrected copy
    """
    tref = data.tref if tref is None else tref
    return data.with_values(data.rv - perspective_drift(data.time,
                                                        acceleration, tref))


def _tap(url: str, query: str, timeout: float) -> list:
    """One synchronous ADQL query, as a list of dict"""
    body = urllib.parse.urlencode(dict(REQUEST='doQuery', LANG='ADQL',
                                       FORMAT='csv', QUERY=query)).encode()
    request = urllib.request.Request(url, data=body)
    with urllib.request.urlopen(request, timeout=timeout) as resp:
        text = resp.read().decode()
    return list(csv.DictReader(io.StringIO(text)))


def gaia_astrometry(name: str, timeout: float = 60.0) -> Dict[str, Any]:
    """
    The Gaia DR3 astrometry of a star, by name

    SIMBAD resolves the name to a Gaia DR3 source; the Gaia archive gives
    its parallax, proper motions, their errors and correlations, its radial
    velocity and its RUWE. Needs the network.

    :param name: str, any name SIMBAD knows (GJ 699, Kepler-21, TOI-2120)
    :param timeout: float, per query [s]

    :return: dict, the Gaia columns (floats) plus name and source_id
    """
    safe = name.replace("'", "''")
    rows = _tap(SIMBAD_TAP,
                "SELECT id2.id FROM ident AS id1 JOIN ident AS id2 "
                f"USING(oidref) WHERE id1.id = '{safe}' "
                "AND id2.id LIKE 'Gaia DR3 %'", timeout)
    if not rows:
        raise ValueError(f'SIMBAD knows no Gaia DR3 source for {name}')
    source = int(rows[0]['id'].split()[-1])
    rows = _tap(GAIA_TAP, f"SELECT {', '.join(GAIA_COLUMNS)} FROM "
                          f"gaiadr3.gaia_source WHERE source_id = {source}",
                timeout)
    if not rows:
        raise ValueError(f'The Gaia archive has no source {source}')
    out: Dict[str, Any] = dict(name=name, source_id=source)
    for key in GAIA_COLUMNS[1:]:
        val = rows[0].get(key, '')
        out[key] = float(val) if val not in ('', None) else float('nan')
    return out


#: the names of the first version, when "secular" meant the perspective term
secular_acceleration = perspective_acceleration
secular_from_astrometry = perspective_from_astrometry
secular_drift = perspective_drift
remove_secular = remove_perspective


# =============================================================================
# End of code
# =============================================================================

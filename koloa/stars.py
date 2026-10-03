#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The star: its spectral type and a rough mass, for the minimum mass of a
companion

The mass of a main-sequence star from its spectral type is that of the
table of Pecaut & Mamajek (2013, ApJS 208, 9), as Mamajek keeps it
(EEM_dwarf_UBVIJHK_colors_Teff.txt, version 2022.04.16), interpolated in
the subclass; it is rough (a dwarf of a type spans some 10 % in mass), and
the NASA Exoplanet Archive's mass of a planet host is used first when
there is one. The masses of the Earth and Jupiter are the IAU 2015 nominal
values; Neptune's, its GM from JPL Horizons (6835099.97 km^3/s^2, the planet
without Triton).

Created on 2026-10-03

@author: artigau
"""
import re
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from koloa import kepler

# =============================================================================
# Define variables
# =============================================================================
#: the mass of a dwarf by spectral type [solar masses]: Pecaut & Mamajek
#: (2013), Mamajek's table of 2022.04.16
SPT_MASS: List[Tuple[str, float]] = [
    ('O3V', 59.0),
    ('O4V', 48.0),
    ('O5V', 43.0),
    ('O5.5V', 38.0),
    ('O6V', 35.0),
    ('O6.5V', 31.0),
    ('O7V', 28.0),
    ('O7.5V', 26.0),
    ('O8V', 23.6),
    ('O8.5V', 21.9),
    ('O9V', 20.2),
    ('O9.5V', 18.7),
    ('B0V', 17.7),
    ('B0.5V', 14.8),
    ('B1V', 11.8),
    ('B1.5V', 9.9),
    ('B2V', 7.3),
    ('B2.5V', 6.1),
    ('B3V', 5.4),
    ('B4V', 5.1),
    ('B5V', 4.7),
    ('B6V', 4.3),
    ('B7V', 3.92),
    ('B8V', 3.38),
    ('B9V', 2.75),
    ('B9.5V', 2.68),
    ('A0V', 2.18),
    ('A1V', 2.05),
    ('A2V', 1.98),
    ('A3V', 1.86),
    ('A4V', 1.93),
    ('A5V', 1.88),
    ('A6V', 1.83),
    ('A7V', 1.77),
    ('A8V', 1.81),
    ('A9V', 1.75),
    ('F0V', 1.61),
    ('F1V', 1.5),
    ('F2V', 1.46),
    ('F3V', 1.44),
    ('F4V', 1.38),
    ('F5V', 1.33),
    ('F6V', 1.25),
    ('F7V', 1.21),
    ('F8V', 1.18),
    ('F9V', 1.13),
    ('F9.5V', 1.08),
    ('G0V', 1.06),
    ('G1V', 1.03),
    ('G2V', 1.0),
    ('G3V', 0.99),
    ('G4V', 0.985),
    ('G5V', 0.98),
    ('G6V', 0.97),
    ('G7V', 0.95),
    ('G8V', 0.94),
    ('G9V', 0.9),
    ('K0V', 0.88),
    ('K1V', 0.86),
    ('K2V', 0.82),
    ('K3V', 0.78),
    ('K4V', 0.73),
    ('K5V', 0.7),
    ('K6V', 0.69),
    ('K7V', 0.64),
    ('K8V', 0.62),
    ('K9V', 0.59),
    ('M0V', 0.57),
    ('M0.5V', 0.54),
    ('M1V', 0.5),
    ('M1.5V', 0.47),
    ('M2V', 0.44),
    ('M2.5V', 0.4),
    ('M3V', 0.37),
    ('M3.5V', 0.27),
    ('M4V', 0.23),
    ('M4.5V', 0.184),
    ('M5V', 0.162),
    ('M5.5V', 0.123),
    ('M6V', 0.102),
    ('M6.5V', 0.093),
    ('M7V', 0.09),
    ('M7.5V', 0.088),
    ('M8V', 0.085),
    ('M8.5V', 0.08),
    ('M9V', 0.079),
    ('M9.5V', 0.078),
    ('L0V', 0.077),
    ('L1V', 0.076),
    ('L2V', 0.075)]
#: the rough relative error of a mass from a spectral type
SPT_MASS_ERR = 0.10
#: the masses of Neptune and Jupiter in Earth masses (GM ratios: Neptune's
#: from JPL Horizons, Jupiter's and the Earth's the IAU 2015 nominal values)
MNEP_MEARTH = 6835099.97e9 / 3.986004e14
MJUP_MEARTH = kepler.MJUP_MSUN / kepler.MEARTH_MSUN
#: the letters of the spectral classes, in order
CLASSES = 'OBAFGKML'


# =============================================================================
# Define functions
# =============================================================================
def spectral_code(sptype: Any) -> Optional[float]:
    """
    A spectral type as a number (O0 = 0, B0 = 10, ... M0 = 60, L0 = 70; the
    subclass after), for a dwarf or a type with no luminosity class; None
    for a giant or a subgiant (I to IV), a white dwarf or what cannot be
    read ('M3.5Ve' 63.5, 'dM2' 62, 'K7V+M0V' 57)

    :param sptype: str, the spectral type (SIMBAD's)

    :return: float or None
    """
    if not sptype:
        return None
    text = str(sptype).strip()
    text = re.sub(r'^(sd|d)', '', text)
    match = re.match(r'([OBAFGKML])\s*([0-9](?:\.[0-9]+)?)?\s*'
                     r'((?:I{1,3}|IV|V|VI)?)', text)
    if not match or text.startswith(('D', 'W')):
        return None
    if match.group(3) and match.group(3) not in ('V', 'VI'):
        return None
    sub = float(match.group(2)) if match.group(2) else 5.0
    return CLASSES.index(match.group(1)) * 10 + sub


def mass_from_spectral_type(sptype: Any) -> Optional[Tuple[float, float]]:
    """
    The rough mass of a dwarf from its spectral type (Pecaut & Mamajek
    2013), interpolated in the subclass

    :param sptype: str, the spectral type

    :return: tuple, the mass and its rough error [solar masses], or None
    """
    code = spectral_code(sptype)
    if code is None:
        return None
    codes = np.array([spectral_code(spt) for spt, _ in SPT_MASS])
    masses = np.array([mass for _, mass in SPT_MASS])
    if code < codes.min() - 1 or code > codes.max() + 1:
        return None
    mass = float(np.interp(code, codes, masses))
    return mass, SPT_MASS_ERR * mass


def stellar(ident: Dict[str, Any], archive_star: Optional[Dict[str, Any]]
            = None) -> Dict[str, Any]:
    """
    The spectral type and the mass of a star: the archive's mass when the
    star hosts known planets, else a rough one from SIMBAD's spectral type

    :param ident: dict, the star (koloa.archive.resolve: sptype)
    :param archive_star: dict or None, the archive's star (known_planets)

    :return: dict, sptype, mass, mass_err [solar masses] and source (None
             where there is none)
    """
    archive_star = archive_star or {}
    sptype = ident.get('sptype') or archive_star.get('spectral_type')
    out = dict(sptype=sptype, mass=None, mass_err=None, source=None)
    if archive_star.get('mass'):
        out.update(mass=float(archive_star['mass']),
                   mass_err=SPT_MASS_ERR * float(archive_star['mass']),
                   source='NASA Exoplanet Archive')
        return out
    found = mass_from_spectral_type(sptype)
    if found is not None:
        out.update(mass=found[0], mass_err=found[1],
                   source=f'its spectral type (Pecaut & Mamajek 2013)')
    return out


def minimum_mass(amp: float, period: float, ecc: float, mstar: float,
                 amp_err: float = 0.0, period_err: float = 0.0,
                 ecc_err: float = 0.0, mstar_err: float = 0.0,
                 ndraw: int = 4000, seed: int = 5) -> Dict[str, Any]:
    """
    m sin i of a companion and its error (draws of K, P, e and the mass of
    the star, gaussian; e kept in [0, 0.99)), in Earth, Neptune and Jupiter
    masses

    :return: dict, earth, neptune, jupiter: each (value, minus, plus), the
             16th, 50th and 84th percentiles
    """
    rng = np.random.default_rng(seed)
    draw = lambda val, err: val + (err or 0.0) * rng.standard_normal(ndraw)
    amps = np.abs(draw(amp, amp_err))
    pers = np.abs(draw(period, period_err))
    eccs = np.clip(np.abs(draw(ecc, ecc_err)), 0.0, 0.99)
    mstars = np.clip(draw(mstar, mstar_err), 0.01, None)
    massfn = (pers * kepler.DAY * amps ** 3 * (1 - eccs ** 2) ** 1.5
              / (2 * np.pi * kepler.GM_SUN))
    msini = (massfn * mstars ** 2) ** (1.0 / 3.0)
    for _ in range(50):
        msini = (massfn * (mstars + msini) ** 2) ** (1.0 / 3.0)
    earth = msini / kepler.MEARTH_MSUN
    best = kepler.minimum_mass(amp, period, ecc, mstar)
    low, high = np.percentile(earth, [15.87, 84.13])
    out = {}
    for key, unit in (('earth', 1.0), ('neptune', MNEP_MEARTH),
                      ('jupiter', MJUP_MEARTH)):
        out[key] = (best / unit, max(best - low, 0.0) / unit,
                    max(high - best, 0.0) / unit)
    return out


# =============================================================================
# End of code
# =============================================================================

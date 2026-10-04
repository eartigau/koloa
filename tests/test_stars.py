#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa.stars: a spectral type read, a rough mass from it, m sin i and its
error in Earth, Neptune and Jupiter masses

Created on 2026-10-03

@author: artigau
"""
import numpy as np
import pytest

from koloa import kepler, stars


def test_a_spectral_type_read():
    assert stars.spectral_code('M3.5Ve') == 63.5
    assert stars.spectral_code('dM2') == 62.0
    assert stars.spectral_code('K7V+M0V') == 57.0
    assert stars.spectral_code('K0III') is None
    assert stars.spectral_code('DA2') is None
    assert stars.spectral_code('') is None


def test_a_mass_from_a_spectral_type():
    # Pecaut & Mamajek: M0V 0.57, G2V 1.00; between subclasses, between
    mass, err = stars.mass_from_spectral_type('M0V')
    assert mass == 0.57 and abs(err - 0.057) < 1e-12
    assert abs(stars.mass_from_spectral_type('G2V')[0] - 1.0) < 1e-9
    low = stars.mass_from_spectral_type('M4V')[0]
    high = stars.mass_from_spectral_type('M3V')[0]
    assert low < stars.mass_from_spectral_type('M3.5V')[0] < high
    # the archive's mass first
    got = stars.stellar(dict(sptype='M0V'), dict(mass=0.6))
    assert got['mass'] == 0.6 and got['source'] == 'NASA Exoplanet Archive'
    assert stars.stellar(dict(sptype='K0III'))['mass'] is None


def test_the_masses_of_the_planets():
    # IAU 2015 nominal GM ratios, and Neptune's from JPL
    assert abs(stars.MJUP_MEARTH - 317.8284) < 1e-3
    assert abs(stars.MNEP_MEARTH - 17.1477) < 1e-3
    got = stars.minimum_mass(17.09, 2.64388, 0.138, 0.45, 0.22, 1e-6, 0.01,
                             0.045)
    best = kepler.minimum_mass(17.09, 2.64388, 0.138, 0.45)
    assert got['earth'][0] == best
    assert abs(got['jupiter'][0] * stars.MJUP_MEARTH - best) < 1e-9
    # its error, mostly the mass of the star's: (2/3) 10 %
    assert abs(0.5 * (got['earth'][1] + got['earth'][2]) / best
               - np.hypot(2 / 3 * 0.1, 0.22 / 17.09)) < 0.01


def test_the_equilibrium_temperature():
    """the Earth about the Sun: 278 K with no albedo, 255 K with 0.3; the
    effective temperature of a dwarf from its type (Pecaut & Mamajek)"""
    from koloa import stars
    assert stars.equilibrium_temperature(5772, 1.0, 1.0, 365.25) == \
        pytest.approx(278, abs=1)
    assert stars.equilibrium_temperature(5772, 1.0, 1.0, 365.25, 0.3) == \
        pytest.approx(255, abs=1)
    assert stars.teff_from_spectral_type('M2V') == 3560
    assert stars.teff_from_spectral_type('G2V') == 5770
    assert stars.stellar_teff(dict(sptype='K5V'), tic=False)['teff'] == 4440
    assert stars.stellar_teff(dict(sptype='K5V'), dict(teff=4500.0),
                              tic=False)['source'] == 'NASA Exoplanet Archive'

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Doppler shifts, relativistic, never first order.

koloa works in velocity and never needs a wavelength itself; these helpers
are for the steps around it (a velocity turned into a shift of a spectrum,
a pixel shift on a log-wavelength grid read back as a velocity). The
first-order (1 + v/c) is off by about (v/c)^2 / 2, 1.5 m/s at a 30 km/s
barycentric correction, which matters at the precision of radial velocity
work.

    ratio = sqrt((1 + v/c) / (1 - v/c))        (v > 0: receding)
    ln ratio = atanh(v/c)
    v = c (r^2 - 1) / (r^2 + 1) = c tanh(ln r)

Created on 2026-09-27

@author: artigau
"""
import numpy as np

# =============================================================================
# Define variables
# =============================================================================
#: the speed of light [m/s]
SPEED_OF_LIGHT = 299792458.0


# =============================================================================
# Define functions
# =============================================================================
def velocity_to_ratio(velocity):
    """
    The wavelength ratio of a velocity (observed over emitted)

    :param velocity: float or np.ndarray, the velocity [m/s], positive when
                     receding

    :return: the ratio lambda_obs / lambda_emit
    """
    beta = np.asarray(velocity, dtype=float) / SPEED_OF_LIGHT
    return np.sqrt((1 + beta) / (1 - beta))


def ratio_to_velocity(ratio):
    """
    The velocity of a wavelength ratio

    :param ratio: float or np.ndarray, lambda_obs / lambda_emit

    :return: the velocity [m/s]
    """
    ratio2 = np.asarray(ratio, dtype=float) ** 2
    return SPEED_OF_LIGHT * (ratio2 - 1) / (ratio2 + 1)


def velocity_to_logshift(velocity):
    """
    The shift in ln(wavelength) of a velocity: atanh(v / c)

    Divided by the ln step of a log-wavelength grid, this is a shift in
    pixels.

    :param velocity: float or np.ndarray, the velocity [m/s]

    :return: the shift in ln(wavelength)
    """
    return np.arctanh(np.asarray(velocity, dtype=float) / SPEED_OF_LIGHT)


def logshift_to_velocity(shift):
    """
    The velocity of a shift in ln(wavelength): c tanh(shift)

    :param shift: float or np.ndarray, the shift in ln(wavelength)

    :return: the velocity [m/s]
    """
    return SPEED_OF_LIGHT * np.tanh(np.asarray(shift, dtype=float))


def shift_wavelength(wavelength, velocity):
    """
    Wavelengths as seen from a source moving at a velocity

    :param wavelength: float or np.ndarray, the emitted wavelengths
    :param velocity: float, the velocity [m/s], positive when receding

    :return: the observed wavelengths
    """
    return np.asarray(wavelength, dtype=float) * velocity_to_ratio(velocity)


# =============================================================================
# End of code
# =============================================================================

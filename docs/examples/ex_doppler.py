#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Relativistic Doppler shifts.

koloa works in velocity. When a velocity has to become a wavelength ratio,
or a pixel shift on a log-wavelength grid has to become a velocity, these
conversions do it relativistically. The first-order (1 + v/c) is off by
1.5 m/s at a 30 km/s barycentric correction.

Created on 2026-09-27

@author: artigau
"""
# numpy, koloa's relativistic conversions (velocities in m/s), and the
#   log (timestamped lines, 'value' for a number)
import numpy as np

from koloa.doppler import (SPEED_OF_LIGHT, logshift_to_velocity,
                           ratio_to_velocity, shift_wavelength,
                           velocity_to_logshift, velocity_to_ratio)
from koloa.log import log

# 30 km/s, the size of a barycentric correction (positive: receding)
velocity = 30.0e3
# lambda_obs / lambda_emit = sqrt((1 + v/c) / (1 - v/c)), exactly
ratio = velocity_to_ratio(velocity)
# the first-order ratio, for comparison (SPEED_OF_LIGHT in m/s)
first = 1 + velocity / SPEED_OF_LIGHT
log(f'30 km/s: ratio {ratio:.12f}, first order {first:.12f}', 'value')
# ratio_to_velocity inverts the exact formula, v = c (r^2 - 1) / (r^2 + 1):
#   the velocity the first-order ratio really stands for
log(f'the first-order ratio is the shift of {ratio_to_velocity(first):.2f} '
    f'm/s, {velocity - ratio_to_velocity(first):.2f} m/s short', 'value')

# the He I line at 1083.330 nm (vacuum), from a source receding at 30 km/s
# (the wavelength is multiplied by the ratio, so it keeps its unit)
log(f'1083.330 nm is observed at '
    f'{shift_wavelength(1083.330, velocity):.4f} nm', 'value')

# a grid of constant step in ln(wavelength), 1e-6 (about 300 m/s a pixel)
step = 1.0e-6
# a velocity moves ln(wavelength) by atanh(v/c); over the step, in pixels
pixels = velocity_to_logshift(velocity) / step
log(f'30 km/s is a shift of {pixels:.5f} pixels; ln(1 + v/c) would give '
    f'{np.log(first) / step:.5f}', 'value')
# and back: a shift of 100 pixels measured on that grid
# (v = c tanh(shift), the shift in ln(wavelength))
log(f'100 pixels is {logshift_to_velocity(100 * step):.2f} m/s', 'value')

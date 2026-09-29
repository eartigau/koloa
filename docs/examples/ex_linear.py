#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The analytic evidence of linear models.

Offsets, trends and the two amplitudes of a circular orbit enter the
velocity linearly; with gaussian priors they integrate out exactly. The
gain in evidence of one more sinusoid is a two-by-two Schur complement, so
a whole grid of periods costs one matrix product. This is what makes the
FIP sampler exact and fast.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's linear models, log (timestamped lines, 'value' for a
#   number) and noise covariance
import numpy as np

from koloa.linear import LinearModel, sinusoid_gain, weighted_lstsq
from koloa.log import log
from koloa.noise import BlockCov

# 60 random times over 300 days, error bars of 1.5 m/s, and the data: an
#   offset of 2 m/s, a sinusoid of 3 m/s at 17.3 days, and the noise
rng = np.random.default_rng(3)
time = np.sort(rng.uniform(0, 300, 60))
err = np.full(len(time), 1.5)
value = (2.0 + 3.0 * np.sin(2 * np.pi * time / 17.3)
         + err * rng.normal(size=len(time)))

# an offset with a wide gaussian prior (100 m/s), integrated out
# cov: the noise covariance V, here the variances alone, on its diagonal
#   (no visit blocks) [(m/s)^2]
cov = BlockCov(err ** 2)
# y = X beta + noise, with X a column of ones; prior_var: the prior
#   variance of each column [(m/s)^2]; logz is the log evidence, ln Z, and
#   mean the posterior mean of beta
base = LinearModel(np.ones((len(time), 1)), value, cov,
                   prior_var=np.array([100.0 ** 2]))
log(f'offset only: log Z = {base.logz:.2f}, offset {base.mean[0]:.2f} m/s',
    'value')

# every period of a grid at once: Delta log Z of a sinusoid whose two
#   amplitudes have a gaussian prior of 3 m/s
# 20000 frequencies from 1/150 to 1/1.5 per day, and the cos and sin of
#   each at every time, (n x G) matrices
freq = np.linspace(1 / 150, 1 / 1.5, 20000)
cos = np.cos(2 * np.pi * np.outer(time, freq))
sin = np.sin(2 * np.pi * np.outer(time, freq))
# the products sinusoid_gain needs (V is symmetric): X^T V^-1 cos and sin
#   (wdesign is V^-1 X), y^T V^-1 cos and sin (wvalue is V^-1 y), then
#   cos^T V^-1 cos, sin^T V^-1 sin and cos^T V^-1 sin; tau2: the prior
#   variances to try [(m/s)^2], one column of the result each
gain = sinusoid_gain(base, base.wdesign.T @ cos, base.wdesign.T @ sin,
                     base.wvalue @ cos, base.wvalue @ sin,
                     cov.quad_diag(cos), cov.quad_diag(sin),
                     cov.quad_diag(cos, sin), tau2=np.array([3.0 ** 2]))[:, 0]
# the frequency whose sinusoid gains the most evidence
best = int(np.argmax(gain))
log(f'{len(freq)} periods: best {1 / freq[best]:.3f} d (injected 17.3), '
    f'Delta log Z = {gain[best]:.4f}', 'value')

# the same number the long way: a model with the two columns in it
# (with the same priors: 100 m/s on the offset, 3 m/s on each amplitude)
design = np.column_stack([np.ones(len(time)), cos[:, best], sin[:, best]])
full = LinearModel(design, value, cov, prior_var=np.array([1e4, 9.0, 9.0]))
log(f'from the two evidences: {full.logz - base.logz:.4f}', 'value')
# plain weighted least squares at that period: the coefficients, their
#   covariance (not used) and the chi^2; K is the hypotenuse of the two
#   amplitudes
coeffs, _, chi2 = weighted_lstsq(design, value, err)
log(f'least squares at that period: K = {np.hypot(*coeffs[1:]):.2f} m/s, '
    f'chi2 = {chi2:.1f} for {len(time) - 3} degrees of freedom', 'value')

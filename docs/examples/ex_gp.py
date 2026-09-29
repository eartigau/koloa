#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gaussian processes for stellar activity.

The kernels of activity work (SHO, rotation, quasi-periodic, Matern-3/2,
squared exponential) and a dense GP that numpy factorises in a millisecond
for a few hundred points: its likelihood, and its prediction between the
points.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's GP tools, log (timestamped lines, 'value' for a
#   number) and simulated campaigns
import numpy as np

from koloa.gp import KERNEL_PARAMS, DenseGP, kernel_matrix
from koloa.log import log
from koloa.simulate import activity_signal, observing_times

# a spotted star: an SHO of 4 m/s at a 23-day rotation, seen 80 times
# observing_times: 80 nights of one exposure, in observing seasons (240
#   days a year) over 400 days; activity_signal: one draw of the GP, with
#   its sigma [m/s], period [days] and quality factor Q
time = observing_times(nvisits=80, per_visit=1, baseline=400, seed=2)
truth = activity_signal(time, kernel='sho', sigma=4.0, period=23.0,
                        quality=3.0, seed=5)
# error bars of 1 m/s, and the data: the activity plus white noise
err = np.full(len(time), 1.0)
rv = truth + err * np.random.default_rng(1).normal(size=len(time))

# the parameters are logs of the natural units (the mixing term: a logit)
# sho: sigma, period, Q; rotation: sigma, period, Q0 = 1, dQ = 1 and a
#   mix of 0.5 (logit 0); qp: sigma, period, a decay of 40 days and a
#   smoothness of 0.5; matern32: sigma and a length of 10 days
pars = dict(sho=np.log([4.0, 23.0, 3.0]),
            rotation=np.array([np.log(4.0), np.log(23.0), 0.0, 0.0, 0.0]),
            qp=np.log([4.0, 23.0, 40.0, 0.5]), matern32=np.log([4.0, 10.0]))
# the lags 0, P/2 and P [days]
lags = np.array([0.0, 11.5, 23.0])
for kernel, par in pars.items():
    # the covariance of the kernel at those lags [(m/s)^2]
    cov = kernel_matrix(kernel, lags, par)
    # the GP plus the noise, factorised once: diag is the noise variance
    #   of each point [(m/s)^2]; each point is its own visit (block) and
    #   there is no visit jitter (blockval=None)
    gp = DenseGP(time, kernel, par, diag=err ** 2, block=np.arange(len(time)),
                 blockval=None)
    # KERNEL_PARAMS: the names of the parameters of each kernel, in order
    log(f'{kernel:9s} ({", ".join(KERNEL_PARAMS[kernel])})', 'info')
    # loglike: ln N(rv | 0, C), the data taken as the residual of a mean
    #   model of zero
    log(f'  k(0, P/2, P) = {", ".join(f"{val:5.1f}" for val in cov)} '
        f'(m/s)^2, ln L = {gp.loglike(rv):7.1f}', 'value')

# the SHO conditioned on the data, at the data and between them
gp = DenseGP(time, 'sho', pars['sho'], diag=err ** 2,
             block=np.arange(len(time)), blockval=None)
# predict: the mean and the 1-sigma of the GP at the times given [m/s]
mean, sigma = gp.predict(rv, time)
log(f'SHO at the data: rms of (prediction - truth) = '
    f'{np.std(mean - truth):.2f} m/s, predicted 1-sigma '
    f'{np.median(sigma):.2f} m/s', 'value')
# 2000 times across the whole campaign, the seasonal gap included
grid = np.linspace(time[0], time[-1], 2000)
_, sigma_grid = gp.predict(rv, grid)
log(f'in the seasonal gaps the 1-sigma goes back up to '
    f'{sigma_grid.max():.1f} m/s, the amplitude of the prior', 'value')

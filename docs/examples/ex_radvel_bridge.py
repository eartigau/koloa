#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa's outliers inside radvel.

A radvel user keeps the model, the priors and the sampler, and swaps one
class: RVLikelihood becomes OutlierRVLikelihood, with two more parameters,
the outlier fraction and width. Here both are maximised from the same start
on a planet (P = 12.3 d, K = 4 m/s) with 8% bad visits. radvel's model has
a white jitter but no visit jitter, which koloa's own fits add.

Created on 2026-09-27

@author: artigau
"""
# numpy, radvel, and koloa's series, log (timestamped lines, 'value' for
#   a number), radvel bridge and simulation
from pathlib import Path

import numpy as np
import radvel

from koloa.data import RVData
from koloa.log import log
from koloa.radvel_bridge import OutlierRVLikelihood, outlier_priors
from koloa.simulate import simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# a real NIRPS sampling with a circular planet (K [m/s], tp its time of
#   periastron [days]), 8 % of the visits moved as a block by 6 median
#   errors or more, and a 1.5 m/s jitter shared by each visit
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
sim = simulate(planets=[dict(P=12.3, K=4.0, tp=60102.0)], template=tpl,
               outliers=[dict(kind='visit', frac=0.08, amplitude=6.0)],
               visit_jitter=1.5, seed=31)
# the series, and which exposures were made outliers (True: bad)
data, truth = sim['data'], sim['outlier_mask']

# radvel's own likelihood, then koloa's, which takes the visit of each
#   point (seq) and makes a whole visit the unit of an outlier
for label, likelihood, extra in (
        ('radvel', radvel.likelihood.RVLikelihood, {}),
        ('radvel + koloa', OutlierRVLikelihood,
         dict(seq=data.seq, unit='sequence'))):
    # one planet in radvel's basis: period, time of conjunction,
    #   sqrt(e) cos(w), sqrt(e) sin(w) and K
    params = radvel.Parameters(1, basis='per tc secosw sesinw k')
    # the same start for both: e held at 0 and no trend (dvdt and curv
    #   held); only per1, tc1 and k1 vary
    for key, val in dict(per1=12.3, tc1=60105.0, secosw1=0.0, sesinw1=0.0,
                         k1=3.0, dvdt=0.0, curv=0.0).items():
        params[key] = radvel.Parameter(value=val, vary=key in ('per1', 'tc1',
                                                                'k1'))
    # the model and the data: time [days], velocity and error [m/s]
    like = likelihood(radvel.RVModel(params), data.time, data.rv, data.err,
                      **extra)
    # the offset gamma, fitted like any parameter (radvel solves it
    #   analytically only with linear=True and vary=False), and the white
    #   jitter, starting at 2 m/s
    like.params['gamma'] = radvel.Parameter(value=0.0, linear=False)
    like.params['jit'] = radvel.Parameter(value=2.0)
    # radvel holds the parameters as a dict and as a vector: sync them
    like.vector.dict_to_vector()
    # the likelihood and its priors: hard bounds keep the jitter and K
    #   between 0 and 50 m/s
    post = radvel.posterior.Posterior(like)
    post.priors += [radvel.prior.HardBounds('jit', 0.0, 50.0),
                    radvel.prior.HardBounds('k1', 0.0, 50.0)]
    # koloa's priors on its two parameters: a beta prior on the outlier
    #   fraction, and hard bounds on the log of the outlier width
    if likelihood is OutlierRVLikelihood:
        post.priors += outlier_priors(like)
    # the maximum a posteriori (scipy's Powell method)
    post = radvel.fitting.maxlike_fitting(post, verbose=False)
    log(f'{label:15s} P = {post.params["per1"].value:.3f} d, '
        f'K = {post.params["k1"].value:.2f} m/s, '
        f'jitter {post.params["jit"].value:.2f} m/s', 'value')
    if likelihood is OutlierRVLikelihood:
        # the probability that each point is in an outlier visit, at the
        #   best parameters
        flagged = like.outlier_probability() > 0.5
        log(f'  flagged: {np.sum(flagged & truth)} of the {truth.sum()} bad '
            f'exposures, {np.sum(flagged & ~truth)} of the {np.sum(~truth)} '
            f'good ones', 'value')

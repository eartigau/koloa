#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Clips, with their accounting.

The classical soft clip returns inflated error bars and says how much of
the weight it removed; above a budget of 5 % it logs a warning and raises a
KoloaWarning. The hard clip is there for comparison only: koloa's own
answers come from the mixture.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's series, log (timestamped lines, 'value' for a
#   number), simulation and weights
from pathlib import Path

import numpy as np

from koloa.data import RVData
from koloa.log import log
from koloa.simulate import simulate
from koloa.weights import (hard_clip, soft_clip, student_weights,
                           weight_summary)

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# a planet (P = 16.65 d, K = 5 m/s, tp its time of periastron [days]) on
#   a real NIRPS sampling, a 2 m/s jitter shared by each visit, and 2 %
#   of the exposures alone moved by 8 median errors or more: a 3-sigma
#   soft clip stays within the budget
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
planet = [dict(P=16.65, K=5.0, tp=60105.3)]
sim = simulate(planets=planet, template=tpl,
               outliers=[dict(kind='spike', frac=0.02, amplitude=8.0)],
               visit_jitter=2.0, seed=10)
# the series, and which exposures were made outliers (True: bad)
data, truth = sim['data'], sim['outlier_mask']
# soft_clip: the residual of a straight line in time, refitted over 5
#   passes; a weight is kept within clip robust sigma and falls as
#   (clip / z)^2 beyond, handed back as inflated error bars (err) with the
#   factor of each point, the fraction of the weight removed and the
#   number of points down-weighted (ndown)
soft = soft_clip(data.time, data.rv, data.err, clip=3.0)
# the exposures whose weight was cut
down = soft['factor'] < 1
log(f'a few spikes: {soft["ndown"]} exposures down-weighted '
    f'({np.sum(down & truth)} of the {truth.sum()} spikes), '
    f'{soft["removed"]:.1%} of the weight removed', 'value')

# the same planet and sampling with 7% bad visits and 3% spikes: over the
#   budget (bad visits of 6 median errors or more, spikes of 8 or more)
sim = simulate(planets=planet, template=tpl,
               outliers=[dict(kind='visit', frac=0.07, amplitude=6.0),
                         dict(kind='spike', frac=0.03, amplitude=8.0)],
               visit_jitter=2.0, seed=10)
sdata, truth = sim['data'], sim['outlier_mask']
# above 5 % of the weight removed, soft_clip logs a warning and raises a
#   KoloaWarning, which Python prints with the line of the call
soft = soft_clip(sdata.time, sdata.rv, sdata.err, clip=3.0)
# the exposures whose weight was cut, and by how much
down = soft['factor'] < 1
log(f'bad visits and spikes: soft clip down-weights {down.sum()} '
    f'exposures, {np.sum(down & truth)} of the {truth.sum()} outliers',
    'value')
log(f'  their weight factors: '
    f'{", ".join(f"{val:.2f}" for val in np.sort(soft["factor"][down]))}',
    'value')
# the hard clip removes points instead: the mask of those kept (a 3-sigma
#   clip against a line, repeated until nothing changes, 10 passes at
#   most); as factors (1 kept, 0 removed), weight_summary gives the
#   fraction of the weight removed and the number of points removed
keep = hard_clip(sdata.time, sdata.rv, sdata.err, clip=3.0)
removed, nout = weight_summary(sdata.err, keep.astype(float))
log(f'bad visits and spikes: hard clip removes {nout} exposures, '
    f'{removed:.1%} of the weight; {np.sum(~keep & truth)} of them are '
    f'outliers', 'value')

# a Student-t of 4 degrees of freedom weights a point by (4 + 1) / (4 + z^2)
# (residuals of 1 to 10 with error bars of 1: z itself)
zval = np.array([1.0, 3.0, 6.0, 10.0])
weights = student_weights(zval, np.ones(len(zval)), dof=4.0)
log('Student-t weight at ' + ', '.join(
    f'{zz:.0f} sigma: {ww:.2f}' for zz, ww in zip(zval, weights)), 'value')

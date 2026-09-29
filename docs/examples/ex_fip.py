#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The outlier-aware false inclusion probability.

oafip samples up to kmax signals, a visit jitter, and point and visit
outliers; fip_comparison gives the one-signal FIP with a fixed noise, with
no clip, a soft clip and a hard clip. Short chains here: a real analysis
uses the defaults (two chains of 1500 sweeps).

Figure: the window and the FIP periodogram of each method. The planet's
period is dashed; the dotted lines are FIP = 1% and 0.1%.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's figures, series, FIPs, log (timestamped lines,
#   'value' for a number), frequency grid and window, and simulation
from pathlib import Path

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fip import fip_comparison, oafip
from koloa.log import log
from koloa.periodogram import frequency_grid, window
from koloa.simulate import simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]
# the period of the injected planet [days]
PERIOD = 16.65

# a real NIRPS sampling lends its times, visits and error bars
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
# a circular planet, K [m/s] and tp its time of periastron [days]; 7 % of
#   the visits moved as a block by 6 median errors or more, 3 % of the
#   exposures alone by 8 or more; a 2 m/s jitter shared by each visit
sim = simulate(planets=[dict(P=PERIOD, K=5.0, tp=60105.3)], template=tpl,
               outliers=[dict(kind='visit', frac=0.07, amplitude=6.0),
                         dict(kind='spike', frac=0.03, amplitude=8.0)],
               visit_jitter=2.0, seed=10)
# the series, and which exposures were made outliers (True: bad)
data, truth = sim['data'], sim['outlier_mask']
# one grid for every method: periods from 1.1 d to twice the baseline,
#   10 points per peak width 1/T
freq = frequency_grid(data.time, pmin=1.1, oversample=10)

# the one-signal FIP with the noise held fixed: no clip, a 3-sigma soft
#   clip and a 3-sigma hard clip (the keys gaussian, soft and hard)
fips = fip_comparison(data, freq=freq)
# koloa's: up to kmax=2 signals, and outliers='both': a lone exposure or
#   a whole visit can be one; one chain of 100 burn-in and 400 recorded
#   sweeps (defaults: two of 300 and 1500); progress=False: no log lines
fips['koloa'] = oafip(data, kmax=2, outliers='both', freq=freq, nsweep=400,
                      nburn=100, nchains=1, progress=False)
# best(): the most significant peak (its period [days], fip, window...);
#   fip_at: the FIP of the interval, about 1/T wide, around a period
for res in fips.values():
    log(f'{res.method:30s} best {res.best()["period"]:7.3f} d, FIP '
        f'{res.best()["fip"]:.1e}; at {PERIOD} d: {res.fip_at(PERIOD):.1e}',
        'value')

koloa = fips['koloa']
# the probability of each number of signals, the best peaks, and the
#   expected number of outlier exposures
for line in koloa.summary().split('\n'):
    log(line)
# outlier_prob: the probability that each exposure is bad (an outlier
#   itself or in an outlier visit), averaged over the chain
flagged = koloa.outlier_prob > 0.5
log(f'flagged: {np.sum(flagged & truth)} of the {truth.sum()} bad exposures, '
    f'{np.sum(flagged & ~truth)} of the {np.sum(~truth)} good ones', 'value')

# koloa's web style; then the window of the sampling and -log10 FIP of
#   each method, with the periods of mark drawn dashed [days]
kplot.set_style('web')
fig = kplot.periodograms(freq, window_power=window(data.time, freq),
                         fips=fips, mark=[PERIOD])
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/fip.svg'))
log(f'figure: {Path(path).relative_to(ROOT)}')

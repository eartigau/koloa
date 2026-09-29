#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Which planets a series could have found.

recovery_map injects circular planets into the public HARPS-N series of
Kepler-21 over a grid of periods and semi-amplitudes, folds the series at
each period, and counts a planet as recovered when a sinusoid improves the
likelihood beyond a false-alarm probability of 1 %. Injections at K = 0
measure the false alarms: the drift and the activity of a real star break
the threshold of white noise, so a slope is fitted beside each sinusoid and
each period bin takes its threshold from its own injections at K = 0. Both
likelihoods are mapped; this series has no strong outliers, and the
outlier-aware map matches the gaussian one.

Figure: the recovered fraction on the grid, gaussian (left) and
outlier-aware (right), with the 50 % (dashed) and 90 % contours.

Created on 2026-09-27

@author: artigau
"""
# numpy for the grids, pyplot for the figure; koloa's figures, recovery
#   maps, series and log (timestamped lines, 'value' for a number)
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from koloa import plotting as kplot
from koloa.completeness import recovery_map
from koloa.data import RVData
from koloa.log import log

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# the real series, the public HARPS-N velocities of Kepler-21: each planet
#   is added to it as observed, its own outliers and signals included
data = RVData.from_csv(ROOT / 'data' / 'kepler21_harpsn_dace_drs3.3.12.csv',
                       name='Kepler-21')
# the edges of 8 bins in period and 8 in K, evenly spaced in log
period_edges = np.geomspace(1.5, 300.0, 9)    # days
amp_edges = np.geomspace(0.8, 12.0, 9)        # m/s
# test='fold': the period is known; test='search' looks for the highest
#   peak of a whole periodogram instead (slower)
# ninj: planets per cell, their period and K drawn log-uniformly in it,
#   at a random phase; one is found when a sinusoid at its period gains
#   more than the threshold in ln L; nnull: injections at K = 0, which
#   count the false alarms; outliers: the outlier-aware likelihood (True)
#   or the gaussian one (False); a visit is one unit (the default
#   unit='sequence' fits the visit means); seed: the draws; quiet: no log
#   line
# the threshold of white noise, Delta ln L > -ln(fap) = 4.6, does not hold
#   on a real star: trend=1 fits a slope beside each sinusoid (and in the
#   null) for the drift of the series, and calibrate=True makes the nnull
#   injections at K = 0 in each period bin, whose 1 - fap quantile becomes
#   the threshold of the bin when it is higher (activity raises it)
maps = {}
for label, outliers in (('gaussian', False), ('koloa', True)):
    maps[label] = recovery_map(data, period_edges, amp_edges, ninj=15,
                               test='fold', outliers=outliers, nnull=200,
                               fap=0.01, trend=1, calibrate=True, seed=1,
                               quiet=True)
# the threshold, the false alarms, and in each period bin the K recovered
#   50 % and 90 % of the time
for line in maps['koloa'].summary().split('\n'):
    log(line, 'value')
# extra holds the threshold of each period bin, and the false alarms that
#   the threshold of white noise alone would have let through
extra = maps['koloa'].extra
log(f'threshold of each period bin: '
    f'{", ".join(f"{val:.1f}" for val in extra["thresholds"])}; at 4.6 '
    f'everywhere, false alarms {extra["false_alarm_nominal"]:.1%}', 'value')
# k_at(0.5): in each period bin, the K recovered half the time (nan when
#   the bin never gets there); false_alarm: the fraction of the K = 0
#   injections counted as found
for label, rmap in maps.items():
    log(f'{label:8s} K recovered half the time: median over the periods '
        f'{np.nanmedian(rmap.k_at(0.5)):.1f} m/s; false alarms '
        f'{rmap.false_alarm:.1%}', 'value')

# koloa's web style, then two panels side by side that share the K axis
#   (figsize in inches)
kplot.set_style('web')
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.3), sharey=True)
# the recovered fraction of each cell, with its 50 % (dashed) and 90 %
#   contours; a single colour bar, beside the right panel
kplot.recovery_map(maps['gaussian'], ax=axes[0], colorbar=False,
                   title='Gaussian')
kplot.recovery_map(maps['koloa'], ax=axes[1], title='outlier-aware')
axes[1].set_ylabel('')
# savefig writes the file, closes the figure and returns the path
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/completeness.svg'))
log(f'figure: {Path(path).relative_to(ROOT)}')

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figures coloured by reliability.

Every point wears its outlier probability: blue is reliable, grey
undecided, red an outlier, and hollow when more likely an outlier than not.
set_style('web') gives SVG files for the site, set_style('paper') PDF
files; both styles are light.

Figure: a simulated planet at 11.2 d with bad visits, visit by visit after
a fit: the residual of each exposure, the mean of each visit, and below,
the probability that the visit is an outlier.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's figures, series, model, log (timestamped lines,
#   'value' for a number) and simulation
from pathlib import Path

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fit import RVModel
from koloa.log import log
from koloa.simulate import simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# a circular planet (P = 11.2 d, K = 5 m/s, tp its time of periastron
#   [days]) on a real NIRPS sampling, 7 % of the visits moved as a block
#   by 6 median errors or more, and a 2 m/s jitter shared by each visit
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
sim = simulate(planets=[dict(P=11.2, K=5.0, tp=60100.0)], template=tpl,
               outliers=[dict(kind='visit', frac=0.07, amplitude=6.0)],
               visit_jitter=2.0, seed=1)
# the series, and which exposures were made outliers (True: bad)
data, truth = sim['data'], sim['outlier_mask']
# a circular orbit started at 11.2 days, the outlier mixture with whole
#   visits as the outlier unit (unit='sequence'); the maximum a posteriori
fit = RVModel(data, [dict(period=11.2)], unit='sequence').fit(quiet=True)
# per visit: the highest outlier probability of its exposures, its mean
#   time [days], and whether it was made bad
visit_prob = np.array([fit.outlier_prob[data.seq == ss].max()
                       for ss in range(data.nseq)])
visit_time = np.array([data.time[data.seq == ss].mean()
                       for ss in range(data.nseq)])
visit_bad = np.array([truth[data.seq == ss].any()
                      for ss in range(data.nseq)])
# per exposure: the probability that it belongs to an outlier visit
prob = fit.outlier_prob
# how many exposures wear each colour of the figures
log(f'exposures: {np.sum(prob < 0.1)} blue (P < 0.1), '
    f'{np.sum((prob >= 0.1) & (prob <= 0.5))} in between, '
    f'{np.sum(prob > 0.5)} red and hollow (P > 0.5)', 'value')
# the flagged visits and the visits made bad, by their mean time [days]
log(f'visits more likely outliers than not: '
    f'{", ".join(f"{tt:.1f}" for tt in visit_time[visit_prob > 0.5])}',
    'value')
log(f'visits made bad: '
    f'{", ".join(f"{tt:.1f}" for tt in visit_time[visit_bad])}', 'value')

# the web style (this site's colours, SVG files); kplot.C holds the
#   colours of the current style
kplot.set_style('web')
log(f'web style: koloa {kplot.C["koloa"]}, gaussian {kplot.C["gaussian"]}, '
    f'reliable {kplot.C["reliable"]}, outlier {kplot.C["outlier"]}')
# also: timeseries, phase, periodograms, jackknife, coherence, corner,
#   recovery_map and detection_limits. highlight=None: no dates on the bars
# fit=fit shows each exposure's residual of the fit, not its velocity;
#   the panel below, the outlier probability of each visit
fig = kplot.sequences(data, fit.outlier_prob, fit=fit, highlight=None,
                      title='Simulation, every visit, after a fit at 11.2 d')
# a path without an extension: the web style adds .svg (paper: .pdf)
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/plotting'))
log(f'figure: {Path(path).relative_to(ROOT)} (the web style adds .svg)')

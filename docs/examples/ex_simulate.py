#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simulate a series with known planets and outliers.

simulate keeps the times, visits and error bars of a real series, adds
orbits, a visit jitter and outliers from ready-made recipes, and returns
the truth next to the data.

Figure: the simulated series on a real NIRPS sampling. The injected
outliers are drawn hollow and red, every other exposure blue.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's figures, series, log (timestamped lines, 'value' for
#   a number) and simulation
from pathlib import Path

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.log import log
from koloa.simulate import REALISTIC, simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# REALISTIC = CLEAR + BORDERLINE: bad visits and spikes, far and near
# each recipe is a dict: kind ('visit': the exposures of a visit move
#   together; 'spike': one exposure alone), frac (of the visits, or of the
#   exposures), amplitude in median errors (a number: that times 1 + E, E
#   exponential of mean 0.5; a (low, high) pair: uniform in it) and label
for rec in REALISTIC:
    amp, unit = rec['amplitude'], ('visits' if rec['kind'] == 'visit'
                                   else 'exposures')
    # np.ndim is 0 for a number, 1 for a (low, high) pair
    size = f'{amp[0]} to {amp[1]}' if np.ndim(amp) else f'{amp}+'
    log(f'recipe {rec["label"]:10s} {rec["kind"]:5s} {rec["frac"]:.0%} of '
        f'the {unit}, off by {size} median errors')

# a real sampling and its error bars make the most honest injection
# template: the times, visits, error bars and instrument of a real NIRPS
#   series are kept, and its velocities (all zero in this file) replaced
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
# an eccentric planet: P [days], K [m/s], e and omega [rad], and with no
#   tp the periastron falls on the first time of the series; the recipes
#   above; a 2 m/s jitter shared by the exposures of each visit; seed: the
#   random draws; name: the name of the series
sim = simulate(planets=[dict(P=16.65, K=5.0, e=0.2, omega=1.0)],
               template=tpl, outliers=REALISTIC, visit_jitter=2.0, seed=10,
               name='NIRPS sampling')
# the truth next to the data: which exposures were made outliers
data, mask = sim['data'], sim['outlier_mask']
log(f'{data.n} exposures in {data.nseq} visits, {mask.sum()} made outliers',
    'value')
# outlier_label: the recipe of each exposure ('' when good), and
#   outlier_offset: what was added to it [m/s]
for label in ('clear', 'borderline'):
    sel = sim['outlier_label'] == label
    log(f'  {label:10s} {sel.sum():2d} exposures in '
        f'{len(np.unique(data.seq[sel]))} visits, median shift '
        f'{np.median(np.abs(sim["outlier_offset"][sel])):.1f} m/s', 'value')
# signal: the velocity of the orbits alone at each exposure [m/s]
log(f'rms of the orbit {np.std(sim["signal"]):.2f} m/s, of the data '
    f'{np.std(data.rv):.2f} m/s', 'value')

# no template: seasons, one visit a night, two exposures back to back,
#   and a flare that decays over the next visits
# err: the error bar [m/s], each drawn between 0.8 and 1.2 times it;
#   nvisits, per_visit and baseline [days] shape the campaign; a flare
#   starts at 5 % of the visits, 6 median errors or more, and fades with
#   an e-folding time of 0.5 to 3 days
free = simulate(planets=[dict(P=5.2, K=3.0)], err=1.5, nvisits=40,
                per_visit=2, baseline=700,
                outliers=[dict(kind='flare', frac=0.05, amplitude=6.0)],
                seed=3)
log(f'drawn campaign: {free["data"].n} exposures in {free["data"].nseq} '
    f'visits over {free["data"].baseline:.0f} days, '
    f'{free["outlier_mask"].sum()} exposures in flares', 'value')

# koloa's web style; prob colours each point by its outlier probability,
#   here the truth (0 or 1): the injected outliers are red and hollow
kplot.set_style('web')
fig = kplot.timeseries(data, prob=mask.astype(float),
                       title='Injected outliers, hollow and red')
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/simulate.svg'))
log(f'figure: {Path(path).relative_to(ROOT)}')

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The duck test.

duck_test asks of a signal every question koloa knows: significance with
outliers accounted for, robustness to single visits, coherence over the
campaign, shape, and power in the activity indicators. Here on a
simulated planet at 11.2 d with bad visits, without the GP absorption
test, which is the slow one (gp=True adds it); the simulation has no
activity indicators, so their check has nothing to find.

Created on 2026-09-27

@author: artigau
"""
# koloa's series, the duck test and one of its checks, the two FIPs, the
#   log (timestamped lines, 'value' for a number) and the simulation
from pathlib import Path

from koloa.data import RVData
from koloa.diagnostics import coherence, duck_test
from koloa.fip import fip_single, oafip
from koloa.log import log
from koloa.simulate import simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]
# the period of the signal put to the test [days]
PERIOD = 11.2

# a circular planet at PERIOD (K = 5 m/s, tp its time of periastron
#   [days]) on a real NIRPS sampling, 7 % of the visits moved as a block
#   by 6 median errors or more, and a 2 m/s jitter shared by each visit
# (a real series also brings its activity indicators, in which the duck
#   test looks for power at the period; ex_data.py reads some)
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
sim = simulate(planets=[dict(P=PERIOD, K=5.0, tp=60100.0)], template=tpl,
               outliers=[dict(kind='visit', frac=0.07, amplitude=6.0)],
               visit_jitter=2.0, seed=1)
data = sim['data']
# the one-signal FIP with a fixed jitter, and koloa's (short chain)
# fip_single: the jitter is the excess of the robust rms over the median
#   error, and nothing is sampled; oafip: at most kmax=1 signal, whole
#   visits as the outliers (outliers='sequence'), one chain of 100 burn-in
#   and 400 recorded sweeps (defaults: two chains of 300 and 1500), no
#   progress lines
gauss = fip_single(data)
koloa = oafip(data, kmax=1, outliers='sequence', nsweep=400, nburn=100,
              nchains=1, progress=False)

# every check at PERIOD: fipres, the outlier-aware FIP, gives the
#   significance and the outlier probability of each exposure; gauss_fip
#   is quoted beside it; gp=False skips the GP absorption test (slow);
#   unit='sequence': whole visits are the outliers of its fits and the
#   units of its jackknife; quiet=True: no log lines, the report is
#   printed below
report = duck_test(data, PERIOD, fipres=koloa, gauss_fip=gauss, gp=False,
                   unit='sequence', quiet=True)
# the report line by line, the verdict in blue
for line in report.text().split('\n'):
    if line:
        log(line, 'value' if line.startswith('VERDICT') else 'info')
# checks: every check made; flags: those that speak against a planet
flagged = [check['name'] for check in report.flags]
log(f'checks that flag the signal: {", ".join(flagged) or "none"} (of '
    f'{len(report.checks)})', 'value')

# one of the checks on its own: amplitude and phase, season by season
# a sinusoid at PERIOD is fitted to each season (a new one after a gap of
#   60 days) on the visit means, each exposure weighted by one minus its
#   outlier probability (prob)
coh = coherence(data, PERIOD, prob=koloa.outlier_prob, split='seasons')
# each chunk: its first time tmin [days], its visits, K and its error sK
#   [m/s]
for chunk in coh['chunks']:
    log(f'  season from {chunk["tmin"]:.0f}: {chunk["nvisits"]:2d} visits, '
        f'K = {chunk["K"]:.1f} +- {chunk["sK"]:.1f} m/s', 'value')
# p_vector: the p-value of a chi^2 test that amplitude and phase are the
#   same in every season (a small p: they differ)
log(f'same amplitude and phase in every season: p = {coh["p_vector"]:.1e}',
    'value')

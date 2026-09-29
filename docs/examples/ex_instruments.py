#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Several instruments: one planet, a noise of their own.

A near-infrared instrument (visits of three exposures) and an optical one
(two a night) see the same planet (P = 9.3 d, K = 4 m/s), and the activity
of the star, weaker in the near infrared: a Matern-3/2 GP of 1.5 m/s there
and of 4 m/s in the optical, the same for the exposures of a visit. Two
near-infrared visits are bad, by 8 m/s.

Three fits of both series together: one visit jitter for every instrument
(the optical activity sets it, and the bad near-infrared visits pass for
good), a visit jitter per instrument (seq_jitter='instrument'), and a GP
per instrument (a list of GPs, each on its own instruments). Maximum a
posteriori fits.

Figure: each instrument folded at the period, its GP removed.

Created on 2026-09-28

@author: artigau
"""
# numpy, and koloa's figures, series, model and priors, log (timestamped
#   lines, 'value' for a number) and simulation
from pathlib import Path

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fit import Prior, RVModel
from koloa.log import log
from koloa.simulate import activity_signal, keplerian_signal, observing_times

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]
# the random numbers of the photon noise
rng = np.random.default_rng(3)
# one circular orbit: P and tp, the time of periastron [days], and K [m/s]
planet = [dict(P=9.3, K=4.0, tp=0.0)]
# each instrument: its name, visits, exposures per visit, error bar
#   [m/s], amplitude of the activity [m/s] and seed
series = []
for name, nvis, nexp, err, sigma, seed in (('NIR', 50, 3, 2.0, 1.5, 1),
                                          ('VIS', 60, 2, 1.0, 4.0, 2)):
    # a drawn campaign over 600 days: seasons, one visit a night, nexp
    #   exposures 15 minutes apart; the times in days
    time = observing_times(nvisits=nvis, per_visit=nexp, baseline=600,
                           seed=seed)
    # one draw of a Matern-3/2 GP, pars the logs of its sigma [m/s] and of
    #   its length [days]; at 8 days, the exposures of a visit share it
    act = activity_signal(time, kernel='matern32', seed=seed + 10,
                          pars=np.log([sigma, 8.0]))
    # the orbit, the activity and the photon noise [m/s]
    rv = keplerian_signal(time, planet) + act + err * rng.normal(size=len(time))
    series.append((time, rv, np.full(len(time), err), np.full(len(time), name)))
# both instruments in one series, in time order; RVData finds the visits
#   of each instrument and takes its median velocity out
time, rv, err, inst = (np.concatenate(parts) for parts in zip(*series))
order = np.argsort(time)
data = RVData(time[order], rv[order], err[order], inst=inst[order])
# two bad near-infrared visits: every exposure moved by 8 m/s
# (the 11th and the 31st visits of the NIR instrument)
nir_visits = np.unique(data.seq[data.inst == 'NIR'])
for visit in nir_visits[[10, 30]]:
    data.rv[data.seq == visit] += 8.0
log(f'{data.n} exposures in {data.nseq} visits: '
    + ', '.join(f'{name} {np.sum(data.inst == name)}'
                for name in data.instruments), 'info')

# the orbit to fit: its period starts at 9.3 days, free from 9.0 to 9.6;
#   bad: the two visits made bad, to score the flags
period = dict(period=9.3, period_range=(9.0, 9.6))
bad = set(nir_visits[[10, 30]])


def report(label, fit):
    """K, the noise, and the near-infrared visits flagged"""
    model = fit.model
    # the jitters in m/s, from their logs: log_jit_<inst>, white, one per
    #   instrument; log_sjit (or log_sjit_<inst>), shared by a visit
    noise = ', '.join(f'{name[4:]} {np.exp(fit.theta[model.index[name]]):.1f}'
                      for name in model.names
                      if name.startswith(('log_sjit', 'log_jit')))
    # a visit is flagged when its exposures are more likely bad than not
    flagged = {visit for visit in nir_visits
               if np.mean(fit.outlier_prob[data.seq == visit]) > 0.5}
    # K as (value, minus, plus) [m/s]
    kval = fit.orbits()[0]['K']
    log(f'{label}: K = {kval[0]:.2f} +- {0.5 * (kval[1] + kval[2]):.2f} m/s;'
        f' jitters [m/s] {noise}; bad NIR visits flagged '
        f'{len(flagged & bad)} of 2, good ones {len(flagged - bad)}', 'value')


# every fit keeps the defaults: the outlier mixture, where a lone exposure
#   or a whole visit can be an outlier, a white jitter per instrument and
#   a straight trend; nstart=2: the number of starts (4 by default)
# seq_jitter=True: one visit jitter for both instruments
one = RVModel(data, [period], seq_jitter=True).fit(nstart=2, quiet=True)
report('one visit jitter', one)
# seq_jitter='instrument': a visit jitter per instrument (sjit_NIR, ...)
per = RVModel(data, [period], seq_jitter='instrument').fit(nstart=2,
                                                           quiet=True)
report('a visit jitter per instrument', per)
# a list of GPs, one per instrument, each with its own parameters
#   (gp_NIR_log_sigma, ...); a Prior on log_length keeps the length from
#   2 to 200 days, uniform in its log
gps = [dict(kernel='matern32', instruments=[name],
            prior={'log_length': Prior('uniform', np.log(2.0),
                                       np.log(200.0))})
       for name in ('NIR', 'VIS')]
# with a GP, the fit alternates the parameters and the outlier flags
#   (the slow part of this example)
withgp = RVModel(data, [period], seq_jitter='instrument',
                 gp=gps).fit(nstart=2, quiet=True)
report('and a GP per instrument', withgp)
# a few words on the model
log(withgp.model.describe(), 'info')
# the GP of each instrument, from its log parameters
for name in ('NIR', 'VIS'):
    sigma = np.exp(withgp.theta[withgp.model.index[f'gp_{name}_log_sigma']])
    length = np.exp(withgp.theta[withgp.model.index[f'gp_{name}_log_length']])
    log(f'  GP {name}: sigma = {sigma:.1f} m/s, length = {length:.0f} d '
        f'(simulated: {1.5 if name == "NIR" else 4.0} m/s, 8 d)', 'value')

# pyplot for the figure (noqa: E402 lets flake8 accept an import below code)
import matplotlib.pyplot as plt    # noqa: E402

# koloa's web style, and two panels side by side that share the velocity
#   axis (figsize in inches)
kplot.set_style('web')
fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.3), sharey=True)
for ax, name in zip(axes, ('NIR', 'VIS')):
    # each exposure (level='point') folded on the joint fit, everything but
    #   the orbit taken out; select: this instrument only; ylim [m/s]:
    #   points beyond it sit at its edge; no legend nor colour bar here
    kplot.phase(withgp, level='point', ax=ax, title=name,
                select=data.inst == name, ylim=(-15.0, 25.0), legend=False,
                colorbar=False)
axes[1].set_ylabel('')
# one colour bar of the outlier probability for both panels (a private
#   helper of koloa.plotting)
kplot._colorbar(fig, list(axes), kplot.reliability_cmap())
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/instruments.svg'))
log(f'figure: {Path(path).relative_to(ROOT)}')

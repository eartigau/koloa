#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Orbital elements by MCMC.

mcmc_orbits finds the maximum a posteriori, then samples the posterior with
emcee; the visits are allowed to be outliers. The result gives every draw
in physical units (posterior), their summaries (orbits), and the
autocorrelation times (diagnostics). Short chains here: the default runs
until the chain holds 50 autocorrelation times.

Figure: the phase fold of the orbit, with binned means and the 1-sigma
envelope of the posterior. The exposures of the flagged visits are hollow
and red.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's figures, series, model and MCMC, log (timestamped
#   lines, 'value' for a number) and simulation
from pathlib import Path

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fit import RVModel, mcmc_orbits
from koloa.log import log
from koloa.simulate import simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]
# the mass of the star, for m sin i
MSTAR = 0.6    # solar masses

# an eccentric planet (P = 11.2 d, K = 6 m/s, e = 0.3) and 7% bad visits
# on a real NIRPS sampling: omega [rad] is the argument of periastron of
#   the star and tp the time of periastron [days]; the bad visits move by
#   6 median errors or more, and a 2 m/s jitter is shared by each visit
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
sim = simulate(planets=[dict(P=11.2, K=6.0, e=0.3, omega=np.radians(60),
                             tp=60100.0)], template=tpl,
               outliers=[dict(kind='visit', frac=0.07, amplitude=6.0)],
               visit_jitter=2.0, seed=4)
data = sim['data']

# the maximum a posteriori, then the posterior with emcee: [11.2] holds
#   the starting period of each orbit [days] (free within +-10 %), each
#   eccentric by default; the outlier mixture, with a whole visit as the
#   outlier unit (unit='sequence'); nsteps and nburn: at least that many
#   steps per walker, and thrown away; the chain grows towards 50
#   autocorrelation times, up to max_steps (by default 20 x nsteps)
post = mcmc_orbits(data, [11.2], unit='sequence', nsteps=1500, nburn=500,
                   max_steps=3000, quiet=True)
# the model, each orbit with its errors, the noise and the sampler
for line in post.summary(mstar=MSTAR).split('\n'):
    log(line, 'value')
# orbits(): per orbit, (median, minus, plus) of P, K, e, omega [deg], tp,
#   tc and, given the stellar mass, m sin i [Earth masses]
orbit = post.orbits(mstar=MSTAR)[0]
log(f'omega = {orbit["omega"][0]:.0f} +{orbit["omega"][2]:.0f}/'
    f'-{orbit["omega"][1]:.0f} deg (injected: P = 11.2 d, K = 6 m/s, '
    f'e = 0.3, omega = 60 deg)', 'value')
# posterior(): every draw kept (one step in 10 of each walker) in physical
#   units, named by orbit (K_0, e_0, ...); entry j of each is one draw
draws = post.posterior(mstar=MSTAR)
log(f'{len(draws["K_0"])} draws; P(e > 0.1) = '
    f'{np.mean(draws["e_0"] > 0.1):.2f}, P(K > 5 m/s) = '
    f'{np.mean(draws["K_0"] > 5):.2f}', 'value')
# diagnostics: the walkers, the steps kept after the burn-in, and the
#   longest autocorrelation time, tau_max [steps]
diag = post.diagnostics
log(f'{diag["nwalkers"]} walkers, {diag["nkept"]} steps kept; 50 '
    f'autocorrelation times would take {50 * diag["tau_max"]:.0f}', 'value')

# the same orbit without outliers in the noise model (maximum a posteriori)
# likelihood='gaussian': no outlier; nstart=1: one start per candidate
#   period (the default 4 also turns the periastron around); K comes as
#   (value, minus, plus) [m/s], the errors from the curvature at the
#   maximum
gauss = RVModel(data, [dict(period=11.2, eccentric=True)],
                likelihood='gaussian').fit(nstart=1, quiet=True)
kval = gauss.orbits()[0]['K']
# the error on K with outliers allowed for: the mean of the two sides of
#   its posterior interval
kerr = 0.5 * (orbit['K'][1] + orbit['K'][2])
log(f'gaussian: K = {kval[0]:.2f} +- {kval[1]:.2f} m/s, e = '
    f'{gauss.orbits()[0]["e"][0]:.2f}: the bad visits pull the orbit, '
    f'and the error on K is {kval[1] / kerr:.1f} times that with outliers '
    f'allowed for', 'value')

# koloa's web style; the phase fold of every exposure (level='point'),
#   with binned means and the 1-sigma envelope over the draws
kplot.set_style('web')
fig = kplot.phase(post, level='point')
fig.axes[0].legend(loc='lower right')    # a bad visit sits upper right
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/fit.svg'))
log(f'figure: {Path(path).relative_to(ROOT)}')

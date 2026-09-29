#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Periodograms that know about outliers.

A planet of K = 5 m/s at 16.65 d behind bad visits and spikes. The GLS and
the gaussian profile periodogram peak elsewhere; the outlier-aware
periodogram (OAP) peaks at the planet. The window, the aliases and a
leave-one-visit-out jackknife say how far to trust a peak.

Figure: the gain in log likelihood of a sinusoid, gaussian (orange) and
outlier-aware (blue), and the window. The planet's period is dashed.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's figures, series, log (timestamped lines, 'value' for
#   a number), periodograms and simulation
from pathlib import Path

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.log import log
from koloa.periodogram import (aliases, find_peaks, frequency_grid, gls,
                               influence, jackknife, oap, window)
from koloa.simulate import simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]
# the period of the injected planet [days]
PERIOD = 16.65

# the times, visits and error bars of a real NIRPS sampling, with a
#   circular planet (K [m/s], tp its time of periastron [days]), 7 % of
#   the visits moved as a block by 6 median errors or more, 3 % of the
#   exposures alone by 8 or more, and a 2 m/s jitter shared by each visit
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
sim = simulate(planets=[dict(P=PERIOD, K=5.0, tp=60105.3)], template=tpl,
               outliers=[dict(kind='visit', frac=0.07, amplitude=6.0),
                         dict(kind='spike', frac=0.03, amplitude=8.0)],
               visit_jitter=2.0, seed=10)
data = sim['data']

# periods from 1.1 d to twice the baseline T, 10 grid points per 1/T
freq = frequency_grid(data.time, pmin=1.1, oversample=10)
# the generalised Lomb-Scargle periodogram, a power from 0 to 1
power = gls(data.time, data.rv, data.err, freq)
# one point per visit, a free jitter, with and without the outlier mixture
# (each gives dlnl, the gain in ln L of an offset and a sinusoid over the
#   offset alone, and the fitted frac, jitter [m/s] and amp [m/s])
gauss = oap(data, freq, unit='sequence', outliers=False)
mixt = oap(data, freq, unit='sequence', outliers=True)
for label, score in (('GLS', power), ('profile, Gaussian', gauss['dlnl']),
                     ('OAP, outlier-aware', mixt['dlnl'])):
    # the grid index of the highest peak
    best = find_peaks(freq, score, npeaks=1)[0]
    log(f'{label:19s} best peak at {1 / freq[best]:7.3f} d', 'value')
# the grid point nearest the planet's frequency
idx = np.argmin(np.abs(freq - 1 / PERIOD))
log(f'at {PERIOD} d the OAP finds an outlier fraction of '
    f'{mixt["frac"][idx]:.2f} and a jitter of {mixt["jitter"][idx]:.1f} m/s',
    'value')

# the spectral window: what the sampling alone puts at each frequency (1
#   at zero frequency)
wpow = window(data.time, freq)
log(f'window power at {PERIOD} d: {wpow[idx]:.3f} (highest '
    f'{wpow.max():.2f}, at {1 / freq[np.argmax(wpow)]:.0f} d)', 'value')
# where else the planet shows: its aliases with the sidereal day, the year
#   and the synodic month, those inside the grid
near = ', '.join(f'{al["period"]:.2f} d ({al["name"]})'
                 for al in aliases(1 / PERIOD, freq[0], freq[-1]))
log(f'aliases of {PERIOD} d: {near}')

# which visits hold up a peak: at the GLS peak and at the planet
# jackknife: the GLS with each visit left out in turn; influence: the
#   power with a visit minus the power without it, here over the power
#   of the peak (positive: the visit holds the peak up)
jack = jackknife(data, freq, unit='sequence')
# the visits that hold an injected outlier
bad = np.unique(data.seq[sim['outlier_mask']])
for period in (1 / freq[np.argmax(power)], PERIOD):
    infl = influence(jack, freq, period) / jack['power'][
        np.argmin(np.abs(freq - 1 / period))]
    log(f'GLS power at {period:6.3f} d held up by the top visit: '
        f'{infl.max():.0%}, by the {len(bad)} with outliers: '
        f'{infl[bad].sum():+.0%}', 'value')

# koloa's web style; the two profile periodograms on the same Delta ln L
#   axis, the window below, and the periods of mark dashed [days]
kplot.set_style('web')
fig = kplot.periodograms(freq, oap_gauss=gauss['dlnl'], oap_mix=mixt['dlnl'],
                         window_power=wpow, mark=[PERIOD],
                         title='A planet at 16.65 d behind bad visits')
fig.axes[0].legend(loc='upper right', ncol=2)    # clear of the peak
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/periodogram.svg'))
log(f'figure: {Path(path).relative_to(ROOT)}')

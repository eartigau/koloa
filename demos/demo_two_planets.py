#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 4: two planets at once, with outliers.

The single-signal FIP looks for one sinusoid: a second planet is noise to
it, and it cannot say how many signals there are. koloa's sampler has up
to kmax signal slots, switched on and off, so the FIP periodogram shows
every planet at once, and the posterior of the number of signals comes
with it (Hara et al. 2022, made outlier-aware).

    python demo_two_planets.py

Created on 2026-09-27

@author: artigau
"""
import numpy as np

from _common import (flag_summary, log_flags, outdir, save_both,
                     save_summary, template)
from koloa import plotting as kplot
from koloa.fip import fip_single, oafip
from koloa.fit import RVModel
from koloa.log import log
from koloa.periodogram import frequency_grid, oap, window
from koloa.simulate import REALISTIC, simulate

DEMO = 'two_planets'
PLANETS = [dict(P=4.31, K=3.5, e=0.0, tp=60100.0),
           dict(P=31.4, K=4.5, e=0.25, omega=1.2, tp=60110.0)]


def main():
    log('Demo: two planets, outliers clear and borderline')
    tpl = template()
    # outliers clear and borderline (koloa.simulate.REALISTIC)
    sim = simulate(planets=PLANETS, template=tpl, outliers=REALISTIC,
                   visit_jitter=1.5, seed=21, name='two planets')
    data = sim['data']
    freq = frequency_grid(data.time, 1.1, None, 10)
    single = fip_single(data, freq=freq)
    koloa = oafip(data, kmax=3, outliers='both', freq=freq, nsweep=2000,
                  nburn=400, nchains=2)
    log('  ' + ', '.join(f'P(k={it}) = {val:.3f}' for it, val in
                         enumerate(koloa.pk)), 'value')
    for pl in PLANETS:
        log(f'  planet {pl["P"]} d: FIP single-signal gaussian '
            f'{single.fip_at(pl["P"]):.1e}, koloa {koloa.fip_at(pl["P"]):.1e}',
            'value')
    fit = RVModel(data, [dict(period=PLANETS[0]['P']),
                         dict(period=PLANETS[1]['P'], eccentric=True)],
                  likelihood='mixture', unit='both').fit()
    for line in fit.summary().split('\n')[:3]:
        log(line, 'value')
    flags = flag_summary(koloa.outlier_prob, sim['outlier_label'])
    log_flags(flags)
    orbits = fit.orbits()
    save_summary(DEMO, dict(
        pk=[float(v) for v in koloa.pk],
        planets=[dict(P=pl['P'], K=pl['K'], fip_single=single.fip_at(pl['P']),
                      fip_koloa=koloa.fip_at(pl['P']),
                      K_fit=orbits[ip]['K'], P_fit=orbits[ip]['P'])
                 for ip, pl in enumerate(PLANETS)], flags=flags))
    og = oap(data, freq, unit='sequence', outliers=False)
    om = oap(data, freq, unit='sequence', outliers=True)
    save_both(lambda: kplot.periodograms(
        freq, oap_gauss=og['dlnl'], oap_mix=om['dlnl'],
        window_power=window(data.time, freq),
        fips=dict(gaussian=single, koloa=koloa),
        mark=[pl['P'] for pl in PLANETS],
        title='Two planets (dashed): one FIP periodogram shows both'),
        'periodograms', DEMO)
    for ip in range(2):
        save_both(lambda ip=ip: kplot.phase(fit, planet=ip), f'phase_{ip}',
                  DEMO)
    log(f'Figures in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

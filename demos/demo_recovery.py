#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 2: a planet that IS there, behind bad visits and spikes.

A planet of K = 5 m/s at P = 16.65 d, on the sampling of a real NIRPS
series, with a visit jitter of 2 m/s, a few bad visits and a few spiked
exposures. The gaussian FIP chases the outliers and puts its best peak
elsewhere; the clips trade the outliers for sensitivity; koloa recovers
the planet, flags every bad exposure, and fits the orbit without them.

    python demo_recovery.py

Created on 2026-09-27

@author: artigau
"""
import numpy as np

from _common import (flag_summary, log_flags, outdir, save_both,
                     save_summary, template, web_fips)
from koloa import plotting as kplot
from koloa.fip import fip_comparison, oafip
from koloa.fit import RVModel
from koloa.log import log
from koloa.periodogram import frequency_grid, oap, window
from koloa.simulate import REALISTIC, simulate

DEMO = 'recovery'
SEED = 10


def main():
    log('Demo: a planet behind bad visits and spikes')
    tpl = template()
    rng = np.random.default_rng(SEED + 500)
    period = float(np.exp(rng.uniform(np.log(3), np.log(60))))
    tperi = float(tpl.time[0] + rng.uniform(0, period))
    sim = simulate(planets=[dict(P=period, K=5.0, tp=tperi)], template=tpl,
                   outliers=REALISTIC, visit_jitter=2.0, seed=SEED,
                   name='planet, outliers clear and borderline')
    data = sim['data']
    truth = sim['outlier_mask']
    log(f'  injected: P = {period:.3f} d, K = 5.0 m/s; '
        f'{int(np.sum(truth))} outlying exposures', 'value')
    freq = frequency_grid(data.time, 1.1, None, 10)
    comp = fip_comparison(data, freq=freq)
    koloa = oafip(data, kmax=2, outliers='both', freq=freq, nsweep=1500,
                  nburn=300, nchains=2)
    for key, res in list(comp.items()) + [('koloa', koloa)]:
        log(f'  {res.method:36s} FIP at the planet {res.fip_at(period):.1e}'
            f'; best {res.best()["period"]:.3f} d', 'value')
    flags = flag_summary(koloa.outlier_prob, sim['outlier_label'])
    log_flags(flags)
    # the orbit, with and without knowing about outliers
    fit_g = RVModel(data, [dict(period=period)], likelihood='gaussian',
                    seq_jitter=True).fit()
    fit_k = RVModel(data, [dict(period=period)], likelihood='mixture',
                    unit='both').fit()
    orbits = {}
    for name, fit in (('gaussian', fit_g), ('koloa', fit_k)):
        orb = fit.orbits()[0]
        orbits[name] = dict(K=orb['K'], P=orb['P'])
        log(f'  {name:9s} K = {orb["K"][0]:.2f} +{orb["K"][2]:.2f}/'
            f'-{orb["K"][1]:.2f} m/s, P = {orb["P"][0]:.3f} d', 'value')
    save_summary(DEMO, dict(
        period=period, K=5.0, flags=flags, orbits=orbits,
        fips={key: dict(at_planet=res.fip_at(period),
                        best_period=res.best()['period'],
                        best_fip=res.best()['fip'])
              for key, res in list(comp.items()) + [('koloa', koloa)]}))
    og = oap(data, freq, unit='sequence', outliers=False)
    om = oap(data, freq, unit='sequence', outliers=True)
    fips = dict(gaussian=comp['gaussian'], soft=comp['soft'],
                hard=comp['hard'], koloa=koloa)
    save_both(lambda: kplot.periodograms(
        freq, oap_gauss=og['dlnl'], oap_mix=om['dlnl'],
        window_power=window(data.time, freq), fips=web_fips(fips), mark=[period],
        title=f'A planet at {period:.2f} d (dashed) behind bad visits and '
              f'spikes'), 'periodograms', DEMO)
    save_both(lambda: kplot.phase(
        fit_g, title='Gaussian fit: K = {:.2f} m/s'.format(
            fit_g.orbits()[0]['K'][0])), 'phase_gaussian', DEMO)
    save_both(lambda: kplot.phase(
        fit_k, title='koloa fit: K = {:.2f} m/s (injected 5.00)'.format(
            fit_k.orbits()[0]['K'][0])), 'phase_koloa', DEMO)
    save_both(lambda: kplot.phase(
        fit_k, level='sequence', title='koloa fit, one point per visit'),
        'phase_koloa_visits', DEMO)
    save_both(lambda: kplot.sequences(
        data, fit_k.outlier_prob, fit=fit_k,
        title='Every visit, after the orbit: bad visits and lone spikes'),
        'visits', DEMO)
    log(f'Figures in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

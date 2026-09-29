#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 1: a planet that is not there.

Series with the sampling and error bars of a real NIRPS series, a visit
jitter of 2 m/s, NO planet, and outliers of both signs: clear ones (whole
visits and single exposures well away from the bulk) and borderline ones
(two to four sigma), as koloa.simulate.REALISTIC makes them. Ten
realisations are drawn, and each is analysed with the single-signal FIP
with a fixed jitter (the usual computation), after a soft and a hard clip,
and with koloa. The demo counts how often each reports a planet (FIP < 1 %)
and shows the first realisation in which the fixed-jitter FIP does; how
often that happens over many more series is measured by the injection
campaign of the paper (paper/).

    python demo_false_alarm.py

Created on 2026-09-27

@author: artigau
"""
import numpy as np

from _common import (flag_summary, log_flags, outdir, save_both,
                     save_summary, template, web_fips)
from koloa import plotting as kplot
from koloa.fip import fip_comparison, oafip
from koloa.log import log
from koloa.periodogram import frequency_grid, oap, window
from koloa.simulate import REALISTIC, simulate

DEMO = 'false_alarm'
#: realisations: how often each method is fooled, and one of them shown
NSEED = 10


def realisation(seed: int):
    """No planet: noise, a visit jitter, outliers clear and borderline"""
    tpl = template()
    sim = simulate(template=tpl, outliers=REALISTIC, visit_jitter=2.0,
                   seed=seed, name='no planet, outliers')
    data = sim['data']
    freq = frequency_grid(data.time, 1.1, None, 10)
    comp = fip_comparison(data, freq=freq)
    koloa = oafip(data, kmax=2, outliers='both', freq=freq, nsweep=1500,
                  nburn=300, nchains=2, progress=False)
    return sim, data, freq, comp, koloa


def main():
    log('Demo: a planet that is not there')
    rows, shown = [], None
    for seed in range(1, NSEED + 1):
        sim, data, freq, comp, koloa = realisation(seed)
        row = dict(seed=seed, noutlier=int(np.sum(sim['outlier_mask'])),
                   flags=flag_summary(koloa.outlier_prob,
                                      sim['outlier_label']))
        for key, res in list(comp.items()) + [('koloa', koloa)]:
            best = res.best()
            row[key] = dict(period=float(best['period']),
                            fip=float(best['fip']))
        rows.append(row)
        log(f'  seed {seed}: best FIP ' + ', '.join(
            f'{key} {row[key]["fip"]:.1e} at {row[key]["period"]:.2f} d'
            for key in ('gaussian', 'soft', 'hard', 'koloa')), 'value')
        # the first realisation that fools the fixed-jitter FIP is shown
        if shown is None and row['gaussian']['fip'] < 0.01:
            shown = (seed, sim, data, freq, comp, koloa)
    if shown is None:
        shown = (1,) + realisation(1)
    seed, sim, data, freq, comp, koloa = shown
    fooled = {key: int(sum(row[key]['fip'] < 0.01 for row in rows))
              for key in ('gaussian', 'soft', 'hard', 'koloa')}
    log('  realisations with a false detection (FIP < 1 %): ' + ', '.join(
        f'{key} {val}/{NSEED}' for key, val in fooled.items()), 'value')
    flags = flag_summary(koloa.outlier_prob, sim['outlier_label'])
    log(f'  shown: seed {seed}', 'value')
    log_flags(flags)
    labels = sim['outlier_label']
    nvisit = {kind: len(np.unique(data.seq[labels == kind]))
              for kind in ('clear', 'borderline')}
    save_summary(DEMO, dict(nseed=NSEED, fooled=fooled, rows=rows,
                            shown=dict(seed=seed, flags=flags,
                                       nvisit=nvisit,
                                       **{key: rows[seed - 1][key] for key in
                                          ('gaussian', 'soft', 'hard',
                                           'koloa')})))
    og = oap(data, freq, unit='sequence', outliers=False)
    om = oap(data, freq, unit='sequence', outliers=True)
    wpow = window(data.time, freq)
    fips = dict(gaussian=comp['gaussian'], soft=comp['soft'],
                hard=comp['hard'], koloa=koloa)
    fake = comp['gaussian'].best()['period']
    save_both(lambda: kplot.periodograms(
        freq, oap_gauss=og['dlnl'], oap_mix=om['dlnl'], window_power=wpow,
        fips=web_fips(fips), mark=[fake],
        title='No planet, outliers clear and borderline: who finds a '
              'planet?'), 'periodograms', DEMO)
    save_both(lambda: kplot.timeseries(
        data, koloa.outlier_prob,
        title='The simulated series, coloured by koloa\'s outlier '
              'probability'), 'timeseries', DEMO)
    save_both(lambda: kplot.sequences(
        data, koloa.outlier_prob, title='Every visit: the bad ones stand '
        'out as a whole, the borderline ones in grey'), 'visits', DEMO)
    log(f'Figures in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

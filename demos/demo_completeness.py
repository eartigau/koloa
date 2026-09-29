#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 10: recovery rates in period and K, with outlier-aware folding.

Circular planets are injected over a grid of periods and semi-amplitudes and
looked for by folding the series at the injected period, with koloa's
outlier-aware likelihood and with a gaussian one (koloa.completeness).

1. The method: a new series for every injection, on a real NIRPS sampling
   with a 2 m/s visit jitter and outliers, clear and borderline
   (koloa.simulate.REALISTIC), and the same series without outliers. The
   fold test (a known period) on the whole grid, and the blind search test
   on a coarser one.
2. Two real series as observed: GL 725 B (SPIRou) and Kepler-21
   (HARPS-N), with a linear trend in the fold and the threshold
   of each period bin calibrated on the series' own empty injections (its
   red noise makes the chi^2 threshold far too low at long periods).

    python demo_completeness.py --workers 8

Created on 2026-09-27

@author: artigau
"""
import argparse
import json
import os
import warnings

import numpy as np

from _common import DATA, ROOT, outdir, save_both
from koloa import plotting as kplot
from koloa.completeness import recovery_map
from koloa.data import RVData
from koloa.log import log
from koloa.simulate import REALISTIC, simulate

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'completeness'
PERIODS = np.geomspace(1.5, 500.0, 13)
AMPS = np.geomspace(0.5, 20.0, 13)
PERIODS_SEARCH = np.geomspace(1.5, 500.0, 7)
AMPS_SEARCH = np.geomspace(0.5, 20.0, 7)
NINJ, NINJ_SEARCH, NNULL = 40, 12, 400
REAL = [('GL 725 B (SPIRou)', os.path.join(ROOT, 'data',
                                          'gl725b_spirou_lbl_2026-09-24.rdb')),
        ('Kepler-21 (HARPS-N)', os.path.join(
            ROOT, 'data', 'kepler21_harpsn_dace_drs3.3.12.csv'))]


# =============================================================================
# Define functions
# =============================================================================
def contaminated(seed: int) -> RVData:
    """A new series on the NIRPS sampling, outliers clear and borderline"""
    tpl = RVData.from_csv(DATA)
    return simulate(template=tpl, outliers=REALISTIC, visit_jitter=2.0,
                    seed=seed)['data']


def clean(seed: int) -> RVData:
    """The same series without its outliers"""
    tpl = RVData.from_csv(DATA)
    return simulate(template=tpl, visit_jitter=2.0, seed=seed)['data']


def describe(rmap) -> dict:
    """What the page and the paper need from a map"""
    return dict(periods=rmap.periods.tolist(), k50=rmap.k_at(0.5).tolist(),
                k90=rmap.k_at(0.9).tolist(), false_alarm=rmap.false_alarm,
                nnull=rmap.nnull, threshold=rmap.threshold,
                thresholds=rmap.extra.get('thresholds'),
                false_alarm_nominal=rmap.extra.get('false_alarm_nominal'),
                calibrated=rmap.extra.get('calibrated', False),
                ninj=int(np.sum(rmap.counts)))


class SavedMap:
    """What the limits figure needs of a map, from the summary"""

    def __init__(self, entry: dict):
        self.periods = np.array(entry['periods'])
        self._k = {0.5: np.array(entry['k50'], dtype=float),
                   0.9: np.array(entry['k90'], dtype=float)}

    def k_at(self, level: float) -> np.ndarray:
        return self._k[level]


def limits_figure(maps: dict):
    """The semi-amplitude recovered half of the time, every map"""
    return kplot.detection_limits(
        [maps['gaussian'], maps['koloa'], maps['gaussian_clean'],
         maps['koloa_clean'], maps['search_gaussian'], maps['search_koloa']],
        ['Gaussian fold', 'koloa fold', 'Gaussian fold (no outliers)',
         'koloa fold (no outliers)', 'Gaussian blind search',
         'koloa blind search'],
        ['gaussian', 'koloa', 'gaussian', 'koloa', 'hard', 'visit'],
        title='K recovered half of the time')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--limits-only', action='store_true',
                        help='only draw the limits figure again, from the '
                             'summary')
    args = parser.parse_args()
    if args.limits_only:
        with open(os.path.join(outdir(DEMO), 'summary.json')) as handle:
            summ = json.load(handle)
        maps = {key: SavedMap(summ[key]) for key in (
            'gaussian', 'koloa', 'gaussian_clean', 'koloa_clean',
            'search_gaussian', 'search_koloa')}
        save_both(lambda: limits_figure(maps), 'limits', DEMO)
        return
    warnings.simplefilter('ignore')
    log('Demo: recovery rates in period and K, outlier-aware folding')
    maps, summary = {}, {}
    # -------------------------------------------------------------------------
    # 1. the method, on simulated series
    # -------------------------------------------------------------------------
    for key, series, outl in (('koloa', contaminated, True),
                              ('gaussian', contaminated, False),
                              ('koloa_clean', clean, True),
                              ('gaussian_clean', clean, False)):
        maps[key] = recovery_map(series, PERIODS, AMPS, ninj=NINJ,
                                 test='fold', outliers=outl, nnull=NNULL,
                                 seed=11, workers=args.workers)
        summary[key] = describe(maps[key])
        log(maps[key].summary())
    for key, outl in (('search_koloa', True), ('search_gaussian', False)):
        maps[key] = recovery_map(contaminated, PERIODS_SEARCH, AMPS_SEARCH,
                                 ninj=NINJ_SEARCH, test='search',
                                 outliers=outl, nnull=NNULL // 4, seed=12,
                                 workers=args.workers)
        summary[key] = describe(maps[key])
        log(maps[key].summary())
    # -------------------------------------------------------------------------
    # 2. real series, as observed
    # -------------------------------------------------------------------------
    for label, path in REAL:
        data = RVData.from_csv(path, name=label)
        key = 'real_' + label.split(' (')[0].replace(' ', '').lower()
        # a real series has drifts and red noise (activity, systematics):
        #   a linear trend in both models, and the threshold of each period
        #   bin set by that bin's own empty injections
        maps[key] = recovery_map(data, PERIODS, AMPS, ninj=NINJ // 2,
                                 test='fold', outliers=True, nnull=NNULL // 2,
                                 seed=13, workers=args.workers, trend=1,
                                 calibrate=True)
        maps[key + '_gauss'] = recovery_map(data, PERIODS, AMPS,
                                            ninj=NINJ // 2, test='fold',
                                            outliers=False, nnull=NNULL // 2,
                                            seed=13, workers=args.workers,
                                            trend=1, calibrate=True)
        summary[key] = dict(describe(maps[key]), label=label,
                            gaussian=describe(maps[key + '_gauss']),
                            npoints=data.n, nvisits=data.nseq,
                            baseline=data.baseline)
        log(f'{label}:\n' + maps[key].summary())
    with open(os.path.join(outdir(DEMO), 'summary.json'), 'w') as handle:
        json.dump(summary, handle, indent=1, default=float)
    # -------------------------------------------------------------------------
    # figures
    # -------------------------------------------------------------------------
    def method_figure():
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.1), sharey=True)
        for ax, (key, title) in zip(axes, (
                ('gaussian', 'Gaussian fold, with outliers'),
                ('koloa', 'koloa fold, with outliers'),
                ('koloa_clean', 'koloa fold, no outliers'))):
            kplot.recovery_map(maps[key], title=title, ax=ax,
                               colorbar=(key == 'koloa_clean'))
        for ax in axes[1:]:
            ax.set_ylabel('')
        fig.tight_layout()
        return fig

    save_both(method_figure, 'maps', DEMO)
    save_both(lambda: limits_figure(maps), 'limits', DEMO)

    def real_figure():
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.1), sharey=True)
        for ax, (label, _) in zip(axes, REAL):
            key = 'real_' + label.split(' (')[0].replace(' ', '').lower()
            kplot.recovery_map(maps[key], title=label, ax=ax,
                               colorbar=(label == REAL[-1][0]))
        for ax in axes[1:]:
            ax.set_ylabel('')
        fig.tight_layout()
        return fig

    save_both(real_figure, 'real_maps', DEMO)
    log(f'Everything in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

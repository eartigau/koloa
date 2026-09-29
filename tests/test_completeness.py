#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Recovery maps, and fits that survive many outliers.

Created on 2026-09-27

@author: artigau
"""
import os

import numpy as np

from koloa import log as klog
from koloa.completeness import fold_test, inject, recovery_map
from koloa.data import RVData
from koloa.fit import RVModel
from koloa.simulate import BORDERLINE, CLEAR, REALISTIC, simulate

klog.VERBOSE = False
#: a real NIRPS sampling: the times and error bars of 181 exposures in
#: 58 visits (its velocities are all zero)
TEMPLATE = os.path.join(os.path.dirname(__file__), '..', 'data',
                        'nirps_template.csv')


def _series(seed: int) -> RVData:
    """A new noise realisation on a real NIRPS sampling"""
    tpl = RVData.from_csv(TEMPLATE, name='NIRPS template')
    return simulate(template=tpl, visit_jitter=2.0, seed=seed)['data']


def test_fold_false_alarms_match_the_threshold():
    # at K = 0, the fold test at fap = 5 % fires about 5 % of the time
    rmap = recovery_map(_series, [2.0, 200.0], [5.0, 6.0], ninj=1,
                        nnull=400, fap=0.05, seed=3, quiet=True)
    assert 0.02 < rmap.false_alarm < 0.09


def test_strong_planets_are_recovered_and_weak_ones_are_not():
    data = _series(5)
    strong = inject(data, 12.3, 15.0, 0.7)
    weak = inject(data, 12.3, 0.05, 0.7)
    assert fold_test(strong, 12.3)['dlnl'] > 20
    assert fold_test(weak, 12.3)['dlnl'] < 4.6
    rmap = recovery_map(data, [3.0, 30.0], [0.2, 0.4, 20.0, 40.0], ninj=10,
                        nnull=0, seed=4, quiet=True)
    assert rmap.rate[0, 0] < 0.3 and rmap.rate[0, 2] == 1.0
    assert np.isfinite(rmap.k_at(0.5)[0])


def test_calibrated_thresholds_hold_the_false_alarms():
    # on white noise the calibrated thresholds stay near the chi^2 one, and
    #   the false alarms near fap in every bin
    rmap = recovery_map(_series, [2.0, 20.0, 200.0], [5.0, 6.0], ninj=1,
                        nnull=300, fap=0.05, seed=7, calibrate=True,
                        quiet=True)
    thr = np.array(rmap.extra['thresholds'])
    assert np.all(thr >= -np.log(0.05) - 1e-9) and np.all(thr < 6.5)
    assert rmap.false_alarm < 0.06


def test_outlier_recipes_are_labelled_and_disjoint():
    tpl = RVData.from_csv(TEMPLATE, name='NIRPS template')
    sim = simulate(template=tpl, outliers=REALISTIC, seed=2)
    labels = sim['outlier_label']
    assert set(np.unique(labels)) == {'', 'clear', 'borderline'}
    assert len(CLEAR) == len(BORDERLINE) == 2
    size = np.abs(sim['outlier_offset']) / np.median(tpl.err)
    assert np.max(size[labels == 'borderline']) < 4.5
    # the two visit recipes never pick the same visit (single exposures of
    #   both kinds may share one)
    visits = simulate(template=tpl, outliers=[CLEAR[0], BORDERLINE[0]],
                      seed=2)
    seq, lab = visits['data'].seq, visits['outlier_label']
    for visit in np.unique(seq[lab != '']):
        assert len(set(lab[seq == visit])) == 1


def test_trimmed_start_survives_half_outliers():
    # half the points are outliers: a start from 'no outliers' stays there;
    #   the trimmed start of polyband does not
    rng = np.random.default_rng(2)
    npts = 200
    xval = np.sort(rng.uniform(0, 10, npts))
    coef = np.array([0.05, -0.8, 3.0, 2.0])
    yval = np.polyval(coef, xval) + rng.normal(size=npts)
    idx = rng.choice(npts, npts // 2, replace=False)
    yval[idx] += rng.choice([-1, 1], len(idx)) * rng.uniform(3, 30, len(idx))
    data = RVData(time=xval, rv=yval, err=np.ones(npts), seq=np.arange(npts),
                  name='cubic')
    fit = RVModel(data, [], trend=3, seq_jitter=False, likelihood='mixture',
                  unit='point').fit(nstart=1, quiet=True)
    flagged = fit.outlier_prob > 0.5
    mask = np.zeros(npts, dtype=bool)
    mask[idx] = True
    assert np.mean(flagged[mask]) > 0.9
    assert np.mean(flagged[~mask]) < 0.1
    curve = fit.model.systematics(fit.theta) + data.zero_point['inst']
    assert np.sqrt(np.mean((curve - np.polyval(coef, xval)) ** 2)) < 0.6


# =============================================================================
# End of code
# =============================================================================

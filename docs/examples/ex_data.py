#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read a series.

RVData reads an LBL .rdb or a csv, groups the exposures into visits from
the gaps in time, and keeps the activity indicators next to the velocities.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's series (RVData) and log (timestamped lines, 'value'
#   for a number)
from pathlib import Path

import numpy as np

from koloa.data import RVData
from koloa.log import log

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# the public HARPS-N velocities of Kepler-21, from DACE (the lines that
#   start with # are skipped)
# rjd, vrad and svrad are found by name; every other column with an error
#   column next to it becomes an indicator (fwhm with sfwhm, smw with
#   ssmw, ...)
# the exposures fall into visits: a gap of more than 0.3 day (the
#   argument sequence_gap) starts a new one; times in days, velocities and
#   errors in m/s
data = RVData.from_csv(ROOT / 'data' / 'kepler21_harpsn_dace_drs3.3.12.csv',
                       name='Kepler-21')
# npoints, nseq (the visits), baseline [days], median_err, rms and
#   robust_rms (1.4826 times the median absolute deviation) [m/s]
summ = data.summary()
log(f'{data.name}: {summ["npoints"]} exposures in {summ["nseq"]} visits '
    f'over {summ["baseline"]:.0f} days', 'value')
log(f'median error {summ["median_err"]:.2f} m/s, rms {summ["rms"]:.2f} m/s, '
    f'robust rms {summ["robust_rms"]:.2f} m/s', 'value')
# the median velocity of each instrument is taken out and kept in
#   zero_point; a file without an instrument column holds one instrument,
#   named 'inst'
log(f'zero point taken out: {data.zero_point["inst"]:.1f} m/s', 'value')
# data.seq: the visit of each exposure (0, 1, ... in time order)
sizes = np.bincount(data.seq)
log(f'exposures per visit: {sizes.min()} to {sizes.max()}, median '
    f'{np.median(sizes):.0f}', 'value')
# the names of the indicators; each holds its values and their errors
log(f'indicators: {", ".join(data.indicators)}')

# one point per visit: the weighted mean of its exposures
# (with the formal error of that mean)
visits = data.binned()
log(f'binned: {visits.n} points, median error {np.median(visits.err):.2f} '
    f'm/s', 'value')

# the first 250 days only (the visits are renumbered)
# select takes a boolean mask or indices and returns a copy
first = data.select(data.time < data.time[0] + 250)
log(f'first 250 days: {first.n} exposures in {first.nseq} visits', 'value')

# an indicator as a series of its own, ready for any koloa tool
# (its median taken out, the points without a value or an error dropped);
#   fwhm: the width of the cross-correlation function [m/s]
fwhm = data.indicator('fwhm')
log(f'{fwhm.name}: {fwhm.n} points, robust rms '
    f'{fwhm.summary()["robust_rms"]:.2f} m/s', 'value')

# the same times with other values: the residual of a straight line
# tref: the reference time of trends and phases, the mean time [days];
#   numpy's polyfit takes its weights as 1/sigma
tt = data.time - data.tref
coeffs = np.polyfit(tt, data.rv, 1, w=1 / data.err)
# with_values: a copy with the same times, visits and errors
resid = data.with_values(data.rv - np.polyval(coeffs, tt))
# the slope is in m/s per day
log(f'trend {365.25 * coeffs[0]:+.2f} m/s per year, residual robust rms '
    f'{resid.summary()["robust_rms"]:.2f} m/s', 'value')

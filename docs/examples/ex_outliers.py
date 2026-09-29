#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Why an outlier is an outlier.

explain_outliers compares every outlier of a fit with the good data, on
everything that was recorded with it: here the columns DACE gives for
HARPS (S/N, exposure time, the drift of the simultaneous reference, the
shape of the CCF) and the activity indicators. A single exposure is
compared with the good exposures of its visit, or, alone in its visit, with
the good visits nearest in time; a key is significantly off beyond the
Bonferroni threshold for the number of keys. The public HARPS velocities of
HD 69830 from 2012 on, its three planets fitted with the mixture.

Figure: each outlier (a row) against each key (a column), coloured by its
deviation from the good data (teal low, orange high); framed cells are
significantly off.

Created on 2026-09-29

@author: artigau
"""
# koloa's figures, public DACE data, the model, the log (timestamped lines,
#   'value' for a number) and the outliers
from pathlib import Path

from koloa import plotting as kplot
from koloa.dace import rvdata
from koloa.fit import RVModel
from koloa.log import log
from koloa.outliers import explain

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# the public velocities of HD 69830 from DACE, one instrument era (HARPS
#   before its 2015 upgrade), with every column of DACE kept in data.meta
data = rvdata(ROOT / 'data' / 'dace' / 'HD69830_dace.csv', name='HD 69830',
              instruments=['HARPS03'])
# the seasons from 2012 on (rjd = BJD - 2400000), to keep the example short
data = data.select(data.time > 56000)
log(f'{data.n} exposures in {data.nseq} visits; columns kept: '
    f'{len(data.meta)}, e.g. {", ".join(sorted(data.meta)[:4])}', 'info')

# the three known planets, each period free within 1 %, with the outlier
#   mixture (unit='both': a single exposure or a whole visit may be bad)
#   and a straight trend; nstart=1: one start per period
planets = [dict(period=per, period_range=(0.99 * per, 1.01 * per))
           for per in (8.669, 31.62, 201.0)]
fit = RVModel(data, planets, likelihood='mixture', unit='both',
              trend=1).fit(nstart=1, quiet=True)

# every outlier against the good data; keys='auto' recognises HARPS from
#   DACE and takes its list (and the indicators); keys=[...] takes any
#   columns by name
report = explain(fit)
for line in report.text().split('\n'):
    log(line, 'value')
# the first outlier: its keys, most deviant first, each with its value, the
#   reference it is compared with, its z and whether it is off
for key in report.units[0]['keys'][:3]:
    log(f'{key["key"]}: {key["value"]:.4g} against {key["reference"]:.4g}, '
        f'z = {key["z"]:+.1f}, significant: {key["significant"]}', 'value')

# koloa's web style; the figure: outliers against keys
kplot.set_style('web')
fig = kplot.outlier_keys(report)
path = kplot.savefig(fig, str(ROOT / 'docs/figures/examples/outliers.svg'))
log(f'figure: {Path(path).relative_to(ROOT)}')

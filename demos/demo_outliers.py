#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 19: why the outliers of real series are outliers.

For each star, the known planets are fitted with koloa's mixture (every
exposure and every visit may be an outlier), and koloa.outliers.explain
compares each outlier with the good data on everything that was recorded
with it: the header keywords that LBL copies into its .rdb files (the
SPIRou and NIRPS lists: S/N, airmass, telluric absorption, the shape of the
image, the Fabry-Perot, the age of the wavelength solution...), the
columns of DACE for HARPS and ESPRESSO, the activity indicators and the
error bar of the velocity. A single exposure is compared with the good
exposures of its visit, a whole visit with the good visits nearest in time.

The stars: GJ 687 (SPIRou) and Proxima (NIRPS), series of the SPIRou and
NIRPS teams that are not in this repository (a star whose file is missing
is skipped), and HD 69830 (HARPS, public, from DACE).

    python demo_outliers.py [--only gl687]

Created on 2026-09-29

@author: artigau
"""
import argparse
import os
import time
import warnings


from _common import ROOT, outdir, save_both, save_summary
from koloa import plotting as kplot
from koloa.dace import rvdata
from koloa.data import RVData
from koloa.fit import RVModel
from koloa.log import log
from koloa.outliers import explain
from koloa.utils import blas_threads

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'outliers'
#: the stars: the series, how to read it, and the known planets (period,
#: and the fraction of it they may move)
TARGETS = [
    dict(key='gl687', name='GJ 687', inst='SPIRou',
         file='gl687_spirou_pca2d_0-7.rdb', planets=[38.142]),
    dict(key='proxima', name='Proxima', inst='NIRPS',
         file='proxima_nirps_pca2d_0-7.rdb', planets=[11.18465]),
    dict(key='hd69830', name='HD 69830', inst='HARPS', dace=True,
         file=os.path.join('dace', 'HD69830_dace.csv'),
         instruments=['HARPS03'], planets=[8.669, 31.62, 201.0])]
PERIOD_FREEDOM = 0.01


# =============================================================================
# Define functions
# =============================================================================
def load(target: dict):
    """the series of a target, or None when its file is not here"""
    path = os.path.join(ROOT, 'data', target['file'])
    if not os.path.exists(path):
        log(f'{target["name"]}: {path} is not here (a series of the '
            f'{target["inst"]} team), skipped', 'warn')
        return None
    if target.get('dace'):
        return rvdata(path, name=target['name'],
                      instruments=target.get('instruments'))
    return RVData.from_csv(path, name=target['name'])


def study(target: dict) -> dict:
    """the fit, the outliers and why, and the figures of a target"""
    data = load(target)
    if data is None:
        return None
    start = time.time()
    planets = [dict(period=per, period_range=((1 - PERIOD_FREEDOM) * per,
                                              (1 + PERIOD_FREEDOM) * per))
               for per in target['planets']]
    with blas_threads(1):
        fit = RVModel(data, planets, likelihood='mixture', unit='both',
                      trend=1).fit(nstart=2, quiet=True)
    report = explain(fit)
    for line in report.text().split('\n'):
        log(line, 'value')
    key = target['key']
    with open(os.path.join(outdir(DEMO), f'{key}.txt'), 'w') as handle:
        handle.write(report.text() + '\n')
    save_both(lambda: kplot.outlier_keys(report), f'{key}_keys', DEMO)
    # the series itself, for the paper only (the site shows the keys)
    kplot.set_style('paper')
    kplot.savefig(kplot.timeseries(
        data, fit.outlier_prob, fit=fit,
        title=f'{target["name"]} ({target["inst"]}): every exposure coloured '
              f'by its outlier probability'),
        os.path.join(outdir(DEMO), f'{key}_timeseries.pdf'))
    return dict(name=target['name'], inst=target['inst'], n=int(data.n),
                nvisits=int(data.nseq), planets=target['planets'],
                orbits=[{name: [float(val) for val in orb[name]]
                         for name in ('P', 'K')} for orb in fit.orbits()],
                report=report.to_dict(), runtime=time.time() - start)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('--only', default=None,
                        help='one target (gl687, proxima or hd69830)')
    args = parser.parse_args()
    warnings.simplefilter('ignore')
    start = time.time()
    log('Demo: why the outliers of real series are outliers')
    path = os.path.join(outdir(DEMO), 'summary.json')
    summary = {}
    if args.only and os.path.exists(path):
        import json
        with open(path) as handle:
            summary = json.load(handle)
    for target in TARGETS:
        if args.only and target['key'] != args.only:
            continue
        log(f'{target["name"]} ({target["inst"]})')
        result = study(target)
        if result is not None:
            summary[target['key']] = result
    summary['runtime_minutes'] = (time.time() - start) / 60
    save_summary(DEMO, summary)
    log(f'done in {summary["runtime_minutes"]:.1f} min; everything in '
        f'{outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

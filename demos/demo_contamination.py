#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 11: how many outliers can a fit take? A Student-t against the mixture.

A Student-t likelihood (the robust option of polyband, Artigau 2026, and
koloa's likelihood='student') gives every point heavy tails; koloa's mixture
says that each point, or each visit, is either good or an outlier. Both are
scale mixtures of gaussians, but they do not weigh an outlier the same: the
t keeps a little of its pull, the mixture none once it is flagged. Both need
a start that the outliers cannot drag, which koloa borrows from polyband:
the fit of the least deviant fraction of the units.

Three tests, with outliers of both signs from 3 to 30 sigma (the borderline
ones included), as their fraction grows:
1. a cubic polynomial, 200 points: least squares, polyband (nu = 4), koloa's
   Student-t (nu = 4) and koloa's mixture;
2. an eccentric orbit (K = 10 m/s) on a real NIRPS sampling, with single
   exposures as outliers;
3. the same orbit with whole visits as outliers and a 2 m/s visit jitter,
   which a Student-t on points cannot describe.

    python demo_contamination.py --workers 8

Created on 2026-09-27

@author: artigau
"""
import argparse
import json
import os
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from _common import DATA, outdir, save_both
from koloa import plotting as kplot
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'contamination'
POLY_FRACS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
ORBIT_FRACS = [0.0, 0.1, 0.2, 0.3, 0.4]
NREAL_POLY, NREAL_ORBIT = 24, 16
#: the cubic, and the noise
COEF = np.array([0.05, -0.8, 3.0, 2.0])
SIGMA = 1.0
#: the orbit (as in demo_posterior.py)
ORBIT = dict(P=11.2, K=10.0, e=0.3, omega=np.radians(60))
POLY_METHODS = [('lsq', 'least squares'), ('polyband', 'polyband (t, nu = 4)'),
                ('koloa_student', 'koloa Student-t (nu = 4)'),
                ('koloa_mixture', 'koloa mixture')]
ORBIT_METHODS = [('gaussian', 'Gaussian'),
                 ('koloa_student', 'koloa Student-t (nu = 4)'),
                 ('koloa_mixture', 'koloa mixture')]


# =============================================================================
# Define functions
# =============================================================================
def poly_case(job):
    """One cubic, one contamination, every method"""
    frac, seed = job
    warnings.simplefilter('ignore')
    from koloa import log as klog
    from koloa.data import RVData
    from koloa.fit import RVModel
    klog.VERBOSE = False
    rng = np.random.default_rng(seed)
    npts = 200
    xval = np.sort(rng.uniform(0, 10, npts))
    truth = np.polyval(COEF, xval)
    yval = truth + SIGMA * rng.normal(size=npts)
    nout = int(round(frac * npts))
    idx = rng.choice(npts, nout, replace=False)
    yval[idx] += rng.choice([-1, 1], nout) * SIGMA * rng.uniform(3, 30, nout)
    grid = np.linspace(0, 10, 400)
    tgrid = np.polyval(COEF, grid)
    out = dict(frac=frac, seed=seed)
    out['lsq'] = float(np.sqrt(np.mean(
        (np.polyval(np.polyfit(xval, yval, 3), grid) - tgrid) ** 2)))
    try:
        from polyband import fit_polyband
        band = fit_polyband(xval, yval, order_mean=3, order_width=0, nu=4)
        out['polyband'] = float(np.sqrt(np.mean((band.predict(grid)
                                                 - tgrid) ** 2)))
    except ImportError:
        out['polyband'] = float('nan')
    # every point is its own visit; the curve is the trend of the model
    data = RVData(time=xval, rv=yval, err=np.full(npts, SIGMA),
                  seq=np.arange(npts), name='cubic')
    gdata = RVData(time=grid, rv=np.zeros(len(grid)), err=np.ones(len(grid)),
                   seq=np.arange(len(grid)), name='grid')
    mask = np.zeros(npts, dtype=bool)
    mask[idx] = True
    for name, kwargs in (('koloa_student', dict(likelihood='student',
                                                dof=4.0)),
                         ('koloa_mixture', dict(likelihood='mixture',
                                                unit='point'))):
        model = RVModel(data, [], trend=3, seq_jitter=False, **kwargs)
        fit = model.fit(nstart=1, quiet=True)
        gmodel = RVModel(gdata, [], trend=3, seq_jitter=False, **kwargs)
        gmodel.tref = model.tref
        gmodel.tnorm = (grid - model.tref) / max(data.baseline, 1e-9)
        # the series is fitted relative to its zero point (its median)
        curve = gmodel.systematics(fit.theta) + data.zero_point['inst']
        out[name] = float(np.sqrt(np.mean((curve - tgrid) ** 2)))
        if name == 'koloa_mixture':
            flagged = fit.outlier_prob > 0.5
            out['recall'] = float(np.mean(flagged[mask])) if nout else None
            out['false'] = float(np.mean(flagged[~mask]))
    return out


def orbit_case(job):
    """One orbit, one contamination, every method"""
    frac, seed, kind = job
    warnings.simplefilter('ignore')
    from koloa import log as klog
    from koloa.data import RVData
    from koloa.fit import RVModel
    from koloa.simulate import simulate
    klog.VERBOSE = False
    tpl = RVData.from_csv(DATA)
    rng = np.random.default_rng(5000 + seed)
    tperi = float(tpl.time[0] + rng.uniform(0, ORBIT['P']))
    outliers = []
    if frac > 0:
        outliers = [dict(kind='spike' if kind == 'point' else 'visit',
                         frac=frac, amplitude=(3.0, 30.0))]
    sim = simulate(planets=[dict(ORBIT, tp=tperi)], template=tpl,
                   outliers=outliers,
                   visit_jitter=0.0 if kind == 'point' else 2.0, seed=seed)
    data = sim['data']
    out = dict(frac=frac, seed=seed, kind=kind)
    visits = kind == 'visit'
    for name, kwargs in (
            ('gaussian', dict(likelihood='gaussian', seq_jitter=visits)),
            ('koloa_student', dict(likelihood='student', dof=4.0,
                                   seq_jitter=False)),
            ('koloa_mixture', dict(likelihood='mixture',
                                   unit='both' if visits else 'point',
                                   seq_jitter=visits))):
        fit = RVModel(data, [dict(period=ORBIT['P'], eccentric=True)],
                      **kwargs).fit(nstart=4, quiet=True)
        orbit = fit.orbits()[0]
        out[name] = dict(K=orbit['K'][0], e=orbit['e'][0], P=orbit['P'][0])
    return out


def summarise(values: np.ndarray) -> list:
    """Median and 16th and 84th percentiles"""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return [float('nan')] * 3
    return [float(v) for v in np.percentile(values, [50, 16, 84])]


def figure(summary: dict):
    """The error of each method as the outliers grow"""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
    colors = dict(lsq=kplot.C['neutral'], polyband=kplot.C['hard'],
                  koloa_student=kplot.C['white'], koloa_mixture=kplot.C['koloa'],
                  gaussian=kplot.C['gaussian'])
    panels = [('poly', POLY_METHODS, 'curve error [$\\sigma$]',
               'a cubic, 200 points'),
              ('point', ORBIT_METHODS, '$|K - K_\\mathrm{true}|$ [m s$^{-1}$]',
               'an orbit, bad exposures'),
              ('visit', ORBIT_METHODS, '$|K - K_\\mathrm{true}|$ [m s$^{-1}$]',
               'an orbit, bad visits')]
    for ax, (key, methods, ylabel, title) in zip(axes, panels):
        if key not in summary:
            ax.set_visible(False)
            continue
        fracs = np.array(summary[key]['fracs'])
        for it, (name, label) in enumerate(methods):
            stats = np.array(summary[key][name])
            if not np.any(np.isfinite(stats)):
                continue
            shift = (it - 1.5) * 0.006
            ax.errorbar(100 * (fracs + shift), stats[:, 0],
                        [stats[:, 0] - stats[:, 1], stats[:, 2] - stats[:, 0]],
                        fmt='-o', color=colors[name], ms=3.5, lw=1.3,
                        elinewidth=0.9, label=label)
        ax.set_yscale('log')
        ax.set_xlabel('outliers [%]')
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc='left')
    # one legend under the three panels, where it hides no curve
    handles, labels = [], []
    for ax in axes:
        for handle, label in zip(*ax.get_legend_handles_labels()):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.legend(handles, labels, loc='lower center', ncol=len(labels),
               fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))
    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--figures-only', action='store_true')
    args = parser.parse_args()
    warnings.simplefilter('ignore')
    path = os.path.join(outdir(DEMO), 'summary.json')
    log('Demo: a Student-t against the mixture, as the outliers grow')
    if not args.figures_only or not os.path.exists(path):
        summary = {}
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            poly = list(pool.map(poly_case, [(f, s) for f in POLY_FRACS
                                             for s in range(NREAL_POLY)]))
            summary['poly'] = dict(fracs=POLY_FRACS, runs=poly)
            for name, _ in POLY_METHODS:
                summary['poly'][name] = [summarise(
                    [r[name] for r in poly if r['frac'] == f])
                    for f in POLY_FRACS]
            for kind in ('point', 'visit'):
                runs = list(pool.map(orbit_case, [
                    (f, s, kind) for f in ORBIT_FRACS
                    for s in range(NREAL_ORBIT)]))
                summary[kind] = dict(fracs=ORBIT_FRACS, runs=runs)
                for name, _ in ORBIT_METHODS:
                    summary[kind][name] = [summarise(
                        [abs(r[name]['K'] - ORBIT['K']) for r in runs
                         if r['frac'] == f]) for f in ORBIT_FRACS]
                    summary[kind][name + '_e'] = [summarise(
                        [abs(r[name]['e'] - ORBIT['e']) for r in runs
                         if r['frac'] == f]) for f in ORBIT_FRACS]
        with open(path, 'w') as handle:
            json.dump(summary, handle, indent=1, default=float)
    with open(path) as handle:
        summary = json.load(handle)
    for key, methods in (('poly', POLY_METHODS), ('point', ORBIT_METHODS),
                         ('visit', ORBIT_METHODS)):
        for it, frac in enumerate(summary[key]['fracs']):
            text = ', '.join(f'{name} {summary[key][name][it][0]:.3f}'
                             for name, _ in methods)
            log(f'  {key:5s} {100 * frac:3.0f} % outliers: {text}', 'value')
    save_both(lambda: figure(summary), 'errors', DEMO)
    log(f'Everything in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

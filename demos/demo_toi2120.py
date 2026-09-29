#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 5: TOI-2120 b under more and more bad visits.

The real SPIRou series of TOI-2120 (316 exposures in 81 visits, most of
them polarimetric sequences of four exposures, two observing seasons with
the 188-day gap of the target's yearly visibility), cleaned by PCA2D
(nominal of 2026-09-26: 0 stellar and 7 observer components, bias anchored
in the observer frame at 10 sigma, cut at 3 sigma) and measured by LBL.

TOI-2120 b transits (Hori et al. 2024): its period and time of conjunction
come from TESS, so the velocities only measure K. The activity and the
slow wander are a Matern-3/2 GP with a correlation length of about 50 days
(a gaussian prior on log length, centred on 50 d, 0.3 wide).

Bad visits are then added to the REAL series, never moving its sampling:
visits moved as a whole, up or down. Two scans:
- the AMPLITUDE of the bad visits grows (10 % of the visits, moved by A to
  1.5 A, A from 0 to 100 m/s; the visits scatter by about 12 m/s on their
  own);
- the FRACTION of bad visits grows (5 to 30 %, moved by 60 to 90 m/s).
At every level the orbit is fitted, always jointly with the GP, with:
- gaussian noise;
- gaussian noise after a 3-sigma hard clip of the exposures;
- koloa's mixture, a visit being an outlier or not;
- and, for reference, gaussian noise with the bad visits removed by hand
  (the oracle, which no real analysis has).

What it shows: bad visits a few times the natural scatter bias the
gaussian K and blow up its error bar, and koloa removes them. Bad visits
barely larger than the scatter cannot be told from it by any noise model;
koloa then falls back on the gaussian answer with a larger jitter, which is
the right answer to the question it is asked, and never a worse one.

    python demo_toi2120.py [--workers 4] [--nreal 6]

Created on 2026-09-27

@author: artigau
"""
import argparse
import os
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from _common import HERE, ROOT, outdir, save_both

DEMO = 'toi2120'
RDB = os.path.join(ROOT, 'data', 'toi2120_spirou_pca2d_0-7_bias10s3.rdb')
#: TOI-2120 b, Hori et al. (2024), from the NASA Exoplanet Archive
PLANET = dict(period=5.7998164, period_err=0.0000035, tc=58795.82368,
              tc_err=0.00041)
GP = dict(kernel='matern32', prior={'log_length': (np.log(50.0), 0.3)},
          init={'log_length': np.log(50.0)})
AMPLITUDES = [0.0, 15.0, 30.0, 60.0, 100.0]
FRACTIONS = [0.05, 0.10, 0.20, 0.30]
AMP_FRACTION = 0.10
FRAC_AMPLITUDE = 60.0


def load():
    """The TOI-2120 series"""
    from koloa.data import RVData
    return RVData.from_csv(RDB, name='TOI-2120 (SPIRou, PCA2D 0-7)')


def corrupt(data, frac, amp, seed):
    """
    Move a fraction of the visits as a whole, up or down, by amp to 1.5 amp

    :return: tuple, the corrupted series and the mask of moved exposures
    """
    rng = np.random.default_rng(seed)
    nbad = int(round(frac * data.nseq)) if amp > 0 else 0
    bad = rng.choice(data.nseq, nbad, replace=False)
    rv = data.rv.copy()
    mask = np.zeros(data.n, dtype=bool)
    for visit in bad:
        sel = data.seq == visit
        rv[sel] += rng.choice([-1, 1]) * amp * rng.uniform(1.0, 1.5)
        mask |= sel
    return data.with_values(rv), mask


def fits(args):
    """Every method at one level and one realisation"""
    scan, frac, amp, seed = args
    from koloa import log as klog
    from koloa.fit import RVModel
    klog.VERBOSE = False
    warnings.simplefilter('ignore')
    data, mask = corrupt(load(), frac, amp, seed)
    out = dict(scan=scan, frac=frac, amp=amp,
               level=frac if scan == 'fraction' else amp, seed=seed,
               nbad=int(np.sum(mask)))
    oracle = RVModel(data.select(~mask), [PLANET], gp=GP,
                     likelihood='gaussian').fit(nstart=1, quiet=True)
    out['oracle'] = oracle.orbits()[0]['K']
    gauss = RVModel(data, [PLANET], gp=GP, likelihood='gaussian').fit(
        nstart=1, quiet=True)
    out['gaussian'] = gauss.orbits()[0]['K']
    # the classical cut: exposures more than 3 robust sigma from the
    #   gaussian model, removed, and the fit done again
    resid = data.rv - gauss.model.mean_model(gauss.theta) - \
        gauss.gp_prediction()[0]
    zval = np.abs(resid - np.median(resid)) / (
        1.4826 * np.median(np.abs(resid - np.median(resid))))
    keep = zval < 3
    clip = RVModel(data.select(keep), [PLANET], gp=GP,
                   likelihood='gaussian').fit(nstart=1, quiet=True)
    out['hard'] = clip.orbits()[0]['K']
    out['nclipped'] = int(np.sum(~keep))
    koloa = RVModel(data, [PLANET], gp=GP, likelihood='mixture',
                    unit='sequence').fit(nstart=1, quiet=True)
    out['koloa'] = koloa.orbits()[0]['K']
    flagged = koloa.outlier_prob > 0.5
    out['caught'] = int(np.sum(flagged & mask))
    out['false'] = int(np.sum(flagged & ~mask))
    return out


def figures(results):
    """K against the outlier level, and the series at a bad level"""
    import matplotlib.pyplot as plt
    from koloa import plotting as kplot
    from koloa.fit import RVModel
    methods = [('gaussian', 'Gaussian + GP', 'gaussian'),
               ('hard', '3-sigma clip + GP', 'hard'),
               ('koloa', 'koloa + GP', 'koloa'),
               ('oracle', 'bad visits removed by hand', 'muted')]

    def k_vs_level():
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
        clean = [res['koloa'][0] for res in results
                 if res['scan'] == 'amplitude' and res['amp'] == 0]
        for ax, scan, levels, xlabel in (
                (axes[0], 'amplitude', AMPLITUDES,
                 f'bad-visit offset A [m/s] ({int(100 * AMP_FRACTION)}% '
                 f'of visits)'),
                (axes[1], 'fraction', FRACTIONS,
                 f'fraction of bad visits (offsets {int(FRAC_AMPLITUDE)} to '
                 f'{int(1.5 * FRAC_AMPLITUDE)} m/s)')):
            ax.axhline(np.mean(clean), color=kplot.C['muted'], lw=0.8,
                       ls=(0, (3, 2)))
            span = (max(levels) - min(levels)) or 1.0
            for im, (key, label, ckey) in enumerate(methods):
                for lev in levels:
                    sub = [res for res in results if res['scan'] == scan
                           and np.isclose(res['level'], lev)]
                    vals = np.array([res[key][0] for res in sub])
                    if len(vals) == 0:
                        continue
                    xpos = lev + (im - 1.5) * 0.028 * span
                    # every realisation as a small dot, the median and the
                    #   16-84 % range of the realisations as the error bar
                    ax.plot(np.full(len(vals), xpos), vals, '.', ms=3,
                            color=kplot.C[ckey], alpha=0.5)
                    low, med, high = np.percentile(vals, [16, 50, 84])
                    ax.errorbar(xpos, med, [[med - low], [high - med]],
                                fmt='D' if key == 'oracle' else 'o', ms=5,
                                color=kplot.C[ckey], elinewidth=1.4,
                                label=label if (lev == levels[0] and
                                                scan == 'amplitude') else
                                None)
            ax.set_xlabel(xlabel)
            ax.set_xticks(levels)
        axes[0].set_ylabel('K of TOI-2120 b [m s$^{-1}$]')
        axes[0].set_title('TOI-2120 b, more and more bad visits', loc='left')
        fig.legend(loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.04),
                   frameon=False)
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        return fig

    save_both(k_vs_level, 'K_vs_level', DEMO)
    # the series with 10 % of the visits moved by 60 to 90 m/s
    data, mask = corrupt(load(), 0.10, 60.0, 100)
    fit = RVModel(data, [PLANET], gp=GP, likelihood='mixture',
                  unit='sequence').fit(nstart=1, quiet=True)
    save_both(lambda: kplot.timeseries(
        data, fit.outlier_prob, fit=fit,
        title='TOI-2120, 10% of the visits moved by 60-90 m/s: GP + orbit, '
              'koloa'), 'timeseries_bad', DEMO)
    save_both(lambda: kplot.phase(
        fit, title='TOI-2120 b with bad visits, koloa: K = {:.2f} m/s'.format(
            fit.orbits()[0]['K'][0])), 'phase_bad', DEMO)
    clean_fit = RVModel(load(), [PLANET], gp=GP, likelihood='mixture',
                        unit='sequence').fit(nstart=1, quiet=True)
    gauss_fit = RVModel(load(), [PLANET], gp=GP,
                        likelihood='gaussian').fit(nstart=1, quiet=True)
    describe_clean(clean_fit, gauss_fit)
    save_both(lambda: kplot.timeseries(
        clean_fit.model.data, clean_fit.outlier_prob, fit=clean_fit,
        title='TOI-2120 as observed: two seasons, GP (Matern-3/2) + orbit'),
        'timeseries_clean', DEMO)
    save_both(lambda: kplot.phase(
        clean_fit, title='TOI-2120 b as observed: K = {:.2f} m/s'.format(
            clean_fit.orbits()[0]['K'][0])), 'phase_clean', DEMO)


def describe_clean(koloa_fit, gauss_fit):
    """
    What koloa finds in the series as observed: the visits it flags, and
    what the GP and the jitters become without them

    Written to clean_series.json, for the paper and the web page.
    """
    import json
    from koloa.log import log
    data = koloa_fit.model.data
    resid = koloa_fit.residuals()
    visits = []
    for visit in range(data.nseq):
        sel = data.seq == visit
        prob = float(koloa_fit.outlier_prob[sel].max())
        if prob < 0.5:
            continue
        wgt = 1 / data.err[sel] ** 2
        visits.append(dict(time=float(np.mean(data.time[sel])),
                           nexp=int(np.sum(sel)), prob=prob,
                           offset=float(np.sum(wgt * resid[sel])
                                        / np.sum(wgt))))
    out = dict(visits=visits)
    for name, fit in (('koloa', koloa_fit), ('gaussian', gauss_fit)):
        model = fit.model
        orbit = fit.orbits()[0]
        out[name] = dict(
            K=orbit['K'], gp_sigma=float(np.exp(fit.theta[
                model.index['gp_log_sigma']])),
            gp_length=float(np.exp(fit.theta[model.index['gp_log_length']])),
            seq_jitter=float(np.exp(fit.theta[model.index['log_sjit']])),
            jitter=float(np.exp(fit.theta[model.index['log_jit_inst']])))
    with open(os.path.join(outdir(DEMO), 'clean_series.json'), 'w') as handle:
        json.dump(out, handle, indent=1)
    for vis in visits:
        log(f'  flagged visit at {vis["time"]:.2f}: {vis["offset"]:+.1f} m/s, '
            f'P(outlier) = {vis["prob"]:.2f}', 'value')
    for name in ('gaussian', 'koloa'):
        val = out[name]
        log(f'  {name:9s} K = {val["K"][0]:.2f} m/s, GP sigma = '
            f'{val["gp_sigma"]:.1f} m/s, length {val["gp_length"]:.0f} d, '
            f'visit jitter {val["seq_jitter"]:.1f} m/s', 'value')


def main():
    from koloa.log import log
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--nreal', type=int, default=6)
    parser.add_argument('--figures-only', action='store_true',
                        help='redraw from results.npy without refitting')
    args = parser.parse_args()
    if args.figures_only:
        results = list(np.load(os.path.join(outdir(DEMO), 'results.npy'),
                               allow_pickle=True))
        figures(results)
        log(f'Figures in {outdir(DEMO)}')
        return
    log('Demo: TOI-2120 b under more and more bad visits')
    jobs = [('amplitude', AMP_FRACTION, amp, 100 + ir) for amp in AMPLITUDES
            for ir in range(args.nreal if amp > 0 else 1)]
    jobs += [('fraction', frac, FRAC_AMPLITUDE, 200 + ir)
             for frac in FRACTIONS for ir in range(args.nreal)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(fits, jobs))
    lines = ['scan       level  nbad  oracle K     gaussian K   clip K (n)'
             '        koloa K   (caught/false)']
    for res in sorted(results, key=lambda rr: (rr['scan'], rr['level'],
                                                  rr['seed'])):
        lines.append(
            f'{res["scan"]:9s} {res["level"]:6.2f} {res["nbad"]:5d}  '
            f'{res["oracle"][0]:5.2f}+-{res["oracle"][1]:4.2f}  '
            f'{res["gaussian"][0]:5.2f}+-{res["gaussian"][1]:4.2f}  '
            f'{res["hard"][0]:5.2f}+-{res["hard"][1]:4.2f} ({res["nclipped"]:3d})'
            f'  {res["koloa"][0]:5.2f}+-{res["koloa"][1]:4.2f} '
            f'({res["caught"]}/{res["false"]})')
    for line in lines:
        log(line, 'value')
    with open(os.path.join(outdir(DEMO), 'K_vs_level.txt'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    np.save(os.path.join(outdir(DEMO), 'results.npy'), results,
            allow_pickle=True)
    figures(results)
    log(f'Figures in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

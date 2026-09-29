#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 12: GL 725 A and B with SPIRou, each pulled by the other.

GL 725 B (GJ 725 B, Struve 2398 B) is an M3.5 dwarf at 3.5 pc, 13 arcsec
from GL 725 A. Its SPIRou velocities (LBL, 959 exposures in 237 visits over
6.7 years) drift, and so do those of A (SPIRou, 2019 to 2022).

Their secular accelerations, the perspective effect of their proper motions
(about 0.4 m/s/yr each from Gaia DR3, koloa.secular), are already out of
these velocities: APERO computes the barycentric correction with
barycorrpy, whose formula (Wright & Eastman 2014) includes the proper
motion and parallax of the star. Adding koloa's secular term to the fit
would count them twice. What drifts is then the orbit of each star around
the other: a free linear trend, fitted with koloa's outlier-aware
likelihood and sampled. The two accelerations must have opposite signs, in
the inverse ratio of the masses. The demo then runs the whole koloa
analysis on GL 725 B with the drift taken out: FIP periodogram, flagged
visits, duck test.

    python demo_gl725b.py [--skip-analysis]

Created on 2026-09-27

@author: artigau
"""
import argparse
import json
import os
import warnings

import numpy as np

from _common import ROOT, outdir, save_both, save_summary
from koloa import plotting as kplot
from koloa.analyze import analyze
from koloa.data import RVData
from koloa.fit import RVModel
from koloa.log import log
from koloa.secular import perspective_from_astrometry, remove_perspective

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'gl725b'
FILES = {'GJ 725 B': os.path.join(ROOT, 'data',
                                  'gl725b_spirou_lbl_2026-09-24.rdb'),
         'GJ 725 A': os.path.join(ROOT, 'data',
                                  'gl725a_spirou_lbl_slinky05.rdb')}
GAIA = os.path.join(ROOT, 'data', 'gaia_dr3_astrometry.json')
#: APERO's barycentric correction (barycorrpy, with the Gaia astrometry of
#:   the star) has already removed the secular acceleration
SECULAR_REMOVED = True


# =============================================================================
# Define functions
# =============================================================================
def drift_fit(data: RVData, secular, trend: int = 1):
    """
    The drift of a star without planets: a free polynomial trend (its
    slope is the acceleration of the star, from the other star), and the
    perspective acceleration as a parameter with its Gaia prior when the
    velocities still hold it (secular None when they do not)

    :return: FitResult, and the orbital acceleration [m/s/yr] as draws
    """
    # whole visits as the outlier units: 30 times faster than 'both' on
    #   nearly a thousand exposures, and the same drift
    model = RVModel(data, [], trend=trend, perspective=secular,
                    likelihood='mixture', unit='sequence')
    res = model.sample(nsteps=6000, nburn=1500, converge=True, quiet=True)
    # the linear trend is per baseline; in m/s/yr at the reference time
    slope = res.chain[:, model.index['trend_1']] / data.baseline * 365.25
    return res, slope


def window_only():
    """The drift of B over the time span of A, and the ratio over it"""
    path = os.path.join(outdir(DEMO), 'summary.json')
    with open(path) as handle:
        summary = json.load(handle)
    data_a = RVData.from_csv(FILES['GJ 725 A'], name='GJ 725 A')
    data_b = RVData.from_csv(FILES['GJ 725 B'], name='GJ 725 B')
    tmin, tmax = data_a.time.min(), data_a.time.max()
    sub = data_b.select((data_b.time >= tmin) & (data_b.time <= tmax))
    res, slope = drift_fit(sub, None, trend=1)
    slope_a = drift_fit(data_a, None, trend=1)[1]
    draws = -np.random.default_rng(1).choice(slope_a, 20000) / \
        np.random.default_rng(2).choice(slope, 20000)
    summary['window'] = dict(
        tmin=float(tmin), tmax=float(tmax), npoints=sub.n, nvisits=sub.nseq,
        orbital_b=[float(v) for v in np.percentile(slope, [16, 50, 84])],
        mass_ratio=[float(v) for v in np.percentile(draws, [16, 50, 84])])
    log(f'  B over the span of A ({sub.nseq} visits): '
        f'{summary["window"]["orbital_b"][1]:.2f} m/s/yr; ratio '
        f'{summary["window"]["mass_ratio"][1]:.2f} '
        f'({summary["window"]["mass_ratio"][0]:.2f} to '
        f'{summary["window"]["mass_ratio"][2]:.2f})', 'value')
    save_summary(DEMO, summary)


def coherence_only():
    """The coherence check of the GL 725 B analysis, done again"""
    from koloa.diagnostics import coherence
    path = os.path.join(outdir(DEMO), 'summary.json')
    with open(path) as handle:
        summary = json.load(handle)
    data = RVData.from_csv(FILES['GJ 725 B'], name='GJ 725 B')
    flat = remove_perspective(data, summary['GJ 725 B']['linear']['total'][1])
    period = summary['analysis']['period']
    # the outlier probabilities of a circular fit at the period
    fit = RVModel(flat, [dict(period=period)], likelihood='mixture',
                  unit='sequence').fit(quiet=True)
    out = {}
    for split in ('halves', 'seasons'):
        coh = coherence(flat, period, prob=fit.outlier_prob, split=split)
        out[split] = dict(p_amplitude=coh['p_amplitude'],
                          p_vector=coh['p_vector'], jitter=coh['jitter'],
                          K=[ch['K'] for ch in coh['chunks']],
                          sK=[ch['sK'] for ch in coh['chunks']])
        log(f'  coherence by {split}: K = '
            + ', '.join(f'{k:.1f}+-{e:.1f}' for k, e in
                        zip(out[split]['K'], out[split]['sK']))
            + f' m/s, p(amplitude) {coh["p_amplitude"]:.2e}, p(vector) '
              f'{coh["p_vector"]:.2e}, jitter {coh["jitter"]:.2f} m/s',
            'value')
    summary['analysis']['coherence'] = out
    # the verdict again, with the same rule as koloa.diagnostics.duck_test
    #   (flag when the p of amplitude and phase, halves or seasons, is below
    #   0.01; one flag is inconclusive, two or more are not a planet), the
    #   checks that failed otherwise kept
    ana = summary['analysis']
    old = ana.get('verdict_before_coherence_fix', ana['verdict'])
    ana['verdict_before_coherence_fix'] = old
    if ': ' in old and not old.startswith('NOT DETECTED'):
        flags = [name for name in old.split(': ', 1)[1].split(', ')
                 if name != 'coherence']
        if min(out['halves']['p_vector'], out['seasons']['p_vector']) < 0.01:
            flags.append('coherence')
        if not flags:
            verdict = ('PLANET CANDIDATE: it looks, swims and quacks like a '
                       'planet')
        elif len(flags) == 1:
            verdict = ('INCONCLUSIVE: one test speaks against a planet ('
                       + flags[0] + ')')
        else:
            verdict = ('NOT A PLANET (activity or systematics more likely): '
                       + ', '.join(flags))
        ana['verdict'] = verdict
        log(f'  verdict with the coherence check done again: {verdict}',
            'value')
    save_summary(DEMO, summary)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-analysis', action='store_true',
                        help='keep the analysis of GL 725 B already in the '
                             'summary (the drift it removes is the same)')
    parser.add_argument('--window-only', action='store_true',
                        help='only fit the drift of B over the time span '
                             'of A, and the ratio there')
    parser.add_argument('--coherence-only', action='store_true',
                        help='only redo the coherence check of the '
                             'analysis (visit means, jitter measured)')
    args = parser.parse_args()
    if args.coherence_only:
        coherence_only()
        return
    if args.window_only:
        window_only()
        return
    warnings.simplefilter('ignore')
    log('Demo: GL 725 A and B with SPIRou, each pulled by the other')
    path_summary = os.path.join(outdir(DEMO), 'summary.json')
    previous = {}
    if os.path.exists(path_summary):
        with open(path_summary) as handle:
            previous = json.load(handle)
    with open(GAIA) as handle:
        gaia = json.load(handle)['stars']
    summary = {}
    fits = {}
    for star, path in FILES.items():
        if not os.path.exists(path):
            log(f'  {path} not found, {star} skipped', 'warn')
            continue
        data = RVData.from_csv(path, name=star)
        sec = perspective_from_astrometry(gaia[star])
        log(f'  {star}: {data.n} exposures in {data.nseq} visits over '
            f'{data.baseline:.0f} d; secular acceleration from Gaia DR3 '
            f'{sec[0]:.4f} +- {sec[1]:.5f} m/s/yr', 'value')
        entry = dict(npoints=data.n, nvisits=data.nseq,
                     baseline=data.baseline, secular=sec,
                     secular_removed=SECULAR_REMOVED)
        # a straight drift, and a curved one to check it is straight enough
        for trend in (1, 2):
            res, slope = drift_fit(data, None if SECULAR_REMOVED else sec,
                                   trend=trend)
            if SECULAR_REMOVED:
                secacc = np.zeros_like(slope)
            else:
                secacc = res.chain[:, res.model.index['secacc']]
            key = 'linear' if trend == 1 else 'quadratic'
            entry[key] = dict(
                orbital=[float(v) for v in np.percentile(slope,
                                                         [16, 50, 84])],
                secacc=[float(v) for v in np.percentile(secacc,
                                                        [16, 50, 84])],
                total=[float(v) for v in np.percentile(slope + secacc,
                                                       [16, 50, 84])],
                bic=res.bic(), tau_max=res.diagnostics['tau_max'],
                converged=res.diagnostics['converged'],
                nflag=int(np.sum(res.outlier_prob > 0.5)))
            if trend == 2:
                curv = res.chain[:, res.model.index['trend_2']]
                entry[key]['curvature'] = [float(v) for v in np.percentile(
                    curv, [16, 50, 84])]
            orb = entry[key]['orbital']
            log(f'    {key} drift: orbital acceleration {orb[1]:.3f} '
                f'(-{orb[1] - orb[0]:.3f} +{orb[2] - orb[1]:.3f}) m/s/yr, '
                f'total {entry[key]["total"][1]:.3f} m/s/yr, BIC '
                f'{entry[key]["bic"]:.1f}, {entry[key]["nflag"]} exposures '
                f'flagged', 'value')
            if trend == 1:
                fits[star] = (data, res, slope)
        summary[star] = entry
    # -------------------------------------------------------------------------
    # the two stars pull on each other: opposite signs, ratio of masses
    # -------------------------------------------------------------------------
    if len(fits) == 2:
        ratio = -fits['GJ 725 A'][2].mean() / fits['GJ 725 B'][2].mean()
        draws = -np.random.default_rng(1).choice(fits['GJ 725 A'][2], 20000) \
            / np.random.default_rng(2).choice(fits['GJ 725 B'][2], 20000)
        summary['mass_ratio'] = [float(v) for v in
                                 np.percentile(draws, [16, 50, 84])]
        # what the ratio would be had the secular accelerations been taken
        #   out a second time
        sec_a = summary['GJ 725 A']['secular'][0]
        sec_b = summary['GJ 725 B']['secular'][0]
        twice = -(np.random.default_rng(1).choice(fits['GJ 725 A'][2], 20000)
                  - sec_a) / (np.random.default_rng(2).choice(
                      fits['GJ 725 B'][2], 20000) - sec_b)
        summary['mass_ratio_twice'] = [float(v) for v in
                                       np.percentile(twice, [16, 50, 84])]
        log(f'  - a_orb(A) / a_orb(B) = M_B / M_A = {ratio:.3f} '
            f'({summary["mass_ratio"][0]:.3f} to '
            f'{summary["mass_ratio"][2]:.3f})', 'value')
    # -------------------------------------------------------------------------
    # the whole analysis of GL 725 B, the drift taken out
    # -------------------------------------------------------------------------
    data, res, slope = fits['GJ 725 B']
    total = summary['GJ 725 B']['linear']['total'][1]
    if args.skip_analysis and 'analysis' in previous:
        summary['analysis'] = previous['analysis']
        log('  analysis of GL 725 B kept from the previous run')
    else:
        flat = remove_perspective(data, total)
        flat.name = 'GL 725 B (SPIRou), drift removed'
        report = analyze(flat, os.path.join(outdir(DEMO), 'analysis'),
                         kmax=3, unit='both')
        kfip = report['fip']['koloa']
        summary['analysis'] = dict(
            period=float(report['period']),
            fip_at_period={key: float(res.fip_at(report['period']))
                           for key, res in report['fip'].items()},
            pk=[float(v) for v in kfip.pk],
            nflag=int(np.sum(kfip.outlier_prob > 0.5)),
            verdict=report['duck'].verdict if 'duck' in report else None)
    log(f'  analysis: best period {summary["analysis"]["period"]:.3f} d, '
        f'koloa FIP {summary["analysis"]["fip_at_period"]["koloa"]:.2e}, '
        f'{summary["analysis"]["nflag"]} exposures flagged, verdict '
        f'{summary["analysis"]["verdict"]}', 'value')
    save_summary(DEMO, summary)

    def drift_figure():
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.8))
        for ax, star in zip(axes, ('GJ 725 B', 'GJ 725 A')):
            if star not in fits:
                ax.set_visible(False)
                continue
            dat, fit, slp = fits[star]
            means = kplot.sequence_means(dat, dat.rv, fit.outlier_prob)
            ax.errorbar(means['time'], means['value'], means['err'],
                        fmt='o', ms=3, color=kplot.C['koloa'],
                        elinewidth=0.8, label='visit means')
            grid = np.linspace(dat.time.min(), dat.time.max(), 200)
            orb = summary[star]['linear']['orbital'][1]
            # the fitted line through its own zero point
            off = np.average(means['value'] - orb * (means['time'] - dat.tref)
                             / 365.25, weights=1 / means['err'] ** 2)
            ax.plot(grid, off + orb * (grid - dat.tref) / 365.25,
                    color=kplot.C['text'], lw=1.4,
                    label=f'fitted drift, {orb:+.2f} m s$^{{-1}}$ yr$^{{-1}}$')
            ax.set_xlabel('time [BJD - 2400000]')
            ax.set_ylabel('RV [m s$^{-1}$]')
            ax.set_title(f'{star} (SPIRou)', loc='left')
            ax.legend(loc='best', fontsize=7)
        fig.tight_layout()
        return fig

    save_both(drift_figure, 'drift', DEMO)
    log(f'Everything in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

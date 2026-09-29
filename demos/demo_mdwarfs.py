#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 15: three M dwarfs with more than 150 visits, cleaned by PCA2D.

Three campaigns reduced in Montreal, each with more than 150 visits:
GJ 687 and Barnard's star (GJ 699) with SPIRou, and Proxima with NIRPS
(team data, not in this repository: the demo skips a star whose series is
missing). Their telluric-corrected spectra were cleaned by pca2d-preclean
at its nominal setting (no star component, seven observer components, the
3-sigma cut and the 10-sigma bias of 2026-09-26) and measured by LBL, the
spectra as delivered beside them.

For each star koloa:
1. runs the whole analysis (FIP periodograms, the duck test on the best
   interval, the outliers of every visit and exposure);
2. fits the known planet it can see, with the Gaussian likelihood and with
   koloa's, sampled by MCMC with a linear trend: GJ 687 b (K = 6.14 m/s,
   Feng et al. 2020) and Proxima b (K = 1.23 m/s, Suarez Mascareno et al.
   2025), in the corrected series and (koloa only) in the delivered one;
3. measures the drift of each star from the velocities alone. The secular
   acceleration from Gaia DR3 (4.54 m/s/yr for Barnard's star) is already
   out of these velocities: APERO computes the barycentric correction with
   barycorrpy, whose formula (Wright & Eastman 2014) includes the proper
   motion and parallax of the star, so it is kept for reference only and
   never added to a fit;
4. maps what each series could have found (recovery rates, outlier-aware
   fold);
5. compares the corrected series with the delivered one.

    python demo_mdwarfs.py --workers 8

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
from koloa.completeness import recovery_map
from koloa.data import RVData
from koloa.fit import mcmc_orbits
from koloa.log import log
from koloa.periodogram import frequency_grid, gls, oap, window
from koloa.secular import (acceleration, gaia_astrometry,
                           perspective_from_astrometry)

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'mdwarfs'
GAIA = os.path.join(ROOT, 'data', 'gaia_dr3_astrometry.json')
TARGETS = [
    dict(key='gl687', name='GJ 687', gaia='GJ 687', inst='SPIRou',
         file='gl687_spirou_pca2d_0-7.rdb',
         delivered='gl687_spirou_delivered.rdb', mstar=0.40, apero=True,
         planets=[dict(name='b', period=38.142, Perr=0.007, K=6.14, Kerr=0.32,
                       ref='Feng et al. 2020')]),
    dict(key='proxima', name='Proxima', gaia='Proxima Cen', inst='NIRPS',
         file='proxima_nirps_pca2d_0-7.rdb',
         delivered='proxima_nirps_delivered.rdb', mstar=0.1221, apero=True,
         gp=dict(indicator='fwhm', guess=83.0, prange=(40.0, 150.0),
                 # the same fits with the prior of the FWHM fit restricted
                 #   to the published rotation, 83 to 92 d (--gp-check)
                 check=dict(guess=88.0, prange=(75.0, 110.0))),
         planets=[dict(name='b', period=11.18465, Perr=0.00053, K=1.226, Kerr=0.062,
                       ref='Suarez Mascareno et al. 2025')]),
    dict(key='gl699', name="Barnard's star", gaia='GJ 699', inst='SPIRou',
         file='gl699_spirou_pca2d_0-7.rdb',
         delivered='gl699_spirou_delivered.rdb', mstar=0.162, apero=True,
         planets=[]),
]
#: a campaign of more exposures than this is analysed on its visit means
NBIN = 1000
PERIODS = np.geomspace(1.5, 500.0, 11)
AMPS = np.geomspace(0.3, 30.0, 12)


# =============================================================================
# Define functions
# =============================================================================
def astrometry(name: str) -> dict:
    """Gaia DR3 astrometry, from the cache or the archive"""
    with open(GAIA) as handle:
        cache = json.load(handle)
    if name not in cache['stars']:
        cache['stars'][name] = gaia_astrometry(name)
        with open(GAIA, 'w') as handle:
            json.dump(cache, handle, indent=1)
    return cache['stars'][name]


def describe(data: RVData) -> dict:
    """The numbers of a series"""
    summ = data.summary()
    means = kplot.sequence_means(data, data.rv - np.median(data.rv))
    summ['nightly_rms'] = float(np.std(means['value']))
    # the same with a straight line removed: Barnard's star drifts by
    #   4.5 m/s/yr, which would otherwise be most of its scatter
    tmid = np.mean(data.time)
    coef = np.polyfit(means['time'] - tmid, means['value'], 1)
    summ['nightly_rms_line'] = float(np.std(
        means['value'] - np.polyval(coef, means['time'] - tmid)))
    coef = np.polyfit(data.time - tmid, data.rv, 1)
    resid = data.rv - np.polyval(coef, data.time - tmid)
    summ['rms_line'] = float(np.std(resid))
    summ['robust_rms_line'] = float(1.4826 * np.median(
        np.abs(resid - np.median(resid))))
    # the robust scatter of the visit means about a line fitted without
    #   the visits that stand out (five robust sigma, three passes)
    keep = np.ones(len(means['value']), dtype=bool)
    for _ in range(3):
        coef = np.polyfit(means['time'][keep] - tmid, means['value'][keep], 1)
        vres = means['value'] - np.polyval(coef, means['time'] - tmid)
        sig = 1.4826 * np.median(np.abs(vres[keep] - np.median(vres[keep])))
        keep = np.abs(vres - np.median(vres[keep])) < 5 * sig
    summ['nightly_robust_line'] = float(sig)
    return summ


def planet_fits(data: RVData, target: dict, sec, workers: int) -> dict:
    """Every known planet, gaussian and koloa, sampled"""
    out = {}
    for planet in target['planets']:
        per = planet['period']
        entry = dict(literature=planet)
        for name, like in (('gaussian', 'gaussian'), ('koloa', 'mixture')):
            res = mcmc_orbits(
                data, [dict(period=per, period_range=(0.98 * per, 1.02 * per),
                            eccentric=True)],
                likelihood=like, unit='sequence', trend=1,
                perspective=None if target.get('apero') else sec,
                nsteps=6000, nburn=1500, max_steps=60000, quiet=True)
            orb = res.orbits(mstar=target['mstar'])[0]
            entry[name] = dict(K=orb['K'], P=orb['P'], e=orb['e'],
                               msini=orb['msini'],
                               tau=res.diagnostics['tau_max'],
                               converged=res.diagnostics['converged'],
                               nflag=int(np.sum(res.outlier_prob > 0.5)))
            entry[name + '_fit'] = res
            log(f'  {target["name"]} {planet["name"]}, {name:9s}: K = '
                f'{orb["K"][0]:.2f} +{orb["K"][2]:.2f}/-{orb["K"][1]:.2f} m/s, '
                f'e = {orb["e"][0]:.2f}, P = {orb["P"][0]:.3f} d '
                f'(published {planet["K"]} +- {planet["Kerr"]} m/s)', 'value')
        out[planet['name']] = entry
    return out


def delivered_fits(target: dict, sec) -> dict:
    """The known planets in the delivered series, koloa only"""
    dpath = os.path.join(ROOT, 'data', target['delivered'])
    if not target['planets'] or not os.path.exists(dpath):
        return {}
    data = RVData.from_csv(dpath, name=f'{target["name"]} (delivered)')
    out = {}
    for planet in target['planets']:
        per = planet['period']
        res = mcmc_orbits(
            data, [dict(period=per, period_range=(0.98 * per, 1.02 * per),
                        eccentric=True)],
            likelihood='mixture', unit='sequence', trend=1,
            perspective=None if target.get('apero') else sec,
            nsteps=6000, nburn=1500, max_steps=60000, quiet=True)
        orb = res.orbits(mstar=target['mstar'])[0]
        out[planet['name']] = dict(K=orb['K'], P=orb['P'], e=orb['e'],
                                   msini=orb['msini'],
                                   tau=res.diagnostics['tau_max'],
                                   converged=res.diagnostics['converged'],
                                   nflag=int(np.sum(res.outlier_prob > 0.5)))
        log(f'  {target["name"]} {planet["name"]}, koloa on the delivered '
            f'series: K = {orb["K"][0]:.2f} +{orb["K"][2]:.2f}/'
            f'-{orb["K"][1]:.2f} m/s', 'value')
    return out


def gp_fits(target: dict, check: bool = False) -> dict:
    """
    The known planets with a GP of the rotation, on the visit means

    The period of the GP has its prior from an activity indicator (a GP
    fitted to it, koloa.fit.period_prior_from_indicator). A dense GP on
    every exposure costs hours of sampling on campaigns of this size, so
    the visit means are fitted (each visit is then an outlier unit).

    :param target: dict, an entry of TARGETS
    :param check: bool, use the guess and range of the indicator fit given
                  in target['gp']['check'] instead of the nominal ones
    """
    conf = target.get('gp')
    if not conf or not target['planets']:
        return {}
    if check:
        conf = dict(conf, **conf['check'])
    from koloa.fit import RVModel, period_prior_from_indicator
    out = {}
    for kind, fname in (('pca2d', target['file']),
                        ('delivered', target['delivered'])):
        path = os.path.join(ROOT, 'data', fname)
        if not os.path.exists(path):
            continue
        data = RVData.from_csv(path).binned()
        prior, ifit = period_prior_from_indicator(
            data, conf['indicator'], kernel='sho', period_guess=conf['guess'],
            nsteps=6000, nburn=2000, period_range=conf['prange'], trend=1,
            quiet=True)
        prot = float(np.exp(prior['log_period'][0]))
        # a chain of a few tens of effective samples cannot pin the period
        #   to a fraction of a per cent: the prior is kept 5 % wide or more
        prior['log_period'] = (prior['log_period'][0],
                               max(prior['log_period'][1], 0.05))
        entry = dict(prot=prot, prot_log_sigma=float(prior['log_period'][1]),
                     npoints=data.n)
        for planet in target['planets']:
            per = planet['period']
            model = RVModel(
                data, [dict(period=per, period_range=(0.98 * per, 1.02 * per),
                            eccentric=False)],
                gp=dict(kernel='sho', prior=prior,
                        init={'log_period': np.log(prot)}),
                likelihood='mixture', unit='point', trend=1)
            start = model.fit(nstart=2, quiet=True)
            res = model.sample(start=start, nsteps=6000, nburn=2000,
                               converge=True, max_steps=40000, quiet=True,
                               indicators='map')
            orb = res.orbits(mstar=target['mstar'])[0]
            diag = res.diagnostics
            entry[planet['name']] = dict(
                K=orb['K'], P=orb['P'], msini=orb['msini'],
                gp_period=res.param('gp_log_period'),
                nflag=int(np.sum(res.outlier_prob > 0.5)),
                tau_max=diag['tau_max'], converged=diag['converged'])
            log(f'  {target["name"]} {planet["name"]} with a GP of the '
                f'rotation ({kind}, visit means, prior {prot:.1f} d from '
                f'{conf["indicator"]}): K = {orb["K"][0]:.2f} '
                f'+{orb["K"][2]:.2f}/-{orb["K"][1]:.2f} m/s', 'value')
        out[kind] = entry
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--only', nargs='*', default=None)
    parser.add_argument('--gp-only', action='store_true',
                        help='only fit the known planets with a GP of the '
                             'rotation (targets that have one)')
    parser.add_argument('--gp-check', action='store_true',
                        help='with --gp-only, the GP fits with the check '
                             'prior of the rotation, kept apart')
    parser.add_argument('--describe-only', action='store_true',
                        help='only recompute the numbers of each series')
    parser.add_argument('--delivered-only', action='store_true',
                        help='only fit the known planets in the delivered '
                             'series, and add them to the summary')
    parser.add_argument('--allow-delivered', action='store_true',
                        help='analyse the delivered series of a star whose '
                             'PCA2D series is not there yet')
    args = parser.parse_args()
    warnings.simplefilter('ignore')
    log('Demo: three M dwarfs with more than 150 visits, PCA2D and koloa')
    summary = {}
    path_summary = os.path.join(outdir(DEMO), 'summary.json')
    if os.path.exists(path_summary):
        with open(path_summary) as handle:
            summary = json.load(handle)
    for target in TARGETS:
        if args.only and target['key'] not in args.only:
            continue
        if args.gp_only:
            key = target['key']
            if key not in summary or not target.get('gp'):
                continue
            fits = gp_fits(target, check=args.gp_check)
            with open(path_summary) as handle:
                summary = json.load(handle)
            summary[key]['gp_check' if args.gp_check else 'gp'] = fits
            save_summary(DEMO, summary)
            continue
        if args.describe_only:
            key = target['key']
            if key not in summary:
                continue
            path = os.path.join(ROOT, 'data', target['file'])
            if summary[key].get('kind', 'PCA2D') != 'PCA2D':
                path = os.path.join(ROOT, 'data', target['delivered'])
            summary[key]['series'] = describe(RVData.from_csv(path))
            dpath = os.path.join(ROOT, 'data', target['delivered'])
            if os.path.exists(dpath):
                summary[key]['delivered'] = describe(RVData.from_csv(dpath))
            save_summary(DEMO, summary)
            continue
        if args.delivered_only:
            key = target['key']
            if key not in summary or not target['planets']:
                continue
            sec = tuple(summary[key]['secular'])
            fits = delivered_fits(target, sec)
            with open(path_summary) as handle:
                summary = json.load(handle)
            for pname, val in fits.items():
                summary[key]['planets'][pname]['koloa_delivered'] = val
            save_summary(DEMO, summary)
            continue
        path = os.path.join(ROOT, 'data', target['file'])
        kind = 'PCA2D'
        if not os.path.exists(path):
            dpath0 = os.path.join(ROOT, 'data', target['delivered'])
            if args.allow_delivered and os.path.exists(dpath0):
                log(f'  {target["name"]}: no PCA2D series yet, the delivered '
                    f'one is analysed', 'warn')
                path, kind = dpath0, 'delivered'
            else:
                log(f'  {path} not there yet, {target["name"]} skipped',
                    'warn')
                continue
        key = target['key']
        data = RVData.from_csv(path, name=f'{target["name"]} '
                                          f'({target["inst"]}, {kind})')
        entry = dict(name=target['name'], inst=target['inst'],
                     series=describe(data), kind=kind)
        # was the PCA2D basis fitted on nightly coadds or on every exposure
        #   (a memory choice of pca2d-preclean; every exposure is corrected)
        source = path.replace('.rdb', '.source')
        if kind == 'PCA2D' and os.path.exists(source):
            with open(source) as handle:
                entry['basis'] = handle.read().strip()
        astro = astrometry(target['gaia'])
        sec = perspective_from_astrometry(astro)
        entry['secular'] = sec
        # APERO's barycentric correction has already removed it
        entry['secular_removed'] = bool(target.get('apero'))
        log(f'{target["name"]}: {data.n} exposures in {data.nseq} visits over '
            f'{data.baseline / 365.25:.1f} yr, rms {entry["series"]["rms"]:.2f}'
            f' m/s; secular acceleration {sec[0]:.4f} +- {sec[1]:.1e} m/s/yr',
            'value')
        dpath = os.path.join(ROOT, 'data', target['delivered'])
        if os.path.exists(dpath):
            entry['delivered'] = describe(RVData.from_csv(dpath))
            log(f'  delivered: rms {entry["delivered"]["rms"]:.2f} m/s, '
                f'robust {entry["delivered"]["robust_rms"]:.2f}; PCA2D: rms '
                f'{entry["series"]["rms"]:.2f}, robust '
                f'{entry["series"]["robust_rms"]:.2f}', 'value')
        # 1. the whole analysis: whole visits as the outlier units (the
        #   exposures of a sequence move together); a campaign of more than
        #   a thousand exposures is analysed on its visit means, and the GP
        #   test of the duck test, which costs an hour on such a series, is
        #   left out
        fdata = data if data.n <= NBIN else data.binned()
        report = analyze(fdata, os.path.join(outdir(DEMO), key), kmax=3,
                         unit='sequence' if fdata.nseq < fdata.n else 'point',
                         nsweep=1200, nburn=300, gp=False)
        entry['fip_on'] = 'exposures' if fdata is data else 'visit means'
        kfip = report['fip']['koloa']
        entry['analysis'] = dict(
            period=float(report['period']),
            fip_at_period={name: float(res.fip_at(report['period']))
                           for name, res in report['fip'].items()},
            pk=[float(val) for val in kfip.pk],
            nflag=int(np.sum(kfip.outlier_prob > 0.5)),
            verdict=report['duck'].verdict if 'duck' in report else None,
            peaks=[dict(period=float(pk['period']), fip=float(pk['fip']))
                   for pk in kfip.peaks[:5]])
        log(f'  analysis: best interval {report["period"]:.3f} d, koloa FIP '
            f'{entry["analysis"]["fip_at_period"]["koloa"]:.2e}, '
            f'{entry["analysis"]["nflag"]} {entry["fip_on"]} flagged', 'value')
        # 2. the known planets
        fits = planet_fits(data, target, sec, args.workers)
        for pname, pent in fits.items():
            for name in ('gaussian', 'koloa'):
                entry.setdefault('planets', {}).setdefault(pname, {})[name] = \
                    pent[name]
            entry['planets'][pname]['literature'] = pent['literature']
            best = pent['koloa_fit']

            def phase_fig(best=best, pname=pname):
                import matplotlib.pyplot as plt
                # the size it is printed at in the paper
                fig, ax = plt.subplots(figsize=(4.0, 3.0))
                kplot.phase(best, ax=ax, title=f'{target["name"]} {pname} '
                                               f'({target["inst"]}, {kind})')
                return fig
            save_both(phase_fig, f'{key}_phase_{pname}', DEMO)
        if kind == 'PCA2D':
            for pname, val in delivered_fits(target, sec).items():
                entry['planets'][pname]['koloa_delivered'] = val
            if target.get('gp'):
                entry['gp'] = gp_fits(target)
        # 3. the drift: the acceleration of the star, measured from the
        #   velocities alone (for APERO velocities, what the barycentric
        #   correction left)
        res = mcmc_orbits(data, [], likelihood='mixture', trend=1,
                          nsteps=4000, nburn=1000, quiet=True,
                          unit='sequence')
        entry['drift_measured'] = list(acceleration(res)['accel'])
        log(f'  drift from the velocities: {entry["drift_measured"][0]:.3f} '
            f'(-{entry["drift_measured"][1]:.3f} '
            f'+{entry["drift_measured"][2]:.3f}) m/s/yr, Gaia '
            f'{sec[0]:.4f}', 'value')
        # 4. what the series could have found (on the visit means, as the
        #   FIPs, for a campaign of more than a thousand exposures)
        entry['map_on'] = entry['fip_on']
        rmap = recovery_map(fdata, PERIODS, AMPS, ninj=12, test='fold',
                            outliers=True, nnull=150, seed=21,
                            workers=args.workers, trend=1, calibrate=True,
                            quiet=True)
        entry['completeness'] = dict(periods=rmap.periods.tolist(),
                                     k50=rmap.k_at(0.5).tolist(),
                                     k90=rmap.k_at(0.9).tolist(),
                                     false_alarm=rmap.false_alarm,
                                     thresholds=rmap.extra['thresholds'])
        save_both(lambda rmap=rmap: kplot.recovery_map(
            rmap, title=f'{target["name"]} ({target["inst"]}): recovery'),
            f'{key}_recovery', DEMO)
        # figures of the series; the outlier probabilities of visit means
        #   are given to each of their exposures
        prob = (kfip.outlier_prob if fdata is data
                else kfip.outlier_prob[data.seq])
        freq = frequency_grid(fdata.time, 1.1, None, 5)
        def series_fig():
            import matplotlib.pyplot as plt
            # the size it is printed at in the paper
            fig, ax = plt.subplots(figsize=(4.6, 2.6))
            kplot.timeseries(data, prob, ax=ax, title=f'{target["name"]} '
                                                      f'({target["inst"]}, '
                                                      f'{kind})')
            return fig
        save_both(series_fig, f'{key}_timeseries', DEMO)
        save_both(lambda: kplot.periodograms(
            freq, oap_gauss=oap(fdata, freq, outliers=False)['dlnl'],
            oap_mix=oap(fdata, freq, outliers=True)['dlnl'],
            window_power=window(fdata.time, freq),
            fips=dict(gaussian=report['fip']['gaussian'], koloa=kfip),
            mark=[pl['period'] for pl in target['planets']],
            title=f'{target["name"]} ({target["inst"]}, {kind})'),
            f'{key}_periodograms', DEMO)
        # another star may have been written meanwhile: merge, do not clobber
        if os.path.exists(path_summary):
            with open(path_summary) as handle:
                summary = json.load(handle)
        summary[key] = entry
        save_summary(DEMO, summary)
    log(f'Everything in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 14: the acceleration of a star, with its error.

"Secular acceleration" is used in two senses. The first is the acceleration
of the star itself, from a companion too far out to complete an orbit during
the series; the second, the perspective acceleration mu^2 d, is geometry.

1. An acceleration under outliers: a series on a real NIRPS sampling with a
   planet (K = 3 m/s at 11.2 d), a companion that accelerates the star by
   ACCEL m/s/yr, a visit jitter and outliers clear and borderline, NREAL
   times. The acceleration is fitted with the planet (a trend of degree one,
   read in m/s/yr by koloa.secular.acceleration) with the Gaussian
   likelihood and with koloa's, and with the Gaussian likelihood on the same
   series without its outliers: bias, width and coverage of its intervals.
2. GL 725 A and B, which pull on each other (demo_gl725b.py measured their
   accelerations): the smallest mass that each acceleration asks of the
   other star at their Gaia DR3 separation (Torres 1999), against the
   masses of the stars, and where the two accelerations and the masses
   place B along the line of sight.
3. The perspective acceleration of every star of the paper from its Gaia DR3
   astrometry (cached in data/gaia_dr3_astrometry.json; --refresh queries
   SIMBAD and the Gaia archive again), with the error propagated from the
   covariance of the proper motions and the parallax, and the drift it
   makes over each series. APERO (SPIRou, NIRPS) has already removed it.

    python demo_secular.py [--refresh] [--workers 8] [--nreal 24]
                           [--skip-ensemble]

Created on 2026-09-27

@author: artigau
"""
import argparse
import json
import os
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from _common import DATA, ROOT, outdir, save_both, save_summary
from koloa.data import RVData
from koloa.log import log
from koloa.secular import (GM_JUP, GM_SUN, companion_acceleration,
                           companion_min_mass, gaia_astrometry,
                           perspective_drift, perspective_from_astrometry)

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'secular'
GAIA = os.path.join(ROOT, 'data', 'gaia_dr3_astrometry.json')
MASSES = os.path.join(ROOT, 'data', 'gl725_masses.json')
#: the stars, and the series koloa analyses for them
STARS = {'GJ 699': None,
         'TOI-2120': 'toi2120_spirou_pca2d_0-7_bias10s3.rdb',
         'GJ 725 A': None,
         'GJ 725 B': 'gl725b_spirou_lbl_2026-09-24.rdb',
         'Kepler-21': 'kepler21_harpsn_dace_drs3.3.12.csv'}
PLANET = dict(P=11.2, K=3.0)
#: the acceleration a companion gives the simulated star [m/s/yr]
ACCEL = 2.0
#: the fits of each realisation: likelihood, and whether the outliers stay
FITS = (('gaussian', 'gaussian', True), ('koloa', 'mixture', True),
        ('gaussian_clean', 'gaussian', False))


# =============================================================================
# Define functions
# =============================================================================
def realisation(seed: int) -> dict:
    """One series with a planet, an accelerating companion and outliers,
    fitted three ways; the acceleration of every fit, and its draws"""
    warnings.simplefilter('ignore')
    from koloa import log as klog
    from koloa.fit import mcmc_orbits
    from koloa.secular import acceleration
    from koloa.simulate import REALISTIC, simulate
    klog.VERBOSE = False
    tpl = RVData.from_csv(DATA)
    rng = np.random.default_rng(seed)
    tperi = float(tpl.time[0] + rng.uniform(0, PLANET['P']))
    sim = simulate(planets=[dict(P=PLANET['P'], K=PLANET['K'], tp=tperi)],
                   template=tpl, outliers=REALISTIC, visit_jitter=2.0,
                   seed=seed)
    data = sim['data']
    drift = perspective_drift(data.time, ACCEL, data.tref)
    series = {True: data.with_values(data.rv + drift),
              False: data.with_values(data.rv + drift
                                      - sim['outlier_offset'])}
    out = dict(seed=seed, nclear=int(np.sum(sim['outlier_label'] == 'clear')),
               nborder=int(np.sum(sim['outlier_label'] == 'borderline')))
    for key, likelihood, with_outliers in FITS:
        res = mcmc_orbits(series[with_outliers],
                          [dict(period=PLANET['P'],
                                period_range=(0.98 * PLANET['P'],
                                              1.02 * PLANET['P']))],
                          likelihood=likelihood, eccentric=False,
                          unit='both', trend=1, nsteps=3000, nburn=1000,
                          max_steps=15000, seed=seed, quiet=True)
        acc = acceleration(res)
        draws = acc['draws']
        out[key] = dict(accel=list(acc['accel']),
                        q95=[float(v) for v in
                             np.percentile(draws, [2.5, 97.5])],
                        K=[float(v) for v in res.orbits()[0]['K']],
                        converged=bool(res.diagnostics['converged']),
                        tau_max=float(res.diagnostics['tau_max']))
    return out


def ensemble_stats(runs: list) -> dict:
    """Bias, width and coverage of the acceleration, per fit"""
    out = {}
    for key, _, _ in FITS:
        val = np.array([run[key]['accel'][0] for run in runs])
        low = np.array([run[key]['accel'][1] for run in runs])
        high = np.array([run[key]['accel'][2] for run in runs])
        q95 = np.array([run[key]['q95'] for run in runs])
        cover68 = int(np.sum((val - low <= ACCEL) & (ACCEL <= val + high)))
        cover95 = int(np.sum((q95[:, 0] <= ACCEL) & (ACCEL <= q95[:, 1])))
        # (median - truth) / sigma, sigma on the side of the truth
        zval = (val - ACCEL) / np.where(val > ACCEL, low, high)
        out[key] = dict(n=len(runs), median_error=float(np.median(val
                                                                  - ACCEL)),
                        median_sigma=float(np.median(0.5 * (low + high))),
                        cover68=cover68, cover95=cover95,
                        rms_z=float(np.sqrt(np.mean(zval ** 2))),
                        converged=int(sum(run[key]['converged']
                                          for run in runs)))
    return out


def ensemble_figure(runs: list):
    """The interval of the acceleration in every realisation, per fit"""
    import matplotlib.pyplot as plt
    from koloa import plotting as kplot
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    colors = dict(gaussian=kplot.C['gaussian'], koloa=kplot.C['koloa'],
                  gaussian_clean=kplot.C['muted'])
    labels = dict(gaussian='Gaussian, with the outliers',
                  koloa='koloa, with the outliers',
                  gaussian_clean='Gaussian, without them')
    for it, (key, _, _) in enumerate(FITS):
        xpos = np.arange(len(runs)) + (it - 1) * 0.25
        val = np.array([run[key]['accel'][0] for run in runs])
        err = np.array([run[key]['accel'][1:] for run in runs]).T
        ax.errorbar(xpos, val, yerr=err, fmt='o', ms=3, lw=0.8,
                    color=colors[key], label=labels[key])
    ax.axhline(ACCEL, color='k', lw=0.8, ls='--', label='true acceleration')
    ax.set_xlabel('realisation')
    ax.set_ylabel('acceleration [m s$^{-1}$ yr$^{-1}$]')
    ax.legend(fontsize=7, ncol=2, loc='upper center',
              bbox_to_anchor=(0.5, -0.2))
    fig.tight_layout()
    return fig


def separation_arcsec(ra1: float, dec1: float, ra2: float,
                      dec2: float) -> float:
    """The angular separation of two positions [deg] in arcsec (haversine)"""
    ra1, dec1, ra2, dec2 = np.radians([ra1, dec1, ra2, dec2])
    hav = (np.sin(0.5 * (dec2 - dec1)) ** 2 + np.cos(dec1) * np.cos(dec2)
           * np.sin(0.5 * (ra2 - ra1)) ** 2)
    return float(np.degrees(2 * np.arcsin(np.sqrt(hav))) * 3600.0)


def gl725(cache: dict) -> dict:
    """What the accelerations of GL 725 A and B say of each other"""
    from scipy.optimize import brentq
    path = os.path.join(outdir('gl725b'), 'summary.json')
    if not (os.path.exists(path) and os.path.exists(MASSES)):
        log('  GL 725: run demo_gl725b.py first', 'warn')
        return {}
    with open(path) as handle:
        summ = json.load(handle)
    with open(MASSES) as handle:
        masses = json.load(handle)
    star_a, star_b = cache['stars']['GJ 725 A'], cache['stars']['GJ 725 B']
    sep = separation_arcsec(star_a['ra'], star_a['dec'], star_b['ra'],
                            star_b['dec'])
    distance = 1000.0 / (0.5 * (star_a['parallax'] + star_b['parallax']))
    rho = sep * distance
    out = dict(separation=sep, distance=distance, rho_au=rho)
    # how far B lies behind A from the Gaia DR3 parallaxes alone [au]
    au_per_pc = 648000.0 / np.pi
    out['gaia_depth'] = [
        float((1000.0 / star_b['parallax'] - 1000.0 / star_a['parallax'])
              * au_per_pc),
        float(np.hypot(1000.0 * star_a['parallax_error']
                       / star_a['parallax'] ** 2,
                       1000.0 * star_b['parallax_error']
                       / star_b['parallax'] ** 2) * au_per_pc)]
    # the orbital (dynamical) acceleration of each star: percentiles
    #   16, 50, 84 of its draws [m/s/yr]; the bound grows with |a|
    for star, other in (('GJ 725 A', 'B'), ('GJ 725 B', 'A')):
        acc = summ[star]['linear']['orbital']
        mins = sorted(companion_min_mass(np.array(acc), distance, sep))
        mass = masses[f'mass_{other}']
        out[star] = dict(accel=acc, min_mass_other=mins, mass_other=mass)
        log(f'  {star}: {acc[1]:+.2f} m/s/yr asks of {other} at least '
            f'{mins[1]:.3f} (+{mins[2] - mins[1]:.3f}/-{mins[1] - mins[0]:.3f})'
            f' Msun ({mins[1] * GM_SUN / GM_JUP:.0f} MJup); {other} has '
            f'{mass[0]:.3f} +- {mass[1]:.3f} Msun', 'value')
    # both feel the same geometry: |a_A| + |a_B| = G (M_A + M_B) z / r^3
    #   gives the depth z of B behind A, two roots around rho / sqrt(2)
    total = abs(summ['GJ 725 A']['linear']['orbital'][1]) + \
        abs(summ['GJ 725 B']['linear']['orbital'][1])
    mtot = masses['mass_A'][0] + masses['mass_B'][0]
    peak = companion_acceleration(mtot, distance, sep)
    out['relative_accel'] = total
    out['relative_accel_max'] = peak
    if total < peak:
        def func(depth):
            return companion_acceleration(mtot, distance, sep, depth) - total
        near = brentq(func, 1e-6 * rho, rho / np.sqrt(2.0))
        far = brentq(func, rho / np.sqrt(2.0), 1e4 * rho)
        out['depth_au'] = [near, far]
        # the same with the accelerations and the masses drawn from their
        #   errors (percentiles 16, 50, 84 of each acceleration)
        rng = np.random.default_rng(3)
        draws = {'near': [], 'far': []}
        acc_a = summ['GJ 725 A']['linear']['orbital']
        acc_b = summ['GJ 725 B']['linear']['orbital']
        for _ in range(2000):
            tot = abs(rng.normal(acc_a[1], 0.5 * (acc_a[2] - acc_a[0]))) + \
                abs(rng.normal(acc_b[1], 0.5 * (acc_b[2] - acc_b[0])))
            mdraw = rng.normal(masses['mass_A'][0], masses['mass_A'][1]) + \
                rng.normal(masses['mass_B'][0], masses['mass_B'][1])
            if tot >= companion_acceleration(mdraw, distance, sep):
                continue

            def fdraw(depth, mdraw=mdraw, tot=tot):
                return companion_acceleration(mdraw, distance, sep,
                                              depth) - tot
            draws['near'].append(brentq(fdraw, 1e-6 * rho, rho / np.sqrt(2)))
            draws['far'].append(brentq(fdraw, rho / np.sqrt(2), 1e4 * rho))
        for key in ('near', 'far'):
            out[f'depth_{key}'] = [float(v) for v in
                                   np.percentile(draws[key], [16, 50, 84])]
        # A in front of B when B's velocity decreases (it falls towards us)
        order = ('A in front of B' if summ['GJ 725 B']['linear']['orbital'][1]
                 < 0 else 'B in front of A')
        out['order'] = order
        dnear, dfar = out['depth_near'], out['depth_far']
        log(f'  separation {sep:.2f} arcsec at {distance:.3f} pc '
            f'({rho:.1f} au): the accelerations and the masses place B '
            f'{near:.1f} (+{dnear[2] - dnear[1]:.1f}/-{dnear[1] - dnear[0]:.1f})'
            f' or {far:.1f} (+{dfar[2] - dfar[1]:.1f}/-{dfar[1] - dfar[0]:.1f})'
            f' au from A along the line of sight ({order})', 'value')
    else:
        log(f'  the accelerations ({total:.2f} m/s/yr together) exceed what '
            f'the masses can give at this separation ({peak:.2f})', 'warn')
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--refresh', action='store_true')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--nreal', type=int, default=24)
    parser.add_argument('--skip-ensemble', action='store_true',
                        help='keep the ensemble of the previous run')
    args = parser.parse_args()
    warnings.simplefilter('ignore')
    log('Demo: the acceleration of a star, with its error')
    path_summary = os.path.join(outdir(DEMO), 'summary.json')
    previous = {}
    if os.path.exists(path_summary):
        with open(path_summary) as handle:
            previous = json.load(handle)
    summary = {}
    # -------------------------------------------------------------------------
    # 1. an acceleration under outliers
    # -------------------------------------------------------------------------
    if args.skip_ensemble and 'ensemble' in previous:
        summary['ensemble'] = previous['ensemble']
    else:
        log(f'  {args.nreal} series with a companion drift of {ACCEL} '
            f'm/s/yr, a planet of K = {PLANET["K"]} m/s at {PLANET["P"]} d '
            f'and outliers, fitted three ways')
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            runs = list(pool.map(realisation, range(1, args.nreal + 1)))
        stats = ensemble_stats(runs)
        summary['ensemble'] = dict(accel=ACCEL, planet=PLANET,
                                   nreal=args.nreal, stats=stats, runs=runs)
        save_both(lambda: ensemble_figure(runs), 'ensemble_accel', DEMO)
    for key, st in summary['ensemble']['stats'].items():
        log(f'  {key:15s} median error {st["median_error"]:+.2f} m/s/yr, '
            f'median sigma {st["median_sigma"]:.2f}, true value in the 68 % '
            f'interval {st["cover68"]}/{st["n"]}, in the 95 % one '
            f'{st["cover95"]}/{st["n"]}, rms z {st["rms_z"]:.2f}', 'value')
    # -------------------------------------------------------------------------
    # 2. GL 725 A and B
    # -------------------------------------------------------------------------
    with open(GAIA) as handle:
        cache = json.load(handle)
    if args.refresh:
        for name in STARS:
            cache['stars'][name] = gaia_astrometry(name)
        with open(GAIA, 'w') as handle:
            json.dump(cache, handle, indent=1)
    summary['gl725'] = gl725(cache)
    # -------------------------------------------------------------------------
    # 3. the perspective acceleration of every star
    # -------------------------------------------------------------------------
    table = []
    for name, fname in STARS.items():
        astro = cache['stars'][name]
        value, error = perspective_from_astrometry(astro)
        _, drawn = perspective_from_astrometry(astro, ndraw=100000)
        row = dict(name=name, source_id=astro['source_id'],
                   parallax=astro['parallax'],
                   pm=float(np.hypot(astro['pmra'], astro['pmdec'])),
                   secular=value, error=error, error_drawn=drawn)
        if fname:
            data = RVData.from_csv(os.path.join(ROOT, 'data', fname))
            row['baseline'] = data.baseline
            row['drift'] = value * data.baseline / 365.25
        table.append(row)
        extra = (f', {row["drift"]:.2f} m/s over the {row["baseline"]:.0f} d'
                 f' of its series' if 'drift' in row else '')
        log(f'  perspective {name:12s} {value:8.4f} +- {error:.5f} m/s/yr '
            f'(draws {drawn:.5f}){extra}', 'value')
    summary['table'] = table
    save_summary(DEMO, summary)
    log(f'Everything in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

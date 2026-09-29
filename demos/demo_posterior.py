#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 9: the posterior of an eccentric orbit, with and without outliers
modelled.

An eccentric planet (P = 11.2 d, K = 10 m/s, e = 0.3, omega = 60 deg;
K = 6 m/s with --kamp 6) on the sampling and error bars of a real NIRPS
series, with a 2 m/s visit
jitter, bad visits (7 %) and spiked exposures (3 %) of both signs.
koloa.mcmc_orbits samples the orbit twice: with the gaussian likelihood
(the same white and visit jitters, but no outliers) and with koloa's
mixture. The same series without its outliers, fitted with the gaussian
likelihood, is the reference: what the data say when nothing went wrong.

The corner plots overlay the posteriors, one colour per likelihood. One
realisation can be lucky, so the demo then repeats the four fits (gaussian
and koloa, with and without outliers) on independent realisations, and
measures the bias and the coverage of the 68 and 95 % intervals.

    python demo_posterior.py                 # one realisation + ensemble
    python demo_posterior.py --kamp 6 --nreal 40 --workers 9
    python demo_posterior.py --figures-only  # redraw from what is saved

Every chain runs until it holds 50 autocorrelation times of every
parameter (koloa.mcmc_orbits does this by default).

Created on 2026-09-27

@author: artigau
"""
import argparse
import glob
import json
import os
import time
import warnings

import numpy as np

from _common import outdir, save_both, template
from koloa import kepler
from koloa import plotting as kplot
from koloa.fit import mcmc_orbits
from koloa.log import log
from koloa.simulate import REALISTIC, simulate

# =============================================================================
# Define variables
# =============================================================================
#: the planet (omega is that of the star, in the RadVel convention); K is
#: set on the command line
TRUTH = dict(P=11.2, K=10.0, e=0.3, omega_deg=60.0)
#: the outliers, both signs, clear and borderline (koloa.simulate)
OUTLIERS = REALISTIC
VISIT_JITTER = 2.0
#: the chains: at least this long, and extended to 50 autocorrelation times
NSTEPS, NBURN = 10000, 2000
#: the fits of every realisation: (name, likelihood, with outliers)
FITS = [('gaussian_clean', 'gaussian', False),
        ('koloa_clean', 'mixture', False),
        ('gaussian', 'gaussian', True),
        ('koloa', 'mixture', True)]
#: what is compared, and how it is labelled
NAMES = ['P', 'K', 'e', 'omega', 'dtc', 'jit', 'sjit']
LABELS = dict(P='$P$ [d]', K='$K$ [m s$^{-1}$]', e='$e$',
              omega=r'$\omega$ [deg]', dtc=r'$T_\mathrm{c} - T_\mathrm{c,true}$'
              ' [d]', jit=r'$\sigma_\mathrm{w}$ [m s$^{-1}$]',
              sjit=r'$\sigma_\mathrm{v}$ [m s$^{-1}$]')
PERCENTILES = [2.5, 16, 50, 84, 97.5]


# =============================================================================
# Define functions
# =============================================================================
def demo_name(kamp: float) -> str:
    """The output folder and figure prefix of one planet amplitude"""
    return 'posterior' if kamp == 10.0 else f'posterior_k{kamp:g}'


def realisation(seed: int, kamp: float):
    """One series: the planet, the noise, and the outliers"""
    tpl = template()
    rng = np.random.default_rng(1000 + seed)
    period = TRUTH['P']
    tperi = float(tpl.time[0] + rng.uniform(0, period))
    omega = np.radians(TRUTH['omega_deg'])
    sim = simulate(planets=[dict(P=period, K=kamp, e=TRUTH['e'],
                                 omega=omega, tp=tperi)],
                   template=tpl, outliers=OUTLIERS,
                   visit_jitter=VISIT_JITTER, seed=seed,
                   name=f'eccentric planet, seed {seed}')
    data = sim['data']
    clean = data.with_values(data.rv - sim['outlier_offset'])
    clean.name = data.name + ', no outliers'
    tc_true = kepler.tp_to_tc(tperi, period, TRUTH['e'], omega)
    return data, clean, tc_true, sim['outlier_label']


def physical(post: dict, tc_true: float) -> dict:
    """The quantities of the corner plot, from FitResult.posterior"""
    period = TRUTH['P']
    # a conjunction is defined modulo the period: the offset from the
    #   nearest true conjunction, in [-P/2, P/2)
    dtc = (post['tc_0'] - tc_true + 0.5 * period) % period - 0.5 * period
    omega = post['omega_0']
    return dict(P=post['P_0'], K=post['K_0'], e=post['e_0'],
                omega=TRUTH['omega_deg'] + (omega - TRUTH['omega_deg'] + 180)
                % 360 - 180,
                dtc=dtc, jit=post['jit_inst'], sjit=post['sjit'])


def truths(kamp: float) -> dict:
    """The true values of the corner quantities"""
    return dict(P=TRUTH['P'], K=kamp, e=TRUTH['e'],
                omega=TRUTH['omega_deg'], dtc=0.0, jit=0.0,
                sjit=VISIT_JITTER)


def run_one(seed: int, kamp: float, keep_chains: bool = False,
            quiet: bool = True):
    """The four fits of one realisation, summarised (and the chains)"""
    warnings.simplefilter('ignore')
    data, clean, tc_true, labels = realisation(seed, kamp)
    mask = labels != ''
    out = dict(seed=seed, noutlier=int(np.sum(mask)),
               nclear=int(np.sum(labels == 'clear')),
               nborder=int(np.sum(labels == 'borderline')), fits={})
    chains = {}
    for name, like, dirty in FITS:
        tstart = time.time()
        res = mcmc_orbits(data if dirty else clean, [TRUTH['P']],
                          likelihood=like, eccentric=True, unit='both',
                          nsteps=NSTEPS, nburn=NBURN, seed=seed, quiet=quiet)
        phys = physical(res.posterior(), tc_true)
        diag = res.diagnostics
        entry = dict(runtime=time.time() - tstart, nsteps=diag['nsteps'],
                     nburn=diag['nburn'],
                     tau_max=diag['tau_max'], tau_worst=diag['tau_worst'],
                     n_eff_min=diag['n_eff_min'],
                     converged=diag['converged'],
                     acceptance=diag['acceptance'])
        for key, val in phys.items():
            entry[key] = [float(v) for v in np.percentile(val, PERCENTILES)]
        if like == 'mixture':
            flagged = res.outlier_prob > 0.5
            entry['nflag_true'] = int(np.sum(flagged & mask))
            entry['nflag_clear'] = int(np.sum(flagged & (labels == 'clear')))
            entry['nflag_border'] = int(np.sum(flagged
                                               & (labels == 'borderline')))
            entry['nflag_false'] = int(np.sum(flagged & ~mask))
        out['fits'][name] = entry
        if keep_chains:
            chains[name] = phys
            chains[name + '_prob'] = res.outlier_prob
    if keep_chains:
        out['chains'] = chains
    return out


def clip_series(data, niter: int = 6, clip: float = 3.0):
    """
    The practitioner's baseline: a Gaussian orbit fit, the points beyond
    three robust sigma of its residuals clipped, the fit done again, until
    the selection stops changing. A clipped point is given an error bar a
    million times larger, so that the model is still evaluated at it.

    :return: np.ndarray, boolean mask of the points kept
    """
    from koloa.fit import RVModel
    keep = np.ones(data.n, dtype=bool)
    for _ in range(niter):
        work = data.select(np.arange(data.n))
        work.err = np.where(keep, data.err, data.err * 1e6)
        model = RVModel(work, [dict(period=TRUTH['P'], eccentric=True)],
                        likelihood='gaussian', unit='both')
        fit = model.fit(nstart=4, quiet=True)
        resid = data.rv - model.mean_model(fit.theta)
        centre = np.median(resid[keep])
        scale = 1.4826 * np.median(np.abs(resid[keep] - centre))
        new = np.abs(resid - centre) < clip * scale
        if np.array_equal(new, keep):
            break
        keep = new
    return keep


def run_clip(seed: int, kamp: float) -> dict:
    """The Gaussian fit after a 3-sigma clip of the residuals"""
    warnings.simplefilter('ignore')
    data, _, tc_true, labels = realisation(seed, kamp)
    mask = labels != ''
    tstart = time.time()
    keep = clip_series(data)
    res = mcmc_orbits(data.select(keep), [TRUTH['P']], likelihood='gaussian',
                      eccentric=True, unit='both', nsteps=NSTEPS,
                      nburn=NBURN, seed=seed, quiet=True)
    phys = physical(res.posterior(), tc_true)
    diag = res.diagnostics
    entry = dict(runtime=time.time() - tstart, nsteps=diag['nsteps'],
                 nburn=diag['nburn'], tau_max=diag['tau_max'],
                 tau_worst=diag['tau_worst'], n_eff_min=diag['n_eff_min'],
                 converged=diag['converged'], acceptance=diag['acceptance'],
                 nclip=int(np.sum(~keep)),
                 nclip_true=int(np.sum(~keep & mask)),
                 nclip_clear=int(np.sum(~keep & (labels == 'clear'))),
                 nclip_border=int(np.sum(~keep & (labels == 'borderline'))))
    for key, val in phys.items():
        entry[key] = [float(v) for v in np.percentile(val, PERCENTILES)]
    return dict(seed=seed, entry=entry)


def ensemble_stats(runs: list, kamp: float) -> dict:
    """Bias and coverage of every fit, over the realisations"""
    stats = {}
    tru = truths(kamp)
    fits = list(FITS)
    if all('gaussian_clip' in run['fits'] for run in runs):
        fits.append(('gaussian_clip', 'gaussian', True))
    for name, _, _ in fits:
        stats[name] = {}
        for key in ('P', 'K', 'e', 'omega', 'dtc'):
            pct = np.array([run['fits'][name][key] for run in runs])
            sig = 0.5 * (pct[:, 3] - pct[:, 1])
            zval = (pct[:, 2] - tru[key]) / sig
            in68 = (pct[:, 1] <= tru[key]) & (tru[key] <= pct[:, 3])
            in95 = (pct[:, 0] <= tru[key]) & (tru[key] <= pct[:, 4])
            stats[name][key] = dict(
                n=len(pct), median_error=float(np.median(pct[:, 2]
                                                         - tru[key])),
                median_sigma=float(np.median(sig)),
                median_z=float(np.median(zval)),
                rms_z=float(np.sqrt(np.mean(zval ** 2))),
                cover68=int(np.sum(in68)), cover95=int(np.sum(in95)))
        conv = [run['fits'][name]['converged'] for run in runs]
        stats[name]['converged'] = int(np.sum(conv))
        stats[name]['tau_max'] = float(np.max(
            [run['fits'][name]['tau_max'] for run in runs]))
        stats[name]['runtime'] = float(np.median(
            [run['fits'][name]['runtime'] for run in runs]))
        stats[name]['nsteps_max'] = int(np.max(
            [run['fits'][name]['nsteps'] for run in runs]))
    return stats


def ensemble_figure(runs: list, kamp: float):
    """The 68 % intervals of K and e in every realisation"""
    import matplotlib.pyplot as plt
    tru = truths(kamp)
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.4), sharex=True)
    order = np.argsort([run['fits']['gaussian_clean']['K'][2]
                        for run in runs])
    xpos = np.arange(len(runs))
    for ax, key in zip(axes, ('K', 'e')):
        for name, color, shift, label in (
                ('gaussian', kplot.C['gaussian'], -0.18,
                 'Gaussian, with outliers'),
                ('koloa', kplot.C['koloa'], 0.18, 'koloa, with outliers')):
            pct = np.array([runs[it]['fits'][name][key] for it in order])
            ax.errorbar(xpos + shift, pct[:, 2],
                        [pct[:, 2] - pct[:, 1], pct[:, 3] - pct[:, 2]],
                        fmt='o', color=color, ms=3.5, elinewidth=1.1,
                        mec=kplot.C['surface'], mew=0.5, label=label)
        ax.axhline(tru[key], color=kplot.C['text'], lw=0.9,
                   label='true value')
        ax.set_ylabel(LABELS[key])
    axes[0].legend(loc='upper left', ncol=3)
    axes[1].set_xlabel('realisation (sorted by the K of the Gaussian fit '
                       'without outliers)')
    axes[1].set_ylim(bottom=0)
    fig.tight_layout()
    return fig


def corner_figures(chains: dict, seed: int, kamp: float):
    """The two corner plots of the displayed realisation"""
    tru = truths(kamp)
    demo = demo_name(kamp)
    save_both(lambda: kplot.corner(
        [dict(samples=chains['gaussian'], label='Gaussian, with outliers',
              color='gaussian'),
         dict(samples=chains['koloa'], label='koloa, with outliers',
              color='koloa'),
         dict(samples=chains['gaussian_clean'],
              label='Gaussian, outliers removed (reference)',
              color='neutral', dashed=True)],
        NAMES, labels=LABELS, truths=tru, size=1.0,
        title=f'K = {kamp:g} m/s, outliers clear and borderline (seed {seed})'),
        'corner_outliers', demo)
    save_both(lambda: kplot.corner(
        [dict(samples=chains['gaussian_clean'], label='Gaussian',
              color='gaussian'),
         dict(samples=chains['koloa_clean'], label='koloa', color='koloa')],
        NAMES, labels=LABELS, truths=tru, size=1.0,
        title=f'the same series without its outliers (seed {seed})'),
        'corner_clean', demo)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kamp', type=float, default=10.0)
    parser.add_argument('--nreal', type=int, default=40)
    parser.add_argument('--workers', type=int, default=9)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--figures-only', action='store_true')
    parser.add_argument('--clip-only', action='store_true',
                        help='add the Gaussian fit after a 3-sigma clip to '
                             'every realisation already run, then the stats')
    args = parser.parse_args()
    warnings.simplefilter('ignore')
    kamp = args.kamp
    demo = demo_name(kamp)
    folder = outdir(demo)
    ensdir = os.path.join(folder, 'ensemble')
    os.makedirs(ensdir, exist_ok=True)
    if args.clip_only:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        todo = []
        for fname in sorted(glob.glob(os.path.join(ensdir, 'seed_*.json'))):
            with open(fname) as handle:
                run = json.load(handle)
            if 'gaussian_clip' not in run['fits']:
                todo.append(run['seed'])
        log(f'  the clipped Gaussian fit: {len(todo)} realisations on '
            f'{args.workers} workers')
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            jobs = [pool.submit(run_clip, seed, kamp) for seed in todo]
            for job in as_completed(jobs):
                res = job.result()
                fname = os.path.join(ensdir, f'seed_{res["seed"]:03d}.json')
                with open(fname) as handle:
                    run = json.load(handle)
                run['fits']['gaussian_clip'] = res['entry']
                with open(fname, 'w') as handle:
                    json.dump(run, handle, indent=1)
                log(f'    seed {res["seed"]}: {res["entry"]["nclip"]} points '
                    f'clipped, {res["entry"]["nclip_true"]} of them outliers')
        args.figures_only = True
    chain_file = os.path.join(folder, 'chains.npz')
    log('Demo: the posterior of an eccentric orbit, gaussian and koloa')
    log(f'  planet: P = {TRUTH["P"]} d, K = {kamp:g} m/s, e = '
        f'{TRUTH["e"]}, omega = {TRUTH["omega_deg"]} deg; visit jitter '
        f'{VISIT_JITTER} m/s', 'value')
    # -------------------------------------------------------------------------
    # 1. the displayed realisation, with its chains
    # -------------------------------------------------------------------------
    if not args.figures_only or not os.path.exists(chain_file):
        one = run_one(args.seed, kamp, keep_chains=True, quiet=False)
        chains = one.pop('chains')
        np.savez(chain_file, **{f'{name}__{key}': val
                                for name, phys in chains.items()
                                if not name.endswith('_prob')
                                for key, val in phys.items()},
                 **{name: val for name, val in chains.items()
                    if name.endswith('_prob')})
        with open(os.path.join(folder, 'one.json'), 'w') as handle:
            json.dump(one, handle, indent=1)
    with open(os.path.join(folder, 'one.json')) as handle:
        one = json.load(handle)
    saved = np.load(chain_file)
    chains = {}
    for key in saved.files:
        if '__' in key:
            name, qty = key.split('__')
            chains.setdefault(name, {})[qty] = saved[key]
    for name, _, _ in FITS:
        fit = one['fits'][name]
        log(f'  {name:15s} K = {fit["K"][2]:.2f} +{fit["K"][3] - fit["K"][2]:.2f}'
            f'/-{fit["K"][2] - fit["K"][1]:.2f} m/s, e = {fit["e"][2]:.3f} '
            f'+{fit["e"][3] - fit["e"][2]:.3f}/-{fit["e"][2] - fit["e"][1]:.3f}'
            f', omega = {fit["omega"][2]:.0f} deg; tau {fit["tau_max"]:.0f} '
            f'steps ({fit["tau_worst"]}), {fit["n_eff_min"]:.0f} effective '
            f'samples', 'value')
    corner_figures(chains, args.seed, kamp)
    # -------------------------------------------------------------------------
    # 2. the ensemble: bias and coverage
    # -------------------------------------------------------------------------
    if not args.figures_only:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        todo = [seed for seed in range(1, args.nreal + 1)
                if not os.path.exists(os.path.join(ensdir,
                                                   f'seed_{seed:03d}.json'))]
        log(f'  ensemble: {args.nreal} realisations, {len(todo)} to run on '
            f'{args.workers} workers')
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            jobs = {pool.submit(run_one, seed, kamp): seed for seed in todo}
            for job in as_completed(jobs):
                res = job.result()
                with open(os.path.join(ensdir, f'seed_{res["seed"]:03d}.json'),
                          'w') as handle:
                    json.dump(res, handle, indent=1)
                log(f'    seed {res["seed"]} done')
    runs = []
    for fname in sorted(glob.glob(os.path.join(ensdir, 'seed_*.json'))):
        with open(fname) as handle:
            runs.append(json.load(handle))
    if not runs:
        return
    stats = ensemble_stats(runs, kamp)
    for name, _, _ in FITS:
        for key in ('K', 'e', 'omega', 'P'):
            st = stats[name][key]
            log(f'  {name:15s} {key:5s}: median error {st["median_error"]:+.3f}'
                f' (median sigma {st["median_sigma"]:.3f}), rms z '
                f'{st["rms_z"]:.2f}, 68 % coverage {st["cover68"]}/{st["n"]},'
                f' 95 % {st["cover95"]}/{st["n"]}', 'value')
        log(f'  {name:15s} converged (50 tau) in {stats[name]["converged"]}'
            f'/{len(runs)}; longest tau {stats[name]["tau_max"]:.0f} steps; '
            f'median runtime {stats[name]["runtime"]:.0f} s', 'value')
    flags = [(run['fits']['koloa']['nflag_true'],
              run['fits']['koloa']['nflag_false'], run['noutlier'],
              run['fits']['koloa']['nflag_clear'], run['nclear'],
              run['fits']['koloa']['nflag_border'], run['nborder'])
             for run in runs]
    summary = dict(truth=dict(TRUTH, K=kamp), visit_jitter=VISIT_JITTER,
                   outliers=OUTLIERS,
                   nsteps=NSTEPS, nburn=NBURN, one=one, stats=stats,
                   nreal=len(runs),
                   flags=dict(true=int(sum(f[0] for f in flags)),
                              false=int(sum(f[1] for f in flags)),
                              outliers=int(sum(f[2] for f in flags)),
                              clear=int(sum(f[3] for f in flags)),
                              nclear=int(sum(f[4] for f in flags)),
                              border=int(sum(f[5] for f in flags)),
                              nborder=int(sum(f[6] for f in flags))))
    with open(os.path.join(folder, 'summary.json'), 'w') as handle:
        json.dump(summary, handle, indent=1)
    save_both(lambda: ensemble_figure(runs, kamp), 'ensemble', demo)
    log(f'Everything in {folder}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

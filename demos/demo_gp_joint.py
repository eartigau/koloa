#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 3: activity, a planet and bad visits, fitted jointly.

A spotted star (an SHO gaussian process, rotation 23 d), a planet of
K = 4 m/s at 9.7 d, a visit jitter and a few bad visits. An activity
indicator follows the spots (the same process, scaled, with its own noise).

1. The indicator gives the rotation period, as a prior for the GP of the
   velocities (koloa.fit.period_prior_from_indicator).
2. The GP and the orbit are fitted TOGETHER, with bad visits in the noise
   model: the orbit comes back.
3. The GP fitted first and the orbit to what is left (sequential_fit): the
   GP eats part of the planet, and koloa says so loudly.

The seed shown is typical: over four other realisations (seeds 5 to 8) the
joint fit gave K = 3.2, 6.2, 4.2 and 2.9 m/s (injected 4.0) and the
sequential fit 1.4, 1.5, 0.7 and 0.8 m/s.

    python demo_gp_joint.py

Created on 2026-09-27

@author: artigau
"""
import warnings

import numpy as np

from _common import (flag_summary, log_flags, outdir, save_both,
                     save_summary, template)
from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fit import RVModel, period_prior_from_indicator, sequential_fit
from koloa.log import KoloaWarning, log
from koloa.simulate import REALISTIC, activity_signal, simulate

DEMO = 'gp_joint'
SEED = 7
PROT, PPLANET, KPLANET = 23.0, 9.7, 4.0


def realisation(seed: int, quiet: bool = False):
    """One series (activity, a planet, outliers clear and borderline), its
    joint fit and its sequential fit"""
    tpl = template()
    act = activity_signal(tpl.time, sigma=5.0, period=PROT, quality=3.0,
                          seed=seed)
    sim = simulate(planets=[dict(P=PPLANET, K=KPLANET, tp=tpl.time[0] + 2)],
                   template=tpl, outliers=REALISTIC, visit_jitter=1.0,
                   seed=seed)
    rng = np.random.default_rng(seed + 40)
    ind_err = np.full(tpl.n, 1.5)
    indicator = 0.8 * act + ind_err * rng.normal(size=tpl.n)
    data = RVData(time=tpl.time, rv=sim['data'].rv + act, err=tpl.err,
                  seq=tpl.seq, name='activity + planet',
                  indicators={'dTemp': (indicator, ind_err)},
                  zero_point={'inst': 0.0})
    # 1. the rotation period from the indicator
    prior, _ = period_prior_from_indicator(data, 'dTemp', nsteps=2500,
                                           nburn=800)
    # 2. the joint fit
    gp = dict(kernel='sho', prior=prior)
    joint = RVModel(data, [dict(period=PPLANET)], gp=gp, likelihood='mixture',
                    unit='both').fit(nstart=1, quiet=quiet)
    alone = RVModel(data, [], gp=gp, likelihood='mixture',
                    unit='both').fit(nstart=1, quiet=quiet)
    # 3. the sequential fit, which should not be done
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', KoloaWarning)
        seq = sequential_fit(data, PPLANET, gp_prior=prior,
                             likelihood='mixture', unit='both')
    return dict(data=data, sim=sim, joint=joint, alone=alone, seq=seq,
                nwarn=len(caught))


def ensemble_one(seed: int) -> dict:
    """One more realisation: the joint and sequential K, and the flags"""
    warnings.simplefilter('ignore')
    from koloa import log as klog
    klog.VERBOSE = False
    more = realisation(seed, quiet=True)
    kj = more['joint'].orbits()[0]['K']
    ks = more['seq']['orbit'].orbits()[0]['K']
    return dict(seed=seed, joint=[float(v) for v in kj],
                sequential=[float(v) for v in ks],
                flags=flag_summary(more['joint'].outlier_prob,
                                   more['sim']['outlier_label']))


def ensemble(nreal: int, workers: int):
    """Many realisations: the bias and coverage of the two fits"""
    import json
    import os
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as pool:
        runs = list(pool.map(ensemble_one, range(100, 100 + nreal)))
    joint = np.array([run['joint'] for run in runs])
    seq = np.array([run['sequential'] for run in runs])
    # K with its lower and upper errors (MAP with Laplace intervals)
    zjoint = (joint[:, 0] - KPLANET) / np.where(joint[:, 0] > KPLANET,
                                                  joint[:, 1], joint[:, 2])
    clear = sum(run['flags']['clear']['n'] for run in runs)
    clear_flag = sum(run['flags']['clear']['flagged'] for run in runs)
    border = sum(run['flags']['borderline']['n'] for run in runs)
    border_flag = sum(run['flags']['borderline']['flagged'] for run in runs)
    stats = dict(n=nreal, joint_median=float(np.median(joint[:, 0])),
                 joint_mean=float(np.mean(joint[:, 0])),
                 joint_std=float(np.std(joint[:, 0])),
                 joint_cover68=int(np.sum(np.abs(zjoint) < 1)),
                 seq_median=float(np.median(seq[:, 0])),
                 seq_mean=float(np.mean(seq[:, 0])),
                 recall_clear=clear_flag / max(clear, 1),
                 recall_border=border_flag / max(border, 1))
    log(f'  {nreal} realisations: joint K median {stats["joint_median"]:.2f}'
        f' (mean {stats["joint_mean"]:.2f}, scatter {stats["joint_std"]:.2f})'
        f', within 1 sigma in {stats["joint_cover68"]}; sequential median '
        f'{stats["seq_median"]:.2f} m/s; clear outliers flagged '
        f'{100 * stats["recall_clear"]:.0f}%, borderline '
        f'{100 * stats["recall_border"]:.0f}%', 'value')
    path = os.path.join(outdir(DEMO), 'summary.json')
    with open(path) as handle:
        summ = json.load(handle)
    summ['ensemble'] = dict(stats=stats, runs=runs)
    save_summary(DEMO, summ)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ensemble', type=int, default=0,
                        help='only run this many more realisations, in '
                             'parallel, for the bias and coverage')
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    if args.ensemble:
        ensemble(args.ensemble, args.workers)
        return
    log('Demo: activity, a planet and outliers, clear and borderline')
    res = realisation(SEED)
    data, sim, joint = res['data'], res['sim'], res['joint']
    log(f'  injected: rotation {PROT} d, planet {PPLANET} d with K = '
        f'{KPLANET} m/s, {int(np.sum(sim["outlier_mask"]))} outlying '
        f'exposures', 'value')
    kj = joint.orbits()[0]['K']
    log(f'  joint fit: K = {kj[0]:.2f} +{kj[2]:.2f}/-{kj[1]:.2f} m/s, '
        f'GP period {np.exp(joint.param("gp_log_period")[0]):.1f} d; '
        f'the orbit gains {joint.logpost - res["alone"].logpost:.1f} in log '
        f'posterior over the GP alone', 'value')
    seq = res['seq']
    ks = seq['orbit'].orbits()[0]['K']
    log(f'  sequential fit: K = {ks[0]:.2f} +{ks[2]:.2f}/-{ks[1]:.2f} m/s '
        f'({res["nwarn"]} KoloaWarning raised)', 'value')
    flags = flag_summary(joint.outlier_prob, sim['outlier_label'])
    log_flags(flags)
    # four more realisations: is this one typical?
    extra = []
    for seed in range(SEED + 1, SEED + 5):
        more = realisation(seed, quiet=True)
        extra.append(dict(seed=seed,
                          joint=more['joint'].orbits()[0]['K'][0],
                          sequential=more['seq']['orbit'].orbits()[0]['K'][0]))
        log(f'  seed {seed}: joint K = {extra[-1]["joint"]:.2f}, sequential '
            f'K = {extra[-1]["sequential"]:.2f} m/s', 'value')
    save_summary(DEMO, dict(K_joint=kj, K_sequential=ks, K_true=KPLANET,
                            flags=flags, extra=extra))
    save_both(lambda: kplot.timeseries(
        data, joint.outlier_prob, fit=joint,
        title='Joint fit: GP (activity) + orbit, outliers flagged'),
        'timeseries', DEMO)
    save_both(lambda: kplot.phase(
        joint, title=f'Joint fit: K = {kj[0]:.2f} m/s (injected '
                     f'{KPLANET:.1f})'), 'phase_joint', DEMO)
    save_both(lambda: kplot.phase(
        seq['orbit'], title=f'Sequential fit: K = {ks[0]:.2f} m/s '
                            f'(injected {KPLANET:.1f})'),
        'phase_sequential', DEMO)
    log(f'Figures in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

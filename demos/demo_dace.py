#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 17: a known planetary system from DACE, every instrument at once.

HD 69830, a quiet K0 dwarf 12.6 pc from the Sun, has three Neptune-mass
planets, found with HARPS (Lovis et al. 2006), at 8.67, 31.6 and 197 days.
DACE serves its public velocities from four instrument eras: HARPS before
and after its 2015 fibre upgrade (HARPS03, HARPS15), and ESPRESSO before and
after its 2019 intervention (ESPRESSO18, ESPRESSO19); each era has its own
zero point, so each is an instrument. The demo:

1. fetches the velocities from DACE, anonymously (public data only), and
   keeps them in data/dace/ for later runs (koloa.dace.fetch; --refresh
   asks DACE again); one instrument per era (koloa.dace.rvdata);
2. fits the noise without planets, a white and a visit jitter per
   instrument, and inflates the errors of each instrument to it
   (koloa.fip.inflate_to_fit): the FIP sampler has one jitter for all;
3. runs the FIP with up to four signals (koloa.oafip) and keeps every
   interval with FIP < 1%;
4. fits those planets together by MCMC (koloa.mcmc_orbits), every
   instrument with its offset, jitters and outliers;
5. inflates the errors again, to the fit with the planets, and runs the
   FIP a second time: the noise of the first pass held the planets;
6. compares the orbits with the NASA Exoplanet Archive (fetched, and kept
   in data/dace/): its default solution of each planet names it, and K is
   set against every published solution, the most recent first (the
   archive's default may be the discovery paper, from far fewer
   velocities).

DACE filters some networks: run the fetch where it answers (the CSV then
serves every later run, anywhere).

    python demo_dace.py [--target HD69830] [--host "HD 69830"] [--refresh]
                        [--kmax 4] [--nsweep 1500] [--nburn 300]
                        [--nsteps 6000] [--eccentric] [--compare-only]

--compare-only compares the orbits of the last run with the archive again
(step 6 alone, after the archive has new solutions), without the FIP and the
MCMC.

Created on 2026-09-29

@author: artigau
"""
import argparse
import json
import os
import time
import warnings

import numpy as np

from _common import ROOT, outdir, save_both, save_summary
from koloa import plotting as kplot
from koloa.archive import compare, known_planets
from koloa.dace import fetch, instruments_table, rvdata
from koloa.fip import inflate_to_fit, oafip
from koloa.fit import RVModel, mcmc_orbits
from koloa.log import log
from koloa.utils import blas_threads

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'dace'
CACHE = os.path.join(ROOT, 'data', 'dace')
#: the detection threshold of the FIP
THRESHOLD = 0.01


# =============================================================================
# Define functions
# =============================================================================
def run_fip(data, fit, kmax, nsweep, nburn, seed, label):
    """
    The FIP of a series, its errors inflated instrument by instrument to a
    fit

    :return: tuple, the FIPResult and the inflation
    """
    inflated, info = inflate_to_fit(fit)
    log(f'{label}: errors inflated by ' + ', '.join(
        f'{inst} {val:.2f}' for inst, val in info['inflation'].items())
        + f' m/s (visit jitter left to the FIP: '
          f'{info["visit_jitter_ref"]:.2f} m/s)', 'value')
    with blas_threads(1):
        res = oafip(inflated, kmax=kmax, outliers='both', nsweep=nsweep,
                    nburn=nburn, nchains=2, seed=seed, progress=False)
    peaks = [p for p in res.peaks if p['fip'] < THRESHOLD]
    log(f'{label}: P(k) = ' + ', '.join(f'{val:.2f}' for val in res.pk)
        + '; intervals with FIP < 1%: ' + (', '.join(
            f'{p["period"]:.3f} d ({p["fip"]:.1e})' for p in peaks)
            or 'none'), 'value')
    return res, info


def fip_summary(res, info):
    """The numbers of a FIP, for the summary"""
    return dict(pk=[float(val) for val in res.pk],
                peaks=[dict(period=float(p['period']), fip=float(p['fip']))
                       for p in res.peaks],
                inflation=info['inflation'],
                visit_jitter_ref=info['visit_jitter_ref'],
                runtime=float(res.runtime))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('--target', default='HD69830',
                        help='the name DACE knows the target by')
    parser.add_argument('--host', default='HD 69830',
                        help='the host name in the NASA Exoplanet Archive')
    parser.add_argument('--refresh', action='store_true',
                        help='ask DACE and the archive again')
    parser.add_argument('--kmax', type=int, default=4)
    parser.add_argument('--nsweep', type=int, default=1500,
                        help='the sweeps of each FIP chain')
    parser.add_argument('--nburn', type=int, default=None,
                        help='the sweeps thrown away first (a fifth of '
                             'nsweep, at least 300, by default)')
    parser.add_argument('--eccentric', action='store_true')
    parser.add_argument('--nsteps', type=int, default=6000)
    parser.add_argument('--compare-only', action='store_true',
                        help='compare the orbits of the last run with the '
                             'archive again, without the FIP and the MCMC')
    args = parser.parse_args()
    warnings.simplefilter('ignore')
    start = time.time()
    log(f'Demo: {args.host} from DACE, every instrument')
    safe = args.target.replace(' ', '_')
    if args.compare_only:
        # 6. alone: the orbits of the last run against the archive
        with open(os.path.join(outdir(DEMO), 'summary.json')) as handle:
            summary = json.load(handle)
        known = known_planets(host=args.host, path=os.path.join(
            CACHE, f'{safe}_archive.json'), refresh=args.refresh)
        summary['archive'] = known
        summary['orbits'] = [
            compare({key: orb[key] for key in ('P', 'K', 'e', 'msini')},
                    known) for orb in summary['orbits']]
        save_summary(DEMO, summary)
        for orb in summary['orbits']:
            log(f'  {orb["P"][0]:.3f} d: K = {orb["K"][0]:.2f} m/s'
                + (f'; {orb["known"]["name"]}: ' + ', '.join(
                    f'{sol["K"]} +- {sol["K_err"]:.2f} ({sol["reference"]}, '
                    f'{sol["z"]:+.1f} sigma)' for sol in orb['comparisons'])
                   if 'known' in orb else '; not in the archive'), 'value')
        return
    # 1. the velocities, one instrument per era
    path = fetch(args.target, os.path.join(CACHE, f'{safe}_dace.csv'),
                 refresh=args.refresh)
    data = rvdata(path, name=args.host)
    known = known_planets(host=args.host, path=os.path.join(
        CACHE, f'{safe}_archive.json'), refresh=args.refresh)
    for pl in known['planets']:
        log(f'  known: {pl["name"]}, P = {pl["P"]:.3f} d, K = {pl["K"]} m/s '
            f'({pl["reference"]})', 'value')
    # 2. the noise without planets, a jitter per instrument
    noise = RVModel(data, [], likelihood='mixture', unit='both', trend=1,
                    seq_jitter='instrument').fit(nstart=2, quiet=True)
    # 3. the FIP, first pass
    nburn = args.nburn if args.nburn is not None else max(300,
                                                          args.nsweep // 5)
    fip1, info1 = run_fip(data, noise, args.kmax, args.nsweep, nburn, 1,
                          'FIP, noise without planets')
    found = sorted(p['period'] for p in fip1.peaks if p['fip'] < THRESHOLD)
    if not found:
        log('no interval with FIP < 1%: nothing to fit', 'error')
        return
    # 4. the planets, every instrument, by MCMC
    planets = [dict(period=per, period_range=(0.98 * per, 1.02 * per))
               for per in found]
    with blas_threads(1):
        fit = mcmc_orbits(data, planets, likelihood='mixture',
                          eccentric=args.eccentric, unit='both', trend=1,
                          seq_jitter='instrument', nsteps=args.nsteps,
                          nburn=args.nsteps // 3, nstart=2, seed=7,
                          quiet=True)
    mstar = known.get('star', {}).get('mass') or 1.0
    post = fit.posterior(mstar=mstar)

    def interval(key):
        low, mid, high = np.percentile(post[key], [16, 50, 84])
        return [float(mid), float(mid - low), float(high - mid)]
    orbits = []
    for ip, per in enumerate(found):
        # the known planet of that period, if any, and its K as published
        orb = compare(dict(P=interval(f'P_{ip}'), K=interval(f'K_{ip}'),
                           e=interval(f'e_{ip}'),
                           msini=interval(f'msini_{ip}')), known)
        orbits.append(orb)
        log(f'  {per:.3f} d: K = {orb["K"][0]:.2f} +{orb["K"][2]:.2f}/'
            f'-{orb["K"][1]:.2f} m/s, m sin i = {orb["msini"][0]:.1f} Me'
            + (f'; {orb["known"]["name"]}: ' + ', '.join(
                f'{sol["K"]} +- {sol["K_err"]:.2f} ({sol["reference"]}, '
                f'{sol["z"]:+.1f} sigma)' for sol in orb['comparisons'])
               if 'known' in orb else '; not in the archive'), 'value')
    # 5. the FIP again, the errors inflated to the fit with the planets
    fip2, info2 = run_fip(data, fit, args.kmax, args.nsweep, nburn, 2,
                          'FIP, noise with the planets')
    inst = data.inst
    flagged = {name: int(np.sum(fit.outlier_prob[inst == name] > 0.5))
               for name in data.instruments}
    theta, model = fit.theta, fit.model

    def jitter(pname):
        return (float(np.exp(theta[model.index[pname]]))
                if pname in model.index else None)
    noise_out = {name: dict(white=jitter(f'log_jit_{name}'),
                            visit=jitter(f'log_sjit_{name}'))
                 for name in data.instruments}
    summary = dict(target=args.target, host=args.host, n=int(data.n),
                   nvisits=int(data.nseq),
                   baseline_years=float(np.ptp(data.time) / 365.25),
                   instruments=instruments_table(data), archive=known,
                   fip_first=fip_summary(fip1, info1),
                   fip_second=fip_summary(fip2, info2),
                   found=[float(per) for per in found], orbits=orbits,
                   noise=noise_out, flagged=flagged,
                   eccentric=bool(args.eccentric),
                   mcmc=dict(converged=bool(fit.diagnostics['converged']),
                             tau_max=float(fit.diagnostics['tau_max'])),
                   runtime_minutes=(time.time() - start) / 60)
    save_summary(DEMO, summary)
    # the figures
    marks = [pl['P'] for pl in known['planets']]
    save_both(lambda: kplot.periodograms(
        fip2.freq, fips=dict(koloa=fip2), mark=marks,
        title=f'{args.host}: FIP of every instrument (known planets '
              f'dashed)'), 'periodograms', DEMO)
    save_both(lambda: kplot.timeseries(
        data, fit.outlier_prob, fit=fit,
        title=f'{args.host} from DACE: ' + ', '.join(data.instruments)),
        'timeseries', DEMO)
    for ip, orb in enumerate(orbits):
        save_both(lambda ip=ip, orb=orb: kplot.phase(
            fit, planet=ip, level='point',
            title=f'{args.host}, {orb["P"][0]:.2f} d: K = '
                  f'{orb["K"][0]:.2f} m/s'), f'phase_{ip}', DEMO)
    log(f'done in {summary["runtime_minutes"]:.0f} min; everything in '
        f'{outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

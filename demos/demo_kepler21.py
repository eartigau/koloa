#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 13: Kepler-21 with HARPS-N, a small transiting planet under activity.

Kepler-21 (HD 179070) is a bright F6 IV subgiant with a transiting
super-Earth, Kepler-21 b (P = 2.786 d; Howell et al. 2012), whose mass comes
from HARPS-N: K = 1.99 +- 0.65 m/s (Lopez-Morales et al. 2016) and
2.70 +- 0.46 m/s (Bonomo et al. 2023). The star rotates in about 13 days and
its activity dominates the velocities.

The public HARPS-N velocities (DACE, DRS 3.3.12, 98 exposures in 77 visits
over 5.4 years) are analysed with koloa:
1. the FIP periodograms of the velocities, gaussian and outlier-aware;
2. the rotation period from the S index, as the prior of a GP;
3. the orbit of b, its ephemeris from the transits (Bonomo et al. 2023),
   fitted jointly with the GP and the secular acceleration (Gaia DR3), with
   a gaussian likelihood and with koloa's, and sampled.

    python demo_kepler21.py

Created on 2026-09-27

@author: artigau
"""
import json
import os
import warnings

import numpy as np

from _common import ROOT, outdir, save_both, save_summary
from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fip import fip_comparison, oafip
from koloa.fit import RVModel, period_prior_from_indicator
from koloa.log import log
from koloa.periodogram import find_peaks, frequency_grid, gls, window
from koloa.secular import perspective_from_astrometry

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'kepler21'
PATH = os.path.join(ROOT, 'data', 'kepler21_harpsn_dace_drs3.3.12.csv')
GAIA = os.path.join(ROOT, 'data', 'gaia_dr3_astrometry.json')
#: Kepler-21 b from the transits (Bonomo et al. 2023; NASA Exoplanet
#: Archive default solution); BJD - 2400000
PLANET_B = dict(period=2.7858212, period_err=0.0000032, tc=55093.83716,
                tc_err=0.00085)
LITERATURE = {'Lopez-Morales et al. 2016': (1.99, 0.65),
              'Bonomo et al. 2023': (2.70, 0.46)}
INDICATORS = ['fwhm', 'smw', 'bis']


def main():
    warnings.simplefilter('ignore')
    log('Demo: Kepler-21 with HARPS-N')
    data = RVData.from_csv(PATH, name='Kepler-21')
    summ = data.summary()
    log(f'  {summ["npoints"]} exposures in {summ["nseq"]} visits over '
        f'{summ["baseline"]:.0f} d; rms {summ["rms"]:.1f} m/s, median '
        f'error {summ["median_err"]:.2f} m/s', 'value')
    with open(GAIA) as handle:
        sec = perspective_from_astrometry(json.load(handle)['stars']['Kepler-21'])
    out = dict(summary=summ, secular=sec, literature=LITERATURE)
    # -------------------------------------------------------------------------
    # 1. the FIP periodograms of the velocities
    # -------------------------------------------------------------------------
    freq = frequency_grid(data.time, 1.1, None, 10)
    comp = fip_comparison(data, freq=freq)
    koloa = oafip(data, kmax=3, outliers='both', freq=freq, nsweep=2000,
                  nburn=400, nchains=2, progress=False)
    out['fip'] = {key: dict(best_period=res.best()['period'],
                            best_fip=res.best()['fip'],
                            at_b=res.fip_at(PLANET_B['period']))
                  for key, res in list(comp.items()) + [('koloa', koloa)]}
    for key, val in out['fip'].items():
        log(f'  FIP {key:9s}: best {val["best_period"]:8.3f} d '
            f'({val["best_fip"]:.1e}); at b {val["at_b"]:.2e}', 'value')
    # -------------------------------------------------------------------------
    # 2. the rotation from the indicators
    # -------------------------------------------------------------------------
    rot = {}
    ifreq = frequency_grid(data.time, 2.0, 100.0, 20)
    for name in INDICATORS:
        series = data.indicator(name)
        power = gls(series.time, series.rv, series.err, ifreq)
        rot[name] = float(1 / ifreq[find_peaks(ifreq, power, 1)[0]])
        log(f'  {name}: highest GLS peak at {rot[name]:.2f} d', 'value')
    out['indicator_peaks'] = rot
    # the indicators of HARPS-N drift over the years (a slow change of the
    #   FWHM takes a free GP to hundreds of days): the rotation is looked
    #   for between 5 and 60 days, with a quadratic trend beside it, on the
    #   chromospheric S index, which an instrument moves least
    prior, ifit = period_prior_from_indicator(
        data, 'smw', kernel='sho', period_guess=rot['smw'], nsteps=10000,
        nburn=3000, period_range=(5.0, 60.0), trend=2)
    prot = float(np.exp(prior['log_period'][0]))
    # a chain of a few tens of effective samples cannot pin the period to a
    #   fraction of a per cent: the prior is kept at least 5 % wide
    prior['log_period'] = (prior['log_period'][0],
                           max(prior['log_period'][1], 0.05))
    out['prot'] = prot
    out['prot_log_sigma'] = float(prior['log_period'][1])
    out['prot_neff'] = ifit.diagnostics['n_eff_min']
    log(f'  rotation prior from the S index: {prot:.2f} d (log sigma '
        f'{prior["log_period"][1]:.3f}, {ifit.diagnostics["n_eff_min"]:.0f} '
        f'effective samples)', 'value')
    # -------------------------------------------------------------------------
    # 3. b, the GP and the secular acceleration, together
    # -------------------------------------------------------------------------
    gp = dict(kernel='sho', prior=prior, init={'log_period': np.log(prot)})
    fits = {}
    for name, like in (('gaussian', 'gaussian'), ('koloa', 'mixture')):
        model = RVModel(data, [dict(PLANET_B)], gp=gp, likelihood=like,
                        unit='both', perspective=sec)
        start = model.fit(nstart=2, quiet=True)
        # the outlier indicators at their most probable state (koloa flags
        #   none here), the rest sampled by emcee: the Gibbs sampler of the
        #   indicators mixes K too slowly on this series
        fits[name] = model.sample(start=start, nsteps=6000, nburn=2000,
                                  converge=True, max_steps=60000, quiet=True,
                                  indicators='map')
        orbit = fits[name].orbits(mstar=1.41)[0]
        diag = fits[name].diagnostics
        out[name] = dict(K=orbit['K'], msini=orbit['msini'],
                         jitter=fits[name].param('log_jit_inst'),
                         gp_sigma=fits[name].param('gp_log_sigma'),
                         nflag=int(np.sum(fits[name].outlier_prob > 0.5)),
                         tau_max=diag['tau_max'],
                         tau_worst=diag['tau_worst'],
                         tau_k=diag['tau']['k_0'],
                         neff_k=diag['n_eff']['k_0'],
                         converged=diag['converged'])
        log(f'  {name:9s} K_b = {orbit["K"][0]:.2f} +{orbit["K"][2]:.2f}/'
            f'-{orbit["K"][1]:.2f} m/s, m sin i = {orbit["msini"][0]:.1f} '
            f'Earth masses; {out[name]["nflag"]} exposures flagged; '
            f'{diag["sampler"]}, longest tau {diag["tau_max"]:.0f} '
            f'({diag["tau_worst"]}), tau of K {diag["tau"]["k_0"]:.0f} '
            f'({diag["n_eff"]["k_0"]:.0f} effective samples)', 'value')
    for ref, (kval, kerr) in LITERATURE.items():
        log(f'  {ref}: K_b = {kval:.2f} +- {kerr:.2f} m/s', 'value')
    save_summary(DEMO, out)
    best = fits['koloa']
    save_both(lambda: kplot.periodograms(
        freq, gls_power=gls(data.time, data.rv, data.err, freq),
        window_power=window(data.time, freq),
        fips=dict(gaussian=comp['gaussian'], koloa=koloa),
        mark=[PLANET_B['period'], prot],
        title='Kepler-21 (HARPS-N): b and the rotation (dashed)'),
        'periodograms', DEMO)
    save_both(lambda: kplot.timeseries(
        data, best.outlier_prob, fit=best,
        title='Kepler-21: GP (activity) + b, koloa'), 'timeseries', DEMO)
    save_both(lambda: kplot.phase(
        best, title=f'Kepler-21 b, koloa: K = {out["koloa"]["K"][0]:.2f} '
                    f'm/s'), 'phase_b', DEMO)
    log(f'Everything in {outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

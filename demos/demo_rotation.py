#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 18: the rotation of a star from its temperature indicator, with a
robust GP.

Spots on a rotating star change its temperature as they cross the disc, and
LBL measures that change spectrum by spectrum: DTEMP3500, in kelvins
(Artigau et al. 2022). GJ 687, an M3 dwarf 4.5 pc away, has 852 SPIRou
exposures in 214 visits over three years (cleaned by PCA2D); photometry
gives it a rotation of about 60 days (Burt et al. 2014). The demo:

1. the outlier-aware periodogram of DTEMP3500 (koloa.oap) and the gaussian
   one, and their highest peak between 5 and 150 d;
2. a GP with an SHO kernel (a damped oscillator: an amplitude sigma, an
   undamped period P0 and a quality factor Q) fitted to DTEMP3500 with
   koloa's likelihood (every exposure and every visit may be an outlier) and
   with a gaussian likelihood, both sampled by MCMC. A rotation is
   quasi-periodic, so Q is held above 1: below Q = 1/sqrt 2 the power of an
   SHO has no peak, and the GP can then trade the rotation for red noise
   (the maximum a posteriori with Q free down to 0.3, the default, is
   written beside, for comparison). The power peaks at
   P0 / sqrt(1 - 1 / (2 Q^2)), the period to compare with a periodogram;
3. the corner plot of the parameters of the SHO and of the noise, the two
   posteriors overlaid;
4. robustness: the same series with 10% of its visits made bad (every
   exposure of a visit moved by 4 to 8 times the dispersion of DTEMP3500),
   fitted with both likelihoods (maximum a posteriori), in several
   realisations;
5. the FWHM of the lines, another indicator of activity, for comparison.

The fits that do not depend on each other (the gaussian MCMC, the fits of
step 4 and the comparison of step 2) run in parallel, one process per core
(--ncpu), beside the MCMC with koloa's likelihood.

    python demo_rotation.py [--nsteps 8000] [--nreal 6] [--ncpu 9]

Created on 2026-09-29

@author: artigau
"""
import argparse
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from _common import ROOT, outdir, save_both, save_summary
from koloa import plotting as kplot
from koloa.data import RVData
from koloa.fit import Prior, RVModel, mcmc_orbits
from koloa.log import log
from koloa.periodogram import find_peaks, frequency_grid, oap, window
from koloa.utils import blas_threads

# =============================================================================
# Define variables
# =============================================================================
DEMO = 'rotation'
STAR = 'GJ 687'
PATH = os.path.join(ROOT, 'data', 'gl687_spirou_pca2d_0-7.rdb')
INSTRUMENT = 'SPIRou'
INDICATOR = 'DTEMP3500'
#: the rotation from photometry (Burt et al. 2014) [days]
LITERATURE = dict(P=60.0, reference='Burt et al. 2014')
#: the range of periods searched and allowed to the GP [days]
PERIODS = (5.0, 150.0)
#: the range of the quality factor of the SHO: at least 1, so that its
#: power has a peak (P_peak at most sqrt 2 P0), up to koloa's default 1000;
#: FREE_QUALITY is koloa's default range, for the comparison
QUALITY = (1.0, 1000.0)
FREE_QUALITY = (0.3, 1000.0)
#: the fraction of visits made bad, and their shift [dispersions]
BAD_FRACTION, BAD_SHIFT = 0.10, (4.0, 8.0)


# =============================================================================
# Define functions
# =============================================================================
def interval(draws):
    """the median and the distances to the 16th and 84th percentiles"""
    draws = np.asarray(draws, dtype=float)
    draws = draws[np.isfinite(draws)]
    low, mid, high = np.percentile(draws, [16, 50, 84])
    return [float(mid), float(mid - low), float(high - mid)]


def peak_period(period, quality):
    """
    The period at which the power of an SHO peaks, P0 / sqrt(1 - 1/(2 Q^2));
    nan where Q <= 1/sqrt(2) (the power then peaks at zero frequency)
    """
    period, quality = np.asarray(period), np.asarray(quality)
    arg = 1.0 - 1.0 / (2.0 * quality ** 2)
    return np.where(arg > 0, period / np.sqrt(np.clip(arg, 1e-12, None)),
                    np.nan)


def gp_model(likelihood, period_guess, quality=QUALITY):
    """the SHO GP of an indicator, its period in PERIODS and its quality
    factor in quality"""
    gp = dict(kernel='sho',
              prior={'log_period': Prior('uniform', np.log(PERIODS[0]),
                                         np.log(PERIODS[1])),
                     'log_quality': Prior('uniform', np.log(quality[0]),
                                          np.log(quality[1]))},
              init={'log_period': np.log(period_guess)})
    return dict(gp=gp, likelihood=likelihood, unit='both', trend=1)


def with_bad_visits(series, rng):
    """the series with a fraction of its visits moved together"""
    rv = series.rv + np.array([series.zero_point[str(inst)]
                               for inst in series.inst])
    visits = np.unique(series.seq)
    nbad = max(1, int(round(BAD_FRACTION * len(visits))))
    bad = rng.choice(visits, nbad, replace=False)
    scale = np.std(series.rv)
    for visit in bad:
        rv[series.seq == visit] += (rng.choice([-1, 1])
                                    * rng.uniform(*BAD_SHIFT) * scale)
    return RVData(series.time.copy(), rv, series.err.copy(),
                  inst=series.inst.copy(), seq=series.seq.copy(),
                  name=series.name), bad


def sample(series, likelihood, period_guess, nsteps):
    """
    The SHO GP of a series by MCMC

    :return: FitResult
    """
    warnings.simplefilter('ignore')
    with blas_threads(1):
        return mcmc_orbits(series, [], nsteps=nsteps, nburn=nsteps // 4,
                           max_steps=2 * nsteps, nstart=2, seed=3,
                           quiet=True, **gp_model(likelihood, period_guess))


def sample_draws(series, likelihood, period_guess, nsteps):
    """
    The SHO GP of a series by MCMC, in a process of its own: what the
    summary and the corner plot need, not the fit

    :return: dict, the draws (posterior) and whether the chain converged
    """
    fit = sample(series, likelihood, period_guess, nsteps)
    return dict(post=fit.posterior(),
                converged=bool(fit.diagnostics['converged']))


def best_fit(series, likelihood, period_guess, quality=QUALITY):
    """
    The maximum a posteriori of the SHO GP of a series

    :return: dict, P0, Q, the period of the peak of the power (nan without
             one) and the log posterior
    """
    warnings.simplefilter('ignore')
    with blas_threads(1):
        fit = RVModel(series, [], **gp_model(
            likelihood, period_guess, quality)).fit(nstart=2, quiet=True)
    index = fit.model.index
    per = float(np.exp(fit.theta[index['gp_log_period']]))
    qual = float(np.exp(fit.theta[index['gp_log_quality']]))
    return dict(P0=per, Q=qual, P_peak=float(peak_period(per, qual)),
                logpost=float(fit.logpost))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('--nsteps', type=int, default=8000)
    parser.add_argument('--nreal', type=int, default=6)
    parser.add_argument('--ncpu', type=int,
                        default=max(1, (os.cpu_count() or 2) - 1),
                        help='the processes beside the main one')
    args = parser.parse_args()
    warnings.simplefilter('ignore')
    start = time.time()
    log(f'Demo: the rotation of {STAR} from {INDICATOR}, with an SHO GP')
    data = RVData.from_csv(PATH, name=STAR)
    # the indicator as a series of its own: its values [K] and errors
    series = data.indicator(INDICATOR)
    log(f'{INDICATOR}: {series.n} values in {series.nseq} visits over '
        f'{series.baseline:.0f} d, dispersion {np.std(series.rv):.2f} K, '
        f'median error {np.median(series.err):.2f} K', 'value')
    # 1. the periodograms
    freq = frequency_grid(series.time, PERIODS[0], PERIODS[1], 10)
    dl_mix = np.asarray(oap(series, freq, unit='point', outliers=True)['dlnl'])
    dl_gau = np.asarray(oap(series, freq, unit='point', outliers=False)['dlnl'])
    best = float(1 / freq[find_peaks(freq, dl_mix, 1)[0]])
    best_gau = float(1 / freq[find_peaks(freq, dl_gau, 1)[0]])
    log(f'highest peak: {best:.2f} d outlier-aware (Delta lnL '
        f'{np.max(dl_mix):.1f}), {best_gau:.2f} d gaussian', 'value')
    # the FWHM of the lines, for comparison
    fwhm = data.indicator('fwhm')
    dl_fwhm = np.asarray(oap(fwhm, freq, unit='point', outliers=True)['dlnl'])
    best_fwhm = float(1 / freq[find_peaks(freq, dl_fwhm, 1)[0]])
    log(f'FWHM: highest peak at {best_fwhm:.2f} d', 'value')
    # 4. the series with bad visits, drawn here so that the realisations do
    #   not depend on the processes
    rng = np.random.default_rng(11)
    bad_series = [with_bad_visits(series, rng) for _ in range(args.nreal)]
    # the fits that do not depend on each other, in parallel: the gaussian
    #   MCMC, the maximum a posteriori of each series with bad visits and
    #   likelihood, and that of the clean series with Q free and held
    log(f'the gaussian MCMC and {2 * args.nreal + 2} fits on {args.ncpu} '
        f'processes, the MCMC with koloa\'s likelihood here')
    with ProcessPoolExecutor(max_workers=args.ncpu) as pool:
        gau_job = pool.submit(sample_draws, series, 'gaussian', best,
                              args.nsteps)
        bad_jobs = [{like: pool.submit(best_fit, bad, like, best)
                     for like in ('mixture', 'gaussian')}
                    for bad, _ in bad_series]
        free_job = pool.submit(best_fit, series, 'mixture', best,
                               FREE_QUALITY)
        held_job = pool.submit(best_fit, series, 'mixture', best)
        # 2. the SHO GP with koloa's likelihood, by MCMC, here: its fit
        #   draws the time series
        fit_mix = sample(series, 'mixture', best, args.nsteps)
        runs = [dict(nbad=int(len(bad)), **{like: job.result()
                                            for like, job in jobs.items()})
                for (_, bad), jobs in zip(bad_series, bad_jobs)]
        gau = gau_job.result()
        free, held_map = free_job.result(), held_job.result()
    posts = dict(mixture=fit_mix.posterior(), gaussian=gau['post'])
    converged = dict(mixture=bool(fit_mix.diagnostics['converged']),
                     gaussian=gau['converged'])
    for like, post in posts.items():
        post['gp_period_peak'] = peak_period(post['gp_period'],
                                             post['gp_quality'])
        post['log10_quality'] = np.log10(post['gp_quality'])
        pp = interval(post['gp_period_peak'])
        log(f'{like}: P0 = {np.median(post["gp_period"]):.2f} d, Q = '
            f'{np.median(post["gp_quality"]):.2f}, peak period = {pp[0]:.2f} '
            f'+{pp[2]:.2f}/-{pp[1]:.2f} d, sigma = '
            f'{np.median(post["gp_sigma"]):.2f} K', 'value')
    log(f'Q free down to {FREE_QUALITY[0]}: the maximum a posteriori has '
        f'P0 = {free["P0"]:.2f} d, Q = {free["Q"]:.2f} (log posterior '
        f'{free["logpost"]:.1f}, against {held_map["logpost"]:.1f} with Q '
        f'above {QUALITY[0]:.0f})', 'value')
    flagged = int(np.sum(fit_mix.outlier_prob > 0.5))
    log(f'koloa flags {flagged} of {series.n} DTEMP values as outliers',
        'value')
    for ireal, row in enumerate(runs):
        log(f'  realisation {ireal + 1}: {row["nbad"]} bad visits; peak '
            f'period {row["mixture"]["P_peak"]:.2f} d (koloa), '
            f'{row["gaussian"]["P_peak"]:.2f} d (gaussian)', 'value')
    ref = interval(posts['mixture']['gp_period_peak'])
    width = 0.5 * (ref[1] + ref[2])
    held = {like: int(sum(abs(row[like]['P_peak'] - ref[0]) <= width
                          for row in runs if np.isfinite(row[like]['P_peak'])))
            for like in ('mixture', 'gaussian')}
    log(f'within 1 sigma of the clean period ({ref[0]:.2f} +- {width:.2f} '
        f'd): koloa {held["mixture"]} of {args.nreal}, gaussian '
        f'{held["gaussian"]} of {args.nreal}', 'value')
    # the numbers
    names = ['gp_period', 'gp_period_peak', 'gp_quality', 'gp_sigma']
    summary = dict(
        star=STAR, instrument=INSTRUMENT, literature=LITERATURE,
        indicator=INDICATOR, n=int(series.n), nvisits=int(series.nseq),
        baseline=float(series.baseline), rms=float(np.std(series.rv)),
        median_err=float(np.median(series.err)), periods=list(PERIODS),
        quality=list(QUALITY), peak_mix=best,
        peak_mix_dlnl=float(np.max(dl_mix)), peak_gau=best_gau,
        peak_fwhm=best_fwhm, flagged=flagged,
        posteriors={like: {name: interval(posts[like][name])
                           for name in names + [key for key in posts[like]
                                                if key.startswith(('jit_',
                                                                   'sjit'))]}
                    for like in posts},
        free_quality=dict(range=list(FREE_QUALITY), map=free,
                          map_held=held_map),
        converged=converged,
        robustness=dict(nreal=args.nreal, fraction=BAD_FRACTION,
                        shift=list(BAD_SHIFT), runs=runs, held=held,
                        reference=ref),
        runtime_minutes=(time.time() - start) / 60)
    save_summary(DEMO, summary)
    # the figures
    peak_med = ref[0]
    save_both(lambda: kplot.periodograms(
        freq, oap_gauss=dl_gau, oap_mix=dl_mix,
        window_power=window(series.time, freq), mark=[peak_med, best_fwhm],
        title=f'{STAR}, {INDICATOR} ({INSTRUMENT}): the GP rotation and the '
              f'FWHM peak (dashed)'), 'periodograms', DEMO)
    save_both(lambda: kplot.timeseries(
        series, fit_mix.outlier_prob, fit=fit_mix,
        title=f'{STAR}, {INDICATOR}: the SHO GP, koloa',
        ylabel=f'{INDICATOR} [K]'), 'timeseries', DEMO)
    labels = {'gp_period': r'$P_0$ [d]', 'gp_period_peak': r'$P_\mathrm{peak}$ [d]',
              'log10_quality': r'$\log_{10} Q$', 'gp_sigma': r'$\sigma$ [K]'}
    jitters = [key for key in posts['mixture']
               if key.startswith(('jit_', 'sjit')) and key in posts['gaussian']]
    for key in jitters:
        labels[key] = ('visit jitter [K]' if key.startswith('sjit')
                       else 'white jitter [K]')
    corner_names = ['gp_period', 'gp_period_peak', 'log10_quality',
                    'gp_sigma'] + jitters
    save_both(lambda: kplot.corner(
        [dict(samples=posts['mixture'], label='koloa (outlier-aware)',
              color='koloa'),
         dict(samples=posts['gaussian'], label='gaussian likelihood',
              color='gaussian', dashed=True)],
        corner_names, labels=labels, size=1.35,
        title=f'{STAR}, {INDICATOR}: SHO GP'), 'corner', DEMO)

    def robustness_figure():
        fig, ax = kplot.plt.subplots(figsize=(6.2, 3.2))
        xval = np.arange(1, len(runs) + 1)
        ax.axhspan(ref[0] - ref[1], ref[0] + ref[2], color=kplot.C['koloa'],
                   alpha=0.12, lw=0, label='clean series, koloa (68%)')
        ax.axhline(ref[0], color=kplot.C['koloa'], lw=0.8)
        ax.plot(xval - 0.08, [row['mixture']['P_peak'] for row in runs], 'o',
                color=kplot.C['koloa'], ms=7, label='koloa, with bad visits')
        ax.plot(xval + 0.08, [row['gaussian']['P_peak'] for row in runs],
                's', color=kplot.C['gaussian'], ms=6,
                label='gaussian, with bad visits')
        ax.set_xlabel('realisation')
        ax.set_ylabel(r'$P_\mathrm{peak}$ of the SHO [d]')
        ax.set_xticks(xval)
        ax.set_title(f'{int(100 * BAD_FRACTION)}% of the visits moved by '
                     f'{BAD_SHIFT[0]:.0f} to {BAD_SHIFT[1]:.0f} dispersions',
                     loc='left', fontsize=9)
        ax.legend(fontsize=7.5, loc='best')
        return fig
    save_both(robustness_figure, 'robustness', DEMO)
    log(f'done in {summary["runtime_minutes"]:.0f} min; everything in '
        f'{outdir(DEMO)}')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

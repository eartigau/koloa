#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The whole analysis of a velocity series, in one call.

`analyze` runs, in order: the periodograms (GLS, the outlier-aware OAP and
the window), the FIPs with every way of handling outliers side by side
(none, soft clip, hard clip, and koloa's mixture, with and without it in
the multi-signal sampler), a fit of the best signal, the duck test, and
every figure, and writes a text report and a JSON summary next to them.

Created on 2026-09-27

@author: artigau
"""
import json
import os
from typing import Any, Dict, Optional

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.diagnostics import coherence, duck_test
from koloa.fip import fip_comparison, oafip
from koloa.fit import RVModel
from koloa.log import log
from koloa.periodogram import frequency_grid, gls, jackknife, oap, window


# =============================================================================
# Define functions
# =============================================================================
def _jsonable(obj: Any) -> Any:
    """Numbers and arrays made JSON-friendly"""
    if isinstance(obj, dict):
        return {str(key): _jsonable(val) for key, val in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(val) for val in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist() if obj.size < 200 else f'array({obj.shape})'
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, (float, int, str, bool)) or obj is None:
        return obj
    return str(obj)


def analyze(data: RVData, outdir: str = 'koloa_output', kmax: int = 3,
            unit: str = 'both', pmin: float = 1.1,
            pmax: Optional[float] = None, nsweep: int = 2000,
            nburn: int = 400, nchains: int = 2, duck: bool = True,
            gp: bool = True, style: str = 'paper', seed: int = 1,
            period: Optional[float] = None,
            nightly: Optional[bool] = None) -> Dict[str, Any]:
    """
    Everything koloa does, on one series

    :param data: RVData, the series
    :param outdir: str, where the figures and reports go
    :param kmax: int, the largest number of signals of the FIP sampler
    :param unit: str, what an outlier is: both (a lone exposure or a whole
                 visit), sequence or point
    :param pmin: float, the shortest period [days]
    :param pmax: float or None, the longest period [days]
    :param nsweep: int, recorded sweeps per chain of the FIP sampler
    :param nburn: int, burn-in sweeps
    :param nchains: int, the number of chains
    :param duck: bool, run the duck test on the best signal
    :param gp: bool, include the GP absorption test in the duck test
    :param style: str, paper (PDF) or web (SVG)
    :param seed: int, the seed
    :param period: float or None, the period to examine (the best of the
                   outlier-aware FIP, or of the gaussian FIP when nothing is
                   significant, when None) [days]
    :param nightly: bool, analyse the nightly means (RVData.nightly, koloa's
                    default: an outlier is then a night); False keeps the
                    exposures (an exposure or a whole visit can be one);
                    None is koloa.data.NIGHTLY

    :return: dict, every result
    """
    os.makedirs(outdir, exist_ok=True)
    kplot.set_style(style)
    if nightly is None:
        from koloa import data as kdata
        nightly = kdata.NIGHTLY
    if nightly:
        means = data.nightly()
        if means is not data:
            log(f'koloa: the nightly means, {data.n} exposures in {means.n} '
                f'nights (nightly=False keeps the exposures)')
            data = means
            if unit == 'both':
                unit = 'point'
    summ = data.summary()
    log(f'koloa: {data.name}, {summ["npoints"]} points in {summ["nseq"]} '
        f'sequences over {summ["baseline"]:.0f} days')
    log(f'  median error {summ["median_err"]:.2f} m/s, rms '
        f'{summ["rms"]:.2f} m/s, robust rms {summ["robust_rms"]:.2f} m/s',
        'value')
    out: Dict[str, Any] = dict(summary=summ, figures=[])
    # -------------------------------------------------------------------------
    # periodograms
    # -------------------------------------------------------------------------
    freq = frequency_grid(data.time, pmin, pmax, 10)
    log('Periodograms: GLS, OAP (gaussian and outlier-aware), window')
    out['gls'] = gls(data.time, data.rv, data.err, freq)
    punit = 'point' if unit == 'point' else 'sequence'
    out['oap_gauss'] = oap(data, freq, unit=punit, outliers=False)
    out['oap_mix'] = oap(data, freq, unit=punit, outliers=True)
    out['window'] = window(data.time, freq)
    # -------------------------------------------------------------------------
    # the FIPs
    # -------------------------------------------------------------------------
    log('FIP, one signal: no clip, soft clip, hard clip')
    comp = fip_comparison(data, pmin=pmin, pmax=pmax, freq=freq)
    for key, res in comp.items():
        best = res.best()
        log(f'  {res.method}: best {best["period"]:.4f} d, FIP '
            f'{best["fip"]:.2e}', 'value')
    log(f'FIP, up to {kmax} signal{"s" if kmax > 1 else ""}, gaussian noise (sampled)')
    gauss_multi = oafip(data, kmax=kmax, outliers=None, pmin=pmin,
                        pmax=pmax, nsweep=nsweep, nburn=nburn,
                        nchains=nchains, seed=seed, freq=freq)
    log(f'FIP, up to {kmax} signal{"s" if kmax > 1 else ""}, outlier-aware ({unit} outliers)')
    koloa_fip = oafip(data, kmax=kmax, outliers=unit, pmin=pmin, pmax=pmax,
                      nsweep=nsweep, nburn=nburn, nchains=nchains, seed=seed,
                      freq=freq)
    out['fip'] = dict(comp, gaussian_multi=gauss_multi, koloa=koloa_fip)
    prob = koloa_fip.outlier_prob
    nout = koloa_fip.settings['expected_outliers']
    log(f'  expected outlier points ({unit} model): {nout:.1f}; with P(outlier) > '
        f'0.5: {int(np.sum(prob > 0.5))}', 'value')
    for key, val in koloa_fip.rhat.items():
        level = 'value' if val < 1.1 else 'warn'
        log(f'  R-hat of {key}: {val:.3f}', level)
    # -------------------------------------------------------------------------
    # the period to examine
    # -------------------------------------------------------------------------
    if period is None:
        best = koloa_fip.best()
        if not best or best['fip'] > 0.5:
            best = comp['gaussian'].best()
            log(f'Nothing significant in the outlier-aware FIP: examining '
                f'the best gaussian peak, {best["period"]:.4f} d', 'warn')
        period = best['period']
    out['period'] = period
    # -------------------------------------------------------------------------
    # a fit at that period
    # -------------------------------------------------------------------------
    log(f'Fit of a circular orbit at {period:.4f} d (mixture, {unit})')
    fit = RVModel(data, [dict(period=period)], likelihood='mixture',
                  unit=unit).fit()
    out['fit'] = fit
    for line in fit.summary().split('\n'):
        log(line, 'value')
    jack = jackknife(data, freq, unit=punit)
    coh = coherence(data, period, prob=fit.outlier_prob)
    out['coherence'] = coh
    # -------------------------------------------------------------------------
    # figures
    # -------------------------------------------------------------------------
    name = data.name.replace(' ', '_')
    figs = []
    fig = kplot.timeseries(data, prob, title=f'{data.name}: velocities, '
                           f'coloured by outlier probability')
    figs.append(kplot.savefig(fig, os.path.join(outdir, f'{name}_rv')))
    fig = kplot.timeseries(data, prob, level='sequence',
                           title=f'{data.name}: one point per visit')
    figs.append(kplot.savefig(fig, os.path.join(outdir,
                                                f'{name}_rv_sequences')))
    fig = kplot.sequences(data, prob, fit=fit,
                          title=f'{data.name}: every visit')
    figs.append(kplot.savefig(fig, os.path.join(outdir, f'{name}_visits')))
    fips = dict(gaussian=comp['gaussian'], soft=comp['soft'],
                hard=comp['hard'], koloa=koloa_fip)
    fig = kplot.periodograms(freq, oap_gauss=out['oap_gauss']['dlnl'],
                             oap_mix=out['oap_mix']['dlnl'],
                             window_power=out['window'], fips=fips,
                             mark=[period], title=data.name)
    figs.append(kplot.savefig(fig, os.path.join(outdir,
                                                f'{name}_periodograms')))
    fig = kplot.jackknife(freq, jack, period, data, unit=punit,
                          title=f'{data.name}: leave one {punit} out')
    figs.append(kplot.savefig(fig, os.path.join(outdir, f'{name}_jackknife')))
    fig = kplot.phase(fit)
    figs.append(kplot.savefig(fig, os.path.join(outdir, f'{name}_phase')))
    fig = kplot.phase(fit, level='sequence',
                      title=f'P = {period:.4f} d, one point per visit')
    figs.append(kplot.savefig(fig, os.path.join(outdir,
                                                f'{name}_phase_sequences')))
    fig = kplot.coherence(coh, title=f'{data.name}: is the signal coherent?')
    figs.append(kplot.savefig(fig, os.path.join(outdir, f'{name}_coherence')))
    out['figures'] = figs
    # -------------------------------------------------------------------------
    # the duck test
    # -------------------------------------------------------------------------
    if duck:
        log('Duck test')
        report = duck_test(data, period, fipres=koloa_fip,
                           gauss_fip=comp['gaussian'], gp=gp, unit=unit)
        out['duck'] = report
    # -------------------------------------------------------------------------
    # reports
    # -------------------------------------------------------------------------
    lines = [f'koloa report: {data.name}', '=' * 60, '']
    lines += [f'{key}: {val}' for key, val in summ.items()]
    lines += ['', 'FIP at the examined period '
              f'({period:.4f} d), every method:']
    for key, res in out['fip'].items():
        lines.append(f'  {res.method:40s} {res.fip_at(period):.3e}')
    lines += ['', koloa_fip.summary(), '', gauss_multi.summary(), '',
              fit.summary()]
    if duck:
        lines += ['', out['duck'].text()]
    text = '\n'.join(lines)
    with open(os.path.join(outdir, f'{name}_report.txt'), 'w') as handle:
        handle.write(text + '\n')
    jsum = dict(summary=summ, period=period,
                fip_at_period={key: res.fip_at(period)
                               for key, res in out['fip'].items()},
                peaks={key: res.peaks for key, res in out['fip'].items()},
                pk=koloa_fip.pk, expected_outliers=nout,
                outlier_prob=prob, rhat=koloa_fip.rhat,
                orbit=fit.orbits()[0],
                coherence={key: coh[key] for key in coh if key != 'chunks'},
                verdict=out['duck'].verdict if duck else None,
                figures=figs)
    with open(os.path.join(outdir, f'{name}_summary.json'), 'w') as handle:
        json.dump(_jsonable(jsum), handle, indent=1)
    log(f'Report and {len(figs)} figures in {outdir}')
    return out


# =============================================================================
# End of code
# =============================================================================

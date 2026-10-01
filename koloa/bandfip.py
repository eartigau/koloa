#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The FIP by period bands, with a GP of the activity only as flexible as the
data ask for, from the longest periods down.

A local GP of scale L has its power at the long periods: at a period P, its
power relative to zero frequency is exp(-2 pi^2 L^2 / P^2), 0.7 % for
L = P/2. A GP free to take a short L therefore explains a slow planet as
well as the planet does, and the FIP never needs the planet (on G 203-42, a
200-day planet of K = 15 m/s had a FIP of 0.66 with such a GP, 1e-19 with a
GP held at L >= 182 d).

So the periods are cut in bands, and each band is decided by a FIP whose GP
cannot reach it: its scale at least alpha times the longest period of the
band. Going from the longest band down, the GP may become more flexible
only when the data ask for it:

1. for the band, two FIPs over every period (only the band is read from
   them, but the true period of a signal must be on the grid for its aliases
   in the band to be counted as aliases): one with the GP kept so far (none
   at first), one with a local GP whose shortest scale is alpha times the
   longest period of the band;
2. the more flexible one is kept when it gains more than `threshold` in the
   marginal likelihood (the 90th percentile of the chains' ln Z, the
   sinusoids of the FIP in both, so that a planet does not ask for a GP);
3. the band takes its FIP from the model kept.

The GP is thus fitted from the top down, and stops where the activity
stops.

Created on 2026-10-01

@author: artigau
"""
from typing import Any, Dict, List, Optional

import numpy as np

from koloa.data import RVData
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: the shortest GP scale of a band, as a fraction of its longest period
ALPHA = 0.5
#: the gain in ln Z that makes the GP more flexible
THRESHOLD = 3.0


# =============================================================================
# Define functions
# =============================================================================
def _lnz(res) -> float:
    """the 90th percentile of the chains' ln Z (a type-II maximum, robust)"""
    logz = np.asarray(res.chains.get('logz', []), dtype=float)
    return float(np.percentile(logz, 90)) if len(logz) else np.nan


def spec(lmin: Optional[float], baseline: float) -> Any:
    """the GP of a band: a local GP whose scale is at least lmin (None: no
    GP)"""
    if lmin is None:
        return None
    return [dict(kind='local', length=(float(lmin),
                                       float(max(3 * baseline, 2 * lmin))))]


def banded_fip(data: RVData, pmin: float = 1.1, pmax: Optional[float] = None,
               nband: int = 6, alpha: float = ALPHA,
               threshold: float = THRESHOLD, kmax: int = 2,
               nsweep: int = 1000, nburn: int = 300, nchains: int = 2,
               seed: int = 1, decided: Optional[List[Optional[float]]] = None,
               quiet: bool = False, label: Optional[str] = None,
               trend: int = 1) -> Dict[str, Any]:
    """
    The FIP by period bands, the GP of each band only as flexible as needed
    (see the module)

    Every FIP covers every period, so one run with a given GP serves every
    band that keeps that GP: the descent costs one FIP per band at most,
    plus the first one without a GP.

    :param data: RVData, the series
    :param pmin: float, the shortest period [days]
    :param pmax: float or None, the longest (twice the baseline when None)
    :param nband: int, the period bands (log-spaced)
    :param alpha: float, the shortest GP scale of a band over its longest
                  period
    :param threshold: float, the gain in ln Z that makes the GP more flexible
    :param kmax: int, the signals of each FIP
    :param nsweep: int, recorded sweeps per chain
    :param nburn: int, burn-in sweeps
    :param nchains: int, chains per FIP
    :param seed: int, the seed
    :param decided: list or None, the shortest GP scale of each band, the
                    longest band first (a previous descent's 'decided'): no
                    descent, one FIP per GP
    :param quiet: bool, no log lines
    :param trend: int, the degree of the polynomial in time of each FIP
    :param label: str or None, what the progress of each FIP is shown as
                  (with which FIP of how many at most, and its GP)

    :return: dict, bands (per band, the longest first: low, high, gp (the
             shortest scale kept, None for no GP), gain, result (the FIP
             kept)), decided, the periodogram put together (freq, fip,
             family, period, width), and the settings
    """
    from koloa.fip import oafip
    pmax = pmax or 2 * data.baseline
    edges = np.geomspace(pmin, pmax, nband + 1)
    runs: Dict[Any, Any] = {}
    # the FIPs at most: one per band and the one without a GP, or one per
    #   GP decided
    nmax = nband + 1 if decided is None else len(set(decided))

    def run(lmin):
        """the FIP with the GP of shortest scale lmin (once)"""
        key = None if lmin is None else round(float(lmin), 6)
        if key not in runs:
            what = 'no GP' if lmin is None else f'GP L >= {lmin:.0f} d'
            runs[key] = oafip(data, kmax=kmax, nsweep=nsweep, nburn=nburn,
                              nchains=nchains, progress=False, pmin=pmin,
                              pmax=pmax, gp=spec(lmin, data.baseline),
                              seed=seed + 7 * len(runs), trend=trend,
                              label=(f'{label}, FIP {len(runs) + 1} of up '
                                     f'to {nmax} ({what})' if label
                                     else None))
        return runs[key]
    kept: Optional[float] = None
    bands: List[Dict[str, Any]] = []
    for ib in range(nband - 1, -1, -1):
        low, high = float(edges[ib]), float(edges[ib + 1])
        lmin = alpha * high
        gain = np.nan
        if decided is not None:
            kept = decided[nband - 1 - ib]
        else:
            gain = _lnz(run(lmin)) - _lnz(run(kept))
            if gain > threshold:
                kept = lmin
        bands.append(dict(low=low, high=high, gp=kept, gain=float(gain),
                          candidate_lmin=lmin, result=run(kept)))
        if not quiet:
            log(f'band {low:7.2f} to {high:7.2f} d: '
                + ('' if decided is not None else
                   f'GP L >= {lmin:.1f} d gains {gain:+.1f} in ln Z -> ')
                + (f'GP L >= {kept:.1f} d' if kept else 'no GP'), 'value')
    # the periodogram put together: each frequency from its band's FIP
    freq, fip, family, density = [], [], [], []
    for band in bands:
        res = band['result']
        per = 1.0 / res.freq
        sel = (per >= band['low']) & (per < band['high'])
        freq.append(res.freq[sel])
        fip.append(res.fip[sel])
        density.append(res.density[sel])
        fam = res.family if res.family is not None else res.fip
        family.append(fam[sel])
    order = np.argsort(np.concatenate(freq))
    freq = np.concatenate(freq)[order]
    return dict(bands=bands, edges=edges, freq=freq, width=1.0 / data.baseline,
                fip=np.concatenate(fip)[order],
                family=np.concatenate(family)[order],
                density=np.concatenate(density)[order], period=1.0 / freq,
                decided=[band['gp'] for band in bands], nrun=len(runs),
                settings=dict(alpha=alpha, threshold=threshold, kmax=kmax,
                              nsweep=nsweep, nburn=nburn, nchains=nchains))


def as_result(banded: Dict[str, Any], npeaks: int = 5):
    """
    The banded FIP as one FIPResult (what the rest of koloa reads): the
    periodogram put together, the peaks of each band from its own FIP, the
    outlier probabilities averaged over the FIPs, and P(k) and the chains of
    the band of the strongest peak

    :return: FIPResult, with settings['bands'] (low, high, gp, gain)
    """
    from koloa.fip import FIPResult
    peaks = []
    for band in banded['bands']:
        for peak in band['result'].peaks:
            if band['low'] <= peak['period'] < band['high']:
                peaks.append(dict(peak))
    peaks.sort(key=lambda peak: peak['fip'])
    peaks = peaks[:npeaks]
    best = band_of(banded, peaks[0]['period'])['result'] if peaks else \
        banded['bands'][-1]['result']
    runs = {id(band['result']): band['result'] for band in banded['bands']}
    probs = [res.outlier_prob for res in runs.values()
             if res.outlier_prob is not None]
    settings = dict(best.settings)
    settings.pop('gp', None)
    settings['bands'] = [dict(low=band['low'], high=band['high'],
                              gp=band['gp'], gain=band['gain'])
                         for band in banded['bands']]
    settings['banded'] = dict(banded['settings'])
    return FIPResult(freq=banded['freq'], fip=banded['fip'],
                     density=banded['density'],
                     method=best.method.split(', GP')[0] + ', GP by band',
                     pk=best.pk, outlier_prob=(np.mean(probs, axis=0)
                                               if probs else None),
                     unit=best.unit, chains=best.chains, peaks=peaks,
                     settings=settings, family=banded['family'])


def gp_text(bands: List[Dict[str, Any]]) -> str:
    """the GP of each band, in words"""
    return '; '.join(f'{band["low"]:.1f}-{band["high"]:.1f} d: '
                     + (f'L >= {band["gp"]:.0f} d' if band['gp'] else 'none')
                     for band in bands)


def band_of(result: Dict[str, Any], period: float) -> Dict[str, Any]:
    """the band of a period in a banded_fip result"""
    for band in result['bands']:
        if band['low'] <= period < band['high']:
            return band
    return result['bands'][0] if period >= result['bands'][0]['high'] \
        else result['bands'][-1]


def fip_at(result: Dict[str, Any], period: float, family: bool = True
           ) -> float:
    """
    The FIP of a period (or of the period OR any of its aliases), from the
    FIP of its band

    :return: float
    """
    res = band_of(result, period)['result']
    width = result['width']
    if family and res.family is not None:
        return res.family_containing(period, width)
    return res.fip_containing(period, width)


# =============================================================================
# End of code
# =============================================================================

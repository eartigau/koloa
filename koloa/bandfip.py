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


def banded_fip(data: RVData, pmin: float = 1.1, pmax: Optional[float] = None,
               nband: int = 6, alpha: float = ALPHA,
               threshold: float = THRESHOLD, kmax: int = 2,
               nsweep: int = 1000, nburn: int = 300, nchains: int = 2,
               seed: int = 1, quiet: bool = False) -> Dict[str, Any]:
    """
    The FIP by period bands, the GP of each band only as flexible as needed
    (see the module)

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
    :param quiet: bool, no log lines

    :return: dict, bands (per band, the longest first: edges, gp (the shortest
             scale kept, None for no GP), gain, the FIP kept), and the
             periodogram put together: freq, fip, family, period
    """
    from koloa.fip import oafip
    pmax = pmax or 2 * data.baseline
    edges = np.geomspace(pmin, pmax, nband + 1)
    scale_max = 3 * data.baseline
    kept: Any = None
    kept_lmin: Optional[float] = None
    bands: List[Dict[str, Any]] = []
    kw = dict(kmax=kmax, nsweep=nsweep, nburn=nburn, nchains=nchains,
              progress=False, pmin=pmin)
    for ib in range(nband - 1, -1, -1):
        low, high = float(edges[ib]), float(edges[ib + 1])
        lmin = alpha * high
        cand = [dict(kind='local', length=(lmin, max(scale_max, 2 * lmin)))]
        with_kept = oafip(data, pmax=pmax, gp=kept, seed=seed + ib, **kw)
        flexible = oafip(data, pmax=pmax, gp=cand, seed=seed + ib, **kw)
        gain = _lnz(flexible) - _lnz(with_kept)
        if gain > threshold:
            kept, kept_lmin, res = cand, lmin, flexible
        else:
            res = with_kept
        bands.append(dict(low=low, high=high, gp=kept_lmin, gain=float(gain),
                          result=res, candidate_lmin=lmin))
        if not quiet:
            log(f'band {low:7.2f} to {high:7.2f} d: GP L >= {lmin:.1f} d gains '
                f'{gain:+.1f} in ln Z -> '
                + (f'GP L >= {kept_lmin:.1f} d' if kept_lmin else 'no GP'),
                'value')
    # the periodogram put together: each frequency from its band's FIP
    freq, fip, family = [], [], []
    for band in bands:
        res = band['result']
        per = 1.0 / res.freq
        sel = (per >= band['low']) & (per < band['high'])
        freq.append(res.freq[sel])
        fip.append(res.fip[sel])
        fam = res.family if res.family is not None else res.fip
        family.append(fam[sel])
    order = np.argsort(np.concatenate(freq))
    freq = np.concatenate(freq)[order]
    return dict(bands=bands, edges=edges, freq=freq, width=1.0 / data.baseline,
                fip=np.concatenate(fip)[order],
                family=np.concatenate(family)[order], period=1.0 / freq,
                settings=dict(alpha=alpha, threshold=threshold, kmax=kmax,
                              nsweep=nsweep, nburn=nburn, nchains=nchains))


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

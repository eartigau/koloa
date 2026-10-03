#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A transit in the TESS light curve at the period of a velocity signal (or
of a TOI): the light curve high-passed, folded, and a box searched where
the transit would be

    lc = transit.from_csv('archives/GJ_436/phot/tess.csv')
    res = transit.search(lc, period=2.6439, t0=54406.86, t0_err=0.01,
                         mstar=0.45, rstar=0.42)
    res['best']['snr'], res['best']['depth'], res['plausible']

- The light curve: each sector's flux [ppt] (koloa.gather, phot/tess.csv;
  or koloa.tess.light_curves), split where it has a gap of more than half
  a day; flares clipped (more than 4 sigma above a running median); the
  rest high-passed, a running median of a window three times the expected
  transit (at least half a day) taken out, so that a transit stays and the
  spots and the systematics of TESS go.
- The expected transit: its duration from the period and the star (a
  central transit, P / pi * asin(R* / a), a from Kepler's third law); its
  time, when there is one (the conjunction of a velocity fit, or the
  ephemeris of a TOI), carried to the epochs of TESS with its error,
  sqrt(t0_err^2 + (n P_err)^2): the search is within 3 sigma of it (and
  at least one duration), or over the whole phase when that is wider.
- The search: boxes of 0.5, 1 and 2 times the expected duration, their
  centres a quarter of the shortest apart; the depth of a box is minus the
  mean flux in it, its error the scatter of the light curve averaged over
  the box's duration (the red noise of TESS at that time scale, measured
  on the light curve itself) over the square root of the number of
  transits it holds, or, when larger, the scatter of the depths of the
  same box over every phase of the fold (where there is no transit: the
  null it is set against; an active star's residuals hold more than the
  red noise of a bin); its signal-to-noise ratio, depth over error.
- The trend again, the transits found left out of it (interpolated across
  them: on the steep slopes of an active star a running median follows a
  part of a transit), and the same box measured again (not searched
  again, which would pick the noise twice).
- Plausible: a box with a signal-to-noise ratio of 7 or more (TESS's own
  threshold is 7.1), seen in 2 transits or more, still 3.5 or more
  without its deepest transit (one event, a flare's dip or a systematic
  of one orbit, is not a transit seen twice), as deep as a planet can
  make it (Rp < 2.5 Jupiter radii: deeper, an eclipsing binary). The
  transits of the other known planets of the star (the archive, the
  TOIs) are left out first, so that they do not line up by chance at
  another period. Last, the null: the same search at 20 periods where
  there is nothing to find (0.7 to 1.4 times the period, the same width
  of window at a random time, the box's own transits left out), each
  judged as the box is, without its deepest transit; a transit beats the
  best box of every one
  (an active star, its flares and spots, makes dips as deep as a planet's
  at any period: then the light curve cannot tell).
- The radius of the planet from the depth, Rp = R* sqrt(depth), and the
  depths of a 1 Earth-radius and a 1 Jupiter-radius planet before the star
  (the IAU 2015 nominal radii) for the plot.

Created on 2026-10-03

@author: artigau
"""
import csv
import math
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

# =============================================================================
# Define variables
# =============================================================================
#: the radii of the Earth and Jupiter (equatorial) and of the Sun [m]: the
#: IAU 2015 nominal values (Resolution B3)
R_EARTH = 6.3781e6
R_JUPITER = 7.1492e7
R_SUN = 6.957e8
#: the gravitational parameter of the Sun [m^3/s^2] (IAU 2015 nominal)
GM_SUN = 1.3271244e20
#: a box this many sigma deep is a plausible transit (TESS: 7.1)
SNR_PLAUSIBLE = 7.0
#: seen this many times at least
MIN_TRANSITS = 2
#: as deep as a planet can make it: Rp below this [Jupiter radii]
RP_MAX_RJ = 2.5
#: a gap that splits a light curve [days]
GAP = 0.5
#: the shortest high-pass window [days]
WINDOW_MIN = 0.5
#: the durations of the boxes, in expected durations
DURATIONS = (0.5, 1.0, 2.0)
#: the periods of the null: a transit beats the best box of each
NTRIAL = 20


# =============================================================================
# Define functions
# =============================================================================
def from_csv(path: str) -> Optional[Dict[str, np.ndarray]]:
    """
    The light curve koloa.gather keeps (phot/tess.csv: rjd, flux [ppt],
    sflux, sector, pipeline)

    :return: dict, time [BJD - 2400000], flux, err [ppt], sector; or None
    """
    cols: Dict[str, List[float]] = dict(time=[], flux=[], err=[], sector=[])
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle):
            try:
                vals = (float(row['rjd']), float(row['flux']),
                        float(row['sflux']), float(row['sector']))
            except (KeyError, TypeError, ValueError):
                continue
            if all(np.isfinite(vals)):
                for key, val in zip(cols, vals):
                    cols[key].append(val)
    if len(cols['time']) < 100:
        return None
    return {key: np.asarray(val) for key, val in cols.items()}


def from_light_curves(lcs: Dict[str, Any]) -> Optional[Dict[str, np.ndarray]]:
    """the light curves of koloa.tess.light_curves as one (its sectors
    side by side), or None"""
    sectors = lcs.get('sectors') or []
    if not sectors:
        return None
    return dict(time=np.concatenate([lc['time'] for lc in sectors]),
                flux=np.concatenate([lc['flux'] for lc in sectors]),
                err=np.concatenate([lc['err'] for lc in sectors]),
                sector=np.concatenate([np.full(len(lc['time']),
                                               float(lc['sector']))
                                       for lc in sectors]))


def depth_of(radius_earth: float, rstar: float) -> float:
    """the depth [ppt] of a planet of a radius [Earth radii] before a star
    [solar radii]"""
    return 1e3 * (radius_earth * R_EARTH / (rstar * R_SUN)) ** 2


def radius_of(depth: float, rstar: float) -> float:
    """the radius [Earth radii] of a planet that makes a depth [ppt] before
    a star [solar radii]"""
    return math.sqrt(max(depth, 0.0) / 1e3) * rstar * R_SUN / R_EARTH


def duration(period: float, mstar: float, rstar: float) -> float:
    """
    The duration [hours] of a central transit: P / pi * asin(R* / a), a
    from Kepler's third law (the planet's mass neglected)

    :param period: float [days]
    :param mstar: float [solar masses]
    :param rstar: float [solar radii]
    """
    sec = period * 86400.0
    axis = (GM_SUN * mstar * sec ** 2 / (4 * math.pi ** 2)) ** (1 / 3)
    ratio = min(rstar * R_SUN / axis, 1.0)
    return period * 24.0 / math.pi * math.asin(ratio)


def _segments(time: np.ndarray, sector: np.ndarray) -> List[np.ndarray]:
    """the indices of each stretch of a light curve: a sector, split at
    its gaps"""
    order = np.lexsort((time, sector))
    cut = np.where((np.diff(time[order]) > GAP)
                   | (np.diff(sector[order]) != 0))[0] + 1
    return [seg for seg in np.split(order, cut) if len(seg)]


def highpass(lc: Dict[str, np.ndarray], window: float,
             mask: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """
    The light curve high-passed: in each stretch (a sector, split at its
    gaps), a running median of a window [days] taken out; flares (more
    than 4 sigma above it) clipped. The points of a mask (the transits)
    are left out of the running median, interpolated across them: on a
    steep slope of an active star, a running median would follow a part
    of the transit and make it shallower

    :param lc: dict, the light curve
    :param window: float, the window [days]
    :param mask: np.ndarray of bool or None, the points out of the trend

    :return: dict, time, flux (high-passed, ppt), err, sector; sorted in
             time
    """
    from scipy.ndimage import median_filter
    keep_t, keep_f, keep_e, keep_s = [], [], [], []
    for seg in _segments(lc['time'], lc['sector']):
        time, flux = lc['time'][seg], lc['flux'][seg]
        if len(time) < 5:
            continue
        step = float(np.median(np.diff(time))) if len(time) > 1 else 1.0
        use = ~mask[seg] if mask is not None else np.ones(len(seg), bool)
        if use.sum() < 5:
            use = np.ones(len(seg), bool)
        # an odd number of points, no more than the stretch has
        size = max(3, int(round(window / max(step, 1e-6))) | 1)
        size = min(size, int(use.sum()) - 1 + int(use.sum()) % 2)
        trend = median_filter(flux[use], size=size, mode='nearest')
        rest = flux - np.interp(time, time[use], trend)
        sig = 1.4826 * np.median(np.abs(rest - np.median(rest)))
        good = rest < 4 * sig if sig > 0 else np.ones(len(rest), bool)
        keep_t.append(time[good])
        keep_f.append(rest[good])
        keep_e.append(lc['err'][seg][good])
        keep_s.append(lc['sector'][seg][good])
    if not keep_t:
        return dict(time=np.array([]), flux=np.array([]), err=np.array([]),
                    sector=np.array([]))
    out = dict(time=np.concatenate(keep_t), flux=np.concatenate(keep_f),
               err=np.concatenate(keep_e), sector=np.concatenate(keep_s))
    order = np.argsort(out['time'])
    return {key: val[order] for key, val in out.items()}


def red_noise(lc: Dict[str, np.ndarray], width: float) -> float:
    """
    The scatter of a high-passed light curve averaged over a time [days]:
    the noise of a box of that duration in one transit (white and red),
    from the means of the light curve over bins of that width (half full
    at least), robustly

    :return: float [ppt]
    """
    means = []
    for seg in _segments(lc['time'], lc['sector']):
        time, flux = lc['time'][seg], lc['flux'][seg]
        if len(time) < 2:
            continue
        step = float(np.median(np.diff(time)))
        idx = np.floor((time - time[0]) / width).astype(int)
        count = np.bincount(idx)
        total = np.bincount(idx, flux)
        full = count >= max(1, 0.5 * width / max(step, 1e-6))
        means.append(total[full] / count[full])
    vals = np.concatenate(means) if means else np.array([])
    if len(vals) < 5:
        return float(np.std(lc['flux']) / max(1.0, math.sqrt(
            len(lc['flux']))))
    return float(1.4826 * np.median(np.abs(vals - np.median(vals))))


def _boxes(hp: Dict[str, np.ndarray], period: float, ref: float,
           half: float, dur: float) -> Optional[Dict[str, Any]]:
    """
    The best box of a high-passed light curve folded at a period: centres
    within half [phase] of ref, durations DURATIONS times dur [days]; the
    depth of a box minus its mean flux, its error the scatter of the light
    curve over its duration over the root of its transits

    :return: dict (centre, phase [h from ref], duration [h], depth,
             depth_err [ppt], snr, ntransits), or None
    """
    time, flux = hp['time'], hp['flux']
    phase = ((time - ref) / period + 0.5) % 1.0 - 0.5
    # the boxes: their centres a quarter of the shortest apart, on a grid
    #   of the folded light curve
    step = DURATIONS[0] * dur / 4.0 / period
    nbin = max(int(math.ceil(1.0 / step)), 8)
    idx = np.floor((phase + 0.5) * nbin).astype(int) % nbin
    count = np.bincount(idx, minlength=nbin).astype(float)
    total = np.bincount(idx, flux, minlength=nbin)
    cadence = float(np.median(np.diff(time)))
    best = None
    for scale in DURATIONS:
        width = scale * dur
        nwide = max(1, int(round(width / period * nbin)))
        # the sums over each box (wrapped about the phase)
        kern = np.ones(nwide)
        cnt = np.convolve(np.concatenate([count, count[:nwide - 1]]), kern,
                          'valid')[:nbin]
        tot = np.convolve(np.concatenate([total, total[:nwide - 1]]), kern,
                          'valid')[:nbin]
        mid = (np.arange(nbin) + 0.5 * nwide) / nbin - 0.5
        mid = (mid + 0.5) % 1.0 - 0.5
        ok = (cnt >= max(3, 0.5 * width / max(cadence, 1e-6))) & (
            np.abs(mid) <= half + 0.5 / nbin)
        if not np.any(ok):
            continue
        depth = np.where(ok, -tot / np.maximum(cnt, 1), -np.inf)
        noise = red_noise(hp, width)
        null = phase_scatter(hp, period, width)
        for k in np.argsort(depth)[::-1][:5]:
            if not ok[k] or depth[k] <= 0:
                break
            centre = ref + mid[k] * period
            # the transits in the box: the epochs with half of a box of
            #   points or more
            sel = np.abs(((time - centre) / period + 0.5) % 1.0 - 0.5) \
                <= 0.5 * width / period
            ep = np.round((time[sel] - centre) / period).astype(int)
            _, num = np.unique(ep, return_counts=True)
            ntr = int(np.sum(num >= 0.5 * width / max(cadence, 1e-6)))
            if ntr < 1:
                continue
            err = max(noise / math.sqrt(ntr), null)
            snr = float(depth[k] / err) if err > 0 else 0.0
            if best is None or snr > best['snr']:
                best = dict(centre=float(centre), phase=float(
                    mid[k] * period * 24.0), duration=float(width * 24.0),
                    depth=float(depth[k]), depth_err=float(err), snr=snr,
                    ntransits=ntr)
            break
    return best


def phase_scatter(hp: Dict[str, np.ndarray], period: float,
                  width: float) -> float:
    """
    The scatter of the depths of a box of a width [days] over every phase
    of the light curve folded at a period: what a box of that width gives
    where there is no transit, the null its depth is set against (robustly:
    a transit, a few phases, does not count). On a quiet star it is the
    noise of the box over its transits; on an active one, whose residuals
    have more than the red noise of a bin in them, more

    :return: float [ppt]
    """
    time, flux = hp['time'], hp['flux']
    nbin = max(int(math.ceil(4.0 * period / width)), 8)
    idx = np.floor(((time / period) % 1.0) * nbin).astype(int) % nbin
    count = np.bincount(idx, minlength=nbin).astype(float)
    total = np.bincount(idx, flux, minlength=nbin)
    kern = np.ones(4)
    cnt = np.convolve(np.concatenate([count, count[:3]]), kern,
                      'valid')[:nbin]
    tot = np.convolve(np.concatenate([total, total[:3]]), kern,
                      'valid')[:nbin]
    cadence = float(np.median(np.diff(time)))
    ok = cnt >= max(3, 0.5 * width / max(cadence, 1e-6))
    if ok.sum() < 8:
        return 0.0
    depth = -tot[ok] / cnt[ok]
    return float(1.4826 * np.median(np.abs(depth - np.median(depth))))


def _measure(hp: Dict[str, np.ndarray], period: float, centre: float,
             width: float) -> Optional[Dict[str, Any]]:
    """the depth, its error, its signal-to-noise ratio and its transits
    of one box (a centre, a width [days]) of a light curve folded at a
    period"""
    time, flux = hp['time'], hp['flux']
    sel = np.abs(((time - centre) / period + 0.5) % 1.0 - 0.5) * period \
        <= 0.5 * width
    if sel.sum() < 3:
        return None
    cadence = float(np.median(np.diff(time)))
    ep = np.round((time[sel] - centre) / period).astype(int)
    _, num = np.unique(ep, return_counts=True)
    ntr = int(np.sum(num >= 0.5 * width / max(cadence, 1e-6)))
    if ntr < 1:
        return None
    depth = -float(np.mean(flux[sel]))
    noise = red_noise(hp, width)
    null = phase_scatter(hp, period, width)
    err = max(noise / math.sqrt(ntr), null)
    # the depth of each transit (those with half a box of points), and the
    #   signal without the deepest: one event (a flare's dip, a
    #   systematic of one orbit) is not a transit seen twice
    full = [val for val, cnt in zip(*np.unique(ep, return_counts=True))
            if cnt >= 0.5 * width / max(cadence, 1e-6)]
    each = [-float(np.mean(flux[sel][ep == val])) for val in full]
    drop = None
    if len(each) >= 2:
        rest = sorted(each)[:-1]
        # the same error as the box's, for one transit fewer
        drop = float(np.mean(rest) / max(noise / math.sqrt(len(rest)),
                                         null * math.sqrt(len(each)
                                                          / len(rest))))
    return dict(depth=depth, depth_err=float(err),
                snr=float(depth / err) if err > 0 else 0.0, ntransits=ntr,
                epochs=each, snr_drop=drop)


def _detect(lc: Dict[str, np.ndarray], hp: Dict[str, np.ndarray],
            period: float, ref: float, half: float, dur: float,
            window: float):
    """the best box of a light curve at a period (its high-pass given),
    in two passes: searched, then the trend again without its transits
    and the same box measured again; (best or None, the high-pass)"""
    best = _boxes(hp, period, ref, half, dur)
    if best is None:
        return None, hp
    best.update(_measure(hp, period, best['centre'], best['duration']
                         / 24.0) or {})
    width = max(best['duration'] / 24.0, dur)
    near = np.abs(((lc['time'] - best['centre']) / period + 0.5) % 1.0
                  - 0.5) * period <= 0.75 * width
    hp = highpass(lc, window, mask=near)
    again = _measure(hp, period, best['centre'], best['duration'] / 24.0)
    if again is not None:
        best = dict(best, **again)
    return best, hp


def null_trials(lc: Dict[str, np.ndarray], period: float, half: float,
                mstar: float, rstar: float, ntrial: int = NTRIAL,
                avoid: Sequence[float] = (), seed: int = 2) -> List[float]:
    """
    The same search at periods where there is no transit to find: ntrial
    periods drawn between 0.7 and 1.4 times the period (none within 2 % of
    it, of a known planet's, or of their harmonics), the same width of
    window about a time drawn at random; the signal-to-noise ratio of the
    best box of each without its deepest transit (what the light curve
    gives by chance, on the statistic the candidate is judged on: one deep
    event is not a transit seen twice)

    :return: list of float
    """
    rng = np.random.default_rng(seed)
    tested = [period] + [float(val) for val in avoid if val]
    out = []
    tries = 0
    while len(out) < ntrial and tries < 20 * ntrial:
        tries += 1
        trial = float(period * np.exp(rng.uniform(np.log(0.7), np.log(1.4))))
        if any(abs(trial / (val * k) - 1) < 0.02 for val in tested
               for k in (0.5, 1.0, 2.0)):
            continue
        dur = max(duration(trial, mstar, rstar) / 24.0, 0.02)
        window = max(WINDOW_MIN, 3 * dur)
        hp = highpass(lc, window)
        ref = float(lc['time'].min() + rng.uniform(0, trial))
        best, _ = _detect(lc, hp, trial, ref, max(half, dur / trial), dur,
                          window)
        out.append(float(best.get('snr_drop') or 0.0) if best else 0.0)
    return out


def search(lc: Dict[str, np.ndarray], period: float,
           t0: Optional[float] = None, t0_err: Optional[float] = None,
           period_err: Optional[float] = None, mstar: float = 1.0,
           rstar: float = 1.0, npoints: int = 40000,
           seed: int = 1, others: Sequence[Dict[str, Any]] = (),
           ntrial: int = 20) -> Dict[str, Any]:
    """
    A transit in a light curve at a period (see the module)

    :param lc: dict, the light curve (from_csv, from_light_curves)
    :param period: float [days]
    :param t0: float or None, the expected time of a transit [BJD -
               2400000] (a conjunction of the velocities, or a TOI's)
    :param t0_err: float or None, its error [days]
    :param period_err: float or None, the error of the period [days]
    :param mstar: float, the mass of the star [solar masses]
    :param rstar: float, its radius [solar radii]
    :param npoints: int, the most points given back for the plot
    :param seed: int, the draw of those points
    :param others: list of dict, the transits of other planets (P, tc,
                   and duration [h] when known): their points left out,
                   so that they do not line up at another period
    :param ntrial: int, the periods of the null (null_trials), searched
                   when the box passes the other tests

    :return: dict: period, duration (expected, hours), window (the half
             width searched about t0, hours; None for the whole phase),
             best (centre: the time of a transit, phase: hours from t0 or
             the reference, duration [h], depth and depth_err [ppt], snr,
             ntransits, radius [Earth radii]), plausible (bool) and why,
             depth_earth and depth_jupiter [ppt], ntime, sectors, and the
             points to plot: hours (from the centre shown), flux; bins
             (hours, median, error)
    """
    dur_h = duration(period, mstar, rstar)
    dur = max(dur_h / 24.0, 0.02)
    window = max(WINDOW_MIN, 3 * dur)
    # the transits of the other planets out (their duration, or the one
    #   their period gives, and a margin)
    skip = np.zeros(len(lc['time']), bool)
    for other in others:
        if not other.get('P') or other.get('tc') is None:
            continue
        # the same planet (a TOI of a known planet): not another
        if abs(other['P'] / period - 1) < 0.01:
            continue
        width = (other.get('duration') or duration(other['P'], mstar,
                                                   rstar)) / 24.0
        skip |= np.abs(((lc['time'] - other['tc']) / other['P'] + 0.5) % 1.0
                       - 0.5) * other['P'] <= 0.75 * width + 0.02
    if skip.any():
        lc = {key: val[~skip] for key, val in lc.items()}
    hp = highpass(lc, window)
    out = dict(period=float(period), duration=float(dur_h), window=None,
               best=None, plausible=False, why='no light curve',
               depth_earth=depth_of(1.0, rstar),
               depth_jupiter=depth_of(R_JUPITER / R_EARTH, rstar),
               rstar=float(rstar), mstar=float(mstar), ntime=int(len(
                   hp['time'])), sectors=sorted({int(val) for val in
                                                 hp['sector']}),
               t0=t0, hours=[], flux=[], bins=[])
    if len(hp['time']) < 100:
        return out
    ref = float(t0) if t0 is not None else float(hp['time'].min())
    # where to search: within 3 sigma of the expected transit carried to
    #   the epochs of TESS (at least a duration), or everywhere
    half = 0.5
    if t0 is not None:
        epochs = float(np.median(np.abs(hp['time'] - ref))) / period
        sig = math.sqrt((t0_err or 0.0) ** 2 + (epochs * (period_err or 0.0))
                        ** 2)
        half = min(0.5, max(3 * sig, dur) / period)
        if half < 0.5:
            out['window'] = float(half * period * 24.0)
    # the best box, in two passes: searched, then the trend again without
    #   the transits found (each epoch's box, and a margin) and the same
    #   box measured again (not searched again: a second search on the
    #   trend interpolated across it would pick the noise twice)
    best, hp = _detect(lc, hp, period, ref, half, dur, window)
    time, flux = hp['time'], hp['flux']
    cadence = float(np.median(np.diff(time)))
    if best is not None:
        best['radius'] = radius_of(best['depth'], rstar)
        if t0 is not None:
            best['phase'] = float(((best['centre'] - t0) / period + 0.5)
                                  % 1.0 - 0.5) * period * 24.0
    out['best'] = best
    if best is None:
        out['why'] = 'no box with data where the transit would be'
    elif best['snr'] < SNR_PLAUSIBLE:
        out['why'] = (f'the deepest box is {best["snr"]:.1f} sigma deep '
                      f'(a transit: {SNR_PLAUSIBLE:.0f} or more)')
    elif best['ntransits'] < MIN_TRANSITS:
        out['why'] = (f'{best["snr"]:.1f} sigma, but in a single transit')
    elif (best.get('snr_drop') or 0.0) < 0.5 * SNR_PLAUSIBLE:
        out['why'] = (f'{best["snr"]:.1f} sigma, but from one transit: '
                      f'{best.get("snr_drop") or 0.0:.1f} sigma without it')
    elif best['radius'] > RP_MAX_RJ * R_JUPITER / R_EARTH:
        out['why'] = (f'{best["snr"]:.1f} sigma, but too deep for a planet '
                      f'(an eclipsing binary?)')
    else:
        # what the light curve gives by chance: the same search at other
        #   periods; a transit beats every one of them
        #   (the light curve without the box's own transits, which would
        #   line up by chance at some of them)
        width = max(best['duration'] / 24.0, dur)
        own = np.abs(((lc['time'] - best['centre']) / period + 0.5) % 1.0
                     - 0.5) * period <= width + 0.02
        null = null_trials({key: val[~own] for key, val in lc.items()},
                           period, half, mstar, rstar, ntrial,
                           [other.get('P') for other in others])
        out['null'] = null
        top = max(null) if null else 0.0
        if null and best['snr_drop'] <= top:
            out['why'] = (f'{best["snr"]:.1f} sigma in {best["ntransits"]} '
                          f'transits, but as strong a box at another period '
                          f'({best["snr_drop"]:.1f} sigma without its '
                          f'deepest transit, {top:.1f} at one of '
                          f'{len(null)}): this light curve makes dips of its '
                          f'own')
        else:
            out['plausible'] = True
            out['why'] = (f'{best["snr"]:.1f} sigma in {best["ntransits"]} '
                          f'transits ({best["snr_drop"]:.1f} without the '
                          f'deepest; at most {top:.1f} at {len(null)} other '
                          f'periods)')
    # what to plot: the hours from the transit found (plausible) or
    #   expected, the points (a draw of them when too many) and their
    #   medians in bins of a third of the duration
    centre = best['centre'] if (best and out['plausible']) else ref
    out['shown_centre'] = float(centre)
    hours = (((time - centre) / period + 0.5) % 1.0 - 0.5) * period * 24.0
    pick = np.arange(len(hours))
    if len(pick) > npoints:
        pick = np.sort(np.random.default_rng(seed).choice(
            len(pick), npoints, replace=False))
    out['hours'] = np.round(hours[pick], 4).tolist()
    out['flux'] = np.round(flux[pick], 4).tolist()
    width = max(dur_h / 3.0, cadence * 24.0 * 3)
    edges = np.arange(-0.5 * period * 24.0, 0.5 * period * 24.0 + width,
                      width)
    which = np.digitize(hours, edges)
    bins = []
    for k in np.unique(which):
        vals = flux[which == k]
        if len(vals) >= 5:
            med = float(np.median(vals))
            mad = 1.4826 * float(np.median(np.abs(vals - med)))
            bins.append([float(np.mean(hours[which == k])), med,
                         1.2533 * mad / math.sqrt(len(vals))])
    out['bins'] = bins
    return out


# =============================================================================
# End of code
# =============================================================================

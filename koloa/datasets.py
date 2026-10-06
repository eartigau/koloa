#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Which datasets of a series are used.

The archives of a star often hold the same spectra more than once (HARPS
through the ESO pipeline on DACE, through SERVAL in Trifonov et al. 2020,
in a paper with a reduction of its own), and old velocities too imprecise
to constrain anything. All of it stays on disk and can be used; two rules
say what is used by default, a dataset being the velocities of one
instrument from one source ('HARPS03', 'HARPS (Trifonov+ 2020)'):

1. The releases of the same spectra. A spectrum two datasets have (the
   same spectrograph, within SAME_FAMILY) is taken from the more precise
   of the two. On the spectra they share the star does the same in both,
   so the difference of the two is noise alone, and the covariance of
   each with that difference is its own noise: no model of the star is
   needed. What the star does only blurs the comparison (by chance it
   lines up a little with that difference), so its strongest signals are
   taken out of both first (star_signal: the same curve from both, which
   leaves their difference as it is). The more precise is the one with
   the smaller noise when the two differ by more than SIGNIFICANT times
   what chance gives, or the one with the smaller errors when they share
   fewer than MIN_COMPARE.
   When neither says: a dataset that names its spectrograph before a
   table that does not, the latest release, then the one with the more
   spectra. A table that holds the spectra of two spectrographs under one
   name (one offset for both) comes last. A dataset left with fewer than
   MIN_POINTS spectra of its own is left out.
2. The datasets that constrain nothing. A line is fitted to the nightly
   means (an offset per dataset, one slope). The error of a night is its
   own with, in quadrature, what the star does about a line (the scatter
   of the quietest well-sampled dataset beyond its errors: planets and
   activity, the same for every instrument, as a fit scaled to a reduced
   chi2 of 1 would have it) and what the dataset scatters beyond that
   (errors given too small). Taking a dataset away makes the error of the
   mean and that of the slope larger: one that changes both by less than
   WEAK is left out, the weakest first, until none is left.

The datasets of a file given are never left out (protect), and a dataset
asked back (include) is preferred to the other releases of its spectra.

    from koloa import datasets
    used, rows = datasets.choose(data, sources={'HARPS03': 'DACE', ...})
    for row in rows:
        print(datasets.told(row))

Created on 2026-10-06

@author: artigau
"""
import re
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from koloa import periodogram
from koloa.data import RVData, robust_std
from koloa.published import MIN_POINTS, SAME_ANY, SAME_FAMILY, family

# =============================================================================
# Define variables
# =============================================================================
#: the fewest spectra two datasets share for their noise there to say
#: which is the more precise (below, their errors say it)
MIN_COMPARE = 5
#: two median errors this close are the same: the dataset with the more
#: spectra is then the one preferred
TIE = 0.05
#: the noises of two releases differ when their difference is this many
#: times its own error (the same precision otherwise: the one with the
#: more spectra is then preferred). One sigma is enough: the less precise
#: release keeps the spectra it alone has, so taking one for the other
#: when they are as good costs nothing
SIGNIFICANT = 1.0
#: a pair of velocities this far from the others (in robust sigma of the
#: differences) is not in the comparison of two releases
CLIP = 5.0
#: a dataset that names no spectrograph holds spectra of another when,
#: of its velocities within NEARBY of one of the other's, this fraction
#: are within SAME_FAMILY (HARPS-TERRA velocities in a paper's table,
#: beside HARPS on DACE: all of them; two spectrographs that looked at
#: the star the same night: a few, by chance)
MOSTLY = 0.5
#: two velocities this close in time are of the same hour [days]
NEARBY = 1.0 / 24
#: the sinusoids taken out of the velocities before two releases are
#: compared, at most; only those that chance would give less than once in
#: FALSE_ALARM series (a sinusoid fitted to the noise of a release would
#: take that noise out of it, and make it look the more precise)
NSIGNALS = 6
FALSE_ALARM = 1e-3
#: their periods [days]: from this to twice the span of the series, the
#: frequencies OVERSAMPLE to a peak's width (MAXFREQ of them at most)
PMIN_SIGNAL = 1.5
OVERSAMPLE = 4
MAXFREQ = 30000
#: a dataset is left out when taking it away changes the error of the mean
#: and that of the slope by less than this
WEAK = 0.01
#: the fewest points for a scatter about a line (about their median below)
MIN_SCATTER = 5
#: the fewest nights of a dataset for its scatter to say what the star
#: does (the best sampled one when none has as many)
MIN_STAR = 20

#: the year of the releases that name none (Ribas et al. 2023)
ARCHIVE_YEARS = {'CARMENES DR1': 2023}

#: how close in time each pair of datasets that may share spectra is to be
Links = Dict[Tuple[str, str], float]


# =============================================================================
# Define functions
# =============================================================================
def about_line(time: np.ndarray, rv: np.ndarray
               ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Velocities with a line taken out (their median for fewer than
    MIN_SCATTER, or at one time), the line fitted without its outliers

    :return: tuple, what is left of each velocity, and which are not
             outliers (within 4 robust sigma)
    """
    time, rv = np.asarray(time, dtype=float), np.asarray(rv, dtype=float)
    keep = np.ones(len(rv), dtype=bool)
    if len(rv) < 2:
        return rv - np.sum(rv), keep
    rest = rv - np.median(rv)
    if len(rv) < MIN_SCATTER or np.ptp(time) == 0:
        return rest, keep
    mid = float(np.mean(time))
    for _ in range(3):
        coef = np.polyfit(time[keep] - mid, rv[keep], 1)
        rest = rv - np.polyval(coef, time - mid)
        sig = robust_std(rest)
        if not sig > 0:
            break
        again = np.abs(rest - np.median(rest)) < 4 * sig
        if np.sum(again) < MIN_SCATTER or np.all(again == keep):
            break
        keep = again
    return rest, keep


def line_scatter(time: np.ndarray, rv: np.ndarray) -> float:
    """
    The scatter of velocities about a line (about their median for fewer
    than MIN_SCATTER, or at one time): the standard deviation of what is
    left (variances add: the star's and the noise's), the outliers out

    :return: float [the unit of rv], nan for fewer than two velocities
    """
    rest, keep = about_line(time, rv)
    if len(rest) < 2:
        return np.nan
    return float(np.std(rest[keep]))


#: what the star does, for the series seen last (asked again at each tick
#: of the page)
_SIGNALS: Dict[Tuple[int, float, float], Callable] = {}


def star_signal(data: RVData, nsin: int = NSIGNALS
                ) -> Callable[[np.ndarray], np.ndarray]:
    """
    What the star does, as a curve in time, for the comparison of two
    releases of the same spectra: a line and the strongest sinusoids of
    the series night by night (the median of every velocity of a night,
    each dataset about its own median: one value a night, however many
    releases hold its spectra, or their common noise would count as many
    times and pass for a signal), found one after the other in a
    periodogram of what is left and fitted together. Not an analysis of
    the star: any curve taken out of both releases leaves their
    difference as it is, and the less of the star is left in them the
    better their noises are told apart

    :param data: RVData, every dataset of the star
    :param nsin: int, the sinusoids, at most

    :return: callable, the curve: times [days] -> velocities [m/s]
    """
    key = (int(data.n), float(np.round(np.sum(data.time), 3)),
           float(np.round(np.sum(data.rv), 3)))
    if key in _SIGNALS:
        return _SIGNALS[key]
    order = np.argsort(data.time, kind='stable')
    night = np.concatenate([[0], np.cumsum(np.diff(data.time[order]) > 0.5)])
    time, value, err = [], [], []
    for it in range(int(night[-1]) + 1 if len(night) else 0):
        sel = order[night == it]
        time.append(float(np.mean(data.time[sel])))
        value.append(float(np.median(data.rv[sel])))
        err.append(float(np.median(data.err[sel])))
    time, value, err = np.array(time), np.array(value), np.array(err)
    # the precise nights weigh the same, the imprecise ones less
    if len(err):
        err = np.sqrt(err ** 2 + np.median(err) ** 2)
    mid = float(np.mean(time))
    span = float(np.ptp(time))
    freqs: List[float] = []

    def design(when):
        cols = [np.ones(len(when)), when - mid]
        for freq in freqs:
            cols += [np.sin(2 * np.pi * freq * (when - mid)),
                     np.cos(2 * np.pi * freq * (when - mid))]
        return np.array(cols).T

    def solve():
        matrix = design(time) / err[:, None]
        return np.linalg.lstsq(matrix, value / err, rcond=None)[0]
    coef = solve() if len(time) > 2 else np.zeros(2)
    if len(time) >= 4 * MIN_SCATTER and span > 2 * PMIN_SIGNAL:
        num = int(min(MAXFREQ, OVERSAMPLE * span / PMIN_SIGNAL))
        grid = np.linspace(0.5 / span, 1.0 / PMIN_SIGNAL, num)
        left = value - design(time) @ coef
        for _ in range(nsin):
            if len(time) < 2 * len(coef) + 6:
                break
            power = periodogram.gls(time, left, err, grid)
            freqs.append(float(grid[int(np.argmax(power))]))
            # between two frequencies of the grid: the best of the fits of
            #   every sinusoid so far, finer each time, the new one then
            #   the others again (a strong signal a little off its period
            #   is left half in, and found a second time)
            step = float(grid[1] - grid[0])
            for which in [len(freqs) - 1] + list(range(len(freqs) - 1)):
                width = step if which == len(freqs) - 1 else step / 10.0
                for _ in range(3):
                    centre = freqs[which]
                    trials = np.linspace(centre - width, centre + width, 21)
                    chi2 = []
                    for trial in trials:
                        freqs[which] = float(trial)
                        chi2.append(float(np.sum(((
                            value - design(time) @ solve()) / err) ** 2)))
                    freqs[which] = float(trials[int(np.argmin(chi2))])
                    width /= 10.0
            again = solve()
            rest = value - design(time) @ again
            # what it takes out against what chance takes out at the best
            #   of that many frequencies: (chi2 before - after) / after,
            #   times half the degrees of freedom left, is exp(-x)
            #   distributed for noise
            before = float(np.sum((left / err) ** 2))
            after = float(np.sum((rest / err) ** 2))
            gain = (before - after) / max(after, 1e-30) \
                * (len(time) - len(again)) / 2.0
            if gain < np.log(num / OVERSAMPLE / FALSE_ALARM):
                freqs.pop()
                break
            coef, left = again, rest
    used = list(freqs)

    def curve(when):
        when = np.asarray(when, dtype=float)
        cols = [np.ones(len(when)), when - mid]
        for freq in used:
            cols += [np.sin(2 * np.pi * freq * (when - mid)),
                     np.cos(2 * np.pi * freq * (when - mid))]
        return np.array(cols).T @ coef
    curve.periods = [1.0 / freq for freq in used]
    if len(_SIGNALS) >= 8:
        _SIGNALS.pop(next(iter(_SIGNALS)))
    _SIGNALS[key] = curve
    return curve


def arm(name: Any) -> str:
    """the arm of a spectrograph in a dataset's name (VIS, NIR), '' when
    none is named"""
    found = re.search(r'(?<![A-Z])(VIS|NIR)(?![A-Z])', str(name).upper())
    return found.group(1) if found else ''


def nearest(times: np.ndarray, others: np.ndarray
            ) -> Tuple[np.ndarray, np.ndarray]:
    """the distance of each time to the nearest of others, and which"""
    order = np.argsort(others, kind='stable')
    srt = others[order]
    pos = np.searchsorted(srt, times)
    low = np.clip(pos - 1, 0, len(srt) - 1)
    high = np.clip(pos, 0, len(srt) - 1)
    dlow, dhigh = np.abs(times - srt[low]), np.abs(times - srt[high])
    return np.minimum(dlow, dhigh), order[np.where(dlow <= dhigh, low, high)]


def tolerance(data: RVData, one: str, other: str) -> float:
    """
    How close in time a spectrum of a dataset is to the same one in
    another: SAME_FAMILY for one spectrograph (and arm), 0 for two
    spectrographs (HARPS and NIRPS observe at the same time: not the same
    spectra). A dataset that names no spectrograph (a paper's table, of
    one spectrograph or of several) holds spectra of the other when, of
    its velocities within NEARBY of one of the other's, MOSTLY are within
    SAME_FAMILY (SAME_FAMILY then: the times of two releases differ by a
    minute or so, UTC for one and TDB for the other); it shares with it
    only what is within SAME_ANY otherwise

    :return: float [days]
    """
    if arm(one) and arm(other) and arm(one) != arm(other):
        return 0.0
    fam_one, fam_other = family(one), family(other)
    if fam_one and fam_other:
        return SAME_FAMILY if fam_one == fam_other else 0.0
    first, second = data.time[data.inst == one], data.time[data.inst == other]
    if len(first) > len(second):
        first, second = second, first
    dist = nearest(first, second)[0]
    close = int(np.sum(dist < SAME_FAMILY))
    if close >= MIN_COMPARE and close >= MOSTLY * np.sum(dist < NEARBY):
        return SAME_FAMILY
    return SAME_ANY


def links(data: RVData, sources: Optional[Dict[str, str]] = None) -> Links:
    """
    The pairs of datasets that may hold the same spectra, and how close
    in time two of them are then (tolerance): not two datasets of one
    source (the eras of an instrument, its arms), nor two spectrographs

    :param sources: dict or None, where each dataset came from (each its
                    own source when not said)

    :return: dict, {(one, other): days}, in the order of the series
    """
    sources = sources or {}
    names = data.instruments
    out: Links = {}
    for it, one in enumerate(names):
        for other in names[it + 1:]:
            if sources.get(one, one) == sources.get(other, other):
                continue
            tol = tolerance(data, one, other)
            if tol > 0:
                out[(one, other)] = tol
    return out


def _tol(tols: Links, one: str, other: str) -> float:
    """how close the spectra two datasets share are (0: they share none)"""
    return tols.get((one, other), tols.get((other, one), 0.0))


def compare(data: RVData, one: str, other: str, tol: float = SAME_FAMILY,
            signal: Optional[Callable] = None) -> Optional[Dict[str, Any]]:
    """
    Two datasets on the spectra they share: how many, and which is the
    more precise there

    With a = star + noise_a and b = star + noise_b on the same spectra,
    a - b is noise alone: cov(a, a - b) is the variance of noise_a and
    cov(b, b - a) that of noise_b, and var(a) - var(b) their difference,
    known to sqrt(var(a + b) var(a - b) / n).

    A release that gives one velocity for several spectra (those of a
    night, binned) is compared with the mean of the other's velocities of
    those spectra: bin against bin, not a bin against one exposure.

    :param tol: float, how close in time the same spectrum is [days]
    :param signal: callable or None, what the star does (star_signal),
                   taken out of both before they are compared

    :return: dict or None (no spectrum shared): n (the spectra of one that
             other has), ratio (the noise of one over that of other: below
             1, one is the more precise; 1 when they do not differ), by
             (noise or errors), the two noises, or median errors (one,
             other) [m/s], and extra (what the less precise scatters more
             than the other, in quadrature [m/s]; None by errors)
    """
    first, second = data.inst == one, data.inst == other
    dist, which = nearest(data.time[first], data.time[second])
    same = dist < tol
    if not np.any(same):
        return None
    num = int(np.sum(same))
    # from the one with the more velocities to the other: its velocities
    #   that are one velocity of the other are meant (finer, coarser)
    swap = int(np.sum(first)) < int(np.sum(second))
    fine, coarse = (second, first) if swap else (first, second)
    if swap:
        dist, which = nearest(data.time[fine], data.time[coarse])
        same = dist < tol
    pick = which[same]
    out = dict(n=num, by='errors', ratio=1.0, extra=None)
    efine = float(np.median(data.err[fine][same]))
    ecoarse = float(np.median(data.err[coarse][pick]))
    out['one'], out['other'] = (ecoarse, efine) if swap else (efine, ecoarse)
    groups, member = np.unique(pick, return_inverse=True)
    if len(groups) >= MIN_COMPARE:
        tfine, vfine = data.time[fine][same], data.rv[fine][same]
        if signal is not None:
            vfine = vfine - signal(tfine)
        size = np.bincount(member)
        tmean = np.bincount(member, weights=tfine) / size
        vmean = np.bincount(member, weights=vfine) / size
        tcoarse, vcoarse = data.time[coarse][groups], data.rv[coarse][groups]
        if signal is not None:
            vcoarse = vcoarse - signal(tcoarse)
        mine, _ = about_line(tmean, vmean)
        theirs, _ = about_line(tcoarse, vcoarse)
        diff = mine - theirs
        sig = robust_std(diff)
        if sig > 0:
            good = np.abs(diff - np.median(diff)) < CLIP * sig
            mine, theirs, diff = mine[good], theirs[good], diff[good]
        # the same velocities in both (one table copied into another):
        #   no noise to tell them by
        level = max(float(np.std(mine)), float(np.std(theirs)))
        if len(diff) >= MIN_COMPARE and np.std(diff) > 1e-6 * level:
            mine, theirs = mine - np.mean(mine), theirs - np.mean(theirs)
            diff = mine - theirs
            var_fine = float(np.mean(mine * diff))
            var_coarse = float(-np.mean(theirs * diff))
            gap = float(np.var(mine) - np.var(theirs))
            error = float(np.sqrt(np.var(mine + theirs) * np.var(diff)
                                  / len(diff)))
            ratio = 1.0
            if abs(gap) > SIGNIFICANT * error:
                # a noise that chance made negative, or smaller than can
                #   be told: no smaller than what the comparison resolves
                low = max(error, 1e-12)
                ratio = float(np.sqrt(max(var_fine, low)
                                      / max(var_coarse, low)))
            nfine = float(np.sqrt(max(var_fine, 0.0)))
            ncoarse = float(np.sqrt(max(var_coarse, 0.0)))
            out.update(by='noise', extra=float(np.sqrt(abs(gap))),
                       one=ncoarse if swap else nfine,
                       other=nfine if swap else ncoarse,
                       ratio=1.0 / ratio if swap else ratio)
            return out
    if out['other'] > 0 and abs(out['one'] / out['other'] - 1.0) >= TIE:
        out['ratio'] = out['one'] / out['other']
    return out


def year(name: Any, source: Any = '') -> int:
    """the year of a release: in the dataset's name ('HIRES (Teklu+
    2025)'), else in where it came from ('HIRES (CLS)' of Rosenthal et al.
    2021), else that of an archive known here (ARCHIVE_YEARS); 9999 for a
    living archive (DACE: as it is now)"""
    for text in (name, source):
        found = re.findall(r'(?<!\d)((?:19|20)\d{2})(?!\d)', str(text or ''))
        if found:
            return int(found[-1])
    return ARCHIVE_YEARS.get(str(source or ''), 9999)


def ratings(names: Sequence[str],
            pairs: Dict[Tuple[str, str], Dict[str, Any]]) -> Dict[str, float]:
    """
    One number per dataset from the comparisons of the pairs that share
    spectra, the smaller the more precise: the log of a noise, up to a
    constant, such that the difference of two of them is the log of the
    ratio measured for that pair (least squares, each pair weighted by the
    spectra it shares). A dataset is not better for having been compared
    with a poor one only: every comparison constrains the same scale

    :param names: list of str, the datasets
    :param pairs: dict, {(one, other): compare()}

    :return: dict, the rating of each dataset (0 for one compared with
             none), rounded so that equal ones are equal
    """
    index = {name: it for it, name in enumerate(names)}
    out = {name: 0.0 for name in names}
    if not pairs:
        return out
    rows, values, weights = [], [], []
    for (one, other), found in pairs.items():
        row = np.zeros(len(names))
        row[index[one]], row[index[other]] = 1.0, -1.0
        rows.append(row)
        values.append(np.log(found['ratio']))
        weights.append(np.sqrt(found['n']))
    # the scale is free: the mean of the ratings is zero (and that of each
    #   group of datasets compared with each other: the least-norm answer)
    design = np.array(rows) * np.array(weights)[:, None]
    solved = np.linalg.lstsq(design, np.array(values) * np.array(weights),
                             rcond=None)[0]
    for name in names:
        out[name] = float(np.round(solved[index[name]], 6)) + 0.0
    return out


def _upper(names: Optional[Sequence[str]]) -> List[str]:
    """names as they are compared: upper case, no space at their ends"""
    return [str(name).strip().upper() for name in (names or [])]


def ranking(data: RVData, tols: Links, protect: Sequence[str] = (),
            prefer: Sequence[str] = (), signal: Optional[Callable] = None,
            sources: Optional[Dict[str, str]] = None
            ) -> Tuple[List[str], Dict[Tuple[str, str], Dict[str, Any]]]:
    """
    The datasets of a series from the one whose spectra are kept first to
    the last: those of the file (protect, in their order), those asked for
    (prefer), then the more precise before the less precise on the spectra
    they share; for the same precision, a dataset that names its
    spectrograph before a table that names none, the latest release, the
    more spectra. Last, a table that names no spectrograph and holds the
    spectra of two (HIRES and HARPS under one name, one offset for both)

    :param tols: dict, the pairs that may share spectra (links)
    :param protect: list of str, the datasets of the file given
    :param prefer: list of str, the datasets asked back
    :param signal: callable or None, what the star does (star_signal)
    :param sources: dict or None, where each dataset came from (the year
                    of a release that does not name it)

    :return: tuple, the names in order, and the comparison of each pair
             that shares spectra ({(one, other): compare()})
    """
    names = data.instruments
    pairs: Dict[Tuple[str, str], Dict[str, Any]] = {}
    # the spectrographs whose spectra a table that names none holds
    holds: Dict[str, set] = {name: set() for name in names}
    for (one, other), tol in tols.items():
        found = compare(data, one, other, tol, signal)
        if found is None:
            continue
        pairs[(one, other)] = found
        if found['n'] >= MIN_COMPARE:
            for mine, theirs in ((one, other), (other, one)):
                if not family(mine) and family(theirs):
                    holds[mine].add(family(theirs))
    score = ratings(names, pairs)
    count = {name: int(np.sum(data.inst == name)) for name in names}
    first, second = _upper(protect), _upper(prefer)

    def key(name):
        upper = name.upper()
        if upper in first:
            return (0, first.index(upper))
        group = 1 if upper in second else 3 if len(holds[name]) > 1 else 2
        return (group, score[name], 0 if family(name) else 1,
                -year(name, (sources or {}).get(name)), -count[name], name)
    return sorted(names, key=key), pairs


def own_spectra(data: RVData, order: Sequence[str], tols: Links
                ) -> Tuple[Dict[str, np.ndarray], Dict[str, Dict[str, int]]]:
    """
    The spectra of each dataset that no dataset before it has, the
    datasets taken in order (a dataset left with fewer than MIN_POINTS of
    its own keeps none: its spectra are then the next one's to keep)

    :param order: list of str, the datasets used, the preferred first
    :param tols: dict, the pairs that may share spectra (links)

    :return: tuple, per dataset the points of the series it keeps
             (indices), and how many of its spectra each dataset before it
             has ({dataset: {other: n}})
    """
    kept: Dict[str, np.ndarray] = {}
    lost: Dict[str, Dict[str, int]] = {}
    for name in order:
        index = np.where(data.inst == name)[0]
        fresh = np.ones(len(index), dtype=bool)
        lost[name] = {}
        for other, theirs in kept.items():
            tol = _tol(tols, name, other)
            if not len(theirs) or not tol > 0:
                continue
            dist, _ = nearest(data.time[index], data.time[theirs])
            same = (dist < tol) & fresh
            if np.any(same):
                lost[name][other] = int(np.sum(same))
                fresh &= ~same
        if lost[name] and np.sum(fresh) < MIN_POINTS:
            fresh[:] = False
        kept[name] = index[fresh]
    return kept, lost


def budget(parts: Dict[str, Tuple[np.ndarray, np.ndarray]]
           ) -> Dict[str, Dict[str, float]]:
    """
    What each dataset is worth to a line through all of them (an offset
    per dataset, one slope): by how much the errors of the mean and of
    the slope grow when it is taken away

    The offset of a dataset costs it one point: n points weigh as n - 1 in
    the mean, and only their spread in time about their own mean weighs in
    the slope (what an old dataset far from the others would add to the
    slope goes into its offset).

    :param parts: dict, per dataset its times [days] and its errors [m/s]

    :return: dict, per dataset: mean and slope (the fractions by which the
             two errors grow without it; inf for the only one that
             constrains them)
    """
    weight, lever = {}, {}
    for name, (time, err) in parts.items():
        time = np.asarray(time, dtype=float)
        wgt = 1.0 / np.asarray(err, dtype=float) ** 2
        if not len(time):
            weight[name], lever[name] = 0.0, 0.0
            continue
        total = float(np.sum(wgt))
        centre = float(np.sum(wgt * time) / total)
        weight[name] = total * (len(time) - 1) / len(time)
        lever[name] = float(np.sum(wgt * (time - centre) ** 2))

    def grows(mine, total):
        if not total > 0 or not mine > 0:
            return 0.0
        left = total - mine
        return float(np.sqrt(total / left) - 1.0) if left > 0 else np.inf
    wsum, lsum = sum(weight.values()), sum(lever.values())
    return {name: dict(mean=grows(weight[name], wsum),
                       slope=grows(lever[name], lsum)) for name in parts}


def star_scatter(excess: Dict[str, Tuple[int, float]]) -> float:
    """
    What the star does about a line, from the scatter of each dataset
    beyond its errors: the smallest among the datasets of MIN_STAR nights
    or more (a dataset of a few nights can scatter little by chance), or
    that of the best sampled one when none has as many

    :param excess: dict, per dataset its nights and the variance of its
                   scatter beyond its errors [m2/s2]

    :return: float [m/s], 0 when no dataset says
    """
    if not excess:
        return 0.0
    many = [val for num, val in excess.values() if num >= MIN_STAR]
    if not many:
        many = [max(excess.values())[1]]
    return float(np.sqrt(min(many)))


def nightly_parts(data: RVData, kept: Dict[str, np.ndarray]
                  ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    The nights of each dataset and their errors, for the budget: the
    formal error of each nightly mean, with what the star does about a
    line (star_scatter: in every dataset) and what the dataset scatters
    beyond that (errors given too small)

    :return: dict, per dataset its times and its errors [m/s]
    """
    nights, excess = {}, {}
    for name, index in kept.items():
        if not len(index):
            continue
        part = data.select(index).nightly()
        nights[name] = part
        if part.n >= MIN_SCATTER:
            scatter = line_scatter(part.time, part.rv)
            excess[name] = max(0.0, scatter ** 2
                               - float(np.median(part.err ** 2)))
    star = star_scatter({name: (nights[name].n, val)
                         for name, val in excess.items()}) ** 2
    return {name: (part.time, np.sqrt(part.err ** 2 + star + max(
        0.0, excess.get(name, star) - star)))
            for name, part in nights.items()}


def rules(data: RVData, sources: Optional[Dict[str, str]] = None,
          protect: Sequence[str] = (), prefer: Sequence[str] = (),
          exclude: Sequence[str] = (), auto: bool = True
          ) -> List[Dict[str, Any]]:
    """
    What is used of each dataset of a series, and why

    :param data: RVData, every dataset of the star, whole
    :param sources: dict or None, where each dataset came from
    :param protect: list of str, the datasets of the file given: never
                    left out by the rules, their spectra kept first
    :param prefer: list of str, datasets asked back: not left out by the
                   rules, their spectra kept before the other releases'
    :param exclude: list of str, datasets left out as asked
    :param auto: bool, False for the rules to leave nothing out (a
                 spectrum several datasets have is still used once)

    :return: list of dict, per dataset in the order of the series: name,
             n (its points), used (those used), keep (their indices in
             the series), status (on; part: some of its spectra are a
             better release's; release: all of them; weak: it constrains
             nothing; asked: left out as asked), better and same (the
             release that has its spectra, how many), by (noise or
             errors), extra (what it scatters more than the better one,
             in quadrature [m/s]), precision and precision_better (the
             noise, or median error, of the two [m/s]), tie (the two were
             not told apart: the other is preferred), mean and slope
             (what it adds to the errors of the mean and of the slope),
             left (the points of its own of a weak one), nights and source
    """
    sources = sources or {}
    tols = links(data, sources)
    order, pairs = ranking(data, tols, protect, prefer,
                           star_signal(data) if tols else None, sources)
    asked = [name for name in order if name.upper() in _upper(exclude)]
    held = _upper(protect) + _upper(prefer)
    used = [name for name in order if name not in asked]
    weak: Dict[str, Dict[str, Any]] = {}
    while True:
        kept, lost = own_spectra(data, used, tols)
        active = {name: index for name, index in kept.items() if len(index)}
        parts = nightly_parts(data, active)
        worth = budget(parts)
        if not auto:
            break
        # the weakest dataset that constrains neither, taken away, and the
        #   others weighed again without it
        low = [name for name in active if name.upper() not in held
               and worth[name]['mean'] < WEAK and worth[name]['slope'] < WEAK]
        if not low or len(active) < 2:
            break
        last = min(low, key=lambda name: max(worth[name]['mean'],
                                             worth[name]['slope']))
        weak[last] = dict(worth[last], left=int(len(kept[last])),
                          lost=lost[last], nights=int(len(parts[last][0])))
        used.remove(last)
    rows = []
    for name in data.instruments:
        num = int(np.sum(data.inst == name))
        row = dict(name=name, n=num, used=0, keep=np.array([], dtype=int),
                   status='on', better=None, same=0, by=None, extra=None,
                   tie=False,
                   precision=None, precision_better=None, mean=None,
                   slope=None, left=None, nights=None,
                   source=sources.get(name, ''))
        gone: Dict[str, int] = {}
        if name in asked:
            row['status'] = 'asked'
        elif name in weak:
            gone = weak[name]['lost']
            row.update(status='weak', mean=weak[name]['mean'],
                       slope=weak[name]['slope'], left=weak[name]['left'],
                       nights=weak[name]['nights'])
        else:
            gone = lost[name]
            row.update(keep=kept[name], used=int(len(kept[name])),
                       **worth.get(name, {}))
            if name in parts:
                row['nights'] = int(len(parts[name][0]))
            if gone:
                row['status'] = 'part' if row['used'] else 'release'
        if gone:
            # the release that has the most of its spectra
            better = max(gone, key=gone.get)
            found = pairs.get((name, better))
            mine, theirs = 'one', 'other'
            if found is None:
                found = pairs.get((better, name))
                mine, theirs = 'other', 'one'
            row.update(better=better, same=sum(gone.values()))
            if found is not None:
                # a tie: the two could not be told apart, and the better
                #   one is so by its name, its year or its size
                row.update(by=found['by'], extra=found['extra'],
                           precision=found[mine],
                           precision_better=found[theirs],
                           tie=found['ratio'] == 1.0)
        rows.append(row)
    return rows


def choose(data: RVData, sources: Optional[Dict[str, str]] = None,
           protect: Sequence[str] = (), include: Sequence[str] = (),
           exclude: Sequence[str] = (), auto: bool = True
           ) -> Tuple[Optional[RVData], List[Dict[str, Any]]]:
    """
    The series used: its datasets left out by the rules or as asked taken
    out, and each spectrum taken from one dataset only

    :param data: RVData, every dataset of the star, whole
    :param include: list of str, datasets used though the rules would
                    leave them out (and before the other releases of
                    their spectra)
    :param exclude: list of str, datasets left out as asked

    :return: tuple, the series used (None when nothing is left) and what
             was done with each dataset (rules)
    """
    rows = rules(data, sources, protect, include, exclude, auto)
    keep = np.concatenate([row['keep'] for row in rows]).astype(int)
    if not len(keep):
        return None, rows
    if len(keep) == data.n:
        return data, rows
    return data.select(np.sort(keep)), rows


def left_out(rows: Sequence[Dict[str, Any]]) -> List[str]:
    """the datasets the rules leave out (not those left out as asked)"""
    return [row['name'] for row in rows
            if row['status'] in ('release', 'weak')]


def _better(row: Dict[str, Any]) -> str:
    """the release that has spectra of a dataset, and how the two compare"""
    if row.get('precision') is None:
        return str(row['better'])
    if row.get('tie'):
        return (f'{row["better"]} (as precise on the spectra they share: '
                f'the one that names its spectrograph, the latest, or the '
                f'one with the more spectra, is preferred)')
    if row['by'] == 'noise' and row.get('extra') is not None:
        if not row['precision'] > row['precision_better']:
            return (f'{row["better"]} (preferred, from the comparisons of '
                    f'all the releases of these spectra)')
        return (f'{row["better"]}, the more precise ({row["name"]} '
                f'scatters {row["extra"]:.2f} m/s more, in quadrature, on '
                f'the spectra they share)')
    return (f'{row["better"]} (median errors of '
            f'{row["precision_better"]:.2f} against {row["precision"]:.2f} '
            f'm/s on the few spectra they share)')


def _count(num: int, word: str) -> str:
    """a number and its noun: 1 point, 3 points"""
    return f'{num} {word}' + ('' if num == 1 else 's')


def told(row: Dict[str, Any]) -> str:
    """what was done with a dataset, in words"""
    name = row['name']
    points = _count(row['n'], 'point')
    if row['status'] == 'asked':
        return f'{name}: left out, as asked ({points})'
    if row['status'] == 'weak':
        words = (f'the mean ({100 * row["mean"]:.2f} % of its error) nor '
                 f'the slope ({100 * row["slope"]:.2f} %)')
        if row['better']:
            return (f'{name}: left out, {row["same"]} of its {points} are '
                    f'spectra of {_better(row)}, and the {row["left"]} '
                    f'others constrain neither {words}')
        if row['n'] == 1:
            return (f'{name}: left out, one point constrains neither '
                    f'{words}')
        return (f'{name}: left out, its {points} '
                f'({_count(row["nights"], "night")}) constrain neither '
                f'{words}')
    if row['status'] == 'release':
        return (f'{name}: left out, its {points} are spectra of '
                f'{_better(row)}')
    if row['status'] == 'part':
        return (f'{name}: {row["used"]} of {points}, the others are '
                f'spectra of {_better(row)}')
    return f'{name}: {points}'


# =============================================================================
# End of code
# =============================================================================

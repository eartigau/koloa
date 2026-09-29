#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Recovery rates in period and semi-amplitude: what a series could have found.

Circular planets are injected into a series over a grid of periods and
semi-amplitudes, at random phases, and each one is looked for by folding the
series at its period: an offset and a sinusoid are fitted there with the
outlier-aware likelihood of koloa's periodogram (each visit, or each
exposure, is good or an outlier, and the jitter is free; the EM of
koloa.periodogram.profile_periodogram), and the planet counts as recovered
when the sinusoid improves the likelihood by more than a threshold. The same
test with a gaussian likelihood is the reference.

Two tests:

- 'fold': the period is known (a transiting planet, a candidate to
  confirm). Under the null, 2 Delta ln L of a sinusoid at one frequency
  follows a chi^2 with two degrees of freedom, so a false-alarm probability
  fap sets the threshold Delta ln L > -ln(fap).
- 'search': the period is not known. The whole periodogram is computed; the
  planet counts when the highest peak lies within 1/T of the injected
  frequency and passes the threshold of the grid, Delta ln L >
  ln(M / fap), with M = T (f_max - f_min) independent frequencies.

Both thresholds are checked, not trusted: injections at K = 0 measure the
false-alarm rate of each test on the very series used. On a real series a
drift (trend=1 or 2 fits it in both models) or strong activity breaks the
chi^2 null, and that false-alarm rate is how one finds out.

The series can be the real one (the planets are added to it as observed,
its outliers and its own signals included) or a function of a seed that
returns a new realisation for every injection, which maps the method rather
than one draw of the noise.

Created on 2026-09-27

@author: artigau
"""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

import numpy as np

from koloa.data import RVData
from koloa.log import log
from koloa.periodogram import frequency_grid, oap

# =============================================================================
# Define variables
# =============================================================================
#: a series, or a function of a seed that makes one
Series = Union[RVData, Callable[[int], RVData]]


# =============================================================================
# The result
# =============================================================================
@dataclass
class RecoveryMap:
    """
    The recovered fraction on a grid of periods and semi-amplitudes
    """
    period_edges: np.ndarray
    amp_edges: np.ndarray
    rate: np.ndarray
    counts: np.ndarray
    injections: List[Dict[str, float]]
    test: str
    outliers: bool
    unit: str
    fap: float
    threshold: float
    false_alarm: float = float('nan')
    nnull: int = 0
    extra: Dict[str, Any] = field(default_factory=dict)

    @property
    def periods(self) -> np.ndarray:
        """The centres of the period bins (geometric) [days]"""
        return np.sqrt(self.period_edges[1:] * self.period_edges[:-1])

    @property
    def amplitudes(self) -> np.ndarray:
        """The centres of the amplitude bins (geometric) [m/s]"""
        return np.sqrt(self.amp_edges[1:] * self.amp_edges[:-1])

    def k_at(self, level: float = 0.5) -> np.ndarray:
        """
        The semi-amplitude recovered with a given probability, per period
        bin: the detection limit (interpolated in log K; nan when the grid
        never reaches that level)

        :param level: float, the recovered fraction
        :return: np.ndarray, (nperiod) [m/s]
        """
        out = np.full(len(self.periods), np.nan)
        logk = np.log(self.amplitudes)
        for ip, row in enumerate(self.rate):
            # the recovered fraction grows with K; smooth the steps a little
            row = np.maximum.accumulate(np.nan_to_num(row))
            above = np.where(row >= level)[0]
            if not len(above):
                continue
            it = above[0]
            if it == 0:
                out[ip] = self.amplitudes[0]
                continue
            lo, hi = row[it - 1], row[it]
            frac = (level - lo) / (hi - lo) if hi > lo else 1.0
            out[ip] = float(np.exp(logk[it - 1] + frac * (logk[it]
                                                          - logk[it - 1])))
        return out

    def summary(self) -> str:
        """A few lines on the map"""
        like = 'outlier-aware' if self.outliers else 'gaussian'
        lines = [f'Recovery map ({self.test} test, {like}, {self.unit} '
                 f'units): {int(np.sum(self.counts))} injections, threshold '
                 f'Delta ln L > {self.threshold:.2f} (fap {self.fap:g})']
        if self.nnull:
            lines.append(f'  false alarms at K = 0: {self.false_alarm:.3f} '
                         f'({self.nnull} injections)')
        k50, k90 = self.k_at(0.5), self.k_at(0.9)
        for per, low, high in zip(self.periods, k50, k90):
            lines.append(f'  P = {per:7.2f} d: K(50 %) = {low:6.2f} m/s, '
                         f'K(90 %) = {high:6.2f} m/s')
        return '\n'.join(lines)


# =============================================================================
# Define functions
# =============================================================================
def inject(data: RVData, period: float, amp: float, phase: float,
           tref: Optional[float] = None) -> RVData:
    """
    A copy of a series with a circular orbit added

    :param data: RVData, the series
    :param period: float, the period [days]
    :param amp: float, the semi-amplitude [m/s]
    :param phase: float, the phase at tref [rad]
    :param tref: float or None, the reference time (the series' when None)

    :return: RVData, the copy
    """
    tref = data.tref if tref is None else tref
    signal = amp * np.sin(2 * np.pi * (data.time - tref) / period + phase)
    return data.with_values(data.rv + signal)


def fold_test(data: RVData, period: float, outliers: bool = True,
              unit: str = 'sequence', niter: int = 30,
              trend: int = 0) -> Dict[str, float]:
    """
    The gain of a sinusoid at one period, outlier-aware or gaussian

    :param data: RVData, the series
    :param period: float, the period [days]
    :param outliers: bool, the mixture (True) or the gaussian (False)
    :param unit: str, sequence (visit means) or point
    :param niter: int, the EM iterations
    :param trend: int, a polynomial of this degree in both models

    :return: dict, dlnl (the gain in ln L) and amp (the fitted amplitude)
    """
    res = oap(data, np.array([1.0 / period]), unit=unit, outliers=outliers,
              niter=niter, trend=trend)
    return dict(dlnl=float(res['dlnl'][0]), amp=float(res['amp'][0]))


def search_test(data: RVData, freq: np.ndarray, period: float,
                outliers: bool = True, unit: str = 'sequence',
                niter: int = 30, trend: int = 0) -> Dict[str, float]:
    """
    A blind search: the highest peak of the periodogram, and whether it is
    at the injected period

    :param data: RVData, the series
    :param freq: np.ndarray, the frequencies searched [1/day]
    :param period: float, the injected period [days]
    :param outliers: bool, the mixture (True) or the gaussian (False)
    :param unit: str, sequence or point
    :param niter: int, the EM iterations

    :return: dict, dlnl (the highest peak), at (whether it is within 1/T
             of the injected frequency), best (its period), amp
    """
    res = oap(data, freq, unit=unit, outliers=outliers, niter=niter,
              trend=trend)
    best = int(np.argmax(res['dlnl']))
    near = abs(freq[best] - 1.0 / period) <= 1.0 / data.baseline
    return dict(dlnl=float(res['dlnl'][best]), at=bool(near),
                best=float(1.0 / freq[best]), amp=float(res['amp'][best]))


def _threshold(test: str, fap: float, data: RVData,
               freq: Optional[np.ndarray]) -> float:
    """The gain a planet must reach for a false-alarm probability fap"""
    if test == 'fold':
        return float(-np.log(fap))
    nindep = max(data.baseline * (np.max(freq) - np.min(freq)), 1.0)
    return float(np.log(nindep / fap))


def _run_chunk(job) -> List[Dict[str, float]]:
    """A chunk of injections (one worker), quiet; run in the calling
    process (one worker), it gives the log and the warnings back after"""
    import warnings
    from koloa import log as klog
    verbose = klog.VERBOSE
    klog.VERBOSE = False
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return _run_chunk_quiet(job)
    finally:
        klog.VERBOSE = verbose


def _run_chunk_quiet(job) -> List[Dict[str, float]]:
    """The injections of a chunk"""
    series, chunk, test, outliers, unit, freq, threshold, trend = job
    out = []
    for inj in chunk:
        data = series(inj['seed']) if callable(series) else series
        data = inject(data, inj['period'], inj['amp'], inj['phase'])
        if test == 'fold':
            res = fold_test(data, inj['period'], outliers=outliers, unit=unit,
                            trend=trend)
            found = res['dlnl'] > threshold
        else:
            res = search_test(data, freq, inj['period'], outliers=outliers,
                              unit=unit, trend=trend)
            found = res['at'] and res['dlnl'] > threshold
        out.append(dict(inj, dlnl=res['dlnl'], amp_fit=res['amp'],
                        at=bool(res.get('at', True)), recovered=bool(found)))
    return out


def recovery_map(data: Series, period_edges: Sequence[float],
                 amp_edges: Sequence[float], ninj: int = 20,
                 test: str = 'fold', outliers: bool = True,
                 unit: str = 'sequence', fap: float = 0.01,
                 nnull: int = 200, seed: int = 1, workers: int = 1,
                 freq: Optional[np.ndarray] = None, trend: int = 0,
                 calibrate: bool = False,
                 quiet: bool = False) -> RecoveryMap:
    """
    The fraction of injected planets recovered, on a grid of periods and
    semi-amplitudes

    Each cell receives ninj circular planets, with a period and an amplitude
    drawn log-uniformly inside it and a random phase. nnull more injections
    at K = 0 measure the false-alarm rate of the test.

    :param data: RVData, or a function of a seed that returns one (a new
                 realisation for every injection)
    :param period_edges: list of float, the edges of the period bins [days]
    :param amp_edges: list of float, the edges of the semi-amplitude bins
                      [m/s]
    :param ninj: int, injections per cell
    :param test: str, 'fold' (known period) or 'search' (blind)
    :param outliers: bool, the outlier-aware likelihood (True) or the
                     gaussian one (False)
    :param unit: str, 'sequence' (each visit a unit, fitted on the visit
                 means) or 'point'
    :param fap: float, the false-alarm probability that sets the threshold
    :param nnull: int, the injections at K = 0
    :param seed: int, the seed of the injections
    :param workers: int, the processes
    :param freq: np.ndarray or None, the grid of the search test (from the
                 shortest period edge to the baseline when None)
    :param trend: int, a polynomial of this degree beside the sinusoid and
                  in the null: a drift (a companion, the secular
                  acceleration) otherwise fires the test at long periods,
                  which the false alarms at K = 0 show
    :param calibrate: bool, set the threshold of each period bin from the
                      series itself: the (1 - fap) quantile of the gains
                      that nnull injections at K = 0 reach in that bin
                      (nnull per bin). Red noise (activity, systematics)
                      breaks the chi^2 null at long periods; this keeps
                      the false alarms at fap in every bin
    :param quiet: bool, no log lines

    :return: RecoveryMap
    """
    if test not in ('fold', 'search'):
        raise ValueError(f'Unknown test: {test}')
    period_edges = np.asarray(period_edges, dtype=float)
    amp_edges = np.asarray(amp_edges, dtype=float)
    first = data(seed) if callable(data) else data
    if test == 'search' and freq is None:
        freq = frequency_grid(first.time, float(period_edges[0]) * 0.9,
                              None, 5)
    threshold = _threshold(test, fap, first, freq)
    rng = np.random.default_rng(seed)
    injections = []
    for ip in range(len(period_edges) - 1):
        for ik in range(len(amp_edges) - 1):
            for _ in range(ninj):
                injections.append(dict(
                    ip=ip, ik=ik,
                    period=float(np.exp(rng.uniform(
                        np.log(period_edges[ip]),
                        np.log(period_edges[ip + 1])))),
                    amp=float(np.exp(rng.uniform(np.log(amp_edges[ik]),
                                                 np.log(amp_edges[ik + 1])))),
                    phase=float(rng.uniform(0, 2 * np.pi)),
                    seed=int(rng.integers(1, 2 ** 31))))
    if calibrate:
        # nnull empty injections in every period bin, for its threshold
        for ip in range(len(period_edges) - 1):
            for _ in range(nnull):
                injections.append(dict(
                    ip=-1, ik=-1, bin=ip,
                    period=float(np.exp(rng.uniform(
                        np.log(period_edges[ip]),
                        np.log(period_edges[ip + 1])))),
                    amp=0.0, phase=0.0, seed=int(rng.integers(1, 2 ** 31))))
    else:
        for _ in range(nnull):
            injections.append(dict(
                ip=-1, ik=-1, bin=-1,
                period=float(np.exp(rng.uniform(np.log(period_edges[0]),
                                                np.log(period_edges[-1])))),
                amp=0.0, phase=0.0, seed=int(rng.integers(1, 2 ** 31))))
    nchunk = max(workers * 8, 1)
    chunks = [injections[it::nchunk] for it in range(nchunk)]
    jobs = [(data, chunk, test, outliers, unit, freq, threshold, trend)
            for chunk in chunks if chunk]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = [row for part in pool.map(_run_chunk, jobs)
                       for row in part]
    else:
        results = [row for job in jobs for row in _run_chunk(job)]
    shape = (len(period_edges) - 1, len(amp_edges) - 1)
    nominal = threshold
    thresholds = np.full(shape[0], threshold)
    if calibrate:
        # each bin's threshold from its own empty injections, and every
        #   injection judged again against it
        for ip in range(shape[0]):
            gains = [row['dlnl'] for row in results
                     if row['ip'] < 0 and row['bin'] == ip]
            if gains:
                thresholds[ip] = max(float(np.quantile(gains, 1 - fap)),
                                     nominal)
        for row in results:
            ip = row['ip'] if row['ip'] >= 0 else row['bin']
            passed = row['dlnl'] > thresholds[ip]
            if test == 'search':
                passed = passed and row.get('at', True)
            row['recovered'] = bool(passed)
    found, counts = np.zeros(shape), np.zeros(shape)
    nulls = []
    for row in results:
        if row['ip'] < 0:
            nulls.append(row['recovered'])
            continue
        counts[row['ip'], row['ik']] += 1
        found[row['ip'], row['ik']] += row['recovered']
    rate = np.where(counts > 0, found / np.maximum(counts, 1), np.nan)
    out = RecoveryMap(period_edges=period_edges, amp_edges=amp_edges,
                      rate=rate, counts=counts, injections=results,
                      test=test, outliers=outliers, unit=unit, fap=fap,
                      threshold=threshold,
                      false_alarm=float(np.mean(nulls)) if nulls
                      else float('nan'), nnull=len(nulls),
                      extra=dict(nominal_threshold=nominal,
                                 thresholds=thresholds.tolist(),
                                 calibrated=calibrate, trend=trend,
                                 false_alarm_nominal=float(np.mean(
                                     [row['dlnl'] > nominal for row in results
                                      if row['ip'] < 0]))
                                 if nulls else float('nan')))
    if not quiet:
        like = 'outlier-aware' if outliers else 'gaussian'
        log(f'Recovery map ({test}, {like}): {len(results)} injections, '
            f'false alarms at K = 0: {out.false_alarm:.3f}', 'value')
    return out


# =============================================================================
# End of code
# =============================================================================

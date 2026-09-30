#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Is it a planet? The duck test.

"If it looks like a planet, swims like a planet and quacks like a planet,
it is probably a planet." A significant FIP only says that a periodic
signal is there; whether it is an orbit is a separate question, and the
velocities hold several ways of asking it. `duck_test` asks all of them and
says which ones the signal fails:

1. significance: the outlier-aware FIP of its interval, next to the
   gaussian one, the window power and the FIP of its aliases;
2. robustness: whether a handful of visits hold the peak up (jackknife);
3. coherence: whether the amplitude and phase stay the same from one part
   of the campaign to the next (an orbit does, a spotted star does not);
4. shape: whether an eccentric orbit is preferred only by running to the
   edge of the prior, which is what a non-keplerian shape does;
5. activity indicators: whether an indicator has power at the period, or
   at two or three times it (the velocity signal is then a harmonic of the
   rotation);
6. absorption: whether a GP with a period set by an indicator describes the
   velocities as well without the orbit as with it.

The verdict counts the flags. It is a summary for a human, not a
replacement for looking at the figures it comes with.

Created on 2026-09-27

@author: artigau
"""
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np
from scipy import stats

from koloa.data import RVData
from koloa.fip import FIPResult, fip_single
from koloa.linear import weighted_lstsq
from koloa.log import log
from koloa.periodogram import frequency_grid, gls, influence, jackknife
from koloa.periodogram import window

# =============================================================================
# Define variables
# =============================================================================
#: indicators tried, in order, for the activity period
ACTIVITY_INDICATORS = ['DTEMP3500', 'DTEMP', 'dW', 'fwhm', 'FWHM', 'bis',
                       'contrast', 'd2v']

#: columns that are velocities, not activity indicators
VELOCITY_PREFIXES = ('vrad', 'rv', 'vel')

#: how far from a multiple of the period an indicator peak still counts:
#:   rotation is not a sharp period (differential rotation, spot evolution)
HARMONIC_TOLERANCE = 0.1


# =============================================================================
# Individual checks
# =============================================================================
def seasons(time: np.ndarray, gap: float = 60.0) -> np.ndarray:
    """
    The observing season of each point: a new season after a long gap

    :param time: np.ndarray, the time [days]
    :param gap: float, the gap that starts a season [days]

    :return: np.ndarray, the season index of each point
    """
    order = np.argsort(time)
    new = np.concatenate([[0], np.diff(time[order]) > gap]).astype(int)
    out = np.zeros(len(time), dtype=int)
    out[order] = np.cumsum(new)
    return out


def _visit_means(data: RVData, rel: np.ndarray) -> RVData:
    """
    The weighted mean of each visit, each exposure weighted by its inverse
    variance times its reliability (1 - its outlier probability)

    :param data: RVData, the series
    :param rel: np.ndarray, the reliability of each exposure

    :return: RVData, one point per visit
    """
    weight = rel / data.err ** 2
    nseq = data.nseq
    wsum = np.bincount(data.seq, weights=weight, minlength=nseq)
    out = data.binned()
    out.time = np.bincount(data.seq, weights=weight * data.time,
                           minlength=nseq) / wsum
    out.rv = np.bincount(data.seq, weights=weight * data.rv,
                         minlength=nseq) / wsum
    out.err = 1 / np.sqrt(wsum)
    return out


def _excess_jitter(design: np.ndarray, value: np.ndarray,
                   err: np.ndarray) -> float:
    """
    The jitter that brings the reduced chi2 of a linear fit to one

    :param design: np.ndarray, the design matrix (n, p)
    :param value: np.ndarray, the data
    :param err: np.ndarray, the errors

    :return: float, the jitter (0 when the chi2 is already below one)
    """
    npts, npar = design.shape
    dof = npts - npar
    if dof < 1:
        return 0.0

    def excess(jit):
        tot = np.sqrt(err ** 2 + jit ** 2)
        _, _, chi2 = weighted_lstsq(design, value, tot)
        return chi2 / dof - 1.0
    if excess(0.0) <= 0:
        return 0.0
    upper = 10 * float(np.std(value)) + 1e-6
    if excess(upper) > 0:
        return upper
    from scipy.optimize import brentq
    return float(brentq(excess, 0.0, upper, xtol=1e-4 * upper))


def coherence(data: RVData, period: float,
              prob: Optional[np.ndarray] = None,
              jitter: Optional[float] = None,
              split: str = 'halves', min_points: int = 12,
              trend: int = 1, visits: bool = True) -> Dict[str, Any]:
    """
    Amplitude and phase of a sinusoid at a fixed period, chunk by chunk

    Each chunk has its own offset; the sinusoid (A cos + B sin) is fitted
    by weighted least squares, the weight of each point being its inverse
    variance times its reliability. The chunks are then compared: the
    amplitudes K by a chi2 around their weighted mean, and the vectors
    (A, B) by a chi2 around theirs, with their covariances.

    The exposures of a visit are not independent, so by default the test
    runs on the visit means (each exposure weighted by its reliability),
    and the noise that the error bars do not hold (a jitter, activity) is
    added to them: the jitter that brings the reduced chi2 of the global
    fit (the sinusoid, an offset per chunk and the trend) to one. Without
    both, the chi2 tests count the jitter as signal, and their p-values
    are far too small.

    :param data: RVData, the series
    :param period: float, the period [days]
    :param prob: np.ndarray or None, the outlier probability of each point
    :param jitter: float or None, a jitter added to the errors [m/s]; None
                   measures it from the global fit
    :param split: str, 'halves' (two halves in time, by number of visits),
                  'seasons', or 'thirds'
    :param min_points: int, the smallest chunk kept (in visits when
                       visits is True)
    :param trend: int, the degree of a global trend taken out first
    :param visits: bool, test the visit means rather than the exposures

    :return: dict, the chunks and the tests (NaN when fewer than two chunks
             are large enough to compare)
    """
    prob = np.zeros(data.n) if prob is None else np.asarray(prob)
    rel = np.clip(1 - prob, 1e-3, 1)
    if visits and data.nseq < data.n:
        data = _visit_means(data, rel)
        rel = np.ones(data.n)
        min_points = min(min_points, 6)
    err = data.err / np.sqrt(rel)
    tnorm = (data.time - data.tref) / max(data.baseline, 1e-9)
    phase = 2 * np.pi * (data.time - data.tref) / period
    if split == 'seasons':
        label = seasons(data.time)
    else:
        nchunk = 2 if split == 'halves' else 3
        # by visit, so a chunk never cuts a visit in two
        seq_rank = data.seq / max(data.nseq, 1)
        label = np.minimum((seq_rank * nchunk).astype(int), nchunk - 1)
    # the whole series: one sinusoid, offsets per chunk and a trend
    chunks = []
    ids = [cc for cc in np.unique(label) if np.sum(label == cc) >= min_points]
    cols = [np.cos(phase), np.sin(phase)]
    cols += [(label == cc).astype(float) for cc in ids]
    cols += [tnorm ** deg for deg in range(1, trend + 1)]
    keep = np.isin(label, ids)
    out = dict(period=period, split=split, chunks=chunks, K_all=np.nan,
               jitter=np.nan, on_visits=bool(visits), min_points=min_points,
               chi2_amplitude=np.nan, p_amplitude=np.nan,
               chi2_vector=np.nan, p_vector=np.nan)
    # a sparse series: no chunk (a season of a few visits) is large enough
    if not ids:
        return out
    try:
        if jitter is None:
            jitter = _excess_jitter(np.array(cols).T[keep], data.rv[keep],
                                    err[keep])
        err = np.sqrt(err ** 2 + jitter ** 2)
        coeffs, cov, _ = weighted_lstsq(np.array(cols).T[keep],
                                        data.rv[keep], err[keep])
    except np.linalg.LinAlgError:  # the period and the trend degenerate
        return out
    k_all = float(np.hypot(coeffs[0], coeffs[1]))
    resid_trend = np.zeros(data.n)
    for deg in range(1, trend + 1):
        resid_trend += coeffs[2 + len(ids) + deg - 1] * tnorm ** deg
    value = data.rv - resid_trend
    vecs, covs = [], []
    for cc in ids:
        sel = label == cc
        design = np.array([np.cos(phase[sel]), np.sin(phase[sel]),
                           np.ones(np.sum(sel))]).T
        try:
            cf, cv, _ = weighted_lstsq(design, value[sel], err[sel])
        except np.linalg.LinAlgError:  # a chunk that sees one phase only
            continue
        amp = float(np.hypot(cf[0], cf[1]))
        grad = np.array([cf[0], cf[1]]) / max(amp, 1e-12)
        samp = float(np.sqrt(grad @ cv[:2, :2] @ grad))
        # the phase of A cos + B sin = K cos(x - phi), phi = atan2(B, A)
        phi = float(np.arctan2(cf[1], cf[0]))
        gphi = np.array([-cf[1], cf[0]]) / max(amp ** 2, 1e-12)
        sphi = float(np.sqrt(gphi @ cv[:2, :2] @ gphi))
        chunks.append(dict(tmid=float(np.mean(data.time[sel])),
                           tmin=float(np.min(data.time[sel])),
                           tmax=float(np.max(data.time[sel])),
                           npoints=int(np.sum(sel)),
                           nvisits=int(len(np.unique(data.seq[sel]))),
                           K=amp, sK=samp, phase=phi, sphase=sphi,
                           vec=cf[:2], cov=cv[:2, :2]))
        vecs.append(cf[:2])
        covs.append(cv[:2, :2])
    out.update(K_all=k_all, jitter=float(jitter))
    if len(chunks) < 2:
        return out
    amps = np.array([ch['K'] for ch in chunks])
    samps = np.array([ch['sK'] for ch in chunks])
    wts = 1 / samps ** 2
    kmean = np.sum(wts * amps) / np.sum(wts)
    chi2_amp = float(np.sum(wts * (amps - kmean) ** 2))
    dof_amp = len(amps) - 1
    # the vectors: chi2 around their inverse-covariance weighted mean
    invs = [np.linalg.inv(cv) for cv in covs]
    vmean = np.linalg.solve(sum(invs), sum(iv @ vc for iv, vc in
                                          zip(invs, vecs)))
    chi2_vec = float(sum((vc - vmean) @ iv @ (vc - vmean)
                         for iv, vc in zip(invs, vecs)))
    dof_vec = 2 * (len(vecs) - 1)
    out.update(chi2_amplitude=chi2_amp, dof_amplitude=dof_amp,
               p_amplitude=float(stats.chi2.sf(chi2_amp, dof_amp)),
               chi2_vector=chi2_vec, dof_vector=dof_vec,
               p_vector=float(stats.chi2.sf(chi2_vec, dof_vec)),
               K_mean=float(kmean))
    return out


def indicator_check(data: RVData, period: float,
                    names: Optional[Sequence[str]] = None,
                    harmonics: Sequence[int] = (1, 2, 3),
                    pmin: float = 1.1) -> List[Dict[str, Any]]:
    """
    Whether the activity indicators have power at the period, or at a
    multiple of it

    A rotating spotted star puts power at its rotation period and at its
    harmonics P/2 and P/3 in the velocities, while the line-shape
    indicators often peak at the rotation period itself. A velocity signal
    at P with an indicator peak at 2P or 3P is therefore suspicious too.

    :param data: RVData, the series with its indicators
    :param period: float, the period of the velocity signal [days]
    :param names: list of str or None, the indicators (all when None)
    :param harmonics: list of int, the multiples of the period checked
    :param pmin: float, the shortest period of the indicator periodograms

    :return: list of dict, one per indicator: its best period, and at each
             multiple of the period its FIP and GLS power
    """
    if names is None:
        names = [nm for nm in data.indicators
                 if not nm.lower().startswith(VELOCITY_PREFIXES)]
    out = []
    for name in names:
        if name not in data.indicators:
            continue
        series = data.indicator(name)
        if series.n < 20:
            continue
        freq = frequency_grid(series.time, pmin, None, 10)
        res = fip_single(series, freq=freq, npeaks=3,
                         label=f'{name}, one signal')
        power = gls(series.time, series.rv, series.err, freq)
        entry = dict(name=name, best_period=res.best().get('period', np.nan),
                     best_fip=res.best().get('fip', np.nan), at={})
        for mult in harmonics:
            ptest = period * mult
            if ptest > 1 / freq[0]:
                continue
            idx = int(np.argmin(np.abs(freq - 1 / ptest)))
            # the best FIP within one peak width of the period itself, and
            #   within HARMONIC_TOLERANCE of its multiples
            width = 1 / series.baseline
            if mult > 1:
                width = max(width, HARMONIC_TOLERANCE / ptest)
            near = np.abs(freq - 1 / ptest) < width
            entry['at'][mult] = dict(period=ptest, fip=float(
                np.min(res.fip[near])), power=float(power[idx]))
        out.append(entry)
    return out


def eccentricity_check(data: RVData, period: float, **kwargs
                       ) -> Dict[str, Any]:
    """
    Circular against eccentric orbit, both outlier-aware

    :param data: RVData, the series
    :param period: float, the period [days]
    :param kwargs: passed to RVModel

    :return: dict, the two fits and what separates them
    """
    from koloa.fit import RVModel
    circ = RVModel(data, [dict(period=period)], **kwargs).fit(quiet=True)
    ecc = RVModel(data, [dict(period=period, eccentric=True)],
                  **kwargs).fit(quiet=True)
    eval_ = ecc.orbits()[0]['e'][0]
    return dict(circular=circ, eccentric=ecc, e=eval_,
                dlogpost=float(ecc.logpost - circ.logpost),
                at_bound=bool(eval_ > 0.9 * ecc.model.emax),
                K_circular=circ.orbits()[0]['K'],
                K_eccentric=ecc.orbits()[0]['K'])


def gp_absorption(data: RVData, period: float,
                  gp_prior: Optional[Dict[str, Any]] = None,
                  kernel: str = 'sho', **kwargs) -> Dict[str, Any]:
    """
    Does a GP explain the velocities as well without the orbit?

    Both models are fitted jointly (no sequential fit): GP alone, and GP
    plus a circular orbit. The difference of log posterior is compared
    with the penalty of three more parameters (BIC-like, 1.5 ln n).

    :param data: RVData, the series
    :param period: float, the period of the orbit [days]
    :param gp_prior: dict or None, a prior on the GP (the period from an
                     indicator, typically)
    :param kernel: str, the GP kernel
    :param kwargs: passed to RVModel

    :return: dict, the two fits and the comparison
    """
    from koloa.fit import RVModel
    gp = dict(kernel=kernel, prior=gp_prior)
    alone = RVModel(data, [], gp=gp, **kwargs).fit(nstart=1, quiet=True)
    joint = RVModel(data, [dict(period=period)], gp=gp,
                    **kwargs).fit(nstart=1, quiet=True)
    gain = float(joint.logpost - alone.logpost)
    penalty = 1.5 * np.log(data.n)
    return dict(gp_only=alone, joint=joint, gain=gain, penalty=penalty,
                absorbed=bool(gain < penalty),
                K_joint=joint.orbits()[0]['K'])


# =============================================================================
# The duck test
# =============================================================================
@dataclass
class DuckReport:
    """What the duck test found, check by check"""
    period: float
    checks: List[Dict[str, Any]] = field(default_factory=list)
    verdict: str = ''
    details: Dict[str, Any] = field(default_factory=dict)

    @property
    def flags(self) -> List[Dict[str, Any]]:
        """The checks that speak against a planet"""
        return [ch for ch in self.checks if ch['status'] == 'flag']

    def pdf(self, path: Optional[str], data: RVData,
            fipres: Optional[FIPResult] = None, target: Optional[str] = None,
            archive: bool = True, unit: str = 'both',
            outdir: Optional[str] = None) -> str:
        """
        The report as a detailed PDF: every check with its numbers, the
        figures behind them, and the known planets of the star (the NASA
        Exoplanet Archive); see koloa.report

        :param path: str, the PDF written
        :param data: RVData, the series tested
        :param fipres: FIPResult or None, the FIP of the series
        :param target: str or None, the star's name for the archive
        :param archive: bool, ask the archive (needs the network)
        :param unit: str, the outlier unit of the orbit fitted for the folds
        :param outdir: str or None, a folder for everything: the report (there
                       when path is None), each figure as a PDF, the text and
                       a JSON summary

        :return: str, the path of the report
        """
        from koloa.report import duck_pdf
        return duck_pdf(self, data, path, fipres=fipres, target=target,
                        archive=archive, unit=unit, outdir=outdir)

    def text(self) -> str:
        """The report, as text"""
        lines = [f'Duck test at P = {self.period:.4f} d', '']
        mark = {'pass': 'quacks like a planet', 'flag': 'DOES NOT QUACK',
                'info': 'note'}
        for ch in self.checks:
            lines.append(f'[{mark[ch["status"]]}] {ch["name"]}: '
                         f'{ch["summary"]}')
        sols = self.details.get('alias_solutions')
        if sols:
            from koloa.aliases import describe, plan_text
            lines += ['', 'The period and its aliases (an orbit fitted at '
                          'each):'] + ['  ' + line for line in describe(sols)]
            if self.details.get('alias_plan'):
                lines += [''] + plan_text(self.details['alias_plan'])
        lines += ['', f'VERDICT: {self.verdict}']
        return '\n'.join(lines)


def duck_test(data: RVData, period: float,
              fipres: Optional[FIPResult] = None,
              gauss_fip: Optional[FIPResult] = None,
              prob: Optional[np.ndarray] = None,
              activity: Optional[str] = 'auto', gp: bool = True,
              gp_kernel: str = 'sho', unit: str = 'both',
              fip_threshold: float = 0.01, quiet: bool = False,
              pdf: Optional[str] = None, target: Optional[str] = None,
              archive: bool = True, outdir: Optional[str] = None,
              tess: Union[None, bool, Dict[str, Any]] = None,
              aliases: Optional[bool] = None, site: Optional[str] = None,
              nightly: Optional[bool] = None) -> DuckReport:
    """
    Every test of planethood koloa knows, at one period

    :param data: RVData, the series
    :param period: float, the period [days]
    :param fipres: FIPResult or None, the outlier-aware FIP
    :param gauss_fip: FIPResult or None, the gaussian FIP, for comparison
    :param prob: np.ndarray or None, the outlier probability of each point
                 (from fipres when None)
    :param activity: str or None, the indicator that sets the GP period
                     prior ('auto': the first of ACTIVITY_INDICATORS found)
    :param gp: bool, run the GP absorption test (slower)
    :param gp_kernel: str, the GP kernel
    :param unit: str, the outlier unit (sequence or point)
    :param fip_threshold: float, the FIP below which a signal is detected
    :param quiet: bool, no log lines
    :param pdf: str or None, also write a detailed PDF report there (every
                check, its figures, and the known planets of the star from
                the NASA Exoplanet Archive; DuckReport.pdf)
    :param target: str or None, the star's name for the archive (the name
                   of the series when None)
    :param archive: bool, ask the archive for the PDF (needs the network)
    :param outdir: str or None, a folder for everything: the PDF report
                   (there when pdf is None), each of its figures as a PDF,
                   the text of the test and a JSON summary
    :param tess: None, bool or dict, whether the TESS light curves of the
                 star (by target, or the name of the series) are asked for
                 a photometric peak at the period, its half, third or double
                 (koloa.tess; needs the network): None asks when a report is
                 written (pdf or outdir) and archive is True; a dict is light
                 curves already fetched (koloa.tess.light_curves)
    :param aliases: bool or None, fit an orbit at the period and at each
                    alias that holds some of the probability, with its share
                    (koloa.aliases), and when the period is ambiguous, plan
                    the observations that lift the alias (needs the star's
                    position: target, archive True); None fits them when a
                    report is written
    :param site: str or None, the observatory of the plan (a key of
                 koloa.aliases.SITES; from the instruments when None)
    :param nightly: bool, test the nightly means (RVData.nightly, koloa's
                    default; the outlier probability of each exposure is
                    averaged into its night, and an outlier is a night);
                    False tests the exposures; None is koloa.data.NIGHTLY

    Planet or no planet is decided on the FIP of the period OR any of its
    aliases (1 day, 1 year, 1 month): whether there is a planet, whichever
    alias it is. Which alias it is, is reported apart (the 'period' check,
    a note that never counts against a planet).

    :return: DuckReport
    """
    from koloa.fit import period_prior_from_indicator
    if prob is None and fipres is not None:
        prob = fipres.outlier_prob
    if nightly is None:
        from koloa import data as kdata
        nightly = kdata.NIGHTLY
    if nightly:
        means = data.nightly()
        if means is not data:
            from koloa.data import night_index
            night = night_index(data)
            if prob is not None and len(prob) == data.n:
                prob = (np.bincount(night, np.asarray(prob, dtype=float))
                        / np.bincount(night))
            data = means
            if unit == 'both':
                unit = 'point'
    report = DuckReport(period=period)
    prob = np.zeros(data.n) if prob is None else np.asarray(prob)
    # the outliers of the model, out of the checks that have no outlier
    #   model of their own (the jackknife's GLS, the coherence's least
    #   squares): a failed exposure thousands of sigma off would be them
    good = prob < 0.5
    clean = data.select(good) if (np.any(~good) and np.sum(good) > 10) \
        else data
    cprob = prob[good] if clean is not data else prob
    freq0 = 1 / period
    kw = dict(likelihood='mixture', unit=unit)
    # -------------------------------------------------------------------------
    # 1. significance
    # -------------------------------------------------------------------------
    if fipres is not None:
        # the best interval that holds the period (its grid point can be a
        #   fraction of an interval off the peak)
        width = 1 / data.baseline
        fipv = fipres.fip_containing(period, width)
        # planet or no planet: the period OR any of its aliases
        fam = fipres.family_containing(period, width)
        decisive = fam if fam is not None else fipv
        wpow = float(window(data.time, freq0)[0])
        summary = (f'planet or no planet (the period or any of its aliases):'
                   f' FIP = {fam:.2e}; the period alone: FIP = {fipv:.2e}'
                   if fam is not None else f'FIP = {fipv:.2e}')
        summary += ' (outlier-aware)'
        if gauss_fip is not None:
            summary += (f', {gauss_fip.fip_containing(period, 1 / data.baseline):.2e}'
                        f' (gaussian)')
        summary += f', window power {wpow:.3f}'
        # the aliases, and whether one of them is as good
        idx = int(np.argmin(np.abs(fipres.freq - freq0)))
        best_alias = None
        for peak in fipres.peaks:
            if abs(peak['freq'] - fipres.freq[idx]) < 2 / data.baseline:
                for alias in peak.get('aliases', []):
                    if best_alias is None or alias['fip'] < best_alias['fip']:
                        best_alias = alias
        # (with a family, the 'period' check below says it all)
        if best_alias is not None and fam is None:
            summary += (f'; best alias {best_alias["period"]:.3f} d '
                        f'({best_alias["name"]}) FIP {best_alias["fip"]:.2e}')
        status = 'pass' if decisive < fip_threshold else 'flag'
        report.checks.append(dict(name='significance', status=status,
                                  summary=summary, fip=fipv,
                                  family_fip=fam, window=wpow))
        report.details['fip'] = fipv
        report.details['family_fip'] = fam
        # which alias: a note, never a flag (the planet is there all the
        #   same, and observations can tell)
        if fam is not None:
            odds = fipres.alias_odds(period, width)
            report.details['alias_odds'] = odds
            shown = [mem for mem in odds if mem['share'] >= 1e-3]
            text = '; '.join(f'{mem["period"]:.4f} d '
                             f'({"the period" if mem["name"] == "P" else mem["name"] + " alias"})'
                             f' {100 * mem["share"]:.1f} %' for mem in shown)
            if len(shown) > 1:
                text = 'ambiguous between its aliases: ' + text
            else:
                text = f'unambiguous: {text} of the probability'
            report.checks.append(dict(name='period', status='info',
                                      summary=text))
        # the orbit at each alias, and when to observe to tell them apart
        if aliases is None:
            aliases = bool(pdf or outdir)
        if aliases:
            from koloa import aliases as kal
            try:
                sols = kal.solutions(data, period, fipres, unit=unit)
                report.details['alias_solutions'] = sols
                if kal.ambiguous(sols) and archive:
                    from koloa.archive import resolve
                    where = site or kal.site_of(data.instruments)
                    ident = resolve(target or data.name)
                    if where is None:
                        log('duck test: no observatory for the plan (give '
                            'site=)', 'warn')
                    elif ident.get('ra') is not None:
                        report.details['alias_plan'] = kal.plan(
                            sols, ident['ra'], ident['dec'], where)
            except Exception as err:  # a help, not a stop
                log(f'duck test: aliases at {period:.4f} d: {err}', 'warn')
    # -------------------------------------------------------------------------
    # 2. robustness: who holds the peak up
    # -------------------------------------------------------------------------
    freq = frequency_grid(clean.time, 1.1, None, 10)
    jack = jackknife(clean, freq, unit='point' if unit == 'point'
                     else 'sequence')
    infl = influence(jack, freq, period)
    idx = int(np.argmin(np.abs(freq - freq0)))
    peakpow = float(jack['power'][idx])
    top = float(np.max(infl)) if len(infl) else 0.0
    frac_top = top / peakpow if peakpow > 0 else 0.0
    status = 'flag' if frac_top > 0.3 else 'pass'
    report.checks.append(dict(
        name='robustness', status=status,
        summary=(f'the most influential '
                 f'{"point" if unit == "point" else "visit"} holds up '
                 f'{100 * frac_top:.0f}% of the GLS power at the period'
                 + (f' ({data.n - clean.n} outlying exposures left out)'
                    if clean is not data else '')),
        fraction=frac_top))
    report.details['jackknife'] = jack
    report.details['jackknife_freq'] = freq
    report.details['jackknife_data'] = clean
    # -------------------------------------------------------------------------
    # 3. coherence in time
    # -------------------------------------------------------------------------
    coh = {}
    for split in ('halves', 'seasons'):
        coh[split] = coherence(clean, period, prob=cprob, split=split)
    main = coh['halves']
    pmin = min(val for val in (coh['halves']['p_vector'],
                               coh['seasons']['p_vector'])
               if np.isfinite(val)) if any(np.isfinite(
                   [coh['halves']['p_vector'], coh['seasons']['p_vector']])) \
        else np.nan
    amps = ', '.join(f'{ch["K"]:.1f}+-{ch["sK"]:.1f}' for ch in
                     main['chunks'])
    parts = []
    if len(main['chunks']) >= 2:
        parts.append(f'K by half: {amps} m/s (p = '
                     f'{main["p_amplitude"]:.1e})')
    if len(coh['seasons']['chunks']) >= 2:
        parts.append(f'amplitude+phase by season p = '
                     f'{coh["seasons"]["p_vector"]:.1e}')
    for split, what in (('halves', 'half'), ('seasons', 'season')):
        if len(coh[split]['chunks']) < 2:
            parts.append(f'too few visits per {what} to compare (at least '
                         f'{coh[split]["min_points"]} needed)')
    if np.isfinite(pmin):
        status = 'flag' if pmin < 0.01 else 'pass'
    else:
        status = 'info'
    report.checks.append(dict(name='coherence', status=status,
                              summary='; '.join(parts), p=pmin))
    report.details['coherence'] = coh
    # -------------------------------------------------------------------------
    # 4. shape: circular or eccentric
    # -------------------------------------------------------------------------
    ecc = eccentricity_check(data, period, **kw)
    if ecc['at_bound']:
        status, note = 'flag', 'the eccentric fit runs to the prior edge'
    elif ecc['dlogpost'] > 5 and ecc['e'] > 0.1:
        status, note = 'pass', 'a well-defined eccentric orbit'
    else:
        status, note = 'pass', 'consistent with a circular orbit'
    report.checks.append(dict(
        name='shape', status=status,
        summary=(f'e = {ecc["e"]:.2f}, Delta ln post = {ecc["dlogpost"]:.1f} '
                 f'over circular; {note}; K = {ecc["K_circular"][0]:.2f} '
                 f'+- {ecc["K_circular"][1]:.2f} m/s (circular)')))
    report.details['eccentricity'] = ecc
    # -------------------------------------------------------------------------
    # 5. activity indicators
    # -------------------------------------------------------------------------
    inds = indicator_check(data, period)
    report.details['indicators'] = inds
    hits = []
    for entry in inds:
        for mult, val in entry['at'].items():
            if val['fip'] < 0.05:
                hits.append(f'{entry["name"]} at {mult}P '
                            f'({val["period"]:.2f} d, FIP {val["fip"]:.1e})')
    near_best = [f'{entry["name"]} peaks at {entry["best_period"]:.2f} d'
                 for entry in inds if entry['best_fip'] < 0.05]
    if hits:
        status = 'flag'
        summary = 'power at the period or a multiple: ' + '; '.join(hits)
    else:
        status = 'pass'
        summary = 'no indicator has significant power at P, 2P or 3P'
    if near_best:
        summary += ' (' + '; '.join(near_best) + ')'
    report.checks.append(dict(name='indicators', status=status,
                              summary=summary))
    # -------------------------------------------------------------------------
    # 5b. the photometry of TESS: spots put the rotation and its harmonics
    #     in the brightness; a planet does not
    # -------------------------------------------------------------------------
    if tess is None:
        tess = bool(pdf or outdir) and archive
    if tess is not False:
        from koloa import tess as ktess
        try:
            lcs = tess if isinstance(tess, dict) else \
                ktess.light_curves(target or data.name)
            phot = ktess.periodicity(lcs, [period])
            chk = phot['checks'][0] if phot['checks'] else dict(
                status='info', summary='no TESS light curve')
            report.details['tess'] = phot
            report.checks.append(dict(name='photometry',
                                      status=chk['status'],
                                      summary='TESS: ' + chk['summary']))
        except Exception as err:  # no network, or an unknown star
            report.checks.append(dict(
                name='photometry', status='info',
                summary=f'TESS could not be asked ({err})'))
    # -------------------------------------------------------------------------
    # 6. absorption by a GP whose period an indicator sets
    # -------------------------------------------------------------------------
    if gp:
        name = None
        if activity == 'auto':
            name = next((nm for nm in ACTIVITY_INDICATORS
                         if nm in data.indicators), None)
        elif activity is not None:
            name = activity
        gp_prior = None
        if name is not None:
            if not quiet:
                log(f'Duck test: activity period from {name}')
            gp_prior, act_fit = period_prior_from_indicator(
                data, name, kernel=gp_kernel, quiet=True)
            report.details['activity_fit'] = act_fit
            report.details['activity_prior'] = gp_prior
            report.details['activity_indicator'] = name
        if not quiet:
            log('Duck test: joint GP + orbit against GP alone')
        # a GP couples every point, so its outliers are whole sequences
        gpkw = dict(kw, unit='sequence' if unit == 'both' else unit)
        absb = gp_absorption(data, period, gp_prior=gp_prior,
                             kernel=gp_kernel, **gpkw)
        report.details['absorption'] = absb
        prot = ''
        if gp_prior is not None:
            mu, sd = gp_prior['log_period']
            prot = (f' (GP period from {name}: {np.exp(mu):.1f} d, '
                    f'+-{100 * sd:.0f}%)')
        kj = absb['K_joint']
        status = 'flag' if absb['absorbed'] else 'pass'
        report.checks.append(dict(
            name='GP absorption', status=status,
            summary=(f'the orbit adds Delta ln post = {absb["gain"]:.1f} '
                     f'(needs > {absb["penalty"]:.1f}) to a GP{prot}; '
                     f'K = {kj[0]:.2f} +- {kj[1]:.2f} m/s in the joint fit')))
    # -------------------------------------------------------------------------
    # the verdict
    # -------------------------------------------------------------------------
    nflag = len(report.flags)
    decisive = report.details.get('family_fip')
    if decisive is None:
        decisive = report.details.get('fip', 0.0)
    significant = decisive < fip_threshold
    if fipres is not None and not significant:
        others = [ch['name'] for ch in report.flags
                  if ch['name'] != 'significance']
        verdict = (f'NOT DETECTED: FIP = {decisive:.2g} for the period or '
                   f'any of its aliases, once the outliers and the noise are '
                   f'accounted for')
        if others:
            verdict += ('; in addition, ' + ', '.join(others)
                        + (' speaks' if len(others) == 1 else ' speak')
                        + ' against a planet')
    elif nflag == 0:
        verdict = 'PLANET CANDIDATE: it looks, swims and quacks like a planet'
    elif nflag == 1:
        verdict = ('INCONCLUSIVE: one test speaks against a planet (' +
                   report.flags[0]['name'] + ')')
    else:
        verdict = ('NOT A PLANET (activity or systematics more likely): '
                   + ', '.join(ch['name'] for ch in report.flags))
    # a detection whose period is ambiguous says so
    odds = report.details.get('alias_odds') or []
    if significant and len(odds) > 1 and odds[1]['share'] >= 0.01:
        verdict += (' (the period is ambiguous: ' + ', '.join(
            f'{mem["period"]:.4f} d {100 * mem["share"]:.1f} %'
            for mem in odds if mem['share'] >= 0.01) + ')')
    report.verdict = verdict
    if not quiet:
        for line in report.text().split('\n'):
            if line:
                log(line, 'value' if line.startswith('VERDICT') else 'info')
    if pdf or outdir:
        report.pdf(pdf, data, fipres=fipres, target=target,
                   archive=archive, unit=unit, outdir=outdir)
    return report


# =============================================================================
# End of code
# =============================================================================

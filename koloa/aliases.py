#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The aliases of a signal: which period it is, and how to find out.

A periodic signal sampled from the ground also shows at its aliases, |f -
fs| and f + fs for every sampling frequency fs (a sidereal day, a year, a
synodic month). The FIP answers two questions, which koloa keeps apart:

- is there a planet? The FIP of the period OR any of its aliases
  (FIPResult.family_containing): the probability that none of them holds a
  signal. Planet or no planet is decided on it: a planet at 6.92 d whose
  1-day alias at 1.165 d holds 2 % of the probability is detected all the
  same;
- which period is it? The share of the probability that each alias holds
  (FIPResult.alias_odds).

For a tentative planet, this module adds:

- solutions(): an orbit fitted at the period and at each alias that holds
  some of the probability (circular, outliers modelled), with its share of
  the probability and its log posterior against the best;
- figure(): the velocities folded at each of them, side by side;
- plan(): the nights, and the hours, when the solutions differ most, to lift
  the alias: one visit where they predict the most different velocities,
  and two visits in a night where the change they predict differs most (a
  pair does not depend on the offset of the instrument); plan_text() and
  plan_figure() show it.

The sky is computed with low-precision formulas (the Sun to about a
hundredth of a degree, the star without precession, to a few tenths of a
degree), enough to plan a night. Times are UT.

Created on 2026-09-30

@author: artigau
"""
import re
import time as _time
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from koloa.data import RVData
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: a period is ambiguous when an alias holds at least this share of the
#:  probability: its alternatives are shown, and a plan to lift it is made
AMBIGUOUS = 0.01
#: the aliases fitted and shown: those that hold at least this share
MIN_SHARE = 1e-3
#: observatories: latitude, longitude (east positive) [degrees]
SITES = {'CFHT': (19.8253, -155.4689), 'La Silla': (-29.2567, -70.7300),
         'Paranal': (-24.6272, -70.4042), 'Calar Alto': (37.2236, -2.5463),
         'La Palma': (28.7540, -17.8890), 'Subaru': (19.8255, -155.4761),
         'Keck': (19.8260, -155.4747), 'OHP': (43.9308, 5.7133),
         'Gemini North': (19.8238, -155.4690)}
#: the site of an instrument, by the start of its name (letters only, upper
#:  case; HARPS-N before HARPS)
INSTRUMENT_SITES = [('HARPSN', 'La Palma'), ('HARPN', 'La Palma'),
                    ('HARPS', 'La Silla'), ('NIRPS', 'La Silla'),
                    ('CORALIE', 'La Silla'), ('ESPRESSO', 'Paranal'),
                    ('SPIROU', 'CFHT'), ('CARMENES', 'Calar Alto'),
                    ('CARM', 'Calar Alto'), ('IRD', 'Subaru'),
                    ('SOPHIE', 'OHP'), ('HIRES', 'Keck'), ('KPF', 'Keck'),
                    ('MAROON', 'Gemini North')]


# =============================================================================
# The family of a period
# =============================================================================
def label(name: str) -> str:
    """'P' is the period tested; the others are aliases"""
    return 'the period' if name == 'P' else f'{name} alias'


def members(period: float, fmin: float = 0.0,
            fmax: float = np.inf) -> List[float]:
    """
    The frequencies of a period and of its aliases

    :param period: float [days]
    :param fmin: float, the lowest frequency kept [1/day]
    :param fmax: float, the highest [1/day]

    :return: list of float [1/day], the period's own first
    """
    from koloa.periodogram import aliases
    return [1.0 / period] + [al['freq'] for al in
                             aliases(1.0 / period, fmin, fmax)]


def same_family(period1: float, period2: float, width: float) -> bool:
    """
    Whether a period is one of the aliases of another (or the same)

    :param period1: float [days]
    :param period2: float [days]
    :param width: float, the tolerance, 1/T [1/day]

    :return: bool
    """
    return any(abs(freq - 1.0 / period2) <= width
               for freq in members(period1))


# =============================================================================
# The orbit at each alias
# =============================================================================
def _fit(data: RVData, period: float, unit: str, width: float):
    """a circular orbit within one interval of the period, outliers
    modelled"""
    from koloa.fit import RVModel
    freq = 1.0 / period
    prange = (1.0 / (freq + 0.5 * width), 1.0 / max(freq - 0.5 * width,
                                                    1e-9))
    try:
        return RVModel(data, [dict(period=period, period_range=prange)],
                       likelihood='mixture', unit=unit).fit(nstart=1,
                                                            quiet=True)
    except Exception as err:  # an alias that cannot be fitted is skipped
        log(f'aliases: no fit at {period:.4f} d ({err})', 'warn')
        return None


def solutions(data: RVData, period: float, fipres: Optional[Any] = None,
              unit: str = 'both', min_share: float = MIN_SHARE,
              width: Optional[float] = None) -> List[Dict[str, Any]]:
    """
    An orbit at the period and at each of its aliases, with how probable
    each is

    :param data: RVData, the series
    :param period: float, the period of the tentative planet [days]
    :param fipres: FIPResult or None, the FIP of the series: the share of
                   the probability of each alias (without it, every alias
                   on the grid is fitted and only the fits compare them)
    :param unit: str, the outlier unit of the fits
    :param min_share: float, the aliases fitted hold at least this share
    :param width: float or None, the width of an interval (1/T when None)

    :return: list of dict, the most probable first: name ('P' or the
             sampling of the alias), period (fitted), share (of the
             probability, from the FIP), fip (of the interval alone), K and
             tc (value, minus, plus), logpost (of the fit), dlogpost (from
             the best fit), likelihood_share (the fits' relative likelihood,
             exp(dlogpost) normalised), fit (FitResult)
    """
    width = width or 1.0 / data.baseline
    if fipres is not None and getattr(fipres, 'family', None) is not None:
        odds = fipres.alias_odds(period, width)
    else:
        from koloa.periodogram import aliases
        odds = [dict(name='P', period=period, fip=np.nan, share=np.nan)]
        odds += [dict(name=al['name'], period=al['period'], fip=np.nan,
                      share=np.nan)
                 for al in aliases(1.0 / period, 0.5 / data.baseline, 2.0)
                 if abs(al['freq'] - 1.0 / period) > 1.5 * width]
    out = []
    for mem in odds:
        if (mem['name'] != 'P' and np.isfinite(mem['share'])
                and mem['share'] < min_share):
            continue
        fit = _fit(data, mem['period'], unit, width)
        if fit is None:
            continue
        orb = fit.orbits()[0]
        out.append(dict(name=mem['name'], period=float(orb['P'][0]),
                        grid_period=float(mem['period']),
                        share=float(mem['share']), fip=float(mem['fip']),
                        K=[float(val) for val in orb['K']],
                        tc=[float(val) for val in orb['tc']],
                        logpost=float(fit.logpost), fit=fit))
    if not out:
        return out
    best = max(sol['logpost'] for sol in out)
    wts = np.array([np.exp(sol['logpost'] - best) for sol in out])
    for sol, wgt in zip(out, wts / wts.sum()):
        sol['dlogpost'] = sol['logpost'] - best
        sol['likelihood_share'] = float(wgt)
    out.sort(key=lambda sol: (-sol['share'] if np.isfinite(sol['share'])
                              else -sol['likelihood_share']))
    return out


def ambiguous(sols: Sequence[Dict[str, Any]],
              level: float = AMBIGUOUS) -> bool:
    """Whether an alias other than the most probable holds at least level
    of the probability (the fits' relative likelihood without a FIP)"""
    if len(sols) < 2:
        return False
    share = sols[1]['share']
    if not np.isfinite(share):
        share = sols[1]['likelihood_share']
    return bool(share >= level)


def light(sols: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """the solutions without their fits (for a JSON file)"""
    return [{key: val for key, val in sol.items() if key != 'fit'}
            for sol in sols]


def describe(sols: Sequence[Dict[str, Any]]) -> List[str]:
    """
    One line per solution: its period, its share of the probability, its
    amplitude and its log posterior against the best

    :return: list of str
    """
    lines = []
    for sol in sols:
        share = (f'{100 * sol["share"]:.2f} % of the probability'
                 if np.isfinite(sol['share']) else
                 f'{100 * sol["likelihood_share"]:.2f} % of the likelihood')
        fipv = (f', FIP of it alone {sol["fip"]:.1e}'
                if np.isfinite(sol['fip']) else '')
        lines.append(f'{sol["period"]:.4f} d ({label(sol["name"])}): '
                     f'{share}{fipv}; K = {sol["K"][0]:.2f} '
                     f'+{sol["K"][2]:.2f}/-{sol["K"][1]:.2f} m/s; '
                     f'ln posterior {sol["dlogpost"]:+.1f} from the best fit')
    return lines


def figure(sols: Sequence[Dict[str, Any]], title: Optional[str] = None,
           nmax: int = 4):
    """
    The velocities folded at the period and at its most probable aliases,
    side by side, each with its share of the probability

    :param sols: list of dict, from solutions()
    :param title: str or None
    :param nmax: int, the most solutions shown

    :return: matplotlib figure
    """
    from koloa import plotting as kplot
    shown = list(sols[:nmax])
    fig, axes = kplot.plt.subplots(1, len(shown) + 1,
                                   figsize=(3.3 * (len(shown) + 1) + 0.6, 3.3),
                                   squeeze=False)
    for ax, sol in zip(axes[0], shown):
        share = (f'{100 * sol["share"]:.2f} %' if np.isfinite(sol['share'])
                 else f'{100 * sol["likelihood_share"]:.2f} % (fits)')
        kplot.phase(sol['fit'], ax=ax, level='point', legend=False,
                    colorbar=False,
                    title=f'{sol["period"]:.4f} d, {label(sol["name"])}\n'
                          f'{share}; $\\Delta\\ln$ post '
                          f'{sol["dlogpost"]:+.1f}')
    for ax in axes[0][1:len(shown)]:
        ax.set_ylabel('')
    _residuals(axes[0][-1], shown)
    if title:
        fig.suptitle(title, x=0.01, ha='left', fontsize=9)
    fig.tight_layout()
    return fig


def residuals(sol: Dict[str, Any]) -> np.ndarray:
    """the residuals of a solution's fit [m/s], its outliers left out"""
    fit = sol['fit']
    data = fit.model.data
    resid = data.rv - fit.model.mean_model(fit.theta)
    rel = 1 - (fit.outlier_prob if fit.outlier_prob is not None
               else np.zeros(data.n))
    return resid[rel > 0.5]


def _residuals(ax, sols: Sequence[Dict[str, Any]]):
    """the histograms of the residuals of every solution, overplotted"""
    from koloa import plotting as kplot
    colours = [kplot.C['koloa'], kplot.C['outlier'], kplot.C['text'],
               kplot.C['muted']]
    allres = [residuals(sol) for sol in sols]
    span = max(float(np.max(np.abs(res))) for res in allres if len(res))
    bins = np.linspace(-span, span, 25)
    for sol, res, colour in zip(sols, allres, colours):
        ax.hist(res, bins=bins, histtype='step', lw=1.4, color=colour,
                label=f'{sol["period"]:.4f} d: rms {np.std(res):.2f}')
    ax.set_xlabel('residual [m s$^{-1}$]')
    ax.set_ylabel('points')
    ax.set_title('residuals of each fit', loc='left')
    ax.legend(loc='upper left', fontsize=6.5)


# =============================================================================
# Lifting the alias
# =============================================================================
def site_of(instruments: Sequence[str]) -> Optional[str]:
    """
    The observatory of the first instrument koloa knows

    :param instruments: list of str

    :return: str (a key of SITES) or None
    """
    for inst in instruments:
        key = re.sub('[^A-Z]', '', str(inst).upper())
        for name, site in INSTRUMENT_SITES:
            if key.startswith(name):
                return site
    return None


def rjd_now() -> float:
    """now, as BJD - 2400000 (to a few minutes)"""
    return _time.time() / 86400.0 + 2440587.5 - 2400000.0


def utc(rjd: float) -> datetime:
    """a time [BJD - 2400000] as a UT date (MJD 0 is 1858-11-17)"""
    return datetime(1858, 11, 17) + timedelta(days=float(rjd) - 0.5)


def _sun(jd: np.ndarray):
    """the right ascension and declination of the Sun [deg]"""
    day = jd - 2451545.0
    anom = np.radians(357.529 + 0.98560028 * day)
    mlon = 280.459 + 0.98564736 * day
    lon = np.radians(mlon + 1.915 * np.sin(anom) + 0.020 * np.sin(2 * anom))
    obl = np.radians(23.439 - 0.00000036 * day)
    ra = np.degrees(np.arctan2(np.cos(obl) * np.sin(lon), np.cos(lon)))
    dec = np.degrees(np.arcsin(np.sin(obl) * np.sin(lon)))
    return ra, dec


def altitude(ra, dec, jd: np.ndarray, lat: float, lon: float) -> np.ndarray:
    """
    The altitude of a position of the sky from a site

    :param ra: float or np.ndarray [deg]
    :param dec: float or np.ndarray [deg]
    :param jd: np.ndarray, the Julian date
    :param lat: float, the site's latitude [deg]
    :param lon: float, its longitude, east positive [deg]

    :return: np.ndarray [deg]
    """
    gmst = (280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360
    hour = np.radians(gmst + lon - ra)
    phi, delta = np.radians(lat), np.radians(dec)
    return np.degrees(np.arcsin(np.sin(phi) * np.sin(delta) + np.cos(phi)
                                * np.cos(delta) * np.cos(hour)))


def _predict(sol: Dict[str, Any], times: np.ndarray, inst: Optional[str],
             ndraw: int = 100):
    """the velocity a solution predicts, and its spread over the draws of
    its fit (period, phase, amplitude and offset)"""
    from koloa import kepler
    fit = sol['fit']
    model = fit.model
    sel = (model.data.inst == inst) if inst in model.data.instruments \
        else np.ones(model.data.n, bool)

    def one(theta):
        period, tperi, ecc, omega, amp = model.orbit(theta, 0)
        offset = float(np.median(model.systematics(theta)[sel]))
        return kepler.rv_keplerian(times, period, tperi, ecc, omega,
                                   amp) + offset
    draws = np.array([one(theta) for theta in fit.samples(ndraw)])
    return one(fit.theta), np.std(draws, axis=0)


def visit_scatter(sol: Dict[str, Any]) -> float:
    """the scatter of the visit means about a fit [m/s] (robust): the
    noise of one more visit, jitter and activity included"""
    from koloa.data import robust_std
    fit = sol['fit']
    data = fit.model.data
    resid = data.rv - fit.model.mean_model(fit.theta)
    rel = 1 - (fit.outlier_prob if fit.outlier_prob is not None
               else np.zeros(data.n))
    wts = rel / data.err ** 2
    wsum = np.bincount(data.seq, wts)
    good = wsum > 0
    means = np.bincount(data.seq, wts * resid)[good] / wsum[good]
    return float(robust_std(means))


def plan(sols: Sequence[Dict[str, Any]], ra: float, dec: float,
         site: str, start: Optional[float] = None, ndays: float = 60.0,
         min_alt: float = 30.0, sun_alt: float = -12.0,
         step: float = 10.0 / 1440, sigma: Optional[float] = None,
         inst: Optional[str] = None, top: int = 5) -> Dict[str, Any]:
    """
    When to observe to lift the alias: for every night the star can be
    observed, the visit where the most probable solution and each
    alternative predict the most different velocities, and the two visits
    whose predicted changes differ most

    A visit that sees a difference dv adds dv^2 / (sigma^2 + the spread of
    the two predictions) to the chi2 of the wrong solution: every unit of
    it multiplies the odds by exp(1/2). A pair in one night measures a
    change, which does not depend on the offset of the instrument.

    :param sols: list of dict, from solutions() (the first is the most
                 probable)
    :param ra: float, the star's right ascension (J2000) [deg]
    :param dec: float, its declination [deg]
    :param site: str, a key of SITES
    :param start: float or None, the first day [BJD - 2400000] (now)
    :param ndays: float, the days planned
    :param min_alt: float, the lowest altitude [deg] (30: airmass 2)
    :param sun_alt: float, the Sun below this altitude [deg] (-12: nautical
                    twilight)
    :param step: float, the time step [days]
    :param sigma: float or None, the noise of one visit [m/s] (the scatter
                  of the visit means about the best fit when None)
    :param inst: str or None, the instrument whose offset the predictions
                 take (the first of the series)
    :param top: int, the best nights listed

    :return: dict, site, start, ndays, sigma, solutions, nights (date, UT
             window, highest altitude, and per alternative the best single
             visit and the best pair), best_single and best_pair (the top
             nights against the main alternative)
    """
    lat, lon = SITES[site]
    start = rjd_now() if start is None else float(start)
    best, alts = sols[0], list(sols[1:])
    inst = inst or best['fit'].model.data.instruments[0]
    sigma = visit_scatter(best) if sigma is None else float(sigma)
    out = dict(site=site, lat=lat, lon=lon, ra=ra, dec=dec, start=start,
               ndays=ndays, min_alt=min_alt, sun_alt=sun_alt, sigma=sigma,
               instrument=inst,
               solutions=[dict(name=sol['name'], period=sol['period'],
                               share=sol['share']) for sol in sols],
               nights=[], best_single=[], best_pair=[])
    if not alts:
        return out
    times = np.arange(start, start + ndays, step)
    jd = times + 2400000.0
    sra, sdec = _sun(jd)
    salt = altitude(sra, sdec, jd, lat, lon)
    talt = altitude(ra, dec, jd, lat, lon)
    idx = np.where((salt < sun_alt) & (talt > min_alt))[0]
    if not len(idx):
        log(f'plan: {site} never sees the star above {min_alt} deg at night '
            f'in the {ndays:.0f} days from {utc(start):%Y-%m-%d}', 'warn')
        return out
    tobs = times[idx]
    vbest, sbest = _predict(best, tobs, inst)
    preds = [_predict(alt, tobs, inst) for alt in alts]
    groups = np.split(np.arange(len(idx)),
                      np.where(np.diff(tobs) > 0.125)[0] + 1)
    for grp in groups:
        tt = tobs[grp]
        night = dict(date=f'{utc(tt[0]):%Y-%m-%d}', start=float(tt[0]),
                     end=float(tt[-1]),
                     window=f'{utc(tt[0]):%H:%M}-{utc(tt[-1]):%H:%M} UT',
                     alt_max=float(np.max(talt[idx][grp])), versus=[])
        for alt, (valt, salt_) in zip(alts, preds):
            dv = vbest[grp] - valt[grp]
            var = sigma ** 2 + sbest[grp] ** 2 + salt_[grp] ** 2
            chi1 = dv ** 2 / var
            i1 = int(np.argmax(chi1))
            # the change over the night: offset-free
            ddv = ((vbest[grp][:, None] - vbest[grp][None, :])
                   - (valt[grp][:, None] - valt[grp][None, :]))
            chi2 = np.triu(ddv ** 2 / (2 * sigma ** 2), 1)
            i2, j2 = np.unravel_index(int(np.argmax(chi2)), chi2.shape)
            night['versus'].append(dict(
                name=alt['name'], period=alt['period'],
                single=dict(time=float(tt[i1]), ut=f'{utc(tt[i1]):%H:%M}',
                            dv=float(dv[i1]), chi2=float(chi1[i1])),
                pair=dict(time1=float(tt[i2]), time2=float(tt[j2]),
                          ut1=f'{utc(tt[i2]):%H:%M}',
                          ut2=f'{utc(tt[j2]):%H:%M}',
                          ddv=float(ddv[i2, j2]),
                          chi2=float(chi2[i2, j2]))))
        out['nights'].append(night)
    main = [night['versus'][0] for night in out['nights']]
    order1 = np.argsort([-vs['single']['chi2'] for vs in main])[:top]
    order2 = np.argsort([-vs['pair']['chi2'] for vs in main])[:top]
    out['best_single'] = [int(ii) for ii in order1]
    out['best_pair'] = [int(ii) for ii in order2]
    return out


def plan_text(pln: Dict[str, Any]) -> List[str]:
    """
    The plan in words: the best nights against the main alternative, one
    visit and a pair

    :return: list of str
    """
    if not pln.get('nights'):
        return [f'{pln.get("site")}: no night with the star above '
                f'{pln.get("min_alt")} deg in the {pln.get("ndays", 0):.0f} '
                f'days planned']
    sols = pln['solutions']
    lines = [f'Lifting the alias between {sols[0]["period"]:.4f} d and '
             f'{sols[1]["period"]:.4f} d from {pln["site"]} '
             f'({pln["instrument"]}), from {utc(pln["start"]):%Y-%m-%d} for '
             f'{pln["ndays"]:.0f} days (star above {pln["min_alt"]:.0f} deg, '
             f'Sun below {pln["sun_alt"]:.0f} deg; one visit '
             f'{pln["sigma"]:.2f} m/s). Each unit of Delta chi2 multiplies '
             f'the odds by exp(1/2), about 1.65 (10 units: 150 times).',
             'One visit, where the two predict the most different '
             'velocities:']
    for ii in pln['best_single']:
        night = pln['nights'][ii]
        vs = night['versus'][0]['single']
        lines.append(f'  {night["date"]} {vs["ut"]} UT: {vs["dv"]:+.1f} m/s '
                     f'apart, Delta chi2 {vs["chi2"]:.1f} (night '
                     f'{night["window"]})')
    lines.append('Two visits in a night, where the changes they predict '
                 'differ most (independent of the offset):')
    for ii in pln['best_pair']:
        night = pln['nights'][ii]
        vs = night['versus'][0]['pair']
        lines.append(f'  {night["date"]} {vs["ut1"]} and {vs["ut2"]} UT: '
                     f'{vs["ddv"]:+.1f} m/s, Delta chi2 {vs["chi2"]:.1f}')
    return lines


def plan_figure(pln: Dict[str, Any], sols: Sequence[Dict[str, Any]],
                title: Optional[str] = None):
    """
    The velocities the two most probable solutions predict over the planned
    days (one sigma shaded), the nights the star can be observed, and what
    one visit and a pair in each night would tell them apart by

    :return: matplotlib figure
    """
    from koloa import plotting as kplot
    plt, C = kplot.plt, kplot.C
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)
    start, ndays = pln['start'], pln['ndays']
    tt = np.arange(start, start + ndays, 0.02)
    colours = [C['koloa'], C['outlier']]
    for sol, colour in zip(sols[:2], colours):
        val, spread = _predict(sol, tt, pln['instrument'])
        ax1.plot(tt, val, color=colour, lw=0.9,
                 label=f'{sol["period"]:.4f} d ({label(sol["name"])})')
        ax1.fill_between(tt, val - spread, val + spread, color=colour,
                         alpha=0.2, lw=0)
    for night in pln['nights']:
        for ax in (ax1, ax2):
            ax.axvspan(night['start'], night['end'], color=C['muted'],
                       alpha=0.15, lw=0)
    ax1.set_ylabel('predicted velocity [m/s]')
    ax1.legend(loc='upper right')
    if pln['nights']:
        mid = [0.5 * (night['start'] + night['end'])
               for night in pln['nights']]
        ax2.plot(mid, [night['versus'][0]['single']['chi2']
                       for night in pln['nights']], 'o', ms=4,
                 color=C['koloa'], label='one visit')
        ax2.plot(mid, [night['versus'][0]['pair']['chi2']
                       for night in pln['nights']], 's', ms=4,
                 color=C['text'], label='two visits in the night')
        ax2.legend(loc='upper right')
    ax2.set_ylabel('$\\Delta\\chi^2$ against the alias')
    ax2.set_xlabel('time [BJD - 2400000] (the nights the star is up shaded)')
    ax2.set_ylim(bottom=0)
    if title:
        ax1.set_title(title, loc='left')
    fig.tight_layout()
    return fig


# =============================================================================
# End of code
# =============================================================================

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
What more data would do for a candidate: more velocities from the
instrument that found it, and its astrometry in Gaia DR4.

    from koloa import prospects
    told = prospects.candidate(result, quick, star)

More velocities. The candidate is taken as real, at the period, amplitude
and phase of its fold. 50, 100, 150 or 200 more nights of an instrument of
the files of the star (SPIRou, NIRPS...) are spread over 6 months, a year
or two, on the nights it can be observed from the site of the instrument:
above airmass 2 for an hour of night, and for SPIRou within the bright
half of the lunation, where the telescope has it (86 % of the public
SPIRou velocities of GJ 876, Moutou et al. 2023, are within a week of the
full Moon). Each new night has the mean error bar of the instrument, and
the noise the quick look added to it. The Fisher information of a sinusoid
with an offset per instrument and a trend then says by how much the error
of K (and of the period) shrinks: the forecast is the error measured
today times that ratio, so it keeps whatever today's scatter holds that
the error bars do not. A FIP is said too, to an order of magnitude: its
logarithm falls as (K / sigma_K)^2 / 2 grows.

Astrometry. A planet of minimum mass m sin i at a from a star of mass M at
the distance r moves it by at least

    alpha = (m sin i / M) (a / r) = 95.4 uas (m / M_J) (a / au)
                                             (M_sun / M) (10 pc / r)

(Lammers & Winn 2025, AJ 171, 18, their equations 1 and 2; the true mass,
and the motion, are larger by 1 / sin i). One crossing of a field of view
of Gaia measures a position along its scan to sigma_fov = 54 uas for
G < 14, and 10^(0.2 (G - 14)) times that beyond (their equation 9, fitted
to the precision of Gaia DR3 of Holl et al. 2023; worse again for stars
brighter than G = 6). Their threshold for a detection in DR4 is
alpha / sigma_fov > 1.5, a chi2 lower by 50 with the planet; the older
one is 3 (Casertano et al. 2008, Perryman et al. 2014). DR4 holds 66
months of data: the crossings of a star are those of Gaia DR3 (34 months)
times 66 / 34, and 72 when they are not known (the median of the M dwarfs
within 7 pc). The significance is that chi2, to an order of magnitude:
N crossings lower it by about N (alpha / sigma_fov)^2 / 3 (a crossing
measures one direction, an orbit is seen at some angle), which is their
22 (alpha / sigma_fov)^2 of DR4 (equation 31) at 66 crossings, and 50 at
their threshold. An orbit longer than 4 years is not covered well enough
to be told from a proper motion.

Both are orders of magnitude for planning, not results.

Created on 2026-10-08

@author: artigau
"""
import textwrap
import zlib
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: the nights added, and the spans they are spread over [days]
MORE = (50, 100, 150, 200)
SPANS = ((183.0, '6 months'), (365.25, '1 year'), (730.5, '2 years'))
#: a star is observed above this altitude [deg] (airmass 2), the Sun below
#: that one (nautical twilight), for at least this long [days]
MIN_ALT, SUN_ALT, MIN_UP = 30.0, -12.0, 1.0 / 24
#: the synodic month [days], and a full Moon (JD; 2000-01-21 04:48 UT)
SYNODIC, FULL_MOON = 29.530588861, 2451564.7
#: the instruments scheduled in bright time (by the start of their name,
#: letters only, upper case): within a quarter of a lunation of full Moon
BRIGHT_TIME = ('SPIROU',)
#: the draws of the nights whose forecasts are averaged
DRAWS = 8
#: the noise of one crossing of a field of view of Gaia for a bright star
#: [uas], and the magnitude where the photon noise takes over (Lammers &
#: Winn 2025, from the precision of Gaia DR3 of Holl et al. 2023)
SIGMA_FOV, G_BREAK, G_BRIGHT = 54.0, 14.0, 6.0
#: the motion of the Sun by Jupiter at 1 au seen from 10 pc [uas], and
#: Jupiter in Earth masses
ALPHA_JUP, MJUP_EARTH = 95.4, 317.83
#: the data of Gaia DR4 and DR3 [months]; the crossings of DR4 when those
#: of DR3 are not known (the median of the M dwarfs within 7 pc: 37 in
#: DR3); the longest orbit DR4 covers well [years]
DR4_MONTHS, DR3_MONTHS, DR4_TRANSITS, DR4_PERIOD = 66.0, 34.0, 72, 4.0
#: alpha / sigma_fov for a detection: in DR4 and in DR5 (Lammers & Winn
#: 2025), and the older one (Casertano et al. 2008)
SNR_DR4, SNR_DR5, SNR_CLASSIC = 1.5, 1.0, 3.0
#: the chi2 a planet must lower for a detection, and the share of
#: (alpha / sigma_fov)^2 a crossing lowers it by (Lammers & Winn 2025:
#: 22 (alpha / sigma_fov)^2 in DR4, at about 66 crossings)
DCHI2_LIMIT, DCHI2_SHARE = 50.0, 1.0 / 3.0


# =============================================================================
# Define functions
# =============================================================================
def nights(ra: float, dec: float, site: str, start: float, span: float,
           bright: bool = False, step: float = 20.0 / 1440) -> np.ndarray:
    """
    The nights a star can be observed from a site: above MIN_ALT while the
    Sun is below SUN_ALT, for MIN_UP at least; one time a night, the
    middle of when it is up

    :param ra: float, its right ascension [deg]
    :param dec: float, its declination [deg]
    :param site: str, a key of koloa.aliases.SITES
    :param start: float, the first day [BJD - 2400000]
    :param span: float, how long [days]
    :param bright: bool, only the bright half of each lunation (within a
                   quarter of it of full Moon)
    :param step: float, the time step [days]

    :return: np.ndarray, a time for each night [BJD - 2400000]
    """
    from koloa.aliases import SITES, _sun, altitude
    lat, lon = SITES[site]
    times = np.arange(start, start + span, step)
    jd = times + 2400000.0
    sra, sdec = _sun(jd)
    up = (altitude(sra, sdec, jd, lat, lon) < SUN_ALT) & (
        altitude(ra, dec, jd, lat, lon) > MIN_ALT)
    if bright:
        up &= np.abs(((jd - FULL_MOON + SYNODIC / 2) % SYNODIC)
                     - SYNODIC / 2) < SYNODIC / 4
    idx = np.where(up)[0]
    if not len(idx):
        return np.zeros(0)
    out = []
    for grp in np.split(idx, np.where(np.diff(idx) > 0.125 / step)[0] + 1):
        if len(grp) * step >= MIN_UP:
            out.append(float(np.mean(times[grp])))
    return np.array(out)


def _errors(time: np.ndarray, sigma: np.ndarray, inst: np.ndarray,
            period: float, amp: float, tconj: float, trend: bool = True):
    """the errors of K and of the period of a sinusoid (an offset per
    instrument, a trend) from its Fisher information on a sampling"""
    tref = float(np.mean(time))
    tt = time - tref
    arg = 2 * np.pi * (time - tconj) / period
    cols = [(inst == name).astype(float) for name in np.unique(inst)]
    if trend:
        cols.append(tt / 1000.0)
    # v = -K sin(arg): its derivatives in K and in the period
    cols.append(-np.sin(arg))
    cols.append(amp * np.cos(arg) * 2 * np.pi * (time - tconj) / period ** 2)
    jac = np.array(cols).T
    info = jac.T @ (jac / sigma[:, None] ** 2)
    try:
        cov = np.linalg.inv(info)
    except np.linalg.LinAlgError:
        return np.nan, np.nan
    if not np.all(np.diag(cov) > 0):
        return np.nan, np.nan
    return float(np.sqrt(cov[-2, -2])), float(np.sqrt(cov[-1, -1]))


def more_velocities(fold: Dict[str, Any], instrument: str, site: str,
                    ra: float, dec: float, added: Optional[Dict[str, float]]
                    = None, fip: Optional[float] = None,
                    start: Optional[float] = None, seed: int = 0,
                    trend: bool = True) -> Optional[Dict[str, Any]]:
    """
    What MORE nights of an instrument would do for a signal, taken as
    real at the period, amplitude and phase of its fold (see the module)

    :param fold: dict, the fold of the signal as a quick look keeps it:
                 period, K, K_err, P_err, tc, and instruments (name, time,
                 err, valid: its nightly means)
    :param instrument: str, the instrument that would observe again (one
                       of the fold's)
    :param site: str, its site (a key of koloa.aliases.SITES)
    :param ra: float, the right ascension of the star [deg]
    :param dec: float, its declination [deg]
    :param added: dict or None, the noise added to the errors of each
                  instrument [m/s] (the quick look's)
    :param fip: float or None, the FIP of the signal today
    :param start: float or None, when the new nights begin [BJD - 2400000]
                  (None: now, or the last night of the series if later)
    :param seed: int, the seed of the draws of the nights
    :param trend: bool, a trend is fitted with the signal

    :return: dict or None (the instrument has no night in the fold):
             instrument, site, bright (only bright time), error (the mean
             error bar of a night [m/s]), added, sigma (the noise of a new
             night), n (its nights so far), K, K_err, snr (K / K_err
             today), start, and table: for each number of nights and each
             span, n, span (its words), days, nights (those the star can
             be observed), stacked (more nights asked than there are:
             several visits a night), K_err, snr, P_err and log10_fip
             (None when no night can be had)
    """
    from koloa.aliases import rjd_now
    added = added or {}
    parts = [one for one in fold.get('instruments') or [] if one.get('time')]
    mine = next((one for one in parts if one['name'] == instrument), None)
    if mine is None or not fold.get('K') or not fold.get('K_err'):
        return None
    time = np.concatenate([np.asarray(one['time'], float) for one in parts])
    inst = np.concatenate([np.full(len(one['time']), one['name'])
                           for one in parts])
    sigma = np.concatenate([np.sqrt(np.asarray(one['err'], float) ** 2 + float(
        added.get(one['name']) or 0.0) ** 2) for one in parts])
    # a night the quick look holds for an outlier weighs as little
    valid = np.concatenate([np.clip(np.asarray(
        one.get('valid') or np.ones(len(one['time'])), float), 1e-3, 1.0)
        for one in parts])
    sigma = sigma / np.sqrt(valid)
    period, amp = float(fold['period']), float(fold['K'])
    tconj = float(fold.get('tc') or time.min())
    error = float(np.mean(np.asarray(mine['err'], float)))
    more = float(added.get(instrument) or 0.0)
    new = float(np.hypot(error, more))
    k_now, p_now = _errors(time, sigma, inst, period, amp, tconj, trend)
    if not np.isfinite(k_now):
        return None
    key = ''.join(char for char in instrument.upper() if char.isalpha())
    bright = key.startswith(BRIGHT_TIME)
    start = max(rjd_now(), float(time.max()) + 1.0) if start is None \
        else float(start)
    snr = amp / float(fold['K_err'])
    rng = np.random.default_rng(seed)
    table = []
    for days, words in SPANS:
        there = nights(ra, dec, site, start, days, bright)
        for num in MORE:
            row = dict(n=num, span=words, days=days, nights=int(len(there)),
                       stacked=bool(num > len(there)), K_err=None, snr=None,
                       P_err=None, log10_fip=None)
            table.append(row)
            if not len(there):
                continue
            kerr, perr = [], []
            for _ in range(DRAWS):
                # a visit a night; more visits than nights: several a night
                count = np.full(len(there), num // len(there))
                count[rng.choice(len(there), num % len(there),
                                 replace=False)] += 1
                keep = count > 0
                one = _errors(
                    np.concatenate([time, there[keep]]),
                    np.concatenate([sigma, new / np.sqrt(count[keep])]),
                    np.concatenate([inst, np.full(int(keep.sum()),
                                                  instrument)]),
                    period, amp, tconj, trend)
                kerr.append(one[0])
                perr.append(one[1])
            # the error of today, shrunk as the information grows
            row['K_err'] = float(fold['K_err'] * np.nanmedian(kerr) / k_now)
            row['snr'] = float(amp / row['K_err'])
            if fold.get('P_err') and np.isfinite(p_now) and p_now > 0:
                row['P_err'] = float(fold['P_err'] * np.nanmedian(perr)
                                     / p_now)
            if fip is not None:
                row['log10_fip'] = float(
                    np.log10(max(float(fip), 1e-300))
                    - (row['snr'] ** 2 - snr ** 2) / (2 * np.log(10)))
    return dict(instrument=instrument, site=site, bright=bright, error=error,
                added=more, sigma=new, n=int(len(mine['time'])), K=amp,
                K_err=float(fold['K_err']), snr=float(snr), start=start,
                table=table)


def sigma_fov(gmag: float) -> float:
    """the noise of one crossing of a field of view of Gaia [uas] at a G
    magnitude (Lammers & Winn 2025, their equation 9)"""
    return SIGMA_FOV * max(1.0, 10 ** (0.2 * (float(gmag) - G_BREAK)))


def astrometry(msini: float, period: float, mstar: float, parallax: float,
               gmag: Optional[float] = None, transits: Optional[float]
               = None, ruwe: Optional[float] = None) -> Dict[str, Any]:
    """
    The astrometric signal of a planet, at least, and what Gaia DR4 makes
    of it, to an order of magnitude (see the module)

    :param msini: float, its minimum mass [Earth masses]
    :param period: float, its period [days]
    :param mstar: float, the mass of the star [solar masses]
    :param parallax: float, its parallax [mas]
    :param gmag: float or None, its G magnitude (None: a bright star)
    :param transits: float or None, its crossings in Gaia DR3 (34 months;
                     None: DR4_TRANSITS in DR4)
    :param ruwe: float or None, its RUWE in Gaia DR3 (said, not used)

    :return: dict, a (the orbit [au]), alpha (the motion of the star, at
             least [uas]), sigma_fov (the noise of a crossing [uas]),
             snr (alpha / sigma_fov), transits (in DR4) and whether they
             are from DR3 (counted), dchi2 (the chi2 lower with the
             planet, about transits snr^2 / 3), dchi2_limit (what a
             detection asks for), period_yr, msini, mstar, distance [pc],
             gmag, ruwe, and verdict (what it means, in words)
    """
    years = float(period) / 365.25
    axis = years ** (2.0 / 3.0) * float(mstar) ** (1.0 / 3.0)
    distance = 1000.0 / float(parallax)
    alpha = (ALPHA_JUP * (float(msini) / MJUP_EARTH) * axis / float(mstar)
             / (distance / 10.0))
    noise = sigma_fov(gmag if gmag is not None else G_BREAK)
    snr = alpha / noise
    counted = transits is not None and np.isfinite(float(transits))
    ncross = (float(transits) * DR4_MONTHS / DR3_MONTHS if counted
              else float(DR4_TRANSITS))
    if years > DR4_PERIOD:
        verdict = (f'an orbit of {years:.1f} yr, longer than the '
                   f'{DR4_PERIOD:g} yr that DR4 covers well: an '
                   f'acceleration of the star at best')
    elif snr >= SNR_CLASSIC:
        verdict = 'detectable in DR4 (above the threshold of 3 a crossing)'
    elif snr >= SNR_DR4:
        verdict = (f'at the threshold of DR4 ({SNR_DR4:g} a crossing): '
                   f'marginal')
    elif snr >= SNR_DR5:
        verdict = (f'below the threshold of DR4 ({SNR_DR4:g}), at that of '
                   f'DR5 ({SNR_DR5:g})')
    elif snr >= 0.1:
        verdict = 'not detectable in DR4 (well below its threshold of 1.5)'
    else:
        verdict = 'not detectable by Gaia: far too small'
    if gmag is not None and float(gmag) < G_BRIGHT:
        verdict += f'; brighter than G = {G_BRIGHT:g}: Gaia is worse there'
    return dict(a=float(axis), alpha=float(alpha), sigma_fov=float(noise),
                snr=float(snr), transits=float(ncross), counted=bool(counted),
                dchi2=float(ncross * snr ** 2 * DCHI2_SHARE),
                dchi2_limit=DCHI2_LIMIT, period_yr=float(years),
                msini=float(msini), mstar=float(mstar),
                distance=float(distance),
                gmag=None if gmag is None else float(gmag),
                ruwe=None if ruwe is None else float(ruwe), verdict=verdict)


def _round(value: float) -> str:
    """a number to two figures, in full (380, 0.0031)"""
    return f'{float(f"{float(value):.2g}"):g}'


def candidate(res: Dict[str, Any], quick: Optional[Dict[str, Any]],
              star: Optional[Dict[str, Any]] = None,
              start: Optional[float] = None) -> Optional[Dict[str, Any]]:
    """
    The prospects of the candidate of a star of a batch: what more
    velocities of each instrument of its files would do, and its
    astrometry in Gaia DR4

    :param res: dict, the result of the star (with its reading:
                koloa.batchpdf.reading)
    :param quick: dict or None, its quick FIP (quick.json)
    :param star: dict or None, the star (of a survey: ra, dec, its
                 parallax or distance, G, and gaia: what Gaia DR3 has of
                 it)
    :param start: float or None, when new nights begin (None: now)

    :return: dict or None (no candidate): peak (id, period, K, K_err, fip,
             msini), velocities (a list, from more_velocities, one for
             each instrument of the files; none for a star with archives
             only, or of unknown position or site), astrometry (from
             astrometry(); None without a mass or a parallax) and notes
             (why something is missing)
    """
    from koloa.aliases import site_of
    star = star or {}
    told = res.get('reading') or {}
    result = (quick or {}).get('result') or {}
    peak = next((one for one in told.get('peaks') or []
                 if one.get('counts') and not one.get('what')), None)
    if told.get('kind') != 'candidate' or peak is None:
        return None
    fold = next((one for one in result.get('folds') or []
                 if one.get('id') == peak['id']), None)
    out = dict(peak=dict(id=peak['id'], period=peak['period'], K=peak.get('K'),
                         K_err=peak.get('K_err'), fip=peak.get('fip'),
                         msini=(peak.get('msini') or [None])[0]),
               velocities=[], astrometry=None, notes=[])
    summ = res.get('summary') or {}
    sources = summ.get('sources') or {}
    files = [name for name, where in sources.items()
             if str(where).startswith('file')]
    if fold is None:
        out['notes'].append('no fold of the candidate was kept')
    elif not files:
        out['notes'].append('no file of the star: no instrument to ask for '
                            'more nights of')
    elif star.get('ra') is None or star.get('dec') is None:
        out['notes'].append('the position of the star is not known: no '
                            'observing window')
    else:
        seed = zlib.crc32(str(res.get('name')).encode())
        for name in files:
            site = site_of([name])
            if site is None:
                out['notes'].append(f'{name}: its site is not known')
                continue
            try:
                made = more_velocities(
                    fold, name, site, float(star['ra']), float(star['dec']),
                    result.get('inflation'), peak.get('fip'), start, seed,
                    trend=(result.get('settings') or {}).get('trend', 1) >= 1)
            except Exception as err:  # a forecast, not a need
                made = None
                log(f'prospects, {name}: {type(err).__name__}: {err}', 'warn')
            if made is not None:
                out['velocities'].append(made)
    gaia = star.get('gaia') or {}
    plx = gaia.get('parallax') or star.get('plx') or (
        1000.0 / star['distance'] if star.get('distance') else None)
    mass = (told.get('mass') or [None])[0]
    if not out['peak']['msini'] or not mass or not plx:
        out['notes'].append('no mass of the star, or no parallax: no '
                            'astrometric signal')
    else:
        gmag = gaia.get('G') if gaia.get('G') is not None else star.get('G')
        out['astrometry'] = astrometry(
            out['peak']['msini'], peak['period'], mass, plx, gmag,
            gaia.get('transits'), gaia.get('ruwe'))
    return out


def lines(told: Optional[Dict[str, Any]]) -> List[str]:
    """the prospects of a candidate, in words: a table of K over its error
    for each instrument (and the FIP that goes with it, to an order of
    magnitude), and what Gaia DR4 makes of its astrometric signal"""
    if not told:
        return []
    out = []
    for one in told.get('velocities') or []:
        out.append(
            f'more {one["instrument"]} from {one["site"]}'
            + (', in bright time' if one['bright'] else '')
            + f': {one["sigma"]:.2f} m/s a night (its mean error bar'
            + (f' {one["error"]:.2f}, and {one["added"]:.2f} added'
               if one['added'] >= 0.005 else '')
            + f'), {one["n"]} nights so far.')
        out.append(f'   K over its error is {one["snr"]:.1f} today; with')
        out.append(f'   {"more nights":>11s}' + ''.join(
            f'{"over " + words:>25s}' for _, words in SPANS))
        for num in MORE:
            cells = []
            for _, words in SPANS:
                row = next(row for row in one['table']
                           if row['n'] == num and row['span'] == words)
                if row['snr'] is None:
                    cells.append(f'{"not observable":>25s}')
                    continue
                fip = row['log10_fip']
                cells.append(
                    f'{row["snr"]:9.1f}{"*" if row["stacked"] else " "}'
                    + (f'(FIP 1e{max(fip, -300):+.0f})' if fip is not None
                       else '').rjust(15))
            out.append(f'   {"+" + str(num):>11s}' + ''.join(cells))
        if any(row['stacked'] for row in one['table']):
            out.append('   * more nights asked than the star has in that '
                       'time: several visits a night')
        out.append('')
    for note in told.get('notes') or []:
        out.append(note)
    ast = told.get('astrometry')
    if ast:
        out.append(
            f'Gaia DR4: the star moves by at least {ast["alpha"]:.3g} uas '
            f'(m sin i = {ast["msini"]:.3g} ME at {ast["a"]:.3g} au from '
            f'{ast["mstar"]:.2f} Msun, at {ast["distance"]:.1f} pc).')
        out.append(
            f'   One crossing measures {ast["sigma_fov"]:.0f} uas'
            + (f' at G = {ast["gmag"]:.1f}' if ast['gmag'] is not None
               else '') + f': S/N {ast["snr"]:.2g} a crossing; about '
            f'{ast["transits"]:.0f} crossings in DR4'
            + (' (those of DR3 x 66/34).' if ast['counted']
               else ' (a typical number).'))
        out.append(
            f'   Significance: a chi2 lower by about {_round(ast["dchi2"])} '
            f'with the planet, where a detection asks for '
            f'{ast.get("dchi2_limit", DCHI2_LIMIT):g}.')
        out.extend(textwrap.wrap(
            ast['verdict'][:1].upper() + ast['verdict'][1:] + '.'
            + (f' RUWE {ast["ruwe"]:.2f} in DR3.' if ast.get('ruwe')
               else ''), 112, initial_indent='   ',
            subsequent_indent='   '))
    return out


def short(told: Optional[Dict[str, Any]]) -> str:
    """the prospects of a candidate in one line: K over its error with 100
    more nights over a year, and the S/N of a crossing of Gaia"""
    if not told:
        return ''
    parts = []
    for one in (told.get('velocities') or [])[:2]:
        row = next((row for row in one['table'] if row['n'] == 100
                    and 300 < row['days'] < 400 and row['snr'] is not None),
                   None)
        if row is not None:
            parts.append(f'100 more nights of {one["instrument"]} over a '
                         f'year: K over its error {one["snr"]:.1f} to '
                         f'{row["snr"]:.1f}')
    ast = told.get('astrometry')
    if ast:
        parts.append(f'Gaia DR4: {ast["alpha"]:.2g} uas at least, S/N '
                     f'{ast["snr"]:.2g} a crossing')
    return '; '.join(parts)

# =============================================================================
# End of code
# =============================================================================

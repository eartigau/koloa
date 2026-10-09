#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The summary of a batch of stars (koloa.survey): what to think of each star
on one page, and one PDF of them all.

    from koloa import batchpdf
    batchpdf.summary('m_dwarfs_15pc')     # m_dwarfs_15pc/results/summary.pdf

A batch keeps, of each star, the line of its table (result.json) and its
quick FIP as the page draws it (quick.json). reading() says what its peaks
are: a known planet (the period of one of the NASA Exoplanet Archive's),
the rotation of the star or one of its harmonics (the periods SIMBAD and
CARMENES DR1 list), a drift (a period about as long as the series), a
year, or none of these: a candidate. figure() draws the page of a star
from the same two files, with nothing computed again: its verdict, its
peaks, its FIP, its folds, its velocities, its datasets and what the rules
left out, its acceleration, the transit looked for, its rotation, and its
detailed report when one was made. A candidate has a second page, its
prospects (koloa.prospects): what more nights of the instrument of its
files would do for it, and its astrometric signal in Gaia DR4. summary()
puts the table of the stars and the pages of each in one PDF, the
candidates first.

A verdict is a reading of a quick look, without a GP: it sorts the stars,
it does not decide on a planet.

Created on 2026-10-07

@author: artigau
"""
import json
import os
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: a peak counts when its FIP (of the period or any of its aliases) is
#: below this
FIP_LIMIT = 0.01
#: a peak is a known planet when within this fraction of its period
KNOWN_MATCH = 0.01
#: a peak sits at the rotation, or at one of its harmonics, when within
#: this fraction of it (as koloa.detailed)
ROTATION_MATCH = 0.03
#: the harmonics of the rotation looked at (P, P/2, P/3)
ROTATION_HARMONICS = (1, 2, 3)
#: a period longer than this fraction of the series is a drift
DRIFT = 0.5
#: a peak is at a year when within this fraction of it
YEAR, YEAR_MATCH = 365.25, 0.05
#: the verdicts, as they are sorted and coloured
KINDS = ('candidate', 'known', 'drift', 'nothing', 'failed')
COLOURS = dict(candidate='#1baf7a', known='#2a78d6', drift='#eda100',
               nothing='#8a8f98', failed='#e34948')
WORDS = dict(candidate='candidate', known='known planet', drift='drift',
             nothing='nothing', failed='failed')
#: the lines of the table of the stars on a page of the summary
ROWS = 52


# =============================================================================
# Define functions
# =============================================================================
def _fip(value: Any) -> str:
    """a FIP as it is said: one of 0 is below what a number holds"""
    if value is None:
        return ''
    return '< 1e-300' if not value else f'{value:.1e}'


def _mass(star: Dict[str, Any], res: Dict[str, Any]):
    """the mass of a star for a minimum mass: rough, from its spectral
    type (None when it has none)"""
    from koloa.stars import mass_from_spectral_type
    sptype = star.get('sptype') or (res.get('star') or {}).get('spt')
    try:
        return mass_from_spectral_type(sptype) if sptype else None
    except (ValueError, TypeError):
        return None


def reading(res: Dict[str, Any], quick: Optional[Dict[str, Any]],
            star: Optional[Dict[str, Any]] = None,
            limit: float = FIP_LIMIT) -> Dict[str, Any]:
    """
    What the quick look of a star says: each of its numbered peaks with
    what it is, and a verdict

    :param res: dict, the result of the star (its result.json)
    :param quick: dict or None, its quick FIP (quick.json, or its 'result')
    :param star: dict or None, the star (of targets.json: its rotation
                 periods, its spectral type)
    :param limit: float, the FIP below which a peak counts

    :return: dict, kind (candidate: a peak that counts and is neither a
             known planet nor a drift; known: the peaks that count are
             known planets; drift: they are as long as the series;
             nothing: no peak counts; failed), line (the verdict in
             words), peaks (id, period, fip, alone, K, K_err, msini, what:
             known, drift or '', notes, counts), limit, mass
    """
    star = star or {}
    result = (quick or {}).get('result') or quick or {}
    summ = res.get('summary') or {}
    if res.get('status') != 'done' or not result.get('peak_list'):
        failed = res.get('status') != 'done'
        return dict(kind='failed' if failed else 'nothing', limit=limit,
                    line=(f'failed: {res.get("error") or res.get("status")}'
                          if failed else 'no peak in its quick FIP'),
                    peaks=[], mass=None)
    base = float(summ.get('baseline') or 0.0)
    known = result.get('known') or []
    folds = {fold['id']: fold for fold in result.get('folds') or []
             if fold.get('id') is not None}
    spins = [one for one in star.get('rotation') or [] if one.get('period')]
    mass = _mass(star, res)
    peaks = []
    for peak in result['peak_list']:
        if not peak.get('named', True):
            continue
        per = float(peak['period'])
        fold = folds.get(peak['id']) or {}
        what, notes = '', []
        for pl in known:
            if pl.get('P') and abs(per / float(pl['P']) - 1) < KNOWN_MATCH:
                what = 'known'
                notes.append(f'{pl["name"]}, {float(pl["P"]):.4g} d')
                break
        if not what and base and per > DRIFT * base:
            what = 'drift'
            notes.append(f'as long as the series ({base:.0f} d): a drift')
        for spin in spins:
            for harm in ROTATION_HARMONICS:
                if abs(float(spin['period']) / harm / per - 1) \
                        < ROTATION_MATCH:
                    notes.append(('at the rotation' if harm == 1 else
                                  f'at the rotation / {harm}')
                                 + f' ({float(spin["period"]):.4g} d)')
                    break
            else:
                continue
            break
        if abs(per / YEAR - 1) < YEAR_MATCH:
            notes.append('at a year')
        msini = None
        if mass and fold.get('K'):
            from koloa.stars import minimum_mass
            try:
                msini = [float(val) for val in minimum_mass(
                    float(fold['K']), per, 0.0, mass[0],
                    amp_err=float(fold.get('K_err') or 0.0),
                    mstar_err=mass[1], ndraw=1000)['earth']]
            except (ValueError, ZeroDivisionError, FloatingPointError):
                msini = None
        peaks.append(dict(
            id=int(peak['id']), period=per, fip=peak.get('family'),
            alone=peak.get('alone'), K=fold.get('K'),
            K_err=fold.get('K_err'), msini=msini, what=what, notes=notes,
            counts=bool(peak.get('family') is not None
                        and peak['family'] < limit)))
    counted = [peak for peak in peaks if peak['counts']]
    fresh = [peak for peak in counted if not peak['what']]
    if fresh:
        best = fresh[0]
        kind = 'candidate'
        line = (f'candidate: #{best["id"]} at {best["period"]:.4f} d, FIP '
                f'{_fip(best["fip"])}'
                + (f', K = {best["K"]:.2f} m/s' if best.get('K') else '')
                + (f' ({"; ".join(best["notes"])})' if best['notes'] else '')
                + (f'; {len(fresh) - 1} more' if len(fresh) > 1 else ''))
    elif any(peak['what'] == 'known' for peak in counted):
        kind = 'known'
        line = 'known: ' + '; '.join(
            f'#{peak["id"]} is {peak["notes"][0]}' for peak in counted
            if peak['what'] == 'known')
    elif counted:
        kind = 'drift'
        line = (f'a drift: #{counted[0]["id"]} at '
                f'{counted[0]["period"]:.0f} d, as long as the series')
    else:
        kind = 'nothing'
        best = peaks[0] if peaks else None
        line = (f'no peak with a FIP below {limit:g}'
                + (f' (the best, {best["period"]:.4f} d: {_fip(best["fip"])})'
                   if best else ''))
    return dict(kind=kind, line=line, peaks=peaks, limit=limit,
                mass=list(mass) if mass else None)


def _load(root: str, name: str):
    """the result and the quick FIP of a star of a batch folder"""
    from koloa import survey
    path = survey.result_path(root, name)
    if not os.path.exists(path):
        return None, None
    with open(path) as handle:
        res = json.load(handle)
    quick = None
    kept = os.path.join(os.path.dirname(path), 'quick.json')
    if os.path.exists(kept):
        with open(kept) as handle:
            quick = json.load(handle)
    return res, quick


def _short(text: Any, width: int) -> str:
    """a line that fits"""
    text = ' '.join(str(text).split())
    return text if len(text) <= width else text[:width - 3] + '...'


def figure(res: Dict[str, Any], quick: Optional[Dict[str, Any]],
           star: Optional[Dict[str, Any]] = None, batch: str = '',
           limit: float = FIP_LIMIT):
    """
    The page of a star: its verdict, its peaks, its FIP, its folds, its
    velocities, and what was done with its data; drawn from its result and
    its quick FIP as a batch keeps them, nothing computed

    :param res: dict, the result of the star (result.json)
    :param quick: dict or None, its quick FIP (quick.json)
    :param star: dict or None, the star (of targets.json)
    :param batch: str, the name of the batch, for the corner of the page
    :param limit: float, the FIP below which a peak counts

    :return: matplotlib Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from koloa import plotting as kplot
    kplot.set_style('paper')
    star = star or {}
    told = res.get('reading') or reading(res, quick, star, limit)
    result = (quick or {}).get('result') or {}
    summ = res.get('summary') or {}
    name = str(res.get('name') or star.get('name') or '')
    fig = plt.figure(figsize=(7.8, 10.6))
    mono = dict(family='monospace', fontsize=7.4, va='top')
    # the star
    fig.text(0.06, 0.972, name, fontsize=18, fontweight='bold', va='top')
    main = star.get('main') or (res.get('star') or {}).get('raw') or ''
    dist = star.get('distance') or res.get('distance')
    parts = [str(main) if main and str(main) != name else '',
             str(star.get('sptype') or res.get('sptype') or ''),
             f'{dist:.2f} pc' if dist else '',
             (f'about {told["mass"][0]:.2f} solar masses (its type)'
              if told.get('mass') else '')]
    fig.text(0.06, 0.942, '  ·  '.join(part for part in parts if part),
             fontsize=9, color='0.25', va='top')
    fig.text(0.94, 0.972, _short(batch, 40), fontsize=8.5, color='0.4',
             ha='right', va='top')
    fig.text(0.94, 0.955, f'koloa, {res.get("made") or ""}', fontsize=7.5,
             color='0.5', ha='right', va='top')
    # the verdict, in its colour
    box = fig.add_axes([0.06, 0.893, 0.88, 0.032])
    box.axis('off')
    box.add_patch(Rectangle((0, 0), 1, 1, transform=box.transAxes,
                            facecolor=COLOURS[told['kind']], alpha=0.22,
                            edgecolor=COLOURS[told['kind']], lw=1.0))
    box.text(0.012, 0.5, _short(told['line'], 118), fontsize=9.6,
             fontweight='bold', va='center', transform=box.transAxes)
    # its peaks
    top = 0.884
    fig.text(0.06, top, f'{"":3s} {"P [d]":>10s} {"FIP (P or alias)":>17s} '
             f'{"FIP (P alone)":>14s} {"K [m/s]":>14s} {"m sin i [ME]":>13s}'
             f'  what', color='0.35', **mono)
    for it, peak in enumerate(told['peaks'][:5]):
        amp = (f'{peak["K"]:.2f} ± {peak["K_err"]:.2f}'
               if peak.get('K') is not None else '')
        mass = f'{peak["msini"][0]:.1f}' if peak.get('msini') else ''
        what = '; '.join(peak['notes']) or (
            'none of these: a candidate' if peak['counts'] else
            f'FIP above {told["limit"]:g}')
        fig.text(0.06, top - 0.0125 * (it + 1),
                 f'#{peak["id"]:<2d} {peak["period"]:10.4f} '
                 f'{_fip(peak["fip"]):>17s} {_fip(peak["alone"]):>14s} '
                 f'{amp:>14s} {mass:>13s}  {_short(what, 46)}',
                 fontweight='bold' if peak['counts'] else 'normal',
                 color=kplot.C['text'] if peak['counts'] else '0.4', **mono)
    insts = list((result.get('instruments') or {}).keys())
    colour = {inst: kplot.INST_COLOURS[it % 8] for it, inst in
              enumerate(insts)}
    marker = {inst: kplot.INST_MARKERS[it % 8] for it, inst in
              enumerate(insts)}
    # its FIP
    if result.get('period'):
        ax = fig.add_axes([0.09, 0.612, 0.85, 0.180])
        per = np.asarray(result['period'], dtype=float)
        ax.plot(per, result['alone'], color='0.65', lw=0.7,
                label='the period alone')
        ax.plot(per, result['family'], color=kplot.C['text'], lw=0.9,
                label='the period or any of its aliases')
        ax.axhline(-np.log10(told['limit']), color='0.3', lw=0.8, ls=':')
        for pl in result.get('known') or []:
            if pl.get('P') and per.min() <= pl['P'] <= per.max():
                ax.axvline(pl['P'], color=kplot.C['outlier'], lw=0.8, ls='--')
                ax.text(pl['P'], 0.97, f' {pl["name"].split()[-1]}',
                        fontsize=6.5, color=kplot.C['outlier'], va='top',
                        transform=ax.get_xaxis_transform())
        for spin in (star.get('rotation') or [])[:1]:
            for harm in (1, 2):
                val = float(spin['period']) / harm
                if per.min() <= val <= per.max():
                    ax.axvline(val, color='#4a3aa7', lw=0.8, ls='-.')
                    ax.text(val, 0.72, ' P_rot' + (f'/{harm}' if harm > 1
                                                   else ''),
                            fontsize=6.5, color='#4a3aa7', va='top',
                            transform=ax.get_xaxis_transform())
        high = max(max(result['family']), 2.4)
        for peak in told['peaks']:
            height = -np.log10(max(peak['fip'] or 1e-300, 1e-300))
            ax.annotate(f'#{peak["id"]}', (peak['period'],
                                           min(height, high)),
                        xytext=(0, 3), textcoords='offset points',
                        ha='center', fontsize=7.5, fontweight='bold')
        ax.set_xscale('log')
        ax.set_xlim(per.min(), per.max())
        ax.set_ylim(0, 1.12 * high)
        ax.set_xlabel('period [d]', fontsize=8)
        ax.set_ylabel('-log10 FIP', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6.5, loc='upper right', frameon=False)
    # its folds, the two best
    folds = [fold for fold in result.get('folds') or []
             if fold.get('instruments')]
    for it, fold in enumerate(folds[:2]):
        ax = fig.add_axes([0.09 + 0.46 * it, 0.418, 0.39, 0.142])
        for one in fold['instruments']:
            ax.errorbar(one['phase'], one['rv'], one['err'], fmt='none',
                        ecolor=colour.get(one['name'], '0.5'), elinewidth=0.4,
                        alpha=0.5, rasterized=True)
            ax.plot(one['phase'], one['rv'], ls='none', ms=2.4,
                    marker=marker.get(one['name'], 'o'),
                    color=colour.get(one['name'], '0.5'), rasterized=True)
        curve = fold.get('curve') or {}
        if curve.get('phase'):
            ax.plot(curve['phase'], curve['rv'], color=kplot.C['text'],
                    lw=1.1)
        ax.set_title(f'#{fold.get("id")}: P = {fold["period"]:.4f} d, K = '
                     f'{fold["K"]:.2f} ± {fold["K_err"]:.2f} m/s',
                     fontsize=7.5)
        ax.set_xlabel('phase (0: conjunction)', fontsize=7.5)
        if it == 0:
            ax.set_ylabel('RV [m/s]', fontsize=8)
        ax.tick_params(labelsize=7)
    # its velocities, the offsets and the trend of the fold taken out
    if folds:
        ax = fig.add_axes([0.09, 0.238, 0.85, 0.125])
        for one in folds[0]['instruments']:
            ax.errorbar(one['time'], one['rv'], one['err'], fmt='none',
                        ecolor=colour.get(one['name'], '0.5'), elinewidth=0.4,
                        alpha=0.5, rasterized=True)
            ax.plot(one['time'], one['rv'], ls='none', ms=2.4,
                    marker=marker.get(one['name'], 'o'),
                    color=colour.get(one['name'], '0.5'),
                    label=_short(one['name'], 26), rasterized=True)
        ax.set_xlabel('BJD - 2400000', fontsize=7.5)
        ax.set_ylabel('RV [m/s]', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=5.8, ncol=min(4, max(1, len(insts))),
                  loc='lower center', bbox_to_anchor=(0.5, 1.0),
                  frameon=False, handletextpad=0.2, columnspacing=0.9)
    # what was done with its data, and what else is known
    lines = []
    nights = summ.get('instruments') or {}
    sets = res.get('datasets') or {}
    if summ:
        lines.append(f'{summ.get("n")} nights over '
                     f'{summ.get("baseline") or 0:.0f} d'
                     + (f', {sets.get("used")} datasets used of '
                        f'{sets.get("all")}' if sets else '') + ': '
                     + ', '.join(f'{inst} {count}' for inst, count
                                 in nights.items()))
    for line in (sets.get('told') or [])[:3]:
        lines.append('  ' + line)
    if len(sets.get('told') or []) > 3:
        lines.append(f'  and {len(sets["told"]) - 3} more datasets left '
                     f'out or cut (result.json)')
    if summ.get('accel') is not None:
        lines.append(f'acceleration of the star: {summ["accel"]:+.3f} ± '
                     f'{summ.get("accel_err") or 0:.3f} m/s/yr '
                     f'({summ.get("accel_sigma") or 0:.1f} sigma)')
    spins = star.get('rotation') or []
    lines.append('rotation: ' + ('; '.join(
        f'{float(one["period"]):.4g} d ({one.get("source")})'
        for one in spins[:4]) if spins else 'none published (SIMBAD, '
        'CARMENES DR1)'))
    trans = summ.get('transit') or {}
    if trans:
        lines.append('transit at the best peak: ' + str(trans.get('status'))
                     + (f', {trans["snr"]:.1f} sigma' if trans.get('snr')
                        is not None else '')
                     + (f' ({_short(trans.get("why") or "", 60)})'
                        if trans.get('why') else ''))
    planets = (res.get('star') or {}).get('planets') or []
    if planets:
        lines.append('known planets: ' + ', '.join(
            f'{pl["name"].split()[-1]} {pl["P"]:.4g} d' for pl in planets))
    rep = res.get('report') or {}
    if rep.get('status') and rep['status'] != 'none':
        lines.append(f'detailed report: {rep["status"]}'
                     + (f', GP: {rep["gp"]}' if rep.get('gp') else '')
                     + (f' ({rep["where"]})' if rep.get('where') else ''))
        for one in (rep.get('signals') or [])[:4]:
            lines.append('  ' + _short(one, 120))
        if rep.get('error'):
            lines.append('  ' + _short(rep['error'], 120))
    elif rep.get('why'):
        lines.append(f'detailed report: none ({rep["why"]})')
    if res.get('prospects'):
        # in short; its page follows
        from koloa import prospects
        lines.append('prospects (next page): ' + prospects.short(
            res['prospects']))
    if res.get('note'):
        lines.append('note: ' + str(res['note']))
    for it, line in enumerate(lines[:13]):
        lead = '  ' if line.startswith('  ') else ''
        fig.text(0.06, 0.196 - 0.0125 * it, lead + _short(line, 126), **mono)
    fig.text(0.06, 0.018, 'A reading of a quick look, with no GP: it sorts '
             'the stars, it does not decide on a planet. Each line is to be '
             'looked at.', fontsize=6.8, color='0.45')
    return fig


def prospects_figure(res: Dict[str, Any], star: Optional[Dict[str, Any]]
                     = None, batch: str = ''):
    """
    The page of the prospects of a candidate (koloa.prospects): how K over
    its error grows with 50 to 200 more nights of each instrument of the
    files of the star, spread over 6 months, a year or two, and its
    astrometric signal against what Gaia DR4 measures

    :param res: dict, the result of the star, with its prospects
    :param star: dict or None, the star (of targets.json)
    :param batch: str, the name of the batch, for the corner of the page

    :return: matplotlib Figure, or None (no prospects)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from koloa import plotting as kplot
    from koloa import prospects
    told = res.get('prospects')
    if not told:
        return None
    kplot.set_style('paper')
    star = star or {}
    name = str(res.get('name') or star.get('name') or '')
    peak = told['peak']
    fig = plt.figure(figsize=(7.8, 10.6))
    fig.text(0.06, 0.972, name, fontsize=18, fontweight='bold', va='top')
    fig.text(0.06, 0.942, f'what more data would do for #{peak["id"]}, '
             f'{peak["period"]:.4f} d'
             + (f', K = {peak["K"]:.2f} ± {peak["K_err"]:.2f} m/s'
                if peak.get('K') is not None else '')
             + ', if it is real', fontsize=9, color='0.25', va='top')
    fig.text(0.94, 0.972, _short(batch, 40), fontsize=8.5, color='0.4',
             ha='right', va='top')
    shown = (told.get('velocities') or [])[:2]
    styles = (('-', 'o'), ('--', 's'), (':', 'D'))
    for it, one in enumerate(shown):
        ax = fig.add_axes([0.10 + 0.47 * it, 0.62, 0.37, 0.26])
        colour = kplot.INST_COLOURS[it % 8]
        for (_, words), (line, mark) in zip(prospects.SPANS, styles):
            rows = [row for row in one['table'] if row['span'] == words
                    and row['snr'] is not None]
            if not rows:
                continue
            ax.plot([0] + [row['n'] for row in rows],
                    [one['snr']] + [row['snr'] for row in rows], ls=line,
                    marker=mark, ms=4, color=colour, lw=1.2,
                    label=f'over {words} ({rows[0]["nights"]} nights there)')
        ax.axhline(one['snr'], color='0.5', lw=0.8, ls=':')
        ax.set_xlabel(f'more nights of {one["instrument"]}', fontsize=8)
        if it == 0:
            ax.set_ylabel('K / its error', fontsize=8)
        ax.set_title(f'{one["instrument"]} ({one["site"]}'
                     + (', bright time' if one['bright'] else '') + f'): '
                     f'{one["sigma"]:.2f} m/s a night, {one["n"]} so far',
                     fontsize=7.5)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6.5, loc='lower right', frameon=False)
    mono = dict(family='monospace', fontsize=7.2, va='top')
    top = 0.555 if shown else 0.90
    # as they are written: their columns are aligned
    for it, line in enumerate(prospects.lines(told)[:36]):
        fig.text(0.06, top - 0.0128 * it, line[:128], **mono)
    fig.text(0.06, 0.030, 'More velocities: the candidate taken as real; the '
             'error of K of today shrunk as the Fisher information of a '
             'sinusoid grows with the new nights', fontsize=6.8, color='0.45')
    fig.text(0.06, 0.018, '(the mean error bar of the instrument, the star '
             'above airmass 2 at night). Gaia: Lammers & Winn 2025. Orders '
             'of magnitude, to plan with.', fontsize=6.8, color='0.45')
    return fig


def _table(rows: Sequence[Dict[str, Any]], held: Dict[str, Any],
           counts: Dict[str, int], first: int, npage: int, page: int):
    """a page of the table of the stars of a batch"""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(7.8, 10.6))
    fig.text(0.06, 0.972, f'koloa: the batch {held.get("name") or ""}',
             fontsize=15, fontweight='bold', va='top')
    fig.text(0.94, 0.972, f'page {page} of {npage}', fontsize=8, color='0.4',
             ha='right', va='top')
    fig.text(0.06, 0.945, (
        f'{counts["n"]} stars done of {len(held["targets"])}: '
        + ', '.join(f'{counts[kind]} {WORDS[kind]}' for kind in KINDS
                    if counts[kind])
        + (f'; {counts["reports"]} detailed reports' if counts['reports']
           else '')), fontsize=8.5, color='0.25', va='top')
    mono = dict(family='monospace', fontsize=6.9, va='top')
    fig.text(0.06, 0.922, f'{"":4s} {"star":20s} {"type":8s} {"d [pc]":>6s} '
             f'{"nights":>6s} {"sets":>5s} {"best P [d]":>11s} {"FIP":>9s} '
             f'{"K [m/s]":>8s}  {"verdict":13s} report', color='0.35', **mono)
    for it, row in enumerate(rows):
        res, told = row['res'], row['told']
        summ = res.get('summary') or {}
        sets = res.get('datasets') or {}
        best = next((peak for peak in told['peaks'] if peak['counts']
                     and (told['kind'] != 'candidate' or not peak['what'])),
                    told['peaks'][0] if told['peaks'] else {})
        dist = res.get('distance')
        fig.text(0.06, 0.908 - 0.0166 * it, (
            f'{first + it + 1:4d} {_short(res["name"], 20):20s} '
            f'{_short(res.get("sptype") or "", 8):8s} '
            f'{f"{dist:6.2f}" if dist else "":>6s} '
            f'{str(summ.get("n") or ""):>6s} '
            f'{f"{sets.get("used")}/{sets.get("all")}" if sets else "":>5s} '
            f'{f"{best["period"]:11.4f}" if best else "":>11s} '
            f'{_fip(best.get("fip")) if best else "":>9s} '
            f'{f"{best["K"]:8.2f}" if best.get("K") is not None else "":>8s}  '
            f'{WORDS[told["kind"]]:13s} '
            f'{(res.get("report") or {}).get("status") or ""}'),
            color=COLOURS[told['kind']] if told['kind'] in (
                'candidate', 'failed') else '0.1',
            fontweight='bold' if told['kind'] == 'candidate' else 'normal',
            **mono)
    fig.text(0.06, 0.018, 'The candidates first, then the known planets, '
             'the drifts and the stars with nothing, each by its FIP; the '
             'page of each star follows, in this order.', fontsize=6.8,
             color='0.45')
    return fig


def summary(root: str, limit: Optional[float] = None,
            path: Optional[str] = None) -> Optional[str]:
    """
    The summary of a batch folder as one PDF: the table of its stars done,
    the candidates first, then the page of each (figure()); from what the
    batch kept, nothing computed

    :param root: str, the batch folder
    :param limit: float or None, the FIP below which a peak counts (None:
                  as each star was read when it was computed, else
                  FIP_LIMIT)
    :param path: str or None, the PDF (None: results/summary.pdf)

    :return: str or None, the PDF (None when no star is done)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from koloa import survey
    root = os.path.abspath(os.path.expanduser(root))
    held = survey.targets_of(root)
    rows = []
    for star in held['targets']:
        res, quick = _load(root, star['name'])
        if res is None:
            continue
        told = (res.get('reading') if limit is None and res.get('reading')
                else reading(res, quick, star, limit or FIP_LIMIT))
        rows.append(dict(star=star, res=dict(res, reading=told), quick=quick,
                         told=told))
    if not rows:
        log('summary: no star done yet', 'warn')
        return None

    def rank(row):
        told = row['told']
        best = next((peak for peak in told['peaks'] if peak['counts']
                     and (told['kind'] != 'candidate' or not peak['what'])),
                    told['peaks'][0] if told['peaks'] else None)
        return (KINDS.index(told['kind']),
                best['fip'] if best and best['fip'] is not None else 1.0)
    rows.sort(key=rank)
    counts = {kind: sum(row['told']['kind'] == kind for row in rows)
              for kind in KINDS}
    counts.update(n=len(rows), reports=sum(
        (row['res'].get('report') or {}).get('status') == 'done'
        for row in rows))
    path = path or os.path.join(root, 'results', 'summary.pdf')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    npage = -(-len(rows) // ROWS)
    with PdfPages(path + '.part') as book:
        for page in range(npage):
            fig = _table(rows[page * ROWS:(page + 1) * ROWS], held, counts,
                         page * ROWS, npage, page + 1)
            book.savefig(fig)
            plt.close(fig)
        for row in rows:
            fig = figure(row['res'], row['quick'], row['star'],
                         batch=str(held.get('name') or ''),
                         limit=row['told']['limit'])
            book.savefig(fig, dpi=160)
            plt.close(fig)
            # a candidate: what more data would do for it
            fig = prospects_figure(row['res'], row['star'],
                                   batch=str(held.get('name') or ''))
            if fig is not None:
                book.savefig(fig, dpi=160)
                plt.close(fig)
    os.replace(path + '.part', path)
    log(f'summary: {len(rows)} stars ('
        + ', '.join(f'{counts[kind]} {WORDS[kind]}' for kind in KINDS
                    if counts[kind]) + f') in {path}', 'value')
    return path

# =============================================================================
# End of code
# =============================================================================

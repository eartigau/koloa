#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The duck test as a PDF: every check with its numbers, the figures behind
them, and what the NASA Exoplanet Archive knows of the star.

    report = koloa.duck_test(data, period=fip.best()['period'], fipres=fip,
                             pdf='duck.pdf')             # or, afterwards,
    report.pdf('duck.pdf', data, fipres=fip, target='GJ 687')

The pages:
1. the verdict and every check (status, summary, its numbers), the orbit at
   the period, and whether a known planet has that period; planet or no
   planet is decided on the FIP of the period OR any of its aliases, and
   which alias it is comes apart: the velocities folded at the period and
   at each alias that holds some of the probability, with its share, and
   when the period is ambiguous, the nights that would lift the alias
   (koloa.aliases);
2. the signal: the orbit folded at the period, the FIP periodogram (the
   period marked, the known planets dashed);
3. robustness and coherence: the jackknife (which visit holds the peak
   up) and the amplitude and phase through the campaign;
4. the activity indicators, at the period and its multiples, and the
   TESS light curves of the star (koloa.tess): the periodogram of each
   sector, with the period, its half, third and double marked;
5. the NASA Exoplanet Archive: the known planets of the star, every
   published solution, and the fitted K against them.

The archive needs the network (the star's name goes to CDS Sesame and the
archive); without it, or with archive=False, that page says so.

With outdir, everything also goes into that folder: the report, each figure
as a PDF of its own, the text of the test (<name>_duck.txt) and a JSON
summary (<name>_duck.json: the checks, the orbit, the archive).

Created on 2026-09-29

@author: artigau
"""
import json
import os
import textwrap
from typing import Any, Dict, List, Optional

import numpy as np

from koloa import plotting as kplot
from koloa.data import RVData
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: an A4 page [inches]
PAGE = (8.27, 11.69)


# =============================================================================
# Define functions
# =============================================================================
def _text_page(title: str, blocks: List[Any]):
    """
    A page of text: a title, then blocks, each a (text, style) pair where
    style is 'body', 'mono', 'head' or a colour key of koloa's palette

    :return: matplotlib figure
    """
    fig = kplot.plt.figure(figsize=PAGE)
    fig.text(0.07, 0.955, title, fontsize=15, weight='bold',
             color=kplot.C['text'], va='top')
    ypos = 0.915
    for text, style in blocks:
        if style == 'head':
            ypos -= 0.012
            fig.text(0.07, ypos, text, fontsize=11.5, weight='bold',
                     color=kplot.C['text'], va='top')
            ypos -= 0.028
            continue
        mono = style == 'mono'
        width = 98 if mono else 105
        colour = kplot.C.get(style, kplot.C['text']) \
            if style not in ('body', 'mono') else kplot.C['text']
        for para in str(text).split('\n'):
            lines = textwrap.wrap(para, width=width) or ['']
            for line in lines:
                if ypos < 0.04:
                    break
                fig.text(0.07, ypos, line, fontsize=8.2 if mono else 9,
                         family='monospace' if mono else None,
                         color=colour, va='top')
                ypos -= 0.0155
        ypos -= 0.006
    return fig


def _fit(data: RVData, period: float, unit: str):
    """a circular orbit at the period, outliers modelled (the phase fold,
    and the K to compare with the archive)"""
    from koloa.fit import RVModel
    try:
        return RVModel(data, [dict(period=period,
                                   period_range=(0.99 * period,
                                                 1.01 * period))],
                       likelihood='mixture', unit=unit).fit(nstart=1,
                                                            quiet=True)
    except Exception as err:  # a report is a help: go on without the fit
        log(f'duck report: no fit at {period:.4f} d ({err})', 'warn')
        return None


def _archive(target: str, orbit: Optional[Dict[str, Any]]):
    """the known planets of the star, and the orbit against them"""
    from koloa.archive import compare, known_planets
    known = known_planets(name=target)
    if orbit is not None:
        orbit = compare(dict(orbit), known)
    return known, orbit


def _json(obj):
    """numbers, arrays and the rest for json.dump"""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    return str(obj)


def duck_pdf(report: Any, data: RVData, path: Optional[str] = None,
             fipres: Optional[Any] = None, target: Optional[str] = None,
             archive: bool = True, unit: str = 'both',
             outdir: Optional[str] = None) -> str:
    """
    A detailed PDF of a duck test (see the module's docstring)

    :param report: DuckReport, from koloa.duck_test
    :param data: RVData, the series tested
    :param path: str, the PDF written
    :param fipres: FIPResult or None, the FIP of the series (its periodogram)
    :param target: str or None, the star's name for the archive (the name
                   of the series when None)
    :param archive: bool, ask the NASA Exoplanet Archive
    :param unit: str, the outlier unit of the fit (both, sequence, point)
    :param outdir: str or None, a folder for everything: the report (there
                   when path is None), each figure as a PDF of its own, the
                   text of the test and a JSON summary

    :return: str, the path of the report
    """
    from matplotlib.backends.backend_pdf import PdfPages
    period = float(report.period)
    target = target or data.name
    safe = ''.join(ch if ch.isalnum() or ch in '-_.' else '_'
                   for ch in target)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        if path is None:
            path = os.path.join(outdir, f'{safe}_duck.pdf')
    if path is None:
        raise ValueError('duck_pdf needs a path or an outdir')
    kplot.set_style('paper')
    fit = _fit(data, period, unit)
    orbit = None
    if fit is not None:
        orb = fit.orbits()[0]
        orbit = {key: [float(val) for val in orb[key]]
                 for key in ('P', 'K', 'e')}
    known, archive_note = None, ''
    if archive:
        try:
            known, orbit = _archive(target, orbit)
        except Exception as err:  # no network, or an unknown name
            archive_note = f'The archive could not be asked ({err}).'
    else:
        archive_note = 'The archive was not asked (archive=False).'
    with PdfPages(path) as book:
        def keep(fig, name):
            """a page of the report, and a PDF of its own in the folder"""
            if outdir:
                fig.savefig(os.path.join(outdir, f'{safe}_{name}.pdf'),
                            bbox_inches='tight')
            book.savefig(fig)

        # 1. the verdict and every check
        insts = ('' if data.instruments == ['inst'] else
                 f' ({", ".join(data.instruments)})')
        blocks = [(f'{target}: {data.n} exposures in {data.nseq} visits over '
                   f'{data.baseline:.0f} d{insts}; the signal tested at '
                   f'P = {period:.5f} d.', 'body'),
                  ('Verdict', 'head'), (report.verdict, 'body'),
                  ('The checks', 'head')]
        for check in report.checks:
            colour = {'pass': 'koloa', 'flag': 'outlier'}.get(
                check['status'], 'muted')
            mark = {'pass': 'PASS', 'flag': 'FLAG'}.get(check['status'],
                                                        'NOTE')
            blocks.append((f'[{mark}] {check["name"]}: {check["summary"]}',
                           colour))
        if orbit is not None:
            blocks += [('The orbit at the period (outliers modelled)', 'head'),
                       (f'P = {orbit["P"][0]:.5f} d, K = {orbit["K"][0]:.2f} '
                        f'+{orbit["K"][2]:.2f}/-{orbit["K"][1]:.2f} m/s '
                        f'(circular fit; e from the shape check '
                        f'above)', 'body')]
        blocks.append(('The NASA Exoplanet Archive', 'head'))
        if known is not None and orbit is not None and 'known' in orbit:
            blocks.append((f'The period is that of {orbit["known"]["name"]}: '
                           + '; '.join(
                               f'K = {sol["K"]} +- {sol["K_err"]:.2f} m/s '
                               f'({sol["reference"]}, {sol["z"]:+.1f} sigma '
                               f'from this fit)'
                               for sol in orbit.get('comparisons', [])),
                           'body'))
        elif known is not None:
            names = ', '.join(f'{pl["name"]} ({pl["P"]:.4f} d)'
                              for pl in known.get('planets', []))
            blocks.append((f'No known planet has this period. Known: '
                           f'{names or "none"} (host '
                           f'{known.get("host") or "not in the archive"}).',
                           'body'))
        else:
            blocks.append((archive_note, 'body'))
        keep(_text_page('Duck test', blocks), 'summary')
        # 1b. which alias: the folds, their shares, and the plan
        sols = report.details.get('alias_solutions')
        if sols:
            from koloa import aliases as kal
            keep(kal.figure(sols, title=f'{target}: the velocities folded '
                                        f'at the period and its aliases'),
                 'aliases')
            fam = report.details.get('family_fip')
            blocks = []
            if fam is not None:
                blocks.append((f'Planet or no planet is decided on the '
                               f'period OR any of its aliases: FIP = '
                               f'{fam:.2e}, the probability that none of '
                               f'them holds a signal. Which of them it is: '
                               f'the share of the probability each holds '
                               f'(the FIP), and the log posterior of an '
                               f'orbit fitted at each.', 'body'))
            blocks += [('The period and its aliases', 'head')]
            blocks += [(line, 'mono') for line in kal.describe(sols)]
            pln = report.details.get('alias_plan')
            if pln:
                blocks += [('Lifting the alias', 'head')]
                blocks += [(line, 'mono' if line.startswith('  ') else 'body')
                           for line in kal.plan_text(pln)]
            keep(_text_page('Which period?', blocks), 'aliases_text')
            if pln and pln.get('nights'):
                keep(kal.plan_figure(pln, sols, title=f'{target}: when to '
                                                      f'observe to lift the '
                                                      f'alias'), 'plan')
        # 2. the signal
        if fit is not None:
            keep(kplot.phase(fit, level='point',
                             title=f'{target}: folded at {period:.5f} d'),
                 'phase')
        if fipres is not None:
            marks = [period] + [pl['P'] for pl in
                                (known or {}).get('planets', [])
                                if pl.get('P')]
            keep(kplot.periodograms(
                fipres.freq, fips=dict(koloa=fipres), mark=marks,
                title=f'{target}: FIP (the tested period, and the known '
                      f'planets, dashed)'), 'fip')
        # 3. robustness and coherence
        jack = report.details.get('jackknife')
        freq = report.details.get('jackknife_freq')
        if jack is not None and freq is not None:
            keep(kplot.jackknife(
                freq, jack, period, report.details.get('jackknife_data', data),
                unit='point' if unit == 'point' else 'sequence',
                title=f'{target}: leave one visit out'), 'jackknife')
        for split, coh in report.details.get('coherence', {}).items():
            if len(coh['chunks']) < 2:  # too few visits to compare
                continue
            keep(kplot.coherence(
                coh, title=f'{target}: is the signal coherent? (by {split})'),
                 f'coherence_{split}')
        # 4. the activity indicators
        inds = report.details.get('indicators') or []
        if inds:
            lines = []
            for entry in inds:
                row = (f'{entry["name"]:>12s}  best {entry["best_period"]:9.2f}'
                       f' d (FIP {entry["best_fip"]:.0e})')
                for mult, val in sorted(entry.get('at', {}).items()):
                    fipv = val.get('fip', np.nan) if isinstance(val, dict) \
                        else val
                    row += f'  {mult}P {fipv:.0e}'
                lines.append((row, 'mono'))
            keep(_text_page(
                'Activity indicators', [(f'The outlier-aware FIP of each '
                                         f'indicator at the period and its '
                                         f'multiples (P = {period:.4f} d); a '
                                         f'low FIP there speaks against a '
                                         f'planet.', 'body')] + lines),
                 'indicators')
        # 4b. the photometry of TESS
        phot = report.details.get('tess')
        if phot is not None and phot.get('sectors'):
            from koloa import tess as ktess
            keep(ktess.figure(phot, [period],
                              title=f'{target}: TESS, the periodogram of '
                                    f'each sector'), 'tess')
            chk = phot['checks'][0] if phot['checks'] else {}
            keep(_text_page('TESS photometry', [
                (f'TIC {phot["tic"]}: {phot["summary"]}.', 'body'),
                (f'At P = {period:.4f} d: {chk.get("summary", "")}.',
                 {'flag': 'outlier', 'pass': 'koloa'}.get(
                     chk.get('status'), 'body')),
                ('Spots rotating with the star put its rotation period and '
                 'its harmonics in the brightness; a planet does not, '
                 'unless it transits. A sector lasts 27 days, so periods '
                 'longer than half of it are not tested, and the orbit of '
                 'TESS (13.7 d) and its half leave systematics.', 'body'),
                ('Sector, source, scatter of the 30-minute bins, strongest '
                 'periods (GLS power, semi-amplitude)', 'head')]
                + [(line, 'mono') for line in ktess.table(phot)]),
                 'tess_sectors')
        # 5. the archive, in full
        if known is not None:
            blocks = [(f'Host: {known.get("host")}; star: '
                       f'{known.get("star")}; asked on '
                       f'{known.get("fetched")}.', 'body')]
            for pl in known.get('planets', []):
                blocks.append((f'{pl["name"]}: P = {pl["P"]} d', 'head'))
                for sol in reversed(pl.get('solutions') or []):
                    blocks.append((f'{sol["reference"]:<28s} P = '
                                   f'{sol["P"]}  K = {sol["K"]} +- '
                                   f'{sol["K_err"]}', 'mono'))
            keep(_text_page('What the NASA Exoplanet Archive knows', blocks),
                 'archive')
        info = book.infodict()
        info['Title'] = f'Duck test of {target} at {period:.5f} d'
        info['Creator'] = 'koloa'
    kplot.plt.close('all')
    if outdir:
        with open(os.path.join(outdir, f'{safe}_duck.txt'), 'w') as handle:
            handle.write(report.text() + '\n')
        summary = dict(target=target, period=period, verdict=report.verdict,
                       checks=report.checks, orbit=orbit, archive=known,
                       archive_note=archive_note, report=os.path.basename(path))
        if report.details.get('tess') is not None:
            from koloa.tess import light
            summary['tess'] = light(report.details['tess'])
        summary['family_fip'] = report.details.get('family_fip')
        summary['alias_odds'] = report.details.get('alias_odds')
        if report.details.get('alias_solutions'):
            from koloa.aliases import light as alight
            summary['aliases'] = alight(report.details['alias_solutions'])
        summary['alias_plan'] = report.details.get('alias_plan')
        with open(os.path.join(outdir, f'{safe}_duck.json'), 'w') as handle:
            json.dump(summary, handle, indent=1, default=_json)
        log(f'duck test: the report, its figures, text and summary in '
            f'{outdir}')
    else:
        log(f'duck test report: {path}')
    return path


# =============================================================================
# End of code
# =============================================================================

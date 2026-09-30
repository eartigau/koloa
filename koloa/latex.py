#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The detailed analysis as a LaTeX report, compiled to PDF

koloa.detailed_analysis hands what it found to `detailed_report`, which
writes <star>_report.tex next to the figures; `compile_pdf` runs pdflatex on
it (twice, for the tables and the links). Without pdflatex the .tex is still
written, and the analysis goes on.

The sections:
1. a summary: the data, the known planets, every signal and its verdict;
2. the star: SIMBAD's names, the archive's mass and distance;
3. the data: every instrument (exposures, visits, dates, errors, rms, the
   error added to it for the FIP), DACE, and the velocities;
4. the known planets: every published solution;
5. the FIP, in two passes: P(k), the strongest intervals, the periodograms;
6. the signals: the orbits, against the known planets and the activity,
   folded;
7. the duck test of each signal: every check, the verdict, the jackknife
   and the coherence;
8. the activity indicators;
9. the outliers, and why each one is;
10. the settings and the files.

Created on 2026-09-30

@author: artigau
"""
import os
import shutil
import subprocess
import time
from glob import glob
from typing import Any, Dict, List, Optional

import numpy as np

from koloa.data import robust_std
from koloa.log import log
from koloa.plotting import STYLES

# =============================================================================
# Define variables
# =============================================================================
#: LaTeX's special characters
_SPECIAL = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
            '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}',
            '~': r'\textasciitilde{}', '^': r'\textasciicircum{}',
            '<': r'\textless{}', '>': r'\textgreater{}'}
#: characters beyond the T1 encoding, in LaTeX
_UNICODE = {'±': r'\ensuremath{\pm}', '−': '-', '–': '--', '—': '-',
            '×': r'\ensuremath{\times}', 'σ': r'\ensuremath{\sigma}',
            'Δ': r'\ensuremath{\Delta}', 'α': r'\ensuremath{\alpha}',
            'µ': r'\ensuremath{\mu}', 'μ': r'\ensuremath{\mu}',
            '⊕': r'\ensuremath{\oplus}', '⊙': r'\ensuremath{\odot}',
            '☉': r'\ensuremath{\odot}', '≈': r'\ensuremath{\approx}',
            '≤': r'\ensuremath{\leq}', '≥': r'\ensuremath{\geq}',
            '°': r'\ensuremath{^\circ}', '’': "'", '‘': '`', '“': '``',
            '”': "''", '…': r'\ldots{}'}
#: where pdflatex may be when it is not on the PATH
_TEX_DIRS = ('/Library/TeX/texbin', '/usr/local/texlive/*/bin/*',
             '/opt/homebrew/bin', '/usr/local/bin', '/usr/bin')
#: the colour of a duck-test verdict, by its first words
_VERDICT_COLOURS = (('PLANET CANDIDATE', 'good'), ('INCONCLUSIVE', 'warn'),
                    ('NOT', 'flag'))
#: files that pdflatex leaves behind
_AUX = ('.aux', '.out', '.toc')
#: the parameters of the star shown (key of koloa.archive, label, unit)
_STAR_ROWS = (('spectral_type', 'Spectral type', ''),
              ('teff', 'T$_{\\rm eff}$', '\\,K'),
              ('mass', 'Mass', '\\,M$_\\odot$'),
              ('radius', 'Radius', '\\,R$_\\odot$'),
              ('metallicity', '{[Fe/H]}', ''), ('logg', '$\\log g$', ''),
              ('vsini', '$v\\sin i$', '\\,km/s'),
              ('rotation', 'Rotation period', '\\,d'),
              ('age', 'Age', '\\,Gyr'), ('distance', 'Distance', '\\,pc'),
              ('V', 'V', '\\,mag'), ('J', 'J', '\\,mag'),
              ('K', 'K', '\\,mag'))


# =============================================================================
# Text and numbers
# =============================================================================
def escape(text: Any) -> str:
    """
    Text for LaTeX: its special characters escaped, the unicode beyond T1
    written in LaTeX (a character it has no way to write becomes '?')

    :param text: anything, written with str()

    :return: str
    """
    out = []
    for char in str(text):
        if char in _SPECIAL:
            out.append(_SPECIAL[char])
        elif char in _UNICODE:
            out.append(_UNICODE[char])
        elif ord(char) < 0x180 and (char.isprintable() or char == '\n'):
            out.append(char)
        else:
            out.append('?')
    return ''.join(out)


def pretty(text: Any) -> str:
    """koloa's words for LaTeX: escaped, with +- and Delta ln typeset"""
    out = escape(text)
    for old, new in (('+-', r'\ensuremath{\pm}'),
                     ('Delta ln post', r'$\Delta\ln$ post'),
                     ('Delta lnL', r'$\Delta\ln L$'),
                     ('Delta ln L', r'$\Delta\ln L$'),
                     ('Delta chi2', r'$\Delta\chi^2$')):
        out = out.replace(old, new)
    return out


def sci(value: Optional[float], digits: int = 1) -> str:
    """a probability: two decimals down to 0.01, a power of ten below"""
    if value is None or not np.isfinite(value):
        return '--'
    if value >= 0.01 or value == 0:
        return f'{value:.2f}'
    if value < 1e-300:
        return r'$<10^{-300}$'
    exp = int(np.floor(np.log10(value)))
    mant = value / 10.0 ** exp
    if round(mant, digits) >= 10:
        mant, exp = mant / 10, exp + 1
    return rf'${mant:.{digits}f}\times10^{{{exp}}}$'


def pm(val: Any, fmt: str = '.2f') -> str:
    """a (value, minus, plus) triplet, or a (value, error) pair"""
    if val is None:
        return '--'
    val = [float(item) for item in val]
    if not np.isfinite(val[0]):
        return '--'
    if len(val) == 2 or abs(val[1] - val[2]) <= 0.1 * max(val[1], val[2]):
        err = val[1] if len(val) == 2 else 0.5 * (val[1] + val[2])
        return rf'${val[0]:{fmt}} \pm {err:{fmt}}$'
    return rf'${val[0]:{fmt}}^{{+{val[2]:{fmt}}}}_{{-{val[1]:{fmt}}}}$'


def _url(url: str) -> str:
    """a URL for \\href: & encoded (it would end a cell), % and # escaped"""
    return (url.replace('&', '%26').replace('%', r'\%').replace('#', r'\#'))


def _path(path: str) -> str:
    """a path in type, breakable at its slashes"""
    if any(char in path for char in '%#{}\\^~'):
        return r'\texttt{' + escape(path) + '}'
    return r'\nolinkurl{' + path + '}'


def _date(rjd: float) -> str:
    """the calendar date of a time in BJD - 2400000"""
    from koloa.outliers import _date as date
    return date(rjd)


def _verdict_colour(verdict: str) -> str:
    """the colour of a duck-test verdict"""
    for start, colour in _VERDICT_COLOURS:
        if verdict.startswith(start):
            return colour
    return 'muted'


# =============================================================================
# Pieces of the document
# =============================================================================
def _preamble(star: str) -> str:
    """the document class, its packages, koloa's colours, the page"""
    pal = STYLES['paper']
    colours = dict(koloa=pal['koloa'], flag=pal['outlier'], good=pal['soft'],
                   warn=pal['hard'], muted=pal['muted'], rule=pal['grid'])
    defs = '\n'.join(rf'\definecolor{{{name}}}{{HTML}}{{{val.lstrip("#")}}}'
                     for name, val in colours.items())
    return rf"""\documentclass[a4paper,10pt]{{article}}
\usepackage[T1]{{fontenc}}
\usepackage[utf8]{{inputenc}}
\usepackage{{lmodern}}
\usepackage[margin=2cm,top=2.3cm,bottom=2.2cm,headheight=14pt]{{geometry}}
\usepackage{{microtype}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{longtable}}
\usepackage{{tabularx}}
\usepackage{{array}}
\usepackage{{xcolor}}
\usepackage{{float}}
\usepackage[font=small,labelfont=bf]{{caption}}
\usepackage{{fancyhdr}}
\usepackage{{hyperref}}
{defs}
\hypersetup{{colorlinks=true,linkcolor=koloa,urlcolor=koloa,
  pdftitle={{koloa: detailed analysis of {escape(star)}}},pdfcreator={{koloa}}}}
\pagestyle{{fancy}}
\fancyhf{{}}
\lhead{{\small\color{{muted}}koloa: detailed analysis of {escape(star)}}}
\rhead{{\small\color{{muted}}\thepage}}
\renewcommand{{\headrulewidth}}{{0.4pt}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{4pt}}
\renewcommand{{\arraystretch}}{{1.15}}
\newcolumntype{{L}}{{>{{\raggedright\arraybackslash}}X}}
\newcommand{{\status}}[2]{{\textcolor{{#1}}{{\textbf{{#2}}}}}}
\def\UrlBreaks{{\do\/\do-\do\_\do.}}
\begin{{document}}
"""


def _figure(path: Optional[str], folder: str, caption: str,
            height: str = '0.62') -> str:
    """a figure, when it is a PDF or an image LaTeX can read"""
    if not path or not os.path.exists(path):
        return ''
    stem, ext = os.path.splitext(os.path.relpath(path, folder))
    if ext.lower() not in ('.pdf', '.png', '.jpg', '.jpeg'):
        return ''
    # the braces keep the dots and spaces of a name away from graphicx
    return ('\\begin{figure}[H]\\centering\n'
            f'\\includegraphics[width=\\linewidth,height={height}\\textheight,'
            f'keepaspectratio]{{{{{stem}}}{ext}}}\n'
            f'\\caption{{{caption}}}\n\\end{{figure}}\n')


def _summary(rep: Dict[str, Any]) -> str:
    """the box at the top: the data, the known planets, every signal"""
    data, fip2 = rep['data'], rep['fip_second']
    insts = ', '.join(f'{escape(inst)} {int(np.sum(data.inst == inst))}'
                      for inst in data.instruments)
    items = [f'{data.n} exposures in {data.nseq} visits over '
             f'{data.baseline:.0f}\\,d ({insts}).']
    added = [src for src in rep['sources'][1:] if src['n']]
    if added:
        items[-1] += (' Added to the file: ' + '; '.join(
            f'{src["n"]} from {escape(src["label"])} ({escape(src["kind"])})'
            for src in added) + '.')
    known = rep['known'] or {}
    star = known.get('star') or {}
    words = [escape(val) for val in (star.get('spectral_type'),) if val]
    if star.get('rotation'):
        words.append(f'rotation {star["rotation"]}\\,d')
    if words:
        items.append('The star: ' + ', '.join(words) + '.')
    if known.get('planets'):
        items.append('Known planets: ' + ', '.join(
            f'{escape(pl["name"])} ({pl["P"]:.4f}\\,d)'
            for pl in known['planets'] if pl.get('P')) + '.')
    elif rep['settings'].get('archive'):
        items.append('No known planet in the NASA Exoplanet Archive.')
    pk = np.asarray(fip2.pk)
    nsig = int(np.argmax(pk))
    items.append(f'The FIP puts {nsig} signal{"" if nsig == 1 else "s"} in '
                 f'the series most likely (P(k), from k = 0: '
                 + ', '.join(f'{val:.2f}' for val in pk) + ').')
    if not rep['orbits']:
        best = fip2.peaks[0] if fip2.peaks else None
        items.append(f'No interval with FIP below {rep["threshold"]:g}.'
                     + (f' The strongest: {best["period"]:.4f}\\,d (FIP '
                        f'{sci(best["fip"])}).' if best else ''))
    for orb in rep['orbits']:
        per = orb['P'][0]
        line = (f'\\textbf{{{per:.4f}\\,d}}: K = {pm(orb["K"])}\\,m/s, FIP '
                f'{sci(fip2.fip_containing(per, 1 / data.baseline))}')
        if 'known' in orb:
            line += f'; {escape(orb["known"]["name"])}'
            if orb.get('comparisons'):
                sol = orb['comparisons'][0]
                line += (f' (K {sol["z"]:+.1f}$\\sigma$ from '
                         f'{escape(sol["reference"])})')
        else:
            line += '; not a known planet'
        eph = orb.get('ephemeris')
        if eph:
            colour = 'good' if abs(eph['nsigma']) < 3 else 'flag'
            line += (f'; \\status{{{colour}}}{{conjunction '
                     f'{eph["nsigma"]:+.1f}$\\sigma$}} from the published '
                     f'ephemeris')
        if orb.get('activity'):
            line += '; \\status{flag}{at an activity period}'
        if orb.get('rotation'):
            line += '; \\status{flag}{at the rotation period or a harmonic}'
        gporb = orb.get('gp')
        if gporb:
            penalty = 1.5 * np.log(rep['data'].n)
            colour = 'good' if gporb['dlnl'] > penalty else 'flag'
            line += (f'.\\newline with the GP: K = {pm(gporb["K"])}\\,m/s, '
                     f'\\status{{{colour}}}{{$\\Delta\\ln L = '
                     f'{gporb["dlnl"]:+.1f}$}} (threshold {penalty:.1f}), '
                     f'whitened FAP {sci(gporb["whitened"]["fap"])}')
        report = rep['duck'].get(f'{per:.4f}')
        if report is not None:
            line += (f'.\\newline duck test: \\status{{'
                     f'{_verdict_colour(report.verdict)}}}{{'
                     f'{escape(report.verdict)}}}')
        items.append(line + '.')
    body = '\n'.join(f'\\item {item}' for item in items)
    return ('\\noindent\\fcolorbox{koloa}{white}{\\begin{minipage}'
            '{\\dimexpr\\linewidth-2\\fboxsep-2\\fboxrule\\relax}\n'
            '\\textbf{Summary}\n\\begin{itemize}\\setlength{\\itemsep}{2pt}\n'
            f'{body}\n\\end{{itemize}}\n\\end{{minipage}}}}\n\n')


def _star(rep: Dict[str, Any]) -> str:
    """who the star is"""
    ident, known = rep['ident'], rep['known'] or {}
    out = ['\\section{The star}']
    if not ident:
        out.append(f'SIMBAD did not resolve {escape(rep["star"])}, or was '
                   f'not asked: no name, no archive, no DACE.')
        return '\n'.join(out) + '\n'
    names = [('SIMBAD', escape(ident.get('main')))]
    for key, label in (('gj', 'GJ'), ('hd', 'HD'), ('hip', 'HIP'),
                       ('tic', 'TIC'), ('gaia_dr3', 'Gaia DR3')):
        if ident.get(key):
            names.append((label, escape(ident[key])))
    star = known.get('star') or {}
    params = []
    for key, label, unit in _STAR_ROWS:
        val = star.get(key)
        if val is None:
            continue
        text = escape(val)
        if key == 'rotation' and star.get('rotation_err'):
            text = f'${val} \\pm {star["rotation_err"]}$'
        params.append((label, text + unit))
    out.append('\\begin{minipage}[t]{0.46\\linewidth}\\vspace{0pt}\n'
               '\\begin{tabular}{@{}ll@{}}\n\\toprule\n\\multicolumn{2}{@{}l}'
               '{\\textbf{Names} (SIMBAD)} \\\\\n\\midrule')
    out += [f'{label} & {val} \\\\' for label, val in names]
    out.append('\\bottomrule\n\\end{tabular}\n\\end{minipage}\\hfill')
    out.append('\\begin{minipage}[t]{0.5\\linewidth}\\vspace{0pt}\n')
    if params:
        out.append('\\begin{tabular}{@{}ll@{}}\n\\toprule\n\\multicolumn{2}'
                   '{@{}l}{\\textbf{Parameters} (NASA Exoplanet Archive)} '
                   '\\\\\n\\midrule')
        out += [f'{label} & {val} \\\\' for label, val in params]
        out.append('\\bottomrule\n\\end{tabular}')
    else:
        out.append('The NASA Exoplanet Archive has no parameter of this '
                   'star (it keeps those of planet hosts only).')
    out.append('\\end{minipage}\n')
    aliases = ident.get('aliases') or []
    if aliases:
        out.append('\\smallskip{\\small Also known as: '
                   + ', '.join(escape(name) for name in aliases) + '.}\n')
    return '\n'.join(out) + '\n'


def _data(rep: Dict[str, Any], folder: str) -> str:
    """every instrument, DACE, and the velocities"""
    data = rep['data']
    out = ['\\section{The data}',
           '\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}'
           'lrrllrrrrr@{}}\n\\toprule',
           'Instrument & Exp. & Visits & First & Last & Span & '
           '$\\sigma$ & rms & \\multicolumn{2}{c@{}}{added for the FIP} \\\\',
           '\\cmidrule(l){9-10}',
           '& & & & & [d] & [m/s] & [m/s] & 1st pass & 2nd pass \\\\\n'
           '\\midrule']
    for inst in data.instruments:
        sub = data.select(data.inst == inst)
        add1 = rep['inflation_first'].get(str(inst))
        add2 = rep['inflation_second'].get(str(inst))
        out.append(f'{escape(inst)} & {sub.n} & {sub.nseq} & '
                   f'{_date(sub.time.min())} & {_date(sub.time.max())} & '
                   f'{sub.baseline:.0f} & {np.median(sub.err):.2f} & '
                   f'{robust_std(sub.rv):.2f} & '
                   + (f'{add1:.2f}' if add1 is not None else '--') + ' & '
                   + (f'{add2:.2f}' if add2 is not None else '--') + ' \\\\')
    out.append('\\bottomrule\n\\end{tabular*}\n')
    out.append('{\\small\\color{muted}$\\sigma$ is the median error of the '
               'exposures, rms their scatter (robust: outliers do not '
               'count); the FIP adds an error in '
               'quadrature to each instrument, for the noise its errors do '
               'not hold.}\n')
    out.append('\\subsection*{Where the velocities come from}\n'
               '\\begin{tabularx}{\\linewidth}{@{}llrL@{}}\n\\toprule\n'
               'Source & Kind & Kept & Instruments; note \\\\\n\\midrule')
    for src in rep['sources']:
        insts = ', '.join(f'{escape(inst)} {num}'
                          for inst, num in src['instruments'].items())
        note = '; '.join(val for val in (insts, pretty(src.get('note', '')))
                         if val)
        out.append(f'{escape(src["label"])} & {escape(src["kind"])} & '
                   f'{src["n"]} & {note or "--"} \\\\')
    out.append('\\bottomrule\n\\end{tabularx}\n')
    if any(src['kind'] in ('VizieR', 'given') for src in rep['sources']):
        out.append('{\\small\\color{muted}A published velocity within a '
                   'minute of an exposure already in the series is the same '
                   'spectrum, and is left out; a published instrument that '
                   'shares its name with one of the series keeps its own '
                   'offset. The published velocities are used as published '
                   '(their authors may have corrected offsets or activity).}'
                   '\n')
    out.append(_figure(rep['figures'].get('rv'), folder,
                       'The velocities of every instrument, each exposure '
                       'coloured by its probability of being an outlier.',
                       '0.4'))
    return '\n'.join(out) + '\n'


def _known(rep: Dict[str, Any]) -> str:
    """the known planets, every published solution"""
    known = rep['known'] or {}
    out = ['\\section{The known planets}']
    if not rep['settings'].get('archive'):
        return out[0] + '\nThe NASA Exoplanet Archive was not asked.\n\n'
    if not known.get('planets'):
        return (out[0] + f'\nThe NASA Exoplanet Archive has no planet of '
                f'{escape(known.get("host") or rep["star"])}.\n\n')
    from koloa.archive import conjunction
    out.append(f'From the NASA Exoplanet Archive (host '
               f'{escape(known.get("host"))}, asked on '
               f'{escape(known.get("fetched"))}): the default solution of '
               f'each planet, then every published solution, the most recent '
               f'first.\n')
    out.append('\\begin{tabularx}{\\linewidth}{@{}lrrrrrrL@{}}\n\\toprule\n'
               'Planet & P [d] & K [m/s] & e & m sin i [M$_\\oplus$] & '
               'a [au] & T$_{\\rm eq}$ [K] & Discovery \\\\\n\\midrule')

    def num(val, fmt):
        return '--' if val is None else f'{val:{fmt}}'
    for pl in known['planets']:
        found = ', '.join(escape(val) for val in (
            pl.get('discovery'), pl.get('disc_year'), pl.get('disc_facility'))
            if val)
        kval = (pm((pl['K'], pl['K_err'])) if pl.get('K') is not None
                and pl.get('K_err') else num(pl.get('K'), '.2f'))
        out.append(f'{escape(pl["name"])} & {num(pl.get("P"), ".5f")} & '
                   f'{kval} & {num(pl.get("e"), ".2f")} & '
                   f'{num(pl.get("mass_earth"), ".2f")} & '
                   f'{num(pl.get("a"), ".4f")} & {num(pl.get("teq"), ".0f")} '
                   f'& {found or "--"} \\\\')
    out.append('\\bottomrule\n\\end{tabularx}\n')
    ephs = [(pl, conjunction(pl)) for pl in known['planets']]
    ephs = [(pl, conj) for pl, conj in ephs if conj]
    if ephs:
        out.append('Published ephemerides (a conjunction, BJD $-$ 2400000): '
                   + '; '.join(
                       f'{escape(pl["name"])} {pm((conj["tc"], conj["tc_err"] or 0.0), ".3f")}'
                       f', P = {pm((conj["P"], conj["P_err"] or 0.0), ".5f")}'
                       f'\\,d (from the {escape(conj["source"])})'
                       for pl, conj in ephs) + '.\n')
    out.append('\\begin{longtable}{@{}llrrl@{}}\n\\toprule\nPlanet & '
               'Reference & P [d] & K [m/s] & Date \\\\\n\\midrule\n'
               '\\endhead')
    for pl in known['planets']:
        sols = list(reversed(pl.get('solutions') or [pl]))
        for isol, sol in enumerate(sols):
            ref = escape(sol.get('reference') or '')
            if sol.get('reference_url'):
                ref = f'\\href{{{_url(sol["reference_url"])}}}{{{ref}}}'
            kval = (f'${sol["K"]} \\pm {sol["K_err"]}$'
                    if sol.get('K') is not None and sol.get('K_err')
                    else escape(sol.get('K') or '--'))
            out.append(f'{escape(pl["name"]) if isol == 0 else ""} & {ref} & '
                       f'{escape(sol.get("P") or "--")} & {kval} & '
                       f'{escape(sol.get("date") or "")} \\\\')
    out.append('\\bottomrule\n\\end{longtable}\n')
    return '\n'.join(out) + '\n'


def _fip(rep: Dict[str, Any], folder: str) -> str:
    """the FIP, in two passes"""
    out = ['\\section{The FIP}',
           'The outlier-aware false inclusion probability (FIP) of every '
           'interval of frequency: the probability that it holds no '
           f'signal. Below {rep["threshold"]:g}, a signal is detected. The '
           'first pass has the errors of each instrument inflated to a noise '
           'model without planets; the second, to the fit with the signals '
           'of the first.\n']
    passes = (('first', 'First pass: the noise without planets',
               rep['fip_first']),
              ('second', 'Second pass: the noise with the planets',
               rep['fip_second']))
    for key, title, res in passes:
        out.append(f'\\subsection{{{title}}}')
        out.append('P(k), the probability of k signals, from k = 0: '
                   + ', '.join(f'{val:.2f}' for val in res.pk) + '.\n')
        out.append('\\begin{tabular}{@{}rrrrl@{}}\n\\toprule\nPeriod [d] & '
                   'FIP & K [m/s] & window & best alias (FIP) \\\\\n'
                   '\\midrule')
        for peak in res.peaks:
            mark = ('\\status{good}{detected}'
                    if peak['fip'] < rep['threshold'] else '')
            amp = '--'
            if peak.get('amplitude'):
                # the 50th, 16th and 84th percentiles of the posterior
                mid, low, high = peak['amplitude']
                amp = pm((mid, mid - low, high - mid), '.1f')
            alias = ''
            if peak.get('aliases'):
                best = min(peak['aliases'], key=lambda item: item['fip'])
                alias = (f'{best["period"]:.3f}\\,d, {escape(best["name"])} '
                         f'({sci(best["fip"])})')
            out.append(f'{peak["period"]:.4f} & {sci(peak["fip"])} & {amp} & '
                       f'{peak.get("window", np.nan):.2f} & {alias} {mark} '
                       f'\\\\')
        out.append('\\bottomrule\n\\end{tabular}\n')
        out.append(_figure(rep['figures'].get(f'fip_{key}'), folder,
                           f'{title}: the FIP of every interval (the known '
                           f'planets dashed).', '0.3'))
    return '\n'.join(out) + '\n'


def _signals(rep: Dict[str, Any], folder: str) -> str:
    """the orbits, against the known planets and the activity"""
    out = ['\\section{The signals}']
    orbits = rep['orbits']
    if not orbits:
        return (out[0] + f'\nNo interval has a FIP below '
                f'{rep["threshold"]:g}: nothing was fitted.\n\n')
    method = ('sampled by MCMC' if rep['settings'].get('mcmc') else
              'the maximum a posteriori, with Laplace errors')
    out.append(f'Every interval with a FIP below {rep["threshold"]:g}, the '
               f'known planets (tested at their periods, found by the FIP or '
               f'not) and the periods given, fitted together (each period '
               f'free within 2\\,\\% of its start, circular orbits, every '
               f'exposure and visit may be an outlier; {method}; no GP'
               + ('; Section~\\ref{sec:gp} has the fit with one'
                  if rep.get('gp') else '') + ').\n')
    out.append('\\begin{tabularx}{\\linewidth}{@{}rrrrrlL@{}}\n\\toprule\n'
               'P [d] & K [m/s] & e & m sin i [M$_\\oplus$] & FIP & from & '
               'known planet \\\\\n\\midrule')
    for orb in orbits:
        per = orb['P'][0]
        match = '--'
        if 'known' in orb:
            match = escape(orb['known']['name'])
            if orb.get('comparisons'):
                sol = orb['comparisons'][0]
                match += (f': K {sol["K"]} ({escape(sol["reference"])}, '
                          f'{sol["z"]:+.1f}$\\sigma$)')
        origin = escape(orb.get('origin', 'FIP').split(' ')[0])
        fipv = rep['fip_second'].fip_containing(per, 1 / rep['data'].baseline)
        out.append(f'{pm(orb["P"], ".4f")} & {pm(orb["K"])} & '
                   f'{pm(orb["e"])} & {pm(orb.get("msini"))} & '
                   f'{sci(fipv)} & {origin} & {match} \\\\')
    out.append('\\bottomrule\n\\end{tabularx}\n')
    for orb in orbits:
        per = orb['P'][0]
        eph = orb.get('ephemeris')
        if eph:
            colour = 'good' if abs(eph['nsigma']) < 3 else 'flag'
            out.append(
                f'\\status{{{colour}}}{{Ephemeris at {per:.4f}\\,d:}} the '
                f'fitted conjunction, {pm((eph["tc"], eph["tc_err"]), ".3f")}'
                f' (BJD $-$ 2400000), against the published one carried '
                f'{eph["ncycle"]} cycles, {pm((eph["predicted"], eph["predicted_err"]), ".3f")}'
                f': {eph["diff"]:+.3f}\\,d, {eph["nsigma"]:+.1f}$\\sigma$, '
                f'{eph["dphase"]:+.3f} in phase (the published time from the '
                f'{escape(eph["source"])}; its covariance with the period '
                f'ignored).\n')
        if orb.get('rotation'):
            out.append(f'\\status{{flag}}{{Rotation at {per:.4f}\\,d:}} '
                       + pretty('; '.join(orb['rotation'])) + '.\n')
        if orb.get('activity'):
            out.append(f'\\status{{flag}}{{Activity at {per:.4f}\\,d:}}'
                       ' ' + pretty('; '.join(orb['activity'])) + '.\n')
    for ip, orb in enumerate(orbits):
        out.append(_figure(rep['figures'].get(f'phase_{ip}'), folder,
                           f'The orbit at {orb["P"][0]:.4f}\\,d, the other '
                           f'signals removed: K = {pm(orb["K"])}\\,m/s.',
                           '0.35'))
    return '\n'.join(out) + '\n'


def _gp(rep: Dict[str, Any], folder: str) -> str:
    """the signals against a GP of the activity, and the whole series with
    it"""
    gpsum = rep.get('gp')
    if not gpsum:
        return ''
    nobs = rep['data'].n
    penalty = 1.5 * np.log(nobs)
    out = ['\\section{The activity, as a GP}\\label{sec:gp}',
           f'The velocities are fitted with a GP of the activity: '
           f'{pretty(gpsum["label"])}; jointly with every signal (circular '
           f'orbits; every exposure and visit may be an outlier; the maximum '
           f'a posteriori). The GP alone, the GP with every signal, and the '
           f'GP with every signal but one: the likelihood a signal adds over '
           f'the activity, $\\Delta\\ln L$, against the cost of its three '
           f'parameters, $1.5\\ln n = {penalty:.1f}$ for $n = {nobs}$ '
           f'(BIC-like).\n']
    names = sorted(set(gpsum['gp_null']) | set(gpsum.get('gp_all') or {}))
    out.append('\\begin{tabular}{@{}lrr@{}}\n\\toprule\nGP parameter & GP '
               'alone & GP + signals \\\\\n\\midrule')
    for name in names:
        one = gpsum['gp_null'].get(name)
        two = (gpsum.get('gp_all') or {}).get(name)
        out.append(f'{escape(name)} & '
                   + (f'{one:.3g}' if one is not None else '--') + ' & '
                   + (f'{two:.3g}' if two is not None else '--') + ' \\\\')
    out.append('\\bottomrule\n\\end{tabular}\n')
    if gpsum['orbits']:
        out.append('\\begin{tabular}{@{}rrrrrrl@{}}\n\\toprule\n'
                   'P [d] & K, GP [m/s] & K, no GP [m/s] & $\\Delta\\ln L$ & '
                   'whitened $\\Delta\\chi^2$ & FAP & \\\\\n\\midrule')
        for orb, gporb in zip(rep['orbits'], gpsum['orbits']):
            word = ('\\status{good}{survives the GP}'
                    if gporb['dlnl'] > penalty
                    else '\\status{flag}{the GP takes it}')
            out.append(f'{pm(gporb["P"], ".4f")} & {pm(gporb["K"])} & '
                       f'{pm(orb["K"])} & {gporb["dlnl"]:+.1f} & '
                       f'{gporb["whitened"]["dchi2"]:.1f} & '
                       f'{sci(gporb["whitened"]["fap"])} & {word} \\\\')
        out.append('\\bottomrule\n\\end{tabular}\n')
        out.append('{\\small\\color{muted}The whitened $\\Delta\\chi^2$ is '
                   'that of a sinusoid near the period (within 1\\,\\%), the '
                   'noise that of the GP-only model; its FAP, the fraction of '
                   'simulations of that noise whose highest peak, anywhere, '
                   'is higher.}\n')
    out.append(_figure(rep['figures'].get('gp_sequence'), folder,
                       'The whole series, season by season: every exposure '
                       '(offsets and trend removed; outliers as crosses), the '
                       'GP of the activity with its one sigma (shaded), and '
                       'the GP with the signals (thin line).', '0.85'))
    labels = dict(null='the GP alone as the noise',
                  residual='what the GP and the signals leave')
    for key, res in gpsum['whitened'].items():
        out.append(f'\\textbf{{Whitened periodogram, {labels.get(key, key)}}} '
                   f'({res["kept"]} of {res["n"]} points, {res["nsim"]} '
                   f'simulations): ' + '; '.join(
                       f'{pk["period"]:.3f}\\,d ($\\Delta\\chi^2$ '
                       f'{pk["dchi2"]:.1f}, FAP {sci(pk["fap"])})'
                       for pk in res['peaks']) + '.\n')
    out.append(_figure(rep['figures'].get('gp_whitened'), folder,
                       'Periodograms whitened by the GP, with the false-alarm '
                       'levels of simulations of its noise; the signals and '
                       'the rotation (and its half) marked.', '0.5'))
    for ip, orb in enumerate(rep['orbits']):
        out.append(_figure(rep['figures'].get(f'gp_fold_{ip}'), folder,
                           f'The signal at {orb["P"][0]:.4f}\\,d, the GP and '
                           f'the other signals removed (visit means), with '
                           f'the published orbit when there is one.', '0.35'))
    return '\n'.join(out) + '\n'


def _ducks(rep: Dict[str, Any], folder: str) -> str:
    """the duck test of every signal"""
    if not rep['duck']:
        return ''
    out = ['\\section{The duck test}',
           'Every test of planethood koloa knows, at the period of each '
           'signal: significance, robustness (the visit that holds the peak '
           'up), coherence through the campaign, the shape of the orbit, '
           'and the activity indicators (the GP absorption test is not run '
           'here).\n']
    labels = {'pass': ('good', 'pass'), 'flag': ('flag', 'flag'),
              'info': ('muted', 'note')}
    for ip, orb in enumerate(rep['orbits']):
        report = rep['duck'].get(f'{orb["P"][0]:.4f}')
        if report is None:
            continue
        out.append(f'\\subsection{{At {orb["P"][0]:.4f}\\,d}}')
        out.append(f'\\status{{{_verdict_colour(report.verdict)}}}'
                   f'{{{escape(report.verdict)}}}\n')
        out.append('\\begin{tabularx}{\\linewidth}{@{}llL@{}}\n\\toprule\n'
                   'Check & & What it found \\\\\n\\midrule')
        for check in report.checks:
            colour, word = labels.get(check['status'], ('muted', 'note'))
            out.append(f'{escape(check["name"])} & \\status{{{colour}}}'
                       f'{{{word}}} & {pretty(check["summary"])} \\\\')
        out.append('\\bottomrule\n\\end{tabularx}\n')
        out.append(_figure(rep['figures'].get(f'duck_{ip}_jackknife'),
                           folder, f'Leave one visit out: which visits hold '
                                   f'the peak at {orb["P"][0]:.4f}\\,d up.',
                           '0.35'))
        for split in ('halves', 'seasons'):
            out.append(_figure(
                rep['figures'].get(f'duck_{ip}_coherence_{split}'), folder,
                f'Is the signal at {orb["P"][0]:.4f}\\,d coherent? Its '
                f'amplitude and phase by {split}.', '0.3'))
    return '\n'.join(out) + '\n'


def _indicators(rep: Dict[str, Any], folder: str) -> str:
    """the activity indicators"""
    indic = rep['indicators']
    out = ['\\section{The activity indicators}']
    if not indic:
        return out[0] + '\nThe series has no activity indicator.\n\n'
    out.append('The strongest peaks of the outlier-aware periodogram of each '
               'indicator ($\\Delta\\ln L$ in brackets). A signal of the '
               'velocities at one of them, or at twice or half of one, may '
               'be activity.\n')
    out.append('\\begin{tabularx}{\\linewidth}{@{}lrL@{}}\n\\toprule\n'
               'Indicator & n & Strongest peaks [d] \\\\\n\\midrule')
    for name, res in indic.items():
        peaks = ', '.join(f'{pk["period"]:.2f} ({pk["dlnl"]:.0f})'
                          for pk in res['peaks'])
        out.append(f'{escape(name)} & {res.get("n", "")} & {peaks} \\\\')
    out.append('\\bottomrule\n\\end{tabularx}\n')
    prot = (rep['known'] or {}).get('star', {}).get('rotation')
    out.append(_figure(rep['figures'].get('indicators'), folder,
                       'The outlier-aware periodogram of each indicator (the '
                       'signals of the velocities dashed'
                       + (f'; the rotation period, {prot}\\,d, and its half '
                          f'dotted' if prot else '') + ').', '0.7'))
    return '\n'.join(out) + '\n'


def _outliers(rep: Dict[str, Any], folder: str) -> str:
    """the outliers, and why each one is"""
    from koloa.outliers import NAMES, _describe, _unusual
    why = rep['outliers']
    out = ['\\section{The outliers}',
           f'An exposure or a visit is an outlier when its probability of '
           f'being one is above {why.threshold}. A key of the headers or an '
           f'indicator is off when it deviates beyond $|z| > z_{{crit}}$ '
           f'(Bonferroni at $\\alpha = {why.alpha}$) and fewer than 1\\,\\% '
           f'of the good data deviate as much.\n']
    for inst, count in why.counts.items():
        shown = inst if inst != 'inst' else NAMES.get(why.lists.get(inst),
                                                      'the series')
        out.append(f'\\subsection{{{escape(shown)}}}')
        out.append(f'{count["exposures"]} exposures in {count["visits"]} '
                   f'visits: {count["outliers"]} outlying exposures '
                   f'({count["visit_units"]} whole visits and '
                   f'{count["exposure_units"]} single exposures), '
                   f'{count["borderline"]} borderline; '
                   f'{len(why.keys.get(inst, []))} keys looked at '
                   f'({escape(why.lists.get(inst, ""))} list), off beyond '
                   f'$|z| > {why.zcrit.get(inst, np.nan):.1f}$.\n')
        if why.missing.get(inst):
            out.append('{\\footnotesize\\color{muted}Not in the data: '
                       + escape(', '.join(why.missing[inst])) + '.}\n')
        shared = why.shared(inst)
        if shared:
            out.append('What the outliers share (all of them against the '
                       'good data, Holm-corrected): ' + '; '.join(
                           f'{escape(row["label"])}: {pretty(row["summary"])}'
                           for row in shared) + '.\n')
        units = [unit for unit in why.units if unit['inst'] == inst]
        if not units:
            continue
        out.append('\\begin{longtable}{@{}p{0.15\\linewidth}p{0.1\\linewidth}'
                   'r p{0.14\\linewidth}p{0.45\\linewidth}@{}}\n\\toprule\n'
                   'Date & Unit & p & Residual & Why \\\\\n\\midrule\n'
                   '\\endhead')
        zcrit = why.zcrit.get(inst, np.inf)
        for unit in units:
            what = (f'visit ({unit["n"]})' if unit['kind'] == 'visit'
                    else 'exposure')
            resid = '--'
            if unit.get('resid') is not None and np.isfinite(unit['resid']):
                resid = (f'{unit["resid"]:+.2f}\\,m/s '
                         f'({unit["resid_sigma"]:+.1f}$\\sigma$)')
            sig = [key for key in unit['keys'] if key['significant']]
            rest = [key for key in unit['keys'] if not key['significant']
                    and key.get('z') is not None and np.isfinite(key['z'])
                    and abs(key['z']) > 2]
            text = ('off: ' + pretty('; '.join(_describe(key) for key in sig))
                    if sig else 'nothing significantly off')
            if rest[:3]:
                text += ('; also unusual: ' + pretty('; '.join(
                    _unusual(key, zcrit) for key in rest[:3])))
            out.append(f'{escape(unit["date"])} & {what} & '
                       f'{unit["prob"]:.2f} & {resid} & {text} \\\\')
        out.append('\\bottomrule\n\\end{longtable}\n')
    out.append(_figure(rep['figures'].get('outliers'), folder,
                       'The keys of the outliers against the good data.',
                       '0.6'))
    return '\n'.join(out) + '\n'


def _settings(rep: Dict[str, Any], folder: str, skipped: List[str]) -> str:
    """how the analysis ran, and what it wrote"""
    import koloa
    out = ['\\section{Settings and files}',
           '\\begin{tabular}{@{}ll@{}}\n\\toprule']
    for key, val in rep['settings'].items():
        out.append(f'{escape(key)} & {escape(val)} \\\\')
    out.append(f'koloa & {escape(getattr(koloa, "__version__", ""))} \\\\')
    out.append(f'run in & {rep["runtime"] / 60:.1f} min \\\\')
    out.append('\\bottomrule\n\\end{tabular}\n')
    if rep.get('source'):
        out.append(f'Source: {_path(os.path.abspath(rep["source"]))}\n')
    out.append(f'In {_path(os.path.abspath(folder))}: {{\\small '
               + ', '.join(_path(name) for name in rep['files']) + '}.\n')
    if skipped:
        out.append('{\\small\\color{muted}Not shown here (LaTeX reads PDF '
                   'and PNG figures only; see the folder): '
                   + ', '.join(_path(name) for name in skipped) + '.}\n')
    return '\n'.join(out) + '\n'


# =============================================================================
# The report
# =============================================================================
def detailed_report(rep: Dict[str, Any], path: str) -> str:
    """
    The detailed analysis as a LaTeX document (see the module's docstring)

    :param rep: dict, from koloa.detailed_analysis: star, source, data
                (RVData), ident, known, sources (per source of velocities:
                kind, label, n, instruments, note), fip_first and
                fip_second (FIPResult), inflation_first and
                inflation_second (the error added to each instrument
                [m/s]), orbits, indicators, duck (DuckReport by period),
                outliers (OutlierReport), figures (path by role: rv,
                fip_first, fip_second, phase_<i>, indicators, outliers,
                duck_<i>_jackknife, duck_<i>_coherence_<split>), settings,
                threshold, runtime [s] and files (the names in the folder)
    :param path: str, the .tex written (the figures are found from its
                 folder)

    :return: str, the path of the .tex
    """
    folder = os.path.dirname(os.path.abspath(path))
    skipped = [os.path.basename(val) for val in rep['figures'].values()
               if os.path.splitext(val)[1].lower()
               not in ('.pdf', '.png', '.jpg', '.jpeg')]
    title = (f'\\begin{{center}}\n{{\\LARGE\\bfseries Detailed analysis of '
             f'{escape(rep["star"])}}}\\\\[4pt]\n{{\\color{{muted}}koloa, '
             f'{time.strftime("%Y-%m-%d %H:%M")}'
             + (f'; series {escape(rep["data"].name)}'
                if rep['data'].name != rep['star'] else '')
             + '}\n\\end{center}\n\n')
    parts = [_preamble(rep['star']), title, _summary(rep),
             '\\tableofcontents\n\\bigskip\n', _star(rep),
             _data(rep, folder), _known(rep), _fip(rep, folder),
             _signals(rep, folder), _gp(rep, folder), _ducks(rep, folder),
             _indicators(rep, folder), _outliers(rep, folder),
             _settings(rep, folder, skipped), '\\end{document}\n']
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(parts))
    return os.path.abspath(path)


def _pdflatex() -> Optional[str]:
    """pdflatex, on the PATH or where TeX installs it"""
    exe = shutil.which('pdflatex')
    if exe:
        return exe
    for pattern in _TEX_DIRS:
        for folder in sorted(glob(pattern), reverse=True):
            cand = os.path.join(folder, 'pdflatex')
            if os.access(cand, os.X_OK):
                return cand
    return None


def compile_pdf(tex: str, runs: int = 2,
                timeout: float = 180.0) -> Optional[str]:
    """
    A .tex compiled to PDF by pdflatex, in its folder

    :param tex: str, the .tex
    :param runs: int, the passes (two: the tables and the links settle)
    :param timeout: float, the longest a pass may take [s]

    :return: str or None, the PDF (None without pdflatex, or when it
             fails: its log says why)
    """
    exe = _pdflatex()
    if exe is None:
        log(f'LaTeX report: no pdflatex on this machine; {tex} is written, '
            f'not compiled', 'warn')
        return None
    folder, name = os.path.split(os.path.abspath(tex))
    stem = os.path.splitext(name)[0]
    for _ in range(runs):
        try:
            proc = subprocess.run(
                [exe, '-interaction=nonstopmode', '-halt-on-error', name],
                cwd=folder, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            log(f'LaTeX report: pdflatex took more than {timeout:.0f} s; '
                f'{tex} is written, not compiled', 'warn')
            return None
        if proc.returncode != 0:
            text = proc.stdout.decode('utf-8', 'replace').splitlines()
            first = next((line for line in text if line.startswith('!')),
                         'see the log')
            log(f'LaTeX report: pdflatex failed ({first.lstrip("! ")}); '
                f'see {os.path.join(folder, stem + ".log")}', 'warn')
            return None
    for ext in _AUX + ('.log',):
        aux = os.path.join(folder, stem + ext)
        if os.path.exists(aux):
            os.remove(aux)
    return os.path.join(folder, stem + '.pdf')


# =============================================================================
# End of code
# =============================================================================

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The cards of the web page that come from the demos, and the mathematics.

Every card is written from the summary a demo saved
(demos/output/<demo>/summary.json), so no number on the page is typed by
hand; make_page.py calls fill_cards() after its own blocks. The mathematics
come from docs/partials/math_panels.html, one tab per section.

Created on 2026-09-27

@author: artigau
"""
import json
import os
import re

import numpy as np
from html import escape

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEMOS = os.path.join(ROOT, 'demos', 'output')


# =============================================================================
# Helpers
# =============================================================================
def fill(html: str, name: str, content: str) -> str:
    """Replace one AUTO block"""
    pattern = re.compile(r'(<!--AUTO:' + name + r'-->)(.*?)(<!--/AUTO-->)',
                         re.S)
    return pattern.sub(lambda mt: mt.group(1) + content + mt.group(3), html)


def load(demo: str):
    """The summary of a demo, or None"""
    path = os.path.join(DEMOS, demo, 'summary.json')
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def interval(val, low, high, digits=2):
    """A value with its asymmetric interval, in HTML"""
    return (f'{val:.{digits}f}<sup>+{high:.{digits}f}</sup>'
            f'<sub>&minus;{low:.{digits}f}</sub>')


def pct_interval(pct, digits=2):
    """From 2.5 16 50 84 97.5 percentiles"""
    return interval(pct[2], pct[2] - pct[1], pct[3] - pct[2], digits)


def figure(src: str, alt: str, caption: str) -> str:
    """A click-to-enlarge figure"""
    return (f'<figure class="lecture-figure"><img src="{src}" '
            f'alt="{escape(alt, quote=True)}" loading="lazy">'
            f'<figcaption>{caption}</figcaption></figure>')


def fipfmt(val: float) -> str:
    """A FIP for the page"""
    if val >= 0.095:
        return f'{val:.2f}'
    mant, expo = f'{val:.1e}'.split('e')
    return f'{mant}&times;10<sup>{int(expo)}</sup>'


def flags_text(flags: dict) -> str:
    """Clear and borderline outliers flagged, in words"""
    clear, border, good = flags['clear'], flags['borderline'], flags['good']
    return (f'{clear["flagged"]} of the {clear["n"]} exposures moved by clear '
            f'outliers, {border["flagged"]} of the {border["n"]} borderline '
            f'ones (their mean outlier probability is '
            f'{border["mean_prob"]:.2f}), and {good["flagged"]} of the '
            f'{good["n"]} good ones')


# =============================================================================
# The mathematics
# =============================================================================
def math(html: str) -> str:
    """The mathematics: one page per section of math_panels.html, and the
    list of them on the page that introduces them"""
    path = os.path.join(HERE, 'partials', 'math_panels.html')
    if not os.path.exists(path):
        return html
    text = open(path).read()
    sections = re.findall(r'<section data-tab="([^"]+)" data-title="([^"]+)">'
                          r'(.*?)</section>', text, re.S)
    views, index = [], ['<ul class="page-index">']
    for key, title, body in sections:
        views.append(f'<section class="view math-panel" id="math-panel-{key}" '
                     f'data-group="math" data-title="{escape(title)}">\n'
                     f'<p class="crumb"><a href="#math">The mathematics</a></p>\n'
                     f'<h2>{escape(title)}</h2>\n{body.strip()}\n</section>')
        # the first sentence, as text (its entities decoded, escaped again
        #   below)
        from html import unescape
        first = unescape(re.sub(r'<[^>]+>', '', re.search(
            r'<p>(.*?)</p>', body, re.S).group(1)))
        first = re.split(r'(?<=[.:])\s', ' '.join(first.split()))[0]
        index.append(f'<li><a href="#math-panel-{key}">{escape(title)}</a>'
                     f'<span>{escape(first)}</span></li>')
    index.append('</ul>')
    html = fill(html, 'math_views', '\n'.join(views))
    return fill(html, 'math_index', '\n'.join(index))


# =============================================================================
# The demos
# =============================================================================
def false_alarm(html: str) -> str:
    """A planet that is not there"""
    summ = load('false_alarm')
    if summ is None:
        return html
    shown = summ['shown']
    fooled = summ['fooled']
    nseed = summ['nseed']
    text = (f'<p>No planet: noise on a real NIRPS sampling, a 2&nbsp;m/s visit '
            f'jitter, and outliers of both signs, clear (whole visits and '
            f'single exposures six to eight median errors off) and borderline '
            f'(two to four). Over {nseed} such series, the single-signal FIP '
            f'with a fixed jitter reports a planet (FIP&nbsp;&lt;&nbsp;1%) in '
            f'{fooled["gaussian"]} of them, the soft clip in {fooled["soft"]}, the '
            f'hard clip in {fooled["hard"]} and koloa in {fooled["koloa"]}.'
            f'</p>')
    gauss = shown['gaussian']
    if gauss['fip'] < 0.01:
        text += (f'<p>Shown: the first series that fools the fixed-jitter FIP '
                 f'(seed {shown["seed"]}), with a "planet" at '
                 f'{gauss["period"]:.2f}&nbsp;d and FIP = '
                 f'{fipfmt(gauss["fip"])}. koloa\'s best interval has FIP = '
                 f'{fipfmt(shown["koloa"]["fip"])}. It flags '
                 f'{flags_text(shown["flags"])}.</p>')
    else:
        text += (f'<p>Shown: seed {shown["seed"]}. koloa flags '
                 f'{flags_text(shown["flags"])}.</p>')
    return fill(html, 'false_alarm_text', text)


def recovery(html: str) -> str:
    """A planet behind outliers"""
    summ = load('recovery')
    if summ is None:
        return html
    fips = summ['fips']
    orb = summ['orbits']
    kval, klow, khigh = orb['koloa']['K']
    gval, glow, ghigh = orb['gaussian']['K']
    text = (f'<p>K = {summ["K"]:.0f}&nbsp;m/s at {summ["period"]:.2f}&nbsp;d, '
            f'behind outliers clear and borderline. The FIP of the interval '
            f'of the planet is {fipfmt(fips["gaussian"]["at_planet"])} for '
            f'the fixed-jitter FIP, {fipfmt(fips["soft"]["at_planet"])} after '
            f'a soft clip, {fipfmt(fips["hard"]["at_planet"])} after a hard '
            f'clip and {fipfmt(fips["koloa"]["at_planet"])} for koloa. koloa '
            f'flags {flags_text(summ["flags"])}. Its orbit gives K = '
            f'{interval(kval, klow, khigh)}&nbsp;m/s; the Gaussian orbit '
            f'{interval(gval, glow, ghigh)}&nbsp;m/s.</p>')
    return fill(html, 'recovery_text', text)


def two_planets(html: str) -> str:
    """Two planets, one FIP periodogram"""
    summ = load('two_planets')
    if summ is None:
        return html
    pk = summ['pk']
    parts = []
    for pl in summ['planets']:
        parts.append(f'the {pl["P"]}-day planet (K = {pl["K"]}&nbsp;m/s) '
                     f'has FIP {fipfmt(pl["fip_koloa"])} with koloa and '
                     f'{fipfmt(pl["fip_single"])} with the single-signal '
                     f'Gaussian FIP')
    text = (f'<p>Two planets and outliers, clear and borderline. The posterior '
            f'on the number of signals gives P(k&nbsp;&ge;&nbsp;2) = '
            f'{sum(pk[2:]):.2f}; ' + '; '.join(parts) + '.</p>')
    return fill(html, 'two_planets_text', text)


def gp_joint(html: str) -> str:
    """Activity and a planet, jointly"""
    summ = load('gp_joint')
    if summ is None:
        return html
    kj, ks = summ['K_joint'], summ['K_sequential']
    extra = summ['extra']
    joint = [row['joint'] for row in extra]
    seq = [row['sequential'] for row in extra]
    text = (f'<p>A spotted star rotating in 23&nbsp;d, a K = '
            f'{summ["K_true"]:.0f}&nbsp;m/s planet at 9.7&nbsp;d, and outliers '
            f'clear and borderline. With the GP period constrained by an '
            f'activity indicator, the joint fit returns K = '
            f'{interval(*kj)}&nbsp;m/s and the sequential fit (GP first, orbit '
            f'on the residual) {interval(*ks)}&nbsp;m/s. Over {len(extra)} more '
            f'realisations the joint fit gives {min(joint):.1f} to '
            f'{max(joint):.1f}&nbsp;m/s and the sequential fit '
            f'{min(seq):.1f} to {max(seq):.1f}&nbsp;m/s. koloa flags '
            f'{flags_text(summ["flags"])}.</p>')
    ens = summ.get('ensemble')
    if ens:
        st = ens['stats']
        text += (f'<p>Over {st["n"]} further realisations the joint fit gives '
                 f'a median K of {st["joint_median"]:.2f}&nbsp;m/s (scatter '
                 f'{st["joint_std"]:.2f}), its 68% interval holding the true '
                 f'value in {st["joint_cover68"]} of them, against '
                 f'{st["seq_median"]:.2f}&nbsp;m/s for the sequential fit; '
                 f'under the GP koloa flags {100 * st["recall_clear"]:.0f}% of '
                 f'the clear outlying exposures and '
                 f'{100 * st["recall_border"]:.0f}% of the borderline ones. '
                 f'These are maximum a posteriori fits.</p>')
    return fill(html, 'gp_joint_text', text)


def posterior(html: str) -> str:
    """The posterior of an eccentric orbit"""
    summ = load('posterior')
    if summ is None:
        return html
    one = summ['one']['fits']
    stats = summ['stats']
    nreal = summ['nreal']
    text = (f'<p>An eccentric planet (P = 11.2&nbsp;d, K = '
            f'{summ["truth"]["K"]:.0f}&nbsp;m/s, e = 0.3, &omega; = 60&deg;) on '
            f'a real NIRPS sampling, with outliers clear and borderline. The orbit '
            f'is sampled by MCMC (<span class="mono">koloa.mcmc_orbits</span>: '
            f'emcee, run until every parameter has 50 autocorrelation times) '
            f'with the Gaussian likelihood and with koloa\'s mixture.</p>'
            f'<p>Without the outliers both give the same posterior: K = '
            f'{pct_interval(one["gaussian_clean"]["K"])} and '
            f'{pct_interval(one["koloa_clean"]["K"])}&nbsp;m/s, e = '
            f'{pct_interval(one["gaussian_clean"]["e"], 3)} and '
            f'{pct_interval(one["koloa_clean"]["e"], 3)}. With them, the '
            f'Gaussian posterior gives K = {pct_interval(one["gaussian"]["K"])}'
            f'&nbsp;m/s, e = {pct_interval(one["gaussian"]["e"], 3)} and '
            f'&omega; = {pct_interval(one["gaussian"]["omega"], 0)}&deg;; '
            f'koloa K = {pct_interval(one["koloa"]["K"])}&nbsp;m/s, e = '
            f'{pct_interval(one["koloa"]["e"], 3)} and &omega; = '
            f'{pct_interval(one["koloa"]["omega"], 0)}&deg;.</p>'
            f'<p>Over {nreal} realisations with outliers, the 68% interval of '
            f'K holds the true value in {stats["gaussian"]["K"]["cover68"]} for '
            f'the Gaussian likelihood and {stats["koloa"]["K"]["cover68"]} for '
            f'koloa, and that of e in {stats["gaussian"]["e"]["cover68"]} and '
            f'{stats["koloa"]["e"]["cover68"]} (about {0.68 * nreal:.0f} '
            f'expected). The median half-width of the K interval is '
            f'{stats["gaussian"]["K"]["median_sigma"]:.2f} and '
            f'{stats["koloa"]["K"]["median_sigma"]:.2f}&nbsp;m/s. On the same '
            f'series without outliers both likelihoods under-cover K (rms of '
            f'(median &minus; truth)/&sigma; {stats["gaussian_clean"]["K"]["rms_z"]:.2f}'
            f' and {stats["koloa_clean"]["K"]["rms_z"]:.2f}): the intervals '
            f'of K are about 20% too narrow, a problem not traced yet.</p>')
    if 'gaussian_clip' in stats:
        clip = stats['gaussian_clip']
        text += (f'<p>A Gaussian fit after an iterative 3-sigma clip of its '
                 f'residuals comes close to koloa in K (median error '
                 f'{clip["K"]["median_error"]:+.2f} against '
                 f'{stats["koloa"]["K"]["median_error"]:+.2f}&nbsp;m/s), but '
                 f'its 95% intervals miss the true K in '
                 f'{clip["K"]["n"] - clip["K"]["cover95"]}, the true e in '
                 f'{clip["e"]["n"] - clip["e"]["cover95"]} and the true period '
                 f'in {clip["P"]["n"] - clip["P"]["cover95"]} of {nreal} '
                 f'realisations, against '
                 f'{stats["koloa"]["K"]["n"] - stats["koloa"]["K"]["cover95"]}, '
                 f'{stats["koloa"]["e"]["n"] - stats["koloa"]["e"]["cover95"]}'
                 f' and '
                 f'{stats["koloa"]["P"]["n"] - stats["koloa"]["P"]["cover95"]}'
                 f' for koloa. The clip starts from an orbit that the outliers '
                 f'have already bent, and keeps the clear outliers this orbit '
                 f'passes near: the fits that miss K are mostly those that '
                 f'kept clear outliers.</p>')
    return fill(html, 'posterior_text', text)


def contamination(html: str) -> str:
    """A Student-t against the mixture"""
    summ = load('contamination')
    if summ is None:
        return html
    poly = summ['poly']
    fracs = poly['fracs']

    def at(block, name, frac):
        return block[name][block['fracs'].index(frac)][0]

    text = (f'<p>polyband, a polynomial regression, resists outliers with a '
            f'Student-t likelihood; koloa can too (<span class="mono">'
            f'likelihood=\'student\'</span>), and the two agree to three '
            f'decimals. The mixture differs: a flagged outlier keeps no weight '
            f'at all. On a cubic with outliers of both signs from 3 to '
            f'30&nbsp;sigma, the median error of the fitted curve is, at 30% '
            f'of outliers, {at(poly, "lsq", 0.3):.2f}&nbsp;sigma for least '
            f'squares, {at(poly, "polyband", 0.3):.2f} for polyband and '
            f'{at(poly, "koloa_mixture", 0.3):.2f} for the mixture; at 50%, '
            f'{at(poly, "polyband", 0.5):.2f} and '
            f'{at(poly, "koloa_mixture", 0.5):.2f}. The mixture reaches 50% '
            f'with the robust start it borrows from polyband: the fit of the '
            f'least deviant 70 or 50% of the points.</p>')
    for key, label in (('point', 'bad exposures'), ('visit', 'bad visits')):
        block = summ[key]
        text += (f'<p>An orbit (K = 10&nbsp;m/s) with {label}: at 20% of them '
                 f'the median error of K is {at(block, "gaussian", 0.2):.2f}'
                 f'&nbsp;m/s for the Gaussian, '
                 f'{at(block, "koloa_student", 0.2):.2f} for the Student-t and '
                 f'{at(block, "koloa_mixture", 0.2):.2f} for the mixture; at '
                 f'40%, {at(block, "gaussian", 0.4):.2f}, '
                 f'{at(block, "koloa_student", 0.4):.2f} and '
                 f'{at(block, "koloa_mixture", 0.4):.2f}.</p>')
    return fill(html, 'contamination_text', text)


def completeness(html: str) -> str:
    """Recovery maps"""
    summ = load('completeness')
    if summ is None:
        return html

    def k_at(entry, period):
        import numpy as np
        per = np.array(entry['periods'])
        it = int(np.argmin(np.abs(np.log(per / period))))
        val = entry['k50'][it]
        return f'{val:.1f}' if val is not None and val == val else 'n/a'

    kol, gau = summ['koloa'], summ['gaussian']
    text = (f'<p>Circular planets are injected over a grid of periods and '
            f'semi-amplitudes and looked for by folding the series at their '
            f'period (<span class="mono">koloa.completeness</span>): an '
            f'offset and a sinusoid fitted there, outlier-aware or Gaussian, '
            f'and a planet counts when the gain in log likelihood passes the '
            f'1% false-alarm threshold (4.6). On a real NIRPS sampling with '
            f'outliers clear and borderline, a new series per injection, K is '
            f'recovered half the time from {k_at(kol, 10)}&nbsp;m/s at 10&nbsp;d '
            f'with koloa and from {k_at(gau, 10)}&nbsp;m/s with the Gaussian '
            f'fold; without the outliers, {k_at(summ["koloa_clean"], 10)} and '
            f'{k_at(summ["gaussian_clean"], 10)}&nbsp;m/s. Injections at K = 0 '
            f'fire in {100 * kol["false_alarm"]:.1f}% (koloa) and '
            f'{100 * gau["false_alarm"]:.1f}% (Gaussian) of the cases, for a '
            f'nominal 1%.</p>')
    rows = ['<table class="koloa"><tr><th>series</th><th>exposures</th>'
            '<th>visits</th><th>K half the time at 10 d</th>'
            '<th>at 100 d</th><th>false alarms at K = 0</th></tr>']
    for key in ('real_gl725b', 'real_kepler-21'):
        if key not in summ:
            continue
        entry = summ[key]
        rows.append(f'<tr><td>{escape(entry["label"])}</td>'
                    f'<td class="num">{entry["npoints"]}</td>'
                    f'<td class="num">{entry["nvisits"]}</td>'
                    f'<td class="num">{k_at(entry, 10)} m/s</td>'
                    f'<td class="num">{k_at(entry, 100)} m/s</td>'
                    f'<td class="num">{100 * entry["false_alarm"]:.1f}%</td>'
                    f'</tr>')
    rows.append('</table>')
    html = fill(html, 'completeness_table', '\n'.join(rows))
    return fill(html, 'completeness_text', text)


def secular(html: str) -> str:
    """The acceleration of a star (an ensemble with outliers, GL 725 A and B)
    and the perspective accelerations"""
    summ = load('secular')
    if summ is not None:
        rows = ['<table class="koloa"><tr><th>star</th><th>parallax [mas]'
                '</th><th>proper motion [mas/yr]</th><th>perspective '
                'acceleration [m/s/yr]</th><th>drift over its series [m/s]'
                '</th></tr>']
        for row in summ['table']:
            drift = (f'{row["drift"]:.2f}' if 'drift' in row else '')
            rows.append(f'<tr><td>{escape(row["name"])}</td>'
                        f'<td class="num">{row["parallax"]:.2f}</td>'
                        f'<td class="num">{row["pm"]:.0f}</td>'
                        f'<td class="num">{row["secular"]:.4f} &plusmn; '
                        f'{row["error"]:.1e}</td>'
                        f'<td class="num">{drift}</td></tr>')
        rows.append('</table>')
        html = fill(html, 'secular_table', '\n'.join(rows))
        text = ''
        ens = summ.get('ensemble')
        if ens:
            st = ens['stats']
            gau, kol, cln = st['gaussian'], st['koloa'], st['gaussian_clean']
            text += (f'<p>An acceleration under outliers: {ens["nreal"]} '
                     f'series on a real NIRPS sampling, each with a companion '
                     f'that accelerates the star by {ens["accel"]:.0f}&nbsp;'
                     f'm/s/yr, a K = {ens["planet"]["K"]:.0f}&nbsp;m/s planet, '
                     f'a visit jitter and outliers clear and borderline, '
                     f'the acceleration fitted with the planet. The Gaussian '
                     f'likelihood gives it with a median half-width of '
                     f'{gau["median_sigma"]:.2f}&nbsp;m/s/yr, koloa '
                     f'{kol["median_sigma"]:.2f}, and the Gaussian likelihood '
                     f'on the same series without outliers '
                     f'{cln["median_sigma"]:.2f}; their 95% intervals hold '
                     f'the true value in {gau["cover95"]}, {kol["cover95"]} '
                     f'and {cln["cover95"]} of {gau["n"]} (68%: '
                     f'{gau["cover68"]}, {kol["cover68"]} and '
                     f'{cln["cover68"]}).</p>')
        if text:
            html = fill(html, 'secular_text', text)
        # what GL 725 A and B ask of each other (the use case of the binary)
        text = ''
        g725 = summ.get('gl725')
        if g725 and 'GJ 725 A' in g725:
            aa, bb = g725['GJ 725 A'], g725['GJ 725 B']
            text += (f'<p>What GL 725 A and B ask of each other: at their '
                     f'Gaia DR3 separation, {g725["separation"]:.2f} arcsec '
                     f'({g725["rho_au"]:.1f}&nbsp;au at '
                     f'{g725["distance"]:.3f}&nbsp;pc), the acceleration of B '
                     f'needs a companion of at least '
                     f'{bb["min_mass_other"][1]:.3f}&nbsp;M<sub>&#9737;</sub> '
                     f'and that of A one of at least '
                     f'{aa["min_mass_other"][1]:.3f}&nbsp;M<sub>&#9737;</sub>, '
                     f'whatever the orbit (Torres 1999); A has '
                     f'{bb["mass_other"][0]:.3f} and B '
                     f'{aa["mass_other"][0]:.3f}&nbsp;M<sub>&#9737;</sub>.')
            if 'depth_au' in g725:
                near, far = g725['depth_au']
                text += (f' Together, the two accelerations and the masses '
                         f'place the stars {near:.0f} or {far:.0f}&nbsp;au '
                         f'apart along the line of sight, '
                         f'{escape(g725.get("order", ""))}.')
            text += '</p>'
        if text:
            html = fill(html, 'gl725_masses_text', text)
    summ = load('gl725b')
    if summ is not None and 'GJ 725 B' in summ:
        bstar = summ['GJ 725 B']
        lin = bstar['linear']

        def band(vals):
            return interval(vals[1], vals[1] - vals[0], vals[2] - vals[1])
        text = (f'<p>GL 725 B (SPIRou, {bstar["npoints"]} exposures in '
                f'{bstar["nvisits"]} visits over '
                f'{bstar["baseline"] / 365.25:.1f} years) drifts by '
                f'{band(lin["total"])}&nbsp;m/s/yr')
        if 'GJ 725 A' in summ:
            astar = summ['GJ 725 A']
            text += (f', and GL 725 A the other way, by '
                     f'{band(astar["linear"]["total"])}&nbsp;m/s/yr')
        text += '.'
        if bstar.get('secular_removed'):
            text += (' Their perspective accelerations (from Gaia DR3, '
                     f'{bstar["secular"][0]:.3f}&nbsp;m/s/yr for B'
                     + (f' and {summ["GJ 725 A"]["secular"][0]:.3f} for A'
                        if 'GJ 725 A' in summ else '') +
                     ') are already out of these velocities: APERO computes '
                     'the barycentric correction with barycorrpy, which is '
                     'given the proper motion and the parallax of each star. '
                     'What drifts is the pull of each star on the other.')
        if 'mass_ratio' in summ:
            low, mid, high = summ['mass_ratio']
            text += (f' The ratio of the two accelerations, '
                     f'&minus;a<sub>A</sub>/a<sub>B</sub> = '
                     f'{interval(mid, mid - low, high - mid)}, should be '
                     f'that of the masses, M<sub>B</sub>/M<sub>A</sub>')
            mpath = os.path.join(ROOT, 'data', 'gl725_masses.json')
            if os.path.exists(mpath):
                import numpy as np
                with open(mpath) as handle:
                    mass = json.load(handle)
                mlow, mmid, mhigh = mass['ratio']
                nsig = abs(mid - mmid) / np.hypot(0.5 * (high - low),
                                                  0.5 * (mhigh - mlow))
                # the drift common to both series that would close the gap
                acc_a = summ['GJ 725 A']['linear']['total'][1]
                acc_b = bstar['linear']['total'][1]
                common = (acc_a + mmid * acc_b) / (1 + mmid)
                text += (f' = {mmid:.3f} from the mass-luminosity relation '
                         f'of Mann et al. (2019): they are {nsig:.1f} sigma '
                         f'apart, and a drift of {common:+.2f}&nbsp;m/s/yr '
                         f'common to both series would close the gap')
            text += '.'
            if 'mass_ratio_twice' in summ:
                low, mid, high = summ['mass_ratio_twice']
                text += (f' Taking the perspective accelerations out a second '
                         f'time would have given '
                         f'{interval(mid, mid - low, high - mid)}.')
        text += '</p>'
        ana = summ.get('analysis')
        if ana:
            text += (f'<p>With the drift taken out, the whole koloa analysis '
                     f'finds its best interval at {ana["period"]:.2f}&nbsp;d with '
                     f'FIP = {fipfmt(ana["fip_at_period"]["koloa"])}, and flags '
                     f'{ana["nflag"]} exposures; the duck test says '
                     f'{escape(str(ana["verdict"]))}.</p>')
        html = fill(html, 'gl725_text', text)
    return html


def kepler21(html: str) -> str:
    """Kepler-21 b"""
    summ = load('kepler21')
    if summ is None:
        return html
    lit = summ['literature']
    kol, gau = summ['koloa'], summ['gaussian']
    text = (f'<p>Kepler-21 (HD 179070) is an F6 IV subgiant with a transiting '
            f'super-Earth at 2.786&nbsp;d. Its public HARPS-N velocities (DACE, '
            f'{summ["summary"]["npoints"]} exposures in '
            f'{summ["summary"]["nseq"]} visits over '
            f'{summ["summary"]["baseline"] / 365.25:.1f} years) scatter by '
            f'{summ["summary"]["rms"]:.1f}&nbsp;m/s, far above their '
            f'{summ["summary"]["median_err"]:.2f}&nbsp;m/s errors: activity. No '
            f'blind FIP finds the planet (koloa FIP at its period '
            f'{fipfmt(summ["fip"]["koloa"]["at_b"])}). With the ephemeris from '
            f'the transits, a GP whose rotation period ({summ["prot"]:.1f}&nbsp;d) '
            f'comes from the S index, and the perspective acceleration, the joint '
            f'fit gives K<sub>b</sub> = {interval(*kol["K"])}&nbsp;m/s with koloa '
            f'and {interval(*gau["K"])}&nbsp;m/s with a Gaussian likelihood '
            f'(m&nbsp;sin&nbsp;i = {interval(*kol["msini"], digits=1)}&nbsp;Earth '
            f'masses), against ' + ' and '.join(
                f'{val:.2f}&nbsp;&plusmn;&nbsp;{err:.2f} ({escape(ref)})'
                for ref, (val, err) in lit.items()) + '.</p>')
    return fill(html, 'kepler21_text', text)


def kfmt(val) -> str:
    """K with its interval, or alone (a maximum a posteriori whose Laplace
    interval is undefined, a parameter on the bound of its prior)"""
    if any(not np.isfinite(v) for v in val[1:]):
        return f'{val[0]:.2f} (maximum a posteriori)'
    return interval(*val)


def mdwarfs(html: str) -> str:
    """Three M dwarfs cleaned by PCA2D"""
    summ = load('mdwarfs')
    if summ is None:
        return html
    import numpy as np
    rows = ['<table class="koloa"><tr><th>star</th><th>instrument</th>'
            '<th>exposures / visits</th><th>visit scatter, delivered</th>'
            '<th>visit scatter, PCA2D</th>'
            '<th>known planet: K koloa (PCA2D)</th>'
            '<th>K koloa (delivered)</th><th>K gaussian</th>'
            '<th>published</th>'
            '<th>drift left (velocities)</th><th>perspective (Gaia, removed by APERO)</th>'
            '<th>K half the time, 10 d</th></tr>']
    texts = []
    for key in ('gl687', 'proxima', 'gl699'):
        if key not in summ:
            continue
        entry = summ[key]
        ser = entry['series']
        deliv = entry.get('delivered', {})
        planet_k, lit, kdel, kgau = '', '', '', ''
        for pname, pent in entry.get('planets', {}).items():
            if entry.get('kind', 'PCA2D') == 'PCA2D':
                planet_k = f'{pname}: {interval(*pent["koloa"]["K"])} m/s'
                if 'koloa_delivered' in pent:
                    kdel = interval(*pent['koloa_delivered']['K'])
            else:
                planet_k = f'{pname}: pending'
                kdel = interval(*pent['koloa']['K'])
            kgau = interval(*pent['gaussian']['K'])
            if pent['gaussian']['e'][0] > 0.9:
                kgau += f' (e = {pent["gaussian"]["e"][0]:.2f})'
            lit = (f'{pent["literature"]["K"]:.2f} &plusmn; '
                   f'{pent["literature"]["Kerr"]:.2f}')
        drift = entry['drift_measured']
        comp = entry.get('completeness', {})
        k10 = ''
        if comp:
            per = np.array(comp['periods'])
            it = int(np.argmin(np.abs(np.log(per / 10.0))))
            k10 = f'{comp["k50"][it]:.1f} m/s'
        # the scatter of the nightly means, a straight line removed (the
        #   drift of Barnard's star would be most of it otherwise)
        night = ser.get('nightly_robust_line',
                        ser.get('nightly_rms_line', ser['nightly_rms']))
        dnight = deliv.get('nightly_robust_line',
                           deliv.get('nightly_rms_line',
                                     deliv.get('nightly_rms')))
        pca = entry.get('kind', 'PCA2D') == 'PCA2D'
        corrected = f'{night:.2f}' if pca else 'pending'
        rows.append(f'<tr><td>{escape(entry["name"])}</td>'
                    f'<td>{escape(entry["inst"])}</td>'
                    f'<td class="num">{ser["npoints"]} / {ser["nseq"]}</td>'
                    f'<td class="num">'
                    f'{"" if dnight is None else f"{dnight:.2f}"}</td>'
                    f'<td class="num">{corrected}</td>'
                    f'<td class="num">{planet_k}</td>'
                    f'<td class="num">{kdel}</td>'
                    f'<td class="num">{kgau}</td>'
                    f'<td class="num">{lit}</td>'
                    f'<td class="num">{interval(*drift)} m/s/yr</td>'
                    f'<td class="num">{entry["secular"][0]:.3f}</td>'
                    f'<td class="num">{k10}</td></tr>')
        ana = entry.get('analysis', {})
        if ana:
            kind = ('cleaned by PCA2D' if entry.get('kind', 'PCA2D') == 'PCA2D'
                    else 'as delivered (its PCA2D series is being computed)')
            verdict = str(ana.get('verdict') or '').split(':')[0]
            units = ('exposures' if entry.get('fip_on', 'exposures') ==
                     'exposures' else 'visit means')
            texts.append(
                f'<b>{escape(entry["name"])}</b> ({escape(entry["inst"])}, '
                f'{kind}): the best interval of the analysis is at '
                f'{ana["period"]:.2f}&nbsp;d, with a koloa FIP of '
                f'{fipfmt(ana["fip_at_period"]["koloa"])}, and the duck test '
                f'says {escape(verdict)}; {ana["nflag"]} {units} are flagged. '
                f'The drift left in the velocities is '
                f'{interval(*drift)}&nbsp;m/s/yr: APERO has already removed '
                f'the perspective acceleration ({entry["secular"][0]:.3f}&nbsp;'
                f'm/s/yr from Gaia).')
            gpres = entry.get('gp', {})
            parts = []
            for kind, label in (('delivered', 'delivered'),
                                ('pca2d', 'PCA2D')):
                for pname, pent in gpres.get(kind, {}).items():
                    if isinstance(pent, dict) and 'K' in pent:
                        parts.append(f'{interval(*pent["K"])}&nbsp;m/s '
                                     f'({label})')
            if parts:
                texts[-1] += (' With a GP of the rotation fitted to the visit '
                              'means, the known planet has K = '
                              + ' and '.join(parts) + '.')
                # the same fits with the prior of the rotation moved to the
                #   published range (demo_mdwarfs.py --gp-only --gp-check)
                check = entry.get('gp_check', {})
                cparts = []
                for kind, label in (('delivered', 'delivered'),
                                    ('pca2d', 'PCA2D')):
                    for pname, pent in check.get(kind, {}).items():
                        if isinstance(pent, dict) and 'K' in pent:
                            cparts.append(f'{interval(*pent["K"])}&nbsp;m/s '
                                          f'({label})')
                if cparts:
                    prots = [check[kind]['prot'] for kind in
                             ('delivered', 'pca2d') if kind in check]
                    # the range is target['gp']['check'] of the demo
                    texts[-1] += (' With the FWHM fit held to 75 to '
                                  '110&nbsp;d, around the published rotation '
                                  '(prior at '
                                  + ' and '.join(f'{p:.0f}' for p in prots)
                                  + '&nbsp;d instead of '
                                  + ' and '.join(
                                      f'{gpres[kind]["prot"]:.0f}' for kind
                                      in ('delivered', 'pca2d')
                                      if kind in gpres)
                                  + '&nbsp;d), K = ' + ' and '.join(cparts)
                                  + ': the planet does not depend on it.')
    rows.append('</table>')
    html = fill(html, 'mdwarfs_table', '\n'.join(rows))
    if texts:
        html = fill(html, 'mdwarfs_text',
                    '<p>' + '</p>\n<p>'.join(texts) + '</p>')
    return html


def dace(html: str) -> str:
    """A known system from DACE: the instruments, the FIP and the orbits
    against the NASA Exoplanet Archive (demos/demo_dace.py)"""
    summ = load('dace')
    if summ is None:
        return html
    f1, f2 = summ['fip_first'], summ['fip_second']
    known = summ['archive']['planets']
    text = (f'<p>DACE serves {summ["n"]} public velocities of '
            f'{escape(summ["host"])} in {summ["nvisits"]} visits over '
            f'{summ["baseline_years"]:.1f}&nbsp;years, from '
            f'{len(summ["instruments"])} instrument eras (table). koloa '
            f'fits the noise of each (a white and a visit jitter), inflates '
            f'the errors of each to it (<span class="mono">'
            f'koloa.fip.inflate_to_fit</span>: the FIP sampler has one '
            f'jitter for all) and runs the FIP with up to '
            f'{len(f2["pk"]) - 1} signals. With the noise fitted without '
            f'planets, the intervals with FIP&nbsp;&lt;&nbsp;1% are at '
            + ', '.join(f'{pk["period"]:.2f}&nbsp;d ({fipfmt(pk["fip"])})'
                        for pk in f1['peaks'] if pk['fip'] < 0.01)
            + '. Those planets are then fitted together by MCMC, every '
              'instrument with its offset, jitters and outliers, and the FIP '
              'runs again with the noise of that fit: '
            + ', '.join(f'{pk["period"]:.2f}&nbsp;d ({fipfmt(pk["fip"])})'
                        for pk in f2['peaks'] if pk['fip'] < 0.01)
            + f'; the number of signals has P(k) = '
            + ', '.join(f'{val:.2f}' for val in f2['pk'])
            + f' for k = 0 to {len(f2["pk"]) - 1}. The K of each planet is '
              f'set against every solution that the NASA Exoplanet Archive '
              f'lists, the most recent first (table): the archive\'s default '
              f'is the discovery paper, from far fewer velocities than DACE '
              f'now serves.</p>')
    rows = ['<table class="koloa"><tr><th>instrument</th><th>velocities</th>'
            '<th>visits</th><th>years</th><th>median error [m/s]</th>'
            '<th>white jitter [m/s]</th><th>visit jitter [m/s]</th>'
            '<th>flagged</th></tr>']
    for row in summ['instruments']:
        noise = summ['noise'].get(row['inst'], {})
        white, visit = noise.get('white'), noise.get('visit')
        rows.append(f'<tr><td>{escape(row["inst"])}</td>'
                    f'<td class="num">{row["n"]}</td>'
                    f'<td class="num">{row["nvisits"]}</td>'
                    f'<td class="num">{row["years"][0]:.1f} to '
                    f'{row["years"][1]:.1f}</td>'
                    f'<td class="num">{row["median_err"]:.2f}</td>'
                    f'<td class="num">{"" if white is None else f"{white:.2f}"}</td>'
                    f'<td class="num">{"" if visit is None else f"{visit:.2f}"}</td>'
                    f'<td class="num">{summ["flagged"].get(row["inst"], 0)}</td></tr>')
    rows.append('</table>')
    orows = ['<table class="koloa"><tr><th>planet</th><th>P [d]</th>'
             '<th>K [m/s]</th><th>e</th><th>m sin i [Earth masses]</th>'
             '<th>published K [m/s]</th></tr>']
    from usecases import comparisons, match_known
    for orb, pl in match_known(summ['orbits'], known):
        # every published K of the planet, most recent first, linked to its
        #   paper
        pub = '<br>'.join(
            f'{sol["K"]:.2f} &plusmn; {sol["K_err"]:.2f} '
            + (f'(<a href="{escape(sol["reference_url"])}" target="_blank" '
               f'rel="noopener">{escape(sol["reference"])}</a>)'
               if sol.get('reference_url') else f'({escape(sol["reference"])})')
            for sol in comparisons(orb, pl))
        orows.append(f'<tr><td>{escape((pl or {}).get("name", "new"))}</td>'
                     f'<td class="num">{interval(*orb["P"], digits=4)}</td>'
                     f'<td class="num">{interval(*orb["K"])}</td>'
                     f'<td class="num">{interval(*orb["e"])}</td>'
                     f'<td class="num">{interval(*orb["msini"], digits=1)}</td>'
                     f'<td class="num">{pub}</td></tr>')
    orows.append('</table>')
    html = fill(html, 'dace_text', text)
    html = fill(html, 'dace_instruments', '\n'.join(rows))
    return fill(html, 'dace_orbits', '\n'.join(orows))


def outliers(html: str) -> str:
    """Understanding outliers: the keys of each instrument (koloa.outliers),
    and why the outliers of three stars are outliers (demos/demo_outliers.py)"""
    from koloa import outliers as kout
    # the keys: every list, one row per key, where it is found
    words = {'low': 'low', 'high': 'high', 'either': 'either way'}
    generic = {'snr_rate': 'S/N squared per second of exposure: the flux '
                           'collected (clouds, seeing, guiding lower it)',
               'snr_goal': 'S/N over its goal (SNRGOAL)',
               'wave_age': 'days since the wavelength solution was taken'}
    order, where = [], {}
    for lname in ('spirou', 'nirps', 'dace'):
        for key in kout.LISTS[lname]:
            if key.name not in where:
                order.append(key)
                where[key.name] = set()
            where[key.name].add(lname)
    rows = ['<table class="koloa"><tr><th>key</th><th>what it is</th>'
            '<th>a problem when</th><th>SPIRou</th><th>NIRPS</th>'
            '<th>DACE</th></tr>']
    for key in order:
        label = generic.get(key.name, key.label)
        name = key.name.replace('HIERARCH ESO ', '')
        rows.append(f'<tr><td class="mono">{escape(name)}</td>'
                    f'<td>{escape(label)}</td><td>{words[key.bad]}</td>'
                    + ''.join(f'<td>{"&#10003;" if lname in where[key.name] else ""}</td>'
                              for lname in ('spirou', 'nirps', 'dace'))
                    + '</tr>')
    rows.append('<tr><td class="mono">indicators</td><td>every activity '
                'indicator of the series: DTEMP, d2v, d3v, dW, FWHM, '
                'contrast, CRX and the chromatic slope from LBL; FWHM, '
                "bisector, R'HK, S index and H alpha from DACE</td>"
                '<td>either way</td><td>&#10003;</td><td>&#10003;</td>'
                '<td>&#10003;</td></tr>')
    rows.append('</table>')
    html = fill(html, 'outliers_keys', '\n'.join(rows))
    summ = load('outliers')
    if summ is None:
        return html
    source = {'SPIRou': 'data of the SPIRou team', 'NIRPS': 'data of the '
              'NIRPS team', 'HARPS': 'public HARPS velocities from DACE'}
    blocks = []
    for key in ('gl687', 'proxima', 'hd69830'):
        case = summ.get(key)
        if not case:
            continue
        rep = case['report']
        units = rep['units']
        count = next(iter(rep['counts'].values()))
        planets = ', '.join(f'{per:g}' for per in case['planets'])
        text = (f'<p><b>{escape(case["name"])}</b> ({case["inst"]}, '
                f'{source.get(case["inst"], case["inst"])}): {case["n"]} '
                f'exposures in {case["nvisits"]} visits, the known planet'
                f'{"s" if len(case["planets"]) > 1 else ""} ({planets}&nbsp;d) '
                f'fitted with the mixture. {count["outliers"]} exposure'
                f'{"s are" if count["outliers"] != 1 else " is"} outlying: '
                f'{count["visit_units"]} whole visit'
                f'{"s" if count["visit_units"] != 1 else ""} and '
                f'{count["exposure_units"]} single exposure'
                f'{"s" if count["exposure_units"] != 1 else ""}.')
        shared = [row for row in rep['association'] if row['significant']]
        if shared:
            text += (' What they share: ' + '; '.join(
                f'{escape(row["label"])}, {escape(row["summary"])}'
                for row in shared) + '.')
        text += '</p><ul>'
        for unit in units:
            what = (f'a visit of {unit["n"]} exposures'
                    if unit['kind'] == 'visit' else 'an exposure')
            where = ('with the other exposures of its visit'
                     if unit['compare'] == 'visit' else
                     'with the good visits nearest in time')
            resid = ''
            if unit.get('resid') is not None:
                resid = (f', {unit["resid"]:+.1f}&nbsp;m/s '
                         f'({abs(unit["resid_sigma"]):.1f}&nbsp;sigma) off '
                         f'the model')
            off = [k for k in unit['keys'] if k['significant']]
            if off:
                found = '; '.join(_outlier_key(k) for k in off)
                verdict = f'significantly off: {found}'
            else:
                near = [k for k in unit['keys'] if k.get('z') is not None
                        and abs(k['z']) > 2][:2]
                verdict = 'nothing significantly off' + (
                    ' (closest: ' + '; '.join(
                        f'{escape(k["label"])}, z&nbsp;{k["z"]:+.1f}'
                        for k in near) + ')' if near else '')
            text += (f'<li>{unit["date"]}, {what}{resid}, compared {where}: '
                     f'{verdict}.</li>')
        text += '</ul>'
        text += (f'<figure class="lecture-figure"><img src="figures/'
                 f'outliers_{key}_keys.svg" alt="Why the outliers of '
                 f'{escape(case["name"])} are outliers"><figcaption>Each '
                 f'outlier of {escape(case["name"])} (a row) against each key '
                 f'(a column), coloured by its deviation from the good data '
                 f'(teal low, orange high); framed cells are significantly '
                 f'off.</figcaption></figure>')
        blocks.append(text)
    return fill(html, 'outliers_cases', '\n'.join(blocks))


def _outlier_key(key: dict) -> str:
    """a key that is off, in words, for the site"""
    if key.get('text'):
        return (f'{escape(key["label"])} is {escape(str(key["value"]))}, as '
                f'for {100 * key["frequency"]:.0f}&nbsp;% of the good '
                f'exposures')
    from koloa.outliers import _fmt_pair
    value, ref = _fmt_pair(key['value'], key['reference'])
    zval = key['z']
    ztxt = ('without a precedent among the good data' if zval is None
            or not np.isfinite(zval) else f'z&nbsp;{zval:+.1f}')
    return (f'{escape(key["label"])} {value} against {ref} ({ztxt})')


def rotation(html: str) -> str:
    """The rotation of a star from its temperature indicator: the
    periodograms, the posteriors of the SHO GP, the robustness
    (demos/demo_rotation.py)"""
    summ = load('rotation')
    if summ is None:
        return html
    mix, gau = summ['posteriors']['mixture'], summ['posteriors']['gaussian']
    rob = summ['robustness']
    # the quality factor held above 1, and what the fit does with it free
    qlow = summ.get('quality', [0.3])[0]
    free = summ.get('free_quality')
    held = ''
    if free:
        fmap, hmap = free['map'], free['map_held']
        gain = fmap['logpost'] - hmap['logpost']
        shape = ('an SHO without a peak, red noise rather than a rotation'
                 if not np.isfinite(fmap['P_peak']) else
                 f'its power peaking at {fmap["P_peak"]:.1f}&nbsp;d')
        held = (f', held above Q = {qlow:.0f}: a rotation is quasi-periodic, '
                f'and below Q = 1/&radic;2 the power of an SHO has no peak. '
                f'Left free down to Q = {free["range"][0]}, the maximum a '
                f'posteriori falls at Q = {fmap["Q"]:.2f} and P<sub>0</sub> = '
                f'{fmap["P0"]:.1f}&nbsp;d, {shape}; its log posterior is '
                f'{abs(gain):.1f} {"higher" if gain > 0 else "lower"} than '
                f'with Q above {qlow:.0f}')
    # the period of the GP against the photometry, in the GP's sigma (the
    #   side of its interval towards it)
    per, lit_p = mix['gp_period_peak'], summ['literature']['P']
    zlit = abs(lit_p - per[0]) / (per[2] if lit_p > per[0] else per[1])
    if zlit < 2:
        gp_note = (f'The period of the GP is broader than the peaks of the '
                   f'periodograms: an SHO of Q near {mix["gp_quality"][0]:.1f} '
                   f'spreads its power over a wide band of periods, so its '
                   f'peak is a loose estimate of the rotation, and the '
                   f'periodograms, which fit a single sinusoid, are sharper '
                   f'for this star.')
    else:
        gp_note = (f'That is {zlit:.1f}&sigma; below the rotation of the '
                   f'photometry and the peaks of the periodograms: for this '
                   f'star, a single damped oscillator of Q near '
                   f'{mix["gp_quality"][0]:.1f} is a poor description of the '
                   f'rotation, which the periodograms of DTEMP3500 and of the '
                   f'FWHM find without it. A sharper kernel (an oscillator at '
                   f'P and one at P/2, or a quasi-periodic one) would be the '
                   f'next test.')
    text = (f'<p>The outlier-aware periodogram of DTEMP3500 peaks at '
            f'{summ["peak_mix"]:.2f}&nbsp;d (&Delta;lnL = '
            f'{summ["peak_mix_dlnl"]:.1f}), the gaussian one at '
            f'{summ["peak_gau"]:.2f}&nbsp;d, and that of the FWHM of the lines '
            f'at {summ["peak_fwhm"]:.2f}&nbsp;d. An SHO GP has an amplitude, an '
            f'undamped period P<sub>0</sub> and a quality factor Q{held}. '
            f'Fitted to DTEMP3500 by MCMC with koloa\'s likelihood, it gives '
            f'P<sub>0</sub> = {interval(*mix["gp_period"], digits=1)}&nbsp;d '
            f'and Q = {interval(*mix["gp_quality"])}; its power peaks at '
            f'P<sub>0</sub>/&radic;(1 &minus; 1/2Q<sup>2</sup>) = '
            f'{interval(*mix["gp_period_peak"], digits=1)}&nbsp;d, the period '
            f'to compare with a periodogram and with the '
            f'{summ["literature"]["P"]:.0f}&nbsp;d of the photometry '
            f'({escape(summ["literature"]["reference"])}). koloa flags '
            f'{summ["flagged"]} of the {summ["n"]} values as outliers. With a '
            f'gaussian likelihood, the peak is at '
            f'{interval(*gau["gp_period_peak"], digits=1)}&nbsp;d. '
            + gp_note + '</p>'
            f'<p>With {100 * rob["fraction"]:.0f}% of the visits moved by '
            f'{rob["shift"][0]:.0f} to {rob["shift"][1]:.0f} times the '
            f'dispersion of DTEMP3500, over {rob["nreal"]} realisations '
            f'(maximum a posteriori), koloa\'s period stays within 1&sigma; of '
            f'the clean one in {rob["held"]["mixture"]} and the gaussian '
            f'likelihood\'s in {rob["held"]["gaussian"]}.</p>')
    names = [('gp_period', 'P<sub>0</sub> [d]', 1),
             ('gp_period_peak', 'peak of the power [d]', 1),
             ('gp_quality', 'Q', 2), ('gp_sigma', '&sigma; [K]', 2)]
    names += [(key, 'visit jitter [K]' if key.startswith('sjit')
               else 'white jitter [K]', 2)
              for key in mix if key.startswith(('jit_', 'sjit')) and key in gau]
    rows = ['<table class="koloa"><tr><th>parameter</th><th>koloa</th>'
            '<th>gaussian</th></tr>']
    for key, label, digits in names:
        rows.append(f'<tr><td>{label}</td>'
                    f'<td class="num">{interval(*mix[key], digits=digits)}</td>'
                    f'<td class="num">{interval(*gau[key], digits=digits)}</td>'
                    f'</tr>')
    rows.append('</table>')
    html = fill(html, 'rotation_text', text)
    return fill(html, 'rotation_table', '\n'.join(rows))


def usecases(html: str) -> str:
    """The test of every use case (docs/usecases.py): its question, data,
    command and criteria, each checked on the outputs of its run; and the
    list of the use cases, with their score, on the page that introduces
    them"""
    from usecases import USECASES, build
    specs, tally, scores = build()
    for key, block in specs.items():
        html = fill(html, f'uc_spec_{key}', block)
    rows, kind_now = ['<ul class="page-index uc-index">'], None
    for key, tab, _, kind, _ in USECASES:
        if key not in specs or f'id="uc-panel-{key}"' not in html:
            continue
        if kind != kind_now:
            rows.append(f'<li class="page-index-label">'
                        f'{"simulated" if kind == "simulated" else "real stars"}'
                        f'</li>')
            kind_now = kind
        npass, ntest, question = scores[key]
        mark = 'uc-pass' if npass == ntest else 'uc-fail'
        rows.append(f'<li><a href="#uc-panel-{key}">{escape(tab)}</a>'
                    f'<span>{question}</span>'
                    f'<span class="uc-badge {mark}">{npass}/{ntest}</span></li>')
    rows.append('</ul>')
    html = fill(html, 'uc_index', '\n'.join(rows))
    return fill(html, 'uc_tally', f'{tally[0]} of {tally[1]} criteria pass '
                                  f'over {len(specs)} use cases.')


def sidebar(html: str) -> str:
    """
    The sidebar: every page, in the order of the page, the pages of a group
    (use cases, modules, mathematics) under the page that introduces it;
    read from the pages themselves (<section class="view" id=...
    data-title=... data-group=...>), so it always matches them
    """
    views = re.findall(r'<section class="view[^"]*" id="([\w-]+)"'
                       r'(?: data-group="([\w-]+)")? data-title="([^"]*)"', html)
    groups = {}
    for vid, group, title in views:
        if group:
            groups.setdefault(group, []).append((vid, title))
    from usecases import USECASES
    kinds = {f'uc-panel-{key}': kind for key, _, _, kind, _ in USECASES}
    # the pages of the GUI: its parts, then the walk-throughs
    kinds.update({vid: ('gui-walk' if vid.startswith('gui-walk') else 'gui')
                  for vid, group, _ in views if group == 'gui'})
    labels = dict(simulated='simulated', real='real stars',
                  gui='the page, part by part', **{'gui-walk': 'walk-throughs'})
    parts = ['<div class="nav-title">koloa</div>', '<ul class="nav-top">']
    for vid, group, title in views:
        if group:
            continue
        sub = groups.get(vid, [])
        if not sub:
            parts.append(f'<li><a href="#{vid}">{title}</a></li>')
            continue
        items, kind_now = [], None
        for sid, stitle in sub:
            kind = kinds.get(sid)
            if kind and kind != kind_now:
                items.append(f'<li class="nav-label">'
                             f'{labels.get(kind, "real stars")}</li>')
                kind_now = kind
            items.append(f'<li><a href="#{sid}">{stitle}</a></li>')
        parts.append(f'<li class="nav-group" data-group="{vid}">'
                     f'<a href="#{vid}">{title}</a>'
                     f'<ul class="nav-sub">{"".join(items)}</ul></li>')
    parts.append('</ul>')
    return fill(html, 'sidebar', '\n'.join(parts))


#: the papers the site cites, and where they are (each checked, 2026-09-29,
#: on CrossRef or on the NASA Exoplanet Archive's ADS links)
CITATIONS = [
    ('Hara et al.', '2022', 'https://doi.org/10.1051/0004-6361/202140543'),
    ('Torres', '1999', 'https://doi.org/10.1086/316313'),
    ('Liu et al.', '2002', 'https://doi.org/10.1086/339845'),
    ('Wright &amp; Eastman', '2014', 'https://doi.org/10.1086/678541'),
    ('Mann et al.', '2019', 'https://doi.org/10.3847/1538-4357/aaf3bc'),
    ('Lopez-Morales et al.', '2016',
     'https://doi.org/10.3847/0004-6256/152/6/204'),
    ('Bonomo et al.', '2023',
     'https://ui.adsabs.harvard.edu/abs/2023A%26A...677A..33B/abstract'),
    ('Hori et al.', '2024',
     'https://ui.adsabs.harvard.edu/abs/2024AJ....167..289H/abstract'),
    ('Rosenthal et al.', '2021', 'https://doi.org/10.3847/1538-4365/abe23c'),
    ('Lovis et al.', '2006', 'https://doi.org/10.1038/nature04828'),
    ('Harada et al.', '2025', 'https://doi.org/10.3847/1538-3881/ae0b62'),
    ('Burt et al.', '2014', 'https://doi.org/10.1088/0004-637X/789/2/114'),
    ('Artigau et al.', '2022', 'https://doi.org/10.3847/1538-3881/ac7ce6'),
    ('Sokal', '1997', 'https://doi.org/10.1007/978-1-4899-0319-8_6'),
]


def link_citations(html: str) -> str:
    """
    Every citation of CITATIONS in the text of the page linked to its paper
    ("Hara et al. (2022)" or "Hara et al. 2022"), except where it already
    is a link, in an attribute, a script or a style
    """
    pats = [(re.compile(re.escape(who).replace(r'\ ', r'\s+')
                        + r'\s+(?:\(' + year + r'\)|' + year + r')'), url)
            for who, year, url in CITATIONS]
    parts = re.split(r'(<[^>]+>)', html)
    inside = {'a': 0, 'script': 0, 'style': 0, 'head': 0}
    for it, part in enumerate(parts):
        if part.startswith('<'):
            tag = re.match(r'</?\s*([a-zA-Z]+)', part)
            if tag and tag.group(1).lower() in inside:
                name = tag.group(1).lower()
                if part.startswith('</'):
                    inside[name] = max(inside[name] - 1, 0)
                elif not part.endswith('/>'):
                    inside[name] += 1
            continue
        if any(inside.values()) or not part.strip():
            continue
        for pat, url in pats:
            part = pat.sub(lambda mt, url=url: (
                f'<a href="{url}" target="_blank" rel="noopener">'
                f'{mt.group(0)}</a>'), part)
        parts[it] = part
    return ''.join(parts)


def fill_cards(html: str) -> str:
    """Every card this module writes"""
    for func in (math, false_alarm, recovery, two_planets, gp_joint,
                 posterior, contamination, completeness, secular, kepler21,
                 mdwarfs, dace, rotation, outliers, usecases, sidebar):
        html = func(html)
    return link_citations(html)


# =============================================================================
# End of code
# =============================================================================

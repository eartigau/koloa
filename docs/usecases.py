#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The use cases of the koloa site, one tab each, written as tests.

Every use case says what it asks, on which data, what a correct method must
do (its criteria), and how to run it. Each criterion is checked on the
outputs of its demo (demos/output/) and shown as
PASS or FAIL with the value measured, or as INFO when there is no truth to
check against. docs/page_cards.py writes the uc_spec_<key> blocks of
index.html from here, so the page never shows a number typed by hand.

Created on 2026-09-28

@author: artigau
"""
import json
import os
from html import escape

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEMOS = os.path.join(ROOT, 'demos', 'output')
#: the binomial distribution, linked
BINOMIAL = ('<a href="https://en.wikipedia.org/wiki/Binomial_distribution" '
            'target="_blank" rel="noopener">binomial</a>')


# =============================================================================
# Helpers
# =============================================================================
def load(*parts):
    """A JSON file under the repository, or None"""
    path = os.path.join(ROOT, *parts)
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def demo(name, fname='summary.json'):
    """The summary of a demo, or None"""
    return load('demos', 'output', name, fname)


def half(val):
    """The half-width of a [median, low, high] interval"""
    return 0.5 * (val[1] + val[2])


def pm(val, digits=2):
    """A [median, low, high] interval, in HTML"""
    return (f'{val[0]:.{digits}f}&nbsp;&plusmn;&nbsp;'
            f'{half(val):.{digits}f}')


def fipfmt(val):
    """A FIP or a probability"""
    if val >= 0.095:
        return f'{val:.2f}'
    mant, expo = f'{val:.1e}'.split('e')
    return f'{mant}&times;10<sup>{int(expo)}</sup>'


def binomial(nval, prob, nsig=2.0):
    """The counts expected within nsig of a binomial"""
    mid, sig = nval * prob, np.sqrt(nval * prob * (1 - prob))
    return mid - nsig * sig, mid + nsig * sig


def template_sampling():
    """The real NIRPS sampling on which the simulations are made (times,
    visits and error bars; data/nirps_template.csv)"""
    path = os.path.join(ROOT, 'data', 'nirps_template.csv')
    if not os.path.exists(path):
        return 'a real NIRPS sampling'
    from koloa.data import RVData
    ser = RVData.from_csv(path, name='NIRPS template')
    return (f'a real NIRPS sampling ({ser.n} exposures in {ser.nseq} visits '
            f'over {np.ptp(ser.time):.0f}&nbsp;d)')


def series_text(ser):
    """exposures, visits and baseline of a series summary"""
    nvis = ser.get('nseq', ser.get('nvisits'))
    return (f'{ser["npoints"]} exposures in {nvis} visits over '
            f'{ser["baseline"] / 365.25:.1f}&nbsp;yr')


# =============================================================================
# The use cases: simulated
# =============================================================================
def uc_nofalse():
    """No planet, outliers of both kinds"""
    summ = demo('false_alarm')
    if summ is None:
        return None
    fooled, nval = summ['fooled'], summ['nseed']
    flags = summ['shown']['flags']
    return dict(
        question='Does koloa report a planet where there is none, when some '
                 'exposures and some whole visits are outliers?',
        data=f'Simulated: {nval} series of noise on {template_sampling()}, with '
             f'a visit jitter and outliers of both signs, clear (whole visits '
             f'and single exposures) and borderline. No planet.',
        criteria=[
            ('koloa reports no planet: no interval reaches FIP&nbsp;&lt;&nbsp;1% '
             'in any series.', fooled['koloa'] == 0,
             f'koloa {fooled["koloa"]} of {nval}; for comparison the '
             f'fixed-jitter FIP {fooled["gaussian"]}, the soft clip '
             f'{fooled["soft"]}, the hard clip {fooled["hard"]}'),
            ('It flags the clear outliers and leaves the good exposures '
             'alone.',
             (flags['clear']['flagged'] == flags['clear']['n']
              and flags['good']['flagged'] <= 0.02 * flags['good']['n']),
             f'clear {flags["clear"]["flagged"]} of {flags["clear"]["n"]}, '
             f'good {flags["good"]["flagged"]} of {flags["good"]["n"]} '
             f'(the series shown)')],
        run='python demos/demo_false_alarm.py')


def uc_recovery():
    """A planet behind outliers"""
    summ = demo('recovery')
    if summ is None:
        return None
    per, kval = summ['period'], summ['K']
    fk = summ['fips']['koloa']
    korb = summ['orbits']['koloa']['K']
    flags = summ['flags']
    return dict(
        question='Does koloa find a planet hidden behind outliers, and measure '
                 'its amplitude?',
        data=f'Simulated: a circular planet, K&nbsp;=&nbsp;{kval:.0f}&nbsp;m/s '
             f'at {per:.2f}&nbsp;d, on {template_sampling()}, with a visit jitter '
             f'and outliers, clear and borderline.',
        criteria=[
            ('The best interval of koloa\'s FIP is the planet.',
             abs(fk['best_period'] / per - 1) < 0.01,
             f'{fk["best_period"]:.2f}&nbsp;d'),
            ('Its FIP passes the detection threshold, FIP&nbsp;&lt;&nbsp;1%.',
             fk['at_planet'] < 0.01,
             f'koloa {fipfmt(fk["at_planet"])}; the fixed-jitter FIP '
             f'{fipfmt(summ["fips"]["gaussian"]["at_planet"])}, the hard clip '
             f'{fipfmt(summ["fips"]["hard"]["at_planet"])}'),
            ('The 68% interval of K holds the true K.',
             korb[0] - korb[1] <= kval <= korb[0] + korb[2],
             f'K&nbsp;=&nbsp;{korb[0]:.2f} +{korb[2]:.2f}/&minus;{korb[1]:.2f}'
             f'&nbsp;m/s; the Gaussian orbit '
             f'{summ["orbits"]["gaussian"]["K"][0]:.2f}&nbsp;m/s'),
            ('The clear outliers are flagged, the good exposures are not.',
             (flags['clear']['flagged'] == flags['clear']['n']
              and flags['good']['flagged'] <= 0.02 * flags['good']['n']),
             f'clear {flags["clear"]["flagged"]} of {flags["clear"]["n"]}, '
             f'good {flags["good"]["flagged"]} of {flags["good"]["n"]}')],
        run='python demos/demo_recovery.py')


def uc_posterior():
    """Posteriors of an eccentric orbit under outliers"""
    summ = demo('posterior')
    if summ is None:
        return None
    tru, nval = summ['truth'], summ['nreal']
    kst, gst = summ['stats']['koloa'], summ['stats']['gaussian']
    low, high = binomial(nval, 0.68)
    zlim = 2.0 / np.sqrt(2 * nval)
    return dict(
        question='Do outliers bias the posterior of an eccentric orbit, and '
                 'are koloa\'s intervals honest?',
        data=f'Simulated: an eccentric planet (P&nbsp;=&nbsp;{tru["P"]}&nbsp;d, '
             f'K&nbsp;=&nbsp;{tru["K"]:.0f}&nbsp;m/s, e&nbsp;=&nbsp;{tru["e"]}, '
             f'&omega;&nbsp;=&nbsp;{tru["omega_deg"]:.0f}&deg;) on '
             f'{template_sampling()}, with a visit jitter and outliers, in '
             f'{nval} realisations, each sampled by <a href="https://en.wikipedia.org/wiki/'
             f'Markov_chain_Monte_Carlo" target="_blank" rel="noopener">MCMC</a> with the Gaussian '
             f'likelihood and with koloa\'s mixture.',
        criteria=[
            ('With the outliers, koloa\'s K is unbiased: its median error is '
             'smaller than its median uncertainty.',
             abs(kst['K']['median_error']) < kst['K']['median_sigma'],
             f'median error {kst["K"]["median_error"]:+.2f}&nbsp;m/s, '
             f'uncertainty {kst["K"]["median_sigma"]:.2f}; the Gaussian '
             f'likelihood {gst["K"]["median_error"]:+.2f} and '
             f'{gst["K"]["median_sigma"]:.2f}'),
            ('Its 68% intervals hold the true K and e as often as they should '
             f'({low:.0f} to {high:.0f} of {nval}, {BINOMIAL} 2&sigma;).',
             all(low <= kst[key]['cover68'] <= high for key in ('K', 'e')),
             f'K {kst["K"]["cover68"]}, e {kst["e"]["cover68"]}'),
            ('The errors scatter as the uncertainties say: the rms of '
             '(median &minus; truth)/&sigma; of K is 1 (within '
             f'{zlim:.2f}).', abs(kst['K']['rms_z'] - 1) <= zlim,
             f'{kst["K"]["rms_z"]:.2f}: the intervals of K are too narrow, '
             f'on clean series too (not traced yet)')],
        run='python demos/demo_posterior.py')


def uc_twoplanets():
    """Two planets and outliers"""
    summ = demo('two_planets')
    if summ is None:
        return None
    pk2 = float(sum(summ['pk'][2:]))
    pls = summ['planets']
    desc = ' and '.join(f'K&nbsp;=&nbsp;{pl["K"]}&nbsp;m/s at {pl["P"]}&nbsp;d'
                        for pl in pls)
    return dict(
        question='Does koloa\'s sampler find two planets at once, under '
                 'outliers?',
        data=f'Simulated: two planets, {desc}, on {template_sampling()}, with '
             f'outliers clear and borderline; the FIP with up to three '
             f'signals.',
        criteria=[
            ('It finds two signals: P(k&nbsp;&ge;&nbsp;2)&nbsp;&gt;&nbsp;0.9.',
             pk2 > 0.9, f'{pk2:.2f}'),
            ('Each planet has a lower FIP with koloa than with the '
             'single-signal Gaussian FIP.',
             all(pl['fip_koloa'] < pl['fip_single'] for pl in pls),
             '; '.join(f'{pl["P"]}&nbsp;d: {fipfmt(pl["fip_koloa"])} against '
                       f'{fipfmt(pl["fip_single"])}' for pl in pls)),
            ('Each planet passes FIP&nbsp;&lt;&nbsp;1%.',
             all(pl['fip_koloa'] < 0.01 for pl in pls),
             '; '.join(f'{pl["P"]}&nbsp;d: {fipfmt(pl["fip_koloa"])}'
                       for pl in pls))],
        run='python demos/demo_two_planets.py')


def uc_gpjoint():
    """Activity and a planet, jointly"""
    summ = demo('gp_joint')
    if summ is None:
        return None
    ens, ktrue = summ['ensemble']['stats'], summ['K_true']
    return dict(
        question='Does a joint fit of the activity (a <a href="https://'
                 'en.wikipedia.org/wiki/Gaussian_process" target="_blank" '
                 'rel="noopener">GP</a>) and the orbit '
                 'recover K, where a GP fitted first and an orbit fitted to '
                 'what it leaves does not?',
        data=f'Simulated: a spotted star rotating in 23&nbsp;d, a planet with '
             f'K&nbsp;=&nbsp;{ktrue:.0f}&nbsp;m/s at 9.7&nbsp;d, and outliers '
             f'clear and borderline, in {ens["n"]} realisations; the GP period '
             f'constrained by an activity indicator.',
        criteria=[
            ('The joint fit recovers K: the median over the realisations is '
             'within 20% of the truth.',
             abs(ens['joint_median'] / ktrue - 1) < 0.2,
             f'{ens["joint_median"]:.2f}&nbsp;m/s; its 68% intervals hold the '
             f'truth in {ens["joint_cover68"]} of {ens["n"]}'),
            ('The sequential fit, for comparison, loses the planet into the '
             'GP (median below 80% of the truth).',
             ens['seq_median'] < 0.8 * ktrue,
             f'{ens["seq_median"]:.2f}&nbsp;m/s'),
            ('Under the GP, koloa flags the clear outliers (90% or more).',
             ens['recall_clear'] >= 0.9,
             f'{100 * ens["recall_clear"]:.0f}% of the clear, '
             f'{100 * ens["recall_border"]:.0f}% of the borderline')],
        run='python demos/demo_gp_joint.py')


def uc_contamination():
    """How many outliers a fit can take"""
    summ = demo('contamination')
    if summ is None:
        return None
    rows = []
    for kind, label in (('visit', 'bad visits'), ('point', 'bad exposures')):
        fracs = summ[kind]['fracs']
        itop = len(fracs) - 1
        top = fracs[itop]
        mix = summ[kind]['koloa_mixture'][itop][0]
        gau = summ[kind]['gaussian'][itop][0]
        stu = summ[kind]['koloa_student'][itop][0]
        rows.append((kind, label, top, mix, gau, stu))
    crit = []
    for kind, label, top, mix, gau, stu in rows:
        crit.append((f'With {100 * top:.0f}% {label}, the error of K stays '
                     f'below 1&nbsp;m/s with the mixture.', mix < 1.0,
                     f'mixture {mix:.2f}&nbsp;m/s; Gaussian {gau:.1f}, '
                     f'Student-t {stu:.1f}'))
    return dict(
        question='How large a fraction of outliers can a fit take? The '
                 'mixture against a <a href="https://en.wikipedia.org/wiki/'
                 'Student%27s_t-distribution" target="_blank" '
                 'rel="noopener">Student-t</a> and least squares.',
        data='Simulated: an eccentric orbit with a growing fraction of bad '
             'exposures, or of bad visits, and a cubic with outliers of both '
             'signs (the test of polyband).',
        criteria=crit,
        run='python demos/demo_contamination.py')


def uc_completeness():
    """Recovery maps"""
    summ = demo('completeness')
    if summ is None:
        return None
    kol, gau, cln = summ['koloa'], summ['gaussian'], summ['koloa_clean']
    low, high = binomial(kol['nnull'], kol['false_alarm_nominal'])
    per = np.array(kol['periods'])
    it = int(np.argmin(np.abs(np.log(per / 10.0))))
    return dict(
        question='Which planets could a series have found, and is the '
                 'false-alarm rate what it claims?',
        data=f'Simulated: {kol["ninj"]} circular planets injected over a grid '
             f'of periods and semi-amplitudes into series on '
             f'{template_sampling()} with outliers, and {kol["nnull"]} series '
             f'without a planet; each searched by folding at the injected '
             f'period.',
        criteria=[
            ('The false alarms at K&nbsp;=&nbsp;0 are at the nominal rate '
             f'({100 * kol["false_alarm_nominal"]:.0f}%, {BINOMIAL} 2&sigma;).',
             low <= kol['false_alarm'] * kol['nnull'] <= high,
             f'{100 * kol["false_alarm"]:.1f}% of {kol["nnull"]}'),
            ('With outliers, koloa needs a smaller K than the Gaussian fold, '
             'at every period.',
             all(kk < gg for kk, gg in zip(kol['k50'], gau['k50'])),
             f'K found half the time at {per[it]:.0f}&nbsp;d: '
             f'{kol["k50"][it]:.1f}&nbsp;m/s against {gau["k50"][it]:.1f} '
             f'(without the outliers {cln["k50"][it]:.1f})')],
        run='python demos/demo_completeness.py')


def uc_campaign():
    """The injection-recovery campaign"""
    summ = load('paper', 'campaign_summary.json')
    if summ is None:
        return None
    stats, counts = summ['stats'], summ['counts']
    fpk = stats['visits_null']['koloa']['fp'][0]
    fpg = stats['visits_null']['gauss']['fp'][0]
    nnull = stats['visits_null']['koloa']['n']
    cal = summ['calibration']['koloa'][2]
    calg = summ['calibration']['gauss'][2]
    detk = stats['clean_planet']['koloa']['det'][0]
    detg = stats['clean_planet']['gauss']['det'][0]
    return dict(
        question='Over many series, how often does each method report a '
                 'false planet or find a real one, and do its probabilities '
                 'mean what they say?',
        data=f'Simulated: {sum(counts.values())} series on '
             f'{template_sampling()}, with its error bars and a visit jitter: no '
             f'planet or a K&nbsp;=&nbsp;3&nbsp;m/s planet, each without '
             f'outliers, with bad visits, with spiked exposures and with '
             f'borderline outliers; the FIP computed seven ways on each.',
        criteria=[
            ('With bad visits and no planet, koloa reports a planet '
             f'(FIP&nbsp;&lt;&nbsp;1%) in at most 1 of {nnull} series.',
             fpk * nnull <= 1,
             f'koloa {100 * fpk:.0f}%, the fixed-jitter FIP '
             f'{100 * fpg:.0f}%'),
            ('Its peaks with a TIP of 0.90 to 0.99 are real at least 90% of '
             'the time.', cal['true'] >= 0.9 * cal['n'],
             f'{cal["true"]} of {cal["n"]}; the fixed-jitter FIP '
             f'{calg["true"]} of {calg["n"]}'),
            ('On clean series, it detects the planet as often as the '
             'fixed-jitter FIP at FIP&nbsp;&lt;&nbsp;1%.', detk >= detg,
             f'koloa {100 * detk:.0f}%, the fixed-jitter FIP '
             f'{100 * detg:.0f}%: the visit jitter makes koloa conservative')],
        run='python paper/campaign.py')


def uc_acceleration():
    """The acceleration of a star, under outliers"""
    summ = demo('secular')
    if summ is None:
        return None
    ens = summ['ensemble']
    kst, gst = ens['stats']['koloa'], ens['stats']['gaussian']
    cst = ens['stats']['gaussian_clean']
    nval = kst['n']
    low, high = binomial(nval, 0.95)
    return dict(
        question='Can the acceleration of a star (the pull of a companion) be '
                 'measured under outliers?',
        data=f'Simulated: {nval} series on {template_sampling()}, each with a '
             f'companion that accelerates the star by '
             f'{ens["accel"]:.0f}&nbsp;m/s/yr, a planet with '
             f'K&nbsp;=&nbsp;{ens["planet"]["K"]:.0f}&nbsp;m/s, a visit '
             f'jitter and outliers; the acceleration fitted with the planet.',
        criteria=[
            ('koloa\'s 95% intervals hold the true acceleration as often as '
             f'they should ({max(low, 0):.0f} to {min(high, nval):.0f} of '
             f'{nval}).', low <= kst['cover95'] <= high,
             f'{kst["cover95"]} of {nval}'),
            ('Its intervals are nearly as narrow as without the outliers '
             '(within 50%).',
             kst['median_sigma'] <= 1.5 * cst['median_sigma'],
             f'median half-width {kst["median_sigma"]:.2f}&nbsp;m/s/yr, '
             f'{cst["median_sigma"]:.2f} without the outliers, '
             f'{gst["median_sigma"]:.2f} for the Gaussian likelihood with '
             f'them')],
        run='python demos/demo_secular.py')


# =============================================================================
# The use cases: real stars
# =============================================================================
def uc_toi2120():
    """TOI-2120 b under added bad visits"""
    clean = demo('toi2120', 'clean_series.json')
    path = os.path.join(DEMOS, 'toi2120', 'results.npy')
    if clean is None or not os.path.exists(path):
        return None
    res = [rr for rr in np.load(path, allow_pickle=True) if 'scan' in rr]
    sub = [rr for rr in res if rr['scan'] == 'fraction'
           and np.isclose(rr['level'], 0.2)]
    ksc = float(np.std([rr['koloa'][0] for rr in sub]))
    gsc = float(np.std([rr['gaussian'][0] for rr in sub]))
    kk, kg = clean['koloa']['K'], clean['gaussian']['K']
    data = 'Real: the SPIRou velocities of TOI-2120 (PCA2D, LBL)'
    rdb = os.path.join(ROOT, 'data', 'toi2120_spirou_pca2d_0-7_bias10s3.rdb')
    if os.path.exists(rdb):
        from koloa.data import RVData
        ser = RVData.from_csv(rdb, name='TOI-2120')
        data += (f', {ser.n} exposures in {ser.nseq} visits over '
                 f'{np.ptp(ser.time):.0f}&nbsp;d')
    data += (', its transiting planet with the ephemeris of the transits as '
             'priors and a Mat&eacute;rn-3/2 GP of the activity; bad visits '
             'added to the real series, growing in amplitude, then in number.')
    return dict(
        question='Does the K of a known transiting planet survive more and '
                 'more bad visits?',
        data=data,
        criteria=[
            ('On the series as observed, koloa and the Gaussian likelihood '
             'agree on K (within their uncertainties).',
             abs(kk[0] - kg[0]) < np.hypot(half(kk), half(kg)),
             f'koloa {pm(kk, 1)}, Gaussian {pm(kg, 1)}&nbsp;m/s'),
            ('With 20% of the visits made bad, koloa\'s K scatters from one '
             'realisation to the next by no more than 1.5 times its '
             'uncertainty.', ksc <= 1.5 * half(kk),
             f'{ksc:.1f}&nbsp;m/s against {half(kk):.1f}; the Gaussian '
             f'{gsc:.1f}')],
        run='python demos/demo_toi2120.py')


#: a signal is a known planet when their periods are within this fraction
MATCH = 0.05


def match_known(orbits, known):
    """each fitted orbit with the known planet of the nearest period within
    MATCH, or None"""
    out = []
    for orb in orbits:
        best = min(known, key=lambda pl: abs(pl['P'] / orb['P'][0] - 1),
                   default=None)
        if best is not None and abs(best['P'] / orb['P'][0] - 1) >= MATCH:
            best = None
        out.append((orb, best))
    return out


def comparisons(orb, pl):
    """
    The published solutions of a known planet against a fitted orbit, the
    most recent first, each with its distance to the fit in sigma (z): those
    the demo wrote, or, from a run before it kept every solution, the
    archive's default alone

    :return: list of dict (reference, reference_url, K, K_err, z)
    """
    if orb.get('comparisons'):
        return orb['comparisons']
    if not pl or pl.get('K') is None or not pl.get('K_err'):
        return []
    sols = pl.get('solutions') or [pl]
    return [dict(sol, z=float((orb['K'][0] - sol['K'])
                              / np.hypot(half(orb['K']), sol['K_err'])))
            for sol in reversed(sols)
            if sol.get('K') is not None and sol.get('K_err')]


def uc_dace():
    """A known planetary system from DACE, every instrument"""
    summ = demo('dace')
    if summ is None:
        return None
    fip = summ['fip_second']
    known = summ['archive']['planets']
    found = [pk for pk in fip['peaks'] if pk['fip'] < 0.01]
    crit = []
    for pl in known:
        hit = [pk for pk in found if abs(pk['period'] / pl['P'] - 1) < MATCH]
        crit.append((f'The FIP finds {escape(pl["name"])} ({pl["P"]:.2f}&nbsp;d): '
                     f'an interval with FIP&nbsp;&lt;&nbsp;1% within '
                     f'{100 * MATCH:.0f}% of its period.', bool(hit),
                     (f'{hit[0]["period"]:.3f}&nbsp;d, FIP {fipfmt(hit[0]["fip"])}'
                      if hit else 'no such interval')))
    extra = [pk for pk in found
             if not any(abs(pk['period'] / pl['P'] - 1) < MATCH for pl in known)]
    crit.append(('Another interval with FIP&nbsp;&lt;&nbsp;1%?', None,
                 ', '.join(f'{pk["period"]:.2f}&nbsp;d (FIP {fipfmt(pk["fip"])})'
                           for pk in extra) or 'none'))
    # K against the most recent published solution (the archive's default
    #   may be the discovery paper, from far fewer velocities), the others
    #   beside it
    for orb, pl in match_known(summ['orbits'], known):
        sols = comparisons(orb, pl)
        if not sols:
            continue
        last = sols[0]
        crit.append((f'The K of {escape(pl["name"])} agrees with its most '
                     f'recent published value, {escape(last["reference"])} '
                     f'(within 2&sigma;).', abs(last['z']) < 2,
                     f'{pm(orb["K"])}&nbsp;m/s against '
                     + '; '.join(('' if it == 0 else
                                  f'{escape(sol["reference"])}: ')
                                 + f'{sol["K"]:.2f}&nbsp;&plusmn;&nbsp;'
                                   f'{sol["K_err"]:.2f} '
                                   f'({abs(sol["z"]):.1f}&sigma;)'
                                 for it, sol in enumerate(sols))))
    insts = ', '.join(f'{row["inst"]} ({row["n"]})' for row in summ['instruments'])
    return dict(
        question='Does koloa find the known planets of a star from the public '
                 'velocities of every instrument that observed it, and measure '
                 'them as published?',
        data=f'Real: the public velocities of {escape(summ["host"])} on DACE, '
             f'fetched by koloa.dace: {summ["n"]} velocities over '
             f'{summ["baseline_years"]:.1f}&nbsp;yr from {insts}; each '
             f'instrument era with its own offset and noise; the NASA Exoplanet '
             f'Archive for the known planets.',
        criteria=crit,
        run='python demos/demo_dace.py --target HD69830 --host "HD 69830"')


def uc_rotation():
    """The rotation of a star from its temperature indicator, robust GP"""
    summ = demo('rotation')
    if summ is None:
        return None
    mix = summ['posteriors']['mixture']
    gau = summ['posteriors']['gaussian']
    per = mix['gp_period_peak']
    lit = summ['literature']
    rob = summ['robustness']
    need = rob['nreal'] - 1
    # the distance of the GP's period from the photometric one, in its own
    #   sigma (the side of the interval towards it)
    side = per[2] if lit['P'] > per[0] else per[1]
    zlit = abs(lit['P'] - per[0]) / side
    return dict(
        question='Can the rotation of a star be measured from its temperature '
                 'indicator, robustly, with a GP?',
        data=f'Real: the {summ["instrument"]} DTEMP3500 of '
             f'{escape(summ["star"])}, {summ["n"]} values in {summ["nvisits"]} '
             f'visits over {summ["baseline"] / 365.25:.1f}&nbsp;yr; an SHO GP '
             + (f'(its quality factor held above {summ["quality"][0]:.0f}, so '
                f'that its power has a peak) ' if 'quality' in summ else '')
             + f'fitted with koloa\'s likelihood and with a gaussian one, both '
               f'sampled by MCMC.',
        criteria=[
            (f'The outlier-aware periodogram of DTEMP3500 finds the rotation of '
             f'the photometry ({lit["reference"]}, {lit["P"]:.0f}&nbsp;d) '
             f'within 5%.', abs(summ['peak_mix'] / lit['P'] - 1) < 0.05,
             f'{summ["peak_mix"]:.2f}&nbsp;d'),
            ('Another indicator, the FWHM of the lines, peaks at the same '
             'period (within 5%).',
             abs(summ['peak_fwhm'] / summ['peak_mix'] - 1) < 0.05,
             f'{summ["peak_fwhm"]:.2f}&nbsp;d'),
            (f'The period at which the power of the SHO peaks agrees with the '
             f'photometry within 2&sigma;.', zlit < 2,
             f'{pm(per, 1)}&nbsp;d, {zlit:.1f}&sigma; from {lit["P"]:.0f}&nbsp;d '
             f'(the SHO has Q = {mix["gp_quality"][0]:.1f})'),
            (f'Robust: with {100 * rob["fraction"]:.0f}% of the visits made bad, '
             f'koloa\'s period stays within 1&sigma; of the clean one in at '
             f'least {need} of {rob["nreal"]} realisations.',
             rob['held']['mixture'] >= need,
             f'koloa {rob["held"]["mixture"]} of {rob["nreal"]}; the gaussian '
             f'likelihood {rob["held"]["gaussian"]}'),
            ('How many DTEMP values does koloa flag?', None,
             f'{summ["flagged"]} of {summ["n"]}; the gaussian posterior of the '
             f'period {pm(gau["gp_period_peak"], 1)}&nbsp;d')],
        run='python demos/demo_rotation.py --nsteps 16000 --nreal 6 --ncpu 9')


def uc_kepler21():
    """Kepler-21 b under activity"""
    summ = demo('kepler21')
    if summ is None:
        return None
    ser, kol = summ['summary'], summ['koloa']
    lit = summ['literature']
    name = max(lit, key=lambda key: int(key.split()[-1]))
    kpub, epub = lit[name]
    others = '; '.join(f'{key}: {val[0]:.2f}&nbsp;&plusmn;&nbsp;{val[1]:.2f}'
                       for key, val in lit.items() if key != name)
    return dict(
        question='Can koloa measure a small transiting planet under stellar '
                 'activity?',
        data=f'Real: the public HARPS-N velocities of Kepler-21 (DACE), '
             f'{series_text(ser)}, scattering by {ser["rms"]:.1f}&nbsp;m/s '
             f'against errors of {ser["median_err"]:.1f}&nbsp;m/s; a transiting '
             f'super-Earth, the ephemeris from the transits, and a GP whose '
             f'rotation period comes from the S index.',
        criteria=[
            ('Is the planet found blind, without its ephemeris?', None,
             f'no: koloa FIP at its period {fipfmt(summ["fip"]["koloa"]["at_b"])}'),
            (f'With the ephemeris and the GP, K agrees with {name} '
             '(within 1&sigma;).',
             abs(kol['K'][0] - kpub) < np.hypot(half(kol['K']), epub),
             f'{pm(kol["K"])}&nbsp;m/s against {kpub:.2f}&nbsp;&plusmn;&nbsp;'
             f'{epub:.2f} ({others})')],
        run='python demos/demo_kepler21.py')


def uc_gl725():
    """GL 725 A and B pull on each other"""
    summ = demo('gl725b')
    sec = demo('secular')
    if summ is None or sec is None:
        return None
    sb, sa = summ['GJ 725 B'], summ['GJ 725 A']
    acb, aca = sb['linear']['orbital'], sa['linear']['orbital']
    ratio = summ['mass_ratio']
    gl = sec['gl725']
    pred = gl['GJ 725 A']['mass_other'][0] / gl['GJ 725 B']['mass_other'][0]
    sig = 0.5 * (ratio[2] - ratio[0])
    zval = (ratio[1] - pred) / sig
    return dict(
        question='Do the two stars of a binary pull on each other as their '
                 'masses say?',
        data=f'Real: the SPIRou velocities of GL&nbsp;725&nbsp;B '
             f'({series_text(sb)}) and GL&nbsp;725&nbsp;A ({series_text(sa)}), '
             f'{gl["separation"]:.1f}&nbsp;arcsec apart; APERO has already '
             f'removed their perspective accelerations.',
        criteria=[
            ('The two stars accelerate in opposite directions.',
             np.sign(acb[1]) != np.sign(aca[1]),
             f'B {acb[1]:+.2f}, A {aca[1]:+.2f}&nbsp;m/s/yr'),
            ('The ratio of the accelerations is that of the masses '
             '(within 2&sigma;).', abs(zval) < 2,
             f'&minus;a<sub>A</sub>/a<sub>B</sub>&nbsp;=&nbsp;{ratio[1]:.2f}'
             f'&nbsp;&plusmn;&nbsp;{sig:.2f} against M<sub>B</sub>/M<sub>A</sub>'
             f'&nbsp;=&nbsp;{pred:.3f} ({abs(zval):.1f}&sigma;)')],
        run='python demos/demo_gl725b.py; python demos/demo_secular.py')


def uc_mdwarfs():
    """Three M dwarfs cleaned by PCA2D"""
    summ = demo('mdwarfs')
    if summ is None:
        return None
    parts, planets, scat = [], [], []
    for key in ('gl687', 'proxima', 'gl699'):
        if key not in summ:
            continue
        ent = summ[key]
        ser, deliv = ent['series'], ent.get('delivered', {})
        parts.append(f'{escape(ent["name"])} ({ent["inst"]}, '
                     f'{ser["npoints"]} exposures in {ser["nseq"]} visits)')
        for pname, pent in ent.get('planets', {}).items():
            lit = pent['literature']
            kk = pent['koloa']['K']
            ok = abs(kk[0] - lit['K']) < 2 * np.hypot(half(kk), lit['Kerr'])
            planets.append((f'{escape(ent["name"])} {pname}', ok,
                            f'{pm(kk)} against {lit["K"]:.2f}&nbsp;&plusmn;'
                            f'&nbsp;{lit["Kerr"]:.2f}'))
        night = ser.get('nightly_robust_line',
                        ser.get('nightly_rms_line', ser.get('nightly_rms')))
        dnight = deliv.get('nightly_robust_line',
                           deliv.get('nightly_rms_line',
                                     deliv.get('nightly_rms')))
        if night is not None and dnight is not None:
            scat.append((escape(ent['name']), dnight, night))
    crit = [('The known planets have their published K, within 2&sigma;.',
             all(ok for _, ok, _ in planets),
             '; '.join(f'{name}: {txt}&nbsp;m/s' for name, _, txt in planets))]
    if scat:
        crit.append(('PCA2D lowers the scatter of the visits of every star.',
                     all(after < before for _, before, after in scat),
                     '; '.join(f'{name} {before:.2f} &rarr; {after:.2f}'
                               for name, before, after in scat)
                     + '&nbsp;m/s'))
    return dict(
        question='Does cleaning the spectra with PCA2D help, and does koloa '
                 'find the known planets of M dwarfs with many visits?',
        data='Real: ' + ', '.join(parts) + ', their telluric-corrected '
             'spectra cleaned by pca2d-preclean and measured by LBL, beside '
             'the spectra as delivered.',
        criteria=crit,
        run='python demos/demo_mdwarfs.py')


#: the use cases, in the order of the tabs: key, tab, title, kind, function
USECASES = [
    ('nofalse', 'No planet', 'A planet that is not there', 'simulated',
     uc_nofalse),
    ('recovery', 'A planet behind outliers',
     'A planet that is there, behind outliers clear and borderline',
     'simulated', uc_recovery),
    ('posterior', 'Orbit posteriors',
     'What outliers do to an orbit: posteriors by MCMC', 'simulated',
     uc_posterior),
    ('twoplanets', 'Two planets', 'Two planets, one FIP periodogram',
     'simulated', uc_twoplanets),
    ('gpjoint', 'Activity and a planet', 'Activity and a planet, fitted jointly',
     'simulated', uc_gpjoint),
    ('contamination', 'How many outliers',
     'How many outliers can a fit take? A Student-t against the mixture',
     'simulated', uc_contamination),
    ('completeness', 'Recovery maps',
     'What a series could have found: recovery rates in period and K',
     'simulated', uc_completeness),
    ('campaign', 'The campaign', 'How often: the injection-recovery campaign',
     'simulated', uc_campaign),
    ('acceleration', 'A star\'s acceleration',
     'The acceleration of a star under outliers', 'simulated',
     uc_acceleration),
    ('toi2120', 'TOI-2120 b', 'TOI-2120 b: the same planet under more and more '
     'bad visits', 'real', uc_toi2120),
    ('dace', 'HD 69830 from DACE', 'HD 69830 from DACE: three known planets, '
     'every instrument, fetched and fitted', 'real', uc_dace),
    ('kepler21', 'Kepler-21 b', 'Kepler-21 b with HARPS-N: a small transiting '
     'planet under activity', 'real', uc_kepler21),
    ('gl725', 'GL 725 A and B', 'GL 725 A and B with SPIRou: each star pulled '
     'by the other', 'real', uc_gl725),
    ('rotation', 'Rotation from DTEMP', 'GJ 687: the rotation from its '
     'temperature indicator, with a robust GP', 'real', uc_rotation),
    ('mdwarfs', 'M dwarfs, PCA2D', 'Three M dwarfs with more than 150 visits, '
     'cleaned by PCA2D', 'real', uc_mdwarfs),
]


# =============================================================================
# The HTML
# =============================================================================
def badge(ok):
    """PASS, FAIL or INFO, with a symbol: never colour alone"""
    if ok is None:
        return '<span class="uc-badge uc-info">&#9432; INFO</span>'
    if ok:
        return '<span class="uc-badge uc-pass">&#10003; PASS</span>'
    return '<span class="uc-badge uc-fail">&#10007; FAIL</span>'


def spec_html(key, kind, case):
    """The test of one use case: the question, the data, the criteria with
    their verdicts, and the command"""
    rows = ''.join(
        f'<li>{badge(None if ok is None else bool(ok))}'
        f'<span class="uc-crit">{text}</span>'
        f'<span class="uc-measured">{measured}</span></li>'
        for text, ok, measured in case['criteria'])
    npass = sum(1 for _, ok, _ in case['criteria']
                if ok is not None and bool(ok))
    ntest = sum(1 for _, ok, _ in case['criteria'] if ok is not None)
    return (f'<div class="uc-spec">'
            f'<div class="uc-row"><span class="uc-key">question</span>'
            f'<span>{case["question"]}</span></div>'
            f'<div class="uc-row"><span class="uc-key">data</span>'
            f'<span><span class="uc-kind uc-{kind}">{kind}</span> '
            f'{case["data"].split(": ", 1)[-1]}</span></div>'
            f'<div class="uc-row"><span class="uc-key">run</span>'
            f'<code>{escape(case["run"])}</code></div>'
            f'<div class="uc-row"><span class="uc-key">criteria</span>'
            f'<span class="uc-score">{npass} of {ntest} pass</span></div>'
            f'<ul class="uc-criteria">{rows}</ul></div>')


def build():
    """
    Every use case that has its outputs

    :return: tuple, dict key: the HTML of its test (question, data,
             criteria, command); the tally of the criteria; dict key: (passed,
             tested, question)
    """
    out, tally, scores = {}, [0, 0], {}
    for key, _, _, kind, func in USECASES:
        case = func()
        if case is None:
            continue
        out[key] = spec_html(key, kind, case)
        tested = [bool(ok) for _, ok, _ in case['criteria'] if ok is not None]
        tally[0] += sum(tested)
        tally[1] += len(tested)
        scores[key] = (sum(tested), len(tested), case['question'])
    return out, tally, scores


if __name__ == '__main__':
    specs, tally, _ = build()
    for key, html in specs.items():
        import re
        text = re.sub('<[^>]+>', ' ', html)
        print(f'== {key}: {re.sub(" +", " ", text)[:600]}')
    print(f'{tally[0]} of {tally[1]} criteria pass')

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa.prospects: what more velocities would do for a candidate, and its
astrometric signal in Gaia DR4.

Created on 2026-10-08

@author: artigau
"""
import numpy as np
import pytest

from koloa import log as klog
from koloa import prospects

klog.VERBOSE = False


def test_the_nights_a_star_can_be_observed():
    # GJ 876 (dec -14) from Maunakea: about 230 nights a year above
    #   airmass 2 at night (as astropy counts them), half of them in the
    #   bright half of the lunations, none of a star that never rises
    start = 61300.0
    year = prospects.nights(343.3197, -14.2637, 'CFHT', start, 365.25)
    assert 225 <= len(year) <= 235
    assert np.all(np.diff(year) > 0.5) and year[0] >= start
    bright = prospects.nights(343.3197, -14.2637, 'CFHT', start, 365.25,
                              bright=True)
    assert 0.45 < len(bright) / len(year) < 0.6
    # within a quarter of a lunation of a full Moon
    off = ((bright + 2400000.0 - prospects.FULL_MOON + prospects.SYNODIC / 2)
           % prospects.SYNODIC) - prospects.SYNODIC / 2
    assert np.max(np.abs(off)) < prospects.SYNODIC / 4 + 1.0
    # GJ 15 A (dec +44) never reaches airmass 2 from La Silla
    assert len(prospects.nights(4.595, 44.023, 'La Silla', start, 365.25)) == 0
    # six months hold about half of the nights of the year, or fewer
    assert len(prospects.nights(343.3197, -14.2637, 'CFHT', start, 183.0)) \
        < 0.6 * len(year)


def _fold(seed=3, period=23.7, amp=3.0):
    """the fold of a signal in two instruments, as a quick look keeps it"""
    rng = np.random.default_rng(seed)
    out = []
    for name, t0, num, err in (('SPIRou', 59700.0, 45, 2.4),
                               ('NIRPS', 60100.0, 35, 1.8)):
        time = np.sort(t0 + rng.uniform(0, 500, num))
        out.append(dict(name=name, time=list(time), err=[err] * num,
                        rv=[0.0] * num, valid=[1.0] * num))
    return dict(id=1, period=period, K=amp, K_err=0.5, P_err=0.05,
                tc=59700.0, instruments=out)


def test_what_more_nights_would_do():
    # the error of K of a sinusoid sampled at random: sigma sqrt(2 / N)
    rng = np.random.default_rng(1)
    time = np.sort(rng.uniform(0, 400, 800))
    kerr, perr = prospects._errors(time, np.full(800, 2.0), np.full(
        800, 'A'), 17.3, 5.0, 3.0)
    assert kerr == pytest.approx(2.0 * np.sqrt(2 / 800), rel=0.05)
    assert 0 < perr < 0.01
    fold = _fold()
    made = prospects.more_velocities(fold, 'NIRPS', 'La Silla', 258.829,
                                     4.964, dict(NIRPS=0.6), fip=1e-3,
                                     start=61300.0, seed=4)
    assert made['n'] == 35 and not made['bright']
    assert made['error'] == pytest.approx(1.8)
    assert made['sigma'] == pytest.approx(np.hypot(1.8, 0.6))
    assert made['snr'] == pytest.approx(6.0)
    assert len(made['table']) == 12
    year = {row['n']: row for row in made['table'] if row['span'] == '1 year'}
    # more nights: a smaller error of K, a larger K over it, a smaller FIP
    snr = [year[num]['snr'] for num in prospects.MORE]
    assert snr == sorted(snr) and snr[0] > made['snr']
    assert year[200]['K_err'] < year[50]['K_err'] < fold['K_err']
    assert year[200]['log10_fip'] < year[50]['log10_fip'] < -3
    assert year[100]['log10_fip'] == pytest.approx(
        -3 - (year[100]['snr'] ** 2 - 36.0) / (2 * np.log(10)))
    # a longer span: a better period
    two = next(row for row in made['table'] if row['n'] == 100
               and row['span'] == '2 years')
    half = next(row for row in made['table'] if row['n'] == 100
                and row['span'] == '6 months')
    assert two['P_err'] < half['P_err'] < fold['P_err']
    # SPIRou has the bright half of the nights: more visits asked than
    #   nights in six months, several a night, and said
    spirou = prospects.more_velocities(fold, 'SPIRou', 'CFHT', 258.829, 4.964,
                                       start=61300.0)
    assert spirou['bright']
    row = next(row for row in spirou['table'] if row['n'] == 200
               and row['span'] == '6 months')
    assert row['stacked'] and row['nights'] < 200 and row['snr'] > 6.0
    # a star the site never sees: no forecast in any window
    never = prospects.more_velocities(fold, 'NIRPS', 'La Silla', 4.595,
                                      44.023, start=61300.0)
    assert all(row['snr'] is None and row['nights'] == 0
               for row in never['table'])
    assert prospects.more_velocities(fold, 'HARPS', 'La Silla', 0, 0) is None


def test_the_astrometric_signal_of_a_planet():
    # Jupiter moves the Sun by 5.2 au / 1047: 497 uas from 10 pc
    jup = prospects.astrometry(317.83, 4332.6, 1.0, 100.0, 5.0)
    assert jup['a'] == pytest.approx(5.2, abs=0.01)
    assert jup['alpha'] == pytest.approx(497.0, rel=0.01)
    assert 'longer than' in jup['verdict']          # 11.9 yr: not in DR4
    assert 'brighter than G = 6' in jup['verdict']
    # GJ 876 b, at least: 274 uas (HST measured 250 +- 60, Benedict et al.
    #   2002), five times the noise of a crossing
    b = prospects.astrometry(723.0, 61.12, 0.37, 214.04, 8.88, 23, 1.34)
    assert b['alpha'] == pytest.approx(274.0, rel=0.02)
    assert b['sigma_fov'] == 54.0 and b['snr'] == pytest.approx(5.07, rel=0.02)
    assert b['transits'] == pytest.approx(23 * 66 / 34) and b['counted']
    # the chi2 of Lammers & Winn 2025: 22 snr^2 at 66 crossings, 50 at 1.5
    assert b['dchi2'] == pytest.approx(b['transits'] * b['snr'] ** 2 / 3)
    edge = prospects.astrometry(1.5 * 54.0 / 95.4 * 317.83, 365.25, 1.0,
                                100.0, 10.0, 34)
    assert edge['dchi2'] == pytest.approx(49.5, rel=0.01)
    assert edge['dchi2_limit'] == 50.0
    assert b['verdict'].startswith('detectable in DR4')
    # a small planet close to an M dwarf: far below the noise
    small = prospects.astrometry(5.0, 23.7, 0.18, 68.3, 13.0)
    assert small['snr'] < 0.02 and not small['counted']
    assert small['transits'] == prospects.DR4_TRANSITS
    assert 'far too small' in small['verdict']
    # the noise of a crossing: flat to G = 14, the photons beyond
    assert prospects.sigma_fov(9.0) == prospects.sigma_fov(14.0) == 54.0
    assert prospects.sigma_fov(16.5) == pytest.approx(54.0 * 10 ** 0.5)
    # between the thresholds of DR5, of DR4 and the older one
    for alpha, word in ((0.6, 'not detectable in DR4'), (1.2, 'DR5'),
                        (2.0, 'marginal'), (4.0, 'detectable in DR4')):
        told = prospects.astrometry(alpha * 54.0 / 95.4 * 317.83, 365.25,
                                    1.0, 100.0, 10.0)
        assert told['snr'] == pytest.approx(alpha) and word in told['verdict']


def test_the_prospects_of_a_candidate():
    fold = _fold()
    res = dict(name='GJ 1214', status='done', summary=dict(
        K=3.0, K_err=0.5, sources=dict(SPIRou='file: lbl.rdb',
                                       NIRPS='file: lbl.rdb',
                                       HARPS03='DACE')),
        reading=dict(kind='candidate', mass=[0.18, 0.02], peaks=[dict(
            id=1, period=23.7, fip=1e-4, K=3.0, K_err=0.5, counts=True,
            what='', notes=[], msini=[5.0, 1.0, 1.0])]))
    quick = dict(result=dict(folds=[fold], inflation=dict(SPIRou=0.0),
                             settings=dict(trend=1)))
    star = dict(ra=258.829, dec=4.964, plx=68.3, G=13.0, gaia=dict(
        G=13.0, parallax=68.3, ruwe=1.01, transits=52.0))
    told = prospects.candidate(res, quick, star, start=61300.0)
    # the instruments of the files, not those of the archives
    assert [one['instrument'] for one in told['velocities']] == [
        'SPIRou', 'NIRPS']
    assert [one['site'] for one in told['velocities']] == ['CFHT', 'La Silla']
    assert told['astrometry']['alpha'] == pytest.approx(0.51, rel=0.05)
    assert told['astrometry']['transits'] == pytest.approx(52 * 66 / 34)
    assert told['notes'] == []
    words = prospects.lines(told)
    assert words[0].startswith('more SPIRou from CFHT, in bright time: '
                               '2.40 m/s a night')
    assert any(line.startswith('Gaia DR4: the star moves by at least')
               for line in words)
    assert any('a detection asks for 50' in line for line in words)
    assert words[-1].strip().startswith('Not detectable by Gaia')
    assert max(len(line) for line in words) <= 124
    assert 'Gaia DR4' in prospects.short(told)
    # a star with archives only: no instrument to ask for more nights of;
    #   no parallax: no astrometric signal; no candidate: nothing
    alone = dict(res, summary=dict(res['summary'],
                                   sources=dict(HARPS03='DACE')))
    told = prospects.candidate(alone, quick, dict(ra=1.0, dec=2.0))
    assert told['velocities'] == [] and told['astrometry'] is None
    assert len(told['notes']) == 2
    known = dict(res, reading=dict(res['reading'], kind='known'))
    assert prospects.candidate(known, quick, star) is None
    assert prospects.lines(None) == [] and prospects.short(None) == ''


def test_what_gaia_has_of_the_stars(monkeypatch):
    from koloa import gather, survey
    asked = []

    def tap(url, query, timeout=120.0):
        asked.append((url, query))
        return [dict(source_id='4393265392167891712',
                     phot_g_mean_mag='12.9968', parallax='68.2986',
                     ruwe='1.0087', astrometric_matched_transits='52',
                     visibility_periods_used='21')]
    monkeypatch.setattr(gather, '_tap', tap)
    stars = [dict(name='GJ 1214', ids={'Gaia DR3': 'Gaia DR3 '
                                       '4393265392167891712'}),
             dict(name='No number', ids={})]
    assert survey.gaia_dr3(stars) == 1
    assert asked[0][0] == survey.GAIA_TAP
    assert '4393265392167891712' in asked[0][1]
    assert stars[0]['gaia'] == dict(G=12.9968, parallax=68.2986, ruwe=1.0087,
                                    transits=52.0, periods=21.0)
    assert 'gaia' not in stars[1]
    # asked once: what a star already has is kept
    assert survey.gaia_dr3(stars) == 0 and len(asked) == 1

# =============================================================================
# End of code
# =============================================================================

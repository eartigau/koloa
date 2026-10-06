#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa's transit search in TESS (koloa.transit), on simulated light curves:
a transit injected found, about its ephemeris and over the whole phase;
none where there is none; one event alone not taken for a transit; the
transits of another planet not lined up at a wrong period

Created on 2026-10-03

@author: artigau
"""
import numpy as np
import pytest

from koloa import transit


def _light_curve(seed=1, transits=(), sectors=((58000.0, 27.0),
                                               (58400.0, 27.0)),
                 spots=2.0, white=0.6):
    """two sectors of 2-minute TESS-like flux [ppt]: white noise, a spotted
    star's modulation (6.3 d), and the transits given (period, t0, depth
    [ppt], duration [h])"""
    rng = np.random.default_rng(seed)
    times, sec = [], []
    for k, (start, length) in enumerate(sectors):
        tt = np.arange(start, start + length, 2.0 / 1440)
        # the gap of the orbit of TESS in the middle of a sector
        tt = tt[np.abs(tt - start - length / 2) > 0.5]
        times.append(tt)
        sec.append(np.full(len(tt), float(k + 1)))
    time = np.concatenate(times)
    flux = spots * np.sin(2 * np.pi * time / 6.3) + rng.normal(0, white,
                                                               len(time))
    for period, t0, depth, dur in transits:
        phase = ((time - t0) / period + 0.5) % 1.0 - 0.5
        flux[np.abs(phase * period * 24.0) < dur / 2] -= depth
    return dict(time=time, flux=flux, err=np.full(len(time), white),
                sector=np.concatenate(sec))


def test_the_duration_and_the_depths():
    """a central transit of the Earth before the Sun (13 h), and the depth
    of the Earth and of Jupiter before the Sun"""
    assert transit.duration(365.25, 1.0, 1.0) == pytest.approx(13.0, abs=0.2)
    assert transit.depth_of(1.0, 1.0) == pytest.approx(0.0840, abs=2e-4)
    assert transit.depth_of(transit.R_JUPITER / transit.R_EARTH, 1.0) == \
        pytest.approx(10.56, abs=0.02)
    assert transit.radius_of(transit.depth_of(2.0, 0.4), 0.4) == \
        pytest.approx(2.0)


def test_a_transit_found():
    """a 1.2 ppt transit every 3.1 d: found about its ephemeris (even one
    off by an hour, with its error) and over the whole phase, plausible,
    at its depth"""
    lc = _light_curve(transits=[(3.1, 58001.3, 1.2, 1.6)])
    for t0, t0_err in ((58001.3, 0.01), (58001.34, 0.05), (None, None)):
        res = transit.search(lc, 3.1, t0, t0_err, 1e-5, mstar=0.4,
                             rstar=0.4, ntrial=10)
        best = res['best']
        assert res['plausible'], res['why']
        # which sentence the verdict is, for a page in another language
        assert res['why_code'] == 'plausible' and res['snr_need'] > 0
        assert best['depth'] == pytest.approx(1.2, abs=0.25)
        # the transit at its time
        off = ((best['centre'] - 58001.3) / 3.1 + 0.5) % 1.0 - 0.5
        assert abs(off * 3.1 * 24) < 1.0
        assert best['ntransits'] >= 12
    assert res['window'] is None and len(res['bins']) > 10
    assert res['depth_jupiter'] > res['depth_earth'] > 0


def test_no_transit_where_there_is_none():
    """the same light curve without a transit: nothing plausible, about an
    ephemeris or over the whole phase"""
    lc = _light_curve(seed=2)
    for t0 in (58001.3, None):
        res = transit.search(lc, 3.1, t0, 0.05, 1e-4, mstar=0.4, rstar=0.4,
                             ntrial=10)
        assert not res['plausible'], res['why']
        assert res['why_code'] in ('none', 'weak', 'single', 'one', 'dips',
                                   'chance')


def test_one_event_is_not_a_transit():
    """one deep dip (a systematic of one orbit): not a transit seen twice"""
    lc = _light_curve(seed=3)
    lc['flux'][np.abs(lc['time'] - 58410.0) < 0.05] -= 4.0
    res = transit.search(lc, 5.0, None, mstar=0.4, rstar=0.4, ntrial=10)
    assert not res['plausible']
    assert 'one transit' in res['why'] or 'single' in res['why']
    assert res['why_code'] in ('one', 'single')


def test_another_planet_left_out():
    """the transits of a known planet (2.2 d, deep) do not line up at a
    wrong period once left out"""
    lc = _light_curve(seed=4, transits=[(2.2, 58000.7, 3.0, 1.5)])
    other = [dict(P=2.2, tc=58000.7, duration=1.5)]
    res = transit.search(lc, 2.9, None, mstar=0.4, rstar=0.4, others=other,
                         ntrial=10)
    assert not res['plausible'], res['why']


def test_a_deep_transit_of_a_faint_m_dwarf():
    """GJ 1214 b as TESS would see it (TESS has not: no sector holds the
    star): 13.5 ppt, 52 minutes, every 1.58 d, on a faint M dwarf (3.5 ppt
    a point, a slow rotation). Found at its ephemeris, about a conjunction
    as the velocities give one (its period then scanned across sectors two
    years apart) and over the whole phase, at its depth and radius"""
    period, tzero, depth, dur = 1.580404531, 59001.1, 13.5, 0.87
    lc = _light_curve(seed=3, sectors=((59000.0, 27.0), (59750.0, 27.0)),
                      spots=0.0, white=3.5,
                      transits=[(period, tzero, depth, dur)])
    lc['flux'] = lc['flux'] + 3.0 * np.sin(2 * np.pi * lc['time'] / 40.0)
    asked = ((tzero, 0.002, 1e-7),
             (tzero + 0.04, transit.conjunction_error(0.02, period),
              period ** 2 / (4 * 3000.0)),
             (None, None, None))
    for t0, t0_err, perr in asked:
        res = transit.search(lc, period, t0, t0_err, perr, mstar=0.18,
                             rstar=0.215, ntrial=10)
        assert res['plausible'], res['why']
        assert res['best']['snr'] > 40 and res['best']['ntransits'] > 25
        fit = res['fit']
        assert fit['depth'] == pytest.approx(depth, abs=1.0)
        assert fit['duration'] == pytest.approx(dur, abs=0.25)
        # 2.7 Earth radii before a 0.215 Rsun star
        assert fit['radius'] == pytest.approx(
            transit.radius_of(depth, 0.215), rel=0.06)
        assert res['fold_period'] == pytest.approx(period, abs=2e-4)
    # no deeper than Jupiter's: nothing rules it out as a planet
    assert res['depth_jupiter'] > depth


def test_a_deep_transit_keeps_its_points():
    """a transit tens of sigma deep a point (a Jupiter before an M dwarf,
    117 ppt; a hot Jupiter before a Sun, 12 ppt at 0.6 ppt a point): only
    the points above the curve (flares) are cut, never its dips; found at
    its depth, a planet of Jupiter's size"""
    for depth, dur, white, mstar, rstar in ((117.0, 1.2, 2.0, 0.3, 0.3),
                                            (12.0, 2.6, 0.6, 1.0, 1.0)):
        lc = _light_curve(seed=5, spots=2.0, white=white,
                          transits=[(3.2, 58001.1, depth, dur)])
        inside = np.abs((((lc['time'] - 58001.1) / 3.2 + 0.5) % 1.0 - 0.5)
                        * 3.2 * 24.0) < dur / 2
        trimmed = transit.trim_edges(lc)
        kept = transit.highpass(trimmed, 3 * dur / 24.0)
        left = np.abs((((kept['time'] - 58001.1) / 3.2 + 0.5) % 1.0 - 0.5)
                      * 3.2 * 24.0) < dur / 2
        asked = np.abs((((trimmed['time'] - 58001.1) / 3.2 + 0.5) % 1.0
                        - 0.5) * 3.2 * 24.0) < dur / 2
        # every point of the transits is still there
        assert inside.sum() > 300 and left.sum() == asked.sum()
        for t0 in (58001.1, None):
            res = transit.search(lc, 3.2, t0, 0.01 if t0 else None,
                                 1e-6 if t0 else None, mstar=mstar,
                                 rstar=rstar, ntrial=5)
            assert res['plausible'], res['why']
            assert res['fit']['depth'] == pytest.approx(depth, rel=0.02)
            assert res['fit']['radius'] == pytest.approx(
                transit.radius_of(depth, rstar), rel=0.02)


def test_the_light_curve_of_the_archives(tmp_path):
    """phot/tess.csv as koloa.gather writes it"""
    lc = _light_curve(seed=5)
    path = tmp_path / 'tess.csv'
    path.write_text('rjd,flux,sflux,sector,pipeline\n' + ''.join(
        f'{tt},{ff},{ee},{int(ss)},SPOC\n' for tt, ff, ee, ss in zip(
            lc['time'][:500], lc['flux'][:500], lc['err'][:500],
            lc['sector'][:500])))
    read = transit.from_csv(str(path))
    assert len(read['time']) == 500
    assert read['flux'][3] == pytest.approx(lc['flux'][3])


def test_the_period_scanned():
    """a period off by enough that, over two sectors 400 days apart, its
    transits do not line up: folded at it, the box is weak; scanned within
    its error, the transits line up again, at the true period"""
    lc = _light_curve(seed=7, transits=[(3.1, 58001.3, 1.2, 1.6)])
    asked = 3.1 + 0.0008
    fixed = transit.search(lc, asked, mstar=0.4, rstar=0.4, ntrial=5)
    scanned = transit.search(lc, asked, period_err=0.0005, mstar=0.4,
                             rstar=0.4, ntrial=5)
    assert 'scan' in scanned and 'scan' not in fixed
    assert scanned['fold_period'] == pytest.approx(3.1, abs=1.5e-4)
    assert scanned['best']['snr'] > 1.3 * fixed['best']['snr']
    assert scanned['plausible'], scanned['why']
    assert scanned['best']['depth'] == pytest.approx(1.2, abs=0.25)
    # a period known well enough is not scanned
    assert 'scan' not in transit.search(lc, 3.1, period_err=1e-6,
                                        mstar=0.4, rstar=0.4, ntrial=5)


def test_the_chance_of_the_null():
    """a Gumbel fitted to the null: a value far beyond it is unlikely, one
    within it is not"""
    null = [3.1, 4.0, 3.5, 2.8, 4.4, 3.9, 3.0, 3.6, 4.1, 3.3]
    assert transit.null_chance(12.0, null) < 1e-6
    assert transit.null_chance(4.0, null) > 0.1


def test_a_conjunction_from_the_velocities():
    """a sinusoid's conjunction: its own error, and the offset an eccentric
    orbit gives it (P e / pi), in quadrature"""
    assert transit.conjunction_error(0.0, 2.644) == pytest.approx(
        2.644 * 0.1 / np.pi)
    assert transit.conjunction_error(0.3, 2.644) == pytest.approx(
        np.hypot(0.3, 2.644 * 0.1 / np.pi))


def test_the_box_fitted_to_the_medians():
    """a 1.2 ppt, 2.4 h transit: the box fitted to the medians gives its
    depth and duration (the search's box, the deepest of many, may not),
    and the gaps of TESS split the light curve, their ramps left out"""
    lc = _light_curve(seed=8, transits=[(4.3, 58002.1, 1.2, 2.4)])
    res = transit.search(lc, 4.3, 58002.1, 0.01, 1e-5, mstar=0.5,
                         rstar=0.5, ntrial=5)
    fit = res['fit']
    assert res['plausible'], res['why']
    assert fit['depth'] == pytest.approx(1.2, abs=0.15)
    assert fit['duration'] == pytest.approx(2.4, abs=0.5)
    assert abs(fit['centre']) < 0.5
    assert fit['radius'] == pytest.approx(transit.radius_of(1.2, 0.5),
                                          rel=0.1)
    trimmed = transit.trim_edges(lc)
    # each sector starts again twice (at its start, after its gap): 0.25 d
    #   left out each time
    assert len(lc['time']) - len(trimmed['time']) == pytest.approx(
        4 * 0.25 * 720, abs=8)

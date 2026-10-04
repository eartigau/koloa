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


def test_one_event_is_not_a_transit():
    """one deep dip (a systematic of one orbit): not a transit seen twice"""
    lc = _light_curve(seed=3)
    lc['flux'][np.abs(lc['time'] - 58410.0) < 0.05] -= 4.0
    res = transit.search(lc, 5.0, None, mstar=0.4, rstar=0.4, ntrial=10)
    assert not res['plausible']
    assert 'one transit' in res['why'] or 'single' in res['why']


def test_another_planet_left_out():
    """the transits of a known planet (2.2 d, deep) do not line up at a
    wrong period once left out"""
    lc = _light_curve(seed=4, transits=[(2.2, 58000.7, 3.0, 1.5)])
    other = [dict(P=2.2, tc=58000.7, duration=1.5)]
    res = transit.search(lc, 2.9, None, mstar=0.4, rstar=0.4, others=other,
                         ntrial=10)
    assert not res['plausible'], res['why']


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

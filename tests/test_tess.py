#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The TESS check, on light curves made up here (no network)

Created on 2026-09-30

@author: artigau
"""
import numpy as np

from koloa import tess


def fake(period, amp=2.0, nsec=3, noise=1.0, seed=1):
    """sectors of 27 days at 2-minute cadence, a sinusoid of amp ppt"""
    rng = np.random.default_rng(seed)
    sectors = []
    for isec in range(nsec):
        time_ = 60000.0 + 400 * isec + np.arange(0, 27, 2 / 1440)
        flux = amp * np.sin(2 * np.pi * time_ / period) \
            + rng.normal(0, noise, time_.size)
        sectors.append(dict(sector=10 + isec, provenance='SPOC',
                            exposure=120.0, column='PDCSAP_FLUX',
                            time=time_, flux=flux,
                            err=np.full(time_.size, noise)))
    return dict(tic=1, sectors=sectors)


def test_rotation_is_found_at_the_period_and_its_double():
    res = tess.periodicity(fake(4.0), [4.0, 8.0, 9.1, 60.0])
    status = {chk['period']: chk['status'] for chk in res['checks']}
    # at P, and at 2P (its half is the rotation): flagged
    assert status[4.0] == 'flag' and status[8.0] == 'flag'
    # 9.1 d: not the flank of the 4 d peak
    assert status[9.1] == 'pass'
    # 60 d and its harmonics are beyond a sector
    assert status[60.0] == 'info'
    assert abs(res['sectors'][0]['peaks'][0]['period'] - 4.0) < 0.2
    assert abs(res['sectors'][0]['peaks'][0]['amplitude'] - 2.0) < 0.2


def test_the_orbit_of_tess_is_a_note():
    res = tess.periodicity(fake(6.85), [6.85])
    assert res['checks'][0]['status'] == 'info'
    assert 'orbit of TESS' in res['checks'][0]['summary']


def test_flares_are_clipped():
    lcs = fake(4.0, nsec=1)
    sec = lcs['sectors'][0]
    sec['flux'][5000:5010] += 200.0
    res = tess.periodicity(lcs, [4.0])
    assert res['checks'][0]['status'] == 'flag'
    assert res['sectors'][0]['rms'] < 5


def test_duck_test_takes_light_curves(tmp_path):
    from koloa.diagnostics import duck_test
    from koloa.simulate import simulate
    sim = simulate(planets=[dict(P=4.0, K=5.0, e=0.0, tp=0.0)], seed=3)
    report = duck_test(sim['data'], 4.0, gp=False, quiet=True,
                       tess=fake(4.0))
    phot = [chk for chk in report.checks if chk['name'] == 'photometry']
    assert phot and phot[0]['status'] == 'flag'
    assert 'photometry' in [chk['name'] for chk in report.flags]
    # and in the report, a figure of its own
    report.pdf(None, sim['data'], archive=False, outdir=str(tmp_path))
    assert any(item.name.endswith('_tess.pdf') for item in tmp_path.iterdir())
    # offline by default without a report
    report = duck_test(sim['data'], 4.0, gp=False, quiet=True)
    assert 'photometry' not in [chk['name'] for chk in report.checks]


def test_a_star_tess_has_not_observed(monkeypatch):
    """GJ 1214: MAST has no light curve of it and no sector holds its
    position; said as such, not as a transit not found. A star in the full
    frames only is told apart (its sectors given)"""
    import io
    import json
    from koloa import archive, gui
    asked = []

    def urlopen(url, timeout=None):
        asked.append(url)
        # TESScut: no sector at GJ 1214, two at the other star
        rows = ([] if 'ra=258.82' in url else
                [dict(sector='0022'), dict(sector='0049')])
        return io.BytesIO(json.dumps(dict(results=rows)).encode())
    monkeypatch.setattr(tess, '_mast', lambda service, params, timeout=60:
                        [])
    monkeypatch.setattr(tess.urllib.request, 'urlopen', urlopen)
    monkeypatch.setattr(archive, 'resolve', lambda name, **kwargs: (
        dict(tic='TIC 467929202', ra=258.8289, dec=4.9639)
        if '1214' in name else dict(tic='TIC 1', ra=175.5, dec=26.7)))
    assert tess.sectors_at(258.8289, 4.9639) == []
    assert tess.sectors_at(None, 4.9639) is None
    out = tess.light_curves('GJ 1214')
    assert out['sectors'] == [] and out['observed'] == []
    assert tess.light_curves('Other star')['observed'] == [22, 49]
    # the page: nothing to look for a transit in, and why
    monkeypatch.setattr(gui, 'TESS_LC', {})
    monkeypatch.setattr(gui, 'TESS_WHY', {})
    res = gui.transit_check(dict(target='GJ 1214', root='nowhere',
                                 period=1.5804, fetch=True, name='#1'))
    assert res['missing'] and res['reason'] == 'unobserved'
    assert 'not observed' in res['lc']
    res = gui.transit_check(dict(target='Other star', root='nowhere',
                                 period=1.5804, fetch=True, name='#1'))
    assert res['reason'] == 'frames' and res['sectors'] == [22, 49]

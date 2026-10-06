#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The light curves of Kepler, K2 and CoRoT (koloa.photometry), on files made
up here (no network): each mission's file read, the star found by its
position, K2SFF taken before the mission's own flux, a transit found in
long cadence, and the mission with a transit shown by the page

Created on 2026-10-06

@author: artigau
"""
import io
import shutil

import numpy as np
import pytest

from koloa import gui, photometry, transit

fits = pytest.importorskip('astropy.io.fits')


def _kepler_file(path, quarter=None, campaign=None, start=120.0, days=30.0,
                 depth=0.0, seed=1):
    """a long-cadence light curve as Kepler and K2 write them (times from
    BJD 2454833), a transit every 3.5 d when asked, a few cadences flagged
    and wild"""
    rng = np.random.default_rng(seed)
    time_ = start + np.arange(0, days, 29.4244 / 1440)
    flux = 1e5 * (1 + rng.normal(0, 2e-4, time_.size))
    phase = ((time_ - 121.0) / 3.5 + 0.5) % 1.0 - 0.5
    flux[np.abs(phase * 3.5 * 24) < 1.5] *= 1 - depth
    qual = np.zeros(time_.size, int)
    qual[::97] = 128
    flux[::97] *= 3.0
    cols = [fits.Column('TIME', 'D', array=time_),
            fits.Column('PDCSAP_FLUX', 'E', array=flux),
            fits.Column('PDCSAP_FLUX_ERR', 'E', array=np.full(time_.size,
                                                             20.0)),
            fits.Column('SAP_QUALITY', 'J', array=qual)]
    head = fits.Header()
    if quarter is not None:
        head['QUARTER'] = quarter
    if campaign is not None:
        head['CAMPAIGN'] = campaign
    head['OBSMODE'] = 'long cadence'
    table = fits.BinTableHDU.from_columns(cols)
    table.header['BJDREFI'], table.header['BJDREFF'] = 2454833, 0.0
    fits.HDUList([fits.PrimaryHDU(header=head), table]).writeto(
        path, overwrite=True)
    return time_


def _sff_file(path, start=1977.0, seed=2):
    """a K2SFF light curve: its best aperture, a few cadences moving"""
    rng = np.random.default_rng(seed)
    time_ = start + np.arange(0, 60.0, 29.4244 / 1440)
    flux = 1 + rng.normal(0, 1e-4, time_.size)
    moving = np.zeros(time_.size, int)
    moving[::50] = 1
    table = fits.BinTableHDU.from_columns([
        fits.Column('T', 'D', array=time_),
        fits.Column('FRAW', 'D', array=flux),
        fits.Column('FCOR', 'D', array=flux),
        fits.Column('MOVING', 'I', array=moving)], name='BESTAPER')
    table.header['BJDREFI'], table.header['BJDREFF'] = 2454833, 0.0
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(path, overwrite=True)
    return time_


def _corot_file(path, bright=False, seed=3):
    """a CoRoT N2 file: a faint star (its flux corrected for systematics,
    512 s) or a bright one (its regular flux, 32 s); a tenth of its points
    flagged"""
    rng = np.random.default_rng(seed)
    step = (32.0 if bright else 512.0) / 86400.0
    time_ = 54400.0 + np.arange(0, 20.0, step)
    flux = 1e6 * (1 + rng.normal(0, 5e-4, time_.size))
    status = np.zeros(time_.size, int)
    status[::10] = 256
    flux[::10] *= 1.5
    if bright:
        table = fits.BinTableHDU.from_columns([
            fits.Column('DATEBARREGTT', 'D', array=time_),
            fits.Column('FLUXBARREG', 'D', array=flux),
            fits.Column('STATUSBARREG', 'J', array=status)], name='BARREG')
    else:
        table = fits.BinTableHDU.from_columns([
            fits.Column('DATEBARTT', 'D', array=time_),
            fits.Column('WHITEFLUXSYS', 'E', array=flux),
            fits.Column('STATUSSYS', 'J', array=status)], name='SYSTEMATIC')
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(path, overwrite=True)
    return time_


def test_the_name_of_a_stretch():
    """a sector of TESS, a quarter of Kepler, a campaign of K2, a run of
    CoRoT: told apart by their number, each with its name"""
    assert photometry.label(22) == 'sector 22'
    assert photometry.label(1003) == 'Kepler Q3'
    assert photometry.label(2005) == 'K2 C5'
    assert photometry.label(3001, 'LRa01') == 'CoRoT LRa01'
    assert [photometry.mission_of(num) for num in (22, 1000, 2018, 3002)] \
        == ['tess', 'kepler', 'k2', 'corot']


def test_each_file_is_read(tmp_path):
    """Kepler's PDCSAP flux (its flagged cadences left out, its times BJD
    - 2400000, ppt about its median), K2SFF's best aperture (not where the
    spacecraft moved), CoRoT's corrected flux, and a bright star's binned
    to 512 s"""
    path = str(tmp_path / 'kplr.fits')
    time_ = _kepler_file(path, quarter=3)
    lc = photometry._read_kepler(path)
    assert lc['quarter'] == 3 and lc['column'] == 'PDCSAP_FLUX'
    assert len(lc['time']) == time_.size - len(time_[::97])
    assert lc['time'][0] == pytest.approx(time_[1] + 54833.0)
    assert abs(np.median(lc['flux'])) < 0.05 and np.std(lc['flux']) < 0.5
    assert lc['err'][0] == pytest.approx(0.2, rel=0.01)
    path = str(tmp_path / 'sff.fits')
    time_ = _sff_file(path)
    lc = photometry._read_k2sff(path)
    assert len(lc['time']) == time_.size - len(time_[::50])
    assert lc['column'].startswith('FCOR') and np.all(lc['err'] > 0)
    path = str(tmp_path / 'faint.fits')
    time_ = _corot_file(path)
    lc = photometry._read_corot(path)
    assert lc['column'] == 'SYSTEMATIC.WHITEFLUXSYS'
    assert len(lc['time']) == time_.size - len(time_[::10])
    assert np.std(lc['flux']) < 1.0 and lc['exposure'] > 500
    path = str(tmp_path / 'bright.fits')
    _corot_file(path, bright=True)
    lc = photometry._read_corot(path)
    assert lc['exposure'] == photometry.COROT_BIN
    assert np.median(np.diff(lc['time'])) * 86400 == pytest.approx(512, 1)
    # means of 16 points: four times less scatter
    assert np.std(lc['flux']) < 0.2


def _archive(tmp_path, monkeypatch):
    """MAST and VizieR made up: a star with two quarters of Kepler (a
    neighbour 6 arcsec away too), one campaign of K2 (the mission's and
    K2SFF's) and one run of CoRoT"""
    src = tmp_path / 'src'
    src.mkdir()
    _kepler_file(str(src / 'kplr001-2009_llc.fits'), quarter=1, start=131.0,
                 depth=1e-3)
    _kepler_file(str(src / 'kplr001-2010_llc.fits'), quarter=5, start=443.0,
                 depth=1e-3, seed=4)
    _kepler_file(str(src / 'kplr002-2009_llc.fits'), quarter=1, start=131.0)
    _kepler_file(str(src / 'ktwo009-c01_llc.fits'), campaign=1, start=1977.0)
    _sff_file(str(src / 'hlsp_k2sff_k2_lightcurve_009-c01_kepler_v1_llc.fits'))
    _corot_file(str(src / 'EN2_STAR_CHR_0000000077_A_B.fits'))
    rows = [
        dict(obsid=1, obs_collection='Kepler', provenance_name='Kepler',
             target_name='kplr001', t_exptime=1800, distance=0.3),
        dict(obsid=2, obs_collection='Kepler', provenance_name='Kepler',
             target_name='kplr002', t_exptime=1800, distance=6.0),
        dict(obsid=3, obs_collection='Kepler', provenance_name='Kepler',
             target_name='kplr001', t_exptime=60, distance=0.3),
        dict(obsid=4, obs_collection='K2', provenance_name='K2',
             target_name='ktwo009', t_exptime=1800, distance=0.5),
        dict(obsid=5, obs_collection='HLSP', provenance_name='K2SFF',
             target_name='ktwo009', t_exptime=1800, distance=0.5),
        dict(obsid=6, obs_collection='HLSP', provenance_name='EVEREST',
             target_name='ktwo009', t_exptime=1800, distance=0.5)]
    products = {
        1: ['kplr001-2009_llc.fits', 'kplr001-2010_llc.fits',
            'kplr001_dvs.pdf'],
        2: ['kplr002-2009_llc.fits'], 3: ['kplr001-2009_slc.fits'],
        4: ['ktwo009-c01_llc.fits'],
        5: ['hlsp_k2sff_k2_lightcurve_009-c01_kepler_v1_llc.fits'],
        6: ['hlsp_everest_k2_llc_009-c01_kepler_v2.0_lc.fits']}
    asked = []

    def mast(service, params, timeout=60):
        asked.append(service)
        if service.endswith('Position'):
            return rows
        ids = [int(num) for num in str(params['obsid']).split(',')]
        return [dict(obsID=num, productFilename=name, dataURI=f'mast:{name}')
                for num in ids for name in products[num]]

    def download(url, path, timeout=300.0):
        name = url.rsplit(':', 1)[-1].rsplit('/', 1)[-1]
        import os
        os.makedirs(os.path.dirname(path), exist_ok=True)
        shutil.copy(str(src / name), path)

    def urlopen(url, timeout=None):
        # VizieR's log of CoRoT: one run in the table of the faint stars
        text = '#\n' if 'Bright' in url else (
            'CoRoT\tRun\tFileName\tdate1\n\t\t\t\n-----\t---\t---\t---\n'
            '77\tLRa01\tN2/EN2_STAR_CHR_0000000077_A_B.fits\t2007-10-23\n')
        return io.BytesIO(text.encode())
    from koloa import tess
    monkeypatch.setattr(tess, '_mast', mast)
    monkeypatch.setattr(photometry, '_download', download)
    monkeypatch.setattr(photometry.urllib.request, 'urlopen', urlopen)
    monkeypatch.setattr(photometry, 'CACHE', str(tmp_path / 'cache'))
    monkeypatch.setattr(photometry, '_PRODUCTS', {})
    monkeypatch.setattr(photometry, '_COROT_LOG', {})
    return asked


def test_a_star_in_their_fields(tmp_path, monkeypatch):
    """asked by position: the nearest target's long cadence (not its
    neighbour's, not its short cadence), K2SFF before the mission's flux,
    the run of CoRoT; which missions have the star, asked once; a star
    without a position is in none"""
    asked = _archive(tmp_path, monkeypatch)
    ident = dict(main='A star', ra=290.0, dec=44.0)
    assert photometry.available(ident) == dict(kepler=2, k2=1, corot=1)
    lcs = photometry.light_curves(ident)
    kep = lcs['kepler']
    assert kep['target'] == 'kplr001'
    assert [lc['label'] for lc in kep['sectors']] == ['Kepler Q1',
                                                      'Kepler Q5']
    assert [lc['sector'] for lc in kep['sectors']] == [1001, 1005]
    assert [lc['provenance'] for lc in lcs['k2']['sectors']] == ['K2SFF']
    assert lcs['k2']['sectors'][0]['label'] == 'K2 C1'
    run = lcs['corot']['sectors'][0]
    assert lcs['corot']['target'] == '77' and run['label'] == 'CoRoT LRa01'
    assert run['sector'] == 3001 and run['file'].endswith('LRa01_77.csv')
    # MAST's position asked once for all of it
    assert asked.count('Mast.Caom.Filtered.Position') == 1
    # the run of CoRoT read back from what was kept of it (not its file)
    monkeypatch.setattr(photometry, '_download', lambda *a, **k: 1 / 0)
    again = photometry.corot_light_curves(ident)['sectors'][0]
    assert np.allclose(again['flux'], run['flux'], atol=1e-4)
    assert again['column'] == run['column']
    none = photometry.light_curves(dict(main='Nowhere', ra=None, dec=None))
    assert all(not val['sectors'] for val in none.values())
    # a transit in two quarters of long cadence: found, at its depth
    lc = transit.from_light_curves(kep)
    res = transit.search(lc, 3.5, 54833.0 + 121.0, 0.01, 1e-5, mstar=1.0,
                         rstar=1.0, ntrial=5)
    assert res['plausible'], res['why']
    assert res['best']['depth'] == pytest.approx(1.0, abs=0.2)
    assert sorted(set(res['sector'])) == [1001, 1005]


def test_the_archives_keep_them(tmp_path, monkeypatch):
    """koloa.gather: a CSV per quarter, campaign and run, one per mission
    (its stretches by their number), read back with their names"""
    from koloa import gather
    _archive(tmp_path, monkeypatch)
    ident = dict(main='A star', ra=290.0, dec=44.0)
    folder = tmp_path / 'phot'
    lcs = gather.space_photometry(ident, 'A star', str(folder))
    assert (folder / 'kepler' / 'Kepler_Q1.csv').exists()
    assert (folder / 'corot' / 'CoRoT_LRa01.csv').exists()
    for key in photometry.MISSIONS:
        assert (folder / f'{key}.csv').exists()
        read = transit.from_csv(str(folder / f'{key}.csv'))
        assert len(read['time']) == sum(len(lc['time'])
                                        for lc in lcs[key]['sectors'])
    assert gui._csv_labels(str(folder / 'kepler.csv')) == {
        1001: 'Kepler Q1', 1005: 'Kepler Q5'}
    assert gui._csv_labels(str(folder / 'corot.csv')) == {
        3001: 'CoRoT LRa01'}


def test_the_mission_with_a_transit_is_shown(monkeypatch):
    """the page: TESS searched first, then the other missions that have
    the star, until one has a plausible transit (Kepler's here: TESS's
    light curve is too noisy for it); a mission asked by its name; none
    has the star"""
    rng = np.random.default_rng(5)

    def curve(sector, white, depth):
        time_ = 59000.0 + np.arange(0, 60.0, 10.0 / 1440)
        time_ = time_[np.abs(time_ - 59030.0) > 1.0]
        flux = rng.normal(0, white, time_.size)
        phase = ((time_ - 59001.0) / 3.5 + 0.5) % 1.0 - 0.5
        flux[np.abs(phase * 3.5 * 24) < 1.2] -= depth
        return dict(time=time_, flux=flux, err=np.full(time_.size, white),
                    sector=np.full(time_.size, float(sector)))
    have = dict(tess=(curve(40, 3.0, 0.3), 'a test', {}),
                kepler=(curve(1003, 0.1, 0.3), 'a test',
                        {1003: 'Kepler Q3'}))

    def light(target, root='', fetch=False, mission='tess'):
        return have.get(mission, (None, 'not in its fields', {}))
    monkeypatch.setattr(gui, 'space_light', light)
    monkeypatch.setattr(gui, 'resolve_star', lambda target, root='',
                        refresh=False: dict(star=dict(
                            mass=1.0, radius=1.0, radius_source='a test')))
    monkeypatch.setattr(gui, 'known_periods', lambda target: [])
    monkeypatch.setattr(gui, 'transits_of', lambda target, known: [])
    monkeypatch.setattr(gui, '_missions_with', lambda target, root, used:
                        ['tess', 'kepler'])
    opts = dict(target='A star', period=3.5, tc=59001.0, tc_err=0.01,
                p_err=1e-5, fetch=True, name='b')
    res = gui.transit_check(opts)
    assert res['mission'] == 'kepler' and res['mission_name'] == 'Kepler'
    assert res['plausible'] and res['stretches'] == 'Kepler Q3'
    assert res['labels'] == {'1003': 'Kepler Q3'}
    assert res['searched']['tess']['plausible'] is False
    assert res['searched']['kepler']['plausible'] is True
    assert res['missions'] == ['tess', 'kepler']
    # asked of TESS alone: its own answer
    alone = gui.transit_check(dict(opts, mission='tess'))
    assert alone['mission'] == 'tess' and not alone['plausible']
    assert alone['stretches'] == 'TESS 40'
    # no mission has the star
    have.clear()
    monkeypatch.setattr(gui, 'tess_why', lambda target, root='':
                        dict(reason='unobserved'))
    none = gui.transit_check(opts)
    assert none['missing'] and none['others_none']
    assert gui._stretches('kepler', list(range(1000, 1018)), {
        num: f'Kepler Q{num - 1000}' for num in range(1000, 1018)}) == \
        'Kepler Q0 to Q17 (18)'

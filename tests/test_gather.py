#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The folder of a star's archives (koloa.gather), without the network

Created on 2026-10-01

@author: artigau
"""
import csv
import json
import os

import numpy as np

from koloa import gather as kg
from koloa.data import RVData


def _series(seed=1):
    rng = np.random.default_rng(seed)
    time = np.sort(rng.uniform(58000, 59000, 30))
    inst = np.where(np.arange(30) < 12, 'HARPS15', 'NIRPS')
    rv = np.where(inst == 'HARPS15', -21000.0, -20950.0) + rng.normal(0, 3, 30)
    fwhm = (rng.normal(6000, 5, 30), np.full(30, 2.0))
    return RVData(time, rv, np.full(30, 1.5), inst=inst,
                  indicators=dict(fwhm=fwhm), name='star')


def test_folder_name():
    assert kg.folder_name('GJ 436') == 'GJ_436'
    assert kg.folder_name(' TOI-700 ') == 'TOI-700'


def test_a_series_written_and_read_back(tmp_path):
    data = _series()
    path = kg.write_rv(data, str(tmp_path / 'rv.csv'), source='dace')
    back = RVData.from_csv(path, inst='inst')
    assert back.instruments == data.instruments
    assert np.allclose(back.time, data.time)
    # the velocities as given: the zero points put back on both sides
    for name in data.instruments:
        assert np.isclose(back.zero_point[name], data.zero_point[name])
    assert np.allclose(back.rv, data.rv) and np.allclose(back.err, data.err)
    assert np.allclose(back.indicators['fwhm'][0], data.indicators['fwhm'][0])


def test_carmenes_keeps_the_corrected_velocities(tmp_path):
    folder = tmp_path / 'carmenes'
    folder.mkdir()
    rows = []
    for it in range(8):
        row = {key: '' for key in kg.CARMENES_COLUMNS}
        row.update(bjd=2458000.5 + it, rv=10.0 + it, err_rv=1.0,
                   dlw=5.0 + it, err_dlw=1.0)
        if it != 3:
            # one exposure without the nightly zero point: left out
            row.update(a_rv=30.0 + it, err_a_rv=1.2)
        rows.append(row)
    with open(folder / 'raw.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=kg.CARMENES_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    data = kg.carmenes_rv(dict(carmenes_id='J00000+000'), 'star',
                          str(folder))
    assert data.n == 7 and data.instruments == ['CARMENES']
    assert np.isclose(data.time[0], 58000.5)
    assert 'dlw' in data.indicators
    assert os.path.exists(folder / 'CARMENES.csv')


def test_load_reads_what_gather_wrote(tmp_path):
    folder = tmp_path / 'star'
    kg.write_rv(_series(), str(folder / 'rv' / 'all_rv.csv'))
    os.makedirs(folder / 'phot' / 'tess')
    with open(folder / 'phot' / 'tess' / 's0022_SPOC.csv', 'w') as handle:
        handle.write('rjd,flux,sflux\n58900.1,0.5,0.2\n58900.2,-0.3,0.2\n')
    manifest = dict(target='star', archives=dict(tess=dict(
        status='ok', sectors=[dict(sector=22, pipeline='SPOC',
                                   file='phot/tess/s0022_SPOC.csv')])))
    with open(folder / 'manifest.json', 'w') as handle:
        json.dump(manifest, handle)
    star = kg.load(str(folder))
    assert star['rv'].n == 30 and len(star['tess']) == 1
    assert np.allclose(star['tess'][0]['flux'], [0.5, -0.3])


def test_the_dace_key_is_looked_for(tmp_path, monkeypatch):
    from koloa import dace
    rc = tmp_path / '.dacerc'
    rc.write_text('[user]\nkey = apiKey:abc123\n')
    monkeypatch.setattr(dace, 'DACERC', str(rc))
    monkeypatch.delenv('DACE_API_KEY', raising=False)
    assert dace.find_key() == ('abc123', '~/.dacerc')
    monkeypatch.setenv('DACE_API_KEY', 'xyz')
    assert dace.find_key() == ('xyz', 'DACE_API_KEY')
    # set to nothing, or False: the public data only
    monkeypatch.setenv('DACE_API_KEY', '')
    assert dace.find_key() == (None, None)
    assert dace.find_key(False) == (None, None)
    assert dace.find_key('given') == ('given', 'given')


def test_a_carmenes_star_found_in_the_list_kept(tmp_path, monkeypatch):
    """the stars of CARMENES DR1 kept on disk: a star found by its position,
    without the network"""
    kept = tmp_path / 'carmenes_objects.json'
    kept.write_text(json.dumps([
        dict(carmenes_id='J11421+267', name='Ross 905', ra='175.54622',
             dec='26.70657', p_rot='44.6'),
        dict(carmenes_id='J00067-075', name='GJ 1002', ra='1.67',
             dec='-7.54', p_rot='')]))
    monkeypatch.setattr(kg, 'CARMENES_CACHE', str(kept))
    monkeypatch.setattr(kg, '_tap', lambda *a, **k: 1 / 0)
    assert kg.carmenes_star(175.5463, 26.7065)['carmenes_id'] == 'J11421+267'
    assert kg.carmenes_star(10.0, 10.0) is None


def test_the_velocities_alone_or_a_sector_gone(tmp_path):
    """what gather put in a folder, read without its photometry, or with a
    sector the manifest lists that the folder no longer has"""
    import json
    import numpy as np
    from koloa.data import RVData
    from koloa.gather import load, write_rv
    write_rv(RVData(np.arange(5.0), np.zeros(5), np.ones(5),
                    inst=np.array(['HARPS'] * 5)),
             str(tmp_path / 'rv' / 'all_rv.csv'))
    (tmp_path / 'manifest.json').write_text(json.dumps(dict(
        target='X', archives=dict(tess=dict(sectors=[
            dict(sector=1, pipeline='SPOC', file='phot/tess/s0001.csv')])))))
    got = load(str(tmp_path), photometry=False)
    assert got['rv'].n == 5 and got['tess'] == []
    assert load(str(tmp_path))['tess'] == []   # the sector is gone


def test_what_an_archive_gave_in_words():
    """the line under a step of a gather in koloa's GUI: the points of an
    archive by instrument, the sectors of TESS, or why nothing"""
    from koloa.gather import _told
    assert _told(dict(status='ok', npoints=260, instruments=dict(
        HARPS03=19, HARPS15=114, NIRPS=127))) == (
        '260 points: HARPS03 19, HARPS15 114, NIRPS 127')
    assert _told(dict(status='ok', sectors=[
        dict(sector=13, pipeline='SPOC', npoints=100),
        dict(sector=66, pipeline='QLP', npoints=50)])) == (
        '2 sectors, 150 points (s13 SPOC, s66 QLP)')
    assert _told(dict(status='none', message='not in DR1')) == (
        'nothing: not in DR1')

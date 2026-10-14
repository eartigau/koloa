#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LBL from the spectra of the archives (koloa.lblrun), without the network
and without LBL: the parameters, the script, the spectra fetched, the rdb
kept, and the tab of the GUI

Created on 2026-10-14

@author: artigau
"""
import gzip
import json
import os
import sys

import pytest

from koloa import lblrun


def test_the_parameters_of_a_run(tmp_path):
    params = lblrun.parameters('GJ 436', 'HARPSN', str(tmp_path), 3420.0,
                               nproc=6)
    assert params['INSTRUMENT'] == 'HARPSN' and params['DATA_SOURCE'] == 'ESO'
    assert params['OBJECT_SCIENCE'] == ['GJ436'] == params['OBJECT_COMPARISON']
    assert params['OBJECT_TEFF'] == [3420] and params['COMPUTE_NPROC'] == 6
    # the gradient of temperature nearest the star's
    assert params['RESPROJ_TABLES'] == {
        'DTEMP3500': 'temperature_gradient_3500.fits'}
    assert params['RUN_LBL_TELLUCLEAN'] is True
    assert lblrun.parameters('GJ 436', 'HARPSN', str(tmp_path), 2500.0,
                             tellu=False)['OBJECT_TEFF'] == [3000]
    # the folder of the instrument, its models shared
    here = params['DATA_DIR']
    assert here.endswith('HARPSN_ESO') and os.path.islink(
        os.path.join(here, 'models'))
    assert lblrun.rdb_of(params).endswith('lblrdb/lbl_GJ436_tc_GJ436_tc.rdb')
    # the script: LBL's wrapper with these parameters
    text = lblrun.wrap(params)
    assert 'from lbl import lbl_wrap' in text and 'lbl_wrap(rparams)' in text
    space = {}
    exec(text.replace('from lbl import lbl_wrap', 'lbl_wrap = print')
         .split('if __name__')[0], space)
    assert space['rparams'] == params


def test_a_run_and_its_rdb_kept(tmp_path):
    """LBL as a process of its own (here a Python that writes the rdb):
    its lines said, its rdb kept in the folder of the star"""
    params = lblrun.parameters('GJ 436', 'HARPSN', str(tmp_path / 'lbl'),
                               3400.0)
    fake = tmp_path / 'fakelbl'
    (fake / 'lbl').mkdir(parents=True)
    (fake / 'lbl' / '__init__.py').write_text(
        'import json, os, sys\n'
        'def lbl_wrap(params):\n'
        '    print("\\x1b[92;1mcomputing\\x1b[0;0m")\n'
        '    obj = params["OBJECT_SCIENCE"][0] + "_tc"\n'
        '    out = os.path.join(params["DATA_DIR"], "lblrdb")\n'
        '    os.makedirs(out, exist_ok=True)\n'
        '    open(os.path.join(out, f"lbl_{obj}_{obj}.rdb"), "w").write(\n'
        '        "rjd\\tvrad\\tsvrad\\n---\\t----\\t-----\\n1\\t2\\t3\\n")\n')
    os.environ['PYTHONPATH'] = str(fake)
    try:
        said = []
        rdb = lblrun.run(params, sys.executable, echo=said.append)
    finally:
        del os.environ['PYTHONPATH']
    assert said == ['computing']
    kept = lblrun.keep(rdb, 'GJ 436', 'HARPSN', str(tmp_path / 'archives'))
    assert kept.endswith('archives/GJ_436/lbl/GJ_436_HARPSN_lbl.rdb')
    assert open(kept).read().splitlines()[2] == '1\t2\t3'
    # a run that fails says why
    os.environ['PYTHONPATH'] = str(tmp_path / 'nothing')
    try:
        with pytest.raises(RuntimeError, match='LBL stopped'):
            lblrun.run(params, sys.executable, echo=said.append)
    finally:
        del os.environ['PYTHONPATH']


def test_the_spectra_fetched(tmp_path, monkeypatch):
    """the S2D of the exposures kept, gunzipped, those there not fetched
    again, a failure said"""
    asked = []

    def urlopen(request, timeout):
        import io
        asked.append(request.full_url)
        if 'bad' in request.full_url:
            raise OSError('gone')
        return io.BytesIO(gzip.compress(b'SIMPLE  =                    T'))
    monkeypatch.setattr(lblrun.urllib.request, 'urlopen', urlopen)
    names = [f'r.HARPN.2016-03-0{it}T00-00-00.000_S2D_A.fits.gz'
             for it in range(1, 4)]
    got = lblrun.fetch('GJ 436', 'HARPSN', str(tmp_path), names=names)
    assert got['n'] == 3 and got['fetched'] == 3 and not got['failed']
    assert asked[0] == lblrun.TNG_FILES + names[0]
    there = os.path.join(got['folder'], names[0][:-3])
    assert open(there, 'rb').read().startswith(b'SIMPLE')
    asked.clear()
    got = lblrun.fetch('GJ 436', 'HARPSN', str(tmp_path),
                       names=names + ['r.HARPN.bad_S2D_A.fits.gz'])
    assert len(asked) == 1 and got['failed'] == 1 and 'gone' in got['why']
    assert not [name for name in os.listdir(got['folder'])
                if name.endswith('.part')]


def test_the_spectra_of_the_exposures_kept(tmp_path):
    """the S2D of the espdr exposures koloa keeps from their CCFs: the
    companion pointed at in the same place left out by its velocity"""
    from koloa import tng
    folder = tmp_path / 'archives' / 'GJ_1' / 'rv' / 'published'
    folder.mkdir(parents=True)
    rows = []
    for it in range(6):
        name = 'GJ1' if it < 4 else 'GJ1B'
        row = {key: '' for key in tng.COLUMNS}
        row.update(file=f'r.HARPN.2016-03-0{it + 1}T00-00-00.000_CCF_A.fits.gz',
                   reduction='espdr', policy='FREE', object=name,
                   bjd=2457451.5 + it, rv=10.0 if it < 4 else 30.0,
                   err=0.001, mask='M2', dpr='STAR,WAVE', qc=1,
                   version='espdr/3.2.0')
        rows.append(row)
    tng._keep(rows, str(folder / 'tng_exposures.csv'))
    names = lblrun.harpsn_spectra('GJ 1', str(tmp_path / 'archives'))
    assert names == [f'r.HARPN.2016-03-0{it + 1}T00-00-00.000_S2D_A.fits.gz'
                     for it in range(4)]


def test_the_python_of_lbl(tmp_path, monkeypatch):
    monkeypatch.setattr(lblrun, 'PYTHONS', [])
    assert lblrun.find_python() == ''
    assert not lblrun.check_python('')['ok']
    lblrun.save_settings(python=sys.executable, nproc=3)
    assert lblrun.find_python() == sys.executable
    assert lblrun.settings()['nproc'] == 3
    # a Python without LBL: why said
    fake = tmp_path / 'python'
    fake.write_text('#!/bin/sh\necho "No module named lbl" >&2\nexit 1\n')
    fake.chmod(0o755)
    state = lblrun.check_python(str(fake))
    assert not state['ok'] and state['error'] == 'No module named lbl'


def test_the_tab_of_the_gui(tmp_path, monkeypatch):
    """the command line of the tab, its state, one LBL at a time"""
    from koloa import gui
    monkeypatch.setattr(lblrun, 'check_python', lambda python: dict(
        ok=True, version='0.67.008', where='/x/lbl', error=''))
    opts = dict(target='GJ 436', root=str(tmp_path / 'archives'),
                lbl_dir=str(tmp_path / 'lbl'), lbl_HARPSN=True, nproc='6',
                teff='3400', tellu=False)
    args = gui.command('lbl', opts)
    assert args == ['GJ 436', '--lbl', 'HARPSN', '--archives',
                    str(tmp_path / 'archives'), '--lbl-dir',
                    str(tmp_path / 'lbl'), '--nproc', '6', '--teff', '3400',
                    '--no-tellu']
    with pytest.raises(ValueError, match='tick an instrument'):
        gui.command('lbl', dict(opts, lbl_HARPSN=False))
    with pytest.raises(ValueError, match='SIMBAD name'):
        gui.command('lbl', dict(opts, target=''))
    # the state: the spectra there, the rdb kept
    science = tmp_path / 'lbl' / 'HARPSN_ESO' / 'science' / 'GJ436'
    science.mkdir(parents=True)
    (science / 'a_S2D_A.fits').write_text('x')
    rdb = tmp_path / 'archives' / 'GJ_436' / 'lbl' / 'GJ_436_HARPSN_lbl.rdb'
    rdb.parent.mkdir(parents=True)
    rdb.write_text('rjd\tvrad\n---\t----\n1\t2\n3\t4\n')
    state = gui.lbl_state('GJ 436', str(tmp_path / 'archives'),
                          str(tmp_path / 'lbl'))
    one = state['instruments'][0]
    assert one['there'] == 1 and one['npoints'] == 2
    assert one['dataset'] == 'HARPS-N (LBL)' and state['check']['ok']
    saved = gui.lbl_save(dict(nproc='5', data=str(tmp_path / 'lbl')))
    assert saved['nproc'] == 5
    assert json.load(open(lblrun.SETTINGS))['data'] == str(tmp_path / 'lbl')

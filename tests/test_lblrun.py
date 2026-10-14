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


def _group(names=('a.fits.gz',)):
    """a run of HARPS-N, as groups() makes it"""
    return dict(key='', obj='GJ436', lbl='HARPSN', source='ESO',
                input='*S2D_A.fits', tellu=True, params={},
                spectra=[dict(id=name, name=name[:-3]) for name in names],
                eras=[('', 'HARPS-N (LBL)') + lblrun.ALWAYS],
                get=lblrun._get_tng)


def test_the_parameters_of_a_run(tmp_path):
    params = lblrun.parameters('GJ 436', _group(), str(tmp_path), 3420.0,
                               nproc=6)
    assert params['INSTRUMENT'] == 'HARPSN' and params['DATA_SOURCE'] == 'ESO'
    assert params['OBJECT_SCIENCE'] == ['GJ436'] == params['OBJECT_COMPARISON']
    assert params['OBJECT_TEFF'] == [3420] and params['COMPUTE_NPROC'] == 6
    # the gradient of temperature nearest the star's
    assert params['RESPROJ_TABLES'] == {
        'DTEMP3500': 'temperature_gradient_3500.fits'}
    assert params['RUN_LBL_TELLUCLEAN'] is True
    assert lblrun.parameters('GJ 436', _group(), str(tmp_path), 2500.0,
                             tellu=False)['OBJECT_TEFF'] == [3000]
    # no cleaning where LBL does none (NIRPS, SPIRou), whatever is asked
    assert lblrun.parameters('GJ 436', dict(_group(), tellu=False),
                             str(tmp_path), 3400.0, tellu=True
                             )['RUN_LBL_TELLUCLEAN'] is False
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
    params = lblrun.parameters('GJ 436', _group(), str(tmp_path / 'lbl'),
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
    kept = lblrun.keep(rdb, 'GJ 436', 'HARPSN', _group(),
                       str(tmp_path / 'archives'))
    assert kept[0][0].endswith('archives/GJ_436/lbl/GJ_436_HARPSN_lbl.rdb')
    assert kept[0][1:] == ('HARPS-N (LBL)', 1)
    assert open(kept[0][0]).read().splitlines()[2] == '1\t2\t3'
    # a file per era (SOPHIE and SOPHIE+), their names beside them
    two = dict(_group(), eras=[('', 'SOPHIE (LBL)', -1e9, 2.0),
                               ('p', 'SOPHIE+ (LBL)', 2.0, 1e9)])
    with open(rdb, 'a') as handle:
        handle.write('5\t6\t7\n')
    kept = lblrun.keep(rdb, 'GJ 436', 'SOPHIE', two,
                       str(tmp_path / 'archives'))
    assert [(os.path.basename(path), label, n) for path, label, n in kept] \
        == [('GJ_436_SOPHIE_lbl.rdb', 'SOPHIE (LBL)', 1),
            ('GJ_436_SOPHIEp_lbl.rdb', 'SOPHIE+ (LBL)', 1)]
    listed = lblrun.kept('GJ 436', str(tmp_path / 'archives'))
    assert [(one['inst'], one['label'], one['npoints']) for one in listed] \
        == [('HARPSN', 'HARPS-N (LBL)', 1), ('SOPHIE', 'SOPHIE (LBL)', 1),
            ('SOPHIE', 'SOPHIE+ (LBL)', 1)]
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
    got = lblrun.fetch('GJ 436', 'HARPSN', _group(names), str(tmp_path))
    assert got['n'] == 3 and got['fetched'] == 3 and not got['failed']
    assert asked[0] == lblrun.TNG_FILES + names[0]
    there = os.path.join(got['folder'], names[0][:-3])
    assert open(there, 'rb').read().startswith(b'SIMPLE')
    asked.clear()
    got = lblrun.fetch('GJ 436', 'HARPSN', _group(
        names + ['r.HARPN.bad_S2D_A.fits.gz']), str(tmp_path))
    # the failure asked once more, then said
    assert len(asked) == 2 and got['failed'] == 1 and 'gone' in got['why']
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
    assert one['there'] == 1 and one['files'][0]['npoints'] == 2
    assert one['files'][0]['label'] == 'HARPS-N (LBL)'
    assert state['check']['ok']
    saved = gui.lbl_save(dict(nproc='5', data=str(tmp_path / 'lbl')))
    assert saved['nproc'] == 5
    assert json.load(open(lblrun.SETTINGS))['data'] == str(tmp_path / 'lbl')


def test_the_spectra_of_sophie(tmp_path, monkeypatch):
    """an e2ds of SOPHIE: its spectrum written with the name of its blaze
    in its header, the blaze beside it; a protected one left out"""
    from astropy.io import fits
    import numpy as np

    def e2ds(berv=-13.1, bjd=2454304.41):
        head = fits.Header()
        head['HIERARCH OHP DRS BJD'] = bjd
        if berv is not None:
            head['HIERARCH OHP DRS BERV'] = berv
        hdul = fits.HDUList([fits.PrimaryHDU(np.ones((2, 5)), head),
                             fits.ImageHDU(np.full((2, 5), 2.0),
                                           name='BLAZE_A')])
        hdul[0].name = 'E2DS_A'
        out = io.BytesIO()
        hdul.writeto(out)
        return out.getvalue()
    import io
    answers = {'1': e2ds(), '2': e2ds(berv=None, bjd=2460547.0), '3': b''}
    monkeypatch.setattr(lblrun.urllib.request, 'urlopen', lambda request,
                        timeout: io.BytesIO(answers[request.full_url.split(
                            ',')[-1].split('%5D')[0]]))
    sci, cal = tmp_path / 'sci', tmp_path / 'cal'
    sci.mkdir()
    cal.mkdir()
    assert lblrun._get_sophie(dict(id='1', name='1_e2ds_A.fits'), str(sci),
                              str(cal)) == ''
    head = fits.getheader(sci / '1_e2ds_A.fits')
    assert head['HIERARCH OHP DRS BLAZE FILE'] == '1_blaze_A.fits'
    assert fits.getdata(cal / '1_blaze_A.fits')[0, 0] == 2.0
    assert 'protected' in lblrun._get_sophie(dict(id='2', name='2.fits'),
                                             str(sci), str(cal))
    assert lblrun._get_sophie(dict(id='3', name='3.fits'), str(sci),
                              str(cal)) == 'not public yet'


def test_a_tar_read_by_parts(tmp_path, monkeypatch):
    """the files of a tar from their headers alone, a long name too, and
    one of them fetched by its bytes"""
    import io
    import tarfile
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode='w', format=tarfile.GNU_FORMAT) as tar:
        for name, data in (('./r.NIRPS.x_S2D_BLAZE_A.fits', b'a' * 700),
                           ('./' + 'long' * 40 + '_BLAZE_A.fits', b'b' * 10)):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    blob = raw.getvalue()

    def urlopen(request, timeout):
        first, last = request.headers['Range'][6:].split('-')
        return io.BytesIO(blob[int(first):int(last) + 1])
    monkeypatch.setattr(lblrun.urllib.request, 'urlopen', urlopen)
    members = lblrun.tar_members('http://x/tar')
    assert [(name, size) for name, _, size in members] == [
        ('r.NIRPS.x_S2D_BLAZE_A.fits', 700),
        ('long' * 40 + '_BLAZE_A.fits', 10)]
    name, start, size = members[1]
    assert lblrun._range('http://x/tar', start, size) == b'b' * 10


def test_the_runs_of_eso(monkeypatch):
    """NIRPS: HE and HA apart, the spectra corrected for the tellurics with
    their blaze; ESPRESSO: a run per mode and binning, ESPRESSO18 and 19"""
    found = [dict(id=f'ADP.{it}', target='X', mjd=58000.0 + 500 * it,
                  mode=mode, bin=binning, object='X')
             for it, (mode, binning) in enumerate(
                 [('HE', '11'), ('HE', '11'), ('HA', '11')])]
    monkeypatch.setattr(lblrun, 'eso_spectra', lambda star, ins: found)
    ident = dict(ra=1.0, dec=2.0, main='X', aliases=[], astrometry=None)
    runs = lblrun.groups('X', 'NIRPS', ident=ident)
    assert [(run['key'], run['lbl'], len(run['spectra'])) for run in runs] \
        == [('', 'NIRPS_HE', 2), ('HA', 'NIRPS_HA', 1)]
    assert runs[0]['params'] == dict(BLAZE_CORRECTED=False)
    assert runs[0]['tellu'] is False and runs[0]['input'] == '*TCORR_A.fits'
    assert runs[1]['eras'][0][1] == 'NIRPS_HA (LBL)'
    found[:] = [dict(id=f'ADP.{it}', target='X', mjd=58000.0 + 500 * it,
                     mode=mode, bin=binning, object='X')
                for it, (mode, binning) in enumerate(
                    [('SINGLEHR', '21'), ('SINGLEHR', '21'),
                     ('SINGLEUHR', '11')])]
    runs = lblrun.groups('X', 'ESPRESSO', ident=ident)
    assert [(run['key'], run['obj'], len(run['spectra'])) for run in runs] \
        == [('', 'X', 2), ('UHR11', 'XUHR11', 1)]
    assert [era[1] for era in runs[1]['eras']] == [
        'ESPRESSO18_UHR11 (LBL)', 'ESPRESSO19_UHR11 (LBL)']
    assert all(':' not in one['name'] for run in runs
               for one in run['spectra'])


def test_a_template_of_another_star(tmp_path, monkeypatch):
    """the rdb of a run with another star's template; that star's spectra
    of the same setting, spread; the BERV of the near infrared read"""
    params = lblrun.parameters('TOI-1452', dict(_group(), obj='TOI1452',
                                                tellu=False), str(tmp_path),
                               3200.0)
    params['OBJECT_COMPARISON'] = ['GJ699']
    assert lblrun.rdb_of(params).endswith('lblrdb/lbl_TOI1452_GJ699.rdb')
    many = [dict(id=str(it), name=f'{it}t.fits') for it in range(300)]
    monkeypatch.setattr(lblrun, 'groups', lambda star, inst, root='archives',
                        ident=None: [dict(_group(), lbl='SPIROU',
                                          spectra=many, obj='GJ699')])
    mine = lblrun.template_group(dict(_group(), lbl='SPIROU'), 'SPIROU',
                                 'GJ 699')
    assert len(mine['spectra']) == 60 and mine['spectra'][1]['id'] == '5'
    assert lblrun.template_group(dict(_group(), lbl='NIRPS_HA'), 'NIRPS',
                                 'GJ 699') is None
    # the span of the BERV of SPIRou's spectra
    from astropy.io import fits
    import numpy as np
    group = dict(_group(), lbl='SPIROU', source='CADC', obj='X')
    science = lblrun.folders(group, str(tmp_path))[1]
    for it, berv in enumerate((1.7, 4.8, 3.0)):
        head = fits.Header()
        head['BERV'] = berv
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(
            np.zeros((2, 2)), head, name='FluxAB')]).writeto(
                os.path.join(science, f'{it}t.fits'))
    assert abs(lblrun.berv_span(group, str(tmp_path)) - 3.1) < 1e-9
    assert lblrun.berv_span(_group(), str(tmp_path)) is None


def test_another_star_pointed_nearby():
    """an exposure named as another star of a catalogue the star is in is
    left out (the companion of a binary); the observers' own names are not"""
    barnard = dict(names=['NAME Barnard\'s Star', 'GJ 699', 'HIP 87937'])
    names = ['GJ699-HA', 'GJ699-HE', 'Gl699', 'K4_02', 'BARNARD', 'Gl 15A']
    kept = lblrun._named([dict(target=one) for one in names], barnard, 'X')
    assert [one['target'] for one in kept] == names[:-1]
    toi = dict(names=['TOI-1452', 'TIC 420112589'])
    assert [one['target'] for one in lblrun._named(
        [dict(target='TOI-1452'), dict(target='TOI-1760')], toi, 'X')] == [
        'TOI-1452']
    gj15a = dict(names=['GJ 15 A', 'HD 1326'])
    assert [one['target'] for one in lblrun._named(
        [dict(target=one) for one in ('GJ15A', 'GJ15B', 'Gl 15', 'HD1326')],
        gj15a, 'X')] == ['GJ15A', 'Gl 15', 'HD1326']

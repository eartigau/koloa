#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa.observatories, without the network: the velocities of the pipelines
from the archives of NExScI (NEID), the OHP (SOPHIE, ELODIE), ESO (HARPS,
ESPRESSO, NIRPS) and the CFHT (SPIRou); the headers kept in a table; and
a file believed before the velocities of a CCF in the near infrared

Created on 2026-10-13

@author: artigau
"""
import json

import numpy as np
import pytest

from koloa import observatories as obs
from koloa import published
from koloa.data import RVData

#: a star that does not move, and one that does (Barnard's)
STAR = dict(ra=100.0, dec=20.0, pmra=0.0, pmdec=0.0)
BARNARD = dict(ra=269.45207696, dec=4.69336497, pmra=-801.551,
               pmdec=10362.394)


def test_neid(monkeypatch):
    asked = []
    rows = [dict(object='TIC 1', ccfjdsum=str(2459300.5 + it), rv=str(
        -110243.5 + it), err='1.6', obsmode='hr', swversion='v1.5.3',
        rejected='0', qrad='100.0', qdecd='20.0') for it in range(5)]
    rows += [dict(rows[0], ccfjdsum='2459320.5', obsmode='he'),
             # rejected by the pipeline; pointed elsewhere; another star
             #   pointed at the same place, by its velocity
             dict(rows[0], ccfjdsum='2459321.5', rejected='1'),
             dict(rows[0], ccfjdsum='2459322.5', qdecd='20.02'),
             dict(rows[0], ccfjdsum='2459323.5', object='TIC 2',
                  rv='-108000.0')]
    monkeypatch.setattr(obs, '_tap', lambda url, query, where, fmt='csv':
                        asked.append((url, query)) or rows)
    parts, words = obs.neid(STAR, 'nowhere', name='a star')
    assert asked[0][0] == obs.NEID_TAP and 'neidl2' in asked[0][1]
    assert words == 'NEID (NExScI): 5; NEID-HE (NExScI): 1'
    data, about = parts[0]
    assert data.n == 5 and data.time[0] == pytest.approx(59300.5)
    assert data.rv[1] - data.rv[0] == pytest.approx(1.0)
    assert np.all(data.err == 1.6)
    # its pipeline is told how the star moves
    assert about['perspective'] is False and about['defer'] is False
    assert '1 as TIC 2, of another star' in about['left']
    assert published.family(about['label']) == 'NEID'
    monkeypatch.setattr(obs, '_tap', lambda *a, **k: [])
    assert obs.neid(STAR, 'nowhere') == (
        [], 'NEID has no public velocity of the star')


def test_the_archives_of_the_ohp(monkeypatch):
    urls = []
    ccf = [dict(seq=str(it), objname='HD1', bjd=str(2455000.5 + 300 * it),
                rv=str(-33.18 + 1e-3 * it), err='0.0010', mask='G2',
                fwhm='7.4', contrast='37.0', span='-0.012')
           for it in range(6)]
    ccf += [dict(ccf[0], mask='K5'),            # the same spectrum, twice
            # a protected programme: no date, a rounded velocity
            dict(seq='20', objname='HD1', bjd='', rv='-33.2', err='0.1',
                 mask='G2')]
    modes = [dict(seq=str(it), mode='HR') for it in range(6)]
    modes[5]['mode'] = 'HE'

    def table(url, where):
        urls.append(url)
        if 'n=sophiecc' in url:
            return ccf
        if 'n=sophie&' in url:
            return modes
        return [dict(objname='HD1', jdb=str(50000.5 + it), imatyp=kind,
                     masque='R37000K0', vfit=str(4.768 + 1e-3 * it),
                     sn='111.0', sigfit='4.0') for it, kind in enumerate(
            ['OBTH', 'OBTHs', 'OBTHd', 'OBJOd'])]
    monkeypatch.setattr(obs, '_table', table)
    parts, words = obs.sophie(STAR, 'nowhere')
    # before and after the fibres of June 2011, and its two resolutions
    assert words == 'SOPHIE (OHP): 3; SOPHIE+ (OHP): 2; SOPHIE+ HE (OHP): 1'
    assert 'p=J06%2040%2000.0%20%2B20%2000%2000' in urls[0]
    data, about = parts[0]
    assert data.time[0] == pytest.approx(55000.5)
    assert data.err[0] == pytest.approx(1.0)
    assert data.rv[1] - data.rv[0] == pytest.approx(1.0)
    assert about['perspective'] and '1 with another mask than G2' in \
        about['left']
    # the indicators of its CCF, as koloa keeps those of DACE: m/s and
    #   per cent, the errors of the pipelines
    assert set(data.indicators) == {'fwhm', 'bis', 'contrast'}
    assert data.indicators['fwhm'][0][0] == pytest.approx(7400.0)
    assert data.indicators['fwhm'][1][0] == pytest.approx(2.0)
    assert data.indicators['bis'][0][0] == pytest.approx(-12.0)
    assert data.indicators['contrast'][1][0] == pytest.approx(
        37.0 * 2.0 / 7400.0)
    assert {published.family(one[1]['label']) for one in parts} == {
        'SOPHIE'}
    assert len({one[1]['key'] for one in parts}) == 3
    # ELODIE: the exposures with the lamp, and the error of Perrier et
    #   al. (2003): 6 m/s and the photon noise, 2.8 m/s at 70 Vir's S/N
    parts, words = obs.elodie(STAR, 'nowhere')
    data, about = parts[0]
    assert words == 'ELODIE (OHP): 3' and data.time[0] == 50000.5
    assert data.err[0] == pytest.approx(np.hypot(6.0, 2.8), abs=0.01)
    assert '1 with no lamp beside the star' in about['left']
    # the width of its CCF alone, from the sigma fitted
    assert list(data.indicators) == ['fwhm']
    assert data.indicators['fwhm'][0][0] == pytest.approx(9419.3, abs=0.1)
    assert published.family(about['label']) == 'ELODIE'


def test_the_headers_kept_in_a_table(tmp_path):
    read = []

    def reader(one):
        read.append(one['file'])
        return dict(error='' if one['file'] != 'c' else 'no', rv=1.0)
    listed = [dict(file=name, public=name != 'd') for name in 'abcd']
    path = str(tmp_path / 'sub' / 'table.csv')
    rows = obs._headers(listed, path, reader, 'a test')
    assert read == ['a', 'b', 'c'] and [row['file'] for row in rows] == [
        'a', 'b']
    # nothing fetched twice, but the one that failed; a file made public
    listed[3]['public'] = True
    rows = obs._headers(listed, path, reader, 'a test')
    assert read == ['a', 'b', 'c', 'c', 'd'] and len(rows) == 3
    assert len(obs._headers(listed, path, reader, 'a test', True)) == 3
    assert len(read) == 9


def _eso_rows(star):
    cat = [dict(tel_object='Gl699', ins_mode=mode, drs_bjd=bjd,
                drs_ccf_rvc=-110.28 + 1e-3 * it, drs_ccf_noise=0.0005,
                drs_mask='M2', ra_simbad=star['ra'], dec_simbad=star['dec'])
           for it, (mode, bjd) in enumerate([
               ('HAM', 2454194.9), ('HAM', 2454195.9), ('EGGS', 2454196.9),
               ('HAM', 2458000.9), ('HAM', 2458001.9)])]
    # a star 20 arcsec away in the catalogue: not it
    cat.append(dict(cat[0], dec_simbad=star['dec'] + 20 / 3600))
    files = []
    for it, (ins, mjd, release, orig) in enumerate([
            ('ESPRESSO', 58600.0, '2020', 'ES_SCCA_1.fits'),
            ('ESPRESSO', 58700.0, '2020', 'ESPRESSO_CCF_A_2.fits'),
            ('ESPRESSO', 58701.0, '2020', 'ESPRESSO_CCF_TELL_CORR_A_2.fits'),
            ('ESPRESSO', 58702.0, '2999', 'ESPRESSO_CCF_A_3.fits'),
            ('ESPRESSO', 60091.0, '2024', 'ESPRESSO_CCF_A_4.fits'),
            ('NIRPS', 60069.0, '2024', 'r.NIRPS.1_CCF_TELL_CORR_A.fits'),
            ('NIRPS', 60069.0, '2024', 'r.NIRPS.1_CCF_A.fits'),
            ('HARPS', 61308.0, '2026', 'r.HARPS.1_CCF_A.fits'),
            ('HARPS', 61308.0, '2026', 'r.HARPS.1_CCF_B.fits')]):
        years = (mjd - 51544.5) / 365.25
        files.append(dict(
            ins=ins, target='GJ699', t_min=mjd, release=release, orig=orig,
            id=f'ADP.{it}', dec=star['dec'] + years * star['pmdec'] / 3.6e6,
            ra=star['ra'] + years * star['pmra'] / 3.6e6
            / np.cos(np.radians(star['dec']))))
    return cat, files


def test_eso(tmp_path, monkeypatch):
    cat, files = _eso_rows(BARNARD)
    monkeypatch.setattr(obs, '_tap', lambda url, query, where, fmt='csv':
                        cat if url == obs.ESO_CAT else files)
    read = []

    def header(one):
        read.append(one['id'])
        mjd = next(row['t_min'] for row in files if row['id'] == one['id'])
        return dict(error='', ins=one['ins'], target=one['target'],
                    bjd=2400000.5 + mjd, rv=-110.18, err=0.0002,
                    mask='M5', mode='SINGLEHR' if one['ins'] == 'ESPRESSO'
                    else 'HE', version='espdr/3.0.0' if mjd > 60000
                    and one['ins'] == 'ESPRESSO' else 'espdr/2.2.1', qc=1)
    monkeypatch.setattr(obs, '_eso_read', header)
    parts, words = obs.eso(BARNARD, str(tmp_path), name='GJ 699')
    got = {about['label']: (data.n, about['defer'], about['perspective'])
           for data, about in parts}
    # HARPS from the catalogue: its high-accuracy mode, before and after
    #   2015; ESPRESSO before and after 2019; one CCF an exposure; the
    #   star's fibre; what is public
    assert got == {'HARPS03 (ESO)': (2, False, True),
                   'HARPS15 (ESO)': (2, False, True),
                   'ESPRESSO18 (ESO)': (1, False, True),
                   'ESPRESSO19 (ESO)': (1, False, True),
                   # a new pipeline, another zero point
                   'ESPRESSO19 (ESO v3)': (1, False, True),
                   'NIRPS (ESO)': (1, True, True),
                   'HARPS-espdr (ESO)': (1, False, True)}
    assert read == ['ADP.0', 'ADP.1', 'ADP.4', 'ADP.5', 'ADP.7']
    assert words.endswith('; 1 files not public yet')
    harps = parts[0][0]
    assert harps.time[0] == pytest.approx(54194.9)
    assert harps.err[0] == pytest.approx(0.5)
    assert (tmp_path / 'eso_exposures.csv').exists()
    # a file of NIRPS is given: its CCFs are not fetched
    read.clear()
    parts, _ = obs.eso(dict(BARNARD, have=['NIRPS']), str(tmp_path / 'two'))
    assert 'ADP.5' not in read
    assert 'NIRPS (ESO)' not in [about['label'] for _, about in parts]


def test_spirou_at_the_cfht(tmp_path, monkeypatch):
    rows = [dict(target='Gl699', ra=BARNARD['ra'] + 23.0 * BARNARD['pmra']
                 / 3.6e6, dec=BARNARD['dec'] + 23.0 * BARNARD['pmdec']
                 / 3.6e6, mjd=59945.0 + it, uri=f'cadc:CFHT/{it}v.fits',
                 release='2020-01-01T00:00:00.000' if it < 3
                 else '2999-01-01T00:00:00.000') for it in range(4)]
    # the place of 2000, on a night of 2023: not where the star was
    rows.append(dict(rows[0], ra=BARNARD['ra'], dec=BARNARD['dec'] + 0.01,
                     uri='cadc:CFHT/9v.fits'))
    monkeypatch.setattr(obs, '_tap', lambda *a, **k: rows)
    monkeypatch.setattr(obs, '_cadc_read', lambda one: dict(
        error='', target=one['target'], bjd=2459945.5 + int(
            one['file'][0]), rv=-110.1297, err=4.5, mask='GL699_neg.fits',
        version='0.7.295', qc='T', dpr='OBJ_FP'))
    parts, words = obs.cfht(BARNARD, str(tmp_path))
    data, about = parts[0]
    assert words == 'SPIRou (CFHT): 3; 1 files not public yet'
    assert data.n == 3 and data.err[0] == 4.5
    assert data.time[0] == pytest.approx(59945.5)
    # APERO is told how the star moves; a CCF in the near infrared
    #   defers to a file
    assert about['perspective'] is False and about['defer'] is True
    assert published.family(about['label']) == 'SPIRou'


def test_a_file_believed_first(tmp_path, monkeypatch):
    """the velocities of a CCF in the near infrared are left out when a
    file of that spectrograph is given, and the archive of the CFHT is
    not asked when one is known to be"""
    def cfht(star, folder, refresh=False, name=None):
        time = 59000.5 + np.arange(4.0)
        return [(RVData(time, np.zeros(4), np.full(4, 4.5),
                        inst=np.full(4, 'SPIRou (CFHT)')),
                 dict(key='spirou_cfht', label='SPIRou (CFHT)', left='',
                      reference='CFHT archive', perspective=False,
                      defer=True))], 'SPIRou (CFHT): 4'
    monkeypatch.setattr(obs, 'SOURCES', dict(obs.SOURCES, cfht=cfht))
    ident = dict(main='GJ 699', aliases=[], ra=1.0, dec=2.0)
    folder = str(tmp_path / 'pub')
    # slow, and of little use: asked only when it is named
    assert 'cfht' in published.OPTIONAL
    assert [note for note in published.fetch(
        ident, str(tmp_path / 'usual'),
        sources=None, papers=False) if note['key'] == 'spirou_cfht'] == []
    notes = published.fetch(ident, folder, sources=['cfht'])
    assert notes[0]['defer'] and not notes[0]['perspective']
    assert notes[0]['where'] == 'CFHT archive'
    assert published.deferring(folder) == ['SPIRou (CFHT)']
    assert published.holding(folder) == []
    data = published.load(folder)
    # with a file of another spectrograph: kept; of SPIRou: left out
    kept, gone = published.left_to_files(data, folder, ['NIRPS'])
    assert kept is data and gone == []
    kept, gone = published.left_to_files(data, folder, ['SPIROU'])
    assert kept is None and gone == ['SPIRou (CFHT)']
    assert published.left_to_files(None, folder, ['SPIROU']) == (None, [])
    # known when the archives are gathered: not asked at all
    notes = published.fetch(ident, str(tmp_path / 'two'), sources=['cfht'],
                            have=['SPIROU'])
    assert notes[0]['n'] == 0 and notes[0]['note'] == (
        'not asked: a file of SPIRou is given')
    # on the page: the file of SPIRou, and not the archive's
    from koloa import gui
    from koloa.gather import write_rv
    star = tmp_path / 'arch' / 'GJ_699'
    (star / 'rv').mkdir(parents=True)
    (star / 'rv' / 'published').mkdir()
    for name in ('published.json', 'spirou_cfht.csv'):
        (star / 'rv' / 'published' / name).write_text(
            open(f'{folder}/{name}').read())
    mine = str(tmp_path / 'lbl.csv')
    write_rv(RVData(59100.5 + np.arange(5.0), np.zeros(5), np.ones(5),
                    inst=np.full(5, 'SPIROU')), mine)
    whole, source, notes, _, _ = gui.whole_series(
        mine, 'GJ 699', str(tmp_path / 'arch'), vizier=True)
    assert whole.instruments == ['SPIROU']
    assert any(note.startswith('SPIRou (CFHT): left out') for note in notes)
    # without the file, the archive's
    whole, source, notes, _, _ = gui.whole_series(
        '', 'GJ 699', str(tmp_path / 'arch'), vizier=True)
    assert whole.instruments == ['SPIRou (CFHT)']
    assert source['SPIRou (CFHT)'] == 'CFHT archive'


def test_the_table_of_the_archives(tmp_path, monkeypatch, capsys):
    """a gather tells koloa's GUI where each archive is, a row each (the
    papers of the star in one), out of its log; and the page reads the
    same rows from what a star gathered earlier left"""
    from koloa import gui, tng
    from koloa import log as klog
    from koloa.gather import space_state, state_of
    # said only to the GUI
    monkeypatch.delenv('KOLOA_PROGRESS', raising=False)
    klog.row('dace', 'asking')
    assert capsys.readouterr().out == ''
    monkeypatch.setenv('KOLOA_PROGRESS', 'gui')
    klog.row('tng', 'asking', 'two\nlines', 3, 10)
    assert capsys.readouterr().out.strip().endswith(
        '| row: tng | asking | 3 | 10 | two lines')
    # the files of an archive as they come: its row, and no line
    monkeypatch.setattr(klog, 'VERBOSE', True)
    assert len(list(tng._steps(iter(range(5)), 5, 'TNG', 'tng'))) == 5
    said = capsys.readouterr().out.splitlines()
    assert said and all('row: tng | asking' in line for line in said)
    assert said[-1].rstrip().endswith('| row: tng | asking | 5 | 5 |')
    # a run keeps the rows apart from its log
    job = gui.Job.__new__(gui.Job)
    job.rows = {}
    assert job._row('row: tng | asking | 3 | 10 | ') is True
    assert job._row('row: neid | ok |  |  | 215 velocities: NEID 215')
    assert job._row('published, NEID: 215 kept') is False
    assert job.rows == dict(
        tng=dict(state='asking', text='', done=3, total=10),
        neid=dict(state='ok', text='215 velocities: NEID 215', done=None,
                  total=None))
    # what a gather asks, from its command line
    usual = [key for key in published.SOURCES if key != 'cfht']
    assert gui.asked_archives(['x', '--gather', 'a']) == (
        ['dace', 'carmenes'] + usual + ['tess', 'space'])
    assert gui.asked_archives(['x', '--gather', 'a', '--no-dace',
                               '--vizier-sources', 'tng', 'cfht',
                               '--no-tess']) == ['carmenes', 'tng', 'cfht']
    assert gui.asked_archives(['x', '--gather', 'a', '--no-vizier',
                               '--no-carmenes']) == ['dace', 'tess', 'space']
    # the states: of an entry of the manifest, of the notes of a source
    assert state_of(dict(status='ok', npoints=3, instruments=dict(A=3))
                    ) == ('ok', '3 velocities: A 3')
    assert state_of(dict(status='ok', sectors=[dict(npoints=5)] * 2)
                    ) == ('ok', '2 sectors, 10 points')
    assert state_of(dict(status='unreachable', message='down')) == (
        'failed', 'down')
    assert state_of(None) == ('none', '')
    assert space_state(dict(kepler=dict(status='none', message='x'))
                       )[0] == 'none'
    notes = [
        dict(kind='survey', key='teklu25', file='t.csv', n=4,
             instruments={'HIRES (Teklu+ 2025)': 4}),
        dict(kind='survey', key='cls21', file=None, n=0,
             note='not in the survey'),
        dict(kind='paper', key='a', file='a.csv', n=7, instruments={'X': 7}),
        dict(kind='paper', key='b', file='b.csv', n=3, instruments={'Y': 3}),
        dict(kind='archive', key='tng_drs37', source='tng', file='x.csv',
             n=5, instruments={'HARPS-N (TNG 3.7)': 5},
             note='5 exposures; 2 files not public yet'),
        dict(kind='archive', key='tng_espdr', source='tng', file='y.csv',
             n=6, instruments={'HARPS-N (TNG espdr)': 6}, note='6'),
        dict(kind='archive', key='neid', source='neid', file=None, n=0,
             note='the NEID archive could not be reached (down)'),
        dict(kind='archive', key='cfht', source='cfht', file=None, n=0,
             note='not asked: a file of SPIRou is given')]
    assert published.papers_state(notes, 40) == (
        'ok', '10 velocities in 2 papers of the 40 of the star in SIMBAD')
    assert published.papers_state([])[0] == 'none'
    assert published.archive_state(notes[4:6]) == (
        'ok', '11 velocities: HARPS-N (TNG 3.7) 5, HARPS-N (TNG espdr) 6; '
        '2 files not public yet')
    assert published.archive_state(notes[6:7])[0] == 'failed'
    assert published.archive_state(notes[7:])[0] == 'off'
    # the table of a star gathered earlier
    star = tmp_path / 'arch' / 'GJ_1'
    (star / 'rv' / 'published').mkdir(parents=True)
    assert gui.archive_rows('GJ 1', str(tmp_path / 'arch')) == []
    (star / 'manifest.json').write_text(json.dumps(dict(archives=dict(
        dace=dict(status='ok', npoints=9, instruments=dict(HARPS03=9)),
        carmenes=dict(status='none', message='not in CARMENES DR1'),
        tess=dict(status='unreachable', message='MAST is down')))))
    (star / 'rv' / 'published' / 'published.json').write_text(
        json.dumps(notes))
    (star / 'rv' / 'published' / 'sources.json').write_text(json.dumps(
        ['teklu25', 'cls21', 'papers', 'tng', 'neid', 'cfht']))
    rows = {row['key']: row for row in gui.archive_rows(
        'GJ 1', str(tmp_path / 'arch'))}
    assert list(rows) == (['dace', 'carmenes'] + list(published.SOURCES)
                          + ['tess', 'space'])
    assert {key: row['state'] for key, row in rows.items()} == dict(
        dace='ok', carmenes='none', teklu25='ok', cls21='none',
        talor19='off', fischer14='off', rvbank20='off', papers='ok',
        tng='ok', eso='off', neid='failed', sophie='off', elodie='off',
        cfht='off', tess='failed', space='off')
    assert rows['papers']['text'] == '10 velocities in 2 papers'
    assert rows['tng']['text'].startswith('11 velocities')
    assert rows['dace']['name'] == 'DACE'
    # while a run gathers it: where each archive is
    job.action, job.returncode = 'gather', None
    job.args, job.outputs = ['GJ 1', '--gather', 'arch'], str(star)
    monkeypatch.setattr(gui, 'JOBS', dict(a=job))
    rows = {row['key']: row for row in gui.archive_rows(
        'GJ 1', str(tmp_path / 'arch'))}
    assert rows['tng'] == dict(key='tng', name=published.LABELS['tng'],
                               state='asking', text='', done=3, total=10)
    assert rows['neid']['state'] == 'ok' and rows['dace']['state'] == 'waiting'
    assert rows['cfht']['state'] == 'off'

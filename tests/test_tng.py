#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa.tng, without the network: the header of a CCF file, a star found by
where it was (and not its companion), the files of a star, their headers
kept in a table, the two reductions as two datasets; then the archive as
a source of published velocities, and the perspective acceleration taken
out of the velocities of the pipelines where a series is put together

Created on 2026-10-13

@author: artigau
"""
import gzip
import json

import numpy as np
import pytest

from koloa import published, secular, tng
from koloa.data import RVData

#: Barnard's star (J2000, Gaia DR3) and its perspective acceleration
BARNARD = dict(ra=269.45207696, dec=4.69336497, pmra=-801.551,
               pmdec=10362.394, parallax=546.976)


def _card(key, value, comment=''):
    """a card of a FITS header"""
    if isinstance(value, str):
        value = f"'{value:8s}'"
    elif isinstance(value, bool):
        value = 'T' if value else 'F'
    head = (f'HIERARCH {key} = {value}' if ' ' in key or len(key) > 8
            else f'{key:8s}= {value!s:>20s}')
    return f'{head} / {comment}'[:80].ljust(80)


def _file(cards):
    """a gzipped FITS file of one header"""
    text = ''.join(_card(*card) for card in cards) + 'END'.ljust(80)
    return gzip.compress(text.ljust(2880).encode('ascii'))


def test_the_header_of_a_ccf_file():
    raw = _file([('SIMPLE', True), ('NAXIS', 0), ('EXPTIME', 300.0),
                 ('TNG QC CCF RV', -110.303515536763, 'Radial velocity'),
                 ('TNG QC CCF RV ERROR', 0.0034593, 'its error [km/s]'),
                 ('TNG QC BJD', 2456476.72238674), ('TNG QC CCF MASK', 'M2'),
                 ('TNG DPR TYPE', 'STAR WAVE'),
                 ('ESO PRO REC1 PIPE ID', 'espdr/3.2.0')])
    head = tng.header(raw)
    assert head['SIMPLE'] is True and head['NAXIS'] == 0
    assert head['TNG QC CCF RV'] == pytest.approx(-110.303515536763)
    assert head['TNG QC BJD'] == pytest.approx(2456476.72238674)
    assert head['TNG QC CCF MASK'] == 'M2'
    assert head['ESO PRO REC1 PIPE ID'] == 'espdr/3.2.0'
    # not gzipped: the same
    assert tng.header(gzip.decompress(raw)) == head
    # the reduction of a file, by its name
    assert tng.reduction_of('HARPN.2016-03-04T05-13-47.777_ccf_M2_A.fits.gz'
                            ) == '3.7'
    assert tng.reduction_of('r.HARPN.2013-07-03T05-09-26.181_CCF_A.fits.gz'
                            ) == 'espdr'
    for name in ('r.HARPN.2013-07-03T05-09-26.181_CCF_RESIDUALS_A.fits.gz',
                 'HARPN.2016-03-04T05-13-47.777_ccf_M2_A.tbl.gz',
                 'r.HARPN.2013-07-03T05-09-26.181_S2D_A.fits.gz'):
        assert tng.reduction_of(name) is None


def test_a_star_found_by_where_it_was():
    star = BARNARD
    # where it was in 2000 (the catalogue of most), and 20 years later:
    #   40 arcsec from where it was at the latest Gaia release
    assert tng.away(star['ra'], star['dec'], star) == pytest.approx(0.0)
    later = star['dec'] + 20 * star['pmdec'] / 3.6e6
    assert tng.away(star['ra'], later, star) > 40
    assert float(tng.away(star['ra'], later, star, tng.EPOCHS + (2020.0,))
                 ) < 20
    # that place is on its path since HARPS-N began, not off it
    assert tng.offpath(star['ra'], later, star, end=2026.0) < 20
    assert tng.offpath(star['ra'] + 0.02, later, star, end=2026.0) > 60
    # a companion 15 arcsec north of a star that moves east lies on the
    #   path of neither: its own place is what counts
    pair = dict(ra=10.0, dec=40.0, pmra=2900.0, pmdec=0.0)
    north = 40.0 + 15.0 / 3600.0
    assert tng.away(10.0, north, pair) == pytest.approx(15.0, abs=0.01)
    # a star whose motion is not known: its place of 2000
    still = dict(ra=10.0, dec=40.0)
    assert tng.away(10.0, 40.0, still) == 0.0
    assert tng.offpath(10.0, 40.001, still) == pytest.approx(3.6, abs=0.01)


def test_the_pointings_of_a_star(monkeypatch):
    monkeypatch.setattr(tng, 'pointings', lambda refresh=False: dict(
        asked=0.0, object=['GJ699', 'K4_02', 'other', 'still private'],
        ra=np.array([269.45207696, 269.45207, 100.0, 269.45207]),
        dec=np.array([4.69336497, 4.69336, 4.7, 4.69337]),
        free=np.array([10, 20, 5, 0]), private=np.array([0, 0, 0, 7])))
    got = tng.near(BARNARD)
    assert got == dict(free=True, private=True, nfree=30, nprivate=7,
                       objects=['GJ699', 'K4_02', 'still private'])
    assert tng.near(dict(ra=0.0, dec=0.0)) == dict(
        free=False, private=False, nfree=0, nprivate=0, objects=[])


def _rows(star, dates, name='GJ699', policy='FREE', shift=0.0, years=0.0):
    """the rows of the archive for exposures of a star: a CCF file of
    each reduction, and a file that is not one; pointed at where the star
    was in 2000, or some years later, or off it [deg]"""
    out = []
    star = dict(star, ra=star['ra'] + years * star['pmra'] / 3.6e6
                / np.cos(np.radians(star['dec'])),
                dec=star['dec'] + years * star['pmdec'] / 3.6e6)
    for date in dates:
        for file in (f'HARPN.{date}_ccf_M2_A.fits.gz',
                     f'r.HARPN.{date}_CCF_A.fits.gz',
                     f'r.HARPN.{date}_CCF_RESIDUALS_A.fits.gz'):
            out.append(dict(
                OBJECT=name, policy=policy, DATE_OBS=date.replace('-', ':', 0),
                RA_RAD=str(np.radians(star['ra'])),
                DEC_RAD=str(np.radians(star['dec'] + shift)),
                file_url=f'http://archive/files/{file}'))
    return out


def test_the_files_of_a_star(monkeypatch):
    asked = []
    rows = (_rows(BARNARD, ['2016-03-04T05-13-47.777'])
            + _rows(BARNARD, ['2019-08-04T22-00-00.000'], 'K4_02', 'PRIV')
            # pointed 30 arcsec away: another star
            + _rows(BARNARD, ['2018-08-02T22-00-00.000'], 'far', shift=30 / 3600)
            # at where the star was that night, 18 years on
            + _rows(BARNARD, ['2018-08-03T22-00-00.000'], 'of the night',
                    years=18.6))
    monkeypatch.setattr(tng, '_ask', lambda query, **k: asked.append(query)
                        or rows)
    got = tng.files(BARNARD)
    assert len(asked) == 1 and "INSTRUMENT = 'HARPN'" in asked[0]
    assert [(one['object'], one['reduction'], one['policy']) for one in got
            ] == [('GJ699', '3.7', 'FREE'), ('GJ699', 'espdr', 'FREE'),
                  ('of the night', '3.7', 'FREE'),
                  ('of the night', 'espdr', 'FREE'),
                  ('K4_02', '3.7', 'PRIV'), ('K4_02', 'espdr', 'PRIV')]


def _exposure(file, reduction, bjd, rv, err=0.001, name='GJ699', mask='M2',
              policy='FREE', dpr='STAR,WAVE', qc=1):
    """a row of the table of the exposures of a star"""
    row = {key: '' for key in tng.COLUMNS}
    row.update(file=file, reduction=reduction, policy=policy, object=name,
               date='2016-03-04', bjd=bjd, rv=rv, err=err, mask=mask,
               version='espdr/3.2.0' if reduction == 'espdr' else 'HARPN_3.7',
               dpr=dpr, qc=qc if reduction == 'espdr' else '')
    return row


def test_the_headers_kept_in_a_table(tmp_path, monkeypatch):
    listed = [dict(file=f'r.HARPN.{it}_CCF_A.fits.gz', url=f'http://a/{it}',
                   reduction='espdr', policy='FREE' if it < 3 else 'PRIV',
                   object='GJ699', date=f'2016-03-0{it + 1}T00:00:00')
              for it in range(4)]
    monkeypatch.setattr(tng, 'files', lambda star: listed)
    fetched = []

    def fetch(one, timeout=120.0):
        fetched.append(one['file'])
        return _exposure(one['file'], 'espdr', 2457451.5 + len(fetched),
                         -110.3)
    monkeypatch.setattr(tng, '_fetch', fetch)
    rows = tng.exposures(BARNARD, str(tmp_path))
    assert len(fetched) == 3 and len(rows) == 4
    # the private file is listed, with nothing read
    assert rows[3]['policy'] == 'PRIV' and rows[3]['rv'] == ''
    # nothing is fetched twice; a file made public since is
    listed[3]['policy'] = 'FREE'
    rows = tng.exposures(BARNARD, str(tmp_path))
    assert len(fetched) == 4 and all(row['rv'] != '' for row in rows)
    assert (tmp_path / 'tng_exposures.csv').exists()
    rows = tng.exposures(BARNARD, str(tmp_path))
    assert len(fetched) == 4
    rows = tng.exposures(BARNARD, str(tmp_path), refresh=True)
    assert len(fetched) == 8


def test_two_reductions_two_datasets():
    rows = []
    for it in range(6):
        rows.append(_exposure(f'a{it}', '3.7', 2457451.5 + it,
                              -110.2465 + 1e-3 * it, 0.0008))
        rows.append(_exposure(f'b{it}', 'espdr', 2457451.5008 + it,
                              -110.3035 + 1e-3 * it, 0.0010))
    # another mask (another zero point), no velocity, a lamp, a failed
    #   check, a velocity far from the others, a private file
    rows += [_exposure('c', 'espdr', 2457460.5, -110.30, mask='K5'),
             _exposure('d', 'espdr', 2457461.5, 'nan'),
             _exposure('e', 'espdr', 2457462.5, -110.30, dpr='WAVE,WAVE'),
             _exposure('f', 'espdr', 2457463.5, -110.30, qc=0),
             _exposure('g', 'espdr', 2457464.5, -90.0),
             _exposure('h', 'espdr', '', '', policy='PRIV')]
    # the companion, pointed at where the star was to be: its velocity
    #   says it is another star
    rows += [_exposure(f'i{it}', 'espdr', 2457470.5 + it, -109.2,
                       name='companion') for it in range(3)]
    parts = tng.series(rows, 'GJ 699')
    assert [note['label'] for _, note in parts] == [
        'HARPS-N (TNG 3.7)', 'HARPS-N (TNG espdr)']
    old, new = parts[0][0], parts[1][0]
    assert old.n == 6 and new.n == 6
    # m/s and BJD - 2400000
    assert old.time[0] == pytest.approx(57451.5)
    assert old.rv[1] - old.rv[0] == pytest.approx(1.0)
    assert new.err[0] == pytest.approx(1.0)
    assert published.family(new.instruments[0]) == 'HARPS-N'
    note = parts[1][1]
    assert note['mask'] == 'M2' and note['version'] == 'espdr/3.2.0'
    assert note['objects'] == ['GJ699']
    for words in ('3 with no velocity', '1 with another mask than M2',
                  '3 as companion, of another star', '1 more than 2 km/s'):
        assert words in note['left']
    assert tng.series([_exposure('h', 'espdr', '', '', policy='PRIV')]) == []


def test_the_indicators_of_a_ccf(tmp_path):
    """the width, the contrast and the bisector span of the CCF come with
    the velocities, as koloa keeps those of DACE, and are kept on disk"""
    rows = []
    for it in range(4):
        row = _exposure(f'b{it}', 'espdr', 2457451.5 + it, -110.3, 0.0035)
        row.update(fwhm=3.1327, contrast=39.512, bis=-0.00382 + 1e-3 * it)
        rows.append(row)
        # the older reduction has no bisector span in its header
        row = _exposure(f'a{it}', '3.7', 2457451.5 + it, -110.25, 0.0008)
        row.update(fwhm=2.9713, contrast=29.516)
        rows.append(row)
    old, new = (part[0] for part in tng.series(rows))
    assert set(new.indicators) == {'fwhm', 'contrast', 'bis'}
    assert set(old.indicators) == {'fwhm', 'contrast'}
    # the errors the pipeline gives: 6.9 m/s and 0.087 per cent here
    assert new.indicators['fwhm'][0][0] == pytest.approx(3132.7)
    assert new.indicators['fwhm'][1][0] == pytest.approx(7.0)
    assert new.indicators['bis'][0][1] == pytest.approx(-2.82)
    assert new.indicators['contrast'][1][0] == pytest.approx(0.0883,
                                                             abs=1e-3)
    # a span of half the width is a measure that failed: blank
    got = tng.ccf_indicators(np.ones(2), fwhm=[3.1, 3.1], bis=[1.8, 0.004])
    assert np.isnan(got['bis'][0][0]) and got['bis'][0][1] == pytest.approx(
        4.0)
    assert np.isnan(got['bis'][1][0]) and got['bis'][1][1] == 2.0
    assert tng.ccf_indicators(np.ones(3)) == {}
    assert tng.ccf_indicators(np.ones(2), fwhm=[np.nan, np.nan]) == {}
    # written with the velocities, and read back
    path = str(tmp_path / 'tng_espdr.csv')
    published._write(new, path)
    back = RVData.from_csv(path, inst='inst')
    assert set(back.indicators) == {'fwhm', 'contrast', 'bis'}
    assert back.indicators['bis'][0][1] == pytest.approx(-2.82)
    # a series without indicators: as before
    plain = RVData(np.arange(3.0), np.zeros(3), np.ones(3),
                   inst=np.full(3, 'x'))
    published._write(plain, path)
    assert open(path).readline().strip() == 'rjd,vrad,svrad,inst'


def test_how_far_the_files_are(monkeypatch, capsys):
    """a bar in a terminal; a line every tenth of the files where the
    output is read line by line (the page, a log), not one at each step
    of a bar"""
    from koloa import log as klog
    # koloa silent: nothing at all
    monkeypatch.setattr(klog, 'VERBOSE', False)
    assert tng._progress(range(5), 5, 'a test') == range(5)
    monkeypatch.setattr(klog, 'VERBOSE', True)
    monkeypatch.setattr(tng.sys.stderr, 'isatty', lambda: False,
                        raising=False)
    assert list(tng._progress(iter(range(374)), 374, 'HARPS-N, TNG archive')
                ) == list(range(374))
    lines = [line for line in capsys.readouterr().out.splitlines()
             if 'HARPS-N, TNG archive' in line]
    assert len(lines) == 10
    assert '38 of 374 files, 0 s so far, about 0 s left' in lines[0]
    assert lines[-1].rstrip('\x1b[0m').endswith('374 of 374 files, 0 s so far')
    # a few files: a line each
    list(tng._progress(iter(range(3)), 3, 'a test'))
    assert len(capsys.readouterr().out.splitlines()) == 3


def test_the_velocities_of_a_star(tmp_path, monkeypatch):
    # never pointed at: nothing is asked beyond the list
    monkeypatch.setattr(tng, 'files', lambda star: 1 / 0)
    assert tng.velocities(BARNARD, str(tmp_path)) == (
        [], 'HARPS-N was not pointed at the star')
    monkeypatch.setattr(tng, 'near', lambda star, refresh=False: dict(
        free=False, private=True, objects=['M108']))
    assert 'nothing public yet' in tng.velocities(BARNARD, str(tmp_path))[1]
    monkeypatch.setattr(tng, 'near', lambda star, refresh=False: dict(
        free=True, private=False, objects=['GJ699']))
    monkeypatch.setattr(tng, 'exposures', lambda star, folder, refresh, told=None: [
        _exposure(f'b{it}', 'espdr', 2457451.5 + it, -110.3)
        for it in range(4)] + [_exposure('p', 'espdr', '', '', policy='PRIV')])
    parts, words = tng.velocities(BARNARD, str(tmp_path), name='GJ 699')
    assert len(parts) == 1 and parts[0][0].n == 4
    assert words == ('4 exposures by the pipeline espdr/3.2.0 (mask M2), as '
                     'GJ699; 1 files not public yet')
    # pointed near it, at another star
    monkeypatch.setattr(tng, 'exposures', lambda star, folder, refresh, told=None: [])
    assert 'not at it' in tng.velocities(BARNARD, str(tmp_path))[1]


def _two(star, folder, refresh=False, name=None):
    """the two reductions of four exposures"""
    time = 57451.5 + np.arange(4.0) * 400.0
    parts = [(RVData(time + 0.0008 * (red == 'espdr'), np.arange(4.0),
                     np.ones(4), inst=np.full(4, about['label']), name=name),
              dict(about, reduction=red, n=4, mask='M2', left=''))
             for red, about in tng.REDUCTIONS.items()]
    return parts, '4 exposures twice, as GJ699'


def test_the_archive_among_the_published_velocities(tmp_path, monkeypatch):
    """the archive of the TNG is a source like a survey: each reduction a
    dataset kept whole, two releases of the same spectra; a folder
    gathered before it was read has it added without VizieR asked again"""
    monkeypatch.setattr(published, 'CACHE', str(tmp_path / 'cache'))
    asked = []

    def survey(item, ra, dec, timeout=60.0):
        asked.append(item['key'])
        return None, 'not in the survey'
    monkeypatch.setattr(published, 'survey_velocities', survey)
    monkeypatch.setattr(tng, 'velocities', _two)
    ident = dict(main='GJ 699', aliases=[], **{key: BARNARD[key]
                                              for key in ('ra', 'dec')})
    folder = str(tmp_path / 'pub')
    assert 'tng' in published.SOURCES
    # as a folder of before: the surveys only
    old = [survey['key'] for survey in published.SURVEYS]
    notes = published.fetch(ident, folder, sources=old)
    assert len(asked) == 5 and not any(note['file'] for note in notes)
    notes = published.fetch(ident, folder, sources=old + ['tng'])
    assert len(asked) == 5
    got = [note for note in notes if note['kind'] == 'archive']
    assert [note['key'] for note in got] == ['tng_drs37', 'tng_espdr']
    assert all(note['perspective'] and note['whole'] and note['n'] == 4
               and note['where'] == 'TNG archive' for note in got)
    assert published.holding(folder) == ['HARPS-N (TNG 3.7)',
                                         'HARPS-N (TNG espdr)']
    assert set(published.shown(folder).values()) == {'TNG archive'}
    # two sources, one spectrograph: they may hold the same spectra
    came = published.origins(folder)
    assert len(set(came.values())) == 2
    data = published.load(folder)
    from koloa import datasets
    assert list(datasets.links(data, came).values()) == [
        published.SAME_FAMILY]
    # read back, nothing asked
    monkeypatch.setattr(tng, 'velocities', lambda *a, **k: 1 / 0)
    assert published.fetch(ident, folder, sources=old + ['tng']) == notes
    # an archive that cannot be reached: said, the rest goes on
    def down(*args, **kwargs):
        raise RuntimeError('the TNG archive could not be reached (down)')
    monkeypatch.setattr(tng, 'velocities', down)
    notes = published.fetch(ident, str(tmp_path / 'other'), sources=['tng'])
    assert notes[0]['n'] == 0 and 'could not be reached' in notes[0]['note']
    # a paper whose velocities hold the perspective acceleration too
    (tmp_path / 'paper').mkdir()
    (tmp_path / 'paper' / 'published.json').write_text(json.dumps([dict(
        kind='paper', key='x', reference='Mignon+ 2024', file='x.csv',
        catalogue='J/A+A/689/A32', instruments={'HARPS (Mignon+ 2024)': 3}),
        dict(kind='survey', key='y', reference='RVBank', file='y.csv',
             catalogue='J/A+A/636/A74', instruments={'HARPS03 (RVBank)': 3})]))
    assert published.holding(str(tmp_path / 'paper')) == [
        'HARPS (Mignon+ 2024)']


def test_the_astrometry_of_a_star(tmp_path, monkeypatch):
    monkeypatch.undo()
    monkeypatch.setattr(secular, 'ASTROMETRY_CACHE', str(tmp_path))
    monkeypatch.setattr(secular, '_FAILED', {})
    asked = []

    def tap(url, query, timeout):
        asked.append(url)
        if url == secular.GAIA_TAP:
            row = {key: '0.1' for key in secular.GAIA_COLUMNS}
            row.update(source_id='4472832130942575872', pmra='-801.551',
                       pmdec='10362.394', parallax='546.976')
            return [row]
        return [dict(pmra='-802.8', pmdec='10362.5', plx_value='547.0')]
    monkeypatch.setattr(secular, '_tap', tap)
    ident = dict(main="NAME Barnard's star",
                 gaia_dr3='Gaia DR3 4472832130942575872')
    astro = secular.star_astrometry(ident)
    assert astro['source'] == 'Gaia DR3 4472832130942575872'
    assert astro['pmdec'] == pytest.approx(10362.394)
    # kept: not asked again
    assert secular.star_astrometry(ident) == astro and len(asked) == 1
    acc = secular.star_perspective(ident)
    assert acc['value'] == pytest.approx(4.539, abs=0.002)
    assert acc['source'].startswith('Gaia DR3')
    # the star carries its astrometry (the folder of a gathered star)
    assert secular.star_perspective(dict(astrometry=dict(
        pmra=-801.551, pmdec=10362.394, parallax=546.976)))['value'] == \
        pytest.approx(4.539, abs=0.002)
    # Gaia stalls: SIMBAD's, and nothing kept
    def stalled(url, query, timeout):
        if url == secular.GAIA_TAP:
            raise TimeoutError('The read operation timed out')
        return tap(url, query, timeout)
    monkeypatch.setattr(secular, '_tap', stalled)
    other = dict(main='GJ 1', gaia_dr3='Gaia DR3 12')
    assert secular.star_astrometry(other)['source'] == 'SIMBAD'
    assert not (tmp_path / '12.json').exists()
    # nothing can be asked: none, and not asked again at once
    def down(url, query, timeout):
        asked.append('down')
        raise OSError('no network')
    monkeypatch.setattr(secular, '_tap', down)
    assert secular.star_astrometry(dict(main='GJ 2')) is None
    count = len(asked)
    assert secular.star_astrometry(dict(main='GJ 2')) is None
    assert len(asked) == count
    assert secular.star_perspective(dict(main='GJ 2')) is None
    assert secular.star_astrometry(dict()) is None


def test_the_perspective_drift_out_of_the_pipelines(tmp_path):
    """the velocities of the pipelines (DACE, the archive of an
    observatory) hold the perspective acceleration of the star, the
    others do not: taken out of the first where the page, a batch or a
    report put the series together"""
    time = secular.PERSPECTIVE_TREF + 365.25 * np.arange(6.0)
    data = RVData(np.tile(time, 2), np.zeros(12), np.ones(12),
                  inst=np.array(['HARPS03'] * 6 + ['CARMENES'] * 6))
    freed = secular.without_perspective(data, ['HARPS03'], 4.5)
    assert np.allclose(freed.rv[freed.inst == 'HARPS03'],
                       -4.5 * np.arange(6.0))
    assert np.all(freed.rv[freed.inst == 'CARMENES'] == 0.0)
    assert secular.without_perspective(data, ['NIRPS'], 4.5) is data
    # a gathered star, as the page reads it
    from koloa import gui
    from koloa.gather import stale, write_rv
    star = tmp_path / 'arch' / 'GJ_699'
    write_rv(data, str(star / 'rv' / 'all_rv.csv'))
    pub = star / 'rv' / 'published'
    pub.mkdir(parents=True)
    lines = ''.join(f'{val + 0.3},0.0,1.0,{name}\n' for name in (
        'HARPS-N (TNG espdr)', 'HIRES (Teklu+ 2025)') for val in time)
    (pub / 'both.csv').write_text('rjd,vrad,svrad,inst\n' + lines)
    (pub / 'published.json').write_text(json.dumps([
        dict(kind='archive', key='tng_espdr', where='TNG archive',
             reference='TNG archive (HARPS-N, pipeline espdr)',
             file='both.csv', n=6, perspective=True,
             instruments={'HARPS-N (TNG espdr)': 6}),
        dict(kind='survey', key='teklu25', reference='Teklu', n=6,
             file=None, instruments={'HIRES (Teklu+ 2025)': 6})]))
    (star / 'target.json').write_text(json.dumps(dict(main='GJ 699')))
    kept = dict(dace=True, carmenes=True, vizier=True)
    root = str(tmp_path / 'arch')
    # no astrometry of the star: as they are, and said
    whole, source, notes, _, _ = gui.whole_series('', 'GJ 699', root, **kept)
    assert np.all(whole.rv == 0.0)
    assert any('not taken out' in note for note in notes)
    assert stale(str(star))
    (star / 'target.json').write_text(json.dumps(dict(
        main='GJ 699', astrometry=dict(source='Gaia DR3 1', **{
            key: BARNARD[key] for key in ('pmra', 'pmdec', 'parallax')}))))
    assert not stale(str(star))
    whole, source, notes, _, _ = gui.whole_series('', 'GJ 699', root, **kept)
    slope = {name: np.polyfit((whole.time[whole.inst == name]
                               - secular.PERSPECTIVE_TREF) / 365.25,
                              whole.rv[whole.inst == name], 1)[0]
             for name in whole.instruments}
    assert slope['HARPS03'] == pytest.approx(-4.539, abs=0.002)
    assert slope['HARPS-N (TNG espdr)'] == pytest.approx(-4.539, abs=0.002)
    assert slope['CARMENES'] == 0.0 and slope['HIRES (Teklu+ 2025)'] == 0.0
    assert source['HARPS-N (TNG espdr)'] == 'TNG archive'
    assert source['HIRES (Teklu+ 2025)'] == 'VizieR'
    assert any(note.startswith('perspective acceleration of the star (4.54 '
                               'm/s/yr, Gaia DR3 1) taken out of HARPS03, '
                               'HARPS-N (TNG espdr)') for note in notes)


def test_what_koloas_users_shared(tmp_path, monkeypatch):
    """the shared table of a star read before the archive: when it has
    every public file the pointings count at the star, the archive is not
    asked; when it has fewer, only the files it lacks are fetched; a table
    written and read back as it was; the star found by where it is"""
    from koloa import shared
    rows = [_exposure(f'r.HARPN.2016-03-0{it + 1}T01-02-03.456_CCF_A.fits.gz',
                      'espdr', 2457451.5 + it, -110.3) for it in range(3)]
    text = shared.write_rdb(rows)
    assert text.splitlines()[0].split('\t') == shared.FIELDS
    back = shared.read_rdb(text)
    assert [row['file'] for row in back] == [row['file'] for row in rows]
    assert back[0]['date'] == '2016-03-01T01:02:03.456'
    assert back[0]['policy'] == 'FREE' and float(back[1]['rv']) == -110.3
    # the star of the index: by where it is (within 2 arcsec)
    entry = dict(id='J175748.5+044137', name='GJ 699', ra=BARNARD['ra'],
                 dec=BARNARD['dec'] + 1.0 / 3600, path='stars/x.rdb',
                 pointed=3)
    monkeypatch.setattr(shared, 'index', lambda refresh=False: [entry])
    assert shared.entry(BARNARD) == entry
    assert shared.entry(dict(BARNARD, dec=BARNARD['dec'] + 3 / 3600)) is None
    monkeypatch.setattr(shared, '_get', lambda path, timeout=20.0:
                        text.encode() if path == 'stars/x.rdb' else None)
    # the pointings count three public files there: the archive not asked
    monkeypatch.setattr(tng, 'pointings', lambda refresh=False: dict(
        asked=0.0, object=['GJ699'], ra=np.array([BARNARD['ra']]),
        dec=np.array([BARNARD['dec']]), free=np.array([3]),
        private=np.array([0])))
    monkeypatch.setattr(tng, 'files', lambda star: 1 / 0)
    told = {}
    got = tng.exposures(BARNARD, str(tmp_path), told=told)
    assert len(got) == 3 and told == dict(shared=3, fetched=0, asked=False)
    parts, words = tng.velocities(BARNARD, str(tmp_path / 'v'), name='GJ 699')
    assert parts[0][0].n == 3 and 'the archive not asked' in words
    # a fourth since: the archive asked, the fourth alone fetched
    monkeypatch.setattr(tng, 'pointings', lambda refresh=False: dict(
        asked=0.0, object=['GJ699'], ra=np.array([BARNARD['ra']]),
        dec=np.array([BARNARD['dec']]), free=np.array([4]),
        private=np.array([0])))
    listed = [dict(file=row['file'], url='http://a/' + row['file'],
                   reduction='espdr', policy='FREE', object='GJ699',
                   date=shared.date_of(row['file'])) for row in rows]
    listed.append(dict(file='r.HARPN.2016-03-09T01-02-03.456_CCF_A.fits.gz',
                       url='http://a/4', reduction='espdr', policy='FREE',
                       object='GJ699', date='2016-03-09T01:02:03.456'))
    monkeypatch.setattr(tng, 'files', lambda star: listed)
    fetched = []

    def fetch(one, timeout=120.0):
        fetched.append(one['file'])
        return _exposure(one['file'], 'espdr', 2457459.5, -110.3)
    monkeypatch.setattr(tng, '_fetch', fetch)
    told = {}
    got = tng.exposures(BARNARD, str(tmp_path / 'w'), told=told)
    assert len(got) == 4 and fetched == [listed[-1]['file']]
    assert told == dict(shared=3, fetched=1, asked=True)
    # refresh: what was shared not read, every file fetched
    fetched.clear()
    tng.exposures(BARNARD, str(tmp_path / 'w'), refresh=True)
    assert len(fetched) == 4


def test_a_star_asked_to_be_shared(monkeypatch):
    """only when sharing is on, gh logged in, and the star not among the
    names never to share: an issue that says where the star is"""
    from koloa import shared
    calls = []

    def gh(*args, timeout=60.0):
        calls.append(args)

        class Done:
            returncode, stdout, stderr = 0, b'', b''
        return Done()
    monkeypatch.setattr(shared, '_gh', gh)
    star = dict(BARNARD, names=['NAME Barnard\'s Star', 'GJ 699', 'Gl 699'])
    # off by default
    assert shared.request(star, 'GJ 699', 5) == '' and not calls
    shared.save_settings(True, ['* alf Cen', 'gl699'])
    assert shared.settings()['share'] is True
    assert shared.request(star, 'GJ 699', 5).startswith('not shared')
    assert not calls
    shared.save_settings(True, ['alf Cen'])
    # gh not logged in: nothing more
    assert shared.request(star, 'GJ 699', 5) == '' and not calls
    monkeypatch.setattr(shared, 'can_share', lambda: True)
    assert shared.request(star, 'GJ 699', 5) == 'asked to be shared'
    args = calls[-1]
    assert args[:4] == ('issue', 'create', '--repo', shared.REPO)
    body = json.loads(args[args.index('--body') + 1])
    assert body['ra'] == BARNARD['ra'] and body['fetched'] == 5
    assert set(body) == {'name', 'ra', 'dec', 'pmra', 'pmdec', 'fetched',
                         'koloa'}
    # KOLOA_SHARE=0 says no, whatever the file says
    monkeypatch.setenv('KOLOA_SHARE', '0')
    assert shared.request(star, 'GJ 699', 5) == ''

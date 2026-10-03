#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa.published, without the network: a survey's star by position, the
same spectrum published twice, an instrument column that is not one, the
strict choice of a paper's table, a star's published velocities kept and
read back, and the page's series with them

Created on 2026-10-03

@author: artigau
"""
import json
import xml.etree.ElementTree as ElementTree

import numpy as np
import pytest

from koloa import literature, published
from koloa.data import RVData


def test_a_star_in_a_survey_list():
    stars = [('HIP57087', 175.54622, 26.70657), ('GJ1', 10.0, -20.0)]
    assert published.nearest_star(stars, 175.5463, 26.7066) == 'HIP57087'
    # 30 arcsec away at most
    assert published.nearest_star(stars, 175.5463, 26.7166) is None
    assert published.nearest_star([], 0.0, 0.0) is None


def test_the_spectrograph_of_an_instrument():
    assert published.family('HIRES-k (CLS)') == 'HIRES'
    assert published.family('HARPS03 (RVBank)') == 'HARPS'
    assert published.family('HARPS-N (Mortier+ 2020)') == 'HARPS-N'
    assert published.family('Mayor+ 2009') == ''


def test_the_same_spectrum_published_twice():
    kept = [(np.array([100.0, 200.0]), np.array(['HIRES', 'HIRES']))]
    data = RVData(np.array([100.0 + 3 / 1440, 150.0, 200.0 + 3 / 1440,
                            300.0]), np.zeros(4), np.ones(4),
                  inst=np.array(['HIRES (CLS)'] * 3 + ['APF (CLS)']))
    # HIRES 3 minutes apart: the same spectra; the others new
    fresh = published._new(data, kept)
    assert fresh.n == 2 and set(fresh.time) == {150.0, 300.0}
    # another spectrograph 3 minutes apart: another spectrum
    other = RVData(np.array([100.0 + 3 / 1440]), np.zeros(1), np.ones(1),
                   inst=np.array(['APF (CLS)']))
    assert published._new(other, kept).n == 1


def test_an_instrument_column_that_is_not_one():
    data = RVData(np.arange(10.0), np.zeros(10), np.ones(10),
                  inst=np.array([f'0.{70 + it}' for it in range(10)]))
    tidy = published._tidy(data, 'Lopez-Morales+ 2014')
    assert tidy.instruments == ['Lopez-Morales+ 2014']
    named = RVData(np.arange(4.0), np.zeros(4), np.ones(4),
                   inst=np.array(['HARPS', 'HARPS', 'HIRES', 'HIRES']))
    assert set(published._tidy(named, 'X+ 2020').instruments) == {
        'HARPS (X+ 2020)', 'HIRES (X+ 2020)'}


def _votable(described: str):
    """a VizieR answer: one catalogue, one table of velocities, no column
    of names"""
    text = f'''<VOTABLE><RESOURCE><DESCRIPTION>{described}</DESCRIPTION>
<TABLE name="J/A+A/1/1/table1"><DESCRIPTION>Radial velocities</DESCRIPTION>
<FIELD name="BJD" ucd="time.epoch" unit="d"/>
<FIELD name="RV" ucd="spect.dopplerVeloc" unit="m/s"/>
<FIELD name="e_RV" ucd="stat.error" unit="m/s"/>
<DATA><TABLEDATA><TR><TD>2459000.5</TD><TD>1.0</TD><TD>1.0</TD></TR>
<TR><TD>2459001.5</TD><TD>2.0</TD><TD>1.0</TD></TR></TABLEDATA></DATA>
</TABLE></RESOURCE></VOTABLE>'''
    return ElementTree.fromstring(text)


def test_a_paper_that_only_mentions_the_star(monkeypatch):
    monkeypatch.setattr(literature, 'has_catalogue', lambda *a: True)
    # the catalogue of another star: not its velocities (strict)
    monkeypatch.setattr(literature, '_votable',
                        lambda *a: _votable('A super-Earth around GJ 876'))
    data, note = literature.vizier_velocities('J/A+A/1/1', ['GJ 436'],
                                              strict=True)
    assert data is None
    # the catalogue of the star: its velocities
    monkeypatch.setattr(literature, '_votable',
                        lambda *a: _votable('Radial velocities of GJ 436'))
    data, note = literature.vizier_velocities('J/A+A/1/1', ['GJ 436'],
                                              strict=True)
    assert data.n == 2
    # the papers of the known planets: the only table is the star's
    monkeypatch.setattr(literature, '_votable',
                        lambda *a: _votable('A super-Earth around GJ 876'))
    assert literature.vizier_velocities('J/A+A/1/1', ['GJ 436'])[0].n == 2


def test_a_star_fetched_and_read_back(tmp_path, monkeypatch):
    """the surveys then the papers; a spectrum twice kept once, an
    imprecise instrument left out; kept as CSVs and read back"""
    monkeypatch.setattr(published, 'CACHE', str(tmp_path / 'cache'))

    def survey(item, ra, dec, timeout=60.0):
        if item['key'] == 'teklu25':
            return RVData(np.array([100.0, 200.0]), np.array([1.0, 2.0]),
                          np.ones(2), inst=np.array(['HIRES (Teklu+ 2025)'] * 2),
                          name='teklu25'), 'HIP1: 2'
        if item['key'] == 'cls21':
            return RVData(np.array([100.0 + 2 / 1440, 300.0]), np.zeros(2),
                          np.ones(2), inst=np.array(['HIRES (CLS)'] * 2),
                          name='cls21'), 'hip1: 2'
        return None, 'not in the survey'
    monkeypatch.setattr(published, 'survey_velocities', survey)
    monkeypatch.setattr(published, 'star_bibcodes',
                        lambda main, timeout=60.0: [
                            ('2020A&A...600A..10X', 2020),
                            ('2019A&A...601A..11Y', 2019)])
    monkeypatch.setattr(published, '_describe', lambda cat, timeout=60.0:
                        dict(title=('X+ 2020' if '600' in cat else 'Y+ 2019'),
                             rv=True))

    def paper(cat, names, title, timeout=60.0, strict=False):
        assert strict
        if title == 'X+ 2020':
            return RVData(np.array([400.0, 500.0]), np.zeros(2),
                          np.ones(2) * 2.0, inst=np.array(['HARPS'] * 2),
                          name=cat), f'{cat}: table1 (2)'
        # a km/s binary survey: left out
        return RVData(np.array([600.0]), np.zeros(1), np.array([300.0]),
                      inst=np.array(['CORAVEL']), name=cat), f'{cat}: t (1)'
    monkeypatch.setattr(literature, 'vizier_velocities', paper)
    ident = dict(main='GJ 1', aliases=['HIP 1'], ra=1.0, dec=2.0)
    notes = published.fetch(ident, str(tmp_path / 'pub'))
    got = {note['key']: note['n'] for note in notes}
    assert got['teklu25'] == 2 and got['cls21'] == 1
    assert got['J_A_A_600_A10'] == 2 and got['J_A_A_601_A11'] == 0
    data = published.load(str(tmp_path / 'pub'))
    assert data.n == 5 and 'HARPS (X+ 2020)' in data.instruments
    # kept: read back, nothing asked
    asked = published.survey_velocities
    monkeypatch.setattr(published, 'survey_velocities', None)
    assert published.fetch(ident, str(tmp_path / 'pub')) == notes
    # some of the sources only: asked again, the others left out (the
    #   spectrum CLS shares with Teklu now its own)
    monkeypatch.setattr(published, 'survey_velocities', asked)
    some = published.fetch(ident, str(tmp_path / 'pub'),
                           sources=['cls21', 'rvbank20'])
    assert [note['key'] for note in some] == ['cls21', 'rvbank20']
    assert some[0]['n'] == 2
    assert published.load(str(tmp_path / 'pub')).n == 2
    # the same sources again: read back
    monkeypatch.setattr(published, 'survey_velocities', None)
    assert published.fetch(ident, str(tmp_path / 'pub'),
                           sources=['rvbank20', 'cls21']) == some
    with pytest.raises(ValueError, match='hires'):
        published.fetch(ident, str(tmp_path / 'pub'), sources=['hires'])
    # what was learnt of the catalogues, for the next star
    with open(tmp_path / 'cache' / 'catalogues.json') as handle:
        assert json.load(handle)['J/A+A/600/A10']['rv']


def test_the_page_with_the_published_velocities(tmp_path):
    """the published velocities of the archives on the page, their spectra
    the archives already have left out"""
    from koloa import gui
    from koloa.gather import write_rv
    star = tmp_path / 'arch' / 'GJ_436'
    write_rv(RVData(np.array([100.0, 200.0]), np.zeros(2), np.ones(2),
                    inst=np.array(['HARPS03'] * 2)),
             str(star / 'rv' / 'all_rv.csv'))
    pub = star / 'rv' / 'published'
    pub.mkdir(parents=True)
    (pub / 'rvbank20.csv').write_text(
        'rjd,vrad,svrad,inst\n100.0001,1.0,1.0,HARPS03 (RVBank)\n'
        '150.0,2.0,1.0,HARPS03 (RVBank)\n160.0,2.0,1.0,HARPS03 (RVBank)\n'
        '170.0,2.0,1.0,HARPS03 (RVBank)\n')
    (pub / 'teklu25.csv').write_text(
        'rjd,vrad,svrad,inst\n300.0,1.0,1.0,HIRES (Teklu+ 2025)\n'
        '310.0,2.0,1.0,HIRES (Teklu+ 2025)\n'
        '320.0,2.0,1.0,HIRES (Teklu+ 2025)\n'
        '330.0,1.0,1.0,APF (CLS)\n')
    (pub / 'published.json').write_text(json.dumps([
        dict(kind='survey', key='rvbank20', reference='RVBank', file=
             'rvbank20.csv', n=2), dict(kind='survey', key='teklu25',
                                        reference='Teklu', file='teklu25.csv',
                                        n=2)]))
    res = gui.velocities('', 'GJ 436', str(tmp_path / 'arch'), dace=True,
                         vizier=True)
    got = {inst['name']: (inst['source'], inst['n'])
           for inst in res['instruments']}
    # the spectrum DACE has left out; the APF of a single velocity too
    assert got == {'HARPS03': ('DACE', 2),
                   'HARPS03 (RVBank)': ('VizieR', 3),
                   'HIRES (Teklu+ 2025)': ('VizieR', 3)}
    res = gui.velocities('', 'GJ 436', str(tmp_path / 'arch'), dace=True)
    assert any('VizieR: 8 points gathered, not ticked' in note
               for note in res['notes'])
    assert gui.command('gather', dict(target='GJ 436', vizier=False))[-1] \
        == '--no-vizier'

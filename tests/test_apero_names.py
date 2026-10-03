#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa's APERO names: names cleaned and looked for as APERO does, the star
of a file, a copy of the database from the assets tarball (all offline)

Created on 2026-10-03

@author: artigau
"""
import io
import json
import tarfile

import pytest

from koloa import apero_names

pytest.importorskip('yaml')

ENTRIES = {
    'verified/STAR_A.yaml': """APERO_NAME: STAR_A
ORIGINAL_NAME: Star-A
SIMBAD_NAME: NAME Star A
STATUS: verified
RA:
  value: 10.5
  source: Gaia
SPT:
  value: M2V
MASS_STAR_MANN15: 0.45
ALIASES:
- GJ   1001
- V* XY Abc
- Shared 1
""",
    'pending/BD_M03_2870.yaml': """APERO_NAME: BD_M03_2870
ORIGINAL_NAME: Gl382
SIMBAD_NAME: BD-03  2870
STATUS: pending
SPT:
  value: M2V
ALIASES:
- GJ 382
- Shared 1
""",
    'rejected/TOI756.yaml': """APERO_NAME: TOI756
SIMBAD_NAME: '* alf Xyz'
STATUS: rejected
ALIASES: []
""",
    # a copy at the root of what a tier has: the tier's wins
    'STAR_A.yaml': """APERO_NAME: STAR_A
SIMBAD_NAME: wrong
ALIASES: []
""",
}


def _database(folder):
    """a small database of APERO names on disk"""
    for name, text in ENTRIES.items():
        path = folder / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return folder


def test_names_cleaned_as_apero_does():
    """+ to P, - to M, punctuation and spaces to _, upper case; the
    variants of spaces, signs and underscores"""
    assert apero_names.clean_object('Gl 382') == 'GL_382'
    assert apero_names.clean_object('BD-03  2870') == 'BDM03_2870'
    assert apero_names.clean_object('* alf Cen C') == 'ALF_CEN_C'
    assert apero_names.clean_object('2MASS J0014+3245') == \
        '2MASS_J0014P3245'
    assert apero_names.clean_object(None) == 'Null'
    assert apero_names.clean_object(' nan ') == 'Null'
    variants = apero_names.name_variants('BD-03 2870')
    assert variants[0] == 'BDM03_2870'
    assert {'BDM032870', 'BD03_2870', 'BD032870'} <= set(variants)
    assert apero_names.name_variants('') == []


def test_names_as_apero_itself_finds_them():
    """the same cleaning and variants as APERO's own (when installed)"""
    astro = pytest.importorskip('apero.core.drs_astrometrics')
    for name in ('Gl 382', 'BD-03  2870', '* alf Cen C', "Barnard's star",
                 'TOI-756.01', 'V* AN Sex', '2MASS J00140580-3245599',
                 '[RHG95] 1595', 'GJ   551', 'Proxima-HE', 'HD 189733 B'):
        assert apero_names.clean_object(name) == astro.clean_object(name)
        assert apero_names.name_variants(name) == \
            astro.name_search_variants(name)


def test_a_name_looked_up(tmp_path, monkeypatch):
    """by its APERO name, original name, SIMBAD name or an alias, in any
    variant; a verified entry first; the index kept beside the copy"""
    monkeypatch.setenv(apero_names.ENV_FOLDER, str(_database(tmp_path)))
    for name in ('STAR_A', 'Star-A', 'star a', 'GJ 1001', 'gj1001',
                 'V* XY Abc'):
        assert apero_names.lookup(name)['apero'] == 'STAR_A', name
    entry = apero_names.lookup('Gl382')
    assert entry['apero'] == 'BD_M03_2870' and entry['status'] == 'pending'
    assert entry['spt'] == 'M2V'
    # an alias two entries share: the verified one's
    assert apero_names.lookup('Shared 1')['apero'] == 'STAR_A'
    # the tier's copy, not the root's
    assert apero_names.lookup('STAR_A')['simbad'] == 'NAME Star A'
    assert apero_names.lookup('nobody') is None
    assert (tmp_path / apero_names.INDEX).exists()
    # SIMBAD's prefixes and spaces left out for the resolver
    assert apero_names.simbad_target(apero_names.lookup('STAR_A')) == \
        'Star A'
    assert apero_names.simbad_target(entry) == 'BD-03 2870'
    assert apero_names.simbad_target(apero_names.lookup('TOI756')) == \
        'alf Xyz'


def test_the_star_of_a_file(tmp_path, monkeypatch):
    """its OBJECT column, else its name (lbl_<OBJECT>_<TEMPLATE>.rdb), in
    the database; the name itself when the database lacks it"""
    monkeypatch.setenv(apero_names.ENV_FOLDER,
                       str(_database(tmp_path / 'db')))
    rdb = tmp_path / 'series.rdb'
    rdb.write_text('rjd\tvrad\tsvrad\tOBJECT\n---\t----\t-----\t------\n'
                   + ''.join(f'{60000 + i}\t1.0\t1.0\tGl382\n'
                             for i in range(5)))
    star = apero_names.star_of_file(str(rdb))
    assert star['source'] == 'OBJECT column' and star['raw'] == 'Gl382'
    assert star['apero'] == 'BD_M03_2870' and star['target'] == 'BD-03 2870'
    lbl = tmp_path / 'lbl_STAR_A_STAR_A_tc.rdb'
    lbl.write_text('rjd\tvrad\tsvrad\n---\t----\t-----\n60000\t1\t1\n')
    star = apero_names.star_of_file(str(lbl))
    assert star['source'] == 'file name' and star['apero'] == 'STAR_A'
    other = tmp_path / 'lbl2_NOBODY_NOBODY.rdb'
    other.write_text('rjd,vrad,svrad\n60000,1,1\n')
    star = apero_names.star_of_file(str(other))
    assert star['apero'] is None and star['target'] == 'NOBODY'
    assert apero_names.name_candidates('lbl_A_B_A_B.rdb')[0] == 'A_B'
    assert apero_names.name_candidates('lbl_TOI756_GL846.rdb') == ['TOI756']
    assert apero_names.name_candidates('AN_SEX_rv_clean.rdb') == ['AN_SEX']
    assert apero_names.name_candidates('plain.csv') == []


def test_a_copy_from_the_assets_tarball(tmp_path, monkeypatch):
    """the newest tarball of the server, its astrometrics streamed and the
    rest left (not even read); a new copy in place of the old"""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w:gz') as tar:
        def add(name, data):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
        add('README.md', b'assets')
        for name, text in ENTRIES.items():
            add(f'astrometrics/{name}', text.encode())
        add('astrometrics/.history/STAR_A.jsonl', b'{}')
        add('astrometrics/.name_index.json', b'{}')
        add('models/big.fits', b'x' * 100000)
    tarball = buf.getvalue()
    asked = []

    class Answer(io.BytesIO):
        """what urlopen gives"""

    def urlopen(url, timeout=None):
        asked.append(url)
        if url.endswith('/'):
            return Answer(b'<a href="1700000000_1_assets.tar.gz">old</a>'
                          b'<a href="1790000000_2_assets.tar.gz">new</a>')
        assert url.endswith('1790000000_2_assets.tar.gz')
        return Answer(tarball)
    monkeypatch.delenv(apero_names.ENV_FOLDER)
    monkeypatch.setattr(apero_names, 'FOLDER', str(tmp_path / 'copy'))
    monkeypatch.setattr(apero_names.urllib.request, 'urlopen', urlopen)
    out = apero_names.refresh('http://server/')
    assert out['tarball'] == '1790000000_2_assets.tar.gz'
    assert out['files'] == 4 and out['objects'] == 3
    assert not (tmp_path / 'copy' / '.history').exists()
    with open(tmp_path / 'copy' / 'source.json') as handle:
        assert json.load(handle)['entries'] == 4
    assert apero_names.lookup('GJ 382', fetch=False)['apero'] == \
        'BD_M03_2870'
    # again: the new copy replaces the old
    apero_names.refresh('http://server/')
    assert not (tmp_path / 'copy.old').exists()
    assert apero_names.lookup('STAR_A', fetch=False) is not None
    # the database read from a folder of one's own is not copied over
    monkeypatch.setenv(apero_names.ENV_FOLDER, str(tmp_path / 'copy'))
    with pytest.raises(ValueError, match=apero_names.ENV_FOLDER):
        apero_names.refresh('http://server/')

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The known planets of a star against a fit (koloa.archive; the queries of
the archive need the network and are not tested here).

Created on 2026-09-29

@author: artigau
"""
import numpy as np
import pytest

from koloa.archive import compare, _reference


KNOWN = dict(host='HD 1', planets=[dict(
    name='HD 1 b', P=10.0, K=3.0, K_err=0.3, reference='Old et al. 2006',
    reference_url='', solutions=[
        dict(reference='Old et al. 2006', reference_url='', date='2006-05',
             P=9.8, K=3.0, K_err=0.3),
        dict(reference='New et al. 2025', reference_url='', date='2025-12',
             P=10.02, K=2.0, K_err=0.2)])])


def test_compare_takes_the_most_recent_solution_first():
    orb = compare(dict(P=[10.01, 0.01, 0.01], K=[2.1, 0.1, 0.1]), KNOWN)
    assert orb['known']['name'] == 'HD 1 b'
    refs = [row['reference'] for row in orb['comparisons']]
    assert refs == ['New et al. 2025', 'Old et al. 2006']
    # (2.1 - 2.0) / hypot(0.1, 0.2) and (2.1 - 3.0) / hypot(0.1, 0.3)
    assert np.isclose(orb['K_z'], 0.1 / np.hypot(0.1, 0.2))
    assert np.isclose(orb['comparisons'][1]['z'], -0.9 / np.hypot(0.1, 0.3))


def test_compare_leaves_an_unknown_period_alone():
    orb = compare(dict(P=[25.0, 0.1, 0.1], K=[2.0, 0.2, 0.2]), KNOWN)
    assert 'known' not in orb and 'comparisons' not in orb


def test_reference_reads_the_archive_link():
    text, url = _reference('<a refstr=X href=https://ui.adsabs.harvard.edu/'
                           'abs/2025AJ....170..343H/abstract target=ref>'
                           'Harada et al. 2025</a>')
    assert text == 'Harada et al. 2025'
    assert url.endswith('2025AJ....170..343H/abstract')


def test_the_archive_is_kept_and_looked_up_locally(tmp_path, monkeypatch):
    """the two tables kept on disk: a star found by its Gaia DR3 number, its
    TIC or a name (any case and spacing), without the network"""
    import json as _json
    from koloa import archive
    row = {col: None for col in archive.PSCOMP_COLUMNS}
    row.update(hostname='GJ 436', gaia_dr3_id='Gaia DR3 1', tic_id='TIC 2',
               pl_name='GJ 436 b', pl_orbper=2.6439, pl_rvamp=17.1,
               st_rotp=44.0)
    sol = dict(hostname='GJ 436', pl_name='GJ 436 b', pl_orbper=2.6439,
               pl_rvamp=17.4, pl_rvamperr1=0.2, pl_rvamperr2=-0.2,
               pl_refname='<a href=x>Someone et al. 2020</a>',
               pl_pubdate='2020-01')
    (tmp_path / 'tables.json').write_text(_json.dumps(dict(
        fetched='2026-10-01 10:00', pscomppars=[row], ps=[sol], toi=[])))
    monkeypatch.setattr(archive, 'CACHE', str(tmp_path))
    monkeypatch.setattr(archive, '_TABLES', None)
    monkeypatch.setattr(archive, '_query', lambda *a, **k: 1 / 0)
    for ident in (dict(gaia_dr3='Gaia DR3 1'), dict(tic='TIC 2'),
                  dict(main='Ross 905', aliases=['gj436'])):
        assert archive.host_name(ident) == 'GJ 436'
    assert archive.host_name(dict(main='Gliese 707', aliases=[])) is None
    known = archive.known_planets(host='GJ 436')
    assert known['fetched'] == '2026-10-01 10:00'
    assert known['star']['rotation'] == 44.0
    assert known['planets'][0]['solutions'][0]['K'] == 17.4


EU_CSV = (
    'name,star_name,star_alternate_names,ra,dec,orbital_period,'
    'orbital_period_error_min,orbital_period_error_max,k,k_error_min,'
    'k_error_max,eccentricity,omega,tperi,tconj,tzero_tr,'
    'tzero_tr_error_min,tzero_tr_error_max,mass,mass_sini,radius,'
    'discovered,updated,detection_type,extra\n'
    'GJ 436 b,GJ 436,"Ross 905, LHS 310",175.5458,26.7064,2.64394,'
    '0.0001,0.0001,18.07,1.03,1.03,0.19,354.4,2452992.1,,2454221.61588,'
    '0.0002,0.0004,0.05836,,0.325,2004,2025-11-24,Radial Velocity,x\n'
    'GJ 436 c,GJ 436,,175.5458,26.7064,5.2,,,1.1,0.2,0.2,,,,,,,,,0.01,,'
    '2008,2010-01-01,Radial Velocity,x\n'
    'GJ 436 d,GJ 436,,175.5458,26.7064,,,,,,,,,,,,,,,,,2012,2012-01-01,'
    'Other,x\n'
    'Far b,Far star,,10.0,-5.0,7.5,,,3.0,,,,,,,,,,1.0,,,2020,2021-01-01,'
    'Primary Transit,x\n')


def test_exoplanet_eu_beside_the_archive(tmp_path, monkeypatch):
    """exoplanet.eu's catalogue kept once fetched; a star found in it by a
    name of its own (SIMBAD's aliases) or by its position; its values set
    beside the archive's planet of the same period, the planet it alone
    has added (one without a period left out)"""
    import io
    import json as _json
    from koloa import archive
    row = {col: None for col in archive.PSCOMP_COLUMNS}
    row.update(hostname='GJ 436', gaia_dr3_id='Gaia DR3 1', tic_id='TIC 2',
               pl_name='GJ 436 b', pl_orbper=2.6439, pl_rvamp=17.1)
    (tmp_path / 'tables.json').write_text(_json.dumps(dict(
        fetched='2026-10-01 10:00', pscomppars=[row], ps=[], toi=[])))
    asked = []

    def urlopen(url, timeout=None):
        asked.append(url)
        return io.BytesIO(EU_CSV.encode())
    monkeypatch.setattr(archive, 'CACHE', str(tmp_path))
    monkeypatch.setattr(archive, '_TABLES', None)
    monkeypatch.setattr(archive, '_EU', None)
    monkeypatch.setattr(archive, '_query', lambda *a, **k: 1 / 0)
    monkeypatch.setattr(archive.urllib.request, 'urlopen', urlopen)
    kept = archive.encyclopaedia()
    assert len(kept['planets']) == 4 and asked == [archive.ENCYCLOPAEDIA]
    assert (tmp_path / 'exoplanet_eu.json').exists()
    # read back from the disk, not asked again
    monkeypatch.setattr(archive, '_EU', None)
    assert len(archive.encyclopaedia()['planets']) == 4 and len(asked) == 1
    assert archive.encyclopaedia_fetched() == kept['fetched']
    # by a name (an alias of SIMBAD's against its alternate names; Gl is
    #   GJ), and by the position when no name is its
    ident = dict(main='Ross 905', name='Gl 436', aliases=['HIP 57087'],
                 ra=175.54622, dec=26.70657, gaia_dr3='Gaia DR3 1',
                 tic='TIC 2')
    found = archive.encyclopaedia_planets(ident)
    assert [pl['name'] for pl in found] == ['GJ 436 b', 'GJ 436 c',
                                            'GJ 436 d']
    assert found[0]['tc'] == pytest.approx(54221.61588)
    assert found[0]['mass_earth'] == pytest.approx(0.05836 * 317.828)
    assert found[0]['K_err'] == pytest.approx(1.03)
    by_place = archive.encyclopaedia_planets(dict(
        main='Nobody', name='Nobody', aliases=[], ra=10.001, dec=-5.001))
    assert [pl['name'] for pl in by_place] == ['Far b']
    assert archive.encyclopaedia_planets(dict(
        main='Nobody', name='Nobody', aliases=[], ra=50.0, dec=5.0)) == []
    # beside the archive's
    known = archive.known_planets(host='GJ 436', ident=ident)
    assert known['eu']['planets'] == 3
    both, only = known['planets']
    assert both['name'] == 'GJ 436 b' and both['K'] == 17.1
    assert both['eu']['K'] == 18.07 and 'source' not in both
    assert only['name'] == 'GJ 436 c' and only['source'] == 'exoplanet.eu'
    assert only['solutions'] == [] and only['tp_err'] is None
    # a star the archive does not have, but exoplanet.eu does
    far = archive.known_planets(ident=dict(
        main='Far star', name='Far star', aliases=[], ra=0.0, dec=0.0))
    assert far['host'] is None
    assert [pl['name'] for pl in far['planets']] == ['Far b']
    # not kept, and not waited for: pending, fetched in the background
    started = []
    monkeypatch.setattr(archive, '_EU', None)
    monkeypatch.setattr(archive, '_eu_background',
                        lambda: started.append(1))
    monkeypatch.setattr(archive, '_EU_THREAD', None)
    (tmp_path / 'exoplanet_eu.json').unlink()
    pending = archive.known_planets(host='GJ 436', ident=ident, wait=False)
    assert pending['eu'] == dict(pending=True)
    assert [pl['name'] for pl in pending['planets']] == ['GJ 436 b']
    archive._EU_THREAD.join(5)
    assert started == [1]

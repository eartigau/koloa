#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The known planets of a star against a fit (koloa.archive; the queries of
the archive need the network and are not tested here).

Created on 2026-09-29

@author: artigau
"""
import numpy as np

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
        fetched='2026-10-01 10:00', pscomppars=[row], ps=[sol])))
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

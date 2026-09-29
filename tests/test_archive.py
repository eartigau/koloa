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

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Published velocities (koloa.literature), the LaTeX report (koloa.latex) and
the GP of the activity (koloa.gpcheck), without the network.

Created on 2026-09-30

@author: artigau
"""
import numpy as np

from koloa.archive import bibcode, conjunction
from koloa.data import RVData
from koloa.latex import escape, pm, pretty, sci
from koloa.literature import (_columns, _name_pattern, add, read, to_rjd,
                              vizier_catalogue)


def test_vizier_catalogue_from_bibcode():
    assert vizier_catalogue('2023A&A...680A..28G') == 'J/A+A/680/A28'
    assert vizier_catalogue('2019MNRAS.482.4231C') == 'J/MNRAS/482/4231'
    assert vizier_catalogue('2020AJ....160..123R') == 'J/AJ/160/123'
    assert vizier_catalogue('2020ApJ...890L..12X') == 'J/ApJ/890/L12'
    # a journal whose tables VizieR does not name after the paper
    assert vizier_catalogue('2016Natur.536..437A') is None
    assert vizier_catalogue('too short') is None


def test_bibcode_and_conjunction_from_the_archive():
    url = 'https://ui.adsabs.harvard.edu/abs/2023A&A...680A..28G/abstract'
    assert bibcode(url) == '2023A&A...680A..28G'
    assert bibcode('') is None
    circular = dict(P=6.9442, P_err=0.001, e=0.0, tc=None, tp=57509.78,
                    tp_err=0.26)
    conj = conjunction(circular)
    assert conj['tc'] == 57509.78 and 'circular' in conj['source']
    # an eccentric orbit with only a time of periastron: no conjunction
    assert conjunction(dict(circular, e=0.3)) is None
    transit = dict(circular, e=0.3, tc=57500.1, tc_err=0.01)
    assert conjunction(transit)['source'] == 'transit'


def test_times_to_bjd_minus_2400000():
    assert np.allclose(to_rjd([2457500.5]), [57500.5])
    assert np.allclose(to_rjd([7500.5], 'BJD-2450000'), [57500.5])
    assert np.allclose(to_rjd([57500.0], 'MJD'), [57500.5])
    assert np.allclose(to_rjd([57500.5], 'rjd'), [57500.5])


def test_read_a_vizier_dat_and_a_csv(tmp_path):
    dat = tmp_path / 'paper_rvs.dat'
    dat.write_text('2457499.55 19.48 3.82 CARM-VIS\n'
                   '2457503.63 -3.23 2.20 CARM-VIS\n'
                   '2458540.10 1.00 2.50 IRD\n')
    data = read(str(dat))
    assert data.n == 3 and data.instruments == ['CARM-VIS', 'IRD']
    assert np.isclose(data.time.min(), 57499.55)
    csv = tmp_path / 'other.csv'
    csv.write_text('bjd,rv,e_rv\n2459000.1,1.0,1.0\n2459001.1,2.0,1.0\n'
                   '2459002.1,0.5,1.0\n')
    other = read(str(csv), inst='HARPS')
    assert other.instruments == ['HARPS']
    assert np.isclose(other.time.min(), 59000.1)


def test_add_leaves_out_the_same_spectra_and_renames_a_shared_instrument():
    ours = RVData(time=np.array([100.0, 101.0, 102.0]),
                  rv=np.array([1.0, 2.0, 3.0]), err=np.ones(3),
                  inst=np.array(['SPIRou'] * 3), name='star')
    theirs = RVData(time=np.array([100.0 + 20 / 86400, 50.0, 60.0]),
                    rv=np.array([1.0, 5.0, 6.0]), err=np.ones(3),
                    inst=np.array(['SPIRou'] * 3), name='paper')
    merged, info = add(ours, [theirs], ['Author 2025'])
    assert info[0]['same'] == 1 and info[0]['n'] == 2
    assert 'SPIRou (Author 2025)' in merged.instruments
    assert merged.n == 5


def test_the_names_of_a_star_in_a_table():
    pattern = _name_pattern(['GJ 3988', 'LHS 3262'])
    assert pattern.search('gj3988rv')
    assert pattern.search('rv timeseries of gj 3988')
    assert not pattern.search('gj39881')
    assert not pattern.search('gj724rv')


def test_the_columns_of_a_vizier_table():
    fields = [dict(name='recno', ucd='meta.record', unit='', description=''),
              dict(name='BJD', ucd='time.epoch', unit='d', description=''),
              dict(name='RV', ucd='phys.veloc;pos.heliocentric', unit='m/s',
                   description=''),
              dict(name='e_RV', ucd='stat.error', unit='m/s',
                   description=''),
              dict(name='BIS', ucd='phys.veloc', unit='m/s', description=''),
              dict(name='Inst', ucd='meta.id;instr', unit='',
                   description='')]
    cols = _columns(fields)
    assert cols['time']['name'] == 'BJD' and cols['rv']['name'] == 'RV'
    assert cols['err']['name'] == 'e_RV' and cols['inst']['name'] == 'Inst'
    assert _columns(fields[:2]) is None


def test_latex_text_and_numbers():
    assert escape('a_b & 5 % #1') == r'a\_b \& 5 \% \#1'
    assert 'ensuremath' in escape('±')
    assert pretty('K = 2 +- 1') == r'K = 2 \ensuremath{\pm} 1'
    assert sci(0.5) == '0.50'
    assert sci(1.2e-5) == r'$1.2\times10^{-5}$'
    assert pm((3.0, 0.5, 0.5)) == r'$3.00 \pm 0.50$'
    assert pm((3.0, 0.2, 0.9)) == r'$3.00^{+0.90}_{-0.20}$'


def test_gp_signals_on_a_simulation():
    from koloa.gpcheck import gp_signals
    from koloa.simulate import simulate
    sim = simulate(planets=[dict(P=7.3, K=8.0)], err=1.5, jitter=0.5,
                   seed=3, nvisits=30, per_visit=2, baseline=300)
    res = gp_signals(sim['data'], [7.3], workers=2, nsim=200)
    assert len(res['orbits']) == 1
    orb = res['orbits'][0]
    assert abs(orb['P'][0] / 7.3 - 1) < 0.01 and orb['dlnl'] > 20
    assert res['whitened']['null']['peaks'][0]['period'] > 7.0

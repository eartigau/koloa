#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Why an outlier is an outlier: koloa.outliers, and the columns RVData keeps.

Created on 2026-09-29

@author: artigau
"""
import numpy as np

from koloa import log as klog
from koloa.data import RVData
from koloa.outliers import explain

klog.VERBOSE = False


def series(seed: int = 2, nvis: int = 60, nexp: int = 3) -> RVData:
    """visits of nexp exposures with a S/N, an airmass, a key of pure noise,
    a text key, an exposure time and a DTEMP indicator, none of them off"""
    rng = np.random.default_rng(seed)
    time = np.concatenate([night + 0.01 * np.arange(nexp)
                           for night in np.sort(rng.uniform(0, 900, nvis))])
    npts = len(time)
    meta = {'EXTSN035': rng.normal(150, 10, npts),
            'AIRMASS': 1.3 + 0.1 * rng.normal(size=npts),
            'NOISEKEY': rng.normal(0, 1, npts),
            'EXPTIME': np.full(npts, 60.0),
            'DPRTYPE': np.array(['OBJ_FP'] * npts, dtype='<U12')}
    dtemp = (rng.normal(0, 1, npts), np.full(npts, 0.5))
    return RVData(time, rng.normal(0, 2, npts), np.full(npts, 2.0),
                  indicators={'DTEMP3500': dtemp}, meta=meta, name='test')


KEYS = ['svrad', 'EXTSN035', 'snr_rate', 'AIRMASS', 'NOISEKEY', 'DPRTYPE',
        'DTEMP3500']


def test_each_outlier_gets_its_own_cause():
    data = series()
    prob = np.zeros(data.n)
    # a bad visit of low S/N, a bad exposure of high airmass inside a good
    #   visit, and a bad visit with a rare DPRTYPE and a DTEMP excursion
    visit = np.where(data.seq == 10)[0]
    data.meta['EXTSN035'][visit] = 60.0
    prob[visit] = 1.0
    single = np.where(data.seq == 30)[0][1]
    data.meta['AIRMASS'][single] = 2.4
    prob[single] = 0.99
    other = np.where(data.seq == 45)[0]
    data.meta['DPRTYPE'][other] = 'OBJ_DARK'
    data.indicators['DTEMP3500'][0][other] += 6.0
    prob[other] = 0.97
    rep = explain(data, keys=KEYS, prob=prob)
    units = {unit['visit']: unit for unit in rep.units}
    assert len(rep.units) == 3

    def off(unit):
        return {key['key'] for key in unit['keys'] if key['significant']}
    # the S/N (and the flux per second it makes) of the bad visit
    assert {'EXTSN035', 'snr_rate'} == off(units[10])
    assert units[10]['kind'] == 'visit' and \
        units[10]['compare'] == 'neighbours'
    # the airmass of the single exposure, against its own visit
    assert off(units[30]) == {'AIRMASS'}
    assert units[30]['kind'] == 'exposure' and \
        units[30]['compare'] == 'visit'
    assert off(units[45]) == {'DPRTYPE', 'DTEMP3500'}
    # a key of pure noise, and the error bars, are never off
    for unit in rep.units:
        assert 'NOISEKEY' not in off(unit) and 'svrad' not in off(unit)
    # one outlier each: nothing that they share
    assert rep.shared() == []
    assert 'why' in rep.text() and rep.to_dict()['units']


def test_a_cause_that_the_outliers_share_is_found():
    rng = np.random.default_rng(5)
    data = series(seed=5, nvis=80)
    prob = np.zeros(data.n)
    for visit in (5, 17, 29, 41, 53, 65):
        members = np.where(data.seq == visit)[0]
        data.meta['EXTSN035'][members] -= rng.uniform(40, 60)
        prob[members] = 1.0
    rep = explain(data, keys=KEYS, prob=prob)
    shared = {row['key']: row for row in rep.shared()}
    assert set(shared) <= {'EXTSN035', 'snr_rate'} and 'EXTSN035' in shared
    assert shared['EXTSN035']['direction'] == 'lower'
    assert shared['EXTSN035']['nagree'] == 6


def test_a_constant_key_that_changes_is_off():
    data = series()
    data.meta['RESET_RV'] = np.zeros(data.n)
    prob = np.zeros(data.n)
    visit = np.where(data.seq == 20)[0]
    data.meta['RESET_RV'][visit] = 1.0
    prob[visit] = 1.0
    rep = explain(data, keys=['RESET_RV', 'NOISEKEY'], prob=prob)
    key = [row for row in rep.units[0]['keys'] if row['key'] == 'RESET_RV'][0]
    assert key['significant'] and np.isinf(key['z']) and key['z'] > 0


def test_auto_keys_follow_the_instrument():
    data = series()
    prob = np.zeros(data.n)
    rep = explain(data, prob=prob)
    # the columns say SPIRou: its list, the keys it has, and DTEMP
    assert rep.lists['inst'] == 'spirou'
    assert {'svrad', 'EXTSN035', 'snr_rate', 'AIRMASS', 'DPRTYPE',
            'DTEMP3500'} <= set(rep.keys['inst'])
    assert 'NOISEKEY' not in rep.keys['inst']
    assert 'TLPEH2O' in rep.missing['inst']


def test_rdb_columns_are_kept_through_selection_and_binning(tmp_path):
    path = tmp_path / 'star.rdb'
    rows = ['rjd\tvrad\tsvrad\tEXTSN035\tDPRTYPE\tBERV',
            '---\t----\t-----\t--------\t-------\t----']
    for it in range(8):
        rows.append(f'{60000 + (it // 2) + 0.01 * (it % 2)}\t{it}\t1.0\t'
                    f'{100 + it}\t{"OBJ_FP" if it else "OBJ_DARK"}\tNaN')
    path.write_text('\n'.join(rows) + '\n')
    data = RVData.from_csv(path)
    assert set(data.meta) == {'EXTSN035', 'DPRTYPE', 'BERV'}
    assert data.meta['EXTSN035'].dtype.kind == 'f'
    assert data.meta['DPRTYPE'].dtype.kind == 'U'
    # a column of nan is numbers, not text
    assert data.meta['BERV'].dtype.kind == 'f'
    sel = data.select(np.arange(data.n) >= 2)
    assert np.allclose(sel.meta['EXTSN035'], 100 + np.arange(2, 8))
    binned = data.binned()
    assert np.allclose(binned.meta['EXTSN035'], [100.5, 102.5, 104.5, 106.5])
    assert list(binned.meta['DPRTYPE']) == ['OBJ_DARK', 'OBJ_FP', 'OBJ_FP',
                                            'OBJ_FP']


def test_merge_keeps_instruments_visits_indicators_and_columns():
    from koloa.data import merge
    one = series(seed=1, nvis=10)
    one.inst[:] = 'A'
    one.zero_point = {'A': 0.0}
    rng = np.random.default_rng(3)
    time = np.sort(rng.uniform(0, 900, 12))
    two = RVData(time, 100 + rng.normal(0, 1, 12), np.full(12, 1.0),
                 inst=np.array(['B'] * 12), meta={'texp': np.full(12, 900.0)},
                 indicators={'fwhm': (np.ones(12), np.full(12, 0.1))})
    both = merge([one, two], name='both')
    assert both.n == one.n + two.n and both.instruments[0] in ('A', 'B')
    assert set(both.instruments) == {'A', 'B'}
    # the visits of each series stay whole, numbered in time order
    assert both.nseq == one.nseq + two.nseq
    first = [both.time[both.seq == it].min() for it in range(both.nseq)]
    assert np.all(np.diff(first) > 0)
    # each instrument keeps its velocities (its zero point out again)
    assert np.allclose(np.sort(both.rv[both.inst == 'B']),
                       np.sort(two.rv))
    # what a series lacks is nan, or empty text
    assert np.all(np.isnan(both.meta['texp'][both.inst == 'A']))
    assert np.all(both.meta['DPRTYPE'][both.inst == 'B'] == '')
    assert np.all(np.isnan(both.indicators['fwhm'][0][both.inst == 'A']))
    assert np.all(np.isnan(both.indicators['DTEMP3500'][0][both.inst == 'B']))


def test_detailed_analysis_runs_offline(tmp_path):
    """the whole path of koloa.detailed on a small simulation, without the
    network (no SIMBAD, archive nor DACE) and with a short FIP"""
    from koloa.detailed import detailed_analysis
    from koloa.simulate import simulate
    sim = simulate(planets=[dict(P=7.3, K=8.0)], err=1.5, jitter=0.5,
                   seed=3, nvisits=30, per_visit=2, baseline=300,
                   outliers=[dict(kind='visit', frac=0.1, amplitude=10.0)])
    # (a GP inside the FIP, koloa's default, needs a few hundred sweeps)
    out = detailed_analysis(sim['data'], outdir=str(tmp_path), archive=False,
                            dace=False, kmax=1, nsweep=300, nburn=150,
                            duck=False, gp=False, detection_map=False)
    assert out['fip_first'].pk is not None
    assert len(list(tmp_path.glob('*_report.txt'))) == 1
    assert len(list(tmp_path.glob('*_summary.json'))) == 1
    # the report in LaTeX, and in PDF where pdflatex is
    from koloa.latex import _pdflatex
    assert len(list(tmp_path.glob('*_report.tex'))) == 1
    if _pdflatex() is not None:
        pdf = list(tmp_path.glob('*_report.pdf'))
        assert len(pdf) == 1 and pdf[0].read_bytes()[:4] == b'%PDF'
    # the planet, found and fitted
    assert any(abs(orb['P'][0] / 7.3 - 1) < 0.01 for orb in out['orbits'])


def test_detection_map_of_a_series():
    """the map of koloa.detailed: injections on a grid, the signal found
    taken out, the limits per period bin"""
    from koloa.detailed import _detection_map
    from koloa.fit import RVModel
    from koloa.simulate import simulate
    sim = simulate(planets=[dict(P=7.3, K=8.0)], err=1.5, seed=3,
                   nvisits=30, per_visit=1, baseline=300)
    data = sim['data']
    fit = RVModel(data, [dict(period=7.3, period_range=(7.1, 7.5))],
                  likelihood='mixture', unit='point').fit(nstart=1,
                                                          quiet=True)
    dmap = _detection_map(data, fit, 1, 1.1, 1, 4, 1, 0.5)
    rmap = dmap['map']
    assert rmap.rate.shape == (12, 10)
    assert len(dmap['K90']) == 12 and dmap['msini90'] is not None
    # a planet of 15 errors is found at short periods
    assert np.nanmax(rmap.rate[:4, -1]) >= 0.5

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Data, periodograms, the simulator and the figures.

Created on 2026-09-27

@author: artigau
"""
import os

import numpy as np
import pytest

from koloa import log as klog
from koloa.data import RVData, find_sequences
from koloa.periodogram import frequency_grid, gls, jackknife, window
from koloa.simulate import add_outliers, simulate

klog.VERBOSE = False
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# =============================================================================
# Data
# =============================================================================
def test_reads_lbl_csv_and_finds_visits(tmp_path):
    # a real NIRPS sampling: 181 exposures in 58 visits
    data = RVData.from_csv(os.path.join(ROOT, 'data', 'nirps_template.csv'),
                           name='NIRPS template')
    assert data.n == 181 and data.nseq == 58
    binned = data.binned()
    assert binned.n == 58
    assert np.all(binned.err < np.max(data.err))
    # the columns of LBL, five exposures on two nights: every column with an
    #   error column (sX or sig_X) is an indicator, BERV has none
    table = dict(rjd=[60100.80, 60100.81, 60100.82, 60102.79, 60102.80],
                 vrad=[-33010.0, -33011.5, -33008.0, -33007.5, -33006.0],
                 svrad=[3.1, 2.7, 3.0, 3.3, 2.9],
                 DTEMP3500=[-0.41, -0.62, -0.20, 0.15, 0.33],
                 sDTEMP3500=[0.45, 0.41, 0.43, 0.47, 0.42],
                 dW=[-9.9e5, -1.0e6, -9.8e5, -9.7e5, -9.6e5],
                 sdW=[5.0e4, 4.4e4, 4.7e4, 5.1e4, 4.6e4],
                 fwhm=[10963.2, 10961.5, 10962.0, 10960.1, 10959.8],
                 sig_fwhm=[4.5, 4.0, 4.2, 4.6, 4.1],
                 BERV=[12.3, 12.3, 12.3, 11.9, 11.9])
    lines = [','.join(table)]
    lines += [','.join(str(col[it]) for col in table.values())
              for it in range(5)]
    path = tmp_path / 'lbl.csv'
    path.write_text('\n'.join(lines) + '\n')
    lbl = RVData.from_csv(path)
    assert lbl.n == 5 and lbl.nseq == 2
    assert set(lbl.indicators) == {'DTEMP3500', 'dW', 'fwhm'}
    value, error = lbl.indicators['DTEMP3500']
    assert np.allclose(value, table['DTEMP3500'])
    assert np.allclose(error, table['sDTEMP3500'])
    assert np.allclose(lbl.indicators['fwhm'][1], table['sig_fwhm'])
    # the instrument zero point is out of the velocities
    assert abs(lbl.zero_point['inst'] + 33008.0) < 1e-9
    assert abs(np.median(lbl.rv)) < 1e-9


def test_reads_lbl_rdb():
    path = os.path.join(ROOT, 'data', 'toi2120_spirou_pca2d_0-7_bias10s3.rdb')
    if not os.path.exists(path):
        pytest.skip('no TOI-2120 file')
    data = RVData.from_csv(path)
    assert data.n == 316 and data.nseq == 81


def test_sequences_split_instruments_and_nights():
    time = np.array([0.0, 0.01, 0.02, 1.0, 1.01, 0.005, 5.0])
    inst = np.array(['A', 'A', 'A', 'A', 'A', 'B', 'A'])
    seq = find_sequences(time, inst, gap=0.3)
    # A: {0, 0.01, 0.02}, {1.0, 1.01}, {5.0}; B: {0.005}
    assert len(np.unique(seq)) == 4
    assert seq[0] == seq[1] == seq[2]
    assert seq[5] != seq[0]


# =============================================================================
# Periodograms
# =============================================================================
def test_gls_matches_astropy():
    astropy = pytest.importorskip('astropy.timeseries')
    rng = np.random.default_rng(3)
    time = np.sort(rng.uniform(0, 200, 80))
    err = rng.uniform(1, 2, 80)
    value = 3 * np.sin(2 * np.pi * time / 13.1) + err * rng.normal(size=80)
    freq = frequency_grid(time, 1.5, None, 5)
    ours = gls(time, value, err, freq)
    # astropy's exact method (its default on a regular grid is an
    #   approximation)
    ref = astropy.LombScargle(time, value, err, fit_mean=True).power(
        freq, normalization='standard', method='cython')
    assert np.allclose(ours, ref, atol=1e-8)


def test_jackknife_leaves_each_unit_out():
    sim = simulate(planets=[dict(P=11.0, K=4.0)], err=1.5, seed=2,
                   nvisits=30, per_visit=2, baseline=300)
    data = sim['data']
    freq = frequency_grid(data.time, 2.0, None, 5)
    jack = jackknife(data, freq, unit='sequence')
    # leaving visit 3 out is the same as computing the GLS without it
    keep = data.seq != 3
    direct = gls(data.time[keep], data.rv[keep], data.err[keep], freq)
    assert np.allclose(jack['loo'][3], direct, atol=1e-8)
    assert np.all(jack['low'] <= jack['high'])


def test_window_is_one_at_zero_frequency():
    time = np.sort(np.random.default_rng(1).uniform(0, 100, 50))
    assert np.isclose(window(time, np.array([0.0]))[0], 1.0)


# =============================================================================
# Simulator
# =============================================================================
def test_outliers_go_both_ways():
    rng = np.random.default_rng(5)
    time = np.arange(400) * 1.0
    seq = np.arange(400) // 2
    err = np.ones(400)
    for kind in ('spike', 'visit'):
        out = add_outliers(time, seq, err, kind=kind, frac=0.2,
                           amplitude=8.0, rng=rng)
        moved = out['offset'][out['mask']]
        assert np.sum(moved > 0) > 5 and np.sum(moved < 0) > 5
    # a bad visit moves all its exposures together
    out = add_outliers(time, seq, err, kind='visit', frac=0.1, amplitude=8,
                       rng=rng)
    for visit in np.unique(seq[out['mask']]):
        vals = out['offset'][seq == visit]
        assert np.all(np.sign(vals) == np.sign(vals[0]))


def test_simulated_campaign_has_seasons():
    sim = simulate(seed=4, nvisits=60, per_visit=3, baseline=1000,
                   season=240)
    time = np.sort(sim['data'].time)
    gaps = np.diff(time)
    # a target seen 240 days a year leaves yearly gaps of about 125 days
    assert np.sum(gaps > 90) >= 2


# =============================================================================
# Figures
# =============================================================================
def test_figures_are_drawn(tmp_path):
    from koloa import plotting as kplot
    from koloa.fip import fip_single
    from koloa.fit import RVModel
    sim = simulate(planets=[dict(P=9.0, K=4.0)], err=1.5, seed=3,
                   nvisits=30, per_visit=3, baseline=400,
                   outliers=[dict(kind='visit', frac=0.1, amplitude=8)])
    data = sim['data']
    fit = RVModel(data, [dict(period=9.0)], likelihood='mixture',
                  unit='sequence').fit(quiet=True)
    freq = frequency_grid(data.time, 2.0, None, 5)
    single = fip_single(data, freq=freq)
    for style in ('paper', 'web'):
        kplot.set_style(style)
        paths = [kplot.savefig(kplot.timeseries(data, fit.outlier_prob,
                                                fit=fit),
                               str(tmp_path / f'ts_{style}')),
                 kplot.savefig(kplot.phase(fit), str(tmp_path / f'ph_{style}')),
                 kplot.savefig(kplot.sequences(data, fit.outlier_prob,
                                               fit=fit),
                               str(tmp_path / f'sq_{style}')),
                 kplot.savefig(kplot.periodograms(
                     freq, gls_power=gls(data.time, data.rv, data.err, freq),
                     window_power=window(data.time, freq),
                     fips=dict(gaussian=single)), str(tmp_path / f'pg_{style}'))]
        for path in paths:
            assert os.path.getsize(path) > 1000
            assert path.endswith('.svg' if style == 'web' else '.pdf')
    kplot.set_style('paper')


# =============================================================================
# Doppler
# =============================================================================
def test_doppler_is_relativistic_and_invertible():
    from koloa import doppler
    vel = np.array([-3.0e4, -10.0, 0.0, 1.0, 3.0e4, 2.0e5])
    ratio = doppler.velocity_to_ratio(vel)
    assert np.allclose(doppler.ratio_to_velocity(ratio), vel, atol=1e-7)
    shift = doppler.velocity_to_logshift(vel)
    assert np.allclose(np.log(ratio), shift, rtol=0, atol=1e-15)
    assert np.allclose(doppler.logshift_to_velocity(shift), vel, atol=1e-7)
    # at a 30 km/s barycentric correction the first-order formula is wrong
    #   by about (v/c)^2 / 2, i.e. 1.5 m/s
    beta = 3.0e4 / doppler.SPEED_OF_LIGHT
    err = (doppler.velocity_to_ratio(3.0e4) - (1 + beta)) * \
        doppler.SPEED_OF_LIGHT
    assert 1.4 < err < 1.6


def test_coherence_is_calibrated_with_a_visit_jitter():
    """A coherent planet with a visit jitter is not flagged as incoherent"""
    import numpy as np
    from koloa.data import RVData
    from koloa.diagnostics import coherence
    from koloa.simulate import simulate
    tpl = RVData.from_csv(os.path.join(ROOT, 'data', 'nirps_template.csv'),
                          name='NIRPS template')
    pvals = []
    for seed in range(20):
        sim = simulate(planets=[dict(P=13.7, K=4.0, e=0.0,
                                     tp=float(tpl.time[0]))],
                       template=tpl, visit_jitter=3.0, seed=seed)
        coh = coherence(sim['data'], 13.7, split='seasons')
        pvals.append(coh['p_vector'])
    # the p-values of a true null are uniform: 2 of 20 below 0.01 at most
    assert np.sum(np.array(pvals) < 0.01) <= 2


def test_coherence_of_a_sparse_series_is_nan():
    """a few visits per season: nothing to compare, and no crash"""
    from koloa.diagnostics import coherence
    rng = np.random.default_rng(3)
    tv = np.concatenate([60000 + 365 * s + rng.uniform(0, 90, 4)
                         for s in range(4)])
    time = np.concatenate([t + np.arange(4) * 0.01 for t in tv])
    data = RVData(time=time, rv=rng.normal(0, 2, time.size),
                  err=np.full(time.size, 2.0))
    coh = coherence(data, 444.0, split='seasons')
    assert coh['chunks'] == [] and np.isnan(coh['p_vector'])
    assert len(coherence(data, 444.0, split='halves')['chunks']) == 2


def test_duck_test_writes_a_pdf_report(tmp_path):
    """duck_test(..., pdf=...) writes the detailed report (offline: no
    archive)"""
    from koloa.diagnostics import duck_test
    from koloa.simulate import simulate
    sim = simulate(planets=[dict(P=6.1, K=6.0)], err=1.5, seed=5, nvisits=30,
                   per_visit=2, baseline=300)
    path = tmp_path / 'duck.pdf'
    report = duck_test(sim['data'], 6.1, gp=False, quiet=True,
                       pdf=str(path), archive=False)
    assert report.verdict
    assert path.exists() and path.read_bytes()[:4] == b'%PDF'
    assert path.stat().st_size > 20000
    # everything in a folder: the report, each figure, the text, the summary
    folder = tmp_path / 'duck'
    report.pdf(None, sim['data'], archive=False, outdir=str(folder))
    names = {item.name for item in folder.iterdir()}
    assert any(name.endswith('_duck.pdf') for name in names)
    assert any(name.endswith('_phase.pdf') for name in names)
    assert any(name.endswith('_duck.txt') for name in names)
    assert any(name.endswith('_duck.json') for name in names)


def test_duck_test_leaves_a_failed_visit_out_of_the_coherence():
    """a visit thousands of m/s off, flagged by the outlier model, does not
    reach the jackknife and the coherence (they have no outlier model)"""
    from koloa.diagnostics import duck_test
    from koloa.simulate import simulate
    sim = simulate(planets=[dict(P=6.1, K=6.0)], err=1.5, seed=5, nvisits=30,
                   per_visit=2, baseline=300)
    data = sim['data']
    rv = data.rv.copy()
    bad = data.seq == data.seq[10]
    rv[bad] -= 25000.0
    data = data.with_values(rv)
    prob = bad.astype(float)
    report = duck_test(data, 6.1, prob=prob, gp=False, quiet=True,
                       nightly=False)
    coh = report.details['coherence']['halves']
    assert all(abs(chunk['K'] - 6.0) < 3.0 for chunk in coh['chunks'])
    assert report.details['jackknife_data'].n == data.n - int(np.sum(bad))


def test_a_merge_keeps_each_instrument_on_its_own_zero_point():
    """a part cut out of a larger series still holds the zero points of
    instruments it no longer has; merged with a series that has one of
    them, it must not shift it (NIRPS of a file against DACE's NIRPS)"""
    import numpy as np
    from koloa.data import RVData, merge
    rng = np.random.default_rng(1)
    archive = RVData(np.arange(20.0), np.r_[np.full(10, -28468.0),
                                            np.full(10, -2.0)]
                     + rng.normal(0, 1, 20), np.ones(20),
                     inst=np.array(['NIRPS'] * 10 + ['CARMENES'] * 10))
    carm = archive.select(archive.inst == 'CARMENES')
    assert set(carm.zero_point) == {'CARMENES'}
    # an older cut that still holds the other's zero point
    carm.zero_point['NIRPS'] = -28468.0
    lbl = RVData(np.arange(30.0, 40.0), -28203.0 + rng.normal(0, 1, 10),
                 np.ones(10), inst=np.full(10, 'NIRPS'))
    both = merge([lbl, carm])
    for inst in both.instruments:
        assert abs(np.median(both.rv[both.inst == inst])) < 2.0


def test_a_csv_names_its_instruments(tmp_path):
    """a column named inst (or instrument) that holds names gives the
    instrument of each velocity, each its own zero point (HARPS before and
    after its change of fibres); numbers there, or no such column: one
    instrument"""
    from koloa.data import RVData
    path = tmp_path / 'two.csv'
    path.write_text('rjd,vrad,svrad,inst,OBJECT\n' + ''.join(
        f'{60000 + k},{10.0 + (15.0 if k >= 6 else 0.0)},1.0,'
        f'{"HARPS15" if k >= 6 else "HARPS03"},GJ 1\n' for k in range(10)))
    data = RVData.from_csv(str(path))
    assert data.instruments == ['HARPS03', 'HARPS15']
    assert data.zero_point['HARPS15'] - data.zero_point['HARPS03'] == 15.0
    assert 'OBJECT' in data.meta and 'inst' not in data.meta
    path.write_text('rjd,vrad,svrad,inst\n' + ''.join(
        f'{60000 + k},10.0,1.0,{k % 2}\n' for k in range(10)))
    assert RVData.from_csv(str(path)).instruments == ['inst']
    path.write_text('rjd,vrad,svrad\n' + ''.join(
        f'{60000 + k},10.0,1.0\n' for k in range(10)))
    assert RVData.from_csv(str(path)).instruments == ['inst']


def test_the_figures_of_a_report_by_instrument():
    """the series of a report with each instrument about the drift of the
    star, not about zero; its points and those of the model and of the
    fold translucent, a colour per instrument, circled when flagged"""
    from koloa import plotting as kplot
    from koloa.fit import RVModel
    rng = np.random.default_rng(4)
    parts = []
    for it, (name, t0, t1, num) in enumerate((('OLD', 51000, 54000, 80),
                                              ('NEW', 57000, 60000, 80))):
        time = np.sort(rng.uniform(t0, t1, num))
        rv = (2.0 * (time - 55000) / 365.25 + 5.0 * np.sin(
            2 * np.pi * time / 17.3) + rng.normal(0, 1.0, num) + 40 * it)
        rv[::20] += 30.0  # a few bad points
        parts.append((time, rv, np.ones(num), np.array([name] * num)))
    time, rv, err, inst = (np.concatenate([part[key] for part in parts])
                           for key in range(4))
    data = RVData(time, rv, err, inst=inst, name='toy')
    fit = RVModel(data, [dict(period=17.3)], trend=1).fit(quiet=True)
    prob = fit.outlier_prob
    assert 0 < np.sum(prob > 0.5) < 20

    def points(ax, name):
        return next(line for line in ax.lines
                    if str(line.get_label()).startswith(name + ' ('))
    fig = kplot.series_instruments(data, prob, fit=fit)
    ax = fig.axes[0]
    levels = {}
    for name in ('OLD', 'NEW'):
        line = points(ax, name)
        xval, yval = line.get_data()
        # along the trend of the fit: its median is of what the trend leaves
        drift = kplot.trend_curve(fit.model, fit.theta, xval)
        assert abs(np.median(yval - drift)) < 1.0
        assert line.get_alpha() == kplot.POINT_ALPHA
        levels[name] = np.median(yval)
    # the drift of 2 m/s/yr over the 16 years between them, not a step
    #   back to zero
    assert 25 < levels['NEW'] - levels['OLD'] < 40
    assert abs(np.median(data.rv[data.inst == 'NEW'])) < 1e-6
    # circled: the points more likely outliers than not, and only those
    rings = next(line for line in ax.lines
                 if line.get_label() == 'P(outlier) > 0.5')
    assert len(rings.get_xdata()) == int(np.sum(prob > 0.5))
    assert rings.get_markerfacecolor() == 'none'
    assert any('drift' in str(line.get_label()) for line in ax.lines)
    kplot.plt.close(fig)
    # without a fit: each instrument about its own median (that of its
    #   points not flagged), and no drift drawn
    fig = kplot.series_instruments(data, prob)
    new = data.inst == 'NEW'
    assert np.allclose(points(fig.axes[0], 'NEW').get_ydata(), data.rv[new]
                       - np.median(data.rv[new & (prob <= 0.5)]))
    assert not any('drift' in str(line.get_label())
                   for line in fig.axes[0].lines)
    kplot.plt.close(fig)
    # the model: translucent points in its three panels, the flagged circled
    fig = kplot.model_series(fit)
    for ax in fig.axes:
        marks = [line for line in ax.lines if line.get_alpha()
                 == kplot.POINT_ALPHA]
        assert len(marks) == 2
        assert any(line.get_markerfacecolor() == 'none' and len(
            line.get_xdata()) == int(np.sum(prob > 0.5)) for line in ax.lines)
    kplot.plt.close(fig)
    # the fold: by instrument when asked (no scale of probability, the
    #   instruments under the axis), by probability otherwise
    fig = kplot.phase(fit, colour='inst')
    assert len(fig.axes) == 1 and len(fig.legends) == 1
    labels = [text.get_text() for text in fig.legends[0].get_texts()]
    assert labels == ['OLD (80)', 'NEW (80)', 'P(valid) < 0.5']
    assert sum(line.get_alpha() == kplot.POINT_ALPHA
               for line in fig.axes[0].lines) == 2
    kplot.plt.close(fig)
    fig = kplot.phase(fit)
    assert len(fig.axes) == 2 and not fig.legends
    kplot.plt.close(fig)

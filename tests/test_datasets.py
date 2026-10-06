"""Which datasets of a series are used (koloa.datasets)."""
import numpy as np
import pytest

from koloa import datasets
from koloa.data import RVData, merge


def _star(time):
    """what the star does: a planet and a slope"""
    return 6.0 * np.sin(2 * np.pi * time / 17.3) + 0.002 * (time - 57000.0)


def _set(name, time, noise, rng, err=None):
    """a dataset of a star: its velocities, a noise and an offset"""
    time = np.asarray(time, dtype=float)
    rv = _star(time) + rng.normal(0.0, noise, len(time)) + rng.normal(0, 50)
    return RVData(time=time, rv=rv, err=np.full(len(time), err or noise),
                  inst=np.array([name] * len(time)), name='star')


def _rows(rows):
    return {row['name']: row for row in rows}


def test_the_better_release_of_the_same_spectra():
    """two reductions of the same spectra: the one that scatters less is
    used, the other left out, whatever their errors say"""
    rng = np.random.default_rng(1)
    time = np.sort(rng.uniform(56000, 58000, 150))
    # the errors of the noisier one are the smaller: not what decides
    data = merge([_set('HARPS03', time, 3.0, rng, err=0.5),
                  _set('HARPS (Trifonov+ 2020)', time + 2e-4, 1.0, rng,
                       err=1.5)])
    used, rows = datasets.choose(data, sources={
        'HARPS03': 'DACE', 'HARPS (Trifonov+ 2020)': 'rvbank'})
    rows = _rows(rows)
    assert used.instruments == ['HARPS (Trifonov+ 2020)']
    assert rows['HARPS03']['status'] == 'release'
    assert rows['HARPS03']['better'] == 'HARPS (Trifonov+ 2020)'
    assert rows['HARPS03']['same'] == 150 and rows['HARPS03']['used'] == 0
    assert rows['HARPS03']['by'] == 'noise'
    # the noise of each, from their difference: what was put in
    assert rows['HARPS03']['precision'] == pytest.approx(3.0, rel=0.2)
    assert rows['HARPS03']['precision_better'] == pytest.approx(1.0, abs=0.5)
    assert 'spectra of HARPS (Trifonov+ 2020)' in datasets.told(
        rows['HARPS03'])
    # what it scatters more than the other: sqrt(3^2 - 1^2)
    assert rows['HARPS03']['extra'] == pytest.approx(2.83, rel=0.25)
    assert datasets.left_out(rows.values()) == ['HARPS03']


def test_asked_back_a_release_is_preferred():
    """the release the rules leave out, asked back: its spectra are the
    ones used, the other release's left out"""
    rng = np.random.default_rng(2)
    time = np.sort(rng.uniform(56000, 58000, 120))
    data = merge([_set('HARPS03', time, 3.0, rng),
                  _set('HARPS (Trifonov+ 2020)', time, 1.0, rng)])
    used, rows = datasets.choose(data, include=['harps03'])
    assert used.instruments == ['HARPS03']
    assert _rows(rows)['HARPS (Trifonov+ 2020)']['status'] == 'release'
    # left out as asked: the other release is then the one used
    used, rows = datasets.choose(data, exclude=['HARPS (Trifonov+ 2020)'],
                                 include=['HARPS03'])
    assert used.instruments == ['HARPS03']
    assert _rows(rows)['HARPS (Trifonov+ 2020)']['status'] == 'asked'


def test_a_release_keeps_the_spectra_of_its_own():
    """the less precise release has spectra the other does not: those are
    used, each spectrum once"""
    rng = np.random.default_rng(3)
    time = np.sort(rng.uniform(56000, 58000, 200))
    data = merge([_set('HARPS03', time, 2.5, rng),
                  _set('HARPS (Trifonov+ 2020)', time[:120], 1.0, rng)])
    used, rows = datasets.choose(data)
    rows = _rows(rows)
    assert rows['HARPS03']['status'] == 'part'
    assert rows['HARPS03']['used'] == 80 and rows['HARPS03']['same'] == 120
    assert used.n == 200
    assert np.sum(used.inst == 'HARPS03') == 80
    assert np.all(used.time[used.inst == 'HARPS03'] > time[119])
    assert 'of 200 points' in datasets.told(rows['HARPS03'])


def test_two_spectrographs_at_the_same_time_are_not_one():
    """HARPS and NIRPS observe together, the two arms of CARMENES too, and
    the two eras of a file are one source: none is another's release"""
    rng = np.random.default_rng(4)
    time = np.sort(rng.uniform(59000, 60000, 60))
    data = merge([_set('HARPS15', time, 1.0, rng),
                  _set('NIRPS', time, 2.0, rng),
                  _set('CARMENES VIS', time + 500, 1.5, rng),
                  _set('CARMENES NIR (Paper+ 2021)', time + 500, 6.0, rng),
                  _set('VIS (Other+ 2019)', time + 900, 1.5, rng),
                  _set('NIR (Other+ 2019)', time + 900, 6.0, rng)])
    used, rows = datasets.choose(data, auto=False)
    assert used.n == data.n
    assert all(row['status'] == 'on' for row in rows)
    # two datasets of one source are not compared
    data = merge([_set('HIRES-j', time, 1.0, rng),
                  _set('HIRES-k', time, 3.0, rng)])
    used, rows = datasets.choose(data, sources={'HIRES-j': 'CLS',
                                                'HIRES-k': 'CLS'}, auto=False)
    assert used.n == data.n


def test_a_dataset_that_constrains_nothing_is_left_out():
    """a few old, imprecise velocities beside many precise ones add
    nothing to the mean nor to the slope: left out; a long old series of
    fair precision holds the slope: kept"""
    rng = np.random.default_rng(5)
    harps = _set('HARPS03', np.sort(rng.uniform(56000, 58000, 200)), 1.0, rng)
    coarse = _set('CORAVEL', np.sort(rng.uniform(50000, 50400, 12)), 300.0,
                  rng)
    hires = _set('HIRES', np.sort(rng.uniform(50500, 56500, 80)), 2.5, rng)
    used, rows = datasets.choose(merge([harps, coarse, hires]))
    rows = _rows(rows)
    assert rows['CORAVEL']['status'] == 'weak'
    assert rows['CORAVEL']['mean'] < datasets.WEAK
    assert rows['CORAVEL']['slope'] < datasets.WEAK
    assert rows['HIRES']['status'] == 'on'
    assert rows['HIRES']['slope'] > 0.5
    assert used.instruments == ['HIRES', 'HARPS03']
    assert 'neither' in datasets.told(rows['CORAVEL'])
    # the dataset of the file given is never left out, nor one asked back
    for kwargs in (dict(protect=['CORAVEL']), dict(include=['coravel'])):
        used, rows = datasets.choose(merge([harps, coarse, hires]), **kwargs)
        assert 'CORAVEL' in used.instruments
    # with the rules off, everything is used
    used, rows = datasets.choose(merge([harps, coarse, hires]), auto=False)
    assert used.n == 292 and all(row['status'] == 'on' for row in rows)


def test_the_budget_is_that_of_a_line():
    """the errors of the slope with and without a dataset, against the
    covariance of the least squares (an offset per dataset, one slope)"""
    rng = np.random.default_rng(6)
    parts = {'A': (np.sort(rng.uniform(0, 1000, 40)), np.full(40, 1.0)),
             'B': (np.sort(rng.uniform(2000, 2600, 25)), np.full(25, 2.0)),
             'C': (np.sort(rng.uniform(500, 3500, 30)), np.full(30, 4.0))}

    def slope_error(names):
        time = np.concatenate([parts[name][0] for name in names])
        err = np.concatenate([parts[name][1] for name in names])
        design = np.zeros((len(time), len(names) + 1))
        start = 0
        for col, name in enumerate(names):
            design[start:start + len(parts[name][0]), col] = 1.0
            start += len(parts[name][0])
        design[:, -1] = time
        cov = np.linalg.inv(design.T @ (design / err[:, None] ** 2))
        return np.sqrt(cov[-1, -1])
    worth = datasets.budget(parts)
    full = slope_error(['A', 'B', 'C'])
    for name in parts:
        others = [key for key in parts if key != name]
        assert worth[name]['slope'] == pytest.approx(
            slope_error(others) / full - 1.0, rel=1e-6)
    # one point has an offset and nothing else
    alone = datasets.budget({'A': parts['A'],
                             'D': (np.array([5000.0]), np.array([0.1]))})
    assert alone['D'] == dict(mean=0.0, slope=0.0)
    assert alone['A']['mean'] == np.inf


def test_the_scatter_about_a_line():
    rng = np.random.default_rng(7)
    time = np.sort(rng.uniform(0, 1000, 300))
    rv = 0.05 * time + rng.normal(0, 2.0, 300)
    rv[::40] += 80.0
    assert datasets.line_scatter(time, rv) == pytest.approx(2.0, rel=0.2)
    assert np.isnan(datasets.line_scatter(time[:1], rv[:1]))


def test_two_releases_as_precise_as_each_other():
    """the same noise in both: whichever is preferred (chance can make one
    look better), each spectrum is used once and none is lost"""
    time = np.sort(np.random.default_rng(8).uniform(56000, 58000, 160))
    for seed in range(5):
        rng = np.random.default_rng(seed)
        data = merge([_set('HARPS (Paper+ 2015)', time[:60], 1.5, rng),
                      _set('HARPS03', time, 1.5, rng)])
        used, rows = datasets.choose(data, auto=False)
        assert used.n == 160
        assert len(np.unique(np.round(used.time, 4))) == 160
        assert sorted(row['status'] for row in rows) in (
            ['on', 'part'], ['on', 'release'])
    # told apart by nothing (the same velocities): the one with the more
    #   spectra is the one used
    rng = np.random.default_rng(9)
    full = _set('HARPS03', time, 1.5, rng)
    copy = RVData(time=full.time[:60], rv=full.rv[:60], err=full.err[:60],
                  inst=np.array(['HARPS (Paper+ 2015)'] * 60), name='star')
    used, rows = datasets.choose(merge([copy, full]), auto=False)
    assert _rows(rows)['HARPS03']['status'] == 'on'
    assert _rows(rows)['HARPS (Paper+ 2015)']['status'] == 'release'


def test_a_release_that_names_no_spectrograph():
    """a paper's table that names no instrument, its times a few minutes
    from those of HARPS on DACE: the same spectra (not two spectrographs
    that looked at the star in the same months)"""
    rng = np.random.default_rng(10)
    time = np.sort(rng.uniform(56000, 58000, 140))
    harps = _set('HARPS03', time, 2.5, rng)
    terra = _set('Paper+ 2012', time[:100] + 3.0 / 1440, 0.8, rng)
    other = _set('Other+ 2015', np.sort(rng.uniform(56000, 58000, 100)),
                 0.8, rng)
    data = merge([harps, terra, other])
    assert datasets.tolerance(data, 'HARPS03', 'Paper+ 2012') \
        == datasets.SAME_FAMILY
    assert datasets.tolerance(data, 'HARPS03', 'Other+ 2015') \
        == datasets.SAME_ANY
    used, rows = datasets.choose(data, auto=False)
    rows = _rows(rows)
    assert rows['HARPS03']['status'] == 'part'
    assert rows['HARPS03']['better'] == 'Paper+ 2012'
    assert rows['HARPS03']['used'] == 40
    assert rows['Other+ 2015']['status'] == 'on'
    assert used.n == 240


def test_a_weak_dataset_that_is_mostly_a_release():
    """a dataset whose spectra are another's but for a few: those few
    constrain nothing, and what is said of it says both"""
    rng = np.random.default_rng(11)
    time = np.sort(rng.uniform(56000, 58000, 150))
    data = merge([_set('HARPS03', time, 3.0, rng),
                  _set('HARPS (Trifonov+ 2020)', time[:146], 1.0, rng)])
    used, rows = datasets.choose(data)
    row = _rows(rows)['HARPS03']
    assert row['status'] == 'weak' and row['left'] == 4
    assert row['same'] == 146
    assert '146 of its 150 points' in datasets.told(row)
    assert used.instruments == ['HARPS (Trifonov+ 2020)']

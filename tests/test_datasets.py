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
    assert not rows['HARPS03']['tie']
    assert 'the more precise' in datasets.told(rows['HARPS03'])
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
    # said as it is: not told apart, the other preferred
    assert _rows(rows)['HARPS (Paper+ 2015)']['tie']
    assert 'as precise' in datasets.told(_rows(rows)['HARPS (Paper+ 2015)'])


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


def test_a_table_of_two_spectrographs_that_names_none():
    """a paper's table that names no instrument and holds the velocities
    of two spectrographs: those that are spectra of HARPS (a minute
    apart: UTC and TDB) are HARPS's, whatever the others are"""
    rng = np.random.default_rng(12)
    time = np.sort(rng.uniform(56000, 58000, 200))
    harps = _set('HARPS03', time, 1.0, rng)
    # 30 of HARPS's spectra and 150 of another spectrograph, the same
    #   nights, hours later
    mixed = _set('Paper+ 2019', np.concatenate([
        time[:30] + 1.08 / 1440, time[40:190] + 0.2]), 2.5, rng)
    data = merge([harps, mixed])
    assert datasets.tolerance(data, 'HARPS03', 'Paper+ 2019') \
        == datasets.SAME_FAMILY
    used, rows = datasets.choose(data, protect=['HARPS03'], auto=False)
    row = _rows(rows)['Paper+ 2019']
    assert row['status'] == 'part' and row['same'] == 30
    assert row['used'] == 150 and used.n == 350
    # with the other spectrograph named too: the table holds the spectra
    #   of two under one name (one offset), and comes last, however
    #   precise its velocities
    hires = _set('HIRES (Survey+ 2021)', time[40:190] + 0.2, 4.0, rng)
    data = merge([harps, mixed, hires])
    order, _ = datasets.ranking(data, datasets.links(data))
    assert order[-1] == 'Paper+ 2019'
    used, rows = datasets.choose(data, auto=False)
    assert _rows(rows)['Paper+ 2019']['status'] == 'release'
    assert used.n == 350


def test_as_precise_the_named_and_the_latest_release_first():
    """the same velocities in three releases: the one that names its
    spectrograph before the table that does not, the latest of two"""
    rng = np.random.default_rng(13)
    time = np.sort(rng.uniform(56000, 58000, 80))
    full = _set('HIRES (Old+ 2017)', time, 2.0, rng)

    def copy(name, num):
        return RVData(time=full.time[:num], rv=full.rv[:num],
                      err=full.err[:num], inst=np.array([name] * num),
                      name='star')
    data = merge([full, copy('Table+ 2022', 80),
                  copy('HIRES (New+ 2019)', 80)])
    order, _ = datasets.ranking(data, datasets.links(data))
    assert order == ['HIRES (New+ 2019)', 'HIRES (Old+ 2017)', 'Table+ 2022']
    assert datasets.year('HIRES (Teklu+ 2025)') == 2025
    assert datasets.year('HARPS03') == 9999
    # a survey that does not name its year: that of where it came from
    assert datasets.year('HIRES (CLS)', 'Rosenthal et al. 2021 (California '
                         'Legacy Survey)') == 2021
    assert datasets.year('CARMENES', 'CARMENES DR1') == 2023
    assert datasets.year('HARPS03', 'DACE') == 9999
    data = merge([copy('HIRES (CLS)', 80), copy('HIRES (New+ 2025)', 80)])
    order, _ = datasets.ranking(data, datasets.links(data), sources={
        'HIRES (CLS)': 'Rosenthal et al. 2021 (California Legacy Survey)'})
    assert order[0] == 'HIRES (New+ 2025)'


def test_the_ratings_of_several_releases():
    """three releases compared two by two, and a fourth compared with the
    worst only: the order of their noises, the fourth not first for having
    met a poor one"""
    pairs = {('A', 'B'): dict(ratio=0.5, n=100), ('B', 'C'): dict(ratio=0.5,
                                                                  n=100),
             ('A', 'C'): dict(ratio=0.25, n=100),
             ('D', 'C'): dict(ratio=0.9, n=100)}
    rate = datasets.ratings(['A', 'B', 'C', 'D', 'E'], pairs)
    assert rate['A'] < rate['B'] < rate['D'] < rate['C']
    assert rate['B'] - rate['A'] == pytest.approx(np.log(2), abs=1e-5)
    assert rate['E'] == 0.0
    # compared and found the same: the same rating, to the last digit
    same = datasets.ratings(['A', 'B', 'C'], {
        ('A', 'B'): dict(ratio=1.0, n=50), ('B', 'C'): dict(ratio=1.0, n=80)})
    assert same['A'] == same['B'] == same['C'] == 0.0


def test_what_the_star_does_is_taken_out_first():
    """two releases of the spectra of a star with a strong planet: told
    apart once its signal is out of both, not before; and the curve taken
    out holds the planet, not the noise"""
    rng = np.random.default_rng(14)
    time = np.sort(rng.uniform(56000, 58000, 180))
    star = 40.0 * np.sin(2 * np.pi * time / 5.37) \
        + 15.0 * np.sin(2 * np.pi * time / 12.9 + 1.0)

    def release(name, noise):
        return RVData(time=time, rv=star + rng.normal(0, noise, len(time)),
                      err=np.full(len(time), 1.2), inst=np.array(
                          [name] * len(time)), name='star')
    data = merge([release('HARPS03', 1.5), release('HARPS03 (RVBank)', 1.0)])
    curve = datasets.star_signal(data)
    assert sorted(np.round(curve.periods, 1)) == [5.4, 12.9]
    assert np.std(star - curve(time) - np.mean(star - curve(time))) < 0.5
    plain = datasets.compare(data, 'HARPS03', 'HARPS03 (RVBank)')
    clean = datasets.compare(data, 'HARPS03', 'HARPS03 (RVBank)',
                             signal=curve)
    # with the planets in, chance decides (the noises do not stand out of
    #   30 m/s of signal); without them, the noisier is the noisier
    assert plain['ratio'] == 1.0
    assert clean['ratio'] > 1.2
    assert clean['one'] == pytest.approx(1.5, rel=0.25)
    assert clean['other'] == pytest.approx(1.0, rel=0.3)
    used, rows = datasets.choose(data, auto=False)
    assert used.instruments == ['HARPS03 (RVBank)']
    # a quiet star: nothing but a line is taken out (no sinusoid fitted
    #   to the noise)
    quiet = merge([RVData(time=time, rv=rng.normal(0, 2.0, len(time)),
                          err=np.full(len(time), 2.0),
                          inst=np.array(['A'] * len(time)), name='star')])
    assert datasets.star_signal(quiet).periods == []


def test_a_release_binned_by_night_is_compared_bin_for_bin():
    """a paper that gives one velocity for the three exposures of a night
    is not the more precise for it: against the mean of the other's three
    exposures, the same noise; and the comparison reads the same from
    either side"""
    rng = np.random.default_rng(15)
    nights = np.sort(rng.choice(np.arange(56000, 58000), 90, replace=False))
    time = np.sort(np.concatenate([nights + 0.5 + it * 4.0 / 1440
                                   for it in range(3)]))
    star = _star(time)
    exposures = RVData(time=time, rv=star + rng.normal(0, 2.0, len(time)),
                       err=np.full(len(time), 2.0),
                       inst=np.array(['HIRES (Survey+ 2021)'] * len(time)),
                       name='star')
    binned = RVData(time=nights + 0.5 + 4.0 / 1440, rv=_star(
        nights + 0.5 + 4.0 / 1440) + rng.normal(0, 2.0 / np.sqrt(3), 90),
        err=np.full(90, 2.0 / np.sqrt(3)),
        inst=np.array(['HIRES (Paper+ 2010)'] * 90), name='star')
    data = merge([exposures, binned])
    one = datasets.compare(data, 'HIRES (Survey+ 2021)',
                           'HIRES (Paper+ 2010)')
    other = datasets.compare(data, 'HIRES (Paper+ 2010)',
                             'HIRES (Survey+ 2021)')
    assert one['by'] == 'noise' and one['ratio'] == 1.0
    assert other['ratio'] == 1.0
    assert one['n'] == 270 and other['n'] == 90
    assert one['one'] == pytest.approx(other['other'])
    # whichever is used (as precise: chance can make one look better),
    #   each night is used once
    used, rows = datasets.choose(data, auto=False)
    assert len(used.instruments) == 1
    assert len(np.unique(np.floor(used.time))) == 90
    # with nothing to tell them by, the latest release
    order, _ = datasets.ranking(data, datasets.links(data))
    assert order == ['HIRES (Survey+ 2021)', 'HIRES (Paper+ 2010)']

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Planet or no planet on the period and its aliases; which alias; the plan

Created on 2026-09-30

@author: artigau
"""
import numpy as np

from koloa import aliases as kal
from koloa.fip import PeriodGrid, oafip
from koloa.periodogram import ALIAS_FREQUENCIES
from koloa.simulate import simulate


def test_family_fip_is_the_union_of_the_aliases():
    rng = np.random.default_rng(1)
    grid = PeriodGrid(np.sort(rng.uniform(0, 364, 150)), 1.1, None, 10)
    pon = rng.random(grid.size) ** 8
    poff, covered = 0.3, [int(grid.size * 0.4)]
    total = poff + pon.sum()
    fam = grid.bin_family(pon, poff, total, covered)
    single = grid.bin_fip(pon, poff, total, covered)
    # a union is never less probable than one of its intervals
    assert np.all(fam <= single + 1e-12)
    step = 0.5 * np.median(np.diff(grid.freq))
    for kk in rng.integers(0, grid.size, 100):
        freq = grid.freq[kk]
        mask = np.zeros(grid.size, bool)
        for mem in [freq] + [val for fs in ALIAS_FREQUENCIES.values()
                             for val in (abs(freq - fs), freq + fs)]:
            if grid.freq[0] - step <= mem <= grid.freq[-1] + step:
                cen = int(np.argmin(np.abs(grid.freq - mem)))
                mask[max(cen - grid.halfbin, 0):cen + grid.halfbin + 1] = True
        ref = 0.0 if mask[covered[0]] else (poff + pon[~mask].sum()) / total
        assert abs(ref - fam[kk]) < 1e-12


def test_planet_or_no_planet_is_the_period_and_its_aliases(tmp_path):
    from koloa.diagnostics import duck_test
    sim = simulate(planets=[dict(P=3.3, K=6.0, e=0.0, tp=0.0)], seed=2,
                   err=2.0)
    data = sim['data']
    fip = oafip(data, kmax=1, nsweep=300, nburn=200, nchains=1,
                progress=False)
    width = 1 / data.baseline
    fam = fip.family_containing(3.3, width)
    assert fam is not None and fam <= fip.fip_containing(3.3, width) + 1e-12
    odds = fip.alias_odds(3.3, width)
    assert abs(sum(mem['share'] for mem in odds) - 1) < 1e-9
    # no yearly alias within 1.5 widths of the period
    assert all(abs(1 / mem['period'] - 1 / 3.3) > 1.5 * width
               for mem in odds if mem['name'] != 'P')
    report = duck_test(data, 3.3, fipres=fip, gp=False, quiet=True,
                       aliases=True, archive=False)
    assert report.details['family_fip'] == fam
    assert [ch for ch in report.checks if ch['name'] == 'period'][0][
        'status'] == 'info'
    sols = report.details['alias_solutions']
    assert sols and abs(sols[0]['period'] - 3.3) < 0.05
    assert 'or any of its aliases' in [ch for ch in report.checks if
                                       ch['name'] == 'significance'][0][
        'summary']
    # the plan, from La Silla, for a star it can see
    if len(sols) > 1:
        pln = kal.plan(sols, ra=100.0, dec=-30.0, site='La Silla',
                       start=61300.0, ndays=20)
        assert pln['nights'] and np.isfinite(
            pln['nights'][0]['versus'][0]['single']['chi2'])
        assert kal.plan_text(pln)


def test_family_and_sites():
    assert kal.same_family(6.92, 1.1652, 1 / 364)
    assert not kal.same_family(6.92, 3.0, 1 / 364)
    assert kal.site_of(['SPIRou']) == 'CFHT'
    assert kal.site_of(['HARPN', 'HARPS15']) == 'La Palma'
    assert kal.site_of(['HARPS03']) == 'La Silla'
    assert kal.site_of(['inst']) is None
    # the Sun at noon in Greenwich on an equinox is high, at midnight low
    jd = np.array([2461120.0, 2461120.5])  # 2026-03-20 12:00 and 24:00 UT
    ra, dec = kal._sun(jd)
    alt = kal.altitude(ra, dec, jd, 51.48, 0.0)
    assert alt[0] > 30 and alt[1] < -30


def test_detection_map_by_the_fip():
    """the adaptive map: rounds of injections per band, a strong planet
    found, the levels ordered"""
    from koloa.fipmap import fip_map, _logistic_fit, _levels
    sim = simulate(seed=5, err=1.0)
    dmap = fip_map(sim['data'], 1.5, 40.0, nband=2, nround=2, per_round=2,
                   kmax=1, nsweep=60, nburn=40, workers=2, nboot=20)
    assert len(dmap['injections']) == 8
    assert dmap['K90'][0] >= dmap['K50'][0] > 0
    # the logistic fit: found above ln K = 1, missed below
    lnk = np.linspace(-1, 3, 40)
    k50, k90 = _levels(*_logistic_fit(lnk, lnk > 1))
    assert abs(np.log(k50) - 1) < 0.2 and k90 > k50

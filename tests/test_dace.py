#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests of koloa.dace: the decoding of DACE's answers and the series built
from them (no network).

Created on 2026-09-29

@author: artigau
"""
import csv

import numpy as np
import pytest

from koloa import dace


def test_columns_undo_the_run_length_encoding():
    payload = {'parameters': [
        {'variableName': 'instrument_name', 'stringValues': ['HARPS03', 'ESPRESSO19'],
         'occurrences': [2, 1]},
        {'variableName': 'rv', 'doubleValues': [1.0, 2.0, 3.0],
         'minErrorValues': [0.5, 0.6, 0.7]}]}
    cols = dace._columns(payload)
    assert cols['instrument_name'] == ['HARPS03', 'HARPS03', 'ESPRESSO19']
    assert cols['rv'] == [1.0, 2.0, 3.0]
    assert cols['rv_err'] == [0.5, 0.6, 0.7]


def test_rvdata_keeps_good_velocities_one_per_spectrum_per_era(tmp_path):
    path = tmp_path / 'star.csv'
    rows = []
    rng = np.random.default_rng(1)
    for inst, mode, num, zero in (('HARPS03', 'HARPS', 8, 100.0),
                                  ('ESPRESSO19', 'SINGLEHR11', 6, 200.0),
                                  ('ESPRESSO19', 'SINGLEUHR', 2, 200.0)):
        for it in range(num):
            rows.append(dict(rjd=55000 + 10 * len(rows), rv=zero + rng.normal(),
                             rv_err=1.0, drs_qc='True', instrument_name=inst,
                             ins_mode=mode, file_rootname=f'{inst}{mode}{it}',
                             ccf_fwhm=7000.0, ccf_fwhm_err=5.0))
    # a spectrum twice, a failed one and one without an error
    rows.append(dict(rows[0]))
    rows.append(dict(rows[1], drs_qc='False', file_rootname='bad1'))
    rows.append(dict(rows[2], rv_err='', file_rootname='bad2'))
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    data = dace.rvdata(str(path), name='star', min_points=3)
    # the main ESPRESSO mode keeps the name of the era; the other, too
    #   small, is dropped
    assert sorted(data.instruments) == ['ESPRESSO19', 'HARPS03']
    assert np.sum(data.inst == 'HARPS03') == 8
    assert np.sum(data.inst == 'ESPRESSO19') == 6
    # each era has its own zero point
    assert abs(data.zero_point['HARPS03'] - 100.0) < 2.0
    assert 'fwhm' in data.indicators


def test_split_sequences_keeps_every_point_in_short_visits():
    from koloa.data import split_sequences
    seq = np.array([0] * 3 + [1] * 23 + [2] * 2)
    out = split_sequences(seq, max_size=10)
    sizes = np.bincount(out)
    assert sizes.max() <= 10 and sizes.sum() == len(seq)
    # the short visits untouched, the long one in three chunks, in order
    assert list(sizes) == [3, 8, 8, 7, 2]
    assert np.all(np.diff(out) >= 0)


def test_a_silent_dace_asked_again(tmp_path, monkeypatch):
    """a question DACE leaves unanswered is asked again (it stalls now and
    then, and answers at once the second time); one it answers that it
    does not know the target is not; and three silences are said"""
    import io
    import json
    import urllib.error

    from koloa import dace
    payload = json.dumps({'parameters': [
        {'variableName': 'rjd', 'doubleValues': [1.0, 2.0]},
        {'variableName': 'rv', 'doubleValues': [3.0, 4.0],
         'minErrorValues': [0.5, 0.6]}]}).encode()
    calls = []

    def answer(script):
        def urlopen(request, timeout):
            calls.append(timeout)
            what = script[len(calls) - 1]
            if isinstance(what, Exception):
                raise what
            return io.BytesIO(what)
        return urlopen
    monkeypatch.setattr(dace.urllib.request, 'urlopen',
                        answer([TimeoutError('timed out'), payload]))
    path = dace.fetch('GJ1', str(tmp_path / 'a.csv'), api_key=False)
    assert len(calls) == 2 and calls[0] == 45.0
    assert open(path).readline().strip() == 'rjd,rv,rv_err'
    calls.clear()
    missing = urllib.error.HTTPError('u', 404, 'Not Found', {}, None)
    monkeypatch.setattr(dace.urllib.request, 'urlopen', answer([missing]))
    with pytest.raises(RuntimeError, match='DACE does not know GJ2'):
        dace.fetch('GJ2', str(tmp_path / 'b.csv'), api_key=False)
    assert len(calls) == 1
    calls.clear()
    monkeypatch.setattr(dace.urllib.request, 'urlopen',
                        answer([TimeoutError('timed out')] * 3))
    with pytest.raises(RuntimeError, match='DACE could not be reached'):
        dace.fetch('GJ3', str(tmp_path / 'c.csv'), api_key=False)
    assert len(calls) == 3

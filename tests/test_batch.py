#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A batch of detailed reports: its script, and its summary

Created on 2026-10-02

@author: artigau
"""
import json
import os
import subprocess

from koloa import batch, gui


def test_the_script_of_a_batch(tmp_path):
    """one line per target, each its folder, the koloa of the machine (an
    environment variable before it exported), a valid bash script"""
    out = batch.script(
        [dict(target='GJ 436', files=[dict(path='/data/a b.rdb', label=''),
                                      dict(path='/data/c.rdb',
                                           label='NIRPS_LBL2')]),
         dict(target='GJ 436', files=[]), dict(target='', files=[]),
         dict(target='', files=[dict(path='/data/TOI-1.rdb', label='')])],
        dict(toi_on=True, curvature=True),
        koloa='PYTHONPATH=/k /py/bin/python -m koloa.cli', batch='b1',
        workdir='/work', jobs=3)
    text = out['script']
    assert out['targets'] == ['GJ_436', 'GJ_436_2', 'TOI-1']
    assert 'export PYTHONPATH=/k' in text and 'KOLOA=(/py/bin/python -m ' \
        'koloa.cli)' in text and 'JOBS=3' in text
    assert ("one 'GJ 436' GJ_436 '/data/a b.rdb' /data/c.rdb --detailed "
            "--instruments auto NIRPS_LBL2 --target 'GJ 436' --toi "
            "--curvature &") in text
    assert '--batch-summary "$BATCH"' in text
    path = tmp_path / out['name']
    path.write_text(text)
    assert subprocess.run(['bash', '-n', str(path)]).returncode == 0


def test_the_summary_of_a_batch(tmp_path):
    """a target done (its summary read) and one failed (its error from its
    log), as text, csv and pdf"""
    os.makedirs(tmp_path / 'GJ_1')
    os.makedirs(tmp_path / 'logs')
    (tmp_path / 'GJ_1' / 'GJ_1_summary.json').write_text(json.dumps(dict(
        star='GJ 1', n=40, nvisits=40, instruments=dict(NIRPS=40),
        baseline=300.0, known=dict(planets=[]), duck={'5.3000': 'PLANET'},
        orbits=[dict(P=[5.3, 0.01, 0.01], K=[3.0, 0.2, 0.3],
                     family_fip=1e-5, origin='FIP')],
        acceleration=dict(accel=[1.0, 0.5, 0.5]))))
    (tmp_path / 'logs' / 'GJ_2.log').write_text(
        'koloa, detailed: GJ 2\nTraceback\nFileNotFoundError: no such file\n')
    out = batch.summary(str(tmp_path))
    text = open(out['txt']).read()
    assert 'GJ 1: 40 points (NIRPS 40)' in text and '[PLANET]' in text
    assert 'GJ_2: FAILED: FileNotFoundError: no such file' in text
    assert os.path.exists(out['pdf']) and os.path.exists(out['csv'])


def test_a_folder_of_this_machine(tmp_path):
    """the file browser of the page, here: folders first, no hidden files"""
    (tmp_path / 'sub').mkdir()
    (tmp_path / 'a.rdb').write_text('x')
    (tmp_path / '.hidden').write_text('x')
    res = gui.listing('', str(tmp_path))
    assert [ent['name'] for ent in res['entries']] == ['sub', 'a.rdb']
    assert res['entries'][0]['dir'] and res['path'] == str(tmp_path)

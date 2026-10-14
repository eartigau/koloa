#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A package of LBL for a server (koloa.lblbatch), without the network: the
plan of its stars, the list of what is fetched, what LBL gave as the files
of koloa's batch, the parts of a job array, the tab of the GUI

Created on 2026-10-14

@author: artigau
"""
import csv
import json
import os

from koloa import lblbatch, lblrun


def _runs(name, inst, root='archives', ident=None):
    """the runs of an instrument, as koloa.lblrun.groups makes them"""
    obj = lblrun.object_name(name)
    return [dict(key='', obj=obj, lbl='HARPSN', source='ESO',
                 input='*S2D_A.fits', tellu=True, params={},
                 spectra=[dict(id=f'{obj}{it}.fits.gz', name=f'{obj}{it}.fits')
                          for it in range(3)],
                 eras=[('', 'HARPS-N (LBL)') + lblrun.ALWAYS],
                 get=lblrun._get_tng, via='tng')]


def test_the_plan_and_its_list(tmp_path, monkeypatch):
    from koloa import shared, stars
    monkeypatch.setattr(lblrun, '_ident', lambda name: dict(
        main=name.upper(), aliases=[name + ' b'], ra=1.0, dec=2.0))
    monkeypatch.setattr(stars, 'stellar_teff', lambda ident: dict(teff=3300))
    monkeypatch.setattr(lblrun, 'groups', _runs)
    shared.save_settings(False, ['GJ 2'])
    plan = lblbatch.plan(['GJ 1', 'GJ 2'], ['HARPSN'])
    assert [star['name'] for star in plan] == ['GJ 1', 'GJ 2']
    assert plan[0]['teff'] == 3300.0 and not plan[0]['never']
    assert plan[1]['never']
    group = plan[0]['groups'][0]
    assert 'get' not in group and group['via'] == 'tng'
    json.dumps(plan)
    live = lblrun.live_group(json.loads(json.dumps(group)))
    assert live['get'] is lblrun._get_tng and live['eras'][0][1] == \
        'HARPS-N (LBL)'
    lblbatch.fetch_list(plan, str(tmp_path / 'list.csv'))
    rows = list(csv.DictReader(open(tmp_path / 'list.csv')))
    assert rows[0]['star'] == 'GJ 1' and rows[0]['spectra'] == '3'
    assert int(rows[0]['MB']) == round(3 * lblrun.INSTRUMENTS['HARPSN'][
        'size'])
    # the parts of a job array
    many = [dict(name=f'S{it}') for it in range(7)]
    assert [one['name'] for one in lblbatch._mine(many, (2, 3), [])] == [
        'S1', 'S4']
    assert [one['name'] for one in lblbatch._mine(many, None, ['S5'])] == [
        'S5']


def test_what_lbl_gave_as_the_files_of_the_batch(tmp_path):
    root = tmp_path / 'pack'
    (root / lblbatch.DONE).mkdir(parents=True)
    (root / 'targets.json').write_text(json.dumps(dict(name='x', targets=[
        dict(name='GJ 1', files=[]), dict(name='GJ 2', files=[])])))
    (root / lblbatch.DONE / 'GJ_1.json').write_text(json.dumps(dict(
        name='GJ 1', errors={}, files=[
            dict(path='archives/GJ_1/lbl/GJ_1_SOPHIE_lbl.rdb',
                 label='SOPHIE (LBL)', n=3),
            dict(path='archives/GJ_1/lbl/GJ_1_SOPHIEp_lbl.rdb',
                 label='SOPHIE+ (LBL)', n=4)])))
    assert lblbatch.merge(str(root)) == 1
    batch = json.load(open(root / 'targets.json'))
    assert batch['targets'][0]['labels'] == ['SOPHIE (LBL)', 'SOPHIE+ (LBL)']
    assert batch['targets'][1]['files'] == []


def test_the_tab_of_the_batch(tmp_path):
    from koloa import gui
    args = gui.command('lblpack', dict(
        stars='GJ 436\nGl 408, GJ 1002', lbl_HARPSN=True, lbl_SOPHIE=True,
        out=str(tmp_path), server_root='/scratch/me/x', nproc='8',
        push=False, analyse=True, root='archives'))
    assert args[:4] == ['GJ 436', 'Gl 408', 'GJ 1002', '--lbl-pack']
    assert args[4:6] == ['HARPSN', 'SOPHIE']
    assert args[args.index('--pack-name') + 1] == 'lbl_GJ436_Gl408_GJ1002'
    assert '--no-push' in args and '--no-analyse' not in args
    assert '--fetch' not in args
    assert args[args.index('--server-root') + 1] == '/scratch/me/x'
    # what a package holds
    folder = tmp_path / 'lbl_x'
    folder.mkdir()
    (folder / lblbatch.PLAN).write_text(json.dumps(dict(fetched=False, stars=[
        dict(name='GJ 1', never=True, errors=dict(SOPHIE='down'), groups=[
            dict(inst='HARPSN', key='', spectra=[dict(id='a', name='a')] * 4)
        ])])))
    got = gui.lbl_package(str(folder))
    assert got['nspectra'] == 4 and got['rows'][0]['never']
    assert got['rows'][1]['error'] == 'down' and got['tar'] is None

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa's GUI: the command lines it shows (and runs), and its server, without
the network

Created on 2026-10-01

@author: artigau
"""
import json
import threading
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from koloa import gui


def test_the_command_lines():
    assert gui.command('gather', dict(target='GJ 436', root='arch',
                                      dace=True, carmenes=False,
                                      tess=True)) == [
        'GJ 436', '--gather', 'arch', '--no-carmenes']
    args = gui.command('detailed', dict(
        target='GJ 436', file='star.rdb', outdir='out', kmax='4',
        nsweep='2000', pmin='1.1', periods='113.46, 27', detection_map='fip',
        fip_gp='none', tess=False, exclude='NIRPS, HARPS03'))
    # the defaults of the command line are not written
    assert args == ['star.rdb', '--detailed', '--target', 'GJ 436',
                    '--outdir', 'out', '--kmax', '4', '--periods', '113.46',
                    '27', '--no-fip-gp', '--exclude', 'NIRPS', 'HARPS03',
                    '--detection-map', '--no-tess']
    assert gui.line(args).startswith("koloa star.rdb --detailed --target "
                                     "'GJ 436'")
    # a rotation that can be trusted: an SHO at it, unless asked otherwise
    args = gui.command('detailed', dict(target='GL 406', file='a.rdb', rotation='2.704',
                                        fip_gp='sho'))
    assert args[-2:] == ['--rotation', '2.704']
    assert gui.command('detailed', dict(target='GL 406', file='a.rdb', rotation='2.704',
                                        fip_gp='banded'))[-2:] == [
        '--fip-gp', 'banded']
    with pytest.raises(ValueError):
        gui.command('detailed', dict(target='GL 406', file='a.rdb', fip_gp='sho'))
    assert gui.command('detailed', dict(target='x', file='a.rdb', fip_gp='none'))[-1] == (
        '--no-fip-gp')
    assert gui.command('gather', dict(target='x', refresh=True))[-1] == (
        '--refresh')
    # the archives only when asked: a file alone fetches nothing, a name
    #   alone needs an archive ticked
    assert gui.command('detailed', dict(file='a.rdb')) == [
        'a.rdb', '--detailed', '--outdir', 'koloa_output']
    args = gui.command('detailed', dict(target='GJ 436', dace=True,
                                        carmenes=True))
    assert args[0] == '--detailed' and '--dace' in args and \
        '--carmenes' in args and '--vizier' not in args
    with pytest.raises(ValueError):
        gui.command('detailed', dict(target='GJ 436'))
    with pytest.raises(ValueError):
        gui.command('detailed', dict())
    with pytest.raises(ValueError):
        gui.command('gather', dict(file='a.rdb'))


def test_the_server_answers():
    server = ThreadingHTTPServer(('127.0.0.1', 0), gui.Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f'http://127.0.0.1:{server.server_address[1]}'
    try:
        with urllib.request.urlopen(base + '/') as resp:
            assert b'koloa' in resp.read()
        with urllib.request.urlopen(base + '/api/info') as resp:
            assert json.loads(resp.read())['defaults']['kmax'] == 3
        body = json.dumps(dict(action='gather', options=dict(
            target='GJ 436'))).encode()
        request = urllib.request.Request(base + '/api/command', data=body)
        with urllib.request.urlopen(request) as resp:
            assert json.loads(resp.read())['line'] == (
                "koloa 'GJ 436' --gather archives")
    finally:
        server.shutdown()
        server.server_close()


def test_the_dialog_that_picks_a_file(monkeypatch, tmp_path):
    """the dialog of the machine (here a stand-in): the path chosen, a
    folder under the one koloa runs from made relative, a cancel"""
    import subprocess
    import types
    answers = iter([
        types.SimpleNamespace(returncode=0, stdout='/data/lbl_GJ436.rdb\n',
                              stderr=''),
        types.SimpleNamespace(returncode=0, stdout=f'{tmp_path}/out/\n',
                              stderr=''),
        types.SimpleNamespace(returncode=1, stdout='',
                              stderr='execution error: User canceled. '
                                     '(-128)')])
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: next(answers))
    monkeypatch.setattr(gui.sys, 'platform', 'darwin')
    monkeypatch.chdir(tmp_path)
    assert gui.pick('file') == dict(path='/data/lbl_GJ436.rdb')
    assert gui.pick('folder', str(tmp_path)) == dict(path='out')
    assert gui.pick('file') == dict(cancelled=True)


def test_the_file_and_the_archives_are_told_apart(tmp_path):
    """NIRPS in the file and on DACE: two instruments, NIRPS (the file) and
    NIRPS_DACE (what DACE has that the file does not), each with its
    source; the same spectra are not counted twice"""
    import numpy as np
    from koloa.data import RVData
    from koloa.gather import write_rv
    rng = np.random.default_rng(1)
    tfile = np.sort(rng.uniform(60000, 60300, 30))
    # an LBL file of NIRPS (its EXTSN060 column says so)
    (tmp_path / 'lbl.csv').write_text('rjd,vrad,svrad,EXTSN060\n' + ''.join(
        f'{tt!r},{rng.normal(0, 3)!r},1.0,100.0\n' for tt in tfile))
    tnew = np.sort(rng.uniform(60400, 60500, 10))
    tharps = np.sort(rng.uniform(58000, 59000, 15))
    tcarm = np.sort(rng.uniform(57400, 57700, 12))
    times = np.concatenate([tfile[:20], tnew, tharps, tcarm])
    insts = np.array(['NIRPS'] * 30 + ['HARPS15'] * 15 + ['CARMENES'] * 12)
    write_rv(RVData(times, rng.normal(0, 3, len(times)),
                    np.full(len(times), 1.0), inst=insts),
             str(tmp_path / 'arch' / 'GJ_436' / 'rv' / 'all_rv.csv'))
    # the archives not asked for: the file alone, and a word of what was
    #   gathered and not ticked
    res = gui.velocities(str(tmp_path / 'lbl.csv'), 'GJ 436',
                         str(tmp_path / 'arch'))
    assert [inst['name'] for inst in res['instruments']] == ['NIRPS']
    assert 'DACE: 45 points gathered, not ticked' in res['notes']
    # no file, nothing ticked: nothing, and why
    res = gui.velocities('', 'GJ 436', str(tmp_path / 'arch'))
    assert res['instruments'] == [] and any(
        'CARMENES DR1: 12 points' in note for note in res['notes'])
    res = gui.velocities(str(tmp_path / 'lbl.csv'), 'GJ 436',
                         str(tmp_path / 'arch'), dace=True, carmenes=True)
    got = {inst['name']: (inst['source'], inst['n'])
           for inst in res['instruments']}
    assert got == {'NIRPS': ('file: lbl.csv', 30), 'NIRPS_DACE': ('DACE', 10),
                   'HARPS15': ('DACE', 15),
                   'CARMENES': ('CARMENES DR1', 12)}


def test_a_star_asked_again_is_read_from_the_disk(tmp_path, monkeypatch):
    """the resolver asks the network once; then the disk answers (the copy
    kept, or the star's archives folder first), unless refreshed"""
    import json
    from koloa import archive
    calls = []

    def ask(name, refresh=False):
        calls.append(name)
        return dict(name=name, main='Ross 905', aliases=[], tic=None,
                    variability=[], carmenes=None)
    monkeypatch.setattr(gui, '_ask_star', ask)
    monkeypatch.setattr(archive, 'CACHE', str(tmp_path / 'cache'))
    monkeypatch.setattr(archive, 'host_name', lambda ident: None)
    first = gui.resolve_star('GJ 436', str(tmp_path / 'arch'))
    again = gui.resolve_star('GJ 436', str(tmp_path / 'arch'))
    assert calls == ['GJ 436'] and 'disk' not in first
    assert again['disk'] == 'the copy kept' and again['main'] == 'Ross 905'
    # the star's archives folder (koloa.gather) comes first
    folder = tmp_path / 'arch' / 'GJ_436'
    folder.mkdir(parents=True)
    (folder / 'target.json').write_text(json.dumps(dict(
        name='GJ 436', main='GJ 436 (gathered)', aliases=[], tic=None)))
    assert gui.resolve_star('GJ 436', str(tmp_path / 'arch'))['main'] == (
        'GJ 436 (gathered)')
    gui.resolve_star('GJ 436', str(tmp_path / 'arch'), refresh=True)
    assert calls == ['GJ 436', 'GJ 436']


def test_the_quick_look(tmp_path, monkeypatch):
    """the quick FIP of what the page shows (a planet in a file), its two
    passes when a signal is found, its curves; the PDF of the page"""
    import time
    from koloa.simulate import simulate
    from koloa.gather import write_rv
    sim = simulate(planets=[dict(P=5.3, K=8.0, e=0.0)], err=1.5, seed=3,
                   nvisits=50, per_visit=1, baseline=300)['data']
    write_rv(sim, str(tmp_path / 'star.csv'))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    opts = dict(files=[dict(path=str(tmp_path / 'star.csv'))])
    state = gui.quick_fip(opts)
    for _ in range(300):
        state = gui.quick_state(state['id'])
        if state['status'] != 'running':
            break
        time.sleep(0.5)
    assert state['status'] == 'done', state['error']
    res = state['result']
    assert len(res['period']) == len(res['family']) == len(res['alone'])
    assert abs(res['peaks'][0]['period'] / 5.3 - 1) < 0.01
    assert res['passes'] == 2 and res['window']['year'] == 365.25
    # the peaks numbered, one per family of aliases, the planet first, and
    #   folded there with its K
    first = res['peak_list'][0]
    assert first['id'] == 1 and first['named']
    assert abs(first['period'] / 5.3 - 1) < 0.01
    assert abs(res['folds'][0]['K'] - 8.0) < 4 * res['folds'][0]['K_err']
    assert state['each'] == []   # one instrument: the joint FIP is its own
    pdf = gui.quicklook_pdf(opts, qid=state['id'],
                            command_line='koloa star.csv --detailed')
    assert pdf[:4] == b'%PDF'


def test_each_instrument_has_its_own_quick_fip(tmp_path, monkeypatch):
    """several instruments: the joint FIP, then each instrument's own (one
    with too few nights is skipped), in the PDF too; a new target forgets
    the quick looks"""
    import time
    from koloa.simulate import simulate
    from koloa.gather import write_rv
    files = []
    for seed, label, nvisits in ((3, 'AAA', 40), (4, 'BBB', 30),
                                 (5, 'CCC', 6)):
        sim = simulate(planets=[dict(P=5.3, K=8.0, e=0.0, tp=60000.0)],
                       err=1.5, seed=seed, nvisits=nvisits, per_visit=1,
                       baseline=300)['data']
        if label == 'AAA':
            sim.rv[5] += 80.0   # an outlier
        write_rv(sim, str(tmp_path / f'{label}.csv'))
        files.append(dict(path=str(tmp_path / f'{label}.csv'), label=label))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    opts = dict(files=files)
    state = gui.quick_fip(opts)
    for _ in range(600):
        state = gui.quick_state(state['id'])
        if state['status'] != 'running':
            break
        time.sleep(0.5)
    assert state['status'] == 'done', state['error']
    each = {one['name']: one for one in state['each']}
    assert set(each) == {'AAA', 'BBB', 'CCC'}
    assert each['CCC']['skipped'] and each['CCC']['n'] == 6
    for name in ('AAA', 'BBB'):
        one = each[name]
        assert not one['skipped']
        assert len(one['period']) == len(one['family']) == len(one['alone'])
        assert one['peaks'] and one['peaks'][0]['period'] > 0
    assert abs(each['AAA']['peaks'][0]['period'] / 5.3 - 1) < 0.01, each['AAA']['peaks']
    assert abs(state['result']['peaks'][0]['period'] / 5.3 - 1) < 0.01
    # each night of the fold with its probability to be valid: the outlier
    #   below an even chance, the others above
    fold = state['result']['folds'][0]
    valid = {one['name']: one['valid'] for one in fold['instruments']}
    assert all(len(valid[one['name']]) == len(one['phase'])
               for one in fold['instruments'])
    assert sum(val < 0.5 for val in valid['AAA']) == 1
    assert sum(val < 0.5 for val in valid['BBB']) == 0
    # the outlier hardly counts in the fit of the fold
    assert fold['K_err'] < 0.6 and abs(fold['K'] - 8.0) < 4 * fold['K_err']
    pdf = gui.quicklook_pdf(opts, qid=state['id'])
    assert pdf[:4] == b'%PDF'
    gui.forget()
    assert state['id'] not in gui.QUICKS


def test_a_result_remembered_and_recalled(tmp_path, monkeypatch):
    """a quick look remembered (the page, its FIP, the velocities, a copy of
    the file and of the archives), listed, recalled as it was: the
    originals while they have not changed, else their copies; forgotten"""
    import os
    import numpy as np
    from koloa.data import RVData
    from koloa.gather import write_rv
    monkeypatch.setattr(gui, 'REMEMBERED', str(tmp_path / 'remembered'))
    rng = np.random.default_rng(2)
    tfile = np.sort(rng.uniform(60000, 60300, 30))
    lbl = tmp_path / 'lbl.csv'
    lbl.write_text('rjd,vrad,svrad,EXTSN060\n' + ''.join(
        f'{tt!r},{rng.normal(0, 3)!r},1.0,100.0\n' for tt in tfile))
    tarch = np.sort(rng.uniform(58000, 59000, 15))
    write_rv(RVData(tarch, rng.normal(0, 3, 15), np.full(15, 1.0),
                    inst=np.array(['HARPS15'] * 15)),
             str(tmp_path / 'arch' / 'GJ_436' / 'rv' / 'all_rv.csv'))
    qid = 'remember1'
    gui.QUICKS[qid] = dict(
        id=qid, status='done', step='done', step_detail='', progress=None,
        error=None, start=0.0, end=12.0, instruments=['NIRPS', 'HARPS15'],
        each=[], result=dict(n=45, instruments={'NIRPS': 30, 'HARPS15': 15},
                             known=[], peak_list=[dict(
                                 id=1, period=2.644, family=1e-5,
                                 alone=1e-4, named=True)]))
    page = dict(target='GJ 436', files=[dict(path=str(lbl), label='')],
                root=str(tmp_path / 'arch'), outdir='',
                detailed=dict(dace=True, carmenes=False, exclude=''),
                clip=False, view=None, periods=None)
    # a quick FIP that has not ended is not remembered
    with pytest.raises(ValueError):
        gui.remember(page, 'nothing')
    out = gui.remember(page, qid, note='a peak at 2.64 d')
    listed = gui.remembered()
    assert [entry['id'] for entry in listed] == [out['id']]
    assert listed[0]['note'] == 'a peak at 2.64 d'
    assert listed[0]['summary']['archives'] == ['DACE']
    assert listed[0]['summary']['peaks'][0]['period'] == 2.644
    # recalled: the originals, unchanged; the velocities and the FIP kept
    rec = gui.recall(out['id'])
    assert rec['page']['files'][0]['path'] == str(lbl)
    assert rec['page']['root'] == str(tmp_path / 'arch')
    assert rec['notes'] == []
    assert {inst['name']: inst['source'] for inst in rec['rv']['instruments']
            } == {'NIRPS': 'file: lbl.csv', 'HARPS15': 'DACE'}
    state = gui.quick_state(rec['quick']['id'])
    assert state['status'] == 'done' and state['result']['n'] == 45
    assert state['elapsed'] == 12.0
    # the file changed and the archives moved since: their copies
    lbl.write_text(lbl.read_text() + '60301.0,0.0,1.0,100.0\n')
    os.rename(tmp_path / 'arch', tmp_path / 'arch_moved')
    rec = gui.recall(out['id'])
    copy = rec['page']['files'][0]['path']
    assert copy.startswith(str(tmp_path / 'remembered')) and len(rec['notes']) == 2
    assert rec['page']['root'].startswith(str(tmp_path / 'remembered'))
    assert gui.velocities(rec['page']['files'], 'GJ 436', rec['page']['root'],
                          dace=True)['n'] == 45
    # only a result remembered, by its name
    with pytest.raises(ValueError):
        gui.recall('../lbl.csv')
    gui.unremember(out['id'])
    assert gui.remembered() == []

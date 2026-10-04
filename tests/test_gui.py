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
    # an instrument with its source in brackets is one name, with commas
    #   or spaces between the names
    assert gui.instrument_names('HARPS03, HIRES (CLS) HIRES  (Teklu+ 2025)'
                                ' NIRPS') == ['HARPS03', 'HIRES (CLS)',
                                              'HIRES (Teklu+ 2025)', 'NIRPS']
    assert gui.command('detailed', dict(target='x', file='a.rdb',
                                        exclude='HIRES (CLS), HARPS03'))[
        -3:] == ['--exclude', 'HIRES (CLS)', 'HARPS03']
    # the sources of VizieR: all by default, some named, none no VizieR
    assert gui.command('gather', dict(target='x', vz_talor19=False,
                                      vz_papers=False))[-5:] == [
        '--vizier-sources', 'teklu25', 'cls21', 'fischer14', 'rvbank20']
    off = {f'vz_{key}': False for key in ('teklu25', 'cls21', 'talor19',
                                          'fischer14', 'rvbank20', 'papers')}
    assert gui.command('gather', dict(target='x', **off)) == [
        'x', '--gather', 'archives', '--no-vizier']
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
    # the acceleration of the star, measured with its errors (none here)
    acc = res['acceleration']
    assert res['settings']['trend'] == 1 and set(acc) >= {'accel', 'tref'}
    assert acc['accel'][1] > 0 and abs(acc['accel'][0]) < 5 * acc['accel'][1]
    assert 'jerk' not in acc
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
    # a TESS sector in the archives too (not copied: the report fetches its
    #   own), listed in the manifest that is
    star = tmp_path / 'arch' / 'GJ_436'
    (star / 'phot' / 'tess').mkdir(parents=True)
    (star / 'phot' / 'tess' / 's0042_SPOC.csv').write_text(
        'rjd,flux,sflux\n60000.0,1.0,0.001\n')
    (star / 'manifest.json').write_text(json.dumps(dict(
        target='GJ 436', archives=dict(tess=dict(status='ok', sectors=[
            dict(sector=42, pipeline='SPOC',
                 file='phot/tess/s0042_SPOC.csv')])))))
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
    copied = tmp_path / 'remembered' / out['id'] / 'archives' / 'GJ_436'
    assert (copied / 'manifest.json').exists()
    assert not (copied / 'phot').exists()
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
    # a fold kept before folds carried their dates, BERV and solution: made
    #   again from the series on recall, and kept so
    qpath = tmp_path / 'remembered' / out['id'] / 'quick.json'
    kept = json.loads(qpath.read_text())
    nights = gui.selection(gui.page_options(page))[0].nightly()
    nn = int(np.sum(nights.inst == 'NIRPS'))
    kept['result']['folds'] = [dict(
        id=1, period=2.644, K=1.0, K_err=0.5, rms=3.0, tc=60000.0,
        curve=dict(phase=[0, 1], rv=[0, 0]),
        instruments=[dict(name='NIRPS', phase=[0.1] * nn, rv=[0.0] * nn,
                          err=[1.0] * nn, valid=[0.9] * (nn - 1) + [0.2])])]
    qpath.write_text(json.dumps(kept))
    rec = gui.recall(out['id'])
    folds = rec['quick']['result']['folds']
    nirps = [inst for inst in folds[0]['instruments']
             if inst['name'] == 'NIRPS'][0]
    assert 'model' in folds[0] and len(nirps['time']) == nn
    assert sum(val < 0.5 for val in nirps['valid']) == 1
    assert any('made again' in note for note in rec['notes'])
    assert gui.recall(out['id'])['notes'] == []
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


def test_the_quick_look_measures_the_acceleration(tmp_path, monkeypatch):
    """a drift in the velocities: the quick look gives it back (and its
    change with the curvature), its trend as the report's boxes say"""
    import time
    import numpy as np
    from koloa.simulate import simulate
    from koloa.gather import write_rv
    sim = simulate(err=1.5, seed=7, nvisits=60, per_visit=1,
                   baseline=1000)['data']
    years = (sim.time - sim.time.mean()) / 365.25
    sim.rv = sim.rv + 3.0 * years + 0.5 * years ** 2
    write_rv(sim, str(tmp_path / 'drift.csv'))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    assert gui.trend_order(dict(trend=False)) == 0
    assert gui.trend_order(dict(trend=True, curvature=True)) == 2
    assert gui.trend_order({}) == 1
    out = {}
    for curv in (False, True):
        state = gui.quick_fip(dict(files=[dict(path=str(tmp_path /
                                                         'drift.csv'))],
                                   trend=True, curvature=curv))
        for _ in range(300):
            state = gui.quick_state(state['id'])
            if state['status'] != 'running':
                break
            time.sleep(0.5)
        assert state['status'] == 'done', state['error']
        out[curv] = state['result']
    acc = out[True]['acceleration']
    assert out[True]['settings']['trend'] == 2
    # the offset of each instrument and the trend, for the series
    assert set(out[True]['offsets']) == {'inst'}
    assert len(out[True]['trend_model']['coefs']) == 2
    # d2v/dt2 = 2 x 0.5 m/s/yr^2; dv/dt = 3 m/s/yr at the middle
    assert abs(acc['jerk'][0] - 1.0) < 4 * acc['jerk'][1]
    assert abs(acc['accel'][0] - (3.0 + 1.0 * (acc['tref'] - sim.time.mean())
                                  / 365.25)) < 4 * acc['accel'][1]
    assert acc['accel_sigma'] > 5
    assert 'jerk' not in out[False]['acceleration']


def test_the_fold_carries_its_solution_dates_and_berv(tmp_path):
    """a fold gives each night its time and BERV, and its solution (offsets,
    trend, sinusoid) that the series is drawn with; the BERV of archives
    comes from their raw files; the PDF colours by date or BERV, draws the
    solution on the series, and marks what is off scale"""
    import numpy as np
    from koloa.data import RVData
    from koloa.gather import archive_berv
    rng = np.random.default_rng(5)
    time = np.sort(rng.uniform(60000, 60400, 50))
    inst = np.array(['AAA'] * 30 + ['BBB'] * 20)
    rv = (8.0 * np.sin(2 * np.pi * time / 5.3) + 2.0 * (time - 60200)
          / 365.25 + np.where(inst == 'AAA', 100.0, -50.0)
          + rng.normal(0, 1.0, 50))
    berv = 25 * np.cos(2 * np.pi * (time - 60000) / 365.25)
    data = RVData(time, rv, np.full(50, 1.0), inst=inst,
                  meta=dict(BERV=np.where(inst == 'AAA', berv, np.nan)))
    out = gui.fold(data, 5.3)
    one = {item['name']: item for item in out['instruments']}
    assert one['AAA']['berv'] is not None and one['BBB']['berv'] is None
    assert len(one['AAA']['time']) == 30
    # the solution, at the times of the points, is the fit itself
    for name in ('AAA', 'BBB'):
        sel = data.inst == name
        mod = gui.model_at(out, name, data.time[sel], 5.3)
        assert np.std(data.rv[sel] - mod) < 1.5
    assert gui.model_at(out, 'CCC', data.time, 5.3) is None
    # the BERV of the archives: DACE's raw file by rjd, CARMENES' by bjd
    (tmp_path / 'rv' / 'dace').mkdir(parents=True)
    (tmp_path / 'rv' / 'carmenes').mkdir()
    (tmp_path / 'rv' / 'dace' / 'raw_X.csv').write_text(
        'rjd,cal_berv\n60001.5,12.5\n60003.25,-3.0\n')
    (tmp_path / 'rv' / 'carmenes' / 'raw.csv').write_text(
        'bjd,berv\n2460010.0,7.0\n')
    got = archive_berv(str(tmp_path), np.array([60003.25, 60001.5, 60010.0,
                                                60020.0]))
    assert got[:3].tolist() == [-3.0, 12.5, 7.0] and np.isnan(got[3])
    # the PDF, coloured by date then by BERV, the solution on the series
    write = tmp_path / 'series.csv'
    write.write_text('rjd,vrad,svrad,BERV\n' + ''.join(
        f'{tt!r},{vv!r},1.0,{bb!r}\n' for tt, vv, bb in
        zip(time[inst == 'AAA'], rv[inst == 'AAA'], berv[inst == 'AAA'])))
    nights = gui.selection(dict(files=[dict(path=str(write))]))[0].nightly()
    fold = dict(gui.fold(nights, 5.3), id=1)
    qid = 'foldpdf'
    gui.QUICKS[qid] = dict(id=qid, status='done', start=0.0, end=1.0,
                           each=[], result=dict(
                               gui._fip_curves(type('R', (), dict(
                                   freq=np.array([0.1, 0.2, 0.3]),
                                   fip=np.array([0.5, 1e-5, 0.9]),
                                   family=np.array([0.5, 1e-5, 0.9])))()),
                               known=[], window=gui.WINDOW, passes=1,
                               planets=[], peak_list=[], folds=[fold],
                               pk=[0.1, 0.9], peaks=[], n=nights.n,
                               nexp=nights.n, instruments={'inst': nights.n},
                               inflation={}, settings=dict(gui.QUICK,
                                                           trend=1)))
    for colour in ('date', 'berv'):
        pdf = gui.quicklook_pdf(dict(files=[dict(path=str(write))]),
                                yr=[-5.0, 5.0], qid=qid, fold_colour=colour,
                                overlay=1, series_colour=colour)
        assert pdf[:4] == b'%PDF'
    # the series sent to the page: the BERV of each point
    shown = gui.velocities([dict(path=str(write))])['instruments'][0]
    assert len(shown['berv']) == shown['n'] and shown['berv'][0] is not None


def test_a_period_refined_within_its_peak():
    """a sinusoid over 20 years folded at a grid point of its peak (off by
    a tenth of a resolution element, a tenth of a cycle over the series):
    refined, its period found again, within its error; as asked, left"""
    import numpy as np
    from koloa.data import RVData
    rng = np.random.default_rng(8)
    true = 2.64390
    time = np.sort(rng.uniform(52000, 59300, 300))
    rv = 15.0 * np.sin(2 * np.pi * time / true) + rng.normal(0, 3.0, 300)
    data = RVData(time, rv, np.full(300, 3.0))
    asked = 1.0 / (1.0 / true - 0.1 / np.ptp(time))
    fixed = gui.fold(data, asked)
    assert fixed['period'] == asked and fixed['P_err'] is None
    refined = gui.fold(data, asked, refine=True)
    assert refined['asked'] == asked
    assert abs(refined['period'] - true) < 4 * refined['P_err']
    assert refined['P_err'] < 0.2 * abs(asked - true)
    # its phase no longer drifts: a smaller scatter about the sinusoid
    assert refined['rms'] < 3.2 < fixed['rms']


def test_a_fold_asked_its_keplerian_and_the_residuals(tmp_path, monkeypatch):
    """a fold at a period clicked (moved to the dip nearest), its Keplerian
    orbit, and the FIP of the residuals once its signal is subtracted: the
    second planet comes first"""
    import time
    import numpy as np
    from koloa.simulate import simulate
    from koloa.gather import write_rv
    sim = simulate(planets=[dict(P=5.3, K=9.0, e=0.0, tp=60000.0),
                            dict(P=13.7, K=4.0, e=0.0, tp=60003.0)],
                   err=1.2, seed=11, nvisits=90, per_visit=1,
                   baseline=400)['data']
    write_rv(sim, str(tmp_path / 'two.csv'))
    # two signals: the second planet is in the FIP too
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=2, nsweep=200, nburn=100))
    opts = dict(files=[dict(path=str(tmp_path / 'two.csv'))])

    def done(state):
        for _ in range(400):
            state = gui.quick_state(state['id'])
            if state['status'] != 'running':
                break
            time.sleep(0.5)
        assert state['status'] == 'done', state['error']
        return state
    state = done(gui.quick_fip(opts))
    qid = state['id']
    ids = [item['id'] for item in state['result']['folds']]
    # a click near the first peak: the fold already there (#1)
    near = gui.fold_request(qid, opts, period=5.32, snap=True)['fold']
    assert near['id'] == 1 and abs(near['period'] / 5.3 - 1) < 0.003
    # a period typed: a fold of its own, numbered next
    got = gui.fold_request(qid, opts, period=5.3)['fold']
    assert got['forced'] and got['id'] == max(ids) + 1
    assert got['period'] == 5.3
    # the Keplerian orbit of that fold, kept with it
    kep = gui.fold_request(qid, opts, fid=got['id'], kind='kepler')['fold']
    orbit = kep['kepler']
    assert orbit['kind'] == 'kepler' and orbit['e'] < 0.3
    assert abs(orbit['K'] - 9.0) < 4 * orbit['K_err']
    # its BIC: far better than no planet, a circular orbit as good (the
    #   simulated one is circular: e and omega do not earn their place)
    assert orbit['bic']['d_none'] > 50
    assert orbit['bic']['d_circular'] < 6
    assert abs(orbit['period'] - 5.3) < 0.01
    # the 1-sigma envelope of each fit, and draws of its solution
    assert len(got['curve']['lo']) == len(got['curve']['rv'])
    assert all(lo <= mid + 1e-9 <= hi + 2e-9 for lo, mid, hi in zip(
        got['curve']['lo'], got['curve']['rv'], got['curve']['hi']))
    assert len(got['draws']) == gui.NDRAW and len(orbit['draws']) > 30
    assert len(orbit['curve']['lo']) == len(orbit['curve']['rv'])
    # the envelope of the orbit: about the best one (the spread of draws
    #   of its full covariance), the best one always in it
    assert all(lo <= mid + 1e-9 <= hi + 2e-9 for lo, mid, hi in zip(
        orbit['curve']['lo'], orbit['curve']['rv'], orbit['curve']['hi']))
    grid = np.linspace(sim.time.min(), sim.time.max(), 500)
    inst = str(sim.instruments[0])
    best = gui.model_at(orbit, inst, grid, orbit['period'])
    many = np.array([gui.model_at(dict(model=draw), inst, grid,
                                  orbit['period']) for draw in orbit['draws']])
    low, high = gui.envelope(best, many)
    assert np.all((low <= best + 1e-9) & (best <= high + 1e-9))
    assert 0 < np.median(high - low) < 5
    # a known planet: at its published period, whatever the FIP says, with
    #   its published orbit
    gui.QUICKS[qid]['result']['known'] = [dict(
        name='Sim b', P=5.3, K=9.0, e=0.0, omega=0.0, tp=60000.0, tc=None,
        reference='a test')]
    kn = gui.fold_request(qid, opts, known='Sim b', kind='kepler')['fold']
    assert kn['period'] == 5.3 and kn['known']['name'] == 'Sim b'
    assert abs(max(kn['published']['rv']) - 9.0) < 0.01
    assert kn['kepler']['published']['model']['K'] == 9.0
    pub = kn['published']['model']
    sel = sim.inst == sim.instruments[0]
    # the published orbit is the simulated one: the points about it
    resid = sim.rv[sel] - gui.model_at(dict(model=pub), str(sim.instruments[0]),
                                       sim.time[sel], 5.3)
    assert np.std(resid - np.median(resid)) < 7
    # a click on the side of a peak folds at its top; again, the same fold
    side = gui.fold_request(qid, opts, period=13.45, snap=True)['fold']
    assert abs(side['period'] / 13.7 - 1) < 0.003
    again = gui.fold_request(qid, opts, period=13.47, snap=True)['fold']
    assert again['id'] in (side['id'], 2)
    # a signal ticked on the page: the other folds are fitted on the series
    #   without it (the 13.7 d one without the 5.3 d one: far less scatter,
    #   its K), again only when what they are without changes
    first = next(item for item in state['result']['folds']
                 if abs(item['period'] / 5.3 - 1) < 0.01)
    second = gui.fold_request(qid, opts, period=13.7)['fold']
    raw = second['rms']
    clean = gui.fold_request(qid, opts, fid=second['id'],
                             minus=[first['model']], minus_tag='1s')['fold']
    assert clean['id'] == second['id'] and clean['minus_tag'] == '1s'
    assert clean['rms'] < 0.6 * raw
    assert abs(clean['K'] - 4.0) < 4 * clean['K_err']
    kept = clean['rms']
    assert gui.fold_request(qid, opts, fid=second['id'],
                            minus=[first['model']],
                            minus_tag='1s')['fold']['rms'] == kept
    back = gui.fold_request(qid, opts, fid=second['id'], minus=[],
                            minus_tag='')['fold']
    assert back['minus_tag'] == '' and abs(back['rms'] - raw) < 1e-9
    # the FIP of what is left, beside the FIP of the series: no FIP of
    #   each instrument, and no other quick FIP stopped for it
    other = gui.quick_fip(opts)
    stage = gui.quick_fip(dict(opts, subtract=[dict(first['model'],
                                                    label='#1')],
                               stage=True))
    assert gui.QUICKS[stage['id']]['stage']
    assert not gui.QUICKS[other['id']].get('cancel')
    left = done(stage)
    assert left['each'] == []
    assert abs(left['result']['peaks'][0]['period'] / 13.7 - 1) < 0.01
    done(other)
    # a quick FIP recalled has no series kept: the page's is read again
    gui._QUICK_DATA.pop(qid)
    assert gui.fold_request(qid, opts, period=13.7)['fold']['period'] == 13.7
    # the residuals: the 5.3 d signal out, the 13.7 d one first
    res = done(gui.quick_fip(dict(opts, subtract=[dict(got['model'],
                                                       label='#x')])))
    assert res['result']['subtracted'][0]['label'] == '#x'
    assert abs(res['result']['peaks'][0]['period'] / 13.7 - 1) < 0.01


def test_a_quick_fip_is_stopped(tmp_path, monkeypatch):
    """a long quick FIP stopped by the next one (out of date), and that one
    stopped as it runs: its chains ended, within seconds"""
    import time
    from koloa.simulate import simulate
    from koloa.gather import write_rv
    sim = simulate(planets=[dict(P=5.3, K=8.0, e=0.0)], err=1.5, seed=3,
                   nvisits=60, per_visit=1, baseline=300)['data']
    write_rv(sim, str(tmp_path / 'long.csv'))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=2, nsweep=50000, nburn=200))
    opts = dict(files=[dict(path=str(tmp_path / 'long.csv'))])

    def until(qid, status, wait=30.0):
        start = time.time()
        while time.time() - start < wait:
            state = gui.quick_state(qid)
            if status(state):
                return state
            time.sleep(0.25)
        raise AssertionError(gui.quick_state(qid))
    first = gui.quick_fip(opts)
    until(first['id'], lambda st: st['step'] == 'fip1')
    second = gui.quick_fip(opts)
    start = time.time()
    assert until(first['id'], lambda st: st['status'] != 'running'
                 )['status'] == 'stopped'
    until(second['id'], lambda st: st['step'] == 'fip1')
    gui.stop_quick(second['id'])
    assert until(second['id'], lambda st: st['status'] != 'running'
                 )['status'] == 'stopped'
    assert time.time() - start < 30


def test_a_fold_on_a_transit_ephemeris(tmp_path, monkeypatch):
    """a fold on a transit ephemeris: phase 0 at the transit carried to the
    velocities, its error from those of tc and P; the conjunction the
    velocities put on their own set against it; a wrong ephemeris shows"""
    import time
    import numpy as np
    from koloa.simulate import simulate
    from koloa.gather import write_rv
    # a transit some periods before the velocities (simulated from day 0)
    tc, per = -30.123, 3.21
    # a circular orbit: the star moves away fastest a quarter period
    #   before the transit (its velocity falls through zero at it)
    sim = simulate(planets=[dict(P=per, K=6.0, e=0.0, tp=tc - per / 4)],
                   err=1.0, seed=21, nvisits=50, per_visit=1,
                   baseline=200)['data']
    write_rv(sim, str(tmp_path / 'transit.csv'))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    opts = dict(files=[dict(path=str(tmp_path / 'transit.csv'))])
    state = gui.quick_fip(opts)
    for _ in range(300):
        state = gui.quick_state(state['id'])
        if state['status'] != 'running':
            break
        time.sleep(0.5)
    qid = state['id']
    good = dict(name='Sim b', source='archive', P=per, P_err=2e-5, tc=tc,
                tc_err=1e-3, reference='a test')
    bad = dict(good, name='Sim b (off)', tc=tc + 0.125)
    gui.QUICKS[qid]['result']['transits'] = [good, bad]
    fold = gui.fold_request(qid, opts, transit='Sim b')['fold']
    eph = fold['transit']
    ncyc = eph['cycles']
    assert abs(eph['t0'] - (tc + ncyc * per)) < 1e-9
    assert abs(eph['t0_err'] - np.hypot(1e-3, ncyc * 2e-5)) < 1e-12
    assert abs(fold['K'] - 6.0) < 4 * fold['K_err']
    assert eph['shift_sigma'] < 3
    # the same, asked again: the same fold
    assert gui.fold_request(qid, opts, transit='Sim b')['fold']['id'] == \
        fold['id']
    # an ephemeris 3 h off: the velocities say so
    off = gui.fold_request(qid, opts, transit='Sim b (off)')['fold']
    assert off['transit']['shift_sigma'] > 3
    assert abs(off['transit']['shift'] + 0.125) < 0.05
    # the Keplerian, the transits as priors, phase 0 at the transit
    kep = gui.fold_request(qid, opts, transit='Sim b',
                           kind='kepler')['fold']['kepler']
    assert kep['tc'] == eph['t0'] and abs(kep['K'] - 6.0) < 4 * kep['K_err']
    pdf = gui.quicklook_pdf(dict(opts, mstar=0.5, mstar_err=0.05), qid=qid,
                            fold_model='kepler', fip_view='joint')
    assert pdf[:4] == b'%PDF'
    # the LaTeX document itself (not the figures alone, its fall-back)
    import shutil
    fitz = pytest.importorskip('fitz')
    if shutil.which('pdflatex'):
        text = ' '.join(page.get_text() for page in
                        fitz.open(stream=pdf, filetype='pdf'))
        assert 'transit ephemeris' in text and 'Jup' in text


def test_a_batch_fip(tmp_path, monkeypatch):
    """the quick FIP of many files, one after the other: their best peaks
    and accelerations; a file that fails does not stop the others; one
    opened for the page; New target leaves the batch alone; a batch
    stopped"""
    import time
    from koloa.simulate import simulate
    from koloa.gather import write_rv
    for seed, period in ((3, 5.3), (4, 11.7)):
        sim = simulate(planets=[dict(P=period, K=9.0, e=0.0)], err=1.0,
                       seed=seed, nvisits=50, per_visit=1,
                       baseline=300)['data']
        write_rv(sim, str(tmp_path / f'star{seed}.csv'))
    (tmp_path / 'empty.csv').write_text('rjd,vrad,svrad\n')
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    assert len(gui.list_files(str(tmp_path), '*.csv')['paths']) == 3
    paths = [str(tmp_path / name) for name in ('star3.csv', 'empty.csv',
                                               'star4.csv')]
    state = gui.batch_fip(paths, dict(trend=True))
    gui.forget()   # New target: the batch goes on
    for _ in range(600):
        state = gui.batch_state(state['id'])
        if state['status'] != 'running':
            break
        time.sleep(0.5)
    assert state['status'] == 'done'
    items = state['items']
    assert [item['status'] for item in items] == ['done', 'failed', 'done']
    for item, period in ((items[0], 5.3), (items[2], 11.7)):
        assert abs(item['summary']['period'] / period - 1) < 0.01
        assert item['summary']['fip'] < 0.01
        assert item['summary']['accel_err'] > 0
    opened = gui.batch_open(state['id'], 2)
    assert opened['page']['files'][0]['path'] == paths[2]
    assert opened['quick']['status'] == 'done'
    assert opened['rv']['instruments'][0]['n'] == 50
    # a batch stopped: the file running and those after it
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=2, nsweep=50000, nburn=100))
    state = gui.batch_fip(paths, {})
    time.sleep(2.0)
    gui.stop_batch(state['id'])
    for _ in range(100):
        state = gui.batch_state(state['id'])
        if state['status'] != 'running':
            break
        time.sleep(0.25)
    assert state['status'] == 'stopped'
    assert all(item['status'] in ('stopped', 'failed')
               for item in state['items'])


def test_a_batch_fip_with_every_archive(tmp_path, monkeypatch):
    """a batch with the archives: the star of each file by its APERO name
    (its OBJECT column), its archives gathered once (read from disk the
    next time) and put with the file, the instruments of each line; a
    file whose star is not known run alone; opened with its star"""
    import time
    import numpy as np
    from koloa import gather as kgather
    from koloa.data import RVData
    from koloa.gather import write_rv
    from koloa.simulate import simulate
    pytest.importorskip('yaml')
    db = tmp_path / 'apero' / 'pending'
    db.mkdir(parents=True)
    (db / 'GL1.yaml').write_text('APERO_NAME: GL1\nSIMBAD_NAME: GJ    1\n'
                                 'STATUS: pending\nSPT:\n  value: M1.5V\n'
                                 'ALIASES:\n- Gl 1\n')
    monkeypatch.setenv('KOLOA_APERO_ASTROMETRICS', str(tmp_path / 'apero'))
    sim = simulate(planets=[dict(P=5.3, K=9.0, e=0.0)], err=1.0, seed=3,
                   nvisits=40, per_visit=1, baseline=300)['data']
    lines = [f'{tt!r},{rv!r},{er!r},Gl 1\n'
             for tt, rv, er in zip(sim.time, sim.rv, sim.err)]
    (tmp_path / 'gl1.csv').write_text('rjd,vrad,svrad,OBJECT\n'
                                      + ''.join(lines))
    write_rv(simulate(planets=[], err=1.0, seed=4, nvisits=20, per_visit=1,
                      baseline=300)['data'], str(tmp_path / 'other.csv'))
    gathered = []

    def gather(target, root, **kwargs):
        gathered.append(target)
        star = tmp_path / 'arch' / kgather.folder_name(target)
        rng = np.random.default_rng(5)
        tarch = np.sort(rng.uniform(sim.time.min() - 600,
                                    sim.time.min() - 100, 25))
        write_rv(RVData(tarch, 9.0 * np.sin(2 * np.pi * tarch / 5.3)
                        + rng.normal(0, 1.5, 25), np.full(25, 1.5),
                        inst=np.array(['HARPS15'] * 25)),
                 str(star / 'rv' / 'all_rv.csv'))
        (star / 'manifest.json').write_text(json.dumps(dict(
            target=target, archives=dict(dace=dict(status='ok')))))
    monkeypatch.setattr(kgather, 'gather', gather)
    monkeypatch.setattr(gui, 'known_periods', lambda target: [])
    monkeypatch.setattr(gui, 'transits_of', lambda target, known: [])
    monkeypatch.setattr(gui, '_star_planets', lambda target: (
        [dict(name='GJ 1 b', P=5.3)], [dict(toi='999.01', P=40.1,
                                             disposition='PC')]))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    # its TESS light curve (not asked of MAST): a transit at the
    #   conjunction of the signal (9 sin(2 pi t / 5.3): falling through zero
    #   at t = 2.65 + 5.3 n)
    rng = np.random.default_rng(6)
    ltime = np.arange(60100.0, 60127.0, 2.0 / 1440)
    lflux = rng.normal(0, 0.5, len(ltime))
    phase = ((ltime - 2.65) / 5.3 + 0.5) % 1.0 - 0.5
    lflux[np.abs(phase * 5.3 * 24) < 0.8] -= 1.5
    light = dict(time=ltime, flux=lflux, err=np.full(len(ltime), 0.5),
                 sector=np.full(len(ltime), 40.0))
    monkeypatch.setattr(gui, 'tess_light', lambda target, root='',
                        fetch=False: (light, 'a test'))
    monkeypatch.setattr(gui, 'resolve_star', lambda target, root='',
                        refresh=False: dict(star=dict(
                            mass=0.4, radius=0.4, radius_source='a test')))
    paths = [str(tmp_path / 'gl1.csv'), str(tmp_path / 'other.csv')]

    def run():
        state = gui.batch_fip(paths, dict(trend=True), archives=True,
                              root=str(tmp_path / 'arch'))
        for _ in range(600):
            state = gui.batch_state(state['id'])
            if state['status'] != 'running':
                return state
            time.sleep(0.5)
        return state
    state = run()
    assert state['status'] == 'done' and state['archives']
    first, other = state['items']
    assert first['star']['apero'] == 'GL1' and first['star']['raw'] == 'Gl 1'
    assert first['star']['target'] == 'GJ 1' and first['star']['spt'] == \
        'M1.5V'
    assert first['star']['planets'][0]['P'] == 5.3
    assert first['star']['tois'][0]['toi'] == '999.01'
    insts = first['summary']['instruments']
    assert len(insts) == 2 and insts['HARPS15'] == 25
    assert first['summary']['sources']['HARPS15'] == 'DACE'
    assert [first['summary']['sources'][name] for name in insts
            if name != 'HARPS15'] == ['file: gl1.csv']
    assert abs(first['summary']['period'] / 5.3 - 1) < 0.01
    # the transit in TESS at the best peak: flagged
    assert first['summary']['transit']['plausible'], \
        first['summary']['transit']['why']
    assert first['summary']['transit']['depth'] == pytest.approx(1.5,
                                                                 abs=0.4)
    # a file with no star: alone, and why
    assert other['status'] == 'done' and other['star']['apero'] is None
    assert 'HARPS15' not in other['summary']['instruments']
    assert gathered == ['GJ 1']
    opened = gui.batch_open(state['id'], 0)
    assert opened['page']['target'] == 'GJ 1'
    assert opened['page']['detailed']['dace'] is True
    assert {inst['name'] for inst in opened['rv']['instruments']} >= {
        'HARPS15'}
    # again: the archives on disk, not gathered anew
    state = run()
    assert state['status'] == 'done' and gathered == ['GJ 1']
    # the star of a file for the page
    star = gui.star_of_file(paths[0])
    assert star['apero'] == 'GL1' and star['target'] == 'GJ 1'
    assert gui.star_of_file(str(tmp_path / 'none.csv'))['target'] == ''


def test_the_zero_of_an_instrument_on_the_series():
    """an instrument with a few discrepant nights and a few nights of many
    exposures: its fitted offset lies between the two (the fit weighs
    nights), its zero on the plot puts its exposures about the model"""
    import numpy as np
    from koloa.data import RVData
    from koloa.fit import RVModel
    rng = np.random.default_rng(31)
    ta = np.sort(rng.uniform(60000, 61000, 60))
    tb = np.concatenate([60300.1 + np.arange(4),
                         np.concatenate([60700.0 + night + np.linspace(
                             0, 0.25, 30) for night in range(4)])])
    slope = 3.0 / 365.25
    rva = 10.0 + slope * (ta - 60500) + rng.normal(0, 1.0, len(ta))
    rvb = -20.0 + slope * (tb - 60500) + rng.normal(0, 1.0, len(tb))
    rvb[:4] -= 15.0   # four discrepant nights
    data = RVData(np.concatenate([ta, tb]), np.concatenate([rva, rvb]),
                  np.ones(len(ta) + len(tb)),
                  inst=np.array(['A'] * len(ta) + ['B'] * len(tb)))
    nights = data.nightly()
    fit = RVModel(nights, [], likelihood='mixture', unit='both', trend=1,
                  seq_jitter='instrument').fit(nstart=2, quiet=True)
    got = gui._levels(fit, data)
    model = got['trend_model']
    sel = data.inst == 'B'
    drift = model['coefs'][0] * (data.time[sel] - model['tref']) / \
        model['tscale']
    about_zero = np.median(data.rv[sel] - got['offsets']['B'] - drift)
    about_fit = np.median(data.rv[sel] - got['fit_offsets']['B'] - drift)
    assert abs(about_zero) < 0.5 < abs(about_fit)
    assert got['trend_model']['signals'] == []
    # draws of the trend for its envelope; the solution the series is
    #   drawn about, the same for both instruments, each about it
    assert len(model['draws']) > 5
    curve = gui.series_curve(dict(trend_model=model))
    for inst in ('A', 'B'):
        sel = data.inst == inst
        rest = data.rv[sel] - curve['at'](data.time[sel])
        assert abs(np.median(rest - np.median(rest))) < 1e-9
    assert np.isinf(curve['shortest']) and len(curve['draws']) > 5
    assert gui.series_curve(None) is None


def test_the_series_about_one_solution(tmp_path):
    """two instruments far apart, a trend and a planet: the series drawn
    about one solution (the fold shown, else the fit of the quick look),
    the same for both, each instrument about it (the median of its
    velocities minus it, zero), the residuals under it, in the PDF too"""
    import matplotlib.pyplot as plt
    import numpy as np
    rng = np.random.default_rng(9)

    def write(name, t0, t1, offset):
        tt = np.sort(rng.uniform(t0, t1, 50))
        rv = offset + 3.0 * (tt - 60000) / 365.25 + 8.0 * np.sin(
            2 * np.pi * tt / 7.3) + rng.normal(0, 1.5, len(tt))
        (tmp_path / name).write_text('rjd,vrad,svrad\n' + ''.join(
            f'{a!r},{b!r},1.5\n' for a, b in zip(tt, rv)))
    write('a.csv', 59000, 60000, 120.0)
    write('b.csv', 59700, 60900, -45.0)
    opts = dict(files=[dict(path=str(tmp_path / 'a.csv'), label='A'),
                       dict(path=str(tmp_path / 'b.csv'), label='B')])
    data, _, _ = gui.selection(opts)
    shown = gui.fold(data.nightly(), 7.3, None, 1)
    quick = dict(folds=[dict(shown, id=1)])
    curve = gui.series_curve(quick, overlay=1)
    assert curve['shortest'] == pytest.approx(7.3, rel=1e-3)
    rests = []
    for inst in data.instruments:
        sel = data.inst == inst
        rest = data.rv[sel] - curve['at'](data.time[sel])
        rests.append(rest - np.median(rest))
    # one solution for both: their residuals about it are the noise
    assert np.std(np.concatenate(rests)) < 2.5
    # the trend alone (no signal) and the draws of the envelope
    assert len(curve['draws']) > 5
    flat = curve['at'](np.array([59000.0, 60900.0]), False)
    assert flat[1] - flat[0] == pytest.approx(3.0 * 1900 / 365.25, abs=2.5)
    colour = dict(A='C0', B='C1')
    marker = dict(A='o', B='s')
    fig = gui._series_figure(data, {}, quick, None, None, 'test', colour,
                             marker, overlay=1)
    assert len(fig.axes) == 2   # the series, the residuals
    assert fig.axes[1].get_ylabel().startswith('residuals')
    plt.close(fig)
    # without the fold and without a fit: each about its median, alone
    fig = gui._series_figure(data, {}, None, None, None, 'test', colour,
                             marker)
    assert len(fig.axes) == 1   # the series
    plt.close(fig)

"""
koloa.survey, without the network: a spectral type as a number, a sample
from SIMBAD's answer, what the archives have of a star, the files put
with the stars, a batch folder packed and run.
"""
import ast
import json
import os
import tarfile

import numpy as np
import pytest

from koloa import survey


def test_a_spectral_type_as_a_number():
    assert survey.sp_number('M3.5V') == 63.5
    assert survey.sp_number('M2.0V+M0.5V') == 62.0
    assert survey.sp_number('dM4e') == 64.0
    assert survey.sp_number('K7V') == 57.0
    assert survey.sp_number('M') == 65.0
    assert survey.sp_number('L1') == 71.0
    assert survey.sp_number('DA3') is None and survey.sp_number('') is None
    assert survey.sp_kind('M3.5V') == 'dwarf'
    assert survey.sp_kind('dM4e') == 'dwarf'
    assert survey.sp_kind('M4') == 'dwarf'
    assert survey.sp_kind('sdM1') == 'subdwarf'
    assert survey.sp_kind('M1VIp') == 'subdwarf'
    assert survey.sp_kind('M3III') == 'giant'
    assert survey.sp_kind('G8IV') == 'subgiant'
    assert survey.sp_kind('DA3') == 'white dwarf'
    assert survey.name_key('NAME Proxima Centauri') == 'PROXIMACENTAURI'
    assert survey.name_key('Gl 876') == survey.name_key('GJ  876') == 'GJ876'


def _simbad(monkeypatch):
    """SIMBAD's answer to the three questions of a sample"""
    from koloa import gather
    asked = []

    def tap(url, query, timeout=120.0):
        asked.append(query)
        if 'FROM ident' in query:
            return [dict(oidref='1', id='GJ 876'), dict(oidref='1',
                                                        id='HIP 113020'),
                    dict(oidref='1', id='Karmn J22532-142'),
                    dict(oidref='2', id='GJ   65 B'),
                    dict(oidref='3', id='GJ 65'),
                    dict(oidref='5', id='HD 1')]
        if 'FROM mesVar' in query:  # the same paper twice: once is enough
            return [dict(oidref='1', period='87.3',
                         bibcode='2015MNRAS.452.2745S'),
                    dict(oidref='1', period='82.8',
                         bibcode='2023A&A...672A..52F'),
                    dict(oidref='1', period='82.8',
                         bibcode='2023A&A...672A..52F'),
                    dict(oidref='2', period='', bibcode='2016ApJ...821...93N')]
        return [
            dict(oid='1', main_id='BD-15  6290', ra='343.3197', dec='-14.2637',
                 plx_value='214.0', plx_err='0.05', sp_type='M3.5V',
                 otype='PM*', rvz_radvel='-1.5', V='10.19', G='8.9', J='5.9',
                 K='5.0'),
            dict(oid='2', main_id='G 272-61B', ra='24.7566', dec='-17.9506',
                 plx_value='373.8', plx_err='0.3', sp_type='M6V', otype='PM*',
                 rvz_radvel='', V='', G='10.8', J='', K=''),
            dict(oid='3', main_id='G 272-61', ra='24.7560', dec='-17.9500',
                 plx_value='373.7', plx_err='2.7', sp_type='M5.5V+M6V',
                 otype='**', rvz_radvel='', V='12.08', G='', J='6.3', K='5.3'),
            # a giant, a white dwarf and an earlier type: not in the sample
            dict(oid='4', main_id='X Giant', ra='10.0', dec='10.0',
                 plx_value='80.0', plx_err='1', sp_type='M3III', otype='*',
                 rvz_radvel='', V='5', G='', J='', K=''),
            dict(oid='5', main_id='HD 1', ra='20.0', dec='20.0',
                 plx_value='90.0', plx_err='1', sp_type='K7V', otype='*',
                 rvz_radvel='', V='8', G='', J='', K=''),
        ]
    monkeypatch.setattr(gather, '_tap', tap)
    return asked


def test_a_sample_from_simbad(monkeypatch):
    asked = _simbad(monkeypatch)
    stars = survey.sample(('M0', 'M9'), 15.0)
    assert "b.plx_value >= 66.666667" in asked[0]
    assert "b.sp_type LIKE 'M%'" in asked[0] and "LIKE 'dM%'" in asked[0]
    # the nearest first; the giant and the K dwarf are not M dwarfs
    assert [star['name'] for star in stars] == ['GJ 65 B', 'GJ 65', 'GJ 876']
    star = stars[-1]
    assert star['main'] == 'BD-15 6290' and star['spnum'] == 63.5
    assert star['distance'] == pytest.approx(1000 / 214.0)
    assert star['ids']['HIP'] == 'HIP 113020' and star['V'] == 10.19
    assert star['near'] is None
    # the rotation periods SIMBAD lists, the latest paper first
    assert "v.vartyp = 'ROT'" in asked[2] and 'plx_value >= 66.6' in asked[2]
    assert star['rotation'] == [
        dict(period=82.8, source='2023A&A...672A..52F'),
        dict(period=87.3, source='2015MNRAS.452.2745S')]
    assert stars[0]['rotation'] == [] and stars[1]['rotation'] == []
    # a system and its component: each told of the other
    assert stars[0]['near'] == 'GJ 65' and stars[1]['near'] == 'GJ 65 B'
    assert stars[0]['V'] is None and stars[0]['G'] == 10.8
    # a brighter limit (G for a star with no V), the giants too, K stars
    assert [star['name'] for star in survey.sample(('M0', 'M9'), 15.0,
                                                   vmax=10.5)] == ['GJ 876']
    assert len(survey.sample(('M0', 'M9'), 15.0, dwarfs=False)) == 4
    assert [star['name'] for star in survey.sample(('K5', 'K9'), 15.0)] \
        == ['HD 1']
    assert survey.ident_of(star)['gj'] == 'GJ 876'
    with pytest.raises(ValueError):
        survey.sample(('X', 'M9'), 15.0)


def test_what_the_archives_have_of_a_star(tmp_path, monkeypatch):
    """CARMENES and the surveys by position, DACE by asking: kept where
    koloa.gather reads it back"""
    from koloa import gather, published
    from koloa.data import RVData
    _simbad(monkeypatch)
    stars = survey.sample(('M0', 'M9'), 15.0)
    monkeypatch.setattr(gather, 'carmenes_objects', lambda refresh=False: [
        dict(ra='343.3197', dec='-14.2637', nobs='69', p_rot='81.0',
             p_rot_source='DA19', carmenes_id='J22532-142')])
    monkeypatch.setattr(published, 'survey_stars',
                        lambda survey, refresh=False, timeout=120.0: (
                            [('GJ876', 343.3197, -14.2637)]
                            if survey['key'] in ('teklu25', 'rvbank20')
                            else []))
    folders = []

    def dace_rv(ident, target, folder, api_key=None, refresh=False):
        folders.append(folder)
        if ident['gj'] == 'GJ 65':  # known, with no velocity: nothing
            return RVData(np.array([]), np.array([]), np.array([]),
                          inst=np.array([], dtype=str))
        if ident['gj'] != 'GJ 876':
            return None
        return RVData(np.arange(12.0), np.zeros(12), np.ones(12),
                      inst=np.array(['HARPS03'] * 8 + ['HARPS15'] * 4))
    monkeypatch.setattr(gather, 'dace_rv', dace_rv)
    seen = []
    survey.check(stars, str(tmp_path / 'arch'), workers=2,
                 progress=lambda done, total: seen.append((done, total)))
    assert sorted(seen) == [(1, 3), (2, 3), (3, 3)]
    arch = stars[-1]['archives']
    assert arch['dace'] == dict(n=12, instruments=dict(HARPS03=8, HARPS15=4))
    assert arch['carmenes'] == 69 and arch['surveys'] == ['teklu25',
                                                         'rvbank20']
    assert arch['n'] == 81
    # the rotation period of CARMENES DR1 after those of SIMBAD, once
    assert stars[-1]['rotation'][-1] == dict(period=81.0,
                                             source='CARMENES DR1 (DA19)')
    survey.check_star(stars[-1], str(tmp_path / 'arch'))
    assert len(stars[-1]['rotation']) == 3
    assert os.path.join('arch', 'GJ_876', 'rv', 'dace') in ''.join(folders)
    assert survey.has_data(stars[-1]) and not survey.has_data(stars[0])
    assert stars[0]['archives']['dace'] is None
    assert stars[1]['archives']['dace'] is None
    assert not survey.has_data(stars[1])
    # spectrograph by spectrograph: the eras of HARPS and its release in
    #   the RVBank are one line, with the velocities DACE counted
    assert stars[-1]['summary'] == [
        dict(name='CARMENES', n=69, where=['DR1']),
        dict(name='HARPS', n=12, where=['DACE', 'RVBank']),
        dict(name='HIRES', n=None, where=['Teklu+ 2025'])]
    assert stars[0]['summary'] == []
    over = survey.overview(stars)
    assert (over['n'], over['checked'], over['data'], over['none']) == (
        3, 3, 1, 2)
    assert over['velocities'] == 81 and over['several'] == 1
    assert over['spectrographs'][0] == dict(name='CARMENES', stars=1,
                                            velocities=69, where=dict(DR1=1))
    assert survey.spectrographs(dict(dace=dict(n=7, instruments=dict(
        HARPN=4, CORALIE14=3)), carmenes=None, surveys=['cls21', 'talor19'])
        ) == [dict(name='HARPS-N', n=4, where=['DACE']),
              dict(name='CORALIE', n=3, where=['DACE']),
              dict(name='HIRES', n=None, where=['CLS', 'Tal-Or+ 2019'])]


def test_the_files_put_with_the_stars(tmp_path, monkeypatch):
    from koloa import apero_names
    _simbad(monkeypatch)
    stars = survey.sample(('M0', 'M9'), 15.0)
    for name in ('lbl_GL876_GL876.rdb', 'lbl_GL876_other.rdb',
                 'lbl_HD40307_HD40307.rdb', 'notes.txt'):
        (tmp_path / name).write_text('rjd\tvrad\tsvrad\n1\t2\t3\n')

    def star_of_file(path, fetch=True):
        base = os.path.basename(path)
        raw = 'GL876' if 'GL876' in base else 'HD40307'
        return dict(raw=raw, source='file name', apero=raw, entry=None,
                    target='BD-15 6290' if raw == 'GL876' else 'HD 40307')
    monkeypatch.setattr(apero_names, 'star_of_file', star_of_file)
    res = survey.match_files(stars, [str(tmp_path)], '*.rdb')
    assert [os.path.basename(path) for path in res['matched']['GJ 876']] == [
        'lbl_GL876_GL876.rdb', 'lbl_GL876_other.rdb']
    assert [row['target'] for row in res['unmatched']] == ['HD 40307']
    assert len(stars[-1]['files']) == 2 and stars[0]['files'] == []
    assert survey.has_data(stars[-1])
    with pytest.raises(ValueError, match='no folder'):
        survey.match_files(stars, [str(tmp_path / 'nowhere')])


def _batch(tmp_path, monkeypatch, nstar=3):
    """a batch folder packed from files and archives made here"""
    from koloa.data import RVData
    from koloa.gather import write_rv
    from koloa.simulate import simulate
    monkeypatch.setenv('KOLOA_CACHE', str(tmp_path / 'cache_here'))
    (tmp_path / 'cache_here' / 'archive').mkdir(parents=True)
    (tmp_path / 'cache_here' / 'archive' / 'kept.json').write_text('{}')
    (tmp_path / 'cache_here' / 'carmenes_objects.json').write_text('[]')
    stars = []
    for it in range(nstar):
        sim = simulate(planets=[dict(P=5.3 + it, K=9.0, e=0.0)], err=1.0,
                       seed=3 + it, nvisits=40, per_visit=1,
                       baseline=300)['data']
        path = tmp_path / 'lbl' / f'star{it}.csv'
        path.parent.mkdir(exist_ok=True)
        write_rv(sim, str(path))
        stars.append(dict(name=f'GJ {it + 1}', main=f'Main {it + 1}',
                          sptype='M3V', distance=5.0 + it, ra=10.0 * it,
                          dec=-5.0, files=[str(path)] if it != 1 else []))
    # the second star has no file: its archives alone
    rng = np.random.default_rng(5)
    time = np.sort(rng.uniform(58000, 58400, 30))
    write_rv(RVData(time, 8.0 * np.sin(2 * np.pi * time / 7.7)
                    + rng.normal(0, 1.0, 30), np.full(30, 1.0),
                    inst=np.array(['HARPS15'] * 30)),
             str(tmp_path / 'arch' / 'GJ_2' / 'rv' / 'all_rv.csv'))
    (tmp_path / 'arch' / 'GJ_2' / 'manifest.json').write_text('{}')
    return stars


def test_a_batch_folder_packed(tmp_path, monkeypatch):
    """the folder and its tar: the files, the archives, the cache, koloa,
    the script with its ROOT, nothing that names a path of this machine"""
    stars = _batch(tmp_path, monkeypatch)
    told = []
    made = survey.pack(stars, 'm dwarfs/15pc', out=str(tmp_path / 'out'),
                       root=str(tmp_path / 'arch'), gather=False,
                       server_root='/scratch/me/m_dwarfs_15pc', jobs=4,
                       progress=lambda done, total, name: told.append(done))
    folder = made['folder']
    assert os.path.basename(folder) == 'm_dwarfs_15pc'
    assert made['n'] == 3 and made['files'] == 2 and made['missing'] == []
    assert told[-1] == 3
    for part in ('run_batch.py', 'submit.sh', 'README.txt', 'targets.json',
                 'files/GJ_1/star0.csv', 'files/GJ_3/star2.csv',
                 'archives/GJ_2/rv/all_rv.csv', 'cache/archive/kept.json',
                 'cache/carmenes_objects.json', 'koloa_src/koloa/survey.py',
                 'koloa_src/koloa/gui_static/gui.js'):
        assert os.path.exists(os.path.join(folder, part)), part
    script = open(os.path.join(folder, 'run_batch.py')).read()
    ast.parse(script)
    assert "ROOT = '/scratch/me/m_dwarfs_15pc'" in script
    assert 'JOBS = 4' in script and "if __name__ == '__main__':" in script
    # no path of this machine in what the batch reads
    held = survey.targets_of(folder)
    assert [star['files'] for star in held['targets']] == [
        [os.path.join('files', 'GJ_1', 'star0.csv')], [],
        [os.path.join('files', 'GJ_3', 'star2.csv')]]
    assert str(tmp_path) not in json.dumps(held) + script
    assert '--array=0-2%4' in open(os.path.join(folder, 'submit.sh')).read()
    with tarfile.open(made['tar']) as handle:
        names = handle.getnames()
    assert 'm_dwarfs_15pc/run_batch.py' in names
    assert 'm_dwarfs_15pc/files/GJ_1/star0.csv' in names
    assert made['size'] == os.path.getsize(made['tar'])
    # a star with nothing is said
    made = survey.pack(stars + [dict(name='GJ 9', files=[])], 'again',
                       out=str(tmp_path / 'out'),
                       root=str(tmp_path / 'arch'), gather=False, tar=False)
    assert made['missing'] == ['GJ 9'] and made['tar'] is None
    with pytest.raises(ValueError, match='not a batch folder'):
        survey.targets_of(str(tmp_path))


def test_a_batch_run_where_it_was_carried(tmp_path, monkeypatch):
    """the folder moved elsewhere (no path of where it was made): each
    star's quick FIP from its files or its archives, kept as it goes, the
    stars done not done again, a table of them all"""
    from koloa import gui
    pytest.importorskip('yaml')
    stars = _batch(tmp_path, monkeypatch, nstar=2)
    made = survey.pack(stars, 'two', out=str(tmp_path / 'out'),
                       root=str(tmp_path / 'arch'), gather=False, tar=False)
    moved = str(tmp_path / 'elsewhere' / 'two')
    os.makedirs(os.path.dirname(moved))
    os.rename(made['folder'], moved)
    for name in ('lbl', 'arch'):
        os.rename(str(tmp_path / name), str(tmp_path / (name + '_gone')))
    monkeypatch.setattr(gui, 'known_periods', lambda target: [])
    monkeypatch.setattr(gui, 'transits_of', lambda target, known: [])
    monkeypatch.setattr(gui, '_star_planets', lambda target: ([], []))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    monkeypatch.setattr(gui, 'tess_light',
                        lambda target, root='', fetch=False: (None, 'none'))
    monkeypatch.setattr(gui, 'space_light', lambda target, root='',
                        fetch=False, mission='tess': (None, 'none', {}))
    rows = survey.run(moved)
    assert [row['name'] for row in rows] == ['GJ 1', 'GJ 2']
    assert [row['status'] for row in rows] == ['done', 'done']
    assert rows[0]['period'] == pytest.approx(5.3, rel=0.02)
    assert rows[1]['period'] == pytest.approx(7.7, rel=0.02)
    assert rows[0]['files'] == 1 and rows[1]['files'] == 0
    for star in ('GJ_1', 'GJ_2'):
        for part in ('result.json', 'quick.json', 'quicklook.pdf'):
            assert os.path.exists(os.path.join(moved, 'results', star, part))
    with open(os.path.join(moved, 'results', 'GJ_2', 'result.json')) as handle:
        res = json.load(handle)
    assert res['summary']['instruments'] == dict(HARPS15=30)
    assert res['datasets']['used'] == 1
    table = open(os.path.join(moved, 'results', 'table.csv')).read()
    assert table.splitlines()[0].startswith('name,sptype,distance,status')
    assert len(table.splitlines()) == 3
    # done: not done again (one star in two asked: nothing left to do)
    done = []
    monkeypatch.setattr(survey, 'run_target', lambda root, star, *args:
                        done.append(star['name']))
    survey.run(moved)
    survey.run(moved, part=(1, 2))
    assert done == []
    survey.run(moved, only=['gj 2'], again=True)
    assert done == ['GJ 2']


def test_the_survey_tab_of_the_page(tmp_path, monkeypatch):
    """the questions of the survey tab: a sample kept by the server, its
    archives checked in a thread, its files, the batch of those ticked
    (as a batch of the page), the results of a batch folder opened"""
    import time
    from koloa import apero_names, gather, gui, gui_survey, published
    _simbad(monkeypatch)
    state = gui_survey.route('/api/survey/sample', dict(
        sp_from='M0', sp_to='M9', dmax='15', dec_min='', dec_max='',
        vmax='', dwarfs=True))
    sid = state['id']
    assert state['n'] == 3 and state['stars'][-1]['name'] == 'GJ 876'
    assert state['asked']['dmax'] == 15.0 and state['check'] is None
    monkeypatch.setattr(gather, 'carmenes_objects', lambda refresh=False: [])
    monkeypatch.setattr(published, 'survey_stars',
                        lambda survey, refresh=False, timeout=120.0: [])
    monkeypatch.setattr(gather, 'dace_rv', lambda *args, **kwargs: None)
    gui_survey.route('/api/survey/check', dict(id=sid, names=['gj 876'],
                                               root=str(tmp_path / 'arch')))
    for _ in range(100):
        now = gui_survey.route('/api/survey', dict(id=sid, stars='0'))
        if now['check']['status'] != 'running':
            break
        time.sleep(0.05)
    assert now['check']['status'] == 'done' and now['check']['total'] == 1
    assert 'stars' not in now
    full = gui_survey.route('/api/survey', dict(id=sid))
    assert full['stars'][-1]['archives']['dace'] is None
    assert full['stars'][0]['archives'] is None
    # the summary of the cross-match: one star asked, with nothing
    assert (full['overview']['checked'], full['overview']['data'],
            full['overview']['none']) == (1, 0, 1)
    assert full['stars'][-1]['summary'] == []
    # the files of a folder, by the names of the page
    (tmp_path / 'lbl_GL876_GL876.rdb').write_text('rjd\tvrad\tsvrad\n')
    monkeypatch.setattr(apero_names, 'star_of_file', lambda path, fetch=True:
                        dict(raw='GL876', target='BD-15 6290', apero='GL876',
                             entry=None))
    state = gui_survey.route('/api/survey/files', dict(
        id=sid, folders=f'{tmp_path}, ', pattern='*.rdb'))
    assert state['matched'] == 1
    assert state['stars'][-1]['files'] == ['lbl_GL876_GL876.rdb']
    # the batch of those ticked: the stars as a batch takes them
    asked = {}

    def batch_fip(paths, opts, archives, regather, root, rules, targets=None):
        asked.update(targets=targets, root=root, rules=rules)
        return dict(id='b1')
    monkeypatch.setattr(gui, 'batch_fip', batch_fip)
    assert gui_survey.route('/api/survey/run', dict(
        id=sid, names=['GJ 876'], root='arch', rules=False)) == dict(id='b1')
    assert asked['targets'][0]['name'] == 'GJ 876' and not asked['rules']
    assert asked['targets'][0]['files'][0].endswith('lbl_GL876_GL876.rdb')
    with pytest.raises(ValueError, match='no star ticked'):
        gui_survey.route('/api/survey/run', dict(id=sid, names=[]))
    with pytest.raises(ValueError, match='a name'):
        gui_survey.route('/api/survey/pack', dict(id=sid, names=['GJ 876']))
    with pytest.raises(ValueError, match='no such question'):
        gui_survey.route('/api/survey/nothing', dict(id=sid))


def test_the_results_of_a_batch_folder_as_a_batch_of_the_page(tmp_path,
                                                              monkeypatch):
    from koloa import gui, gui_survey
    stars = _batch(tmp_path, monkeypatch, nstar=2)
    made = survey.pack(stars, 'back', out=str(tmp_path / 'out'),
                       root=str(tmp_path / 'arch'), gather=False, tar=False)
    with pytest.raises(ValueError, match='no star done'):
        gui_survey.results(dict(root=made['folder']))
    # one star done: its result and its FIP, as the batch keeps them
    out = os.path.dirname(survey.result_path(made['folder'], 'GJ 1'))
    os.makedirs(out)
    with open(os.path.join(out, 'result.json'), 'w') as handle:
        json.dump(dict(name='GJ 1', status='done', error=None, note=None,
                       star=dict(target='GJ 1'), summary=dict(n=40,
                                                              period=5.3),
                       datasets=dict(used=1, all=1, told=[]), elapsed=12.0),
                  handle)
    with open(os.path.join(out, 'quick.json'), 'w') as handle:
        json.dump(dict(status='done', result=dict(period=[1.0, 2.0]),
                       elapsed=12.0), handle)
    state = gui_survey.results(dict(root=made['folder']))
    assert state['status'] == 'done' and state['todo'] == 1
    assert [item['name'] for item in state['items']] == ['GJ 1']
    item = state['items'][0]
    assert item['summary']['period'] == 5.3 and item['datasets']['used'] == 1
    assert gui.QUICKS[item['qid']]['result'] == dict(period=[1.0, 2.0])
    assert gui.BATCHES[state['id']]['root'] == os.path.join(made['folder'],
                                                           'archives')


def test_a_terminal_and_its_routes(tmp_path, monkeypatch):
    """a shell answers who has the key only; what it showed as it was
    typed is proposed for a route, not what was typed with the echo off;
    a route kept is typed again"""
    import base64
    import time
    from koloa import terminal
    if not terminal.available():
        pytest.skip('no pseudo-terminal here')
    monkeypatch.setattr(terminal, 'ROUTES', str(tmp_path / 'routes.json'))
    assert terminal.route('/api/term/state', {}, None) == dict(
        available=True, allowed=False)
    assert terminal.route('/api/term/state', {}, terminal.KEY)['allowed']
    for key in (None, '', 'not-the-key'):
        with pytest.raises(PermissionError):
            terminal.route('/api/term/open', {}, key)
    sid = terminal.route('/api/term/open', dict(cols=120, rows=24,
                                                cwd=str(tmp_path)),
                         terminal.KEY)['id']

    def typed(text, pause=1.0):
        for char in text + '\r':
            terminal.route('/api/term/write', dict(id=sid, data=char),
                           terminal.KEY)
            time.sleep(0.01)
        time.sleep(pause)

    def shown(since=0):
        out = terminal.route('/api/term/read', dict(id=sid, since=since,
                                                    wait=0.1), terminal.KEY)
        return base64.b64decode(out['data']).decode('utf-8', 'replace'), out
    try:
        time.sleep(1.5)
        typed('echo koloa-$((6*7))')
        typed('python3 -c "import getpass; getpass.getpass()"')
        typed('not-to-be-kept')
        typed('echo after')
        text, out = shown()
        assert 'koloa-42' in text and out['alive']
        assert 'not-to-be-kept' not in text
        lines = terminal.route('/api/term/typed', dict(id=sid),
                               terminal.KEY)['lines']
        assert 'echo koloa-$((6*7))' in lines and 'echo after' in lines
        assert 'not-to-be-kept' not in lines
        kept = terminal.route('/api/term/route_save', dict(
            name='there', lines='echo route-$((2+3))\n\n', host='me@there',
            folder='/data/batches'), terminal.KEY)['routes']
        assert kept['there']['lines'] == ['echo route-$((2+3))']
        assert json.load(open(tmp_path / 'routes.json'))['there']['host'] \
            == 'me@there'
        terminal.route('/api/term/go', dict(id=sid, name='there'),
                       terminal.KEY)
        time.sleep(3.0)
        assert 'route-5' in shown(out['next'])[0]
        terminal.route('/api/term/resize', dict(id=sid, cols=90, rows=30),
                       terminal.KEY)
        assert terminal.route('/api/term/route_delete', dict(name='there'),
                              terminal.KEY)['routes'] == {}
    finally:
        terminal.route('/api/term/close', dict(id=sid), terminal.KEY)
    with pytest.raises(ValueError, match='no such terminal'):
        terminal.route('/api/term/read', dict(id=sid), terminal.KEY)

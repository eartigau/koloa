"""
koloa.survey, without the network: a spectral type as a number, a sample
from SIMBAD's answer, what the archives have of a star, the files put
with the stars, a batch folder packed and run.
"""
import ast
import glob
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
    from koloa import archive
    monkeypatch.setattr(archive, 'resolve', lambda name, **kw: dict(
        name=name, main=name, aliases=[]))
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
    (tmp_path / 'arch' / 'GJ_2' / 'manifest.json').write_text(json.dumps(
        dict(folder=str(tmp_path / 'arch' / 'GJ_2'), archives=dict(
            dace=dict(status='ok', npoints=30)))))
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
    for part in ('run_batch.py', 'submit.sh', 'README.txt', 'CLAUDE.md',
                 'targets.json',
                 'files/GJ_1/star0.csv', 'files/GJ_3/star2.csv',
                 'archives/GJ_2/rv/all_rv.csv', 'cache/archive/kept.json',
                 'cache/carmenes_objects.json', 'koloa_src/koloa/survey.py',
                 'koloa_src/koloa/gui_static/gui.js'):
        assert os.path.exists(os.path.join(folder, part)), part
    script = open(os.path.join(folder, 'run_batch.py')).read()
    ast.parse(script)
    assert "ROOT = '/scratch/me/m_dwarfs_15pc'" in script
    assert 'JOBS = 4' in script and "if __name__ == '__main__':" in script
    # the stars at once given where it runs, in place of the JOBS packed
    assert "if '--jobs' in args:" in script and 'jobs=jobs' in script
    assert 'N = (2/3 x free cores) / 3' in open(os.path.join(
        folder, 'README.txt')).read()
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
    # what the batch writes in is there, empty, once unpacked
    assert 'm_dwarfs_15pc/logs' in names and 'm_dwarfs_15pc/results' in names
    assert made['size'] == os.path.getsize(made['tar'])
    # its README: how to run it there, for a person or an agent, its
    #   stars, and no path of this machine
    readme = open(os.path.join(folder, 'README.txt')).read()
    assert readme.startswith('koloa: the batch m_dwarfs_15pc\n====')
    assert "ROOT = '/scratch/me/m_dwarfs_15pc'" in readme
    assert 'FOR A CLAUDE SESSION' in readme and '--check' in readme
    # for the Claude of that machine, in the file it reads on its own: the
    #   steps of the README written for it, among them the update of the
    #   koloa of the machine, which may be out of date
    claude = open(os.path.join(folder, 'CLAUDE.md')).read()
    assert claude.startswith('# koloa: the batch m_dwarfs_15pc\n')
    assert '  a. Work from this folder' in claude and '  h. Then read' in claude
    assert 'c. Update the koloa of this machine' in claude
    assert ('pip install --upgrade --force-reinstall --no-deps \\\n'
            '             git+https://github.com/eartigau/koloa.git') in claude
    assert 'N = (2/3 x free cores) / 3' in claude
    assert '1. COPY IT' not in claude and str(tmp_path) not in claude
    assert 'm_dwarfs_15pc/CLAUDE.md' in tarfile.open(made['tar']).getnames()
    assert 'koloa_commit' in survey.targets_of(folder)
    assert 'about\n       12 cores' in readme or '12 cores' in readme
    assert ('rsync -av me@server:/scratch/me/m_dwarfs_15pc/results/ \\\n'
            '          m_dwarfs_15pc/results/') in readme
    lines = readme[readme.index('9. THE STARS'):].splitlines()
    assert lines[4].split()[:4] == ['GJ', '1', 'M3V', '5.00']
    assert lines[3].endswith('P rot [d]  its archives')
    assert lines[4].split()[4:] == ['1', 'none']
    assert lines[5].split()[4:] == ['0', 'DACE', '30']
    assert str(tmp_path) not in readme and '{' not in readme
    # (the lines of its stars are as long as their archives are many)
    assert max(len(line) for line in readme.splitlines()
               if not line.startswith('    GJ ')) <= 78
    # nor in the manifest that goes with the archives of a star
    with open(os.path.join(folder, 'archives', 'GJ_2',
                           'manifest.json')) as handle:
        assert json.load(handle)['folder'] == os.path.join('archives', 'GJ_2')
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
    asked = []
    monkeypatch.setattr(gui, 'tess_light',
                        lambda target, root='', fetch=False: asked.append(
                            fetch) or (None, 'none'))
    monkeypatch.setattr(gui, 'space_light', lambda target, root='',
                        fetch=False, mission='tess': asked.append(fetch) or (
                            None, 'none', {}))
    # before: everything is there, nothing computed
    told = survey.status(moved, jobs=2)
    assert told['state'] == 'ready' and told['todo'] == 2
    assert told['cores']['jobs'] >= 1 and told['cores']['each'] == 3
    assert told['missing'] == [] and not os.path.exists(
        os.path.join(moved, 'results', 'table.csv'))
    rows = survey.run(moved)
    told = survey.status(moved)
    assert told['state'] == 'done' and told['done'] == 2
    assert told['stars'][0]['told'].startswith('P = 5.3')
    assert glob.glob(os.path.join(moved, 'logs', '*.pid')) == []
    assert glob.glob(os.path.join(moved, 'results', '*', 'running.json')) \
        == []
    assert [row['name'] for row in rows] == ['GJ 1', 'GJ 2']
    assert [row['status'] for row in rows] == ['done', 'done']
    assert rows[0]['period'] == pytest.approx(5.3, rel=0.02)
    assert rows[1]['period'] == pytest.approx(7.7, rel=0.02)
    assert rows[0]['files'] == 1 and rows[1]['files'] == 0
    # nothing is asked of the network: a light curve that did not come
    #   with the archives of a star is not fetched, and the table says so
    assert asked and not any(asked)
    assert rows[0]['transit'] == 'no light curve'
    # how each star reads (no known planet here: a candidate), no report
    #   unless asked, the summary of them all, and its page first in the
    #   PDF of each star
    assert [row['verdict'] for row in rows] == ['candidate', 'candidate']
    assert [row['report'] for row in rows] == ['none', 'none']
    assert os.path.getsize(os.path.join(moved, 'results',
                                        'summary.pdf')) > 20000
    for star in ('GJ_1', 'GJ_2'):
        # its PDF under the name of the star: those of every star can be
        #   put in one folder
        for part in ('result.json', 'quick.json', f'{star}_quicklook.pdf'):
            assert os.path.exists(os.path.join(moved, 'results', star, part))
    assert survey.quicklook_path(moved, 'GJ 2') == os.path.join(
        moved, 'results', 'GJ_2', 'GJ_2_quicklook.pdf')
    assert len({os.path.basename(path) for path in glob.glob(os.path.join(
        moved, 'results', '*', '*_quicklook.pdf'))}) == 2
    with open(os.path.join(moved, 'results', 'GJ_2', 'result.json')) as handle:
        res = json.load(handle)
    assert res['summary']['instruments'] == dict(HARPS15=30)
    assert res['datasets']['used'] == 1
    assert res['reading']['kind'] == 'candidate'
    assert res['reading']['peaks'][0]['counts']
    # a candidate has its prospects: its astrometric signal (its mass from
    #   its type, its parallax from its distance), in its result, in the
    #   table, and on a page of its own in the summary
    assert res['prospects']['astrometry']['alpha'] > 0
    assert res['prospects']['astrometry']['distance'] == pytest.approx(6.0)
    assert rows[1]['gaia_alpha'] == res['prospects']['astrometry']['alpha']
    assert rows[1]['gaia_snr'] < 1 and rows[1]['snr'] > 5
    assert rows[1]['gaia_chi2'] == res['prospects']['astrometry']['dchi2']
    assert res['reading']['line'].startswith('candidate: #1 at 7.')
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

    def batch_fip(paths, opts, archives, regather, root, rules, targets=None,
                  report=None):
        asked.update(targets=targets, root=root, rules=rules, report=report)
        return dict(id='b1')
    monkeypatch.setattr(gui, 'batch_fip', batch_fip)
    assert gui_survey.route('/api/survey/run', dict(
        id=sid, names=['GJ 876'], root='arch', rules=False)) == dict(id='b1')
    assert asked['targets'][0]['name'] == 'GJ 876' and not asked['rules']
    # its published rotation goes with the star; no report unless asked
    assert asked['targets'][0]['rotation'][0]['period'] == 82.8
    assert asked['report'] is None
    gui_survey.route('/api/survey/run', dict(
        id=sid, names=['GJ 876'], report=True, report_fip='0.001'))
    assert asked['report'] == dict(fip=0.001, folder='reports')
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


def test_where_a_batch_is_and_its_stop(tmp_path, monkeypatch):
    """the state of a batch folder without computing (what lacks, the
    stars being computed and how far), a batch started twice, and the
    stop of its processes and of no other"""
    import subprocess
    import sys
    stars = _batch(tmp_path, monkeypatch, nstar=2)
    made = survey.pack(stars, 'two', out=str(tmp_path / 'out'),
                       root=str(tmp_path / 'arch'), gather=False, tar=False)
    root = made['folder']
    # a file that did not come with the folder: not ready, and why
    kept = os.path.join(root, 'files', 'GJ_1', 'star0.csv')
    os.rename(kept, kept + '.gone')
    told = survey.status(root)
    assert told['state'] == 'not ready'
    assert told['missing'] == ['GJ 1: no ' + os.path.join(
        'files', 'GJ_1', 'star0.csv')]
    os.rename(kept + '.gone', kept)
    # a part that computes a star, as another process of this machine
    part = subprocess.Popen([sys.executable, '-c', 'import time; '
                             'time.sleep(600)', 'koloa.survey'],
                            start_new_session=True)
    other = subprocess.Popen([sys.executable, '-c', 'import time; '
                              'time.sleep(600)'], start_new_session=True)
    try:
        import socket
        os.makedirs(os.path.join(root, 'logs'), exist_ok=True)
        for name, proc in (('part_0_of_2', part), ('stale', other)):
            with open(os.path.join(root, 'logs', name + '.pid'), 'w') as out:
                json.dump(dict(pid=proc.pid, host=socket.gethostname()), out)
        mark = os.path.join(root, 'results', 'GJ_1', 'running.json')
        os.makedirs(os.path.dirname(mark))
        with open(mark, 'w') as out:
            json.dump(dict(pid=part.pid, host=socket.gethostname()), out)
        with open(os.path.join(root, 'logs', 'part_0.log'), 'w') as out:
            out.write('261007 10:00:00.00 | batch: GJ 1 (1 of 1)\n'
                      '261007 10:05:00.00 | quick FIP, first pass: 40 % of '
                      'the sweeps, 5 min so far, about 7 min left\n')
        told = survey.status(root)
        assert told['state'] == 'running' and told['running'] == 1
        assert told['stars'][0] == dict(
            name='GJ 1', state='running', told='quick FIP, first pass: '
            '40 % of the sweeps, 5 min so far, about 7 min left')
        assert told['stars'][1]['state'] == 'to do'
        # started a second time: the star being computed is left to it
        done = []
        monkeypatch.setattr(survey, 'run_target', lambda root, star, *args:
                            done.append(star['name']) or dict(
                                status='done', elapsed=0.0, summary={}))
        survey.run(root)
        assert done == ['GJ 2']
        # stopped: the part and what it leads, not the process that only
        #   has a number in the folder (another's by now, say)
        assert survey.stop(root) == 1
        assert part.wait(timeout=20) != 0 and other.poll() is None
        assert glob.glob(os.path.join(root, 'logs', '*.pid')) == []
        assert not os.path.exists(mark)
        told = survey.status(root)
        assert told['state'] == 'ready' and told['todo'] == 2
    finally:
        for proc in (part, other):
            if proc.poll() is None:
                proc.kill()


def test_how_a_star_of_a_batch_reads():
    """the peaks of a quick look: a known planet, the rotation, a drift, a
    candidate; and the verdict of the star"""
    from koloa import batchpdf

    def told(peaks, known=(), spins=(), baseline=1000.0, status='done'):
        res = dict(status=status, error='boom', summary=dict(
            baseline=baseline))
        quick = dict(result=dict(
            peak_list=[dict(id=it + 1, period=per, family=fip, alone=fip,
                            named=True) for it, (per, fip) in
                       enumerate(peaks)],
            folds=[dict(id=it + 1, K=5.0, K_err=0.5) for it in
                   range(len(peaks))],
            known=[dict(name=f'Star {name}', P=per) for name, per in known]))
        return batchpdf.reading(res, quick, dict(
            sptype='M3V', rotation=[dict(period=per, source='a paper')
                                    for per in spins]))
    # a peak that is nothing known: a candidate, with its minimum mass
    out = told([(12.3, 1e-9), (40.0, 0.4)])
    assert out['kind'] == 'candidate' and out['peaks'][0]['counts']
    assert not out['peaks'][1]['counts']
    assert out['line'].startswith('candidate: #1 at 12.3000 d, FIP 1.0e-09')
    assert 5 < out['peaks'][0]['msini'][0] < 30
    # the same peak when a planet is known there (within 1 %)
    out = told([(12.3, 1e-9)], known=[('b', 12.25)])
    assert out['kind'] == 'known' and 'Star b' in out['line']
    assert told([(12.3, 1e-9)], known=[('b', 13.0)])['kind'] == 'candidate'
    # a known planet and another peak: the other is the candidate
    out = told([(12.3, 1e-30), (33.0, 1e-5)], known=[('b', 12.3)])
    assert out['kind'] == 'candidate' and '#2 at 33.0000 d' in out['line']
    # as long as the series: a drift
    assert told([(900.0, 1e-9)])['kind'] == 'drift'
    # at the rotation, or its half: still a candidate, and said
    out = told([(50.5, 1e-9)], spins=[101.0])
    assert out['kind'] == 'candidate'
    assert out['peaks'][0]['notes'] == ['at the rotation / 2 (101 d)']
    assert 'at the rotation / 2' in out['line']
    # nothing below the limit; a star that failed
    out = told([(12.3, 0.2)])
    assert out['kind'] == 'nothing' and '12.3000 d' in out['line']
    assert told([(12.3, 1e-9)], status='failed')['kind'] == 'failed'
    assert batchpdf.reading(dict(status='done'), None)['kind'] == 'nothing'


def test_the_detailed_report_of_a_candidate(tmp_path, monkeypatch):
    """a star with a candidate has its detailed report on the series its
    quick look used: an SHO GP at its published rotation, by band without
    one; kept in its result, its table and its summary; a report that was
    stopped is taken again without the quick look"""
    import subprocess
    from koloa import batchpdf, gui
    stars = _batch(tmp_path, monkeypatch, nstar=2)
    stars[0]['rotation'] = [dict(period=41.5, source='2023A&A...672A..52F'),
                            dict(period=40.0, source='CARMENES DR1')]
    made = survey.pack(stars, 'two', out=str(tmp_path / 'out'),
                       root=str(tmp_path / 'arch'), gather=False, tar=False,
                       report=True, report_fip=0.001)
    root = made['folder']
    script = open(os.path.join(root, 'run_batch.py')).read()
    assert 'REPORT = True' in script and 'REPORT_FIP = 0.001' in script
    assert '--summary' in script
    held = survey.targets_of(root)
    assert held['targets'][0]['rotation'][0]['period'] == 41.5
    assert held['options']['report'] is True
    monkeypatch.setattr(gui, 'known_periods', lambda target: [])
    monkeypatch.setattr(gui, 'transits_of', lambda target, known: [])
    monkeypatch.setattr(gui, '_star_planets', lambda target: ([], []))
    monkeypatch.setattr(gui, 'QUICK', dict(kmax=1, nsweep=150, nburn=80))
    monkeypatch.setattr(gui, 'space_light', lambda target, root='',
                        fetch=False, mission='tess': (None, 'none', {}))
    # the report itself: its command line, and what it leaves
    ran = []

    def fake(cmd, env=None, stdout=None, stderr=None):
        ran.append(cmd)
        out = cmd[cmd.index('--outdir') + 1]
        stem = cmd[cmd.index('--name') + 1].replace(' ', '_')
        with open(os.path.join(out, f'{stem}_summary.json'), 'w') as handle:
            json.dump(dict(known=dict(star={}), orbits=[dict(
                P=[5.3, 0.1, 0.1], K=[9.0, 0.3, 0.3], msini=[11.0, 1, 1])],
                duck={'5.3000': 'PLANET CANDIDATE: it quacks'}), handle)
        open(os.path.join(out, f'{stem}_report.txt'), 'w').write('a report')
        return subprocess.CompletedProcess(cmd, 0)
    monkeypatch.setattr(subprocess, 'run', fake)
    rows = survey.run(root, report=True, report_fip=0.001)
    assert [row['report'] for row in rows] == ['done', 'done']
    first, second = ran
    assert first[1:3] == ['-m', 'koloa.cli'] and '--detailed' in first
    assert first[3].endswith(os.path.join('GJ_1', 'report', 'velocities.csv'))
    # the first rotation published: an SHO at it; none: the archive's, else
    #   by band; nothing asked of the network; the rules already applied
    assert first[first.index('--rotation') + 1] == '41.5'
    assert '--rotation' not in second
    assert second[second.index('--fip-gp') + 1] == 'sho'
    assert '--no-tess' in first and '--no-rules' in first
    assert first[first.index('--target') + 1] == 'GJ 1'
    with open(survey.result_path(root, 'GJ 1')) as handle:
        res = json.load(handle)
    rep = res['report']
    assert rep['gp'] == 'SHO at the rotation, 41.5 d (2023A&A...672A..52F)'
    assert rep['signals'] == ['5.3000 d, K = 9.00 m/s, m sin i = 11.0 ME: '
                              'PLANET CANDIDATE: it quacks']
    assert rep['text'] == os.path.join('results', 'GJ_1', 'report',
                                       'GJ_1_report.txt')
    assert str(tmp_path) not in json.dumps(rep)
    with open(survey.result_path(root, 'GJ 2')) as handle:
        assert json.load(handle)['report']['gp'] == (
            'local, by period band (no rotation published)')
    # the series the report was made on: the one of the quick look
    kept = open(os.path.join(root, 'results', 'GJ_2', 'report',
                             'velocities.csv')).read().splitlines()
    assert kept[0].startswith('rjd,vrad,svrad,inst') and len(kept) == 31
    told = survey.status(root)
    assert told['state'] == 'done'
    assert told['stars'][0]['told'].endswith('candidate; report done')
    # the page of a star, from what is kept
    fig = batchpdf.figure(res, None, held['targets'][0])
    assert any('detailed report: done' in text.get_text()
               for text in fig.texts)
    # a report stopped half way: taken again, the quick look kept
    res['report'] = dict(status='running')
    survey._keep(root, res)
    told = survey.status(root)
    assert told['state'] == 'ready' and told['todo'] == 1
    assert 'its detailed report was stopped' in told['stars'][0]['told']
    done = []
    monkeypatch.setattr(survey, 'run_target', lambda *args, **kw:
                        done.append(args))
    del ran[:]
    survey.run(root, report=True)
    assert done == [] and len(ran) == 1 and '--rotation' in ran[0]
    with open(survey.result_path(root, 'GJ 1')) as handle:
        assert json.load(handle)['report']['status'] == 'done'
    # a report that fails: said, and the star is kept
    monkeypatch.setattr(subprocess, 'run', lambda cmd, **kw:
                        subprocess.CompletedProcess(cmd, 3))
    out = survey.report_star(str(tmp_path / 'again'), None, 'GJ 1')
    assert out['status'] == 'failed' and 'code 3' in out['error']


def test_the_stars_at_once_on_two_thirds_of_the_free_cores(monkeypatch):
    """the cores of the machine that runs a batch, those that are free,
    and the stars at once that two thirds of them allow"""
    monkeypatch.setattr(os, 'cpu_count', lambda: 96)
    monkeypatch.setattr(os, 'sched_getaffinity', lambda pid: set(range(96)),
                        raising=False)
    monkeypatch.setattr(os, 'getloadavg', lambda: (6.2, 5.0, 4.0))
    have = survey.cores()
    assert (have['total'], have['mine'], have['free']) == (96, 96, 90)
    # two thirds of 90 free cores are 60: 20 stars of three cores each,
    #   15 of four with their detailed reports
    assert have['share'] == 60 and have['jobs'] == 20
    assert survey.cores(report=True)['jobs'] == 15
    # fewer cores for this process than the machine has (a scheduler)
    monkeypatch.setattr(os, 'sched_getaffinity', lambda pid: set(range(8)),
                        raising=False)
    monkeypatch.setattr(os, 'getloadavg', lambda: (0.2, 0.2, 0.2))
    have = survey.cores()
    assert have['mine'] == 8 and have['free'] == 8 and have['jobs'] == 1
    # a machine busier than it has cores, or that does not say: one star

    def unknown():
        raise OSError('no load average here')
    monkeypatch.setattr(os, 'getloadavg', lambda: (40.0, 40.0, 40.0))
    assert survey.cores()['jobs'] == 1 and survey.cores()['free'] == 1
    monkeypatch.setattr(os, 'getloadavg', unknown)
    assert survey.cores()['load'] is None and survey.cores()['free'] == 8

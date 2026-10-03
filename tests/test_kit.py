#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa's analysis kit: the script it writes (both languages), the YAML of the
star, and a kit that runs on its own, offline

Created on 2026-10-03

@author: artigau
"""
import io
import json
import os
import subprocess
import sys
import tarfile

import numpy as np
import pytest

from koloa import kit

SETTINGS = dict(TITLE='a test', MADE='made here', STAR="'GJ 1'",
                FILES="[('data/files/0_a.csv', None)]",
                ARCHIVE_FOLDER="'GJ_1'", USE_DACE=True, USE_CARMENES=False,
                USE_VIZIER=False, EXCLUDE='[]', TREND=1, KMAX=2, NSWEEP=200,
                NBURN=50, PMIN=1.1, PMAX='None', FIP_GP="'none'",
                ROTATION='None', PERIODS='[]', NIGHTLY=True, MSTAR='0.5',
                MSTAR_ERR='0.05')


@pytest.mark.parametrize('lang', ['en', 'fr'])
def test_the_script_in_both_languages(lang):
    """every step explained in the language asked, every setting filled,
    Python that compiles, and no long dash"""
    text = kit.render(SETTINGS, lang)
    compile(text, 'analysis.py', 'exec')
    assert '{{' not in text and '#>' not in text
    assert '\u2014' not in text and '\u2013' not in text
    word = 'Les réglages' if lang == 'fr' else 'The settings'
    assert word in text
    for key, texts in kit.COMMENTS.items():
        assert set(texts) == {'en', 'fr'}, key
    # every comment is used by the script
    used = {line.strip()[3:] for line in kit.SCRIPT.split('\n')
            if line.strip().startswith('#> ')}
    assert used == set(kit.COMMENTS)
    # a setting missing: said, not left in the script
    with pytest.raises(ValueError, match='KMAX'):
        kit.render({key: val for key, val in SETTINGS.items()
                    if key != 'KMAX'}, lang)


def test_the_yaml_of_the_star():
    """star.yaml, written without PyYAML, read by it as the dict it was"""
    yaml = pytest.importorskip('yaml')
    star = dict(name='AN Sex', main='BD-03 2870', hd=None, ra=153.07,
                plx=129.75, aliases=['V* AN Sex', '[RHG95] 1595', 'GJ 382'],
                variability=[dict(type='ROT', period=21.56,
                                  bibcode='2012AcA....62...67K')],
                planets=[dict(name='b', P=3.5, K=None, reference='a: b "c"')],
                tois=[], carmenes=dict(karmn='J10122-037', nobs='77'),
                flag=True, empty={})
    assert yaml.safe_load(kit._yaml(star)) == star


def _series(tmp_path):
    """a file of SPIRou (a planet of 7.3 days, K = 9 m/s) and the archives
    of the star (HARPS15 on DACE, the same planet)"""
    from koloa.data import RVData
    from koloa.gather import write_rv
    rng = np.random.default_rng(3)
    nights = np.sort(rng.choice(np.arange(60000, 60400), 70, replace=False))
    tfile = nights + rng.uniform(0.2, 0.4, len(nights))
    signal = lambda tt: 9.0 * np.sin(2 * np.pi * tt / 7.3)
    path = tmp_path / 'spirou.csv'
    path.write_text('rjd,vrad,svrad\n' + ''.join(
        f'{tt!r},{signal(tt) + rng.normal(0, 2)!r},2.0\n' for tt in tfile))
    tarch = np.sort(rng.uniform(59000, 59400, 30))
    star = tmp_path / 'arch' / 'GJ_1'
    write_rv(RVData(tarch, signal(tarch) + rng.normal(0, 2.5, 30),
                    np.full(30, 2.5), inst=np.array(['HARPS15'] * 30)),
             str(star / 'rv' / 'all_rv.csv'))
    (star / 'manifest.json').write_text(json.dumps(dict(
        target='GJ 1', archives=dict(dace=dict(
            status='ok', instruments=dict(HARPS15=30), npoints=30,
            key='the key of DACE_API_KEY: may hold data that are not '
                'public')))))
    return path


def test_a_kit_that_runs_on_its_own(tmp_path, monkeypatch):
    """the kit of a page: the file, the archives, the star (YAML and JSON),
    the script with the page's settings; extracted elsewhere, it runs
    offline, finds the planet, and saves its figures as PDF"""
    path = _series(tmp_path)
    monkeypatch.setattr(kit, 'star_info', lambda target, root='': dict(
        name=target, main=target, sptype='M2V', mass=0.44, mass_err=0.044,
        planets=[], tois=[], variability=[]))
    opts = dict(target='GJ 1', root=str(tmp_path / 'arch'),
                files=[dict(path=str(path), label='SPIRou')], dace=True,
                kmax='1', nsweep='250', nburn='80', fip_gp='none',
                trend=True, exclude='')
    data = kit.build(opts, 'fr')
    out = tmp_path / 'elsewhere'
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
        names = tar.getnames()
        tar.extractall(out)
    top = 'koloa_GJ_1_analysis'
    for name in ('analysis.py', 'README.md', 'star.yaml', 'star.json',
                 'data/files/0_spirou.csv',
                 'data/archives/GJ_1/rv/all_rv.csv',
                 'data/archives/GJ_1/manifest.json'):
        assert f'{top}/{name}' in names, name
    folder = out / top
    script = (folder / 'analysis.py').read_text()
    assert "('data/files/0_spirou.csv', 'SPIRou')" in script
    assert 'USE_DACE = True' in script and 'NSWEEP = 250' in script
    assert 'MSTAR = 0.44' in script
    readme = (folder / 'README.md').read_text()
    assert 'HARPS15 30' in readme and 'non publiques' in readme
    # run where it is, koloa from this tree
    env = dict(os.environ, PYTHONPATH=os.path.dirname(os.path.dirname(
        os.path.abspath(kit.__file__))), MPLBACKEND='Agg')
    run = subprocess.run([sys.executable, 'analysis.py'], cwd=folder,
                         env=env, capture_output=True, text=True,
                         timeout=600)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    for fig in ('series.pdf', 'fip.pdf', 'folds.pdf', 'fip_residuals.pdf'):
        assert (folder / 'figures' / fig).exists(), fig
    res = json.loads((folder / 'results.json').read_text())
    assert res['instruments'] == ['SPIRou', 'HARPS15'] or set(
        res['instruments']) == {'SPIRou', 'HARPS15'}
    periods = [sig['P'] for sig in res['signals']]
    assert any(abs(per / 7.3 - 1) < 0.01 for per in periods), periods
    sig = [sig for sig in res['signals'] if abs(sig['P'] / 7.3 - 1) < 0.01][0]
    assert abs(sig['K'] - 9.0) < 2.0 and sig['masses']['earth'][0] > 0


def test_a_kit_needs_velocities(tmp_path):
    """no file and nothing gathered: no kit, and why"""
    with pytest.raises(ValueError, match='nothing to put in the kit'):
        kit.build(dict(target='GJ 1', root=str(tmp_path / 'none')))

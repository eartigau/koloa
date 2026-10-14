#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LBL velocities of a star from the public spectra of the archives: the
spectra fetched, and lbl-rapide (the branch of LBL made faster; Artigau
et al. 2022) run on them.

koloa never imports LBL: LBL needs versions of its own (Python 3.12,
numpy 2.3, numba, astropy 7), in an environment of its own (lbl-rapide by
default). koloa writes the parameters of the run in a short script, as one
writes a wrap file, and runs it with that Python, in a process of its own.
The rdb LBL makes is copied into the folder of the star, under its name
(archives/GJ_436/lbl/GJ_436_HARPSN_lbl.rdb), to be given as a file.

    koloa "GJ 436" --lbl HARPSN --archives archives --lbl-dir lbl_data

Each instrument has a folder of LBL of its own in lbl-dir (its templates,
masks, calibrations), with the models of LBL shared between them
(lbl-dir/models, a few hundred MB, fetched by LBL the first time).

The instruments:

HARPSN: the espdr S2D_A files of the archive of the TNG (blaze corrected,
on a barycentric wavelength grid: LBL's HARPSN ESO source, no calibration
file), 10 MB each, of the exposures koloa keeps of the star from their
CCFs (koloa.tng.series: those on the star, with its mask, the companion of
a binary left out by its velocity); the telluric lines cleaned by LBL
first.

Created on 2026-10-14

@author: artigau
"""
import gzip
import json
import os
import re
import shutil
import subprocess
import time as _time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Sequence

from koloa.log import log, outcome, step
from koloa.paths import config

# =============================================================================
# Define variables
# =============================================================================
#: your choices: the Python of LBL, the folder of its data, its processes
SETTINGS = config('lbl.json')
#: where an environment called lbl-rapide may be, when none is chosen
PYTHONS = [os.path.join(os.path.expanduser(base), 'envs', 'lbl-rapide',
                        'bin', 'python')
           for base in ('/opt/anaconda3', '/opt/miniconda3', '~/anaconda3',
                        '~/miniconda3', '~/miniforge3', '~/mambaforge',
                        '/opt/homebrew/Caskroom/miniforge/base')]
#: the files of the archive of the TNG
TNG_FILES = 'http://archives.ia2.inaf.it/files/tng/'
#: the instruments: their name in LBL, its data source, the files it reads,
#: whether the telluric lines are cleaned first, where the spectra are from
INSTRUMENTS: Dict[str, Dict[str, Any]] = {
    'HARPSN': dict(label='HARPS-N', lbl='HARPSN', source='ESO',
                   input='*S2D_A.fits', tellu=True, size=9.8,
                   dataset='HARPS-N (LBL)',
                   archive='the archive of the TNG (espdr S2D_A)'),
}
#: the temperatures of LBL's gradients of temperature [K]
GRADIENTS = (3000, 3500, 4000, 4500, 5000, 5500, 6000)
#: the files fetched at once
WORKERS = 4


# =============================================================================
# Define functions
# =============================================================================
def settings() -> Dict[str, Any]:
    """your choices (SETTINGS): python, data (the folder of LBL's data),
    nproc (its processes)"""
    out: Dict[str, Any] = dict(python='', data='', nproc=0)
    if os.path.exists(SETTINGS):
        try:
            with open(SETTINGS) as handle:
                out.update(json.load(handle))
        except (OSError, ValueError):
            pass
    return out


def save_settings(**kwargs: Any) -> Dict[str, Any]:
    """your choices, kept (those given replace those kept)"""
    out = settings()
    out.update({key: val for key, val in kwargs.items()
                if key in ('python', 'data', 'nproc')})
    os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
    with open(SETTINGS, 'w') as handle:
        json.dump(out, handle, indent=1)
    return out


def find_python(given: Optional[str] = None) -> str:
    """the Python of LBL: the one given, else the one chosen, else that of
    an environment called lbl-rapide where conda keeps them ('' when none)"""
    for one in [given, settings().get('python')] + PYTHONS:
        if one and os.path.isfile(os.path.expanduser(one)):
            return os.path.expanduser(one)
    return ''


def check_python(python: str, timeout: float = 120.0) -> Dict[str, Any]:
    """
    Whether a Python runs LBL: its version of LBL, and where it is

    :return: dict, ok, version, where, error
    """
    if not python or not os.path.isfile(python):
        return dict(ok=False, version='', where='',
                    error='no Python of LBL (an environment with lbl-rapide)')
    code = ('import os, lbl; from lbl.recipes import lbl_wrap; '
            'here = os.path.dirname(lbl.__file__); '
            'print(open(os.path.join(here, "version.txt")).read().split()[0]);'
            ' print(here)')
    try:
        done = subprocess.run([python, '-c', code], capture_output=True,
                              text=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError) as err:
        return dict(ok=False, version='', where='', error=str(err))
    lines = done.stdout.split()
    if done.returncode != 0 or len(lines) < 2:
        why = (done.stderr.strip().splitlines() or ['?'])[-1]
        return dict(ok=False, version='', where='', error=why)
    return dict(ok=True, version=lines[0], where=lines[1], error='')


def object_name(target: str) -> str:
    """the name of the star for LBL (letters and digits: GJ436)"""
    return re.sub(r'[^A-Za-z0-9]+', '', str(target)) or 'star'


def data_dir(lbl_dir: str, inst: str) -> str:
    """the folder of LBL of an instrument, its models shared"""
    about = INSTRUMENTS[inst]
    here = os.path.join(os.path.abspath(os.path.expanduser(lbl_dir)),
                        f'{about["lbl"]}_{about["source"]}')
    models = os.path.join(os.path.dirname(here), 'models')
    os.makedirs(models, exist_ok=True)
    os.makedirs(here, exist_ok=True)
    link = os.path.join(here, 'models')
    if not os.path.lexists(link):
        os.symlink(models, link)
    return here


def _ident(target: str) -> Dict[str, Any]:
    """the star as SIMBAD knows it, with how it moves"""
    from koloa.archive import resolve
    from koloa.secular import star_astrometry
    ident = resolve(target)
    try:
        ident['astrometry'] = star_astrometry(ident)
    except Exception as err:  # its place then is that of 2000
        log(f'astrometry: {type(err).__name__}: {err}', 'warn')
        ident['astrometry'] = None
    return ident


def harpsn_spectra(target: str, root: str = 'archives',
                   ident: Optional[Dict[str, Any]] = None) -> List[str]:
    """
    The espdr S2D_A files of HARPS-N of a star: those of the exposures
    koloa keeps from their CCFs (koloa.tng.series), read from the folder of
    the star when its archives were gathered, else from what koloa's users
    shared or from the archive

    :return: list of str, the names of the files (gzipped, as the archive
             has them), by date
    """
    from koloa import tng
    from koloa.gather import folder_name
    folder = os.path.join(root, folder_name(target), 'rv', 'published')
    table = os.path.join(folder, 'tng_exposures.csv')
    if os.path.exists(table):
        import csv
        with open(table) as handle:
            rows = list(csv.DictReader(handle))
    else:
        ident = ident or _ident(target)
        astro = ident.get('astrometry') or {}
        star = dict(ra=ident['ra'], dec=ident['dec'], pmra=astro.get('pmra'),
                    pmdec=astro.get('pmdec'),
                    names=[ident.get('main') or ''] + list(
                        ident.get('aliases') or []))
        os.makedirs(folder, exist_ok=True)
        rows = tng.exposures(star, folder)
    kept = [about for _, about in tng.series(rows, target)
            if about['reduction'] == 'espdr']
    if not kept:
        return []
    return [name.replace('_CCF_A.fits', '_S2D_A.fits')
            for name in kept[0]['files']]


def _download(url: str, path: str, timeout: float = 300.0) -> str:
    """a file of an archive, gunzipped where it is gzipped, written whole
    or not at all; '' or why not"""
    part = path + '.part'
    try:
        with urllib.request.urlopen(urllib.request.Request(
                url, headers={'User-Agent': 'koloa'}),
                timeout=timeout) as response, open(part, 'wb') as handle:
            source = (gzip.GzipFile(fileobj=response)
                      if url.endswith('.gz') else response)
            shutil.copyfileobj(source, handle, 1 << 20)
        os.replace(part, path)
        return ''
    except (urllib.error.URLError, TimeoutError, OSError, EOFError) as err:
        if os.path.exists(part):
            os.remove(part)
        return f'{type(err).__name__}: {err}'


def fetch(target: str, inst: str, lbl_dir: str, root: str = 'archives',
          names: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    """
    The spectra of a star for LBL, in the science folder of its instrument
    (lbl-dir/<INSTRUMENT>_<SOURCE>/science/<star>/): only those not there
    yet are fetched, a few at a time

    :param names: the files to fetch (as the archive has them), None for
                  those of the instrument (harpsn_spectra)

    :return: dict, folder, n (the spectra there), fetched, failed (and why,
             the first)
    """
    if inst != 'HARPSN':
        raise ValueError(f'no spectra of {inst} yet: {", ".join(INSTRUMENTS)}')
    names = list(harpsn_spectra(target, root) if names is None else names)
    folder = os.path.join(data_dir(lbl_dir, inst), 'science',
                          object_name(target))
    os.makedirs(folder, exist_ok=True)
    todo = [name for name in names if not os.path.exists(
        os.path.join(folder, name[:-3] if name.endswith('.gz') else name))]
    log(f'{INSTRUMENTS[inst]["label"]}: {len(names)} spectra of {target} in '
        f'{INSTRUMENTS[inst]["archive"]}, {len(todo)} to fetch (about '
        f'{len(todo) * INSTRUMENTS[inst]["size"]:.0f} MB)', 'info')

    def one(name):
        path = os.path.join(folder, name[:-3] if name.endswith('.gz')
                            else name)
        return _download(TNG_FILES + name, path)
    errors: List[str] = []
    if todo:
        from koloa.tng import _steps
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            for why in _steps(pool.map(one, todo), len(todo),
                              f'{INSTRUMENTS[inst]["label"]} spectra',
                              f'lbl_{inst}'):
                if why:
                    errors.append(why)
    have = len([name for name in os.listdir(folder)
                if not name.endswith('.part')])
    if errors:
        log(f'{len(errors)} spectra could not be fetched ({errors[0]})',
            'warn')
    return dict(folder=folder, n=have, fetched=len(todo) - len(errors),
                failed=len(errors), why=errors[0] if errors else '')


def parameters(target: str, inst: str, lbl_dir: str, teff: float,
               nproc: int = 4, tellu: Optional[bool] = None
               ) -> Dict[str, Any]:
    """the parameters of LBL for a star and an instrument (as in a wrap
    file): its template and masks made from its own spectra, the telluric
    lines cleaned first (optical spectra), the gradient of temperature of
    LBL nearest its temperature projected out (DTEMP)"""
    about = INSTRUMENTS[inst]
    obj = object_name(target)
    near = min(GRADIENTS, key=lambda one: abs(one - float(teff)))
    tellu = about['tellu'] if tellu is None else bool(tellu)
    return dict(
        INSTRUMENT=about['lbl'], DATA_SOURCE=about['source'],
        DATA_DIR=data_dir(lbl_dir, inst), INPUT_FILE=about['input'],
        DATA_TYPES=['SCIENCE'], OBJECT_SCIENCE=[obj],
        OBJECT_COMPARISON=[obj],
        OBJECT_TEFF=[int(round(min(max(float(teff), 3000.0), 6000.0)))],
        RUN_LBL_RESET=False, RUN_LBL_TELLUCLEAN=tellu,
        RUN_LBL_TEMPLATE=True, RUN_LBL_MASK=True, RUN_LBL_COMPUTE=True,
        RUN_LBL_COMPILE=True, SKIP_LBL_TEMPLATE=True, SKIP_LBL_MASK=True,
        SKIP_LBL_COMPUTE=True, SKIP_LBL_COMPILE=False,
        RESPROJ_TABLES={f'DTEMP{near}': f'temperature_gradient_{near}.fits'},
        COMPUTE_NPROC=max(1, int(nproc)), PLOT=False)


def wrap(params: Dict[str, Any]) -> str:
    """the script that runs LBL with its parameters (a wrap file)"""
    return ('#!/usr/bin/env python\n'
            '"""LBL, run by koloa (koloa.lblrun): its parameters as koloa '
            'wrote them"""\n'
            'import json\n\n'
            'from lbl import lbl_wrap\n\n'
            f'rparams = json.loads({json.dumps(json.dumps(params, indent=1))})'
            '\n\n'
            'if __name__ == "__main__":\n'
            '    lbl_wrap(rparams)\n')


def rdb_of(params: Dict[str, Any]) -> str:
    """where LBL writes the rdb of a run (with the telluric cleaning, the
    star's name ends in _tc)"""
    obj = params['OBJECT_SCIENCE'][0] + (
        '_tc' if params.get('RUN_LBL_TELLUCLEAN') else '')
    return os.path.join(params['DATA_DIR'], 'lblrdb',
                        f'lbl_{obj}_{obj}.rdb')


def run(params: Dict[str, Any], python: str,
        echo: Optional[Callable[[str], None]] = None) -> str:
    """
    LBL, with its Python, in a process of its own: its lines said as they
    come (echo, else printed), the script kept beside its data

    :return: str, the rdb it made
    """
    script = os.path.join(params['DATA_DIR'],
                          f'koloa_wrap_{params["OBJECT_SCIENCE"][0]}.py')
    with open(script, 'w') as handle:
        handle.write(wrap(params))
    echo = echo or (lambda line: print(line, flush=True))
    colour = re.compile(r'\x1b\[[0-9;]*m')
    last: List[str] = []
    env = dict(os.environ, PYTHONUNBUFFERED='1')
    with subprocess.Popen([python, script], cwd=params['DATA_DIR'],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, bufsize=1, env=env) as proc:
        for line in proc.stdout:
            # its words, without the colours of a terminal
            line = colour.sub('', line).rstrip()
            if line:
                echo(line)
                last = (last + [line])[-5:]
        code = proc.wait()
    if code != 0:
        raise RuntimeError(f'LBL stopped (exit {code}): '
                           + ' / '.join(last[-2:]))
    rdb = rdb_of(params)
    if not os.path.exists(rdb):
        raise RuntimeError(f'LBL made no rdb ({rdb})')
    return rdb


def keep(rdb: str, target: str, inst: str, root: str = 'archives') -> str:
    """the rdb of LBL in the folder of the star, under its name
    (archives/GJ_436/lbl/GJ_436_HARPSN_lbl.rdb)"""
    from koloa.gather import folder_name
    name = folder_name(target)
    path = os.path.join(root, name, 'lbl', f'{name}_{inst}_lbl.rdb')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    shutil.copyfile(rdb, path)
    return path


def lbl_star(target: str, instruments: Sequence[str], lbl_dir: str,
             root: str = 'archives', python: Optional[str] = None,
             nproc: Optional[int] = None, teff: Optional[float] = None,
             tellu: Optional[bool] = None, get: bool = True,
             compute: bool = True) -> Dict[str, Any]:
    """
    The LBL velocities of a star from the public spectra of the archives:
    for each instrument, its spectra fetched (get), LBL run on them
    (compute), its rdb kept in the folder of the star

    :param target: str, a name SIMBAD knows
    :param instruments: the instruments (keys of INSTRUMENTS)
    :param lbl_dir: str, the folder of LBL's data (spectra, templates...)
    :param root: str, the folder of the folders of the stars
    :param python: str or None, the Python of LBL (find_python)
    :param nproc: int or None, the processes of LBL's compute (half the
                  cores of the machine when not said)
    :param teff: float or None, the temperature of the star [K]
                 (koloa.stars.stellar_teff when not said)
    :param tellu: bool or None, clean the telluric lines (the instrument's
                  default when None)

    :return: dict, per instrument: n (spectra), rdb (the file kept), or
             error
    """
    unknown = [inst for inst in instruments if inst not in INSTRUMENTS]
    if unknown:
        raise ValueError(f'no instrument {unknown} for LBL: '
                         f'{", ".join(INSTRUMENTS)}')
    python = find_python(python)
    if compute:
        state = check_python(python)
        if not state['ok']:
            raise RuntimeError(f'LBL cannot be run: {state["error"]}')
        log(f'LBL {state["version"]} ({state["where"]}), with {python}',
            'info')
    if teff is None and compute:
        from koloa.stars import stellar_teff
        ident = _ident(target)
        guess = stellar_teff(ident)
        teff = guess.get('teff') or 4000.0
        log(f'{target}: Teff {teff:.0f} K ({guess.get("source") or "a guess"})',
            'value')
    nproc = int(nproc or settings().get('nproc') or
                max(1, (os.cpu_count() or 2) // 2))
    out: Dict[str, Any] = {}
    for inst in instruments:
        label = INSTRUMENTS[inst]['label']
        try:
            if get:
                step(f'{label}: the spectra')
                got = fetch(target, inst, lbl_dir, root)
                outcome(f'{got["n"]} spectra ({got["fetched"]} fetched'
                        + (f', {got["failed"]} not' if got['failed'] else '')
                        + ')')
                if not got['n']:
                    out[inst] = dict(n=0, error='no spectrum')
                    continue
            if compute:
                step(f'{label}: LBL')
                params = parameters(target, inst, lbl_dir, teff, nproc, tellu)
                start = _time.time()
                rdb = run(params, python)
                path = keep(rdb, target, inst, root)
                outcome(f'{path} ({_time.time() - start:.0f} s)')
                out[inst] = dict(rdb=path)
        except (RuntimeError, ValueError, OSError) as err:
            log(f'{label}: {err}', 'warn')
            out[inst] = dict(error=str(err))
    return out

# =============================================================================
# End of code
# =============================================================================

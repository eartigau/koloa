#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa's GUI: a page in the browser, served from this machine.

    koloa --gui                 # or python -m koloa.gui (--port 8765)

A star by its SIMBAD name (the resolver: identifiers, position, TIC, the
periods of variability SIMBAD lists, CARMENES DR1), a file of velocities or
none, the velocities by instrument, and the two long actions of koloa:
gathering the archives of the star (koloa.gather) and the detailed report
(koloa.detailed). Each action runs the koloa command line the page shows,
word for word, in a process of its own: what the page does, a batch does
with the same line. The steps of a run (the 'step:' lines of koloa's log)
show as they go, with the time each took, and the log as it comes.

The server listens on 127.0.0.1 only, and runs nothing but koloa, its
arguments passed as such (no shell).

Created on 2026-10-01

@author: artigau
"""
import argparse
import json
import mimetypes
import os
import shlex
import signal
import subprocess
import sys
import threading
import time
import urllib.parse
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, List, Optional

import numpy as np

# =============================================================================
# Define variables
# =============================================================================
#: the page and its files
STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'gui_static')
#: the defaults of koloa's command line (a flag is written only when its
#: value differs)
DEFAULTS = dict(kmax=3, nsweep=2000, nburn=400, pmin=1.1)
#: the runs of this session
JOBS: Dict[str, 'Job'] = {}


# =============================================================================
# The command lines
# =============================================================================
def _number(value: Any, kind=float) -> Optional[float]:
    """a number from a field of the page, None when empty"""
    if value is None or str(value).strip() == '':
        return None
    return kind(str(value).strip())


def _files(opts: Dict[str, Any]):
    """the files of the page (files: path and instrument of each; or file,
    one path) and their instruments ('' for the name the file gives)"""
    files = opts.get('files')
    if files is None:
        files = [dict(path=opts.get('file'))] if opts.get('file') else []
    files = [item for item in files if str(item.get('path') or '').strip()]
    return ([str(item['path']).strip() for item in files],
            [str(item.get('label') or '').strip() for item in files])


def command(action: str, opts: Dict[str, Any]) -> List[str]:
    """
    The arguments of koloa's command line for an action of the page

    :param action: str, gather or detailed
    :param opts: dict, the fields of the page (target, file, root, outdir,
                 the switches and numbers)

    :return: list of str, the arguments after 'koloa'
    """
    target = str(opts.get('target') or '').strip()
    paths, labels = _files(opts)
    off = {key: not opts.get(key, True) for key in
           ('dace', 'carmenes', 'tess', 'vizier', 'archive', 'gpcheck',
            'duck', 'latex', 'trend')}
    # the GP of the FIP: banded, sho (at the rotation given) or none
    fipgp = opts.get('fip_gp', 'banded')
    fipgp = {True: 'banded', False: 'none'}.get(fipgp, fipgp) or 'banded'
    rotation = _number(opts.get('rotation'))
    if action == 'archive':
        return ['--refresh-archive']
    if action == 'gather':
        if not target:
            raise ValueError('a SIMBAD name to gather the archives of')
        args = [target, '--gather', str(opts.get('root') or 'archives')]
        args += ['--no-dace'] * off['dace'] + ['--no-carmenes'] * off[
            'carmenes'] + ['--no-tess'] * off['tess']
        if opts.get('refresh'):
            args.append('--refresh')
        return args
    if action != 'detailed':
        raise ValueError(f'no action {action}')
    if not target and not paths:
        raise ValueError('a file, a SIMBAD name, or (best) both')
    # the velocities of the archives: only when asked
    ask = {key: bool(opts.get(key, False)) for key in
           ('dace', 'carmenes', 'vizier')}
    if not paths and not any(ask.values()):
        raise ValueError('no file: tick DACE, CARMENES DR1 or VizieR for '
                         'the velocities of the archives')
    if fipgp == 'sho' and not rotation:
        raise ValueError('an SHO GP needs the rotation period')
    args = paths + ['--detailed']
    if any(labels):
        args += ['--instruments'] + [lab or 'auto' for lab in labels]
    if target:
        args += ['--target', target]
    args += ['--outdir', str(opts.get('outdir') or 'koloa_output')]
    for key, kind in (('kmax', int), ('nsweep', int), ('nburn', int),
                      ('pmin', float), ('pmax', float)):
        val = _number(opts.get(key), kind)
        if val is not None and val != DEFAULTS.get(key):
            args += [f'--{key}', f'{val:g}' if kind is float else str(val)]
    periods = str(opts.get('periods') or '').replace(',', ' ').split()
    if periods:
        args += ['--periods'] + [f'{float(per):g}' for per in periods]
    if rotation:
        args += ['--rotation', f'{rotation:g}']
    if fipgp == 'none':
        args.append('--no-fip-gp')
    elif rotation and fipgp == 'banded':
        args += ['--fip-gp', 'banded']
    exclude = str(opts.get('exclude') or '').replace(',', ' ').split()
    if exclude:
        args += ['--exclude'] + exclude
    dmap = opts.get('detection_map') or 'none'
    if dmap == 'fip':
        args.append('--detection-map')
    elif dmap == 'search':
        args.append('--search-map')
    if opts.get('exposures'):
        args.append('--exposures')
    if opts.get('mcmc'):
        args.append('--mcmc')
    if opts.get('toi_on'):
        args += ['--toi'] + str(opts.get('tois') or '').replace(
            ',', ' ').replace('TOI-', '').split()
    if opts.get('curvature'):
        args.append('--curvature')
    elif off['trend']:
        args.append('--no-trend')
    args += [f'--{key}' for key, val in ask.items() if val]
    for key, flag in (('tess', '--no-tess'), ('archive', '--no-archive'),
                      ('gpcheck', '--no-gp'), ('duck', '--no-duck'),
                      ('latex', '--no-latex')):
        if off[key]:
            args.append(flag)
    return args


def line(args: List[str]) -> str:
    """the command line as it would be typed"""
    return shlex.join(['koloa'] + args)


# =============================================================================
# The runs
# =============================================================================
class Job:
    """
    One run of koloa's command line, its log read as it comes
    """

    def __init__(self, action: str, args: List[str], outputs: str):
        self.id = uuid.uuid4().hex[:8]
        self.action, self.args = action, args
        self.line = line(args)
        self.outputs = os.path.abspath(outputs)
        self.lines: List[str] = []
        self.steps: List[Dict[str, Any]] = []
        self.start, self.end = time.time(), None
        self.returncode: Optional[int] = None
        # the koloa of this page, whatever the folder it runs from
        env = dict(os.environ, PYTHONUNBUFFERED='1', KOLOA_PROGRESS='gui')
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env['PYTHONPATH'] = root + os.pathsep + env.get('PYTHONPATH', '')
        self.proc = subprocess.Popen(
            [sys.executable, '-W', 'ignore', '-m', 'koloa.cli'] + args,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            bufsize=1, env=env)
        threading.Thread(target=self._read, daemon=True).start()

    def _read(self):
        """the log, line by line, and the steps it marks"""
        for text in self.proc.stdout:
            text = text.rstrip('\n')
            message = text.split(' | ', 1)[-1]
            if message.startswith('progress: '):
                # the sweeps of a FIP: a bar on the page, not a log line
                part = message[10:].rsplit(' | ', 3)
                if self.steps and len(part) == 4:
                    self.steps[-1]['progress'] = dict(
                        label=part[0], done=int(part[1]), total=int(part[2]),
                        seconds=float(part[3]), at=time.time())
                continue
            self.lines.append(text)
            if message.startswith('step: '):
                now = time.time()
                if self.steps and self.steps[-1]['end'] is None:
                    self.steps[-1]['end'] = now
                self.steps.append(dict(name=message[6:], start=now, end=None,
                                       detail=''))
            elif self.steps and text.strip():
                self.steps[-1]['detail'] = message[:160]
        self.returncode = self.proc.wait()
        self.end = time.time()
        if self.steps and self.steps[-1]['end'] is None:
            self.steps[-1]['end'] = self.end

    @property
    def status(self) -> str:
        """running, done, failed or stopped"""
        if self.returncode is None:
            return 'running'
        if self.returncode == 0:
            return 'done'
        return 'stopped' if self.returncode < 0 else 'failed'

    def files(self) -> List[str]:
        """what the run wrote: the PDFs, text and JSON files of its folder"""
        if not os.path.isdir(self.outputs):
            return []
        out = []
        for base, _, names in os.walk(self.outputs):
            for name in names:
                if name.startswith('.') or not name.endswith(
                        ('.pdf', '.txt', '.json', '.csv')):
                    continue
                out.append(os.path.relpath(os.path.join(base, name),
                                           self.outputs))
        return sorted(out)

    def state(self, since: int = 0) -> Dict[str, Any]:
        """what the page shows of the run"""
        now = time.time()
        report = [name for name in self.files()
                  if name.endswith('_report.pdf')]
        return dict(id=self.id, action=self.action, line=self.line,
                    status=self.status, returncode=self.returncode,
                    elapsed=(self.end or now) - self.start,
                    steps=[dict(step, elapsed=(step['end'] or now)
                                - step['start']) for step in self.steps],
                    nlines=len(self.lines), lines=self.lines[since:],
                    outputs=self.outputs,
                    files=self.files() if self.status != 'running' else [],
                    report=report[0] if report else None)


# =============================================================================
# What the page asks
# =============================================================================
def resolve_star(name: str, root: str = '', refresh: bool = False
                 ) -> Dict[str, Any]:
    """the SIMBAD resolver of the page: identifiers, position, TIC, the
    periods SIMBAD lists, CARMENES DR1, and the planets the NASA Exoplanet
    Archive knows (its copy kept here, koloa.archive)"""
    from koloa.archive import CACHE, host_name, known_planets, resolve, tois
    from koloa.gather import carmenes_star, folder_name, variability
    folder = folder_name(name)
    # the disk first: the star's archives folder (koloa.gather), then what
    #   an earlier question kept; the network only when neither has it
    gathered = os.path.join(root or 'archives', folder, 'target.json')
    kept = os.path.join(CACHE, 'stars', folder + '.json')
    out = None
    for path, where in ((gathered, 'the archives folder'),
                        (kept, 'the copy kept')):
        if not refresh and os.path.exists(path):
            with open(path) as handle:
                out = json.load(handle)
            out.update(disk=where, disk_path=path, disk_date=time.strftime(
                '%Y-%m-%d %H:%M', time.localtime(os.path.getmtime(path))))
            break
    if out is None:
        out = _ask_star(name, refresh)
        os.makedirs(os.path.dirname(kept), exist_ok=True)
        with open(kept, 'w') as handle:
            json.dump(out, handle, default=str)
    out.update(folder=folder, planets=[], tois=[])
    try:
        host = host_name(out)
        known = known_planets(host=host) if host else {}
        out['archive_host'] = host
        out['archive_rotation'] = (known.get('star') or {}).get('rotation')
        out['planets'] = [
            dict(name=pl['name'], P=pl.get('P'), K=pl.get('K'),
                 mass_earth=pl.get('mass_earth'), e=pl.get('e'),
                 reference=pl.get('reference'),
                 reference_url=pl.get('reference_url'),
                 solutions=len(pl.get('solutions') or []),
                 discovery=pl.get('discovery'), year=pl.get('disc_year'))
            for pl in known.get('planets', [])]
    except Exception as err:  # a help, not a need
        out['planets_error'] = str(err)
    if out.get('tic'):
        try:
            out['tois'] = tois(out['tic'])
        except Exception as err:
            out['tois_error'] = str(err)
    return out


def _ask_star(name: str, refresh: bool = False) -> Dict[str, Any]:
    """what the network says of a star: SIMBAD (its names, its periods of
    variability) and CARMENES DR1"""
    from koloa.archive import resolve
    from koloa.gather import carmenes_star, variability
    out = dict(resolve(name, refresh=refresh), variability=[],
               carmenes=None)
    try:
        out['variability'] = variability(out['main'])
    except Exception as err:  # a help, not a need
        out['variability_error'] = str(err)
    if out.get('ra') is not None:
        try:
            out['carmenes'] = carmenes_star(out['ra'], out['dec'])
        except Exception as err:
            out['carmenes_error'] = str(err)
    return out


def archives(target: str, root: str = '') -> Dict[str, Any]:
    """what koloa.gather already put on disk for a star (its manifest), so
    that the page offers to refresh it rather than to gather it again"""
    from koloa.gather import folder_name
    folder = os.path.join(root or 'archives', folder_name(target))
    path = os.path.join(folder, 'manifest.json')
    if not target.strip() or not os.path.exists(path):
        return dict(exists=False, folder=folder)
    with open(path) as handle:
        manifest = json.load(handle)
    return dict(exists=True, folder=folder, created=manifest.get('created'),
                archives={key: val.get('status') for key, val in
                          manifest.get('archives', {}).items()})


#: what a dialog of this machine says, by what it picks
PROMPTS = dict(file='A file of velocities (LBL .rdb, csv, DACE csv)',
               folder='A folder')


def pick(kind: str = 'file', start: str = '') -> Dict[str, Any]:
    """
    A dialog of this machine to choose a file or a folder (the page cannot:
    a browser never gives it the path of a file): Finder's on macOS
    (osascript), Tk's elsewhere, in a process of its own

    :param kind: str, file or folder
    :param start: str, where the dialog opens (the folder of a path given,
                  or the folder koloa runs from)

    :return: dict, path, or cancelled
    """
    start = os.path.abspath(os.path.expanduser(start)) if start else ''
    while start and not os.path.isdir(start):
        start = os.path.dirname(start)
    start = start or os.getcwd()
    what = 'folder' if kind == 'folder' else 'file'
    if sys.platform == 'darwin':
        place = start.replace('\\', '\\\\').replace('"', '\\"')
        proc = subprocess.run(
            ['osascript', '-e', 'activate', '-e',
             f'POSIX path of (choose {what} with prompt "{PROMPTS[what]}" '
             f'default location (POSIX file "{place}"))'],
            capture_output=True, text=True)
        if proc.returncode != 0:
            if '-128' in proc.stderr:
                return dict(cancelled=True)
            raise RuntimeError(proc.stderr.strip() or 'osascript failed')
        path = proc.stdout.strip()
    else:
        code = ('import sys, tkinter as tk\n'
                'from tkinter import filedialog\n'
                'root = tk.Tk(); root.withdraw()\n'
                "root.attributes('-topmost', True)\n"
                + ("print(filedialog.askdirectory(initialdir=sys.argv[1], "
                   "title=sys.argv[2]))" if what == 'folder' else
                   "print(filedialog.askopenfilename(initialdir=sys.argv[1], "
                   "title=sys.argv[2], filetypes=[('velocities', '*.rdb "
                   "*.csv *.dat *.txt'), ('all', '*')]))"))
        proc = subprocess.run([sys.executable, '-c', code, start,
                               PROMPTS[what]], capture_output=True,
                              text=True)
        if proc.returncode != 0:
            raise RuntimeError('no dialog on this machine (Tk): type the '
                               'path')
        path = proc.stdout.strip()
    if not path:
        return dict(cancelled=True)
    if what == 'folder':
        path = path.rstrip('/') or '/'
        # a folder under the one koloa runs from: written relative to it
        rel = os.path.relpath(path, os.getcwd())
        path = rel if not rel.startswith('..') else path
    return dict(path=path)


def gathering(target: str, root: str = '') -> bool:
    """whether the archives of a star are being gathered by a run"""
    from koloa.gather import folder_name
    folder = os.path.abspath(os.path.join(root or 'archives',
                                          folder_name(target)))
    return any(job.action == 'gather' and job.returncode is None
               and job.outputs == folder for job in JOBS.values())


def series_of(files: Any = '', target: str = '', root: str = '',
              dace: bool = False, carmenes: bool = False):
    """
    The series of the page: its files, and the archives gathered for the
    star that are asked for (set apart from the files as the report does)

    :return: tuple, RVData or None, the source of each instrument, notes
    """
    from koloa.data import merge
    from koloa.detailed import distinct, read_files
    from koloa.gather import folder_name, load
    series, notes, source = [], [], {}
    filedata = None
    if isinstance(files, str):
        files = [dict(path=files)] if files else []
    paths, labels = _files(dict(files=files))
    if paths:
        parts = read_files(paths, [lab or None for lab in labels])
        for path, part in zip(paths, parts):
            source.update({name: f'file: {os.path.basename(path)}'
                           for name in part.instruments})
            notes.append(f'{os.path.basename(path)}: {part.n} points')
        series += parts
        filedata = parts[0] if len(parts) == 1 else merge(parts)
    if target:
        folder = os.path.join(root or 'archives', folder_name(target))
        if os.path.exists(os.path.join(folder, 'rv', 'all_rv.csv')):
            gathered = load(folder)['rv']
            # each archive set apart from the file, as the report does
            for arch, tag, sel, asked in (
                    ('DACE', 'DACE',
                     ~np.char.startswith(gathered.inst.astype(str), 'CARM'),
                     dace),
                    ('CARMENES DR1', 'DR1',
                     np.char.startswith(gathered.inst.astype(str), 'CARM'),
                     carmenes)):
                # only the archives the report will use
                if not asked or not np.any(sel):
                    continue
                part = gathered.select(sel)
                if filedata is not None:
                    part = distinct(part, filedata.instruments,
                                    filedata.time, tag)
                if part is None:
                    continue
                series.append(part)
                source.update({name: arch for name in part.instruments})
                notes.append(f'{arch} ({folder}): {part.n} points')
        else:
            notes.append(f'nothing gathered in {folder} yet')
    if not series:
        return None, source, notes
    return (series[0] if len(series) == 1 else merge(series)), source, notes


def velocities(files: Any = '', target: str = '', root: str = '',
               dace: bool = False, carmenes: bool = False
               ) -> Dict[str, Any]:
    """
    The velocities of a file and of a star's gathered archives, by
    instrument (each with its median taken out), for the plot of the page
    """
    data, source, notes = series_of(files, target, root, dace, carmenes)
    if data is None:
        return dict(instruments=[], notes=notes)
    out = []
    for name in data.instruments:
        sel = data.inst == name
        # each instrument about its own median, whatever came before
        rv = data.rv[sel] - np.median(data.rv[sel])
        out.append(dict(name=name, n=int(sel.sum()),
                        source=source.get(name, ''),
                        time=np.round(data.time[sel], 6).tolist(),
                        rv=np.round(rv, 3).tolist(),
                        err=np.round(data.err[sel], 3).tolist(),
                        rms=float(np.std(data.rv[sel]))))
    return dict(instruments=out, notes=notes, n=int(data.n),
                baseline=float(data.baseline))


# =============================================================================
# The quick look: a FIP of what the page shows, before the report
# =============================================================================
#: a FIP without a GP, short: a look at the series before the report
QUICK = dict(kmax=2, nsweep=500, nburn=200)
#: the periods of the window of a series from the ground [days]
WINDOW = dict(day=0.99727, month=29.5306, year=365.25)
#: the quick looks of this session
QUICKS: Dict[str, Dict[str, Any]] = {}
#: one quick FIP at a time (they share the progress hook of koloa.fip)
_QUICK_LOCK = threading.Lock()


def selection(opts: Dict[str, Any]):
    """the series the page shows, its instruments left out taken out"""
    data, source, notes = series_of(
        opts.get('files') or opts.get('file') or '', opts.get('target', ''),
        opts.get('root', ''), bool(opts.get('dace')),
        bool(opts.get('carmenes')))
    left = {name.upper() for name in str(opts.get('exclude') or '')
            .replace(',', ' ').split()}
    if data is not None and left:
        keep = ~np.isin(np.char.upper(data.inst.astype(str)), list(left))
        data = data.select(keep) if keep.any() else None
    return data, source, notes


def known_periods(target: str) -> List[Dict[str, Any]]:
    """the known planets of a star (the archive's copy kept here)"""
    if not target.strip():
        return []
    from koloa.archive import host_name, known_planets, resolve
    try:
        host = host_name(resolve(target))
        return [dict(name=pl['name'], P=float(pl['P']))
                for pl in (known_planets(host=host)['planets'] if host
                           else []) if pl.get('P')]
    except Exception:  # a help, not a need
        return []


def _fip_curves(res, nbin: int = 3000) -> Dict[str, Any]:
    """the FIP of the period alone and of the period or any of its aliases,
    -log10, kept at their peaks in nbin bins of log period"""
    per = 1.0 / np.asarray(res.freq)
    alone = np.asarray(res.fip)
    fam = np.asarray(res.family if res.family is not None else res.fip)
    logp = np.log10(per)
    edges = np.linspace(logp.min(), logp.max() + 1e-9, nbin + 1)
    idx = np.clip(np.digitize(logp, edges) - 1, 0, nbin - 1)
    out = {}
    for key, val in (('alone', alone), ('family', fam)):
        best = np.ones(nbin)
        np.minimum.at(best, idx, np.clip(val, 1e-15, 1.0))
        out[key] = np.round(-np.log10(best), 4).tolist()
    out['period'] = np.round(10 ** (0.5 * (edges[1:] + edges[:-1])),
                             6).tolist()
    return out


def _run_quick(qid: str, data, target: str):
    """the quick FIP, in a thread, in the two passes of the report without
    its GP: the errors of each instrument inflated to the noise of a fit
    without planets, a first FIP; then to the noise of a fit with the
    signals it found and the known planets (their variance is not noise),
    the FIP again"""
    from koloa import fip as kfip
    from koloa.detailed import _fip
    from koloa.fit import RVModel
    job = QUICKS[qid]
    with _QUICK_LOCK:
        try:
            nights = data.nightly()
            seq_jitter = ('instrument' if len(nights.instruments) > 1
                          and nights.nseq < nights.n else None)
            job['step'] = 'noise'
            noise = RVModel(nights, [], likelihood='mixture', unit='both',
                            trend=1, seq_jitter=seq_jitter).fit(
                nstart=2, quiet=True)
            job['step'] = 'fip'

            def hook(label, done, total, seconds):
                job['progress'] = dict(done=done, total=total,
                                       seconds=seconds)
            kfip.PROGRESS_HOOK = hook
            job['step'] = 'fip1'
            res, info = _fip(nights, noise, QUICK['kmax'], QUICK['nsweep'],
                             QUICK['nburn'], 1, 'quick FIP, first pass',
                             gp=None, trend=1)
            width = 1.0 / nights.baseline
            known = known_periods(target)
            # the second pass: the noise of a fit with the signals found and
            #   the known planets
            pers = []
            for per in sorted([float(pk['period']) for pk in res.peaks
                               if res.family_containing(pk['period'], width)
                               < 0.01]
                              + [pl['P'] for pl in known
                                 if 1.0 < pl['P'] < nights.baseline]):
                if all(abs(per / old - 1) > 0.02 for old in pers):
                    pers.append(per)
            passes = 1
            if pers:
                job['step'] = 'planets'
                fit = RVModel(nights, [dict(period=per, period_range=(
                    0.98 * per, 1.02 * per)) for per in pers],
                    likelihood='mixture', unit='both', trend=1,
                    seq_jitter=seq_jitter).fit(nstart=2, quiet=True)
                job['step'] = 'fip2'
                job['progress'] = None
                res, info = _fip(nights, fit, QUICK['kmax'], QUICK['nsweep'],
                                 QUICK['nburn'], 2, 'quick FIP, second pass',
                                 gp=None, trend=1)
                passes = 2
            job['result'] = dict(
                _fip_curves(res), known=known, window=WINDOW, passes=passes,
                planets=pers,
                peaks=[dict(period=float(pk['period']), fip=float(pk['fip']),
                            family=float(res.family_containing(
                                pk['period'], width)))
                       for pk in res.peaks[:5]],
                n=int(nights.n), nexp=int(data.n),
                instruments={str(inst): int(np.sum(nights.inst == inst))
                             for inst in nights.instruments},
                inflation={str(key): float(val) for key, val in
                           info.get('inflation', {}).items()},
                settings=dict(QUICK, gp='none', trend=1))
            job['status'] = 'done'
        except Exception as err:
            job['status'] = 'failed'
            job['error'] = f'{type(err).__name__}: {err}'
        finally:
            kfip.PROGRESS_HOOK = None
            job['end'] = time.time()


def quick_fip(opts: Dict[str, Any]) -> Dict[str, Any]:
    """a quick FIP of what the page shows, started in a thread"""
    data, _, _ = selection(opts)
    if data is None:
        raise ValueError('no velocities to look at')
    qid = uuid.uuid4().hex[:8]
    QUICKS[qid] = dict(id=qid, status='running', step='waiting',
                       progress=None, result=None, error=None,
                       start=time.time(), end=None,
                       instruments=list(data.instruments))
    threading.Thread(target=_run_quick, args=(qid, data,
                                              opts.get('target', '')),
                     daemon=True).start()
    return quick_state(qid)


def quick_state(qid: str) -> Dict[str, Any]:
    """what the page shows of a quick FIP"""
    job = QUICKS[qid]
    return dict(job, elapsed=(job['end'] or time.time()) - job['start'])


def quicklook_pdf(opts: Dict[str, Any], xr=None, yr=None, pr=None,
                  qid: str = '') -> bytes:
    """
    The quick look as a PDF: the velocities shown (each instrument about
    its median, the ranges of the page) over the quick FIP

    :return: bytes, the PDF
    """
    import io
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from koloa import plotting as kplot
    data, source, _ = selection(opts)
    if data is None:
        raise ValueError('no velocities to draw')
    kplot.set_style('paper')
    quick = (QUICKS.get(qid) or {}).get('result')
    fig, axes = plt.subplots(2 if quick else 1, 1,
                             figsize=(10, 7.5 if quick else 4.2),
                             squeeze=False)
    ax = axes[0, 0]
    for it, inst in enumerate(data.instruments):
        sel = data.inst == inst
        ax.errorbar(data.time[sel], data.rv[sel] - np.median(data.rv[sel]),
                    data.err[sel], fmt=kplot.INST_MARKERS[it % 8], ms=3.5,
                    lw=0.6, color=kplot.INST_COLOURS[it % 8],
                    label=f'{inst} ({source.get(inst, "")}, {int(sel.sum())})')
    ax.axhline(0, color='0.6', lw=0.6, ls=':')
    if xr:
        ax.set_xlim(*xr)
    if yr:
        ax.set_ylim(*yr)
    ax.set_xlabel('BJD - 2400000')
    ax.set_ylabel('RV - median [m s$^{-1}$]')
    ax.legend(fontsize=7, ncol=3, frameon=False, loc='upper left')
    title = (opts.get('target') or '').strip() or 'the series'
    ax.set_title(f'{title}: the velocities shown', fontsize=10)
    if quick:
        ax = axes[1, 0]
        per = np.asarray(quick['period'])
        ax.plot(per, quick['alone'], color='0.6', lw=0.8,
                label='the period alone')
        ax.plot(per, quick['family'], color=kplot.C['koloa'], lw=1.0,
                label='the period or any of its aliases')
        ax.axhline(2, color='0.3', lw=0.8, ls=':', label='FIP = 1 %')
        lo, hi = pr if pr else (per.min(), per.max())
        # the window and the known planets, those within the periods shown
        for name, val in quick['window'].items():
            if lo <= val <= hi:
                ax.axvline(val, color='0.5', lw=0.8, ls=':')
                ax.text(val, 0.02, f' {name}', fontsize=7, color='0.4',
                        transform=ax.get_xaxis_transform(), va='bottom')
        for pl in quick['known']:
            if lo <= pl['P'] <= hi:
                ax.axvline(pl['P'], color=kplot.C['outlier'], lw=0.9,
                           ls='--')
                ax.text(pl['P'], 0.80, f' {pl["name"]}', fontsize=7,
                        transform=ax.get_xaxis_transform(), va='top',
                        color=kplot.C['outlier'])
        ax.set_xscale('log')
        kplot.plain_log_ticks(ax, 'x')
        ax.set_xlim(lo, hi)
        ax.set_xlabel('period [d]')
        ax.set_ylabel('-log$_{10}$ FIP')
        ax.legend(fontsize=7, frameon=False, loc='upper right')
        sett = quick['settings']
        ax.set_title(f'quick FIP: {quick["n"]} nights, '
                     f'{", ".join(quick["instruments"])}; no GP, '
                     f'{sett["kmax"]} signals, {sett["nsweep"]} sweeps, '
                     f'{quick["passes"]} pass(es)', fontsize=9)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='pdf')
    plt.close(fig)
    return buf.getvalue()


# =============================================================================
# The server
# =============================================================================
class Handler(BaseHTTPRequestHandler):
    """the page, its files, and its questions (/api/...)"""

    def log_message(self, *args):  # quiet: the page is the log
        return

    def _send(self, code: int, body: bytes, kind: str):
        self.send_response(code)
        self.send_header('Content-Type', kind)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def _json(self, value: Any, code: int = 200):
        self._send(code, json.dumps(value, default=str).encode(),
                   'application/json')

    def _file(self, path: str):
        if not os.path.isfile(path):
            return self._json(dict(error='no such file'), 404)
        kind = mimetypes.guess_type(path)[0] or 'application/octet-stream'
        if kind.startswith('text/') or kind.endswith(('javascript', 'json')):
            kind += '; charset=utf-8'
        with open(path, 'rb') as handle:
            self._send(200, handle.read(), kind)

    def do_GET(self):
        url = urllib.parse.urlparse(self.path)
        query = {key: val[-1] for key, val in
                 urllib.parse.parse_qs(url.query).items()}
        try:
            if url.path in ('/', '/index.html'):
                return self._file(os.path.join(STATIC, 'index.html'))
            if url.path.startswith('/static/'):
                name = os.path.basename(url.path)
                return self._file(os.path.join(STATIC, name))
            if url.path == '/api/info':
                from koloa.archive import fetched
                return self._json(dict(cwd=os.getcwd(),
                                       python=sys.executable,
                                       defaults=DEFAULTS,
                                       archive=fetched()))
            if url.path == '/api/resolve':
                return self._json(resolve_star(
                    query.get('name', ''), query.get('root', ''),
                    query.get('refresh', '') == '1'))
            if url.path == '/api/pick':
                return self._json(pick(query.get('kind', 'file'),
                                       query.get('start', '')))
            if url.path == '/api/quickfip':
                return self._json(quick_state(query['id']))
            if url.path == '/api/archives':
                state = archives(query.get('target', ''),
                                 query.get('root', ''))
                state['busy'] = gathering(query.get('target', ''),
                                          query.get('root', ''))
                return self._json(state)
            if url.path == '/api/rv':
                files = (json.loads(query['files']) if query.get('files')
                         else query.get('file', ''))
                return self._json(velocities(
                    files, query.get('target', ''), query.get('root', ''),
                    dace=query.get('dace') == '1',
                    carmenes=query.get('carmenes') == '1'))
            if url.path == '/api/jobs':
                return self._json([job.state(0) for job in JOBS.values()])
            if url.path == '/api/job':
                job = JOBS[query['id']]
                return self._json(job.state(int(query.get('since', 0))))
            if url.path == '/api/output':
                # a file a run wrote, and only that
                job = JOBS[query['id']]
                path = os.path.abspath(os.path.join(job.outputs,
                                                    query['name']))
                if not path.startswith(job.outputs + os.sep):
                    return self._json(dict(error='outside the run'), 403)
                return self._file(path)
            return self._json(dict(error='not found'), 404)
        except Exception as err:
            return self._json(dict(error=f'{type(err).__name__}: {err}'),
                              400)

    def do_POST(self):
        size = int(self.headers.get('Content-Length', 0))
        try:
            body = json.loads(self.rfile.read(size) or b'{}')
            path = urllib.parse.urlparse(self.path).path
            if path == '/api/command':
                args = command(body['action'], body.get('options', {}))
                return self._json(dict(args=args, line=line(args)))
            if path == '/api/run':
                opts = body.get('options', {})
                args = command(body['action'], opts)
                if body['action'] == 'gather':
                    from koloa.gather import folder_name
                    outputs = os.path.join(opts.get('root') or 'archives',
                                           folder_name(opts['target']))
                    busy = [job for job in JOBS.values()
                            if job.action == 'gather'
                            and job.returncode is None and job.outputs
                            == os.path.abspath(outputs)]
                    if busy:
                        raise ValueError('the archives of this star are '
                                         'being gathered already')
                elif body['action'] == 'archive':
                    from koloa.archive import CACHE
                    outputs = CACHE
                else:
                    outputs = opts.get('outdir') or 'koloa_output'
                job = Job(body['action'], args, outputs)
                JOBS[job.id] = job
                return self._json(job.state())
            if path == '/api/quickfip':
                return self._json(quick_fip(body.get('options', {})))
            if path == '/api/quicklook_pdf':
                pdf = quicklook_pdf(body.get('options', {}), body.get('x'),
                                    body.get('y'), body.get('p'),
                                    body.get('quick', ''))
                self.send_response(200)
                self.send_header('Content-Type', 'application/pdf')
                self.send_header('Content-Length', str(len(pdf)))
                self.send_header('Content-Disposition',
                                 'attachment; filename="koloa_quicklook.pdf"')
                self.end_headers()
                self.wfile.write(pdf)
                return None
            if path == '/api/stop':
                job = JOBS[body['id']]
                if job.returncode is None:
                    job.proc.terminate()
                return self._json(dict(id=job.id))
            return self._json(dict(error='not found'), 404)
        except Exception as err:
            return self._json(dict(error=f'{type(err).__name__}: {err}'),
                              400)


def serve(port: int = 8765, browser: bool = True):
    """
    The GUI: a server on 127.0.0.1 and the page in the browser (the next
    free port when this one is taken); Ctrl-C stops it, and the runs with it

    :param port: int, the first port tried
    :param browser: bool, open the page in the browser
    """
    for trial in range(port, port + 20):
        try:
            server = ThreadingHTTPServer(('127.0.0.1', trial), Handler)
            break
        except OSError:
            continue
    else:
        raise OSError(f'no free port from {port} to {port + 19}')
    url = f'http://127.0.0.1:{server.server_address[1]}/'
    print(f'koloa GUI: {url} (runs from {os.getcwd()}; Ctrl-C to stop)',
          flush=True)
    if browser:
        threading.Timer(0.5, webbrowser.open, args=(url,)).start()

    def stop(*_):
        raise KeyboardInterrupt
    # a kill stops the runs too, as Ctrl-C does
    signal.signal(signal.SIGTERM, stop)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        for job in JOBS.values():
            if job.returncode is None:
                job.proc.terminate()
        server.server_close()


def main(argv=None):
    """python -m koloa.gui [--port 8765] [--no-browser]"""
    parser = argparse.ArgumentParser(prog='koloa.gui',
                                     description="koloa's GUI, in the browser")
    parser.add_argument('--port', type=int, default=8765)
    parser.add_argument('--no-browser', action='store_true')
    args = parser.parse_args(argv)
    serve(args.port, browser=not args.no_browser)


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

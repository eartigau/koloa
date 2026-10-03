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

A quick look worth keeping is remembered (REMEMBERED: the page, its quick
FIP, the velocities plotted, and a copy of the files and of the archives'
velocities), listed in a tab of its own, and recalled as it was.

The server listens on 127.0.0.1 only, and runs nothing but koloa, its
arguments passed as such (no shell).

Created on 2026-10-01

@author: artigau
"""
import argparse
import hashlib
import json
import mimetypes
import os
import re
import shlex
import shutil
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
    # SIMBAD's periods of variability: asked in the background (its TAP
    #   service is at times very slow), kept with the star when they come
    if out.get('variability') is None:
        _start_variability(out.get('main') or name, kept)
        out['variability'], out['variability_pending'] = [], True
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


#: the questions to SIMBAD's TAP service under way, by star
_VARIABILITY: Dict[str, threading.Thread] = {}


def _variability_later(main: str, kept: str):
    """SIMBAD's periods of variability of a star, asked in the background,
    written into the copy kept of the star when they come"""
    from koloa.gather import variability
    try:
        found = variability(main, timeout=180.0)
        with open(kept) as handle:
            star = json.load(handle)
        star['variability'] = found
        with open(kept, 'w') as handle:
            json.dump(star, handle, default=str)
    except Exception:  # no answer: asked again the next time
        return


def _start_variability(main: str, kept: str):
    """the periods of variability of a star asked, once at a time"""
    thread = _VARIABILITY.get(kept)
    if thread is not None and thread.is_alive():
        return
    thread = threading.Thread(target=_variability_later, args=(main, kept),
                              daemon=True)
    thread.start()
    _VARIABILITY[kept] = thread


def _ask_star(name: str, refresh: bool = False) -> Dict[str, Any]:
    """what the network says of a star at once: SIMBAD's names (Sesame)
    and CARMENES DR1 (its list kept here); the periods of variability come
    later (variability None until then)"""
    from koloa.archive import resolve
    from koloa.gather import carmenes_star
    out = dict(resolve(name, refresh=refresh), variability=None,
               carmenes=None)
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
                          manifest.get('archives', {}).items()},
                points={key: val.get('npoints') for key, val in
                        manifest.get('archives', {}).items()
                        if val.get('npoints')})


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
            # the velocities only (the photometry is not plotted)
            gathered = load(folder, photometry=False)['rv']
            # each archive set apart from the file, as the report does
            for arch, tag, sel, asked in (
                    ('DACE', 'DACE',
                     ~np.char.startswith(gathered.inst.astype(str), 'CARM'),
                     dace),
                    ('CARMENES DR1', 'DR1',
                     np.char.startswith(gathered.inst.astype(str), 'CARM'),
                     carmenes)):
                # only the archives the report will use
                if not np.any(sel):
                    continue
                if not asked:
                    notes.append(f'{arch}: {int(np.sum(sel))} points '
                                 f'gathered, not ticked')
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
#: the fewest nights for an instrument to have a quick FIP of its own
QUICK_MIN_NIGHTS = 10
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


#: a FIP of the period or any of its aliases below which a family is
#: certain: families below it are told apart by the FIP of the period alone
#: (1e-17 and 1e-76 say the same thing)
FAMILY_FLOOR = 1e-6


def _distinct_peaks(res, width: float, nmax: int = 8
                    ) -> List[Dict[str, Any]]:
    """the peaks of the FIP, one per family of aliases, the best first:
    within a family (all at the same FIP of the period or any of its
    aliases), the period whose own FIP is the lowest, the others its
    aliases; a peak is a dip of the FIP of the period alone (not a point of
    a flat stretch at 1, where the period itself has nothing)"""
    from koloa.aliases import same_family
    fam = np.asarray(res.family if res.family is not None else res.fip)
    alone = np.asarray(res.fip)
    freq = np.asarray(res.freq)
    low = np.where((alone[1:-1] <= alone[:-2]) & (alone[1:-1] <= alone[2:])
                   & (alone[1:-1] < 1.0 - 1e-6)
                   & (fam[1:-1] < 1.0))[0] + 1
    out: List[Dict[str, Any]] = []
    order = np.lexsort((alone[low], np.maximum(fam[low], FAMILY_FLOOR)))
    for idx in low[order]:
        per = float(1.0 / freq[idx])
        if any(same_family(old['period'], per, width) for old in out):
            continue
        out.append(dict(period=per,
                        family=float(res.family_containing(per, width)),
                        alone=float(res.fip_containing(per, width))))
        if len(out) >= nmax:
            break
    for rank, peak in enumerate(out):
        peak['id'] = rank + 1
    return out


def fold(data, period: float, valid: Optional[np.ndarray] = None,
         trend: int = 1) -> Dict[str, Any]:
    """
    The series folded at a period: a sinusoid fitted with an offset per
    instrument and a trend (weighted least squares, the errors of K scaled
    by the reduced chi^2 when it is above one), phase 0 at the conjunction
    (the velocity falling through zero, as in the report)

    :param valid: np.ndarray or None, the probability of each point to be
                  valid (not an outlier), from the FIP: each point weighs
                  that much more in the fit (an outlier hardly counts), and
                  it is given back with the points
    :param trend: int, the order of the trend fitted with it (0: none, 1: an
                  acceleration, 2: and its change)

    :return: dict, period, K, K_err, tc, rms, the points of each instrument
             (about its offset and the trend, with their probability to be
             valid when given) and the curve
    """
    time_ = data.time
    tref = float(np.median(time_))
    insts = list(data.instruments)
    cols = [(data.inst == inst).astype(float) for inst in insts]
    cols += [((time_ - tref) / 365.25) ** order
             for order in range(1, trend + 1)]
    arg = 2 * np.pi * (time_ - tref) / period
    cols += [np.cos(arg), np.sin(arg)]
    design = np.column_stack(cols)
    good = (np.ones(data.n) if valid is None
            else np.clip(np.asarray(valid, float), 1e-6, 1.0))
    wgt = good / data.err ** 2
    amat = design.T @ (design * wgt[:, None])
    coef = np.linalg.solve(amat, design.T @ (data.rv * wgt))
    cov = np.linalg.inv(amat)
    resid = data.rv - design @ coef
    chi2 = float(np.sum(resid ** 2 * wgt)) / max(np.sum(good) - len(coef),
                                                 1)
    cov *= max(chi2, 1.0)
    acos, asin = coef[-2], coef[-1]
    amp = float(np.hypot(acos, asin))
    grad = np.array([acos, asin]) / max(amp, 1e-12)
    kerr = float(np.sqrt(grad @ cov[-2:, -2:] @ grad))
    # the model is K cos(arg - phi0): its maximum at phi0, the conjunction
    #   a quarter of a period later
    phi0 = float(np.arctan2(asin, acos))
    tc = tref + period * (phi0 / (2 * np.pi) + 0.25)
    phase = ((time_ - tc) / period) % 1.0
    shown = data.rv - design[:, :-2] @ coef[:-2]
    grid = np.linspace(0, 1, 201)
    out = dict(period=float(period), K=amp, K_err=kerr, tc=float(tc),
               rms=float(np.std(resid[good >= 0.5])), chi2=chi2,
               curve=dict(phase=grid.tolist(),
                          rv=np.round(-amp * np.sin(2 * np.pi * grid),
                                      4).tolist()), instruments=[])
    for inst in insts:
        sel = data.inst == inst
        out['instruments'].append(dict(
            name=str(inst), phase=np.round(phase[sel], 5).tolist(),
            rv=np.round(shown[sel], 3).tolist(),
            err=np.round(data.err[sel], 3).tolist()))
        if valid is not None:
            out['instruments'][-1]['valid'] = np.round(valid[sel],
                                                       4).tolist()
    return out


def trend_order(opts: Dict[str, Any]) -> int:
    """the order of the trend the report fits, from its boxes: 0 without
    a trend, 1 an acceleration (the default), 2 with its change"""
    if opts.get('curvature'):
        return 2
    return 0 if opts.get('trend') is False else 1


def _quick_acceleration(fit, trend: int) -> Optional[Dict[str, Any]]:
    """the acceleration of the star (dv/dt) the fit of the quick look
    measured, and its change when fitted: (value, minus, plus) from the
    Laplace covariance of the maximum a posteriori, and its significance"""
    if trend < 1:
        return None
    from koloa.secular import acceleration
    try:
        acc = acceleration(fit)
    except (ValueError, np.linalg.LinAlgError):
        return None
    out = dict(tref=float(acc['tref']))
    for key in ('accel', 'jerk'):
        if key in acc:
            val, low, high = (float(x) for x in acc[key])
            out[key] = [val, low, high]
            out[key + '_sigma'] = abs(val) / max(0.5 * (low + high), 1e-30)
    return out


def _run_quick(qid: str, data, target: str, trend: int = 1):
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
                            trend=trend, seq_jitter=seq_jitter).fit(
                nstart=2, quiet=True)
            job['step'] = 'fip'

            def hook(label, done, total, seconds):
                job['progress'] = dict(done=done, total=total,
                                       seconds=seconds)
            kfip.PROGRESS_HOOK = hook
            job['step'] = 'fip1'
            res, info = _fip(nights, noise, QUICK['kmax'], QUICK['nsweep'],
                             QUICK['nburn'], 1, 'quick FIP, first pass',
                             gp=None, trend=trend)
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
                    likelihood='mixture', unit='both', trend=trend,
                    seq_jitter=seq_jitter).fit(nstart=2, quiet=True)
                job['step'] = 'fip2'
                job['progress'] = None
                res, info = _fip(nights, fit, QUICK['kmax'], QUICK['nsweep'],
                                 QUICK['nburn'], 2, 'quick FIP, second pass',
                                 gp=None, trend=trend)
                passes = 2
            # the peaks: those below a FIP of 10 %, or the best three, named
            #   and folded
            peaks = _distinct_peaks(res, width)
            named = [pk for pk in peaks if pk['family'] < 0.1]
            if len(named) < 3:
                named = peaks[:3]
            for pk in peaks:
                pk['named'] = pk in named
            # the probability of each night to be valid (not an outlier),
            #   as the FIP saw it
            valid = res.reliability
            if valid is not None and len(valid) != nights.n:
                valid = None
            folds = [dict(fold(nights, pk['period'], valid, trend),
                          id=pk['id'])
                     for pk in named]
            job['result'] = dict(
                _fip_curves(res), known=known, window=WINDOW, passes=passes,
                planets=pers, peak_list=peaks, folds=folds,
                pk=[float(val) for val in res.pk],
                peaks=[dict(period=float(pk['period']), fip=float(pk['fip']),
                            family=float(res.family_containing(
                                pk['period'], width)))
                       for pk in res.peaks[:5]],
                n=int(nights.n), nexp=int(data.n),
                instruments={str(inst): int(np.sum(nights.inst == inst))
                             for inst in nights.instruments},
                inflation={str(key): float(val) for key, val in
                           info.get('inflation', {}).items()},
                settings=dict(QUICK, gp='none', trend=trend),
                acceleration=_quick_acceleration(fit if pers else noise,
                                                 trend))
            # the FIP of each instrument on its own (its jitter sampled in
            #   the FIP, no inflation needed), when there are several
            insts = list(nights.instruments)
            job['each'] = []
            if len(insts) > 1:
                from koloa.fip import oafip
                from koloa.utils import blas_threads
                for rank, inst in enumerate(insts):
                    sub = nights.select(nights.inst == inst)
                    if sub.n < QUICK_MIN_NIGHTS:
                        job['each'].append(dict(name=str(inst), n=int(sub.n),
                                                skipped=True))
                        continue
                    job['step'] = 'inst'
                    job['step_detail'] = f'{rank + 1}/{len(insts)}: {inst}'
                    job['progress'] = None
                    with blas_threads(1):
                        one = oafip(sub, kmax=QUICK['kmax'], outliers='both',
                                    nsweep=QUICK['nsweep'],
                                    nburn=QUICK['nburn'], nchains=2,
                                    seed=3 + rank, progress=False, gp=None,
                                    trend=trend, nightly=False,
                                    label=f'quick FIP of {inst}')
                    job['each'].append(dict(
                        _fip_curves(one), name=str(inst), n=int(sub.n),
                        skipped=False,
                        peaks=_distinct_peaks(one, 1.0 / max(sub.baseline,
                                                             1.0), nmax=3)))
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
                       step_detail='', each=[],
                       progress=None, result=None, error=None,
                       start=time.time(), end=None,
                       instruments=list(data.instruments))
    threading.Thread(target=_run_quick, args=(qid, data,
                                              opts.get('target', ''),
                                              trend_order(opts)),
                     daemon=True).start()
    return quick_state(qid)


def forget() -> Dict[str, Any]:
    """
    A new target: the server forgets what the page did, as a new session
    would (its runs that ended, its quick looks, the copy of the NASA
    Exoplanet Archive read into memory, read again from disk when asked);
    a run still running is kept (a process of its own: stopping it is the
    page's Stop), and nothing on disk is touched

    :return: dict, the runs kept
    """
    from koloa import archive as karchive
    for jid in [jid for jid, job in JOBS.items()
                if job.returncode is not None]:
        del JOBS[jid]
    # a quick FIP that runs ends on its own, unseen
    QUICKS.clear()
    karchive._TABLES = None
    return dict(kept=[job.id for job in JOBS.values()])


def quick_state(qid: str) -> Dict[str, Any]:
    """what the page shows of a quick FIP"""
    job = QUICKS[qid]
    return dict(job, elapsed=(job['end'] or time.time()) - job['start'])


def _quicklook_figures(data, source, quick, xr, yr, pr, title, each=()):
    """the figures of the quick look: the velocities shown, the quick FIP
    with its peaks named, the folds at them"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from koloa import plotting as kplot
    kplot.set_style('paper')
    colour = {inst: kplot.INST_COLOURS[it % 8]
              for it, inst in enumerate(data.instruments)}
    marker = {inst: kplot.INST_MARKERS[it % 8]
              for it, inst in enumerate(data.instruments)}
    figs = []
    fig, ax = plt.subplots(figsize=(10, 3.6))
    for inst in data.instruments:
        sel = data.inst == inst
        ax.errorbar(data.time[sel], data.rv[sel] - np.median(data.rv[sel]),
                    data.err[sel], fmt=marker[inst], ms=3.5, lw=0.6,
                    color=colour[inst],
                    label=f'{inst} ({source.get(inst, "")}, {int(sel.sum())})')
    ax.axhline(0, color='0.6', lw=0.6, ls=':')
    if xr:
        ax.set_xlim(*xr)
    if yr:
        ax.set_ylim(*yr)
    ax.set_xlabel('BJD - 2400000')
    ax.set_ylabel('RV - median [m s$^{-1}$]')
    ax.legend(fontsize=7, ncol=3, frameon=False, loc='upper left')
    ax.set_title(f'{title}: the velocities shown', fontsize=10)
    fig.tight_layout()
    figs.append(('series', fig))
    if not quick:
        return figs
    fig, ax = plt.subplots(figsize=(10, 3.8))
    per = np.asarray(quick['period'])
    ax.plot(per, quick['alone'], color='0.6', lw=0.8,
            label='the period alone')
    # each instrument in its colour, under the joint FIP (in black: the
    #   colour of the first instrument is koloa's blue)
    valid = [item for item in each if not item.get('skipped')]
    for item in valid:
        ax.plot(item['period'], item['family'], lw=0.7, alpha=0.85,
                color=colour.get(item['name'], '0.5'),
                label=f'{item["name"]} alone (period or alias)')
    ax.plot(per, quick['family'], color=kplot.C['text'], lw=1.0,
            label='the period or any of its aliases')
    ax.axhline(2, color='0.3', lw=0.8, ls=':', label='FIP = 1 %')
    lo, hi = pr if pr else (per.min(), per.max())
    for name, val in quick['window'].items():
        if lo <= val <= hi:
            ax.axvline(val, color='0.5', lw=0.8, ls=':')
            ax.text(val, 0.02, f' {name}', fontsize=7, color='0.4',
                    transform=ax.get_xaxis_transform(), va='bottom')
    for pl in quick['known']:
        if lo <= pl['P'] <= hi:
            ax.axvline(pl['P'], color=kplot.C['outlier'], lw=0.9, ls='--')
            ax.text(pl['P'], 0.80, f' {pl["name"]}', fontsize=7,
                    transform=ax.get_xaxis_transform(), va='top',
                    color=kplot.C['outlier'])
    top = max(max(quick['family']), 2.2)
    for pk in quick['peak_list']:
        if pk['named'] and lo <= pk['period'] <= hi:
            height = -np.log10(max(pk['family'], 1e-15))
            ax.annotate(f'#{pk["id"]}', (pk['period'], height),
                        xytext=(0, 8), textcoords='offset points',
                        ha='center', fontsize=8, fontweight='bold',
                        arrowprops=dict(arrowstyle='-', lw=0.6))
    ax.set_ylim(0, 1.15 * top)
    ax.set_xscale('log')
    kplot.plain_log_ticks(ax, 'x')
    ax.set_xlim(lo, hi)
    ax.set_xlabel('period [d]')
    ax.set_ylabel('-log$_{10}$ FIP')
    ax.legend(fontsize=7, frameon=False, loc='upper right')
    sett = quick['settings']
    ax.set_title(f'quick FIP: {quick["n"]} nights; no GP, {sett["kmax"]} '
                 f'signals, {sett["nsweep"]} sweeps, {quick["passes"]} '
                 f'pass(es)', fontsize=9)
    fig.tight_layout()
    figs.append(('fip', fig))
    if valid:
        fig, axes = plt.subplots(len(valid), 1, sharex=True, squeeze=False,
                                 figsize=(10, 0.6 + 1.5 * len(valid)))
        for ax, item in zip(axes[:, 0], valid):
            ax.plot(item['period'], item['alone'], color='0.6', lw=0.7)
            ax.plot(item['period'], item['family'], lw=0.9,
                    color=colour.get(item['name'], 'k'))
            ax.axhline(2, color='0.3', lw=0.7, ls=':')
            for pl in quick['known']:
                if lo <= pl['P'] <= hi:
                    ax.axvline(pl['P'], color=kplot.C['outlier'], lw=0.8,
                               ls='--')
            best = item['peaks'][0] if item['peaks'] else None
            ax.text(0.005, 0.92, f'{item["name"]}, {item["n"]} nights'
                    + (f'; best {best["period"]:.4f} d, FIP '
                       f'{best["family"]:.1e}' if best else ''),
                    transform=ax.transAxes, fontsize=8, va='top')
            # room above the curve for the label of the panel
            ax.set_ylim(0, max(2.4, 1.4 * max(item['family'])))
            ax.set_ylabel('-log$_{10}$ FIP', fontsize=7)
        ax = axes[-1, 0]
        ax.set_xscale('log')
        kplot.plain_log_ticks(ax, 'x')
        ax.set_xlim(lo, hi)
        ax.set_xlabel('period [d]')
        fig.tight_layout(h_pad=0.3)
        figs.append(('each', fig))
    folds = quick.get('folds') or []
    if folds:
        ncol = min(3, len(folds))
        nrow = int(np.ceil(len(folds) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(10, 3.1 * nrow),
                                 squeeze=False)
        for ax in axes.ravel()[len(folds):]:
            ax.set_visible(False)
        for ax, item in zip(axes.ravel(), folds):
            for inst in item['instruments']:
                ax.errorbar(inst['phase'], inst['rv'], inst['err'],
                            fmt=marker.get(inst['name'], 'o'), ms=3,
                            lw=0.5, color=colour.get(inst['name'], 'k'))
                # circled: less than an even chance of being valid
                low = np.asarray(inst.get('valid', []), float) < 0.5
                if low.any():
                    ax.plot(np.asarray(inst['phase'])[low],
                            np.asarray(inst['rv'])[low], 'o', mfc='none',
                            mec='k', ms=7, mew=0.8, zorder=5)
            ax.plot(item['curve']['phase'], item['curve']['rv'], color='k',
                    lw=1.0)
            ax.set_title(f'#{item["id"]}: P = {item["period"]:.4f} d, '
                         f'K = {item["K"]:.2f} $\\pm$ {item["K_err"]:.2f} '
                         f'm/s', fontsize=8.5)
            ax.set_xlabel('phase (0 = conjunction)')
            ax.set_ylabel('RV [m s$^{-1}$]')
            # the velocities shown as the series is
            if yr:
                ax.set_ylim(*yr)
        fig.tight_layout()
        figs.append(('folds', fig))
    return figs


def _quicklook_tex(data, source, quick, opts, xr, yr, pr, figs,
                   command_line: str, each=()) -> str:
    """the quick look in LaTeX: its figures, then a page of numbers"""
    import platform
    import textwrap
    import koloa
    from koloa.latex import _path, _preamble, escape
    title = (opts.get('target') or '').strip() or 'the series'
    out = [_preamble(title).replace('detailed analysis of', 'quick look at'),
           f'\\section*{{koloa: a quick look at {escape(title)}}}',
           f'{time.strftime("%Y-%m-%d %H:%M")}, before a detailed report: '
           f'the velocities shown, the quick FIP of what is shown (no GP), '
           f'and the folds at its strongest peaks.\n']
    captions = dict(
        series='The velocities shown, each instrument about its median, in '
               'the ranges of the page.',
        fip='The quick FIP of the series shown: of the period or any of its '
            'aliases (black, what decides on a planet) and of the period '
            'alone (grey); FIP = 1 \\% dotted, the window (a day, a synodic '
            'month, a year) dotted, the known planets dashed; the peaks '
            'below a FIP of 10 \\% (or the best three) numbered.',
        folds='The nightly means folded at each numbered peak: a sinusoid '
              'fitted with an offset per instrument and the trend, phase 0 '
              'at the conjunction (each night weighted by its probability '
              'of being valid), on the velocity range of the series; '
              'circled, a night with less than a 50 \\% probability of '
              'being valid (an outlier, as the FIP saw it).',
        each='The quick FIP of each instrument on its own (its nightly '
             'means, its jitter sampled in the FIP): of the period or any '
             'of its aliases (colour) and of the period alone (grey); '
             'FIP = 1 \\% dotted, the known planets dashed.')
    for name in figs:
        out.append('\\begin{figure}[H]\\centering\n'
                   f'\\includegraphics[width=\\linewidth,height=0.27'
                   f'\\textheight,keepaspectratio]{{{name}.pdf}}\n'
                   f'\\caption{{{captions[name]}}}\n\\end{{figure}}')
    out.append('\\newpage\n\\section*{The numbers}')
    # the star
    star = {}
    if (opts.get('target') or '').strip():
        try:
            from koloa.archive import resolve
            star = resolve(opts['target'])
        except Exception:
            star = {}
    out.append('\\subsection*{The star}\n\\begin{tabular}{@{}ll@{}}')
    for key, val in (('Name given', opts.get('target') or '--'),
                     ('SIMBAD', star.get('main')), ('TIC', star.get('tic')),
                     ('Gaia DR3', star.get('gaia_dr3'))):
        out.append(f'{key} & {escape(val or "--")} \\\\')
    out.append('\\end{tabular}\n')
    # the data, instrument by instrument
    paths = {os.path.basename(str(item.get('path'))): str(item.get('path'))
             for item in (opts.get('files') or []) if item.get('path')}
    nights = data.nightly()
    infl = (quick or {}).get('inflation', {})
    out.append('\\subsection*{The data}\n{\\small\\begin{tabularx}'
               '{\\linewidth}{@{}lLrrrrrrr@{}}\n\\toprule\n'
               'Instrument & Source & N & Nights & First & Last & '
               '$\\sigma$ & rms & Added \\\\\n'
               '& & & & [BJD$-$2400000] & & [m/s] & [m/s] & [m/s] '
               '\\\\\n\\midrule')
    for inst in data.instruments:
        sel = data.inst == inst
        src = source.get(inst, '')
        where = (_path(paths[src[6:]]) if src.startswith('file: ')
                 and src[6:] in paths else escape(src))
        rms = float(np.std(data.rv[sel]))
        added = infl.get(str(inst))
        out.append(f'{escape(inst)} & {where} & {int(sel.sum())} & '
                   f'{int(np.sum(nights.inst == inst))} & '
                   f'{data.time[sel].min():.2f} & {data.time[sel].max():.2f}'
                   f' & {np.median(data.err[sel]):.2f} & {rms:.2f} & '
                   + (f'{added:.2f}' if added is not None else '--')
                   + ' \\\\')
    out.append('\\bottomrule\n\\end{tabularx}}\n')
    out.append(f'{data.n} exposures, {nights.n} nights over '
               f'{data.baseline:.0f}\\,d. Added: the error added in '
               f'quadrature to each instrument for the FIP (its second '
               f'pass), to the noise its errors do not hold (one instrument '
               f'is the reference, its jitter sampled in the FIP).\n')
    # what is shown
    asked = [name for key, name in (('dace', 'DACE'),
                                    ('carmenes', 'CARMENES DR1'))
             if opts.get(key)]
    left = str(opts.get('exclude') or '').split()
    rows = [('Archives', ', '.join(asked) or 'none'),
            ('Left out', ', '.join(left) or 'none')]
    if xr:
        rows.append(('Time shown', f'{xr[0]:.1f} to {xr[1]:.1f}'))
    if yr:
        rows.append(('Velocities shown', f'{yr[0]:.1f} to {yr[1]:.1f} m/s'))
    if pr:
        rows.append(('Periods shown', f'{pr[0]:.3g} to {pr[1]:.4g} d'))
    out.append('\\subsection*{What is shown}\n\\begin{tabular}{@{}ll@{}}')
    out += [f'{key} & {escape(val)} \\\\' for key, val in rows]
    out.append('\\end{tabular}\n')
    if quick:
        sett = quick['settings']
        pkk = ', '.join(f'P(k={it}) = {val:.2f}'
                        for it, val in enumerate(quick.get('pk', [])))
        order = sett.get('trend', 1)
        drift = ('no trend' if order < 1 else 'a trend (the acceleration of '
                 'the star)' if order == 1 else 'a trend and its curvature')
        out.append('\\subsection*{The quick FIP}\n'
                   f'Outlier-aware, nightly means, no GP, {drift}; '
                   f'{sett["kmax"]} signals at most, {sett["nsweep"]} sweeps '
                   f'after {sett["nburn"]} (two chains), {quick["passes"]} '
                   f'pass(es)'
                   + (f' (the second with the noise of a fit with '
                      + ', '.join(f'{per:.4f}' for per in quick["planets"])
                      + '\\,d)' if quick['passes'] > 1 else '')
                   + f'. {pkk}.\n')
        acc = quick.get('acceleration')
        if acc:
            for key, what, unit in (('accel', 'Acceleration of the star, '
                                     '$dv/dt$', 'm\\,s$^{-1}$\\,yr$^{-1}$'),
                                    ('jerk', 'Its change, $d^2v/dt^2$',
                                     'm\\,s$^{-1}$\\,yr$^{-2}$')):
                if key in acc:
                    val, low, high = acc[key]
                    out.append(f'{what} $= {val:+.3g}_{{-{low:.2g}}}'
                               f'^{{+{high:.2g}}}$ {unit} '
                               f'({acc[key + "_sigma"]:.1f}\\,$\\sigma$)'
                               + (f', at BJD$-$2400000 = {acc["tref"]:.1f}'
                                  if key == 'accel' else '') + '.\n')
            out.append('From the fit of the quick look (its signals and '
                       'the known planets, no GP; Laplace errors), the '
                       'perspective acceleration included.\n')
        folds = {item['id']: item for item in quick.get('folds') or []}
        out.append('{\\small\\begin{tabular}{@{}rrrrrrr@{}}\n\\toprule\n'
                   'ID & P [d] & FIP (P or alias) & FIP (P alone) & '
                   'K [m/s] & rms [m/s] & conjunction \\\\\n\\midrule')
        for pk in quick['peak_list']:
            item = folds.get(pk['id'])
            kk = (f'${item["K"]:.2f} \\pm {item["K_err"]:.2f}$'
                  if item else '--')
            out.append(f'\\#{pk["id"]} & {pk["period"]:.4f} & '
                       f'{pk["family"]:.1e} & {pk["alone"]:.1e} & {kk} & '
                       + (f'{item["rms"]:.2f} & {item["tc"]:.4f}' if item
                          else '-- & --') + ' \\\\')
        out.append('\\bottomrule\n\\end{tabular}}\n')
        if quick['known']:
            from koloa.aliases import same_family
            width = 1.0 / max(data.baseline, 1.0)
            out.append('Known planets (NASA Exoplanet Archive): '
                       + '; '.join(
                           f'{escape(pl["name"])}, {pl["P"]:.5g}\\,d'
                           + next((f' (\\#{pk["id"]})'
                                   for pk in quick['peak_list']
                                   if same_family(pl['P'], pk['period'],
                                                  width)), ' (no peak)')
                           for pl in quick['known']) + '.\n')
    if each:
        out.append('\\subsection*{The FIP of each instrument}\n'
                   '{\\small\\begin{tabular}{@{}lrl@{}}\n\\toprule\n'
                   'Instrument & Nights & Its best peaks: P [d] (FIP, P or '
                   'alias) \\\\\n\\midrule')
        for item in each:
            what = ('fewer than ' + str(QUICK_MIN_NIGHTS) + ' nights: no '
                    'FIP of its own' if item.get('skipped') else
                    '; '.join(f'{pk["period"]:.4f} ({pk["family"]:.1e})'
                              for pk in item['peaks']) or 'none')
            out.append(f'{escape(item["name"])} & {item["n"]} & '
                       f'{escape(what)} \\\\')
        out.append('\\bottomrule\n\\end{tabular}}\n')
    if command_line:
        out.append('\\subsection*{The detailed report}\n'
                   'The command line of the detailed report, as the page '
                   'set it:\n\\begin{verbatim}\n'
                   + ' \\\n    '.join(textwrap.wrap(
                       command_line, 88, break_long_words=True,
                       break_on_hyphens=False))
                   + '\n\\end{verbatim}')
    out.append('\\vfill{\\small\\color{muted}koloa '
               f'{escape(getattr(koloa, "__version__", ""))}, Python '
               f'{escape(platform.python_version())}, '
               f'{escape(platform.node())}, '
               f'{time.strftime("%Y-%m-%d %H:%M")}.}}\n')
    out.append('\\end{document}\n')
    return '\n'.join(out)


def quicklook_pdf(opts: Dict[str, Any], xr=None, yr=None, pr=None,
                  qid: str = '', command_line: str = '') -> bytes:
    """
    The quick look as a PDF, a LaTeX document: the velocities shown (the
    ranges of the page), the quick FIP with its peaks named, the folds at
    them, then a page of numbers (the star, every instrument and its file,
    what is shown, the FIP, its peaks and the known planets, the command
    line of the report); the figures alone (matplotlib) where there is no
    pdflatex

    :return: bytes, the PDF
    """
    import io
    import tempfile
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from koloa.latex import compile_pdf
    data, source, _ = selection(opts)
    if data is None:
        raise ValueError('no velocities to draw')
    job = QUICKS.get(qid) or {}
    quick, each = job.get('result'), job.get('each') or []
    title = (opts.get('target') or '').strip() or 'the series'
    figs = _quicklook_figures(data, source, quick, xr, yr, pr, title, each)
    tmp = tempfile.mkdtemp(prefix='koloa_quicklook_')
    try:
        for name, fig in figs:
            fig.savefig(os.path.join(tmp, f'{name}.pdf'))
        tex = os.path.join(tmp, 'quicklook.tex')
        with open(tex, 'w') as handle:
            handle.write(_quicklook_tex(data, source, quick, opts, xr, yr, pr,
                                        [name for name, _ in figs],
                                        command_line, each))
        pdf = compile_pdf(tex)
        if pdf and os.path.exists(pdf):
            with open(pdf, 'rb') as handle:
                return handle.read()
        # no pdflatex: the figures, one per page
        buf = io.BytesIO()
        with PdfPages(buf) as book:
            for _, fig in figs:
                book.savefig(fig)
        return buf.getvalue()
    finally:
        for _, fig in figs:
            plt.close(fig)
        shutil.rmtree(tmp, ignore_errors=True)


# =============================================================================
# The results remembered: a quick look kept, to be recalled as it was
# =============================================================================
#: where the results remembered are kept, one folder each
REMEMBERED = os.path.join(os.path.expanduser('~'), '.cache', 'koloa',
                          'remembered')


def _sha1(path: str) -> str:
    """the SHA-1 of a file, to know whether it changed"""
    digest = hashlib.sha1()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _entry_folder(rid: str) -> str:
    """the folder of a result remembered, its name checked (no path)"""
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9+\-._]*', rid or ''):
        raise ValueError(f'not a result remembered: {rid!r}')
    folder = os.path.join(REMEMBERED, rid)
    if not os.path.exists(os.path.join(folder, 'entry.json')):
        raise ValueError(f'no result remembered as {rid}')
    return folder


def _summary(page: Dict[str, Any], result: Dict[str, Any],
             each: List[Dict[str, Any]]) -> Dict[str, Any]:
    """what the list of the results remembered shows of one"""
    named = [pk for pk in result.get('peak_list', []) if pk.get('named')]
    detailed = page.get('detailed') or {}
    return dict(
        n=result.get('n'), instruments=result.get('instruments', {}),
        peaks=[dict(id=pk['id'], period=pk['period'], family=pk['family'])
               for pk in named],
        known=result.get('known', []),
        each=[dict(name=one['name'], n=one['n'], skipped=one['skipped'],
                   best=(one['peaks'][0] if one.get('peaks') else None))
              for one in each],
        files=[os.path.basename(row['path'])
               for row in page.get('files') or []],
        archives=[name for key, name in (('dace', 'DACE'),
                                         ('carmenes', 'CARMENES DR1'))
                  if detailed.get(key)],
        exclude=str(detailed.get('exclude') or ''),
        acceleration=result.get('acceleration'))


def remember(page: Dict[str, Any], qid: str, note: str = ''
             ) -> Dict[str, Any]:
    """
    A quick look remembered: the page (its fields, its ranges), its quick
    FIP (joint and of each instrument), the velocities plotted, and a copy
    of its files and of the velocities of the star's archives, so that it
    is recalled as it was, the files moved or the archives refreshed since

    :param page: dict, the state of the page (target, files, root, outdir,
                 detailed: the options of the report, clip, view, periods)
    :param qid: str, the quick FIP shown (ended)
    :param note: str, a word of why it is worth keeping

    :return: dict, the entry (its id, and its summary)
    """
    from koloa.gather import folder_name
    job = QUICKS.get(qid)
    if not job or job.get('status') != 'done' or not job.get('result'):
        raise ValueError('nothing to remember yet: the quick FIP has not '
                         'ended')
    target = str(page.get('target') or '').strip()
    rid = (f'{folder_name(target) or "series"}_'
           f'{time.strftime("%Y%m%d-%H%M%S")}')
    folder = os.path.join(REMEMBERED, rid)
    os.makedirs(os.path.join(folder, 'files'))
    try:
        return _remember_into(folder, rid, target, page, job, qid, note)
    except BaseException:
        # nothing half written: the list would not show it anyway
        shutil.rmtree(folder, ignore_errors=True)
        raise


def _remember_into(folder: str, rid: str, target: str,
                   page: Dict[str, Any], job: Dict[str, Any], qid: str,
                   note: str) -> Dict[str, Any]:
    """the copies and the files of a result remembered (remember)"""
    from koloa.gather import folder_name
    files = []
    for rank, row in enumerate(page.get('files') or []):
        path = os.path.abspath(os.path.expanduser(row['path']))
        copy = os.path.join(folder, 'files',
                            f'{rank}_{os.path.basename(path)}')
        shutil.copy2(path, copy)
        files.append(dict(path=path, label=row.get('label') or '',
                          copy=copy, sha1=_sha1(copy)))
    # the velocities of the archives (not the photometry: the report
    #   fetches its own)
    root = str(page.get('root') or 'archives')
    arch = None
    src = os.path.join(root, folder_name(target)) if target else ''
    if src and os.path.isdir(os.path.join(src, 'rv')):
        dst = os.path.join(folder, 'archives', folder_name(target))
        shutil.copytree(os.path.join(src, 'rv'), os.path.join(dst, 'rv'))
        for name in ('target.json', 'manifest.json'):
            if os.path.exists(os.path.join(src, name)):
                shutil.copy2(os.path.join(src, name), os.path.join(dst, name))
        allrv = os.path.join(dst, 'rv', 'all_rv.csv')
        arch = dict(root=os.path.abspath(root),
                    copy=os.path.join(folder, 'archives'),
                    sha1=_sha1(allrv) if os.path.exists(allrv) else '')
    detailed = page.get('detailed') or {}
    rv = velocities([dict(path=row['copy'], label=row['label'])
                     for row in files], target,
                    arch['copy'] if arch else root,
                    dace=bool(detailed.get('dace')),
                    carmenes=bool(detailed.get('carmenes')))
    # the files by their own names, not their copies'
    names = {os.path.basename(row['copy']): os.path.basename(row['path'])
             for row in files}
    for inst in rv.get('instruments', []):
        src = inst.get('source', '')
        if src.startswith('file: ') and src[6:] in names:
            inst['source'] = 'file: ' + names[src[6:]]
    rv['notes'] = [next((note.replace(copy, orig, 1)
                         for copy, orig in names.items()
                         if note.startswith(copy + ':')), note)
                   for note in rv.get('notes', [])]
    if arch:
        rv['notes'] = [note.replace(arch['copy'], root)
                       for note in rv['notes']]
    entry = dict(id=rid, target=target, note=str(note or ''),
                 created=time.strftime('%Y-%m-%d %H:%M'),
                 page=page, files=files, archives=arch,
                 summary=_summary(page, job['result'], job.get('each') or []))
    for name, value in (('entry', entry), ('rv', rv),
                        ('quick', dict(result=job['result'],
                                       each=job.get('each') or [],
                                       elapsed=quick_state(qid)['elapsed']))):
        with open(os.path.join(folder, f'{name}.json'), 'w') as handle:
            json.dump(value, handle)
    return dict(id=rid, created=entry['created'], summary=entry['summary'])


def remembered() -> List[Dict[str, Any]]:
    """the results remembered, the last first (their entries)"""
    out = []
    if os.path.isdir(REMEMBERED):
        for rid in os.listdir(REMEMBERED):
            path = os.path.join(REMEMBERED, rid, 'entry.json')
            if not os.path.exists(path):
                continue
            try:
                with open(path) as handle:
                    entry = json.load(handle)
            except ValueError:   # a folder half written
                continue
            out.append({key: entry.get(key) for key in
                        ('id', 'target', 'note', 'created', 'summary')})
    return sorted(out, key=lambda entry: entry['id'].rsplit('_', 1)[-1],
                  reverse=True)


def recall(rid: str) -> Dict[str, Any]:
    """
    A result remembered, as it was: the page, the velocities plotted, the
    quick FIP (ready for its PDF); each file is the original when it has
    not changed, else its copy, and the archives likewise

    :return: dict, page, rv, quick (its state), entry, notes
    """
    folder = _entry_folder(rid)
    loaded = {}
    for name in ('entry', 'rv', 'quick'):
        with open(os.path.join(folder, f'{name}.json')) as handle:
            loaded[name] = json.load(handle)
    entry = loaded['entry']
    page, notes = dict(entry['page']), []
    rows = []
    for row in entry['files']:
        same = (os.path.exists(row['path'])
                and _sha1(row['path']) == row['sha1'])
        rows.append(dict(path=row['path'] if same else row['copy'],
                         label=row['label']))
        if not same:
            notes.append(f'{os.path.basename(row["path"])}: changed or '
                         f'moved since, its copy used')
    page['files'] = rows
    arch = entry.get('archives')
    if arch:
        from koloa.gather import folder_name
        now = os.path.join(arch['root'], folder_name(entry['target']), 'rv',
                           'all_rv.csv')
        same = os.path.exists(now) and _sha1(now) == arch['sha1']
        page['root'] = arch['root'] if same else arch['copy']
        if not same:
            notes.append('the archives refreshed or moved since: their '
                         'copy used')
    qid = uuid.uuid4().hex[:8]
    quick = loaded['quick']
    QUICKS[qid] = dict(id=qid, status='done', step='done', step_detail='',
                       each=quick['each'], progress=None,
                       result=quick['result'], error=None, start=0.0,
                       end=float(quick.get('elapsed') or 0.0),
                       instruments=list(quick['result']['instruments']))
    return dict(page=page, rv=loaded['rv'], quick=quick_state(qid),
                entry={key: entry.get(key) for key in
                       ('id', 'target', 'note', 'created', 'summary')},
                notes=notes)


def unremember(rid: str) -> Dict[str, Any]:
    """a result remembered, forgotten (its folder and its copies)"""
    folder = _entry_folder(rid)
    shutil.rmtree(folder)
    return dict(id=rid)


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
            if url.path == '/api/remembered':
                return self._json(remembered())
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
                                    body.get('quick', ''),
                                    body.get('command', ''))
                self.send_response(200)
                self.send_header('Content-Type', 'application/pdf')
                self.send_header('Content-Length', str(len(pdf)))
                self.send_header('Content-Disposition',
                                 'attachment; filename="koloa_quicklook.pdf"')
                self.end_headers()
                self.wfile.write(pdf)
                return None
            if path == '/api/forget':
                return self._json(forget())
            if path == '/api/remember':
                return self._json(remember(body.get('page', {}),
                                           body.get('quick', ''),
                                           body.get('note', '')))
            if path == '/api/recall':
                return self._json(recall(body.get('id', '')))
            if path == '/api/unremember':
                return self._json(unremember(body.get('id', '')))
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

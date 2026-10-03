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
from typing import Any, Dict, List, Optional, Tuple

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
            'carmenes'] + ['--no-vizier'] * off['vizier'] + ['--no-tess'] * \
            off['tess']
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
            elif message.startswith('result: ') and self.steps:
                # what the step came to: shown once it is over too
                self.steps[-1]['result'] = message[8:][:400]
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
    known = {}
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
    # its spectral type and mass (the archive's, or rough from the type)
    try:
        from koloa.stars import stellar
        if 'sptype' not in out:
            out['sptype'] = resolve(out.get('main') or name).get('sptype')
        out['star'] = stellar(out, known.get('star'))
    except Exception as err:  # a help, not a need
        out['star'] = dict(sptype=None, mass=None, mass_err=None,
                           source=None, error=str(err))
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


def _pick_files(start: str) -> Dict[str, Any]:
    """the dialog of this machine for several files at once (the batch
    FIP): paths, or cancelled"""
    if sys.platform == 'darwin':
        place = start.replace('\\', '\\\\').replace('"', '\\"')
        script = ('set chosen to choose file with prompt "Files of '
                  f'velocities (LBL .rdb, csv)" default location (POSIX file '
                  f'"{place}") with multiple selections allowed\n'
                  'set out to ""\n'
                  'repeat with one in chosen\n'
                  'set out to out & POSIX path of one & linefeed\n'
                  'end repeat\n'
                  'return out')
        proc = subprocess.run(['osascript', '-e', 'activate', '-e', script],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            if '-128' in proc.stderr:
                return dict(cancelled=True)
            raise RuntimeError(proc.stderr.strip() or 'osascript failed')
        paths = [line for line in proc.stdout.splitlines() if line.strip()]
    else:
        code = ('import sys, tkinter as tk\n'
                'from tkinter import filedialog\n'
                'root = tk.Tk(); root.withdraw()\n'
                "root.attributes('-topmost', True)\n"
                "print('\\n'.join(filedialog.askopenfilenames("
                "initialdir=sys.argv[1], filetypes=[('velocities', '*.rdb "
                "*.csv *.dat *.txt'), ('all', '*')])))")
        proc = subprocess.run([sys.executable, '-c', code, start],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError('no dialog on this machine (Tk): give a '
                               'folder')
        paths = [line for line in proc.stdout.splitlines() if line.strip()]
    return dict(paths=paths) if paths else dict(cancelled=True)


def list_files(folder: str, pattern: str = '*.rdb') -> Dict[str, Any]:
    """the files of a folder that match a pattern (the batch FIP)"""
    import glob
    folder = os.path.expanduser(folder or '.')
    if not os.path.isdir(folder):
        raise ValueError(f'no folder {folder}')
    found = sorted(path for path in glob.glob(os.path.join(
        folder, pattern or '*.rdb')) if os.path.isfile(path))
    return dict(paths=[os.path.abspath(path) for path in found])


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
    if kind == 'files':
        return _pick_files(start)
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
              dace: bool = False, carmenes: bool = False,
              vizier: bool = False):
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
            # the velocities only (the photometry is not plotted), with the
            #   BERV of each (the raw files of DACE and CARMENES)
            from koloa.gather import archive_berv
            gathered = load(folder, photometry=False)['rv']
            gathered.meta['BERV'] = archive_berv(folder, gathered.time)
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
        # the velocities published on VizieR (koloa.published), their
        #   spectra already in the series left out
        pdir = os.path.join(folder, 'rv', 'published')
        from koloa import published as kpub
        pub = kpub.load(pdir) if os.path.isdir(pdir) else None
        if pub is not None and not vizier:
            notes.append(f'VizieR: {pub.n} points gathered, not ticked')
        elif pub is not None:
            # the spectra the series has left out, and then an instrument
            #   with too few velocities for an offset of its own
            fresh = kpub.enough(kpub.new_spectra(pub, series))
            if fresh is not None:
                series.append(fresh)
                source.update({name: 'VizieR' for name in fresh.instruments})
                notes.append(f'VizieR ({pdir}): {fresh.n} points'
                             + (f' ({pub.n - fresh.n} the same spectra as '
                                f'the others)' if fresh.n < pub.n else ''))
            else:
                notes.append('VizieR: only spectra the others have')
    if not series:
        return None, source, notes
    return (series[0] if len(series) == 1 else merge(series)), source, notes


def velocities(files: Any = '', target: str = '', root: str = '',
               dace: bool = False, carmenes: bool = False,
               vizier: bool = False) -> Dict[str, Any]:
    """
    The velocities of a file and of a star's gathered archives, by
    instrument (each with its median taken out), for the plot of the page
    """
    data, source, notes = series_of(files, target, root, dace, carmenes,
                                    vizier)
    if data is None:
        return dict(instruments=[], notes=notes)
    out = []
    berv = _berv(data)
    for name in data.instruments:
        sel = data.inst == name
        # each instrument about its own median, whatever came before
        rv = data.rv[sel] - np.median(data.rv[sel])
        out.append(dict(name=name, n=int(sel.sum()),
                        source=source.get(name, ''),
                        median=float(np.median(data.rv[sel])),
                        time=np.round(data.time[sel], 6).tolist(),
                        rv=np.round(rv, 3).tolist(),
                        err=np.round(data.err[sel], 3).tolist(),
                        berv=(_listed(berv[sel], 3) if berv is not None
                              and np.any(np.isfinite(berv[sel])) else None),
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
#: what each quick FIP ran on, for the folds asked later: its nightly means,
#: the probability of each night to be valid, the order of its trend, the
#: FIP itself
_QUICK_DATA: Dict[str, Dict[str, Any]] = {}


def selection(opts: Dict[str, Any]):
    """the series the page shows, its instruments left out taken out"""
    data, source, notes = series_of(
        opts.get('files') or opts.get('file') or '', opts.get('target', ''),
        opts.get('root', ''), bool(opts.get('dace')),
        bool(opts.get('carmenes')), bool(opts.get('vizier')))
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
    def number(value):
        return float(value) if value is not None else None
    try:
        host = host_name(resolve(target))
        return [dict(name=pl['name'], P=float(pl['P']), K=number(pl.get('K')),
                     e=number(pl.get('e')), omega=number(pl.get('omega')),
                     tp=number(pl.get('tp')), tc=number(pl.get('tc')),
                     P_err=number(pl.get('P_err')),
                     tc_err=number(pl.get('tc_err')),
                     reference=pl.get('reference') or '')
                for pl in (known_planets(host=host)['planets'] if host
                           else []) if pl.get('P')]
    except Exception:  # a help, not a need
        return []


def transits_of(target: str, known: List[Dict[str, Any]]
                ) -> List[Dict[str, Any]]:
    """
    The transit ephemerides of a star: of its known planets that have a
    time of transit (the archive), and of its TOIs (TESS; not the false
    positives)

    :return: list of dict: name, source, P, P_err, tc, tc_err [BJD -
             2400000], reference
    """
    out = [dict(name=pl['name'], source='archive', P=pl['P'],
                P_err=pl.get('P_err'), tc=pl['tc'], tc_err=pl.get('tc_err'),
                reference=pl.get('reference', ''))
           for pl in known if pl.get('tc') is not None]
    try:
        from koloa.archive import TOI_NOT_PLANETS, resolve, tois
        tic = resolve(target).get('tic') if target.strip() else None
        for item in tois(tic) if tic else []:
            if item['tc'] is None or item['disposition'] in TOI_NOT_PLANETS:
                continue
            out.append(dict(name=f'TOI-{item["toi"]}', source='TESS',
                            P=item['P'], P_err=item['P_err'], tc=item['tc'],
                            tc_err=item['tc_err'],
                            reference=f'TESS ({item["disposition"]})'))
    except Exception:  # a help, not a need
        pass
    return out


def fold_transit(data, transit: Dict[str, Any],
                 valid: Optional[np.ndarray] = None, trend: int = 1
                 ) -> Dict[str, Any]:
    """
    The series folded on a transit ephemeris: phase 0 at the transit (not
    at a conjunction the velocities put), a circular orbit with its phase
    fixed by the transits (its K, each instrument's offset and the trend
    fitted, each night weighted by its probability of being valid); the
    ephemeris carried to the velocities with its error; and the
    conjunction the velocities put on their own, set against it

    :param data: RVData, the nightly means
    :param transit: dict, name, P, P_err, tc, tc_err (transits_of)
    :param valid: np.ndarray or None, each night's probability of being
                  valid
    :param trend: int, the order of the trend

    :return: dict, as fold() gives, with transit (the ephemeris at the
             velocities: t0, its error in days and in phase, the cycles
             since the transit; the conjunction of the velocities alone,
             its offset from t0 and in sigma)
    """
    per, tc = float(transit['P']), float(transit['tc'])
    sper = float(transit.get('P_err') or 0.0)
    stc = float(transit.get('tc_err') or 0.0)
    time_ = data.time
    # the transit nearest the middle of the velocities, and its error
    ncyc = int(np.round((np.median(time_) - tc) / per))
    t0 = tc + ncyc * per
    st0 = float(np.hypot(stc, ncyc * sper))
    insts = list(data.instruments)
    cols = [(data.inst == inst).astype(float) for inst in insts]
    tref = float(np.median(time_))
    cols += [((time_ - tref) / 365.25) ** order
             for order in range(1, trend + 1)]
    cols.append(-np.sin(2 * np.pi * (time_ - t0) / per))
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
    amp, kerr = float(coef[-1]), float(np.sqrt(cov[-1, -1]))
    shown = data.rv - design[:, :-1] @ coef[:-1]
    phase = ((time_ - t0) / per) % 1.0
    grid = np.linspace(0, 1, 201)
    curve = -amp * np.sin(2 * np.pi * grid)
    band = kerr * np.abs(np.sin(2 * np.pi * grid))
    # the conjunction of the velocities alone (a sinusoid, its phase free)
    free = fold(data, per, valid, trend)
    shift = (free['tc'] - t0 + 0.5 * per) % per - 0.5 * per
    sshift = float(np.hypot(free.get('tc_err') or 0.0, st0))
    out = dict(kind='sine', period=per, K=amp, K_err=kerr, tc=float(t0),
               rms=float(np.std(resid[good >= 0.5])), chi2=chi2,
               curve=dict(phase=grid.tolist(), rv=np.round(curve, 4).tolist(),
                          lo=np.round(curve - band, 4).tolist(),
                          hi=np.round(curve + band, 4).tolist()),
               instruments=[],
               model=dict(kind='sine', period=per, tref=float(t0),
                          offsets={str(inst): float(coef[it])
                                   for it, inst in enumerate(insts)},
                          trend=[float(val) for val in
                                 coef[len(insts):len(insts) + trend]],
                          ttrend=tref, tscale=365.25, cos=0.0,
                          sin=-amp),
               transit=dict(transit, t0=float(t0), t0_err=st0,
                            phase_err=st0 / per, cycles=ncyc,
                            conj=float(free['tc']),
                            conj_err=free.get('tc_err'),
                            shift=float(shift), shift_err=sshift,
                            shift_sigma=(abs(shift) / sshift if sshift > 0
                                         else None),
                            K_free=free['K'], K_free_err=free['K_err']))
    out['draws'] = []
    try:
        rng = np.random.default_rng(17)
        for draw in rng.multivariate_normal(coef, cov, NDRAW):
            out['draws'].append(dict(
                out['model'], offsets={str(inst): float(draw[it])
                                       for it, inst in enumerate(insts)},
                trend=[float(val) for val in
                       draw[len(insts):len(insts) + trend]],
                sin=-float(draw[-1])))
    except (ValueError, np.linalg.LinAlgError):
        out['draws'] = []
    berv = _berv(data)
    for inst in insts:
        sel = data.inst == inst
        out['instruments'].append(dict(
            name=str(inst), phase=np.round(phase[sel], 5).tolist(),
            rv=np.round(shown[sel], 3).tolist(),
            err=np.round(data.err[sel], 3).tolist(),
            time=np.round(time_[sel], 5).tolist(),
            berv=(_listed(berv[sel], 3) if berv is not None
                  and np.any(np.isfinite(berv[sel])) else None)))
        if valid is not None:
            out['instruments'][-1]['valid'] = np.round(valid[sel],
                                                       4).tolist()
    return out


def published(planet: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """the published orbit of a known planet as a Keplerian signal (the
    archive's default solution: P, K, e, omega [deg], and its phase from tc
    when there is one, the transits being far more precise than an RV tp,
    which the archive may take from another paper; else from tp), or None
    when it lacks K or an epoch"""
    from koloa import kepler
    if not planet.get('K') or (planet.get('tp') is None
                               and planet.get('tc') is None):
        return None
    ecc = float(planet.get('e') or 0.0)
    # an omega of 0 is one (not a missing one: 90 deg, a circular orbit's)
    omega = np.radians(float(planet['omega'] if planet.get('omega')
                             is not None else 90.0))
    tperi = (kepler.tc_to_tp(float(planet['tc']), planet['P'], ecc, omega)
             if planet.get('tc') is not None else float(planet['tp']))
    return dict(kind='kepler', period=float(planet['P']), tp=tperi, e=ecc,
                omega=float(omega), K=float(planet['K']),
                name=planet['name'], reference=planet.get('reference', ''))


def _attach_published(item: Dict[str, Any]) -> None:
    """a fold of a known planet given its published orbit: on the fold's
    phases (its curve) and, with the fold's offsets and trend, as a model
    of the series"""
    from koloa import kepler
    pub = published(item['known'])
    if pub is None:
        return
    for one in [item] + ([item['kepler']] if 'kepler' in item else []):
        grid = np.linspace(0, 1, 401)
        rv = kepler.rv_keplerian(one['tc'] + grid * one['period'],
                                 pub['period'], pub['tp'], pub['e'],
                                 pub['omega'], pub['K'])
        one['published'] = dict(pub, phase=grid.tolist(),
                                rv=np.round(rv, 4).tolist(),
                                model=dict(one['model'], **{
                                    key: pub[key] for key in
                                    ('kind', 'period', 'tp', 'e', 'omega',
                                     'K')}))


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


#: draws of a fold's solution kept for its 1-sigma envelope on the series
#: (from the full covariance of its fit: the orbit, the offsets and the
#: trend together)
NDRAW = 200


def envelope(best: np.ndarray, many: np.ndarray) -> Tuple[np.ndarray,
                                                            np.ndarray]:
    """
    The 1-sigma envelope of a model from draws of its parameters (their
    full covariance): the spread of the draws at each point, its 16th and
    84th percentiles about their median, set about the best model (a
    sharp, eccentric orbit's draws, their peaks at different times, have a
    median below its peak: the envelope is the spread, not that median)

    :param best: np.ndarray, the best model
    :param many: np.ndarray, (ndraw, len(best)), the model of each draw

    :return: tuple, the low and high sides of the envelope
    """
    low, mid, high = np.percentile(many, [15.87, 50.0, 84.13], axis=0)
    return best - (mid - low), best + (high - mid)
#: the names of the BERV column of a series (LBL: BERV; DACE: cal_berv)
BERV_NAMES = ('BERV', 'berv', 'cal_berv')


def _berv(data) -> Optional[np.ndarray]:
    """the BERV of each point of a series [km/s], or None"""
    for key in BERV_NAMES:
        val = (getattr(data, 'meta', None) or {}).get(key)
        if val is not None and np.asarray(val).dtype.kind in 'fiu':
            val = np.asarray(val, dtype=float)
            if np.any(np.isfinite(val)):
                return val
    return None


def _listed(values: np.ndarray, digits: int) -> List[Optional[float]]:
    """numbers for the page, nan as null (JSON has no nan)"""
    return [round(float(val), digits) if np.isfinite(val) else None
            for val in values]


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
             (about its offset and the trend, with their time, their BERV
             when the series has it, and their probability to be valid when
             given), the curve, and the model (its offsets, trend and
             sinusoid, to draw it on the series)
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
    # the error of the conjunction: of phi0 = atan2(sin, cos)
    gphi = np.array([-asin, acos]) / max(amp, 1e-12) ** 2
    tc_err = float(period / (2 * np.pi)
                   * np.sqrt(max(gphi @ cov[-2:, -2:] @ gphi, 0.0)))
    phase = ((time_ - tc) / period) % 1.0
    shown = data.rv - design[:, :-2] @ coef[:-2]
    grid = np.linspace(0, 1, 201)
    # the 1-sigma envelope of the curve: (cos, sin) of each phase through
    #   the covariance of the two amplitudes
    gvec = np.array([np.cos(2 * np.pi * (tc + grid * period - tref) / period),
                     np.sin(2 * np.pi * (tc + grid * period - tref) / period)])
    curve = -amp * np.sin(2 * np.pi * grid)
    sig = np.sqrt(np.einsum('ik,ij,jk->k', gvec, cov[-2:, -2:], gvec))
    out = dict(period=float(period), K=amp, K_err=kerr, tc=float(tc),
               tc_err=tc_err,
               rms=float(np.std(resid[good >= 0.5])), chi2=chi2,
               curve=dict(phase=grid.tolist(), rv=np.round(curve, 4).tolist(),
                          lo=np.round(curve - sig, 4).tolist(),
                          hi=np.round(curve + sig, 4).tolist()),
               instruments=[],
               # v(t) = offset + sum_k trend_k ((t - tref) / yr)^k
               #        + cos cos(2 pi (t - tref) / P) + sin sin(...)
               model=dict(kind='sine', period=float(period), tref=tref,
                          offsets={str(inst): float(coef[it])
                                   for it, inst in enumerate(insts)},
                          trend=[float(val) for val in
                                 coef[len(insts):len(insts) + trend]],
                          ttrend=tref, tscale=365.25,
                          cos=float(acos), sin=float(asin)))
    # draws of the solution (its covariance), for its envelope on the series
    out['draws'] = []
    try:
        rng = np.random.default_rng(13)
        for draw in rng.multivariate_normal(coef, cov, NDRAW):
            out['draws'].append(dict(
                out['model'], offsets={str(inst): float(draw[it])
                                       for it, inst in enumerate(insts)},
                trend=[float(val) for val in
                       draw[len(insts):len(insts) + trend]],
                cos=float(draw[-2]), sin=float(draw[-1])))
    except (ValueError, np.linalg.LinAlgError):
        out['draws'] = []
    berv = _berv(data)
    for inst in insts:
        sel = data.inst == inst
        out['instruments'].append(dict(
            name=str(inst), phase=np.round(phase[sel], 5).tolist(),
            rv=np.round(shown[sel], 3).tolist(),
            err=np.round(data.err[sel], 3).tolist(),
            time=np.round(time_[sel], 5).tolist(),
            berv=(_listed(berv[sel], 3) if berv is not None
                  and np.any(np.isfinite(berv[sel])) else None)))
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


def _levels(fit, data=None) -> Dict[str, Any]:
    """
    The zero of each instrument on the plot of the series, and the trend
    of the fit of the quick look: the series drawn about them shows one
    trend across the instruments, not each about its own median

    The fit weighs nights (its nightly means, each with the jitter of its
    instrument), the plot shows exposures: an instrument with a few
    discrepant nights and a few nights of many exposures has its fitted
    offset between the two, its exposures off the trend. Its zero is its
    fitted offset moved by the median of the residuals of its exposures to
    the fitted model (the trend and the signals): its exposures about the
    model.

    :param fit: FitResult, the fit of the quick look
    :param data: RVData or None, the exposures (the series shown)

    :return: dict, offsets (the zeros of the plot), fit_offsets (as
             fitted), trend_model (tref, tscale, coefs)
    """
    model = fit.model
    offsets = {str(inst): float(fit.theta[model.index[f'offset_{inst}']])
               for inst in model.data.instruments
               if f'offset_{inst}' in model.index}
    trend = [float(fit.theta[model.index[f'trend_{deg}']])
             for deg in range(1, model.trend + 1)]
    tscale = float(max(model.data.baseline, 1e-9))
    zeros = dict(offsets)
    if data is not None:
        span = (data.time - model.tref) / tscale
        rest = sum((val * span ** (deg + 1) for deg, val in enumerate(trend)),
                   np.zeros(data.n))
        for ip in range(len(model.planets)):
            rest = rest + model.planet_rv(fit.theta, ip, data.time)
        for inst in data.instruments:
            sel = data.inst == inst
            if str(inst) in offsets and np.any(sel):
                zeros[str(inst)] = offsets[str(inst)] + float(np.median(
                    data.rv[sel] - offsets[str(inst)] - rest[sel]))
    # the signals of the fit (its second pass: those found and the known
    #   planets), for the model drawn with the trend
    signals = []
    for ip in range(len(model.planets)):
        per, tperi, ecc, omega, amp = model.orbit(fit.theta, ip)
        signals.append(dict(kind='kepler', period=float(per), tp=float(tperi),
                            e=float(ecc), omega=float(omega), K=float(amp)))
    return dict(offsets=zeros, fit_offsets=offsets, trend_model=dict(
        tref=float(model.tref), tscale=tscale, coefs=trend, signals=signals))


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


def _go_on(job: Dict[str, Any]) -> None:
    """between the steps of a quick FIP: on, unless it was stopped"""
    if job.get('cancel'):
        from koloa.fip import FIPCancelled
        raise FIPCancelled('the quick FIP was stopped')


def stop_quick(qid: str) -> Dict[str, Any]:
    """a quick FIP stopped: while it waits for its turn, or as it runs
    (its chains ended at once, the step it is at left)"""
    job = QUICKS[qid]
    if job['status'] == 'running':
        job['cancel'] = True
    return quick_state(qid)


def _run_quick(qid: str, data, target: str, trend: int = 1,
               each: bool = True):
    """the quick FIP, in a thread, in the two passes of the report without
    its GP: the errors of each instrument inflated to the noise of a fit
    without planets, a first FIP; then to the noise of a fit with the
    signals it found and the known planets (their variance is not noise),
    the FIP again"""
    from koloa import fip as kfip
    from koloa.detailed import _fip
    from koloa.fit import RVModel
    job = QUICKS[qid]
    # its turn (one at a time), unless it is stopped while it waits
    while not _QUICK_LOCK.acquire(timeout=0.5):
        if job.get('cancel'):
            job['status'], job['end'] = 'stopped', time.time()
            return
    try:
        try:
            _go_on(job)
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
            # a stop ends the chains of the FIP running at once
            kfip.CANCEL_HOOK = lambda: bool(job.get('cancel'))
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
                _go_on(job)
                job['step'] = 'planets'
                fit = RVModel(nights, [dict(period=per, period_range=(
                    0.98 * per, 1.02 * per)) for per in pers],
                    likelihood='mixture', unit='both', trend=trend,
                    seq_jitter=seq_jitter).fit(nstart=2, quiet=True)
                _go_on(job)
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
            _QUICK_DATA[qid] = dict(nights=nights, valid=valid, trend=trend,
                                    res=res)
            folds = [dict(fold(nights, pk['period'], valid, trend),
                          id=pk['id'])
                     for pk in named]
            job['result'] = dict(
                _fip_curves(res), known=known, window=WINDOW, passes=passes,
                transits=transits_of(target, known),
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
                subtracted=job.get('subtracted') or [],
                acceleration=_quick_acceleration(fit if pers else noise,
                                                 trend),
                **_levels(fit if pers else noise, data))
            # the FIP of each instrument on its own (its jitter sampled in
            #   the FIP, no inflation needed), when there are several
            insts = list(nights.instruments)
            job['each'] = []
            if len(insts) > 1 and each:
                from koloa.fip import oafip
                from koloa.utils import blas_threads
                for rank, inst in enumerate(insts):
                    sub = nights.select(nights.inst == inst)
                    if sub.n < QUICK_MIN_NIGHTS:
                        job['each'].append(dict(name=str(inst), n=int(sub.n),
                                                skipped=True))
                        continue
                    _go_on(job)
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
        except kfip.FIPCancelled:
            job['status'] = 'stopped'
        except Exception as err:
            job['status'] = 'failed'
            job['error'] = f'{type(err).__name__}: {err}'
        finally:
            kfip.PROGRESS_HOOK = None
            kfip.CANCEL_HOOK = None
            job['end'] = time.time()
    finally:
        _QUICK_LOCK.release()


def subtracted(data, opts: Dict[str, Any]):
    """the series without the signals the page subtracted (the solution
    of a fold each: its sinusoid or Keplerian orbit, not its offsets and
    trend, which the FIP fits again)"""
    for mod in opts.get('subtract') or []:
        data.rv = data.rv - signal_at(mod, data.time)
    return data


def _quick_series(qid: str, opts: Dict[str, Any]) -> Dict[str, Any]:
    """what a quick FIP ran on (or, for one recalled, the series of the
    page again, without the probabilities of its nights)"""
    if qid not in _QUICK_DATA:
        data, _, _ = selection(opts)
        if data is None:
            raise ValueError('no velocities to fold')
        _QUICK_DATA[qid] = dict(nights=subtracted(data, opts).nightly(),
                                valid=None, trend=trend_order(opts), res=None)
    return _QUICK_DATA[qid]


def _snap(res, period: float, baseline: float) -> float:
    """the period of the deepest dip of the FIP of the period alone near
    a period clicked, within a peak's width in frequency (1 / the
    baseline) each side: a click anywhere on a peak folds at its top"""
    freq = np.asarray(res.freq)
    near = np.abs(freq - 1.0 / period) <= 1.0 / max(baseline, 1.0)
    if not np.any(near):
        return period
    idx = np.flatnonzero(near)
    return float(1.0 / freq[idx[np.argmin(np.asarray(res.fip)[idx])]])


def fold_request(qid: str, opts: Dict[str, Any], period: Optional[float] = None,
                 fid: Optional[int] = None, kind: str = 'sine',
                 snap: bool = False, known: str = '',
                 transit: str = '') -> Dict[str, Any]:
    """
    A fold the page asks for: at a period of its own (clicked on the FIP,
    then moved to the dip nearest; or typed, as it is), or the Keplerian
    orbit of a fold already there; kept with the quick FIP (its PDF has it)

    :return: dict, fold (with its Keplerian, when asked)
    """
    job = QUICKS.get(qid)
    if not job or not job.get('result'):
        raise ValueError('no quick FIP to fold')
    series = _quick_series(qid, opts)
    nights = series['nights']
    folds = job['result'].setdefault('folds', [])
    if transit:
        # a transit ephemeris: phase 0 at its transit, its period
        eph = next((item for item in job['result'].get('transits') or []
                    if item['name'] == transit), None)
        if eph is None:
            raise ValueError(f'no transit ephemeris {transit}')
        series = _quick_series(qid, opts)
        item = next((item for item in folds
                     if (item.get('transit') or {}).get('name') == transit),
                    None)
        if item is None:
            item = dict(fold_transit(series['nights'], eph, series['valid'],
                                     series['trend']),
                        id=max([item['id'] for item in folds] + [0]) + 1,
                        forced=True)
            folds.append(item)
        if kind == 'kepler' and 'kepler' not in item:
            item['kepler'] = fold_kepler(series['nights'], eph['P'],
                                         series['trend'],
                                         transit=item['transit'])
            item['kepler']['transit'] = item['transit']
        return dict(fold=item)
    planet = None
    if known:
        # a known planet: at its published period, whatever the FIP says
        planet = next((pl for pl in job['result'].get('known') or []
                       if pl['name'] == known), None)
        if planet is None:
            raise ValueError(f'no known planet {known}')
        period, snap = planet['P'], False
    if fid is not None:
        item = next((item for item in folds if item['id'] == int(fid)), None)
        if item is None:
            raise ValueError(f'no fold #{fid}')
    else:
        period = float(period)
        if not 0.05 < period < 1e6:
            raise ValueError(f'not a period to fold at: {period}')
        if snap and series['res'] is not None:
            period = _snap(series['res'], period, nights.baseline)
        # a period already folded: that fold
        item = next((item for item in folds
                     if abs(item['period'] / period - 1) < 1e-7), None)
        if item is None:
            item = dict(fold(nights, period, series['valid'],
                             series['trend']),
                        id=max([item['id'] for item in folds] + [0]) + 1,
                        forced=True)
            folds.append(item)
    if planet is not None:
        item['known'] = planet
    if kind == 'kepler' and 'kepler' not in item:
        item['kepler'] = fold_kepler(nights, item['period'], series['trend'])
    if item.get('known'):
        _attach_published(item)
    return dict(fold=item)


def quick_fip(opts: Dict[str, Any]) -> Dict[str, Any]:
    """a quick FIP of what the page shows (without the signals it
    subtracted), started in a thread"""
    data, _, _ = selection(opts)
    if data is None:
        raise ValueError('no velocities to look at')
    data = subtracted(data, opts)
    # the quick FIPs before it are out of date: stopped, not waited for
    for other in QUICKS.values():
        if other.get('status') == 'running' and not other.get('batch'):
            other['cancel'] = True
    qid = uuid.uuid4().hex[:8]
    QUICKS[qid] = dict(id=qid, status='running', step='waiting',
                       step_detail='', each=[],
                       subtracted=[dict(kind=mod.get('kind', 'sine'),
                                        period=mod.get('period'),
                                        label=mod.get('label', ''))
                                   for mod in opts.get('subtract') or []],
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
    # a quick FIP that runs is stopped, and forgotten (a batch's are kept)
    for qid, job in list(QUICKS.items()):
        if job.get('batch'):
            continue
        if job.get('status') == 'running':
            job['cancel'] = True
        QUICKS.pop(qid)
        _QUICK_DATA.pop(qid, None)
    karchive._TABLES = None
    return dict(kept=[job.id for job in JOBS.values()])


#: the batch FIPs of this session
BATCHES: Dict[str, Dict[str, Any]] = {}


def batch_fip(paths: List[str], opts: Dict[str, Any]) -> Dict[str, Any]:
    """
    The quick FIP of many files, one after the other, in a thread: each
    file on its own, with no SIMBAD name (no archive, no known planet),
    the trend as the report's boxes say; each kept as a quick FIP of its
    own, to be opened in the page

    :param paths: list of str, the files of velocities
    :param opts: dict, the options of the page (trend, curvature)

    :return: dict, the state of the batch (batch_state)
    """
    paths = [os.path.abspath(os.path.expanduser(str(path)))
             for path in paths if str(path).strip()]
    if not paths:
        raise ValueError('no file to run the FIP of')
    bid = uuid.uuid4().hex[:8]
    BATCHES[bid] = dict(id=bid, status='running', start=time.time(),
                        end=None, cancel=False, trend=trend_order(opts),
                        items=[dict(path=path, name=os.path.basename(path),
                                    status='waiting', qid=None, error=None,
                                    summary=None) for path in paths])
    threading.Thread(target=_run_batch, args=(bid,), daemon=True).start()
    return batch_state(bid)


def _batch_summary(result: Dict[str, Any], data) -> Dict[str, Any]:
    """one line of the table of a batch: the best peak, its fold, the
    acceleration"""
    peak = (result.get('peak_list') or [None])[0]
    folds = result.get('folds') or []
    best = next((item for item in folds if peak and item['id'] == peak['id']),
                None)
    acc = (result.get('acceleration') or {}).get('accel')
    return dict(
        n=result.get('n'), nexp=result.get('nexp'),
        instruments=result.get('instruments'),
        baseline=float(data.baseline),
        period=peak['period'] if peak else None,
        fip=peak['family'] if peak else None,
        fip_alone=peak['alone'] if peak else None,
        K=best['K'] if best else None, K_err=best['K_err'] if best else None,
        rms=best['rms'] if best else None,
        accel=acc[0] if acc else None,
        accel_err=0.5 * (acc[1] + acc[2]) if acc else None,
        accel_sigma=(result.get('acceleration') or {}).get('accel_sigma'))


def _run_batch(bid: str) -> None:
    """the files of a batch, one quick FIP after the other"""
    batch = BATCHES[bid]
    for item in batch['items']:
        if batch['cancel']:
            item['status'] = 'stopped'
            continue
        item['status'] = 'running'
        try:
            data, _, _ = selection(dict(files=[dict(path=item['path'])]))
            if data is None:
                raise ValueError('no velocity in the file')
            qid = uuid.uuid4().hex[:8]
            QUICKS[qid] = dict(id=qid, status='running', step='waiting',
                               step_detail='', each=[], subtracted=[],
                               progress=None, result=None, error=None,
                               start=time.time(), end=None, batch=bid,
                               instruments=list(data.instruments))
            item['qid'] = qid
            _run_quick(qid, data, '', batch['trend'], each=False)
            job = QUICKS[qid]
            item['status'] = job['status']
            item['error'] = job.get('error')
            if job.get('result'):
                item['summary'] = _batch_summary(job['result'], data)
        except Exception as err:
            item['status'] = 'failed'
            item['error'] = f'{type(err).__name__}: {err}'
    batch['status'] = 'stopped' if batch['cancel'] else 'done'
    batch['end'] = time.time()


def batch_state(bid: str) -> Dict[str, Any]:
    """what the page shows of a batch: each file, the one running with
    its step and its sweeps"""
    batch = BATCHES[bid]
    items = []
    for item in batch['items']:
        one = {key: item[key] for key in ('path', 'name', 'status', 'qid',
                                          'error', 'summary')}
        job = QUICKS.get(item['qid']) if item['qid'] else None
        if job is not None and item['status'] == 'running':
            one.update(step=job.get('step'), progress=job.get('progress'))
        items.append(one)
    return dict(id=bid, status=batch['status'], items=items,
                elapsed=(batch['end'] or time.time()) - batch['start'])


def stop_batch(bid: str) -> Dict[str, Any]:
    """a batch stopped: the file running and those after it"""
    batch = BATCHES[bid]
    batch['cancel'] = True
    for item in batch['items']:
        job = QUICKS.get(item['qid']) if item['qid'] else None
        if job is not None and job.get('status') == 'running':
            job['cancel'] = True
    return batch_state(bid)


def batch_open(bid: str, index: int) -> Dict[str, Any]:
    """one file of a batch for the page, as a result recalled: the page
    (the file, no star), the velocities, its quick FIP"""
    batch = BATCHES[bid]
    item = batch['items'][int(index)]
    if not item['qid'] or item['qid'] not in QUICKS:
        raise ValueError(f'{item["name"]}: no quick FIP to open')
    detailed = dict(trend=batch['trend'] >= 1, curvature=batch['trend'] >= 2,
                    dace=False, carmenes=False, exclude='')
    page = dict(target='', files=[dict(path=item['path'], label='')],
                root='', outdir='', detailed=detailed, clip=False, view=None,
                periods=None, subtract=[])
    return dict(page=page, rv=velocities([dict(path=item['path'])]),
                quick=quick_state(item['qid']), notes=[],
                entry=dict(id='', target=item['name'], note='',
                           created=time.strftime('%Y-%m-%d %H:%M')))


def quick_state(qid: str) -> Dict[str, Any]:
    """what the page shows of a quick FIP"""
    job = QUICKS[qid]
    return dict(job, elapsed=(job['end'] or time.time()) - job['start'])


def signal_at(mod: Dict[str, Any], time: np.ndarray,
              period: Optional[float] = None) -> np.ndarray:
    """the signal of a fold's solution alone (its sinusoid, or its
    Keplerian orbit) at some times [m/s]"""
    from koloa import kepler
    time = np.asarray(time, dtype=float)
    period = float(mod.get('period') or period)
    if mod.get('kind') == 'kepler':
        return kepler.rv_keplerian(time, period, mod['tp'], mod['e'],
                                   mod['omega'], mod['K'])
    arg = 2 * np.pi * (time - mod['tref']) / period
    return mod['cos'] * np.cos(arg) + mod['sin'] * np.sin(arg)


def model_at(item: Dict[str, Any], inst: str, time: np.ndarray,
             period: float) -> Optional[np.ndarray]:
    """the solution of a fold (its offset for the instrument, its trend,
    its sinusoid or Keplerian orbit) at some times [m/s], None for an
    instrument it lacks"""
    mod = item.get('model') or {}
    if inst not in mod.get('offsets', {}):
        return None
    time = np.asarray(time, dtype=float)
    span = (time - mod.get('ttrend', mod.get('tref', 0.0))) / \
        mod.get('tscale', 365.25)
    out = mod['offsets'][inst] + signal_at(mod, time, period)
    for order, val in enumerate(mod.get('trend', []), start=1):
        out = out + val * span ** order
    return out


def fold_kepler(data, period: float, trend: int = 1,
                width: Optional[float] = None,
                transit: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    The series folded on a Keplerian orbit, its eccentricity free: koloa's
    fit (RVModel: outliers, a jitter per instrument, the trend), its period
    free within the peak (half its width in frequency each side), phase 0
    at the conjunction; errors from draws of the Laplace covariance

    :param data: RVData, the nightly means
    :param period: float, the period of the peak [days]
    :param trend: int, the order of the trend
    :param width: float or None, the width of a peak in frequency [1/days]
                  (1 / the baseline)
    :param transit: dict or None, a transit ephemeris (fold_transit's
                    transit): its period and time of transit as gaussian
                    priors, and phase 0 at its transit nearest the data

    :return: dict, as fold() gives, with kind='kepler', the orbit (P, e,
             omega, tp) and the errors of P, K and e
    """
    from koloa import kepler
    from koloa.fit import RVModel
    width = width or 1.0 / max(data.baseline, 1.0)
    freq = 1.0 / period
    plo = max(1.0 / (freq + 0.5 * width), 0.9 * period)
    phi = min(1.0 / max(freq - 0.5 * width, 1e-9), 1.1 * period)
    seq_jitter = ('instrument' if len(data.instruments) > 1
                  and data.nseq < data.n else None)
    orbit = dict(period=period, eccentric=True, period_range=(plo, phi))
    if transit is not None:
        orbit.update(period=transit['P'], tc=transit['t0'],
                     tc_err=max(transit['t0_err'], 1e-5),
                     period_err=max(float(transit.get('P_err') or 0.0),
                                    1e-7 * transit['P']))
    model = RVModel(data, [orbit],
                    likelihood='mixture', unit='both', trend=trend,
                    seq_jitter=seq_jitter)
    res = model.fit(nstart=4, quiet=True)
    per, tperi, ecc, omega, amp = model.orbit(res.theta, 0)
    # the errors: the orbit of draws of the Laplace covariance
    errs = dict(P=np.nan, K=np.nan, e=np.nan)
    draws, orbits = None, None
    if res.cov is not None and np.all(np.isfinite(res.cov)):
        rng = np.random.default_rng(11)
        try:
            draws = rng.multivariate_normal(res.theta, res.cov, 400)
            orbits = np.array([model.orbit(theta, 0) for theta in draws])
            for key, col in (('P', 0), ('e', 2), ('K', 4)):
                low, high = np.percentile(orbits[:, col], [15.87, 84.13])
                errs[key] = float(0.5 * (high - low))
        except (ValueError, np.linalg.LinAlgError):
            draws, orbits = None, None
    shown = data.rv - model.systematics(res.theta)
    # phase 0: the conjunction of the orbit, or the transit
    tconj = (transit['t0'] if transit is not None
             else kepler.tp_to_tc(tperi, per, ecc, omega))
    phase = ((data.time - tconj) / per) % 1.0
    valid = res.reliability
    sig = kepler.rv_keplerian(data.time, per, tperi, ecc, omega, amp)
    resid = shown - sig
    good = valid >= 0.5
    grid = np.linspace(0, 1, 401)
    insts = list(data.instruments)
    # the 1-sigma envelope of the curve: the orbits of the draws, at the
    #   times of the phases of the best one
    times = tconj + grid * per
    best = kepler.rv_keplerian(times, per, tperi, ecc, omega, amp)
    if orbits is not None:
        many = np.array([kepler.rv_keplerian(times, *orb) for orb in orbits])
        low, high = envelope(best, many)
    else:
        low = high = best
    out = dict(kind='kepler', period=float(per), P_err=errs['P'],
               K=float(amp), K_err=errs['K'], e=float(ecc),
               e_err=errs['e'], omega=float(np.degrees(omega) % 360),
               tp=float(tperi), tc=float(tconj),
               rms=float(np.std(resid[good] if good.any() else resid)),
               curve=dict(phase=grid.tolist(), rv=np.round(best, 4).tolist(),
                          lo=np.round(low, 4).tolist(),
                          hi=np.round(high, 4).tolist()),
               instruments=[],
               model=dict(kind='kepler', period=float(per), tp=float(tperi),
                          e=float(ecc), omega=float(omega), K=float(amp),
                          offsets={str(inst): float(res.theta[
                              model.index[f'offset_{inst}']])
                              for inst in insts},
                          trend=[float(res.theta[model.index[f'trend_{deg}']])
                                 for deg in range(1, trend + 1)],
                          ttrend=float(model.tref),
                          tscale=float(max(data.baseline, 1e-9))))
    # draws of the solution, for its envelope on the series
    out['draws'] = []
    for theta in (draws[:NDRAW] if draws is not None else []):
        try:
            dper, dtp, decc, dom, damp = model.orbit(theta, 0)
        except (ValueError, FloatingPointError):
            continue
        if not 0 <= decc < 1:
            continue
        out['draws'].append(dict(
            out['model'], period=float(dper), tp=float(dtp), e=float(decc),
            omega=float(dom), K=float(damp),
            offsets={str(inst): float(theta[model.index[f'offset_{inst}']])
                     for inst in insts},
            trend=[float(theta[model.index[f'trend_{deg}']])
                   for deg in range(1, trend + 1)]))
    berv = _berv(data)
    for inst in insts:
        sel = data.inst == inst
        out['instruments'].append(dict(
            name=str(inst), phase=np.round(phase[sel], 5).tolist(),
            rv=np.round(shown[sel], 3).tolist(),
            err=np.round(data.err[sel], 3).tolist(),
            time=np.round(data.time[sel], 5).tolist(),
            valid=np.round(valid[sel], 4).tolist(),
            berv=(_listed(berv[sel], 3) if berv is not None
                  and np.any(np.isfinite(berv[sel])) else None)))
    return out


def _mark_offscale(ax, x, y, yr) -> None:
    """red triangles at the edge of a plot, pointing to the points beyond
    its velocity range"""
    if not yr:
        return
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    pad = 0.025 * (yr[1] - yr[0])
    for sel, edge, mark in ((y > yr[1], yr[1] - pad, '^'),
                            (y < yr[0], yr[0] + pad, 'v')):
        if np.any(sel):
            ax.plot(x[sel], np.full(int(sel.sum()), edge), mark, ms=5,
                    color='#d62728', mec='none', clip_on=False, zorder=6)


def _date_ticks(t0: float, t1: float, most: int = 6):
    """calendar dates between two times [rjd]: their times and labels"""
    import matplotlib.dates as mdates
    loc = mdates.AutoDateLocator(minticks=2, maxticks=most)
    nums = loc.tick_values(mdates.num2date(t0 - 40587.5),
                           mdates.num2date(t1 - 40587.5))
    fmt = '%Y' if t1 - t0 > 3 * 365.25 else '%Y-%m' if t1 - t0 > 75 \
        else '%Y-%m-%d'
    ticks = [num + 40587.5 for num in nums if t0 <= num + 40587.5 <= t1]
    return ticks, [mdates.num2date(tt - 40587.5).strftime(fmt)
                   for tt in ticks]


def _quicklook_figures(data, source, quick, xr, yr, pr, title, each=(),
                       fold_colour: str = 'inst',
                       overlay: Optional[int] = None,
                       fold_model: str = 'sine',
                       series_colour: str = 'inst',
                       fip_view: str = 'each', series_zero: str = 'fit'):
    """the figures of the quick look: the velocities shown (with the
    solution of a fold on them when asked), the quick FIP with its peaks
    named, the folds at them (coloured by instrument, date or BERV)"""
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
    # the points by instrument, or by their BERV (on one scale)
    berv = _berv(data) if series_colour == 'berv' else None
    bnorm = None
    if berv is not None and np.any(np.isfinite(berv)):
        top = max(float(np.nanmax(np.abs(berv))), 1e-3)
        bnorm = matplotlib.colors.Normalize(-top, top)
    # each instrument about the offset the quick look fitted (one trend
    #   across them), or about its median
    offsets = ((quick or {}).get('offsets') or {}) if series_zero == 'fit' \
        else {}
    zeros = {inst: offsets.get(str(inst), np.median(data.rv[data.inst == inst]))
             for inst in data.instruments}
    for inst in data.instruments:
        sel = data.inst == inst
        label = f'{inst} ({source.get(inst, "")}, {int(sel.sum())})'
        rel = data.rv[sel] - zeros[inst]
        if bnorm is None:
            ax.errorbar(data.time[sel], rel, data.err[sel], fmt=marker[inst],
                        ms=3.5, lw=0.6, color=colour[inst], label=label)
            continue
        ax.errorbar(data.time[sel], rel, data.err[sel], fmt='none', lw=0.5,
                    color='0.7', zorder=1)
        if np.any(np.isfinite(berv[sel])):
            # a thin edge: a BERV near zero is near white
            ax.scatter(data.time[sel], rel, c=berv[sel], cmap='RdBu_r',
                       norm=bnorm, s=14, marker=marker[inst], zorder=2,
                       edgecolors='0.35', linewidths=0.3, label=label,
                       plotnonfinite=True)
        else:
            ax.plot(data.time[sel], rel, marker[inst], ms=3.5, color='0.6',
                    zorder=2, label=label)
    if bnorm is not None:
        bar = fig.colorbar(matplotlib.cm.ScalarMappable(bnorm, 'RdBu_r'),
                           ax=ax, pad=0.01, fraction=0.03)
        bar.set_label('BERV [km s$^{-1}$]')
    ax.axhline(0, color='0.6', lw=0.6, ls=':')
    # the trend the quick look fitted, the instruments about their offsets
    trend_model = (quick or {}).get('trend_model') or {}
    if offsets and (trend_model.get('coefs') or trend_model.get('signals')):
        from koloa import kepler
        lo, hi = xr if xr else (data.time.min(), data.time.max())
        shortest = min([sig['period'] for sig in
                        trend_model.get('signals') or []] + [hi - lo])
        grid = np.linspace(lo, hi, int(min(20000, max(
            400, 30 * (hi - lo) / max(shortest, 1e-3)))))
        span = (grid - trend_model['tref']) / trend_model['tscale']
        drift = sum((val * span ** (deg + 1) for deg, val in
                     enumerate(trend_model.get('coefs') or [])),
                    np.zeros(len(grid)))
        if trend_model.get('coefs'):
            ax.plot(grid, drift, color='0.35', lw=0.8, ls='--', zorder=1)
        # the model with its signals, where they are resolved (a few tens
        #   of cycles at most: beyond, a band that would hide the points)
        if trend_model.get('signals') and (hi - lo) / shortest <= 60:
            full = drift + sum(kepler.rv_keplerian(
                grid, sig['period'], sig['tp'], sig['e'], sig['omega'],
                sig['K']) for sig in trend_model['signals'])
            ax.plot(grid, full, color='0.55', lw=0.5, zorder=1)
    # the folds as the page shows them: their sinusoid, or their Keplerian
    #   orbit when it was fitted
    if quick and fold_model == 'kepler':
        quick = dict(quick, folds=[dict(item['kepler'], id=item['id'])
                                   if 'kepler' in item else item
                                   for item in quick.get('folds') or []])
    # the solution of a fold on the series, each instrument about its median
    shown = next((item for item in (quick or {}).get('folds') or []
                  if overlay is not None and item['id'] == overlay), None)
    if shown:
        for inst in data.instruments:
            sel = data.inst == inst
            lo, hi = data.time[sel].min(), data.time[sel].max()
            npts = int(min(20000, max(400, 30 * (hi - lo) / shown['period'])))
            grid = np.linspace(lo - 0.01 * (hi - lo + 1),
                               hi + 0.01 * (hi - lo + 1), npts)
            mod = model_at(shown, str(inst), grid, shown['period'])
            if mod is None:
                continue
            zero = zeros[inst]
            # its 1-sigma envelope, from the draws of the solution
            many = [model_at(dict(model=draw), str(inst), grid,
                             shown['period'])
                    for draw in shown.get('draws') or []]
            many = np.array([val for val in many if val is not None])
            if len(many) > 5:
                low, high = envelope(mod, many)
                ax.fill_between(grid, low - zero, high - zero, lw=0,
                                color='0.85', zorder=0)
            ax.plot(grid, mod - zero, lw=0.6, alpha=0.8, color=colour[inst],
                    zorder=1)
            pub = (shown.get('published') or {}).get('model')
            if pub:
                ax.plot(grid, model_at(dict(model=pub), str(inst), grid,
                                       pub['period']) - zero, lw=0.7,
                        ls='--', color='#d97706', zorder=1)
        ax.text(0.99, 0.02, f'the solution of #{shown["id"]} '
                f'({shown["period"]:.4f} d)', transform=ax.transAxes,
                ha='right', fontsize=7, color='0.3')
    if xr:
        ax.set_xlim(*xr)
    if yr:
        ax.set_ylim(*yr)
        rel = np.concatenate([data.rv[data.inst == inst] - zeros[inst]
                              for inst in data.instruments])
        times = np.concatenate([data.time[data.inst == inst]
                                for inst in data.instruments])
        inside = ((times >= xr[0]) & (times <= xr[1]) if xr
                  else np.ones(len(times), bool))
        _mark_offscale(ax, times[inside], rel[inside], yr)
    ax.set_xlabel('BJD - 2400000')
    ax.set_ylabel('RV - median [m s$^{-1}$]')
    ax.legend(fontsize=7, ncol=3, frameon=False, loc='upper left')
    # the calendar dates on top (rjd = JD - 2400000 is 40587.5 at
    #   1970-01-01 0h, matplotlib's date 0)
    import matplotlib.dates as mdates
    top = ax.secondary_xaxis('top', functions=(lambda rjd: rjd - 40587.5,
                                               lambda day: day + 40587.5))
    locator = mdates.AutoDateLocator(minticks=3, maxticks=9)
    top.xaxis.set_major_locator(locator)
    top.xaxis.set_major_formatter(mdates.AutoDateFormatter(locator))
    top.tick_params(labelsize=8)
    ax.set_title(f'{title}: the velocities shown', fontsize=10, pad=20)
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
    for item in (valid if fip_view == 'each' else []):
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
        # the colour of the points: their instrument, their date or their
        #   BERV (on one scale for every fold)
        alls = [inst for item in folds for inst in item['instruments']]
        if fold_colour == 'berv' and not any(inst.get('berv')
                                             for inst in alls):
            fold_colour = 'inst'
        norm, cmap, label = None, None, ''
        if fold_colour == 'date':
            times = np.concatenate([inst['time'] for inst in alls])
            norm = matplotlib.colors.Normalize(times.min(), times.max())
            cmap, label = plt.get_cmap('viridis'), 'date'
        elif fold_colour == 'berv':
            vals = np.array([val for inst in alls for val in
                             (inst.get('berv') or []) if val is not None])
            top = max(float(np.max(np.abs(vals))), 1e-3)
            norm = matplotlib.colors.Normalize(-top, top)
            cmap, label = plt.get_cmap('RdBu_r'), 'BERV [km s$^{-1}$]'
        for ax, item in zip(axes.ravel(), folds):
            for inst in item['instruments']:
                if norm is None:
                    ax.errorbar(inst['phase'], inst['rv'], inst['err'],
                                fmt=marker.get(inst['name'], 'o'), ms=3,
                                lw=0.5, color=colour.get(inst['name'], 'k'))
                    continue
                ax.errorbar(inst['phase'], inst['rv'], inst['err'],
                            fmt='none', lw=0.5, color='0.7', zorder=1)
                vals = (inst['time'] if fold_colour == 'date'
                        else inst.get('berv'))
                if vals is None:
                    ax.plot(inst['phase'], inst['rv'],
                            marker.get(inst['name'], 'o'), ms=3,
                            color='0.6', zorder=2)
                    continue
                vals = np.array([np.nan if val is None else val
                                 for val in vals], dtype=float)
                ax.scatter(inst['phase'], inst['rv'], c=vals, cmap=cmap,
                           norm=norm, s=12, zorder=2,
                           marker=marker.get(inst['name'], 'o'),
                           edgecolors='0.35', linewidths=0.3,
                           plotnonfinite=True)
                # circled: less than an even chance of being valid
                low = np.asarray(inst.get('valid', []), float) < 0.5
                if low.any():
                    ax.plot(np.asarray(inst['phase'])[low],
                            np.asarray(inst['rv'])[low], 'o', mfc='none',
                            mec='k', ms=7, mew=0.8, zorder=5)
            # the 1-sigma envelope of the fit, light grey
            if item['curve'].get('lo'):
                ax.fill_between(item['curve']['phase'], item['curve']['lo'],
                                item['curve']['hi'], color='0.85', lw=0,
                                zorder=0)
            ax.plot(item['curve']['phase'], item['curve']['rv'], color='k',
                    lw=1.0)
            pub = item.get('published')
            if pub:
                ax.plot(pub['phase'], pub['rv'], ls='--', color='#d97706',
                        lw=1.0, label=f'{pub["name"]}, {pub["reference"]}')
                ax.legend(fontsize=6, frameon=False, loc='lower left')
            named = (item.get('published') or {}).get('name')
            ax.set_title(f'#{item["id"]}' + (f' {named}' if named else '')
                         + f': P = {item["period"]:.4f} d, '
                         f'K = {item["K"]:.2f} $\\pm$ {item["K_err"]:.2f} '
                         f'm/s' + (f', e = {item["e"]:.2f} $\\pm$ '
                                   f'{item["e_err"]:.2f}'
                                   if item.get('kind') == 'kepler' else ''),
                         fontsize=8.5)
            ax.set_xlabel('phase (0 = conjunction)')
            ax.set_ylabel('RV [m s$^{-1}$]')
            # the velocities shown as the series is
            if yr:
                ax.set_ylim(*yr)
                _mark_offscale(ax, [ph for inst in item['instruments']
                                    for ph in inst['phase']],
                               [val for inst in item['instruments']
                                for val in inst['rv']], yr)
        if norm is None:
            fig.tight_layout()
        else:
            # the scale of the colours in a strip of its own, on the right
            fig.tight_layout(rect=(0, 0, 0.92, 1))
            cax = fig.add_axes((0.935, 0.15, 0.012, 0.7))
            bar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap),
                               cax=cax)
            bar.set_label(label)
            if fold_colour == 'date':
                ticks, labels = _date_ticks(norm.vmin, norm.vmax)
                bar.set_ticks(ticks)
                bar.set_ticklabels(labels)
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
        series='The velocities shown, each instrument about the offset the '
               'quick look fitted (one trend across them, dashed) or about '
               'its median, as on the page, in its ranges (the calendar '
               'dates on top), with the '
               'solution of a fold when the page shows it; a red triangle '
               'at an edge points to a velocity beyond the range.',
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
              'being valid (an outlier, as the FIP saw it); coloured as on '
              'the page (by instrument, date or BERV); a red triangle at an '
              'edge points to a night beyond the range; light grey, the '
              '1-$\\sigma$ envelope of the fit; dashed orange, the published '
              'orbit of a known planet.',
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
        if quick.get('subtracted'):
            out.append('The FIP of the residuals: the series without '
                       + '; '.join(escape(item.get('label') or
                                          f'{item["period"]:.4f} d')
                                   for item in quick['subtracted'])
                       + ' (each fold\'s signal subtracted; its offsets '
                       'and trend fitted again).\n')
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
        # the folds at a period of the page's own, and the Keplerian orbits
        forced = [item for item in quick.get('folds') or []
                  if item.get('forced')]
        if forced:
            out.append('Folded at a period asked: ' + '; '.join(
                f'\\#{item["id"]} '
                + (f'{escape(item["known"]["name"])} (known) '
                   if item.get('known') else '')
                + f'{item["period"]:.4f}\\,d (K = '
                f'${item["K"]:.2f} \\pm {item["K_err"]:.2f}$\\,m/s)'
                for item in forced) + '.\n')
        for item in quick.get('folds') or []:
            eph = item.get('transit')
            if not eph:
                continue
            out.append(
                f'\\#{item["id"]} {escape(eph["name"])} on its transit '
                f'ephemeris ({escape(eph.get("reference") or "")}): '
                f'$P = {eph["P"]:.7f}$\\,d, the transit nearest the '
                f'velocities $T_0 = {eph["t0"]:.5f} \\pm '
                f'{eph["t0_err"] * 1440:.1f}$\\,min ({eph["cycles"]} '
                f'periods from {eph["tc"]:.5f}; phase $\\pm '
                f'{eph["phase_err"]:.4f}$); $K = {item["K"]:.2f} \\pm '
                f'{item["K_err"]:.2f}$\\,m/s with its phase fixed; the '
                f'velocities alone put the conjunction '
                f'${eph["shift"] * 24:+.2f} \\pm {eph["shift_err"] * 24:.2f}'
                f'$\\,h from it'
                + (f' ({eph["shift_sigma"]:.1f}\\,$\\sigma$)'
                   if eph.get('shift_sigma') is not None else '') + '.\n')
        # the minimum masses: the mass of the star the page gives
        mstar, mstar_err = opts.get('mstar'), opts.get('mstar_err') or 0.0

        def masses(item, transit=False):
            if not mstar:
                return ''
            from koloa.stars import minimum_mass
            mm = minimum_mass(item['K'], item['period'],
                              item.get('e') or 0.0, float(mstar),
                              item.get('K_err') or 0.0,
                              item.get('P_err') or 0.0,
                              item.get('e_err') or 0.0, float(mstar_err))
            return ((r'mass ($\sin i \approx 1$) ' if transit
                     else r'$m \sin i$ ') + ' = '.join(
                f'${val:.3g}_{{-{low:.2g}}}^{{+{high:.2g}}}$\\,{unit}'
                for (val, low, high), unit in (
                    (mm['earth'], r'M$_\oplus$'),
                    (mm['neptune'], r'M$_\mathrm{Nep}$'),
                    (mm['jupiter'], r'M$_\mathrm{Jup}$')))
                + f' (M$_\\star$ = {float(mstar):.3f} $\\pm$ '
                  f'{float(mstar_err):.3f}\\,M$_\\odot$)')
        for item in quick.get('folds') or []:
            if item.get('transit') and mstar:
                out.append(f'\\#{item["id"]} {escape(item["transit"]["name"])}'
                           f': {masses(item, transit=True)}.\n')
        keps = [dict(item['kepler'], id=item['id'])
                for item in quick.get('folds') or [] if 'kepler' in item]
        if keps:
            out.append('{\\small\\begin{tabular}{@{}rrrrrr@{}}\n\\toprule\n'
                       'Keplerian & P [d] & K [m/s] & $e$ & $\\omega$ [deg] '
                       '& rms [m/s] \\\\\n\\midrule')
            for item in keps:
                out.append(f'\\#{item["id"]} & ${item["period"]:.4f} \\pm '
                           f'{item["P_err"]:.4f}$ & ${item["K"]:.2f} \\pm '
                           f'{item["K_err"]:.2f}$ & ${item["e"]:.2f} \\pm '
                           f'{item["e_err"]:.2f}$ & {item["omega"]:.0f} & '
                           f'{item["rms"]:.2f} \\\\')
            out.append('\\bottomrule\n\\end{tabular}}\n\n'
                       'Keplerian orbits: koloa\'s fit (outliers, a jitter '
                       'per instrument, the trend), the eccentricity free, '
                       'the period free within the peak; errors from the '
                       'Laplace covariance.\n')
            for item in keps if mstar else []:
                out.append(f'\\#{item["id"]}: '
                           f'{masses(item, transit=bool(item.get("transit")))}'
                           f'.\n')
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
                  qid: str = '', command_line: str = '',
                  fold_colour: str = 'inst',
                  overlay: Optional[int] = None,
                  fold_model: str = 'sine',
                  series_colour: str = 'inst',
                  fip_view: str = 'each', series_zero: str = 'fit') -> bytes:
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
    figs = _quicklook_figures(data, source, quick, xr, yr, pr, title, each,
                              fold_colour, overlay, fold_model,
                              series_colour, fip_view, series_zero)
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
                                         ('carmenes', 'CARMENES DR1'),
                                         ('vizier', 'VizieR'))
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
        # rv/ whole: DACE, CARMENES and the published velocities
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
                    carmenes=bool(detailed.get('carmenes')),
                    vizier=bool(detailed.get('vizier')))
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
    # a result kept before the quick FIP listed the star's transits
    result = QUICKS[qid]['result']
    if 'transits' not in result:
        result['transits'] = transits_of(page.get('target', ''),
                                         result.get('known') or [])
    try:
        if _refresh_folds(qid, page):
            notes.append('folds kept before their dates, BERV and solution '
                         'were: made again from the series')
            # kept so, once and for all
            with open(os.path.join(folder, 'quick.json'), 'w') as handle:
                json.dump(dict(quick, result=QUICKS[qid]['result']), handle)
    except (ValueError, OSError) as err:  # a help, not a need
        notes.append(f'the folds as they were kept ({err})')
    return dict(page=page, rv=loaded['rv'], quick=quick_state(qid),
                entry={key: entry.get(key) for key in
                       ('id', 'target', 'note', 'created', 'summary')},
                notes=notes)


def page_options(page: Dict[str, Any]) -> Dict[str, Any]:
    """the options of the quick look (the page's shownOptions) from the
    state of a page remembered"""
    detailed = page.get('detailed') or {}
    return dict(files=page.get('files') or [], target=page.get('target', ''),
                root=page.get('root', ''), dace=bool(detailed.get('dace')),
                carmenes=bool(detailed.get('carmenes')),
                vizier=bool(detailed.get('vizier')),
                exclude=str(detailed.get('exclude') or ''),
                trend=detailed.get('trend', True),
                curvature=bool(detailed.get('curvature')),
                subtract=page.get('subtract') or [])


def _refresh_folds(qid: str, page: Dict[str, Any]) -> int:
    """the folds of a result remembered before they carried their dates,
    BERV and solution, made again from its series (each night's
    probability of being valid as it was kept)

    :return: int, how many were made again
    """
    result = QUICKS[qid]['result']
    old = [item for item in result.get('folds') or []
           if 'model' not in item or not all('time' in inst
                                             for inst in item['instruments'])]
    if not old:
        return 0
    series = _quick_series(qid, page_options(page))
    nights = series['nights']
    for item in old:
        # the probability of each night, as the fold kept it, where the
        #   nights are the same (elsewhere, valid)
        valid, found = np.ones(nights.n), False
        for inst in item['instruments']:
            sel = np.flatnonzero(nights.inst == inst['name'])
            if inst.get('valid') and len(inst['valid']) == len(sel):
                valid[sel], found = inst['valid'], True
        new = fold(nights, item['period'], valid if found else None,
                   series['trend'])
        keep = {key: item[key] for key in ('id', 'forced', 'kepler')
                if key in item}
        item.clear()
        item.update(new, **keep)
    return len(old)


def unremember(rid: str) -> Dict[str, Any]:
    """a result remembered, forgotten (its folder and its copies)"""
    folder = _entry_folder(rid)
    shutil.rmtree(folder)
    return dict(id=rid)


# =============================================================================
# The server
# =============================================================================
def _finite(value: Any) -> Any:
    """a value for JSON: its nans and infinities None"""
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, dict):
        return {key: _finite(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(val) for val in value]
    return value


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
        try:
            text = json.dumps(value, default=str, allow_nan=False)
        except ValueError:
            # a nan or an infinity: null (JSON has neither)
            text = json.dumps(_finite(value), default=str)
        self._send(code, text.encode(), 'application/json')

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
            if url.path == '/api/batch':
                return self._json(batch_state(query['id']))
            if url.path == '/api/listfiles':
                return self._json(list_files(query.get('folder', ''),
                                             query.get('pattern', '*.rdb')))
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
                    carmenes=query.get('carmenes') == '1',
                    vizier=query.get('vizier') == '1'))
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
                pdf = quicklook_pdf(dict(body.get('options', {}),
                                         mstar=body.get('mstar'),
                                         mstar_err=body.get('mstar_err')),
                                    body.get('x'),
                                    body.get('y'), body.get('p'),
                                    body.get('quick', ''),
                                    body.get('command', ''),
                                    body.get('fold_colour') or 'inst',
                                    body.get('overlay'),
                                    body.get('fold_model') or 'sine',
                                    body.get('series_colour') or 'inst',
                                    body.get('fip_view') or 'each',
                                    body.get('series_zero') or 'fit')
                self.send_response(200)
                self.send_header('Content-Type', 'application/pdf')
                self.send_header('Content-Length', str(len(pdf)))
                self.send_header('Content-Disposition',
                                 'attachment; filename="koloa_quicklook.pdf"')
                self.end_headers()
                self.wfile.write(pdf)
                return None
            if path == '/api/analysis_script':
                # the analysis kit: the script of the analysis with the
                #   page's settings, its data, its star (koloa.kit)
                from koloa.gather import folder_name
                from koloa.kit import build
                opts = dict(body.get('options', {}),
                            mstar=body.get('mstar'),
                            mstar_err=body.get('mstar_err'))
                kit = build(opts, body.get('lang') or 'en')
                target = str(opts.get('target') or '').strip()
                name = (f'koloa_{folder_name(target) if target else "series"}'
                        f'_analysis.tar.gz')
                self.send_response(200)
                self.send_header('Content-Type', 'application/gzip')
                self.send_header('Content-Length', str(len(kit)))
                self.send_header('Content-Disposition',
                                 f'attachment; filename="{name}"')
                self.end_headers()
                self.wfile.write(kit)
                return None
            if path == '/api/forget':
                return self._json(forget())
            if path == '/api/quickstop':
                return self._json(stop_quick(body.get('id', '')))
            if path == '/api/batch':
                return self._json(batch_fip(body.get('paths') or [],
                                            body.get('options') or {}))
            if path == '/api/batchstop':
                return self._json(stop_batch(body.get('id', '')))
            if path == '/api/batch_open':
                return self._json(batch_open(body.get('id', ''),
                                             body.get('index', 0)))
            if path == '/api/fold':
                return self._json(fold_request(
                    body.get('quick', ''), body.get('options', {}),
                    body.get('period'), body.get('id'),
                    body.get('kind') or 'sine', bool(body.get('snap')),
                    body.get('known') or '', body.get('transit') or ''))
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
        # a quick FIP running: its chains ended, not left behind
        running = [job for job in QUICKS.values()
                   if job.get('status') == 'running']
        for job in running:
            job['cancel'] = True
        if running:
            time.sleep(1.5)
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

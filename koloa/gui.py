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


def command(action: str, opts: Dict[str, Any]) -> List[str]:
    """
    The arguments of koloa's command line for an action of the page

    :param action: str, gather or detailed
    :param opts: dict, the fields of the page (target, file, root, outdir,
                 the switches and numbers)

    :return: list of str, the arguments after 'koloa'
    """
    target = str(opts.get('target') or '').strip()
    rvfile = str(opts.get('file') or '').strip()
    off = {key: not opts.get(key, True) for key in
           ('dace', 'carmenes', 'tess', 'vizier', 'archive', 'gpcheck',
            'fip_gp', 'duck', 'latex')}
    if action == 'gather':
        if not target:
            raise ValueError('a SIMBAD name to gather the archives of')
        args = [target, '--gather', str(opts.get('root') or 'archives')]
        args += ['--no-dace'] * off['dace'] + ['--no-carmenes'] * off[
            'carmenes'] + ['--no-tess'] * off['tess']
        return args
    if action != 'detailed':
        raise ValueError(f'no action {action}')
    if not target and not rvfile:
        raise ValueError('a file, a SIMBAD name, or (best) both')
    args = ([rvfile] if rvfile else []) + ['--detailed']
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
    dmap = opts.get('detection_map') or 'none'
    if dmap == 'fip':
        args.append('--detection-map')
    elif dmap == 'search':
        args.append('--search-map')
    if opts.get('exposures'):
        args.append('--exposures')
    if opts.get('mcmc'):
        args.append('--mcmc')
    for key, flag in (('fip_gp', '--no-fip-gp'), ('dace', '--no-dace'),
                      ('carmenes', '--no-carmenes'), ('tess', '--no-tess'),
                      ('vizier', '--no-vizier'), ('archive', '--no-archive'),
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
        env = dict(os.environ, PYTHONUNBUFFERED='1')
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
            self.lines.append(text)
            message = text.split(' | ', 1)[-1]
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
def resolve_star(name: str) -> Dict[str, Any]:
    """the SIMBAD resolver of the page: identifiers, position, TIC, the
    periods SIMBAD lists, CARMENES DR1"""
    from koloa.archive import resolve
    from koloa.gather import carmenes_star, folder_name, variability
    ident = resolve(name)
    out = dict(ident, folder=folder_name(name), variability=[],
               carmenes=None)
    try:
        out['variability'] = variability(ident['main'])
    except Exception as err:  # a help, not a need
        out['variability_error'] = str(err)
    if ident.get('ra') is not None:
        try:
            out['carmenes'] = carmenes_star(ident['ra'], ident['dec'])
        except Exception as err:
            out['carmenes_error'] = str(err)
    return out


def velocities(rvfile: str = '', target: str = '', root: str = ''
               ) -> Dict[str, Any]:
    """
    The velocities of a file and of a star's gathered archives, by
    instrument (each with its median taken out), for the plot of the page
    """
    from koloa.data import merge
    from koloa.detailed import _read
    from koloa.gather import folder_name, load
    series, notes = [], []
    if rvfile:
        series.append(_read(rvfile, None))
        notes.append(f'{os.path.basename(rvfile)}: {series[-1].n} points')
    if target:
        folder = os.path.join(root or 'archives', folder_name(target))
        if os.path.exists(os.path.join(folder, 'rv', 'all_rv.csv')):
            gathered = load(folder)['rv']
            if series:
                keep = ~np.isin(gathered.inst, series[0].instruments)
                gathered = gathered.select(keep) if keep.any() else None
            if gathered is not None:
                series.append(gathered)
                notes.append(f'{folder}: {gathered.n} points')
        else:
            notes.append(f'nothing gathered in {folder} yet')
    if not series:
        return dict(instruments=[], notes=notes)
    data = series[0] if len(series) == 1 else merge(series)
    out = []
    for name in data.instruments:
        sel = data.inst == name
        out.append(dict(name=name, n=int(sel.sum()),
                        time=np.round(data.time[sel], 6).tolist(),
                        rv=np.round(data.rv[sel], 3).tolist(),
                        err=np.round(data.err[sel], 3).tolist(),
                        rms=float(np.std(data.rv[sel]))))
    return dict(instruments=out, notes=notes, n=int(data.n),
                baseline=float(data.baseline))


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
                return self._json(dict(cwd=os.getcwd(),
                                       python=sys.executable,
                                       defaults=DEFAULTS))
            if url.path == '/api/resolve':
                return self._json(resolve_star(query.get('name', '')))
            if url.path == '/api/rv':
                return self._json(velocities(query.get('file', ''),
                                             query.get('target', ''),
                                             query.get('root', '')))
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
                else:
                    outputs = opts.get('outdir') or 'koloa_output'
                job = Job(body['action'], args, outputs)
                JOBS[job.id] = job
                return self._json(job.state())
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

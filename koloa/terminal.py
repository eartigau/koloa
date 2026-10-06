#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A terminal in koloa's GUI: a shell of this machine in a pseudo-terminal,
shown and typed in by the page, to go to the machine that will run a batch
(ssh, its password or its second factor typed as in any terminal), unpack
the batch there and launch it.

How one got there is remembered as a route, under a name: the lines typed
(ssh maestria, cd /data/me/batches, what loads Python), where the batches
go there, and the host as ssh and rsync name it. One button types them
again. Only what the terminal showed as it was typed is proposed for a
route (a password is typed with the echo off: it is never kept), and the
lines are shown, to be corrected, before they are kept
(~/.config/koloa/routes.json).

A shell is more than koloa's page is otherwise trusted with (it runs
nothing but koloa): the terminal answers only to the page opened with the
key koloa printed when it started (its address ends with ?key=...), as a
notebook server does. Another user of this machine cannot reach it.

Not on Windows (no pseudo-terminal).

Created on 2026-10-06

@author: artigau
"""
import base64
import json
import os
import re
import secrets
import signal
import struct
import threading
import time
from typing import Any, Dict, List, Optional

# =============================================================================
# Define variables
# =============================================================================
#: the key of this server: the terminal answers only to who has it
KEY = secrets.token_urlsafe(18)
#: the terminals of this session: id -> Session
SESSIONS: Dict[str, 'Session'] = {}
#: where the routes are kept
ROUTES = os.path.join(os.path.expanduser('~'), '.config', 'koloa',
                      'routes.json')
#: the most output kept of a terminal [bytes] (the page has the rest)
KEPT = 2_000_000
#: how long a question for output waits for some [s]
WAIT = 20.0
#: a line typed again by a route: the next one when the terminal has been
#: quiet this long, or after this long at most [s]
QUIET, PATIENCE = 0.8, 30.0


# =============================================================================
# Define classes
# =============================================================================
class Session:
    """a shell in a pseudo-terminal: what it wrote, what was typed"""

    def __init__(self, cols: int = 100, rows: int = 28,
                 cwd: Optional[str] = None):
        import pty
        self.id = secrets.token_hex(4)
        self.start = time.time()
        self.buffer = b''
        self.base = 0          # the offset of buffer[0] in all it wrote
        self.cond = threading.Condition()
        self.alive = True
        self.last = time.time()   # when it last wrote
        self.typed: List[str] = []   # the lines typed that it showed
        self._line = ''        # the line being typed
        self._from = 0         # where the output was when it began
        shell = os.environ.get('SHELL') or '/bin/sh'
        self.pid, self.fd = pty.fork()
        if self.pid == 0:  # the shell
            os.environ['TERM'] = 'xterm-256color'
            try:
                os.chdir(cwd or os.path.expanduser('~'))
            except OSError:
                pass
            os.execvp(shell, [shell, '-l'])
        self.resize(cols, rows)
        threading.Thread(target=self._read, daemon=True).start()

    def _read(self) -> None:
        """what the shell writes, kept as it comes"""
        while True:
            try:
                data = os.read(self.fd, 65536)
            except OSError:
                data = b''
            with self.cond:
                if not data:
                    self.alive = False
                    self.cond.notify_all()
                    return
                self.buffer += data
                if len(self.buffer) > KEPT:
                    cut = len(self.buffer) - KEPT
                    self.buffer, self.base = self.buffer[cut:], self.base + cut
                self.last = time.time()
                self.cond.notify_all()

    @property
    def end(self) -> int:
        """how much it wrote so far [bytes]"""
        return self.base + len(self.buffer)

    def read(self, since: int = 0, wait: float = WAIT) -> Dict[str, Any]:
        """what it wrote from an offset on, waited for when there is none
        yet: data (base64), next (the offset after it), alive"""
        with self.cond:
            if since >= self.end and self.alive:
                self.cond.wait(timeout=wait)
            start = max(since, self.base) - self.base
            data = self.buffer[start:]
            return dict(data=base64.b64encode(data).decode(), next=self.end,
                        alive=self.alive)

    def write(self, text: str) -> None:
        """keys typed, passed to the shell; the lines they make are kept
        when the terminal showed them"""
        if not self.alive:
            raise ValueError('this terminal has ended')
        # where its output was before these keys: their echo comes after
        before = self.end
        os.write(self.fd, text.encode())
        for char in text:
            if char in '\r\n':
                self._keep()
            elif char in '\x7f\x08':
                self._line = self._line[:-1]
            elif char == '\x15':   # ctrl-U: the line erased
                self._line = ''
            elif char >= ' ' and char != '\x1b':
                if not self._line:
                    self._from = before
                self._line += char

    def _keep(self) -> None:
        """the line just typed, kept if the terminal showed it: with the
        echo off (a password, a code) it showed nothing of it"""
        line, self._line = self._line.strip(), ''
        if not line:
            return
        # the echo may come a moment after the last key
        time.sleep(0.15)
        with self.cond:
            shown = self.buffer[max(self._from, self.base) - self.base:]
        text = re.sub(rb'\x1b\[[0-9;?]*[A-Za-z]', b'', shown).decode(
            'utf-8', 'replace')
        # a long line is shown on several rows: compared without spaces
        if re.sub(r'\s+', '', line) in re.sub(r'\s+', '', text):
            self.typed.append(line)

    def resize(self, cols: int, rows: int) -> None:
        """the size of the terminal of the page"""
        import fcntl
        import termios
        try:
            fcntl.ioctl(self.fd, termios.TIOCSWINSZ,
                        struct.pack('HHHH', int(rows), int(cols), 0, 0))
        except OSError:
            pass

    def type_lines(self, lines: List[str]) -> None:
        """lines typed one after the other (a route), each once the
        terminal has been quiet for a moment (a login takes its time; a
        password asked on the way is typed by hand, and the lines go on
        after it)"""
        def work():
            for line in lines:
                began = time.time()
                while self.alive and time.time() - began < PATIENCE:
                    if time.time() - self.last > QUIET:
                        break
                    time.sleep(0.1)
                if not self.alive:
                    return
                os.write(self.fd, (line + '\r').encode())
                self.last = time.time()
        threading.Thread(target=work, daemon=True).start()

    def close(self) -> None:
        """the shell ended, with what it runs here (not what it started
        on another machine)"""
        try:
            os.kill(self.pid, signal.SIGHUP)
        except OSError:
            pass
        try:
            os.close(self.fd)
        except OSError:
            pass
        self.alive = False


# =============================================================================
# Define functions
# =============================================================================
def available() -> bool:
    """whether this machine has pseudo-terminals (not Windows)"""
    try:
        import pty  # noqa: F401
        import termios  # noqa: F401
    except ImportError:
        return False
    return hasattr(os, 'fork')


def routes() -> Dict[str, Dict[str, Any]]:
    """the routes kept: name -> lines (typed to get there), host (as ssh
    and rsync name it), folder (where the batches go there), made"""
    if not os.path.exists(ROUTES):
        return {}
    try:
        with open(ROUTES) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return {}


def save_route(name: str, lines: Any, host: str = '', folder: str = ''
               ) -> Dict[str, Dict[str, Any]]:
    """a route kept under a name (its lines: a text, a line each, or a
    list), in place of one of that name"""
    name = str(name).strip()
    if not name:
        raise ValueError('a name for the route')
    if isinstance(lines, str):
        lines = lines.splitlines()
    kept = routes()
    kept[name] = dict(lines=[str(line).strip() for line in lines
                             if str(line).strip()],
                      host=str(host).strip(), folder=str(folder).strip(),
                      made=time.strftime('%Y-%m-%d %H:%M'))
    os.makedirs(os.path.dirname(ROUTES), exist_ok=True)
    with open(ROUTES + '.part', 'w') as handle:
        json.dump(kept, handle, indent=1)
    os.replace(ROUTES + '.part', ROUTES)
    return kept


def delete_route(name: str) -> Dict[str, Dict[str, Any]]:
    """a route forgotten"""
    kept = routes()
    kept.pop(str(name).strip(), None)
    if os.path.exists(ROUTES):
        with open(ROUTES + '.part', 'w') as handle:
            json.dump(kept, handle, indent=1)
        os.replace(ROUTES + '.part', ROUTES)
    return kept


def _session(sid: Any) -> Session:
    if sid not in SESSIONS:
        raise ValueError('no such terminal (koloa was started again?)')
    return SESSIONS[sid]


def route(path: str, body: Dict[str, Any], key: Optional[str]
          ) -> Dict[str, Any]:
    """
    The questions of the terminal (/api/term/...), for who has the key

    :param path: str, the path asked
    :param body: dict, what is asked
    :param key: str or None, the key given (the X-Koloa-Key header)
    """
    what = path.rsplit('/', 1)[-1]
    if what == 'state':
        # no key needed to know whether there is a terminal to ask for
        return dict(available=available(),
                    allowed=bool(key) and secrets.compare_digest(
                        str(key), KEY))
    if not key or not secrets.compare_digest(str(key), KEY):
        raise PermissionError(
            'the terminal answers only to the page opened from the address '
            'koloa printed when it started (it ends with ?key=...)')
    if not available():
        raise ValueError('no terminal on this system')
    if what == 'open':
        session = Session(int(body.get('cols') or 100),
                          int(body.get('rows') or 28), body.get('cwd'))
        SESSIONS[session.id] = session
        return dict(id=session.id)
    if what == 'routes':
        return dict(routes=routes())
    if what == 'route_save':
        return dict(routes=save_route(body.get('name', ''),
                                      body.get('lines') or [],
                                      body.get('host', ''),
                                      body.get('folder', '')))
    if what == 'route_delete':
        return dict(routes=delete_route(body.get('name', '')))
    session = _session(body.get('id'))
    if what == 'read':
        return session.read(int(body.get('since') or 0),
                            float(body.get('wait') or WAIT))
    if what == 'write':
        session.write(str(body.get('data') or ''))
        return dict(ok=True)
    if what == 'resize':
        session.resize(int(body.get('cols') or 100),
                       int(body.get('rows') or 28))
        return dict(ok=True)
    if what == 'typed':
        return dict(lines=list(session.typed))
    if what == 'go':
        kept = routes().get(str(body.get('name') or ''))
        if kept is None:
            raise ValueError('no such route')
        session.type_lines(kept['lines'])
        return dict(ok=True, lines=kept['lines'])
    if what == 'close':
        session.close()
        SESSIONS.pop(session.id, None)
        return dict(ok=True)
    raise ValueError(f'no such question: {path}')


def close_all() -> None:
    """every terminal ended (the server stops)"""
    for session in list(SESSIONS.values()):
        session.close()
    SESSIONS.clear()


# =============================================================================
# End of code
# =============================================================================

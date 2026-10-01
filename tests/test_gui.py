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
    args = gui.command('detailed', dict(target='GL 406', rotation='2.704',
                                        fip_gp='sho'))
    assert args[-2:] == ['--rotation', '2.704']
    assert gui.command('detailed', dict(target='GL 406', rotation='2.704',
                                        fip_gp='banded'))[-2:] == [
        '--fip-gp', 'banded']
    with pytest.raises(ValueError):
        gui.command('detailed', dict(target='GL 406', fip_gp='sho'))
    assert gui.command('detailed', dict(target='x', fip_gp='none'))[-1] == (
        '--no-fip-gp')
    assert gui.command('gather', dict(target='x', refresh=True))[-1] == (
        '--refresh')
    # the name alone, or the file alone, will do; neither will not
    assert gui.command('detailed', dict(target='GJ 436'))[0] == '--detailed'
    assert gui.command('detailed', dict(file='a.rdb'))[0] == 'a.rdb'
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

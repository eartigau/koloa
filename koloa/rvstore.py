#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The velocities koloa gathers, kept for a team in a private repository
(github.com/eartigau/koloa-rvs): its collaborators read them, and those
who choose to add theirs.

A star is a folder of the repository named by SIMBAD's main identifier
(GJ 436 and Ross 905 are stars/Ross_905/), with what koloa keeps of it:

    stars/<star>/target.json      the star (its names, its astrometry)
    stars/<star>/manifest.json    what each archive gave (no local path)
    stars/<star>/rv/...           the velocities: dace/, carmenes/,
                                  published/ (VizieR, the observatories)
    stars/<star>/lbl/...          the rdb of LBL (koloa.lblrun), labels.json

The light curves are not kept (phot/, several MB a star).

Reading: when the archives of a star are gathered, what the repository has
of it and the folder here does not is copied into the folder first, so
that it is read back as if it had been fetched here (the velocities of
DACE asked with the key of a collaborator, say). Nothing here is ever
replaced by what the repository has.

Writing: when you keep your velocities there (rvs.json among koloa's
settings, or the box of the GUI), the files of a star that the repository
does not have, or has otherwise, are added to it after a gather or a run
of LBL, in one commit. A star whose names are among those never to share
(koloa.shared) is never written.

Both go through the command line of GitHub (gh, logged in) and the
repository's API, with no copy of it here; without gh, or without access
to the repository, koloa does as if there were none.

Created on 2026-10-14

@author: artigau
"""
import base64
import hashlib
import json
import os
import tempfile
import time as _time
from typing import Any, Dict, List, Optional, Sequence, Tuple

from koloa.log import log
from koloa.paths import config

# =============================================================================
# Define variables
# =============================================================================
#: the repository, and its branch
REPO = 'eartigau/koloa-rvs'
BRANCH = 'main'
#: your choice: whether to add your velocities to it
SETTINGS = config('rvs.json')
#: what of the folder of a star is kept: its files, and those of these
#: folders (all their files)
FILES = ('target.json', 'manifest.json')
FOLDERS = ('rv', 'lbl')
#: the largest file kept [bytes] (GitHub's API takes up to 100 MB)
LARGEST = 50 * 1024 * 1024
#: the tries of a commit another one got ahead of
TRIES = 5
#: what is known of the repository in this process (asked once)
_ACCESS: Dict[str, Any] = {}


# =============================================================================
# Define functions
# =============================================================================
def settings() -> Dict[str, Any]:
    """your choice (SETTINGS): push, whether to add your velocities"""
    out: Dict[str, Any] = dict(push=False)
    if os.path.exists(SETTINGS):
        try:
            with open(SETTINGS) as handle:
                out.update(json.load(handle))
        except (OSError, ValueError):
            pass
    return out


def save_settings(push: bool) -> Dict[str, Any]:
    """your choice, kept"""
    out = dict(settings(), push=bool(push))
    os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
    with open(SETTINGS, 'w') as handle:
        json.dump(out, handle, indent=1)
    return out


def _gh(*args: str, body: Optional[Dict[str, Any]] = None,
        timeout: float = 120.0) -> Tuple[int, Any]:
    """the API of GitHub through gh: its code and its answer (JSON)"""
    from koloa.shared import _gh as run
    path = None
    if body is not None:
        handle, path = tempfile.mkstemp(suffix='.json')
        with os.fdopen(handle, 'w') as out:
            json.dump(body, out)
        args = args + ('--input', path)
    try:
        done = run('api', *args, timeout=timeout)
    finally:
        if path:
            os.remove(path)
    if done is None:
        return -1, None
    try:
        answer = json.loads(done.stdout or b'null')
    except ValueError:
        answer = None
    return done.returncode, answer


def access() -> Dict[str, bool]:
    """whether the repository can be read here, and written (gh logged
    in, its collaborators); asked once a process"""
    if not _ACCESS:
        code, answer = _gh(f'repos/{REPO}', timeout=30.0)
        perms = (answer or {}).get('permissions') or {} if code == 0 else {}
        _ACCESS.update(read=bool(perms.get('pull')),
                       write=bool(perms.get('push')))
    return dict(_ACCESS)


def key_of(ident: Dict[str, Any]) -> str:
    """the folder of a star in the repository: its main identifier"""
    from koloa.gather import folder_name
    return folder_name(str(ident.get('main') or ''))


def _tree() -> Tuple[Optional[str], Dict[str, str]]:
    """the last commit of the branch, and the files of the repository with
    the identifiers of their contents"""
    code, ref = _gh(f'repos/{REPO}/git/ref/heads/{BRANCH}')
    if code != 0 or not ref:
        return None, {}
    head = ref['object']['sha']
    code, tree = _gh(f'repos/{REPO}/git/trees/{head}?recursive=1')
    files = {one['path']: one['sha'] for one in (tree or {}).get('tree', [])
             if one.get('type') == 'blob'} if code == 0 else {}
    return head, files


def _blob_sha(data: bytes) -> str:
    """the identifier git gives a content"""
    return hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest()


def kept(folder: str) -> Dict[str, bytes]:
    """what of the folder of a star the repository keeps, by its path
    there (the local path of the manifest taken out)"""
    out: Dict[str, bytes] = {}
    for name in FILES:
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            with open(path, 'rb') as handle:
                data = handle.read()
            if name == 'manifest.json':
                try:
                    manifest = json.loads(data)
                    manifest.pop('folder', None)
                    data = json.dumps(manifest, indent=1).encode()
                except ValueError:
                    pass
            out[name] = data
    for sub in FOLDERS:
        for here, _, files in os.walk(os.path.join(folder, sub)):
            for name in files:
                path = os.path.join(here, name)
                if name.startswith('.') or name.endswith('.part') or \
                        os.path.getsize(path) > LARGEST:
                    continue
                with open(path, 'rb') as handle:
                    out[os.path.relpath(path, folder).replace(os.sep, '/')
                        ] = handle.read()
    return out


def pull(ident: Dict[str, Any], folder: str) -> int:
    """
    What the repository has of a star and the folder does not, copied into
    the folder (nothing there replaced)

    :return: int, the files copied
    """
    if not access().get('read'):
        return 0
    prefix = f'stars/{key_of(ident)}/'
    _, files = _tree()
    mine = {path[len(prefix):]: sha for path, sha in files.items()
            if path.startswith(prefix)}
    copied = 0
    for rel, sha in mine.items():
        path = os.path.join(folder, *rel.split('/'))
        if os.path.exists(path):
            continue
        code, blob = _gh(f'repos/{REPO}/git/blobs/{sha}')
        if code != 0 or not blob:
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'wb') as handle:
            handle.write(base64.b64decode(blob.get('content') or ''))
        copied += 1
    if copied:
        log(f'koloa-rvs: {copied} files of {ident.get("main")} copied from '
            f'the repository', 'value')
    return copied


def push(ident: Dict[str, Any], folder: str,
         names: Sequence[str] = (), always: Optional[bool] = None) -> int:
    """
    The files of the folder of a star the repository does not have, or
    has otherwise, added to it in one commit, when you keep your
    velocities there; never for a star among those never to share

    :param names: more names of the star (the one typed)
    :param always: True to write whatever your choice (a batch packed to
                   do so), False never to, None as you chose

    :return: int, the files added or changed
    """
    from koloa.shared import never
    if not (settings().get('push') if always is None else always):
        return 0
    every = [ident.get('main') or ''] + list(ident.get('aliases') or []) + \
        list(names)
    if never(every):
        log(f'koloa-rvs: {ident.get("main")} is among the names never to '
            f'share: not written', 'info')
        return 0
    if not access().get('write'):
        return 0
    prefix = f'stars/{key_of(ident)}/'
    ours = kept(folder)
    for _ in range(TRIES):
        head, files = _tree()
        if head is None:
            return 0
        new = {rel: data for rel, data in ours.items()
               if files.get(prefix + rel) != _blob_sha(data)}
        if not new:
            return 0
        entries = []
        for rel, data in new.items():
            code, blob = _gh(f'repos/{REPO}/git/blobs', '--method', 'POST',
                             body=dict(content=base64.b64encode(
                                 data).decode(), encoding='base64'))
            if code != 0 or not blob:
                log(f'koloa-rvs: {rel} not written', 'warn')
                return 0
            entries.append(dict(path=prefix + rel, mode='100644', type='blob',
                                sha=blob['sha']))
        code, commit = _gh(f'repos/{REPO}/git/commits/{head}')
        base = (commit or {}).get('tree', {}).get('sha')
        code, tree = _gh(f'repos/{REPO}/git/trees', '--method', 'POST',
                         body=dict(base_tree=base, tree=entries))
        if code != 0 or not tree:
            return 0
        code, made = _gh(f'repos/{REPO}/git/commits', '--method', 'POST',
                         body=dict(message=f'{ident.get("main")}: '
                                   f'{len(new)} files',
                                   tree=tree['sha'], parents=[head]))
        if code != 0 or not made:
            return 0
        code, _ = _gh(f'repos/{REPO}/git/refs/heads/{BRANCH}', '--method',
                      'PATCH', body=dict(sha=made['sha'], force=False))
        if code == 0:
            log(f'koloa-rvs: {len(new)} files of {ident.get("main")} kept in '
                f'the repository', 'value')
            return len(new)
        # another commit got there first: asked again from it
        _time.sleep(1.0)
    log(f'koloa-rvs: {ident.get("main")} not written (the repository moved '
        f'on {TRIES} times)', 'warn')
    return 0


def push_all(root: str) -> List[Tuple[str, int]]:
    """every star gathered in a folder (a folder of folders with their
    target.json), its velocities added to the repository

    :return: list of (the folder, the files added)
    """
    out = []
    for name in sorted(os.listdir(root)):
        folder = os.path.join(root, name)
        target = os.path.join(folder, 'target.json')
        if not os.path.isfile(target):
            continue
        with open(target) as handle:
            ident = json.load(handle)
        if not ident.get('main'):
            continue
        out.append((name, push(ident, folder, [name.replace('_', ' ')])))
    return out

# =============================================================================
# End of code
# =============================================================================

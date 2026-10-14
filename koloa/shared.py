#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The velocities of HARPS-N that koloa's users fetched, shared: a public
repository (github.com/eartigau/koloa-harpsn) with, for each star, what
the headers of its public CCF files in the archive of the TNG say (one
small table, rdb), and an index of where the stars are.

Reading needs nothing: before it asks the archive of the TNG, koloa looks
for the star in the index (kept a day), by where it is (MATCH), and reads
its table. The list of the pointings of HARPS-N (koloa.tng.pointings)
says how many public files there are at the star: when the table has as
many, nothing is asked of the archive; when it has fewer, the archive is
asked as before, and only the files the table lacks are fetched.

Sharing is yours to choose, and needs the command line of GitHub logged
in (gh auth login): without it, the velocities are fetched and nothing
more. With it, and sharing on (share.json among koloa's settings, or
KOLOA_SHARE=1), a star whose files were fetched from the archive is asked
to be added: an issue of the repository, under your name on GitHub, that
says where the star is and nothing else. The repository's Action then
reads the archive itself and adds what it finds: nothing a user sends is
taken as it is. A star of your list of names never to share is never
asked for.

Created on 2026-10-14

@author: artigau
"""
import json
import os
import re
import shutil
import subprocess
import time as _time
import urllib.error
import urllib.request
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from koloa.log import log
from koloa.paths import cache, config

# =============================================================================
# Define variables
# =============================================================================
#: the repository, and where anyone reads its files
REPO = 'eartigau/koloa-harpsn'
RAW = f'https://raw.githubusercontent.com/{REPO}/main/'
#: the index, as it was last read, and for how long it is kept [days]
INDEX = cache('tng', 'shared_index.json')
INDEX_AGE = 1.0
#: how far from a star of the index a star may be, where both were in 2000
#: [arcsec] (the coordinates of SIMBAD, which differ little from one
#: user to the next)
MATCH = 2.0
#: the columns of the table of a star: what koloa.tng.series reads of the
#: headers, and the file each line is of
FIELDS = ['file', 'reduction', 'object', 'bjd', 'rv', 'err', 'fwhm',
          'contrast', 'bis', 'mask', 'version', 'berv', 'dpr', 'qc']
#: your choices: whether to share, and the names never to share
SETTINGS = config('share.json')
#: the date of an exposure, from the name of its file
_DATE = re.compile(r'HARPN\.(\d{4}-\d{2}-\d{2})T(\d{2})-(\d{2})-([\d.]+)')


# =============================================================================
# Define functions
# =============================================================================
def _gh(*args: str, timeout: float = 60.0
        ) -> Optional[subprocess.CompletedProcess]:
    """the command line of GitHub (gh), None when it is not there"""
    exe = shutil.which('gh')
    if not exe:
        return None
    try:
        return subprocess.run([exe, *args], capture_output=True,
                              timeout=timeout)
    except (OSError, subprocess.SubprocessError):
        return None


def _get(path: str, timeout: float = 20.0) -> Optional[bytes]:
    """
    A file of the repository: as anyone reads it, else through gh (for a
    repository not public, read by those it is shared with); None when it
    cannot be read
    """
    try:
        with urllib.request.urlopen(urllib.request.Request(
                RAW + path, headers={'User-Agent': 'koloa'}),
                timeout=timeout) as response:
            return response.read()
    except urllib.error.HTTPError as err:
        if err.code != 404:
            return None
    except (urllib.error.URLError, TimeoutError, OSError):
        return None
    done = _gh('api', f'repos/{REPO}/contents/{path}', '-H',
               'Accept: application/vnd.github.raw', timeout=timeout)
    return done.stdout if done is not None and done.returncode == 0 \
        else None


def index(refresh: bool = False) -> List[Dict[str, Any]]:
    """
    The stars of the repository (index.json), kept INDEX_AGE days; the
    last one read when it cannot be read again, none when never read

    :return: list of dict: id, name, ra, dec [deg, 2000], pmra, pmdec
             [mas/yr], path (its table), files (its lines), pointed (the
             public files at the star in the list of the pointings when
             it was made), updated
    """
    kept = None
    if os.path.exists(INDEX):
        try:
            with open(INDEX) as handle:
                kept = json.load(handle)
        except (OSError, ValueError):
            kept = None
    fresh = kept is not None and not refresh and \
        _time.time() - kept.get('asked', 0) < INDEX_AGE * 86400
    if fresh:
        return list(kept.get('stars') or [])
    raw = _get('index.json')
    try:
        stars = list(json.loads(raw)['stars']) if raw else None
    except (ValueError, KeyError, TypeError):
        stars = None
    if stars is None:
        # not read: the one kept, or none (asked again in a day)
        stars = list((kept or {}).get('stars') or [])
    os.makedirs(os.path.dirname(INDEX), exist_ok=True)
    with open(INDEX, 'w') as handle:
        json.dump(dict(asked=_time.time(), stars=stars), handle)
    return stars


def entry(star: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """the star of the index where a star is (within MATCH), the nearest;
    None when there is none"""
    from koloa.tng import _offsets
    stars = index()
    if not stars or star.get('ra') is None:
        return None
    east, north = _offsets([one['ra'] for one in stars],
                           [one['dec'] for one in stars], star)
    far = np.hypot(east, north)
    best = int(np.argmin(far))
    return stars[best] if far[best] < MATCH else None


def date_of(name: str) -> str:
    """the start of an exposure (UTC), from the name of its file"""
    found = _DATE.search(name)
    return (f'{found[1]}T{found[2]}:{found[3]}:{found[4]}' if found
            else '')


def read_rdb(text: str) -> List[Dict[str, Any]]:
    """the lines of a table of the repository, as koloa.tng.exposures
    gives them (public files: policy FREE)"""
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        return []
    names = lines[0].split('\t')
    out = []
    for line in lines[2:]:
        row = dict(zip(names, line.split('\t')))
        if not row.get('file'):
            continue
        out.append(dict(row, policy='FREE', date=date_of(row['file']),
                        error=''))
    return out


def write_rdb(rows: Sequence[Dict[str, Any]]) -> str:
    """a table of the repository: FIELDS, tab-separated, a line of dashes
    under their names (rdb), a line per file, by date"""
    def cell(value: Any) -> str:
        return ' '.join(str('' if value is None else value).split())
    lines = ['\t'.join(FIELDS), '\t'.join('-' * len(name) for name in FIELDS)]
    for row in sorted(rows, key=lambda one: (date_of(one['file']),
                                             one['file'])):
        lines.append('\t'.join(cell(row.get(name)) for name in FIELDS))
    return '\n'.join(lines) + '\n'


def rows(star: Dict[str, Any]) -> Tuple[Optional[Dict[str, Any]],
                                         List[Dict[str, Any]]]:
    """
    What the repository has of a star: its entry in the index, and the
    lines of its table (none when it has none, or it cannot be read)
    """
    one = entry(star)
    if one is None:
        return None, []
    raw = _get(one['path'])
    if not raw:
        return None, []
    return one, read_rdb(raw.decode('utf8', 'replace'))


def settings() -> Dict[str, Any]:
    """your choices (SETTINGS): share (bool), never (names); KOLOA_SHARE
    (1 or 0) says whether to share, whatever the file says"""
    out: Dict[str, Any] = dict(share=False, never=[])
    if os.path.exists(SETTINGS):
        try:
            with open(SETTINGS) as handle:
                out.update(json.load(handle))
        except (OSError, ValueError):
            pass
    env = os.environ.get('KOLOA_SHARE', '').strip().lower()
    if env:
        out['share'] = env in ('1', 'true', 'yes', 'on')
    return out


def save_settings(share: bool, never: Sequence[str]) -> Dict[str, Any]:
    """your choices, kept (SETTINGS)"""
    out = dict(share=bool(share),
               never=[str(name).strip() for name in never
                      if str(name).strip()])
    os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
    with open(SETTINGS, 'w') as handle:
        json.dump(out, handle, indent=1)
    return out


def _bare(name: str) -> str:
    """a name as names are compared: no SIMBAD mark, no space, upper case
    (* alf Cen and ALFCEN are one)"""
    return re.sub(r'^(NAME|\*\*|\*|V\*)\s+', '', str(name).strip()).replace(
        ' ', '').upper()


def never(names: Sequence[str]) -> bool:
    """whether one of the names of a star is among those never to share"""
    mine = {_bare(name) for name in settings().get('never') or []} - {''}
    return bool(mine & {_bare(name) for name in names if name})


def can_share() -> bool:
    """whether gh is there and logged in"""
    done = _gh('auth', 'status', timeout=20.0)
    return done is not None and done.returncode == 0


def request(star: Dict[str, Any], name: Optional[str] = None,
            fetched: int = 0) -> str:
    """
    Ask the repository to add what the archive has of a star (an issue:
    where the star is, and its name), when sharing is on, gh logged in,
    and the star not among those never to share

    :param star: dict, ra, dec [deg, 2000], pmra, pmdec [mas/yr], and its
                 names (names)
    :param name: str or None, its name
    :param fetched: int, the files fetched from the archive, which the
                    repository did not have

    :return: str, what was done, in words ('' when nothing)
    """
    from koloa import __version__
    if not settings().get('share'):
        return ''
    names = [name or ''] + list(star.get('names') or [])
    if never(names):
        return 'not shared: among the names never to share'
    if not can_share():
        return ''
    body = json.dumps(dict(
        name=name or '', ra=float(star['ra']), dec=float(star['dec']),
        pmra=float(star.get('pmra') or 0.0),
        pmdec=float(star.get('pmdec') or 0.0), fetched=int(fetched),
        koloa=__version__))
    done = _gh('issue', 'create', '--repo', REPO, '--title',
               f'star: {name or "?"}', '--body', body)
    if done is None or done.returncode != 0:
        why = (done.stderr.decode('utf8', 'replace').strip().splitlines()
               or ['?'])[-1] if done is not None else 'gh not found'
        log(f'shared HARPS-N velocities: not asked ({why})', 'warn')
        return ''
    return 'asked to be shared'

# =============================================================================
# End of code
# =============================================================================

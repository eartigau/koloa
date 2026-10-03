#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The APERO names of stars: the astrometric database of APERO (the pipeline
of SPIRou and NIRPS), one YAML file per object, keyed by its APERO name (a
uniform name without the characters that are trouble in a file name:
GL699, PROXIMA, TOI_756, 2MASS_J00140580M3245599), with its SIMBAD name,
its aliases, its position and its properties

    from koloa import apero_names
    apero_names.lookup('Gl382')              # the entry of AN Sex
    apero_names.star_of_file('lbl_GL382_GL382.rdb')
    apero_names.refresh()                    # a new copy of the database

The database is APERO's own: the astrometrics folder of its assets, which
APERO downloads from its assets server as one tarball. Its YAML files come
first in the tarball, so a copy of them streams the first few MB of it and
stops (no APERO, no profile, no 900 MB); refresh() makes a new copy, and
the first lookup makes one when there is none. KOLOA_APERO_ASTROMETRICS
points to a folder of YAML files instead (an export made by APERO itself,
say), read as it is.

A name is found the way APERO finds it (drs_astrometrics.find_by_name):
cleaned (+ to P, - to M, any other punctuation and spaces to _, upper
case), in every variant of its spaces, signs and underscores, against the
APERO name, the original name, the SIMBAD name and every alias of each
entry; the verified entries first, then the pending, then the rejected.

The star of a file of velocities: its OBJECT column (LBL keeps the header
key of each spectrum), else its name: lbl_<OBJECT>_<TEMPLATE>.rdb (or
lbl2_) as LBL writes it, or <OBJECT>_rv_<...>.rdb; then its SIMBAD name
from the database, for the resolver and the archives.

Created on 2026-10-03

@author: artigau
"""
import csv
import json
import os
import re
import string
import tarfile
import time
import urllib.request
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: where the copy of the database is kept
FOLDER = os.path.join(os.path.expanduser('~'), '.cache', 'koloa',
                      'apero_astrometrics')
#: a folder of YAML files to read instead (an export of APERO's)
ENV_FOLDER = 'KOLOA_APERO_ASTROMETRICS'
#: the assets server of APERO (its DRS.ASSETS_URLS)
SERVER = 'http://206.12.93.77/ari/data/apero/assets/'
#: the folder of the database in the assets
SUBDIR = 'astrometrics'
#: its tiers, the first one found wins
TIERS = ('verified', 'pending', 'rejected')
#: the index of names kept beside the copy
INDEX = 'koloa_index.json'
#: what is kept of each entry
KEEP = dict(apero='APERO_NAME', original='ORIGINAL_NAME',
            simbad='SIMBAD_NAME', status='STATUS', klass='APERO_CLASS',
            gaia='GAIA_SOURCE_ID', mass_mann15='MASS_STAR_MANN15',
            mass_delfosse00='MASS_STAR_DELFOSSE00', teff_gaia='TEFF_GAIA')
#: the values kept of the entries that have a value and its source
VALUES = dict(ra='RA', dec='DEC', plx='PLX', spt='SPT', teff='TEFF',
              rv='RV', vsini='VSINI')
#: characters APERO replaces (all punctuation but _, and the space)
BAD_CHARS = [' '] + list(string.punctuation.replace('_', ''))
NULL_TEXT = ('', 'None', 'Null', 'NULL', 'null', 'nan', 'NaN', 'inf')
#: the index in memory: (signature of the folder, index)
_CACHE: Dict[str, Tuple[str, Dict[str, Any]]] = {}


# =============================================================================
# Define functions
# =============================================================================
def clean_object(raw: Any) -> str:
    """
    A name as APERO cleans it (drs_astrometrics.clean_object): stripped,
    + to p and - to m, any other punctuation and the spaces to _, upper
    case, _ collapsed and stripped; 'Null' for nothing

    :param raw: the name

    :return: str, the cleaned name
    """
    if raw is None or str(raw).strip() in NULL_TEXT:
        return 'Null'
    name = str(raw).strip().replace('+', 'p').replace('-', 'm')
    for char in BAD_CHARS:
        name = name.replace(char, '_')
    name = name.upper()
    while '__' in name:
        name = name.replace('__', '_')
    return name.strip('_')


def name_variants(raw: Any) -> List[str]:
    """
    The variants of a name APERO looks for (drs_astrometrics.
    name_search_variants): the cleaned name first, then every choice of
    the spaces (_ or none), + (P or none) and - (M or none), each also
    without its underscores

    :param raw: the name

    :return: list of str
    """
    if raw is None:
        return []
    text = str(raw).strip()
    if text in NULL_TEXT:
        return []
    keep = {' ', '\t', '+', '-', '_'}
    norm = ''.join('_' if (char in BAD_CHARS and char not in keep) else char
                   for char in text)
    seen: List[str] = []
    canonical = clean_object(raw)
    if canonical and canonical != 'Null':
        seen.append(canonical)
    for space in ('_', ''):
        for plus in ('P', ''):
            for minus in ('M', ''):
                val = norm.replace(' ', space).replace('\t', space)
                val = val.replace('+', plus).replace('-', minus).upper()
                while '__' in val:
                    val = val.replace('__', '_')
                val = val.strip('_')
                for one in (val, val.replace('_', '')):
                    if one and one not in seen:
                        seen.append(one)
    return seen


def folder() -> str:
    """the folder of the database read: KOLOA_APERO_ASTROMETRICS, else the
    copy kept here"""
    return os.path.expanduser(os.environ.get(ENV_FOLDER) or FOLDER)


def yaml_files(path: str) -> List[str]:
    """
    The YAML files of the database, as APERO lists them (iter_yaml_files):
    those of its tiers (verified, pending, rejected), then those at its
    root that no tier has

    :return: list of str
    """
    if not os.path.isdir(path):
        return []
    out, seen = [], set()
    for sub in list(TIERS) + ['']:
        here = os.path.join(path, sub)
        if not os.path.isdir(here):
            continue
        for name in sorted(os.listdir(here)):
            if name.endswith('.yaml') and name not in seen:
                seen.add(name)
                out.append(os.path.join(here, name))
    return out


def newest_tarball(server: str = SERVER, timeout: float = 30.0) -> str:
    """the newest assets tarball of APERO's server (its index lists them,
    <unix time>_<...>_assets.tar.gz)"""
    with urllib.request.urlopen(server, timeout=timeout) as resp:
        page = resp.read().decode('utf-8', 'replace')
    names = sorted(set(re.findall(r'href="(\d+_\d+_assets\.tar\.gz)"',
                                  page)),
                   key=lambda name: float(name.split('_')[0]))
    if not names:
        raise ValueError(f'no assets tarball listed at {server}')
    return names[-1]


def refresh(server: str = SERVER, timeout: float = 60.0) -> Dict[str, Any]:
    """
    A new copy of the database: the YAML files of the astrometrics folder
    of APERO's newest assets tarball, streamed (they come first in it, a
    few MB) and the rest of it left; then the index of their names

    :param server: str, the assets server of APERO
    :param timeout: float [s]

    :return: dict, the tarball, the number of entries, the folder
    """
    if os.environ.get(ENV_FOLDER):
        raise ValueError(f'{ENV_FOLDER} is set: the database is read from '
                         f'{folder()}, not copied here')
    name = newest_tarball(server, timeout)
    start = time.time()
    part = FOLDER + '.part'
    os.makedirs(part, exist_ok=True)
    count = 0
    with urllib.request.urlopen(server + name, timeout=timeout) as resp, \
            tarfile.open(fileobj=resp, mode='r|gz') as tar:
        inside = False
        for member in tar:
            parts = member.name.strip('./').split('/')
            if parts[0] != SUBDIR:
                # the folder is one stretch of the tarball: once left, done
                if inside:
                    break
                continue
            inside = True
            if not (member.isfile() and member.name.endswith('.yaml')):
                continue
            rest = parts[1:]
            # an entry at its root or in a tier, nothing else (its history)
            if len(rest) == 1 or (len(rest) == 2 and rest[0] in TIERS):
                target = os.path.join(part, *rest)
                os.makedirs(os.path.dirname(target), exist_ok=True)
                with tar.extractfile(member) as src, \
                        open(target, 'wb') as dst:
                    dst.write(src.read())
                count += 1
    if not count:
        raise ValueError(f'no astrometric entry in {name}')
    with open(os.path.join(part, 'source.json'), 'w') as handle:
        json.dump(dict(server=server, tarball=name, entries=count,
                       fetched=time.strftime('%Y-%m-%d %H:%M')), handle)
    # the new copy in place of the old one, whole
    if os.path.isdir(FOLDER):
        old = FOLDER + '.old'
        _remove(old)
        os.replace(FOLDER, old)
        os.replace(part, FOLDER)
        _remove(old)
    else:
        os.replace(part, FOLDER)
    _CACHE.pop(FOLDER, None)
    index = load_index(FOLDER)
    log(f'APERO names: {count} files of {name} ({len(index["entries"])} '
        f'objects) in {time.time() - start:.0f} s', 'value')
    return dict(tarball=name, files=count, objects=len(index['entries']),
                folder=FOLDER)


def _remove(path: str) -> None:
    """a copy of the database removed (its files, its tiers, itself)"""
    if not os.path.isdir(path):
        return
    for here, dirs, files in os.walk(path, topdown=False):
        for name in files:
            os.remove(os.path.join(here, name))
        for name in dirs:
            os.rmdir(os.path.join(here, name))
    os.rmdir(path)


def _signature(path: str) -> str:
    """what changes when a file of the database does"""
    files = yaml_files(path)
    stamp = max((os.path.getmtime(name) for name in files), default=0.0)
    return f'{len(files)}:{stamp:.3f}'


def _value(entry: Dict[str, Any], key: str) -> Any:
    """a value of an entry ({value: ..., source: ...} or as it is)"""
    val = entry.get(key)
    return val.get('value') if isinstance(val, dict) else val


def _summary(entry: Dict[str, Any], status: str) -> Dict[str, Any]:
    """what is kept of an entry"""
    out = {key: entry.get(name) for key, name in KEEP.items()}
    out.update({key: _value(entry, name) for key, name in VALUES.items()})
    out['status'] = out.get('status') or status
    for key, val in list(out.items()):
        if isinstance(val, str) and val.strip() in NULL_TEXT:
            out[key] = None
    return out


def load_index(path: Optional[str] = None) -> Dict[str, Any]:
    """
    The index of the database: each variant of each name to its APERO
    name, and each entry (kept beside the copy, built again when a file
    changes)

    :param path: str or None, the folder (folder() by default)

    :return: dict, names (variant: APERO name), entries (APERO name:
             summary)
    """
    path = path or folder()
    sig = _signature(path)
    cached = _CACHE.get(path)
    if cached is not None and cached[0] == sig:
        return cached[1]
    keep = os.path.join(path, INDEX)
    if os.path.exists(keep):
        with open(keep) as handle:
            index = json.load(handle)
        if index.get('signature') == sig:
            _CACHE[path] = (sig, index)
            return index
    try:
        import yaml
    except ImportError as err:
        raise ImportError('the APERO names need PyYAML: pip install pyyaml'
                          ) from err
    loader = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)
    names: Dict[str, str] = {}
    entries: Dict[str, Dict[str, Any]] = {}
    for name in yaml_files(path):
        try:
            with open(name) as handle:
                entry = yaml.load(handle, Loader=loader)
        except Exception:  # a malformed file: left, as APERO does
            continue
        if not isinstance(entry, dict):
            continue
        tier = os.path.basename(os.path.dirname(name))
        apero = str(entry.get('APERO_NAME') or os.path.splitext(
            os.path.basename(name))[0])
        entries.setdefault(apero, _summary(entry, tier if tier in TIERS
                                           else ''))
        aliases = entry.get('ALIASES') or []
        if isinstance(aliases, str):
            aliases = aliases.split('|')
        for value in [apero, entry.get('ORIGINAL_NAME'),
                      entry.get('SIMBAD_NAME')] + list(aliases):
            for variant in name_variants(value):
                names.setdefault(variant, apero)
    index = dict(signature=sig, names=names, entries=entries)
    try:
        with open(keep, 'w') as handle:
            json.dump(index, handle)
    except OSError:  # a folder that cannot be written: kept in memory
        pass
    _CACHE[path] = (sig, index)
    return index


def available(fetch: bool = True) -> bool:
    """whether there is a database to read (a copy made when there is none
    and fetch)"""
    if yaml_files(folder()):
        return True
    if not fetch or os.environ.get(ENV_FOLDER):
        return False
    try:
        refresh()
    except (OSError, ValueError, tarfile.TarError) as err:
        log(f'APERO names: no copy of the database ({err})', 'warn')
        return False
    return True


def lookup(name: Any, fetch: bool = True) -> Optional[Dict[str, Any]]:
    """
    The entry of a name in the database, as APERO finds it

    :param name: the name (raw: Gl382, GJ 382, AN SEX, GL382...)
    :param fetch: bool, a copy of the database made when there is none

    :return: dict (apero, original, simbad, status, ra, dec, plx, spt,
             teff, masses...) or None
    """
    if not available(fetch):
        return None
    index = load_index()
    for variant in name_variants(name):
        apero = index['names'].get(variant)
        if apero is not None:
            return dict(index['entries'][apero])
    return None


def simbad_target(entry: Dict[str, Any]) -> str:
    """the name of an entry for koloa's resolver: its SIMBAD name without
    SIMBAD's prefixes (NAME, *, V*) and repeated spaces, else its APERO
    name"""
    name = str(entry.get('simbad') or '').strip()
    name = re.sub(r'^(NAME|\*\*?|V\*)\s+', '', name)
    name = re.sub(r'\s+', ' ', name).strip()
    return name or str(entry.get('apero') or '')


def file_object(path: str) -> Optional[str]:
    """the OBJECT column of a file of velocities (LBL keeps the header key
    of each spectrum): its most common value, None without one"""
    try:
        with open(os.path.expanduser(path), newline='') as handle:
            first = handle.readline()
            delim = '\t' if '\t' in first else ','
            head = next(csv.reader([first], delimiter=delim))
            cols = [col.strip() for col in head]
            col = next((cols.index(key) for key in ('OBJECT', 'object')
                        if key in cols), None)
            if col is None:
                return None
            values = []
            for row in csv.reader(handle, delimiter=delim):
                if len(row) > col and row[col].strip() and \
                        not set(row[col].strip()) <= {'-'}:
                    values.append(row[col].strip())
    except (OSError, UnicodeDecodeError, csv.Error):
        return None
    if not values:
        return None
    names, counts = np.unique(values, return_counts=True)
    return str(names[np.argmax(counts)])


def name_candidates(path: str) -> List[str]:
    """
    The object a file's name gives: lbl_<OBJECT>_<TEMPLATE>.rdb (or lbl2_)
    as LBL writes it, an object and a template that may hold _ (each split
    of the name: the two halves first when they are the same name, then
    the splits whose template starts with the object, a template being
    named after its star, then the others, the longest object first), or
    <OBJECT>_rv_<...>.rdb

    :return: list of str, the candidates, the likeliest first
    """
    stem = os.path.splitext(os.path.basename(path))[0]
    match = re.match(r'^lbl2?_(.+)$', stem, re.IGNORECASE)
    if match:
        parts = match.group(1).split('_')
        splits = [('_'.join(parts[:cut]), '_'.join(parts[cut:]))
                  for cut in range(len(parts) - 1, 0, -1)]
        out = [obj for obj, tpl in splits if obj == tpl]
        out += [obj for obj, tpl in splits if tpl.startswith(obj)]
        out += [obj for obj, _ in splits]
        return list(dict.fromkeys(out)) or [match.group(1)]
    match = re.match(r'^(.+?)_rv(_.*)?$', stem, re.IGNORECASE)
    if match:
        return [match.group(1)]
    return []


def star_of_file(path: str, fetch: bool = True) -> Dict[str, Any]:
    """
    The star of a file of velocities: its OBJECT column, else its name,
    found in the database of APERO names (its SIMBAD name for the
    resolver); the name itself when the database does not have it

    :param path: str, the file
    :param fetch: bool, a copy of the database made when there is none

    :return: dict: raw (the name read), source (OBJECT column or file
             name), apero (the APERO name, None when not found), target
             (the name to resolve), entry (the database's, or None)
    """
    raw = file_object(path)
    tries = ([(raw, 'OBJECT column')] if raw else []) + [
        (name, 'file name') for name in name_candidates(path)]
    for name, source in tries:
        entry = lookup(name, fetch)
        if entry is not None:
            return dict(raw=name, source=source, apero=entry['apero'],
                        target=simbad_target(entry), entry=entry)
        fetch = False
    if not tries:
        return dict(raw=None, source=None, apero=None, target='',
                    entry=None)
    name, source = tries[0]
    return dict(raw=name, source=source, apero=None, target=name,
                entry=None)


# =============================================================================
# End of code
# =============================================================================

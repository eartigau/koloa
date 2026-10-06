#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A survey: the stars that meet a few constraints, what the archives have of
each, the files of velocities one has of them.

    from koloa import survey
    stars = survey.sample(sptype=('M0', 'M9'), dmax=15.0)   # asked of SIMBAD
    survey.check(stars, 'archives')       # which have velocities, and where
    survey.match_files(stars, ['/data/nirps/lbl', '/data/spirou/lbl'])
    with_data = [star for star in stars if survey.has_data(star)]

The sample is SIMBAD's: every object with a parallax above 1000 / dmax and
a spectral type in the range asked (M3.5V is 63.5: O0 = 0, B0 = 10 ... M0 =
60, L0 = 70; the type of the primary for 'M2V+M3V'), dwarfs unless said
otherwise, with its magnitudes and its names in the catalogues an archive
may know it by (GJ, HD, HIP, TIC, Gaia DR3, Karmn...). A system and its
components are objects of their own in SIMBAD: each says which other
object of the sample is within NEAR arcsec.

Then a batch of them all, here or on another machine:

    made = survey.pack(with_data, 'm_dwarfs_15pc', root='archives',
                       server_root='/scratch/me/m_dwarfs_15pc')
    # made['tar']: the files, the archives, koloa itself, run_batch.py
    survey.run('m_dwarfs_15pc', jobs=6)        # or, there: python run_batch.py

pack() writes one folder, and its .tar.gz, with everything the batch
needs and no path of this machine in it: files/ (the files of each star),
archives/ (gathered here, so that the batch asks nothing of the network),
cache/ (what koloa fetched once: the NASA Exoplanet Archive, the lists of
the surveys), koloa_src/ (koloa as it is here), targets.json and
run_batch.py, whose first setting is ROOT, where the folder is on the
machine that runs it. run() takes the stars one after the other, a few at
once (each a Python of its own), the quick FIP of each with its files and
its archives chosen by koloa.datasets; each star done is kept
(results/<star>/: result.json, quick.json, quicklook.pdf), so that a batch
stopped goes on where it was, and results/table.csv gathers them.

The check asks no more than is needed to tell which stars have
velocities: the lists of the stars of CARMENES DR1 and of the surveys on
VizieR (kept on this machine: a position looked up), and DACE, star by
star, a few at a time; what DACE answers is kept in the archives folder
of the star, where koloa.gather reads it back. The tables of a star's
papers and its light curves are for when its archives are gathered.

Created on 2026-10-06

@author: artigau
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: the spectral classes, ten subtypes each: O0 is 0, M0 60, L0 70
CLASSES = 'OBAFGKMLTY'
#: two objects of a sample this close are a system and one of its
#: components, or two names of one star [arcsec]
NEAR = 5.0
#: the catalogues whose names of a star are kept (an archive, a file, may
#: know it by any of them)
CATALOGUES = ('GJ', 'HD', 'HIP', 'TIC', 'Gaia DR3', 'Karmn', 'LHS', 'Wolf',
              'Ross', 'TOI', 'NAME', 'LP', 'G', 'L', 'BD', 'CD', '2MASS')
#: the stars of a sample DACE is asked of at once
WORKERS = 6
#: the surveys on VizieR, as an archive of a spectrograph: the
#: spectrograph, and the name of the archive in a few letters (the
#: California Legacy Survey has HIRES for every star, the APF and Lick for
#: some: known once its velocities are gathered)
SURVEY_ARCHIVES = {'teklu25': ('HIRES', 'Teklu+ 2025'),
                   'cls21': ('HIRES', 'CLS'),
                   'talor19': ('HIRES', 'Tal-Or+ 2019'),
                   'fischer14': ('Lick', 'Fischer+ 2014'),
                   'rvbank20': ('HARPS', 'RVBank')}
#: the names DACE gives a spectrograph that koloa.published.family does
#: not read
DACE_NAMES = {'HARPN': 'HARPS-N'}
#: the most objects asked of SIMBAD
MAXREC = 50000
#: the file of a batch folder that lists its stars
TARGETS = 'targets.json'
#: what of koloa's cache goes with a batch (not the light curves)
CACHE_PARTS = ('archive', 'published', 'tic', 'apero_astrometrics',
               'carmenes_objects.json')
#: the columns of the table of a batch
COLUMNS = ('name', 'sptype', 'distance', 'status', 'files', 'datasets',
           'nights', 'baseline', 'period', 'fip', 'fip_alone', 'K', 'K_err',
           'rms', 'accel', 'accel_err', 'accel_sigma', 'transit',
           'transit_snr', 'elapsed', 'error')


# =============================================================================
# Define functions
# =============================================================================
def sp_number(sptype: Any) -> Optional[float]:
    """
    A spectral type as a number: M3.5V is 63.5 (O0 = 0, B0 = 10, A0 = 20,
    F0 = 30, G0 = 40, K0 = 50, M0 = 60, L0 = 70, T0 = 80, Y0 = 90), the
    type of the primary of 'M2V+M3V', the first of 'M2/3V'; a class
    without a subtype is its middle (M: 65)

    :return: float, or None when there is no class to read (a white
             dwarf's DA, an empty type)
    """
    text = str(sptype or '').strip()
    found = re.match(r'^(?:esd|usd|sd|d|g|c)?\s*([OBAFGKMLTY])\s*'
                     r'(\d+(?:\.\d+)?)?', text)
    if not found or text.upper().startswith('D') and not text.startswith('d'):
        return None
    base = 10.0 * CLASSES.index(found.group(1))
    return base + (float(found.group(2)) if found.group(2) else 5.0)


def sp_kind(sptype: Any) -> str:
    """
    What kind of star a spectral type says: dwarf (V, a 'd' prefix, or no
    luminosity class: within a few tens of parsecs, a dwarf), subdwarf
    (VI, sd), subgiant (IV), giant (III, II, I, a 'g' prefix) or white
    dwarf (DA, DC...)
    """
    text = str(sptype or '').strip()
    if re.match(r'^D[ABCOQZX]?', text) and not text.startswith('Dw'):
        return 'white dwarf'
    if re.match(r'^(esd|usd|sd)', text):
        return 'subdwarf'
    if text.startswith('g') or text.startswith('c'):
        return 'giant'
    # the luminosity class after the type of the primary
    first = re.split(r'[+/]', re.sub(r'^[a-z]*\s*[OBAFGKMLTY]\s*'
                                     r'[\d.]*(?:[-/][\d.]+)?', '', text))[0]
    if re.match(r'^-?VI(?!I)', first):
        return 'subdwarf'
    if re.match(r'^-?IV', first):
        return 'subgiant'
    if re.match(r'^-?(III|II|I)(?![V])', first) or re.match(
            r'^-?I[ab]', first):
        return 'giant'
    return 'dwarf'


def _number(value: Any) -> Optional[float]:
    """a number of a table, None when there is none"""
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _clean(name: Any) -> str:
    """a name with single spaces (SIMBAD pads: 'L  356-105')"""
    return re.sub(r'\s+', ' ', str(name or '')).strip()


def name_key(name: Any) -> str:
    """a name as two catalogues may write it, to compare them: upper case,
    no space, no SIMBAD prefix (NAME, *, V*), GL and GLIESE as GJ"""
    text = re.sub(r'^(NAME|\*\*?|V\*)\s+', '', _clean(name))
    text = re.sub(r'[\s_]+', '', text).upper()
    return re.sub(r'^(GLIESE|GL)(?=\d)', 'GJ', text)


def sample(sptype: Sequence[str] = ('M0', 'M9'), dmax: float = 15.0,
           dec: Sequence[float] = (-90.0, 90.0),
           vmax: Optional[float] = None, dwarfs: bool = True,
           timeout: float = 180.0) -> List[Dict[str, Any]]:
    """
    The stars of SIMBAD within a distance, in a range of spectral types

    :param sptype: (str, str), the earliest and the latest type kept (M0
                   to M9: every M; M3 to M5.5), by sp_number
    :param dmax: float, the largest distance, from the parallax [pc]
    :param dec: (float, float), the declinations kept [deg] (a telescope
                does not see the whole sky)
    :param vmax: float or None, the faintest V kept (G when there is no
                 V); None for any, and a star with neither is kept
    :param dwarfs: bool, dwarfs only (sp_kind: not the giants, subgiants,
                   subdwarfs and white dwarfs)
    :param timeout: float [s]

    :return: list of dict, the nearest first: name (its GJ name when it
             has one, SIMBAD's otherwise), main (SIMBAD's), ra, dec [deg],
             distance [pc], plx and plx_err [mas], sptype, spnum, kind,
             otype, V, G, J, K, rv [km/s], ids ({catalogue: name}), names
             (every name kept) and near (another object of the sample
             within NEAR arcsec, or None)
    """
    from koloa.gather import SIMBAD_TAP, _tap
    low, high = sp_number(sptype[0]), sp_number(sptype[1])
    if low is None or high is None or not dmax > 0:
        raise ValueError(f'a range of spectral types (M0 to M9) and a '
                         f'distance: {sptype}, {dmax}')
    low, high = min(low, high), max(low, high)
    # the classes of the range, and their dwarfs and subdwarfs written
    #   with a prefix (dM3, sdM1)
    letters = CLASSES[int(low // 10):int(high // 10) + 1]
    like = ' OR '.join(f"b.sp_type LIKE '{pre}{letter}%'" for letter in letters
                       for pre in ('', 'd', 'sd', 'esd', 'usd', 'g'))
    where = (f'b.plx_value >= {1000.0 / dmax:.6f} AND ({like}) AND '
             f'b.dec >= {min(dec):.5f} AND b.dec <= {max(dec):.5f}')
    rows = _tap(SIMBAD_TAP, (
        f'SELECT TOP {MAXREC} b.oid, b.main_id, b.ra, b.dec, b.plx_value, '
        f'b.plx_err, b.sp_type, b.otype, b.rvz_radvel, f.V, f.G, f.J, f.K '
        f'FROM basic AS b LEFT JOIN allfluxes AS f ON f.oidref = b.oid '
        f'WHERE {where}'), timeout)
    cats = ' OR '.join(f"i.id LIKE '{cat} %'" if cat != 'TOI'
                       else "i.id LIKE 'TOI-%'" for cat in CATALOGUES)
    idents = _tap(SIMBAD_TAP, (
        f'SELECT TOP {10 * MAXREC} i.oidref, i.id FROM ident AS i JOIN '
        f'basic AS b ON b.oid = i.oidref WHERE {where} AND ({cats})'),
        timeout)
    names: Dict[str, List[str]] = {}
    for row in idents:
        names.setdefault(row['oidref'], []).append(_clean(row['id']))
    out = []
    for row in rows:
        num, kind = sp_number(row['sp_type']), sp_kind(row['sp_type'])
        plx = _number(row['plx_value'])
        if num is None or not low <= num <= high or not plx:
            continue
        if dwarfs and kind != 'dwarf':
            continue
        mags = {key: _number(row.get(key)) for key in ('V', 'G', 'J', 'K')}
        bright = mags['V'] if mags['V'] is not None else mags['G']
        if vmax is not None and bright is not None and bright > vmax:
            continue
        main = _clean(row['main_id'])
        mine = [main] + [name for name in names.get(row['oid'], [])
                         if name != main]
        ids = {}
        for cat in CATALOGUES:
            found = next((name for name in mine if name.startswith(
                'TOI-' if cat == 'TOI' else cat + ' ')), None)
            if found:
                ids[cat] = found
        out.append(dict(
            name=ids.get('GJ') or re.sub(r'^(NAME|\*\*?|V\*)\s+', '', main),
            main=main, ra=float(row['ra']), dec=float(row['dec']),
            distance=1000.0 / plx, plx=plx, plx_err=_number(row['plx_err']),
            sptype=_clean(row['sp_type']), spnum=num, kind=kind,
            otype=_clean(row['otype']), rv=_number(row['rvz_radvel']),
            ids=ids, names=mine, near=None, **mags))
    out.sort(key=lambda star: star['distance'])
    _neighbours(out)
    log(f'SIMBAD: {len(out)} stars from {sptype[0]} to {sptype[1]} within '
        f'{dmax:g} pc' + (' (dwarfs)' if dwarfs else ''), 'value')
    return out


def _neighbours(stars: List[Dict[str, Any]]) -> None:
    """each star told which other one of the sample is within NEAR arcsec
    (a system and its component, the same velocities for both)"""
    if len(stars) < 2:
        return
    ra = np.radians([star['ra'] for star in stars])
    dec = np.radians([star['dec'] for star in stars])
    for it, star in enumerate(stars):
        cosd = (np.sin(dec[it]) * np.sin(dec) + np.cos(dec[it]) * np.cos(dec)
                * np.cos(ra - ra[it]))
        sep = np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0))) * 3600.0
        sep[it] = np.inf
        best = int(np.argmin(sep))
        if sep[best] <= NEAR:
            star['near'] = stars[best]['name']


def ident_of(star: Dict[str, Any]) -> Dict[str, Any]:
    """a star of a sample as koloa.archive.resolve gives one (main, its
    names by catalogue, aliases, position): what koloa.dace.names and
    koloa.gather take"""
    ids = star.get('ids') or {}
    return dict(main=star['main'], name=star['name'], ra=star['ra'],
                dec=star['dec'], hd=ids.get('HD'), gj=ids.get('GJ'),
                hip=ids.get('HIP'), tic=ids.get('TIC'),
                aliases=list(star.get('names') or []))


def check_star(star: Dict[str, Any], root: str = 'archives',
               dace: bool = True, api_key: Any = None,
               refresh: bool = False) -> Dict[str, Any]:
    """
    What the archives have of one star of a sample, without gathering
    them: CARMENES DR1 and the surveys on VizieR by its position in their
    lists, DACE by asking it (its answer kept in the star's archives
    folder, where koloa.gather reads it back)

    :param star: dict, a star of sample()
    :param root: str, the folder of the archives
    :param dace: bool, ask DACE (the network, a few seconds a star)
    :param api_key: str, None or False: a DACE API key (koloa.dace)
    :param refresh: bool, ask DACE even when its answer is on disk

    :return: dict, set as star['archives'] too: dace (n, its instruments,
             or None; error when it could not be asked), carmenes (its
             velocities in DR1, or None), surveys (the keys of those that
             have the star), n (the velocities known so far)
    """
    from koloa import published as kpub
    from koloa.gather import carmenes_star, dace_rv, folder_name
    out: Dict[str, Any] = dict(dace=None, carmenes=None, surveys=[], n=0)
    try:
        found = carmenes_star(star['ra'], star['dec'])
        if found is not None:
            out['carmenes'] = int(float(found.get('nobs') or 0)) or 1
    except (OSError, ValueError) as err:
        out['carmenes_error'] = str(err)
    for survey in kpub.SURVEYS:
        try:
            if kpub.nearest_star(kpub.survey_stars(survey), star['ra'],
                                 star['dec']):
                out['surveys'].append(survey['key'])
        except (OSError, ValueError) as err:
            out['surveys_error'] = str(err)
    if dace:
        folder = os.path.join(root, folder_name(star['name']), 'rv', 'dace')
        try:
            data = dace_rv(ident_of(star), star['name'], folder,
                           api_key=api_key, refresh=refresh)
            if data is not None:
                out['dace'] = dict(n=int(data.n), instruments={
                    str(inst): int(np.sum(data.inst == inst))
                    for inst in data.instruments})
        except RuntimeError as err:
            out['dace_error'] = str(err).split(';')[0]
    out['n'] = ((out['dace'] or {}).get('n', 0) + (out['carmenes'] or 0))
    star['archives'] = out
    star['summary'] = spectrographs(out)
    return out


def spectrographs(archives: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    What the archives have of a star, spectrograph by spectrograph (the
    eras of an instrument and its releases put together): HARPS on DACE
    and in the RVBank is one line, with the velocities DACE has

    :param archives: dict, from check_star

    :return: list of dict, the most velocities first: name (the
             spectrograph), n (its velocities where an archive counted
             them: DACE, CARMENES DR1; None when only surveys have it,
             whose velocities are counted when they are gathered), where
             (the archives that have it)
    """
    from koloa.published import family
    found: Dict[str, Dict[str, Any]] = {}

    def add(name, where, num=None):
        one = found.setdefault(name, dict(name=name, n=None, where=[]))
        if where not in one['where']:
            one['where'].append(where)
        if num:
            one['n'] = (one['n'] or 0) + int(num)
    for inst, num in ((archives.get('dace') or {}).get('instruments')
                      or {}).items():
        base = re.sub(r'[_\d].*$', '', str(inst))
        add(family(inst) or DACE_NAMES.get(base, base), 'DACE', num)
    if archives.get('carmenes'):
        add('CARMENES', 'DR1', archives['carmenes'])
    for key in archives.get('surveys') or []:
        name, where = SURVEY_ARCHIVES.get(key, (key, key))
        add(name, where)
    return sorted(found.values(), key=lambda one: (-(one['n'] or 0),
                                                   one['name']))


def overview(stars: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """
    What the archives have of a sample, once its stars are checked: how
    many have velocities, and spectrograph by spectrograph

    :return: dict: n (stars), checked, data (those with velocities in an
             archive), none (checked, with nothing), files (those with
             files of one's own), nfiles, velocities (those counted),
             spectrographs (list of dict, the most stars first: name,
             stars, velocities, where: {archive: stars}), several (stars
             with two spectrographs or more) and single (with one)
    """
    checked = [star for star in stars if star.get('archives') is not None]
    table: Dict[str, Dict[str, Any]] = {}
    several = single = 0
    for star in checked:
        mine = star.get('summary') or []
        several += len(mine) > 1
        single += len(mine) == 1
        for one in mine:
            row = table.setdefault(one['name'], dict(
                name=one['name'], stars=0, velocities=0, where={}))
            row['stars'] += 1
            row['velocities'] += one['n'] or 0
            for where in one['where']:
                row['where'][where] = row['where'].get(where, 0) + 1
    data = sum(bool(star.get('summary')) for star in checked)
    return dict(
        n=len(stars), checked=len(checked), data=data,
        none=len(checked) - data,
        files=sum(bool(star.get('files')) for star in stars),
        nfiles=sum(len(star.get('files') or []) for star in stars),
        velocities=sum(row['velocities'] for row in table.values()),
        spectrographs=sorted(table.values(), key=lambda row: (
            -row['stars'], row['name'])),
        several=several, single=single)


def check(stars: Sequence[Dict[str, Any]], root: str = 'archives',
          dace: bool = True, api_key: Any = None, refresh: bool = False,
          workers: int = WORKERS,
          progress: Optional[Callable[[int, int], None]] = None
          ) -> List[Dict[str, Any]]:
    """
    What the archives have of each star of a sample (check_star), a few
    stars at a time

    :param progress: callable or None, told (done, total) after each star

    :return: list of dict, the stars, each with its 'archives'
    """
    from koloa import published as kpub
    from koloa.gather import carmenes_objects
    # the lists, once, before the stars are asked side by side
    try:
        carmenes_objects()
        for survey in kpub.SURVEYS:
            kpub.survey_stars(survey)
    except (OSError, ValueError) as err:
        log(f'survey: a list could not be fetched ({err})', 'warn')
    done = [0]

    def one(star):
        check_star(star, root, dace, api_key, refresh)
        done[0] += 1
        if progress is not None:
            progress(done[0], len(stars))
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
        list(pool.map(one, stars))
    have = sum(has_data(star) for star in stars)
    log(f'survey: {have} of {len(stars)} stars have velocities in the '
        f'archives', 'value')
    return list(stars)


def has_data(star: Dict[str, Any]) -> bool:
    """whether a star has velocities: in an archive checked, or in a file
    matched to it"""
    arch = star.get('archives') or {}
    return bool(arch.get('dace') or arch.get('carmenes')
                or arch.get('surveys') or star.get('files'))


def list_files(folders: Sequence[str], pattern: str = '*.rdb'
               ) -> List[str]:
    """the files of folders that match a pattern, or any of several
    (separated by spaces or commas), each file once"""
    found = set()
    for folder in folders:
        folder = os.path.expanduser(str(folder).strip())
        if not folder:
            continue
        if not os.path.isdir(folder):
            raise ValueError(f'no folder {folder}')
        for one in re.split(r'[\s,;]+', (pattern or '').strip()) or ['*.rdb']:
            found.update(path for path in glob.glob(os.path.join(
                folder, one or '*.rdb')) if os.path.isfile(path))
    return [os.path.abspath(path) for path in sorted(found)]


def match_files(stars: Sequence[Dict[str, Any]], folders: Sequence[str],
                pattern: str = '*.rdb', fetch: bool = True
                ) -> Dict[str, Any]:
    """
    The files of velocities of folders (LBL's of SPIRou, of NIRPS...) put
    with the stars of a sample: the star of each file is its OBJECT column,
    else its name, through APERO's names (koloa.apero_names), and it is a
    star of the sample when one of its names is one of that star's

    :param stars: list of dict, from sample(): each gets 'files' (its
                  paths, [] when none)
    :param folders: list of str, the folders looked in
    :param pattern: str, the files taken ('*.rdb', or several patterns)
    :param fetch: bool, a copy of APERO's names made when there is none

    :return: dict, matched ({star name: paths}), unmatched (the files
             whose star is not in the sample: path, the name read, the
             SIMBAD name APERO gives it)
    """
    from koloa import apero_names
    keys: Dict[str, Dict[str, Any]] = {}
    for star in stars:
        star['files'] = []
        for name in [star['name'], star['main']] + list(star.get('names')
                                                         or []):
            keys.setdefault(name_key(name), star)
    matched: Dict[str, List[str]] = {}
    unmatched = []
    for path in list_files(folders, pattern):
        try:
            found = apero_names.star_of_file(path, fetch)
        except (ImportError, OSError, ValueError):
            raw = apero_names.file_object(path) or next(iter(
                apero_names.name_candidates(path)), None)
            found = dict(raw=raw, target=raw or '', entry=None, apero=None)
        tries = [found.get('target'), found.get('raw'), found.get('apero')]
        entry = found.get('entry') or {}
        tries += list(entry.get('aliases') or [])
        star = next((keys[name_key(name)] for name in tries
                     if name and name_key(name) in keys), None)
        if star is None:
            unmatched.append(dict(path=path, raw=found.get('raw'),
                                  target=found.get('target')))
            continue
        star['files'].append(path)
        matched.setdefault(star['name'], []).append(path)
    log(f'survey: {sum(len(val) for val in matched.values())} files of '
        f'{len(matched)} stars of the sample; {len(unmatched)} files of '
        f'other stars', 'value')
    return dict(matched=matched, unmatched=unmatched)


# =============================================================================
# A batch of stars, here or on another machine
# =============================================================================
def _folder(name: str) -> str:
    """the folder of a star, as koloa.gather names it"""
    from koloa.gather import folder_name
    return folder_name(name)


def targets_of(root: str) -> Dict[str, Any]:
    """what a batch folder holds (its targets.json): name, options, and
    its targets (name, sptype, distance, files: paths from the folder)"""
    path = os.path.join(root, TARGETS)
    if not os.path.exists(path):
        raise ValueError(f'no {TARGETS} in {root}: not a batch folder (or '
                         f'ROOT, at the top of run_batch.py, is not where '
                         f'the folder is on this machine)')
    with open(path) as handle:
        return json.load(handle)


def result_path(root: str, name: str) -> str:
    """where the result of a star of a batch is kept"""
    return os.path.join(root, 'results', _folder(name), 'result.json')


def run_target(root: str, star: Dict[str, Any], rules: bool = True,
               gather: bool = False, trend: int = 1) -> Dict[str, Any]:
    """
    One star of a batch folder: the quick FIP of its files and of its
    archives (those of the folder; gathered when asked), its datasets
    chosen by koloa.datasets, a transit looked for at its best peak when
    its light curve is with its archives; kept in results/<star>/:
    result.json (its line of the table, what was done with each dataset),
    quick.json (the FIP, as the page draws it) and quicklook.pdf

    :param root: str, the batch folder
    :param star: dict, a target of targets.json
    :param rules: bool, False for the rules to leave no dataset out
    :param gather: bool, gather the archives the folder lacks (the network)
    :param trend: int, the order of the trend fitted (1: an acceleration)

    :return: dict, the result (also written)
    """
    from koloa import gui
    name = str(star['name'])
    out = os.path.dirname(result_path(root, name))
    os.makedirs(out, exist_ok=True)
    files = [path if os.path.isabs(path) else os.path.join(root, path)
             for path in star.get('files') or []]
    item = dict(path=files[0] if files else '', files=files, name=name,
                given=dict(name=name, sptype=star.get('sptype')),
                status='waiting', qid=None, error=None, summary=None,
                star=None, stage=None, note=None)
    bid = uuid.uuid4().hex[:8]
    start = time.time()
    gui.BATCHES[bid] = dict(
        id=bid, status='running', start=start, end=None, cancel=False,
        trend=int(trend), archives=True, regather=False, gather=bool(gather),
        rules=bool(rules), root=os.path.join(root, 'archives'), items=[item])
    # in this process, to its end (the page runs it in a thread)
    gui._run_batch(bid)
    state = gui.batch_state(bid)['items'][0]
    result = dict(name=name, sptype=star.get('sptype'),
                  distance=star.get('distance'), files=star.get('files') or [],
                  status=state['status'], error=state['error'],
                  note=state['note'], star=state['star'],
                  summary=state['summary'], datasets=state['datasets'],
                  elapsed=time.time() - start, made=time.strftime(
                      '%Y-%m-%d %H:%M:%S'))
    qid = item.get('qid')
    if qid and (gui.QUICKS.get(qid) or {}).get('result'):
        with open(os.path.join(out, 'quick.json'), 'w') as handle:
            json.dump(gui._finite(gui.quick_state(qid)), handle,
                      default=_plain)
        try:
            pdf = gui.quicklook_pdf(dict(
                files=[dict(path=path) for path in files], target=name,
                root=os.path.join(root, 'archives'), dace=True,
                carmenes=True, vizier=True, rules='on' if rules else 'off',
                exclude=', '.join(item.get('dace_copies') or []),
                trend=trend >= 1, curvature=trend >= 2), qid=qid)
            with open(os.path.join(out, 'quicklook.pdf'), 'wb') as handle:
                handle.write(pdf)
        except Exception as err:  # the numbers are kept all the same
            result['pdf_error'] = f'{type(err).__name__}: {err}'
    path = result_path(root, name)
    with open(path + '.part', 'w') as handle:
        json.dump(gui._finite(result), handle, indent=1, default=_plain)
    os.replace(path + '.part', path)
    return result


def _plain(value: Any) -> Any:
    """a numpy value for JSON"""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def results(root: str) -> List[Dict[str, Any]]:
    """the results of a batch folder so far, in the order of its targets
    (None for a star not done)"""
    out = []
    for star in targets_of(root)['targets']:
        path = result_path(root, star['name'])
        if os.path.exists(path):
            with open(path) as handle:
                out.append(json.load(handle))
        else:
            out.append(None)
    return out


def table(root: str) -> List[Dict[str, Any]]:
    """
    The table of a batch folder, one line per star done, written as
    results/table.csv and results/table.json: the star, its datasets used
    of those it has, the best peak of its quick FIP (its period, its FIP
    with and without its aliases, K, the rms), the acceleration of the
    star, the transit looked for

    :return: list of dict, the lines (COLUMNS)
    """
    import csv
    rows = []
    for res in results(root):
        if res is None:
            continue
        summ, sets = res.get('summary') or {}, res.get('datasets') or {}
        trans = summ.get('transit') or {}
        rows.append(dict(
            name=res['name'], sptype=res.get('sptype'),
            distance=res.get('distance'), status=res['status'],
            files=len(res.get('files') or []),
            datasets=(f'{sets.get("used")}/{sets.get("all")}' if sets
                      else ''),
            nights=summ.get('n'), baseline=summ.get('baseline'),
            period=summ.get('period'), fip=summ.get('fip'),
            fip_alone=summ.get('fip_alone'), K=summ.get('K'),
            K_err=summ.get('K_err'), rms=summ.get('rms'),
            accel=summ.get('accel'), accel_err=summ.get('accel_err'),
            accel_sigma=summ.get('accel_sigma'),
            transit=trans.get('status'), transit_snr=trans.get('snr'),
            elapsed=res.get('elapsed'), error=res.get('error')))
    os.makedirs(os.path.join(root, 'results'), exist_ok=True)
    with open(os.path.join(root, 'results', 'table.csv'), 'w',
              newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ('' if val is None else f'{val:.6g}'
                                   if isinstance(val, float) else val)
                             for key, val in row.items()})
    with open(os.path.join(root, 'results', 'table.json'), 'w') as handle:
        json.dump(rows, handle, indent=1)
    return rows


def run(root: str, jobs: int = 1, rules: bool = True, gather: bool = False,
        only: Optional[Sequence[str]] = None, again: bool = False,
        part: Optional[Tuple[int, int]] = None, trend: int = 1
        ) -> List[Dict[str, Any]]:
    """
    The batch of a folder: each star not done yet (run_target), a few at
    once

    Each star takes a Python of its own kind: with jobs above 1 this one
    starts that many more (python -m koloa.survey ROOT --part i/n), each
    with one star in n and its log in logs/, and waits for them. A job
    array of a cluster asks for its part itself.

    :param root: str, the batch folder
    :param jobs: int, the stars at once (about three cores each)
    :param rules: bool, False for the rules to leave no dataset out
    :param gather: bool, gather the archives the folder lacks (the network)
    :param only: list of str or None, these stars only
    :param again: bool, the stars already done too
    :param part: (int, int) or None, the stars i, i + n, i + 2n... only
    :param trend: int, the order of the trend fitted

    :return: list of dict, the table (table())
    """
    root = os.path.abspath(os.path.expanduser(root))
    stars = targets_of(root)['targets']
    if only:
        keep = {name_key(name) for name in only}
        stars = [star for star in stars if name_key(star['name']) in keep]
    os.makedirs(os.path.join(root, 'logs'), exist_ok=True)
    jobs = max(1, int(jobs))
    if part is None and jobs > 1 and len(stars) > 1:
        jobs = min(jobs, len(stars))
        # koloa as this process has it, for the others
        src = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = dict(os.environ, PYTHONPATH=src + os.pathsep + os.environ.get(
            'PYTHONPATH', ''))
        base = [sys.executable, '-m', 'koloa.survey', root,
                '--trend', str(trend)]
        base += ([] if rules else ['--no-rules']) + (
            ['--gather'] if gather else []) + (['--again'] if again else [])
        if only:
            base += ['--only'] + [str(name) for name in only]
        procs = []
        for it in range(jobs):
            handle = open(os.path.join(root, 'logs', f'part_{it}.log'), 'a')
            procs.append((subprocess.Popen(
                base + ['--part', f'{it}/{jobs}'], env=env, stdout=handle,
                stderr=subprocess.STDOUT), handle))
        log(f'batch: {len(stars)} stars, {jobs} at once (their logs in '
            f'{os.path.join(root, "logs")})', 'info')
        seen = -1
        while any(proc.poll() is None for proc, _ in procs):
            done = sum(os.path.exists(result_path(root, star['name']))
                       for star in stars)
            if done != seen:
                seen = done
                log(f'batch: {done} of {len(stars)} stars done', 'value')
            time.sleep(5)
        for proc, handle in procs:
            handle.close()
            if proc.returncode:
                log(f'batch: a part ended with code {proc.returncode} (see '
                    f'its log)', 'warn')
    else:
        mine = stars if part is None else stars[part[0]::part[1]]
        for it, star in enumerate(mine):
            if os.path.exists(result_path(root, star['name'])) and not again:
                continue
            log(f'batch: {star["name"]} ({it + 1} of {len(mine)})', 'info')
            try:
                res = run_target(root, star, rules, gather, trend)
                summ = res.get('summary') or {}
                log(f'batch: {star["name"]} {res["status"]} in '
                    f'{res["elapsed"]:.0f} s' + (
                        f', P = {summ["period"]:.4f} d, FIP '
                        f'{summ["fip"]:.2g}' if summ.get('period') else '')
                    + (f' ({res["error"]})' if res.get('error') else ''),
                    'value')
            except Exception as err:  # the next star all the same
                log(f'batch: {star["name"]} failed '
                    f'({type(err).__name__}: {err})', 'warn')
    rows = table(root)
    log(f'batch: {len(rows)} of {len(targets_of(root)["targets"])} stars in '
        f'{os.path.join(root, "results", "table.csv")}', 'value')
    return rows


SCRIPT = '''#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa: the batch {name}, made {made}

{nstar} stars, each with its files (files/) and its archives (archives/):
the quick FIP of each, its datasets chosen by koloa.datasets, kept in
results/ as it goes (a batch stopped goes on where it was).

    python run_batch.py                 # every star not done yet
    python run_batch.py --part 3/20     # the stars 3, 23, 43... (a job array)
    python run_batch.py --root /elsewhere/{name}

Everything is found from ROOT: nothing else of this script names a path.
"""
import os
import sys

# =============================================================================
# The settings
# =============================================================================
# where this folder is ON THE MACHINE THAT RUNS THE BATCH. Change it when
#   the folder is copied elsewhere (a server, a cluster): files/, archives/,
#   cache/, koloa_src/ and results/ are found from it.
ROOT = {root!r}

# the stars at once (each takes about three cores)
JOBS = {jobs}

# the datasets of each star chosen by koloa.datasets (the best release of
#   the same spectra, not those that constrain nothing); False: all of them
RULES = {rules}

# ask the archives for what archives/ lacks (needs the network: not on the
#   nodes of most clusters). False: what is in archives/, as packed.
GATHER = {gather}

# these stars only ([]: all of them), and those already done again
ONLY = []
AGAIN = False

# the trend fitted with the signals: 1 an acceleration, 2 its change too
TREND = {trend}


# =============================================================================
# The batch (under main: koloa computes in processes of its own, which read
#   this file again as they start and must not start the batch themselves)
# =============================================================================
def main():
    root, part = ROOT, None
    args = sys.argv[1:]
    if '--root' in args:
        root = args[args.index('--root') + 1]
    if '--part' in args:
        one, of = args[args.index('--part') + 1].split('/')
        part = (int(one), int(of))
    root = os.path.abspath(os.path.expanduser(root))
    if not os.path.exists(os.path.join(root, 'targets.json')):
        here = os.path.dirname(os.path.abspath(__file__))
        raise SystemExit(
            'ROOT is ' + root + ': there is no targets.json there.\\n'
            'Set ROOT, at the top of this script, to where the folder is '
            'on this machine'
            + (' (this script is in ' + here + ').' if os.path.exists(
                os.path.join(here, 'targets.json')) else '.'))
    # what koloa fetched once, and koloa itself, as they were packed
    os.environ['KOLOA_CACHE'] = os.path.join(root, 'cache')
    src = os.path.join(root, 'koloa_src')
    if os.path.isdir(src):
        sys.path.insert(0, src)
        os.environ['PYTHONPATH'] = src + os.pathsep + os.environ.get(
            'PYTHONPATH', '')
    from koloa import survey
    survey.run(root, jobs=JOBS, rules=RULES, gather=GATHER, only=ONLY,
               again=AGAIN, part=part, trend=TREND)


if __name__ == '__main__':
    main()
'''

SUBMIT = '''#!/bin/bash
# koloa: the batch {name} as a job array of SLURM (the clusters of the
# Alliance, formerly Compute Canada): one star per task, {nstar} tasks, at
# most {jobs} at once. Set your account and what loads Python, then:
#     sbatch submit.sh
#SBATCH --job-name=koloa_{name}
#SBATCH --account=def-CHANGE_ME
#SBATCH --array=0-{last}%{jobs}
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/slurm_%A_%a.log

module load python scipy-stack
cd {root}
python run_batch.py --part ${{SLURM_ARRAY_TASK_ID}}/{nstar}
'''

README = '''koloa: the batch {name}, made {made}

{nstar} stars. In this folder:

  run_batch.py    the script: its first setting, ROOT, is where this folder
                  is on the machine that runs it ({root} as packed)
  submit.sh       the same as a SLURM job array (a star per task)
  targets.json    the stars, and the files of each
  files/          the files of velocities of each star
  archives/       the archives of each star, gathered {when}
  cache/          what koloa fetched once (the NASA Exoplanet Archive, the
                  lists of the surveys, APERO's names)
  koloa_src/      koloa as it was when the batch was packed
  results/        made by the batch: a folder per star (result.json,
                  quick.json, quicklook.pdf) and table.csv
  logs/           made by the batch

Elsewhere (a server, a cluster):

  tar xzf {name}.tar.gz
  cd {name}
  # set ROOT at the top of run_batch.py to this folder, then
  python run_batch.py

Python 3.9 or later with numpy, scipy and matplotlib; nothing is asked of
the network (GATHER = False). A batch stopped goes on where it was: each
star done is kept. Bring results/ back to look at them.
'''


def pack(stars: Sequence[Dict[str, Any]], name: str, out: str = '.',
         root: str = 'archives', server_root: Optional[str] = None,
         gather: bool = True, tess: bool = False, jobs: int = 6,
         rules: bool = True, trend: int = 1, refresh: bool = False,
         api_key: Any = None, tar: bool = True,
         progress: Optional[Callable[[int, int, str], None]] = None,
         workers: int = 3) -> Dict[str, Any]:
    """
    A batch folder, and its .tar.gz: everything the batch of some stars
    needs, to run it here or on another machine (see the module)

    :param stars: list of dict, the stars: name, and files (their paths
                  here), sptype, distance (from sample() and match_files())
    :param name: str, the name of the batch: its folder, its tar
    :param out: str, where the folder is made
    :param root: str, the folder of the archives here (each star's are
                 copied from it)
    :param server_root: str or None, where the folder will be on the
                        machine that runs it: the ROOT of run_batch.py
                        (None: where it is made)
    :param gather: bool, gather first the archives of the stars that have
                   none here (the network: DACE, CARMENES DR1, VizieR)
    :param tess: bool, their light curves too (larger; for the transit
                 looked for at the best peak)
    :param jobs: int, the stars at once, as the script is written
    :param rules: bool, the RULES of the script
    :param trend: int, its TREND
    :param refresh: bool, gather again the archives already here
    :param api_key: str, None or False: a DACE API key (koloa.dace)
    :param tar: bool, make the .tar.gz too
    :param progress: callable or None, told (done, total, star)
    :param workers: int, the stars whose archives are gathered at once
                    (the tables of a star's papers take a minute or two)

    :return: dict, folder, tar (None without), n (stars), files, size (of
             the tar [bytes]), and missing (the stars with nothing: no
             file and no archive)
    """
    from koloa import __version__
    from koloa.gather import gather as gather_star
    from koloa.paths import cache
    name = re.sub(r'[^A-Za-z0-9+\-.]+', '_', str(name).strip()).strip('_')
    if not name:
        raise ValueError('a name for the batch')
    folder = os.path.join(os.path.abspath(os.path.expanduser(out)), name)
    os.makedirs(folder, exist_ok=True)
    targets, missing, nfile = [], [], 0
    # the archives first, a few stars at once: those that have none here
    todo = [star for star in stars if gather and (
        refresh or not os.path.exists(os.path.join(
            root, _folder(star['name']), 'manifest.json')))]
    done = [0]

    def fetch(star):
        try:
            gather_star(star['name'], root, dace=True, carmenes=True,
                        tess=tess, vizier=True, refresh=refresh,
                        api_key=api_key)
        except Exception as err:  # the star with what there is
            log(f'pack: the archives of {star["name"]} not gathered '
                f'({type(err).__name__}: {err})', 'warn')
        done[0] += 1
        if progress is not None:
            progress(done[0], len(todo), star['name'])
    if todo:
        # what every star asks for once (the tables of the NASA Exoplanet
        #   Archive, the lists of the surveys), before the stars are
        #   gathered side by side
        try:
            from koloa import archive, published as kpub
            from koloa.gather import carmenes_objects
            archive.tables()
            carmenes_objects()
            for one in kpub.SURVEYS:
                kpub.survey_stars(one)
        except Exception as err:  # each star asks again, and says
            log(f'pack: a table could not be fetched first ({err})', 'warn')
        with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
            list(pool.map(fetch, todo))
    for it, star in enumerate(stars):
        tag = _folder(star['name'])
        if progress is not None and not todo:
            progress(it, len(stars), star['name'])
        here = os.path.join(root, tag)
        if os.path.isdir(here):
            shutil.copytree(here, os.path.join(folder, 'archives', tag),
                            dirs_exist_ok=True)
        files = []
        for path in star.get('files') or []:
            dest = os.path.join(folder, 'files', tag, os.path.basename(path))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copy2(os.path.expanduser(path), dest)
            files.append(os.path.join('files', tag, os.path.basename(path)))
        nfile += len(files)
        if not files and not os.path.isdir(here):
            missing.append(star['name'])
        targets.append(dict(
            name=star['name'], main=star.get('main'),
            sptype=star.get('sptype'), distance=star.get('distance'),
            ra=star.get('ra'), dec=star.get('dec'), files=files,
            archives=star.get('archives')))
    # what koloa fetched once, and koloa itself
    for part in CACHE_PARTS:
        src = cache(part)
        dest = os.path.join(folder, 'cache', part)
        if os.path.isdir(src):
            shutil.copytree(src, dest, dirs_exist_ok=True)
        elif os.path.exists(src):
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copy2(src, dest)
    package = os.path.dirname(os.path.abspath(__file__))
    shutil.copytree(package, os.path.join(folder, 'koloa_src', 'koloa'),
                    dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    made = time.strftime('%Y-%m-%d %H:%M')
    words = dict(name=name, made=made, nstar=len(targets),
                 root=server_root or folder, jobs=int(jobs),
                 rules=bool(rules), gather=False, trend=int(trend),
                 last=max(len(targets) - 1, 0),
                 when=made if gather else 'before')
    with open(os.path.join(folder, TARGETS), 'w') as handle:
        json.dump(dict(name=name, made=made, koloa=__version__,
                       options=dict(rules=bool(rules), trend=int(trend)),
                       targets=targets), handle, indent=1, default=_plain)
    for file, text in (('run_batch.py', SCRIPT), ('submit.sh', SUBMIT),
                       ('README.txt', README)):
        with open(os.path.join(folder, file), 'w') as handle:
            handle.write(text.format(**words))
    made_tar = None
    if tar:
        made_tar = folder + '.tar.gz'
        with tarfile.open(made_tar, 'w:gz') as handle:
            handle.add(folder, arcname=name, filter=lambda info: (
                None if os.path.basename(info.name) in ('results', 'logs')
                and info.isdir() else info))
    if progress is not None:
        progress(len(stars), len(stars), '')
    log(f'pack: {len(targets)} stars, {nfile} files in {folder}'
        + (f', {os.path.getsize(made_tar) / 1e6:.1f} MB in {made_tar}'
           if made_tar else ''), 'value')
    return dict(folder=folder, tar=made_tar, n=len(targets), files=nfile,
                size=os.path.getsize(made_tar) if made_tar else None,
                missing=missing, root=words['root'])


# =============================================================================
# Start of code
# =============================================================================
def main(argv: Optional[Sequence[str]] = None) -> None:
    """python -m koloa.survey ROOT: the batch of a folder (run)"""
    import argparse
    parser = argparse.ArgumentParser(
        prog='python -m koloa.survey',
        description='The batch of a folder made by koloa.survey.pack: the '
                    'quick FIP of each of its stars not done yet.')
    parser.add_argument('root', help='the batch folder (its targets.json)')
    parser.add_argument('--jobs', type=int, default=1,
                        help='the stars at once (about three cores each)')
    parser.add_argument('--part', default=None, metavar='I/N',
                        help='the stars I, I + N, I + 2N... only (a task '
                             'of a job array)')
    parser.add_argument('--only', nargs='+', default=None, metavar='STAR')
    parser.add_argument('--again', action='store_true',
                        help='the stars already done too')
    parser.add_argument('--no-rules', action='store_true',
                        help='every dataset used (koloa.datasets leaves '
                             'none out)')
    parser.add_argument('--gather', action='store_true',
                        help='gather the archives the folder lacks')
    parser.add_argument('--trend', type=int, default=1)
    args = parser.parse_args(argv)
    part = None
    if args.part:
        one, of = args.part.split('/')
        part = (int(one), int(of))
    run(args.root, jobs=args.jobs, rules=not args.no_rules,
        gather=args.gather, only=args.only, again=args.again, part=part,
        trend=args.trend)


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

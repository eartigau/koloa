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

The folder is run outside koloa, and its README.txt says how, for a person
or for a Claude session on the machine that runs it: the copy, what it
needs of Python, the launch, a job array, the results brought back, each
column of the table. Two more of its script, and here:

    survey.status('m_dwarfs_15pc')    # nothing computed: is all there, and
                                      #   where is each star (--check)
    survey.stop('m_dwarfs_15pc')      # the processes of this batch on this
                                      #   machine, and no other (--stop)

It asks nothing of the network unless told to gather: a light curve that
is not with the archives of a star is not fetched.

Each star done has a page of summary (koloa.batchpdf: how its quick look
reads, a candidate, a known planet, a drift, nothing), the first of its
PDF, and the batch one PDF of them all (results/summary.pdf). With
report=True (REPORT in the script), a star that has a candidate has its
detailed report too (report_star: koloa --detailed on the series its
quick look used, whatever it came from; an SHO GP at its published
rotation period, a local GP by period band without one).

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
import signal
import socket
import subprocess
import sys
import tarfile
import textwrap
import threading
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
#: the mark of a star being computed, in its folder of the results
RUNNING = 'running.json'
#: what the batch needs of this Python
NEEDS = ('numpy', 'scipy', 'matplotlib')
#: the cores a star takes, about (its chains, and what starts them), and
#: with its detailed report (the fits of its GP check, side by side)
CORES = 3
CORES_REPORT = 4
#: the share of the free cores of a machine that a batch takes: the rest is
#: left to the others, and to what the batch itself needs beside its stars
SHARE = 2.0 / 3.0
#: the columns of the table of a batch
COLUMNS = ('name', 'sptype', 'distance', 'status', 'files', 'datasets',
           'nights', 'baseline', 'period', 'fip', 'fip_alone', 'K', 'K_err',
           'rms', 'accel', 'accel_err', 'accel_sigma', 'transit',
           'transit_snr', 'verdict', 'report', 'elapsed', 'error')
#: the FIP below which a peak counts: a star with one that is neither a
#: known planet nor a drift has a candidate, and a detailed report when
#: the batch makes them
REPORT_FIP = 0.01


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
             (every name kept), near (another object of the sample within
             NEAR arcsec, or None) and rotation (the rotation periods
             SIMBAD lists, the latest first: period [days] and source;
             check() adds that of CARMENES DR1)
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
    # the rotation periods SIMBAD lists (its table of variability: the
    #   measurements of type ROT, each with its paper), the latest first
    spins: Dict[str, List[Dict[str, Any]]] = {}
    try:
        for row in _tap(SIMBAD_TAP, (
                f'SELECT TOP {10 * MAXREC} v.oidref, v.period, v.bibcode '
                f'FROM mesVar AS v JOIN basic AS b ON b.oid = v.oidref '
                f"WHERE {where} AND v.vartyp = 'ROT' AND v.period IS NOT "
                f'NULL'), timeout):
            period = _number(row['period'])
            one = dict(period=period, source=_clean(row['bibcode']))
            if period and one not in spins.setdefault(row['oidref'], []):
                spins[row['oidref']].append(one)
    except (OSError, ValueError) as err:  # the sample without them
        log(f'SIMBAD: its rotation periods could not be asked ({err})',
            'warn')
    for found in spins.values():
        found.sort(key=lambda one: one['source'][:4], reverse=True)
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
            ids=ids, names=mine, near=None,
            rotation=list(spins.get(row['oid'], [])), **mags))
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
               refresh: bool = False, carmenes: bool = True,
               vizier: bool = True) -> Dict[str, Any]:
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
    :param carmenes: bool, look the star up in CARMENES DR1
    :param vizier: bool, and in the lists of the surveys on VizieR

    :return: dict, set as star['archives'] too: dace (n, its instruments,
             or None; error when it could not be asked), carmenes (its
             velocities in DR1, or None), surveys (the keys of those that
             have the star), n (the velocities known so far)
    """
    from koloa import published as kpub
    from koloa.gather import carmenes_star, dace_rv, folder_name
    out: Dict[str, Any] = dict(dace=None, carmenes=None, surveys=[], n=0)
    try:
        found = carmenes_star(star['ra'], star['dec']) if carmenes else None
        if found is not None:
            out['carmenes'] = int(float(found.get('nobs') or 0)) or 1
            # its rotation period, beside those SIMBAD lists
            spin = _number(found.get('p_rot'))
            mine = star.setdefault('rotation', [])
            if spin and not any(one['source'].startswith('CARMENES')
                                for one in mine):
                mine.append(dict(period=spin, source='CARMENES DR1' + (
                    f' ({_clean(found.get("p_rot_source"))})'
                    if _clean(found.get('p_rot_source')) else '')))
    except (OSError, ValueError) as err:
        out['carmenes_error'] = str(err)
    for survey in kpub.SURVEYS if vizier else ():
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
            # a star DACE knows with no velocity of it has nothing there
            if data is not None and data.n:
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
          progress: Optional[Callable[[int, int], None]] = None,
          carmenes: bool = True, vizier: bool = True
          ) -> List[Dict[str, Any]]:
    """
    What the archives have of each star of a sample (check_star), a few
    stars at a time

    :param progress: callable or None, told (done, total) after each star
    :param carmenes: bool, CARMENES DR1 among them
    :param vizier: bool, and the surveys on VizieR

    :return: list of dict, the stars, each with its 'archives'
    """
    from koloa import published as kpub
    from koloa.gather import carmenes_objects
    # the lists, once, before the stars are asked side by side
    try:
        if carmenes:
            carmenes_objects()
        for survey in kpub.SURVEYS if vizier else ():
            kpub.survey_stars(survey)
    except (OSError, ValueError) as err:
        log(f'survey: a list could not be fetched ({err})', 'warn')
    done = [0]

    def one(star):
        check_star(star, root, dace, api_key, refresh, carmenes, vizier)
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


def report_star(out: str, data: Any, name: str,
                rotation: Optional[Dict[str, Any]] = None, trend: int = 1,
                network: bool = False) -> Dict[str, Any]:
    """
    The detailed report of a star (koloa --detailed) on the series its
    quick look used, whatever it came from (files, archives, or both):
    the series is written beside the report (velocities.csv, each dataset
    its instrument), and the report is run on that file by a process of
    its own, its log beside it (run.log)

    The GP of its FIP is an SHO at the rotation period of the star when
    one is published (the one given: SIMBAD's or CARMENES DR1's; else the
    NASA Exoplanet Archive's), a local GP by period band otherwise.

    :param out: str, the folder of the report
    :param data: RVData or None, the series (None: the velocities.csv
                 already there, of a report that was stopped)
    :param name: str, the SIMBAD name of the star
    :param rotation: dict or None, a published rotation: period [days]
                     and source
    :param trend: int, the order of the trend fitted (1: an acceleration)
    :param network: bool, let the report ask for the light curves of TESS
                    (False: nothing is asked of the network)

    :return: dict, status (done, failed), gp (the GP of its FIP, in
             words), signals (each signal it fitted, with its verdict),
             pdf and text (its report, None when not written), folder,
             command, elapsed, error
    """
    from koloa.gather import write_rv
    start = time.time()
    os.makedirs(out, exist_ok=True)
    series = os.path.join(out, 'velocities.csv')
    if data is not None:
        write_rv(data, series)
    # its files named after the star, as those of any report
    args = [series, '--detailed', '--target', name, '--name', name,
            '--outdir', out, '--no-rules']
    if rotation and rotation.get('period'):
        args += ['--rotation', f'{float(rotation["period"]):g}']
    else:
        args += ['--fip-gp', 'sho']  # the archive's rotation, else by band
    if trend >= 2:
        args.append('--curvature')
    elif trend < 1:
        args.append('--no-trend')
    if not network:
        args.append('--no-tess')
    # koloa as this process has it
    src = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env = dict(os.environ, PYTHONPATH=src + os.pathsep + os.environ.get(
        'PYTHONPATH', ''))
    told = dict(status='failed', gp=None, signals=[], pdf=None, text=None,
                folder=out, command='koloa ' + ' '.join(
                    f"'{arg}'" if ' ' in arg else arg for arg in args),
                error=None)
    with open(os.path.join(out, 'run.log'), 'w') as handle:
        code = subprocess.run([sys.executable, '-m', 'koloa.cli'] + args,
                              env=env, stdout=handle,
                              stderr=subprocess.STDOUT).returncode
    told['elapsed'] = time.time() - start
    for key, pattern in (('pdf', '*_report.pdf'), ('text', '*_report.txt')):
        found = sorted(glob.glob(os.path.join(out, pattern)))
        told[key] = found[0] if found else None
    kept = sorted(glob.glob(os.path.join(out, '*_summary.json')))
    if code != 0 or not kept:
        last = ''
        with open(os.path.join(out, 'run.log'), errors='replace') as handle:
            lines = [line.strip() for line in handle if line.strip()]
        if lines:
            last = lines[-1].split(' | ', 1)[-1]
        told['error'] = f'ended with code {code}' + (f': {last}' if last
                                                     else '')
        return told
    with open(kept[0]) as handle:
        summ = json.load(handle)
    spin = ((summ.get('known') or {}).get('star') or {}).get('rotation')
    if rotation and rotation.get('period'):
        told['gp'] = (f'SHO at the rotation, {float(rotation["period"]):g} d'
                      + (f' ({rotation["source"]})' if rotation.get('source')
                         else ''))
    elif spin:
        told['gp'] = (f'SHO at the rotation, {float(spin):g} d (NASA '
                      f'Exoplanet Archive)')
    else:
        told['gp'] = 'local, by period band (no rotation published)'
    ducks = summ.get('duck') or {}
    for orbit in summ.get('orbits') or []:
        per, amp = orbit['P'][0], orbit['K'][0]
        verdict = next((val for key, val in ducks.items()
                        if abs(float(key) / per - 1) < 1e-3), '')
        told['signals'].append(
            f'{per:.4f} d, K = {amp:.2f} m/s'
            + (f', m sin i = {orbit["msini"][0]:.1f} ME' if orbit.get('msini')
               and np.isfinite(orbit['msini'][0]) else '')
            + (f': {verdict}' if verdict else ''))
    told['status'] = 'done'
    return told


def _report_of(root: str, star: Dict[str, Any], result: Dict[str, Any],
               data: Any, gather: bool, trend: int) -> None:
    """the detailed report of a star of a batch folder, in
    results/<star>/report/, and what came of it in its result (kept)"""
    name = str(star['name'])
    out = os.path.join(os.path.dirname(result_path(root, name)), 'report')
    spins = [one for one in star.get('rotation') or [] if one.get('period')]
    log(f'batch: {name}, a candidate: its detailed report (GP: '
        + (f'SHO at the rotation, {float(spins[0]["period"]):g} d' if spins
           else 'SHO at the rotation of the archive, else by band')
        + f'), its log in {os.path.join(out, "run.log")}', 'info')
    try:
        told = report_star(out, data, name, spins[0] if spins else None,
                           trend, network=gather)
    except Exception as err:  # the quick look is kept all the same
        told = dict(status='failed', error=f'{type(err).__name__}: {err}',
                    signals=[], gp=None, pdf=None, text=None, folder=out)
    # as the folder holds them: no path of this machine
    for key in ('pdf', 'text', 'folder'):
        if told.get(key):
            told[key] = os.path.relpath(told[key], root)
    if told.get('command'):
        told['command'] = told['command'].replace(root + os.sep, '')
    told['where'] = told.get('pdf') or told.get('text') or told.get('folder')
    result['report'] = told
    log(f'batch: {name}, its detailed report {told["status"]}'
        + (f' in {told.get("elapsed", 0):.0f} s' if told.get('elapsed')
           else '') + (f' ({told["error"]})' if told.get('error') else ''),
        'value' if told['status'] == 'done' else 'warn')


def _keep(root: str, result: Dict[str, Any]) -> None:
    """the result of a star of a batch folder, written whole"""
    from koloa import gui
    path = result_path(root, str(result['name']))
    with open(path + '.part', 'w') as handle:
        json.dump(gui._finite(result), handle, indent=1, default=_plain)
    os.replace(path + '.part', path)


def run_target(root: str, star: Dict[str, Any], rules: bool = True,
               gather: bool = False, trend: int = 1, report: bool = False,
               report_fip: float = REPORT_FIP) -> Dict[str, Any]:
    """
    One star of a batch folder: the quick FIP of its files and of its
    archives (those of the folder; gathered when asked), its datasets
    chosen by koloa.datasets, a transit looked for at its best peak when
    its light curve is with its archives (fetched only with gather); kept
    in results/<star>/: result.json (its line of the table, what was done
    with each dataset, how its peaks read: koloa.batchpdf.reading),
    quick.json (the FIP, as the page draws it) and quicklook.pdf, whose
    first page is the summary of the star (koloa.batchpdf.figure).

    With report, a star that has a candidate (a peak with a FIP below
    report_fip that is neither a known planet nor a drift) has its
    detailed report too, in results/<star>/report/ (report_star), on the
    same series, whatever it came from.

    :param root: str, the batch folder
    :param star: dict, a target of targets.json
    :param rules: bool, False for the rules to leave no dataset out
    :param gather: bool, gather the archives the folder lacks (the network)
    :param trend: int, the order of the trend fitted (1: an acceleration)
    :param report: bool, the detailed report of a star with a candidate
    :param report_fip: float, the FIP below which a peak counts

    :return: dict, the result (also written)
    """
    from koloa import batchpdf, gui
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
        # the network only when asked for (GATHER): a light curve that is
        #   not with the archives of the star is not fetched
        fetch=bool(gather), keep=True,
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
    quick = None
    if qid and (gui.QUICKS.get(qid) or {}).get('result'):
        quick = gui._finite(gui.quick_state(qid))
        with open(os.path.join(out, 'quick.json'), 'w') as handle:
            json.dump(quick, handle, default=_plain)
    # how its peaks read: a known planet, the rotation, a drift, a candidate
    result['reading'] = batchpdf.reading(result, quick, star, report_fip)
    result['report'] = dict(status='none', why=(
        'not asked for' if not report else 'no candidate'))
    if report and result['reading']['kind'] == 'candidate' \
            and item.get('data') is not None:
        # kept before it starts: a report stopped (the end of a job) is
        #   taken again on its own, the quick look is not computed again
        result['report'] = dict(status='running')
        _keep(root, result)
        _report_of(root, star, result, item['data'], gather, trend)
        result['elapsed'] = time.time() - start
    if quick is not None:
        try:
            front = batchpdf.figure(result, quick, star,
                                    batch=str(targets_of(root).get('name')
                                              or ''), limit=report_fip)
            pdf = gui.quicklook_pdf(dict(
                files=[dict(path=path) for path in files], target=name,
                root=os.path.join(root, 'archives'), dace=True,
                carmenes=True, vizier=True, rules='on' if rules else 'off',
                exclude=', '.join(item.get('dace_copies') or []),
                trend=trend >= 1, curvature=trend >= 2), qid=qid,
                front=front)
            with open(os.path.join(out, 'quicklook.pdf'), 'wb') as handle:
                handle.write(pdf)
        except Exception as err:  # the numbers are kept all the same
            result['pdf_error'] = f'{type(err).__name__}: {err}'
    _keep(root, result)
    return result


def resume_report(root: str, star: Dict[str, Any], gather: bool = False,
                  trend: int = 1) -> Optional[Dict[str, Any]]:
    """
    The detailed report of a star whose report was stopped (the end of a
    job, a kill) taken again on the series kept beside it; its quick look
    is not computed again

    :return: dict or None, the result (None when no report was left
             unfinished)
    """
    path = result_path(root, str(star['name']))
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        result = json.load(handle)
    if (result.get('report') or {}).get('status') != 'running':
        return None
    series = os.path.join(os.path.dirname(path), 'report', 'velocities.csv')
    if not os.path.exists(series):
        result['report'] = dict(status='failed', error='stopped before its '
                                'series was written', signals=[])
    else:
        _report_of(root, star, result, None, gather, trend)
    _keep(root, result)
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
            verdict=(res.get('reading') or {}).get('kind'),
            report=(res.get('report') or {}).get('status'),
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


def _fip(value: Any) -> str:
    """a FIP as it is said: one of 0 is below what a number holds"""
    return '< 1e-300' if not value else f'{value:.2g}'


def _mark(path: str, **more: Any) -> None:
    """a process of a batch says it is there: its number, its machine,
    since when"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        json.dump(dict(pid=os.getpid(), host=socket.gethostname(),
                       since=time.time(), **more), handle)


def _unmark(path: str) -> None:
    """a mark taken away (it may be gone already)"""
    try:
        os.remove(path)
    except OSError:
        pass


def _marked(path: str) -> Optional[Dict[str, Any]]:
    """a mark, and whether its process is still there: alive is True,
    False, or None when the mark is of another machine (a node of a
    cluster: not known from here)"""
    try:
        with open(path) as handle:
            mark = json.load(handle)
        pid = int(mark['pid'])
    except (OSError, ValueError, KeyError, TypeError):
        return None
    if mark.get('host') != socket.gethostname():
        return dict(mark, alive=None)
    try:
        os.kill(pid, 0)
        alive = True
    except PermissionError:  # there, and another's
        alive = True
    except OSError:
        alive = False
    return dict(mark, alive=alive)


def _ours(pid: int) -> bool:
    """whether a process is one of a batch, by its command line: a number
    kept from before a restart of the machine may be another's by now"""
    try:
        told = subprocess.run(['ps', '-p', str(pid), '-o', 'command='],
                              capture_output=True, text=True,
                              timeout=10).stdout
    except (OSError, subprocess.SubprocessError):
        return False
    return 'koloa.survey' in told or 'run_batch.py' in told


def _end(pid: int) -> bool:
    """a process of a batch ended, with those it started when it leads
    them (a part does, and so does a batch started from a shell)"""
    try:
        if hasattr(os, 'killpg') and os.getpgid(pid) == pid:
            os.killpg(pid, signal.SIGTERM)
        else:
            os.kill(pid, signal.SIGTERM)
        return True
    except OSError:
        return False


def _progress(root: str) -> Dict[str, str]:
    """the last line of the logs of a batch about each star being
    computed: how far its FIP is"""
    out = {}
    for path in glob.glob(os.path.join(root, 'logs', '*.log')) + [
            os.path.join(root, 'run.log')]:
        try:
            with open(path, errors='replace') as handle:
                lines = [line.split(' | ', 1)[-1].strip()
                         for line in handle if line.strip()]
        except OSError:
            continue
        star, last = None, ''
        for line in lines:
            start = re.match(r'batch: (.+) \(\d+ of \d+\)$', line)
            if start:
                star, last = start.group(1), 'started'
            elif star and line.startswith(f'batch: {star} '):
                star = None  # done, or failed
            elif star and not line.startswith('batch: '):
                last = line
        if star:
            out[star] = last
    return out


def cores(report: bool = False) -> Dict[str, Any]:
    """
    The cores of this machine, and how many stars a batch may compute at
    once on two thirds of those that are free (SHARE)

    :param report: bool, the batch makes the detailed report of its
                   candidates (a star then takes four cores, not three)

    :return: dict, total (the cores of the machine), mine (those this
             process may use: fewer under a scheduler or a taskset), load
             (the processes running, averaged over the last minute; None
             when the system does not say), free (mine minus the load, at
             least 1), share (two thirds of them, in cores), each (the
             cores of a star) and jobs (the stars at once: at least 1)
    """
    total = os.cpu_count() or 1
    mine = (len(os.sched_getaffinity(0))
            if hasattr(os, 'sched_getaffinity') else total)
    try:
        load = float(os.getloadavg()[0])
    except (OSError, AttributeError):
        load = None
    free = max(1, mine - int(round(load or 0.0)))
    each = CORES_REPORT if report else CORES
    return dict(total=total, mine=mine, load=load, free=free,
                share=int(SHARE * free), each=each,
                jobs=max(1, int(SHARE * free) // each))


def status(root: str, jobs: Optional[int] = None,
           report: Optional[bool] = None) -> Dict[str, Any]:
    """
    Where the batch of a folder is, said line by line and returned;
    nothing is computed. This Python and what it has of what the batch
    needs, the cores of the machine and the stars at once that two thirds
    of the free ones allow (cores()), then each star: done (its best
    peak), failed (why), being computed (how far), or to do; and whether
    the batch runs.

    :param root: str, the batch folder
    :param jobs: int or None, the stars at once as the batch is set
    :param report: bool or None, the batch makes detailed reports (None:
                   as it was packed)

    :return: dict, state ('not ready': something lacks, see missing;
             'ready': nothing runs, stars are left; 'running'; 'done':
             every star has its result), n, done, failed, running, todo,
             missing (list of str), stars (name, state, told) and cores
             (cores(): its jobs is what to give as --jobs)
    """
    import importlib
    root = os.path.abspath(os.path.expanduser(root))
    held = targets_of(root)
    missing = []
    found = []
    for name in NEEDS:
        try:
            found.append(f'{name} {importlib.import_module(name).__version__}')
        except ImportError:
            missing.append(f'this Python has no {name}')
    if sys.version_info < (3, 9):
        missing.append('Python 3.9 or later is needed')
    log(f'batch {held.get("name")} (koloa {held.get("koloa")}, made '
        f'{held.get("made")}) in {root}', 'info')
    log(f'Python {sys.version.split()[0]} ({sys.executable}), '
        + ', '.join(found), 'value')
    if report is None:
        report = bool((held.get('options') or {}).get('report'))
    have = cores(report)
    log(f'{have["total"]} cores on {socket.gethostname()}'
        + (f', {have["mine"]} for this process' if have['mine']
           != have['total'] else '')
        + (f', {have["load"]:.1f} busy over the last minute'
           if have['load'] is not None else '') + f': {have["free"]} free',
        'value')
    # two thirds of the free cores, a star about three of them (four with
    #   its detailed report): what to start the batch with
    log(f'two thirds of them, {have["share"]} cores, at about '
        f'{have["each"]} a star'
        + (' (with its detailed report)' if report else '')
        + f': {have["jobs"]} star' + ('s' if have['jobs'] != 1 else '')
        + f' at once (--jobs {have["jobs"]})', 'value')
    lines = _progress(root)
    stars = []
    for star in held['targets']:
        name = str(star['name'])
        gone = [path for path in star.get('files') or []
                if not os.path.exists(os.path.join(root, path))]
        missing += [f'{name}: no {path}' for path in gone]
        has = os.path.isdir(os.path.join(root, 'archives', _folder(name)))
        out = os.path.dirname(result_path(root, name))
        mark = _marked(os.path.join(out, RUNNING))
        left = False
        if os.path.exists(result_path(root, name)):
            with open(result_path(root, name)) as handle:
                res = json.load(handle)
            left = (res.get('report') or {}).get('status') == 'running'
        if left:
            # its quick look is done, its detailed report is not
            alive = mark is not None and mark['alive'] is not False
            state = 'running' if alive else 'to do'
            told = ('its detailed report' + (
                f' (on {mark.get("host")})' if alive and mark['alive'] is None
                else '') if alive else 'its detailed report was stopped: '
                'it starts over (its quick look is kept)')
        elif os.path.exists(result_path(root, name)):
            summ = res.get('summary') or {}
            state = 'done' if res.get('status') == 'done' else 'failed'
            told = (f'P = {summ["period"]:.4f} d, FIP {_fip(summ["fip"])}, '
                    if summ.get('period') else '') + (
                f'{res.get("elapsed") or 0:.0f} s' if state == 'done'
                else str(res.get('error') or res.get('status'))) + (
                f'; {res["reading"]["kind"]}' if res.get('reading') else ''
                ) + (f'; report {res["report"]["status"]}' if (
                    res.get('report') or {}).get('status') in (
                        'done', 'failed') else '')
        elif mark is not None and mark['alive'] is not False:
            state = 'running'
            told = (lines.get(name) or 'started') + (
                f' (on {mark.get("host")})' if mark['alive'] is None else '')
        else:
            state = 'to do'
            told = ('stopped as it was being computed: it starts over'
                    if mark is not None else '')
            if not star.get('files') and not has:
                told = 'no file and no archive of it here: it will fail'
        stars.append(dict(name=name, state=state, told=told))
        log(f'  {state:8s} {name:22s} {told}'.rstrip(),
            'warn' if state == 'failed' or gone else 'value')
    count = {key: sum(star['state'] == key for star in stars)
             for key in ('done', 'failed', 'running', 'to do')}
    marks = [_marked(path) for path in glob.glob(
        os.path.join(root, 'logs', '*.pid'))]
    runs = any(mark is not None and mark['alive'] is not False
               for mark in marks) or count['running'] > 0
    for line in missing:
        log(line, 'error')
    state = ('not ready' if missing else 'running' if runs else
             'done' if not count['to do'] else 'ready')
    if jobs and state == 'ready' and int(jobs) > have['jobs']:
        # the stars at once as the batch is set, when that is too many here
        log(f'as set, {jobs} stars at once (about '
            f'{have["each"] * int(jobs)} cores): more than two thirds of '
            f'the free ones; start it with --jobs {have["jobs"]}', 'warn')
    log(f'state: {state}; {count["done"] + count["failed"]} of '
        f'{len(stars)} stars have their result'
        + (f' ({count["failed"]} failed)' if count['failed'] else '')
        + f', {count["running"]} being computed, {count["to do"]} to do',
        'error' if missing else 'info')
    return dict(state=state, n=len(stars), done=count['done'],
                failed=count['failed'], running=count['running'],
                todo=count['to do'], missing=missing, stars=stars,
                cores=have)


def stop(root: str) -> int:
    """
    Stop the batch of a folder that runs on this machine: each of its
    processes is ended, with the chains of the FIP it started. The stars
    done are kept; those that were being computed start over the next
    time. (The tasks of a job array are stopped by the scheduler: scancel.)

    :param root: str, the batch folder

    :return: int, the processes ended
    """
    root = os.path.abspath(os.path.expanduser(root))
    ended, elsewhere = 0, set()
    for path in sorted(glob.glob(os.path.join(root, 'logs', '*.pid'))):
        mark = _marked(path)
        if mark is not None and mark['alive'] is None:
            elsewhere.add(str(mark.get('host')))
            continue
        if mark is not None and mark['alive'] and int(
                mark['pid']) != os.getpid() and _ours(int(mark['pid'])):
            ended += _end(int(mark['pid']))
        _unmark(path)
    for path in glob.glob(os.path.join(root, 'results', '*', RUNNING)):
        mark = _marked(path)
        if mark is None or mark['alive'] is not None:
            _unmark(path)
    log(f'batch: {ended} processes stopped; the stars done are kept'
        + (f'. Others run on {", ".join(sorted(elsewhere))}: stop them '
           f'there (scancel for a job array)' if elsewhere else ''),
        'warn' if elsewhere else 'value')
    return ended


def run(root: str, jobs: int = 1, rules: bool = True, gather: bool = False,
        only: Optional[Sequence[str]] = None, again: bool = False,
        part: Optional[Tuple[int, int]] = None, trend: int = 1,
        report: bool = False, report_fip: float = REPORT_FIP
        ) -> List[Dict[str, Any]]:
    """
    The batch of a folder: each star not done yet (run_target), a few at
    once; then its table, and the summary of its stars as one PDF
    (koloa.batchpdf.summary: results/summary.pdf)

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
    :param report: bool, the detailed report of each star with a candidate
                   (run_target; a report takes far longer than a quick
                   look, and about four cores)
    :param report_fip: float, the FIP below which a peak counts

    :return: list of dict, the table (table())
    """
    root = os.path.abspath(os.path.expanduser(root))
    stars = targets_of(root)['targets']
    if only:
        keep = {name_key(name) for name in only}
        stars = [star for star in stars if name_key(star['name']) in keep]
    os.makedirs(os.path.join(root, 'logs'), exist_ok=True)
    jobs = max(1, int(jobs))
    # this process said to be there (status() and stop() read it)
    mine = os.path.join(root, 'logs', ('run' if part is None else
                                       f'part_{part[0]}_of_{part[1]}')
                        + '.pid')
    _mark(mine)
    try:
        _run(root, stars, jobs, rules, gather, only, again, part, trend,
             report, report_fip)
    finally:
        _unmark(mine)
    rows = table(root)
    log(f'batch: {len(rows)} of {len(targets_of(root)["targets"])} stars in '
        f'{os.path.join(root, "results", "table.csv")}', 'value')
    if part is None:
        # the summary of them all (a part leaves it to who started it; the
        #   tasks of a job array, to the batch started once they ended)
        try:
            from koloa import batchpdf
            batchpdf.summary(root)
        except Exception as err:  # the table and the stars are kept
            log(f'batch: its summary could not be made '
                f'({type(err).__name__}: {err})', 'warn')
    return rows


def _run(root: str, stars: List[Dict[str, Any]], jobs: int, rules: bool,
         gather: bool, only: Optional[Sequence[str]], again: bool,
         part: Optional[Tuple[int, int]], trend: int, report: bool = False,
         report_fip: float = REPORT_FIP) -> None:
    """the stars of run(): its parts started and waited for, or the stars
    of this part one after the other"""
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
        if report:
            base += ['--report', '--report-fip', f'{report_fip:g}']
        if only:
            base += ['--only'] + [str(name) for name in only]
        procs = []
        for it in range(jobs):
            handle = open(os.path.join(root, 'logs', f'part_{it}.log'), 'a')
            # each part leads what it starts (the chains of its FIPs):
            #   stopped together (stop())
            procs.append((subprocess.Popen(
                base + ['--part', f'{it}/{jobs}'], env=env, stdout=handle,
                stderr=subprocess.STDOUT, start_new_session=True), handle))
        log(f'batch: {len(stars)} stars, {jobs} at once (their logs in '
            f'{os.path.join(root, "logs")})', 'info')
        # a kill of this process stops its parts, as Ctrl-C does
        kept = None
        if threading.current_thread() is threading.main_thread():
            def ended(*_):
                raise KeyboardInterrupt
            kept = signal.signal(signal.SIGTERM, ended)
        seen = -1
        try:
            while any(proc.poll() is None for proc, _ in procs):
                # done: its result there, and no report of it going on
                done = sum(os.path.exists(result_path(root, star['name']))
                           and not os.path.exists(os.path.join(
                               os.path.dirname(result_path(
                                   root, star['name'])), RUNNING))
                           for star in stars)
                if done != seen:
                    seen = done
                    log(f'batch: {done} of {len(stars)} stars done', 'value')
                time.sleep(5)
        except KeyboardInterrupt:
            for proc, _ in procs:
                _end(proc.pid)
            log('batch: stopped, the stars done are kept', 'warn')
            raise
        finally:
            if kept is not None:
                signal.signal(signal.SIGTERM, kept)
        for proc, handle in procs:
            handle.close()
            if proc.returncode:
                log(f'batch: a part ended with code {proc.returncode} (see '
                    f'its log)', 'warn')
    else:
        mine = stars if part is None else stars[part[0]::part[1]]
        for it, star in enumerate(mine):
            # done; or done but for its detailed report, which was stopped
            left = False
            if os.path.exists(result_path(root, star['name'])) and not again:
                with open(result_path(root, star['name'])) as handle:
                    left = (json.load(handle).get('report') or {}).get(
                        'status') == 'running'
                if not left:
                    continue
            mark = os.path.join(os.path.dirname(result_path(
                root, star['name'])), RUNNING)
            other = _marked(mark)
            if other is not None and other['alive'] and int(
                    other['pid']) != os.getpid():
                # the batch started twice: not computed twice
                log(f'batch: {star["name"]} is being computed by another '
                    f'process of this batch ({other["pid"]}): left to it',
                    'warn')
                continue
            log(f'batch: {star["name"]} ({it + 1} of {len(mine)})'
                + (': its detailed report, which was stopped' if left
                   else ''), 'info')
            _mark(mark)
            try:
                res = (resume_report(root, star, gather, trend) if left else
                       run_target(root, star, rules, gather, trend, report,
                                  report_fip))
                summ = res.get('summary') or {}
                log(f'batch: {star["name"]} {res["status"]} in '
                    f'{res["elapsed"]:.0f} s' + (
                        f', P = {summ["period"]:.4f} d, FIP '
                        f'{_fip(summ["fip"])}' if summ.get('period') else '')
                    + (f' ({res["error"]})' if res.get('error') else '')
                    + (f'; {res["reading"]["kind"]}' if res.get('reading')
                       else ''), 'value')
            except Exception as err:  # the next star all the same
                log(f'batch: {star["name"]} failed '
                    f'({type(err).__name__}: {err})', 'warn')
            finally:
                _unmark(mark)


SCRIPT = '''#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa: the batch {name}, made {made}

{nstar} stars, each with its files (files/) and its archives (archives/):
the quick FIP of each, its datasets chosen by koloa.datasets, kept in
results/ as it goes (a batch stopped goes on where it was). README.txt,
beside this script, says the whole of it.

    python run_batch.py --check         # nothing computed: is all there,
                                        #   and where is the batch?
    python run_batch.py                 # every star not done yet
    python run_batch.py --jobs 20       # the same, 20 stars at once
                                        #   (--check says how many fit)
    python run_batch.py --part 3/20     # the stars 3, 23, 43... (a job array)
    python run_batch.py --stop          # stop the batch that runs here
    python run_batch.py --summary       # results/summary.pdf again, from
                                        #   what is done (nothing computed)
    python run_batch.py --root /elsewhere/{name} ...

Everything is found from ROOT: nothing else of this script names a path.
"""
import os
import sys

# =============================================================================
# The settings
# =============================================================================
# where this folder is ON THE MACHINE THAT RUNS THE BATCH. Change it when
#   the folder is copied elsewhere (a server, a cluster), or give it as
#   --root: files/, archives/, cache/, koloa_src/ and results/ are found
#   from it.
ROOT = {root!r}

# the stars at once (each takes about three cores, four with its detailed
#   report). --check says how many two thirds of the free cores of this
#   machine allow; --jobs N gives it without changing this line.
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

# the detailed report (koloa --detailed) of each star that has a candidate:
#   a peak with a FIP below REPORT_FIP that is neither a known planet nor a
#   drift as long as the series. The GP of its FIP is an SHO at the
#   rotation period when one is published, by period band otherwise. A
#   report takes far longer than a quick look, and about four cores.
REPORT = {report}
REPORT_FIP = {report_fip}


# =============================================================================
# The batch (under main: koloa computes in processes of its own, which read
#   this file again as they start and must not start the batch themselves)
# =============================================================================
def main():
    root, part, jobs = ROOT, None, JOBS
    args = sys.argv[1:]
    if '--root' in args:
        root = args[args.index('--root') + 1]
    if '--jobs' in args:
        jobs = int(args[args.index('--jobs') + 1])
    if '--part' in args:
        one, of = args[args.index('--part') + 1].split('/')
        part = (int(one), int(of))
    root = os.path.abspath(os.path.expanduser(root))
    if not os.path.exists(os.path.join(root, 'targets.json')):
        here = os.path.dirname(os.path.abspath(__file__))
        raise SystemExit(
            'ROOT is ' + root + ': there is no targets.json there.\\n'
            'Set ROOT, at the top of this script, to where the folder is '
            'on this machine, or give it: --root PATH'
            + (' (this script is in ' + here + ').' if os.path.exists(
                os.path.join(here, 'targets.json')) else '.'))
    # what the batch needs of this Python, said before koloa asks for it
    lacking = []
    for module in {needs!r}:
        try:
            __import__(module)
        except ImportError:
            lacking.append(module)
    if lacking or sys.version_info < (3, 9):
        raise SystemExit(
            'This Python (' + sys.executable + ', '
            + sys.version.split()[0] + ') '
            + ('has no ' + ', '.join(lacking) if lacking
               else 'is older than 3.9')
            + ': the batch needs Python 3.9 or later with '
            + ', '.join({needs!r}) + ' (README.txt, step 3).')
    # what koloa fetched once, and koloa itself, as they were packed
    os.environ['KOLOA_CACHE'] = os.path.join(root, 'cache')
    src = os.path.join(root, 'koloa_src')
    if os.path.isdir(src):
        sys.path.insert(0, src)
        os.environ['PYTHONPATH'] = src + os.pathsep + os.environ.get(
            'PYTHONPATH', '')
    from koloa import survey
    if '--stop' in args:
        survey.stop(root)
    elif '--summary' in args:
        from koloa import batchpdf
        survey.table(root)
        batchpdf.summary(root)
    elif '--check' in args:
        # 0 when the batch can run (or runs, or is done), 1 when not
        told = survey.status(root, jobs=jobs, report=REPORT)
        raise SystemExit(1 if told['state'] == 'not ready' else 0)
    else:
        try:
            survey.run(root, jobs=jobs, rules=RULES, gather=GATHER,
                       only=ONLY, again=AGAIN, part=part, trend=TREND,
                       report=REPORT, report_fip=REPORT_FIP)
        except KeyboardInterrupt:  # stopped (--stop, Ctrl-C): said, and all
            raise SystemExit(130)


if __name__ == '__main__':
    main()
'''

SUBMIT = '''#!/bin/bash
# koloa: the batch {name} as a job array of SLURM (the clusters of the
# Alliance, formerly Compute Canada): one star per task, {nstar} tasks, at
# most {jobs} at once. Set your account and what loads Python, then, from
# this folder (the tasks run where they were submitted from):
#     sbatch submit.sh
#SBATCH --job-name=koloa_{name}
#SBATCH --account=def-CHANGE_ME
#SBATCH --array=0-{last}%{jobs}
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/slurm_%A_%a.log

module load python scipy-stack
cd "$SLURM_SUBMIT_DIR"
python run_batch.py --root "$PWD" --part ${{SLURM_ARRAY_TASK_ID}}/{nstar}
'''

README = '''koloa: the batch {name}
{rule}

{intro}

This folder runs on its own. koloa itself is in it (koloa_src/), with the
files, the archives and the catalogues the batch reads. Nothing is
installed and nothing is asked of the network. It needs Python 3.9 or
later with numpy, scipy and matplotlib, and that is all.

This file is for whoever runs the batch: a person, or a Claude session (or
another agent) on the machine that runs it. An agent reads the next two
sections first; the numbered steps are the detail of each. CLAUDE.md,
beside this file, is the section written for the agent, where a Claude
session started in this folder reads it on its own.


IN SHORT
--------

    tar xzf {name}.tar.gz
    cd {name}
    python run_batch.py --root "$PWD" --check     # seconds; computes nothing
    nohup python -u run_batch.py --root "$PWD" --jobs N > run.log 2>&1 &
    python run_batch.py --root "$PWD" --check     # again, to see where it is

N is the stars computed at once: the number --check gives, two thirds of
the free cores of the machine at about three cores a star.

python is the Python 3.9 or later that has numpy, scipy and matplotlib:
python3 on many machines. When it is over, results/summary.pdf is the
batch in one PDF: the table of the stars, the candidates first, then a
page for each. Then results/ is brought back (step 5) and opened in
koloa's page (step 6).


FOR A CLAUDE SESSION (OR ANOTHER AGENT) ASKED TO RUN THIS BATCH
---------------------------------------------------------------

The task: compute every star of this folder with run_batch.py, on this
machine, then report its table. Everything needed is in this folder. In
this order:

  a. Work from this folder, the one that holds this README.txt and
     targets.json, and give its path every time: --root "$PWD". Do not
     edit ROOT in the script instead unless you are asked to.

  b. Find a Python 3.9 or later with numpy, scipy and matplotlib:

         python3 -c "import sys, numpy, scipy, matplotlib; print(sys.version)"

     If it fails, look at what the machine already offers (module avail
     python, conda env list, an environment the person named) and use
     that. Do not install anything and do not create an environment
     without asking the person: stop, and say what is missing. In the
     lines below, python stands for the Python found here (python3, or
     its full path).

  c. Update the koloa of this machine, if it has one: it may be out of
     date. This is the one thing to fetch, and it needs the network.

         python -m pip show koloa        # is it installed for this Python?

     If it is:

         python -m pip install --upgrade --force-reinstall --no-deps \\
             git+https://github.com/eartigau/koloa.git

     If the person keeps a clone of it here instead (a folder with .git
     and koloa/ in it, on PYTHONPATH or installed with pip -e), update
     that one: git -C <the clone> pull --ff-only. If koloa is not on
     this machine at all, install nothing: this batch carries its own.
     If the update fails (no network, no right to write), say so and
     go on: the batch does not wait for it.

     Know what it changes. This batch does not run with the koloa of
     the machine: it runs with koloa_src/, koloa as it was when the
     batch was packed,

         koloa {version}{commit}, {made}

     so that every star is computed by the code run_batch.py was
     written for. The update is for everything else run here with
     koloa (koloa, koloanui, a report by hand), and it tells you
     whether this folder is behind: if the koloa you now have is newer
     than the one packed (its last commit, git log -1 in a clone, or
     the date of its install), say so in your report. The person may
     then want to pack the batch again; do not swap koloa_src/
     yourself.

  d. Check. This computes nothing and takes seconds:

         python run_batch.py --root "$PWD" --check

     Its last line is the state of the batch, after the time:
     "... | state: ready; 0 of {nstar} stars have their result, ...".
     "state: ready" (exit code 0): go on. "state: not ready" (exit code
     1): the lines above it say what is missing; report them and stop.
     "state: running" or "state: done": the batch was already started
     here; do not start it again, go to g or h.

  e. Choose how to run it, from what the machine is.
     - No scheduler (a server reached by ssh): step 3, on two thirds of
       the cores that are free, no more and not many fewer. Find how
       many cores there are and how many are busy yourself:

           nproc           # the cores this shell may use
           uptime          # the load: the cores busy, over 1, 5, 15 min

       free cores = nproc - the 1 minute load (rounded). A star takes
       about three cores (four when REPORT is True: {report} as packed).
       So the stars at once are

           N = (2/3 x free cores) / 3, rounded down, at least 1

       (/ 4 with REPORT). --check does the same sum and prints it:
       "... : N stars at once (--jobs N)". Check that your number and
       its number agree; if they differ by more than one, look again
       (a load that is changing, cores this process may not use) before
       starting. The JOBS of the script ({jobs} as packed) was set on
       another machine: do not use it, give N as --jobs. Look at the
       memory too once it runs (free -g): if it runs short, stop the
       batch (--stop) and start it again with a smaller N; nothing
       done is lost.
     - A cluster (sbatch exists): step 4, never on a login node. The
       scheduler shares the cores there: the two thirds do not apply.
       submit.sh needs an account (--account); if you were not given
       one, ask.

  f. Start it once, detached, so that it outlives your shell, with the
     N of step e:

         nohup python -u run_batch.py --root "$PWD" --jobs N > run.log 2>&1 &

     It runs for minutes to many hours (step 3 has real timings). Do
     not wait for it in the foreground, and do not use tail -f, which
     never returns.

  g. Follow it by asking again, every few minutes at first, then less
     often:

         python run_batch.py --root "$PWD" --check

     It lists each star: done (with its best peak), failed (with why),
     running (with how far its FIP is and the time left), or to do. The
     batch is over when the last line reads "state: done". If it
     reads "state: ready" again while stars are still to do, the batch
     stopped (the machine restarted, say): start it again as in f.
     The batch does not need you: if your session may end before it
     does, tell the person so, and give them this --check line.

  h. Then read results/table.csv (step 6 says what each column is) and
     report: how many stars are done and how many failed, each failure
     with its error as the table gives it; how the stars read (the
     verdict column: candidate, known, drift, nothing) and the
     candidates with their best peak (period, fip, K), without calling
     any a planet (step 6: the table is a list of questions); which
     detailed reports were made and which failed (the report column);
     how many cores the machine has, how many were free, and the --jobs
     you gave; whether the koloa of the machine was updated, and whether
     it is newer than the one packed; how long it took; where
     results/summary.pdf is; and the
     line that brings the results back (step 5) with the real path of
     this folder in it.

What not to do, however convenient it looks:

  - Fetch or install nothing but the update of koloa of step c: no other
    pip, no conda install, no git clone, and GATHER stays False. The
    batch needs no network.
  - files/, archives/, cache/, koloa_src/ and targets.json are to be
    read, not changed: do not edit, move or delete them.
  - results/ is the work done. Never delete it to start clean: a star is
    computed again with ONLY and AGAIN (step 7).
  - Do not start the batch a second time while it runs.
  - If a star fails, or takes far longer than the others, leave it and
    let the others end. Report it. Do not change the code in koloa_src/
    to get past it.
  - To stop the batch: python run_batch.py --root "$PWD" --stop. Do not
    kill Pythons on a guess: others on this machine are not yours.
  - The velocities in files/ and archives/ may not be public. Do not
    copy this folder elsewhere, upload it, or quote its data; give the
    numbers of the table to the person who asked, and to no one else.


1. COPY IT TO THE MACHINE THAT RUNS IT
--------------------------------------

From the machine that made it, where the tar is beside this folder:

    rsync -av --progress {name}.tar.gz me@server:/where/the/batches/go/

(scp does the same; rsync goes on where it was if the link drops.) There:

    cd /where/the/batches/go
    tar xzf {name}.tar.gz
    cd {name}

The tar holds everything: nothing else is copied.


2. TELL IT WHERE IT IS: ROOT
----------------------------

Everything is found from one path, ROOT, the first setting of
run_batch.py: where this folder is on the machine that runs it. As packed:

    ROOT = {root!r}

If the folder is somewhere else, either change that line, or give the path
when starting the script, which leaves it as packed:

    python run_batch.py --root "$PWD" ...

A wrong ROOT stops at once and says so. Nothing is computed elsewhere.


3. CHECK, THEN RUN
------------------

Load what gives Python 3.9 or later with numpy, scipy and matplotlib: a
conda environment, or on a cluster of the Alliance
"module load python scipy-stack". Then:

    python run_batch.py --root "$PWD" --check

It computes nothing. It says which Python it is and which numpy, scipy and
matplotlib it has, how many cores the machine has, and of each star
whether it is done, failed, running or to do. Its last line is the state
of the batch: not ready (something is missing, named in the lines above;
exit code 1), ready, running or done (exit code 0).

Then the batch, left running after you log out, N stars at once:

    nohup python -u run_batch.py --root "$PWD" --jobs N > run.log 2>&1 &

or the same without nohup inside screen or tmux. To see where it is:

    python run_batch.py --root "$PWD" --check
    tail run.log                    # the stars done, as they end
    tail logs/part_0.log            # one of the stars being computed
    ls results/                     # a folder for each star done

Cores: each star is a Python of its own that takes about three cores.
The batch is meant to take two thirds of the cores that are free on the
machine, and leave the rest: N = (2/3 x free cores) / 3 stars at once,
rounded down, where the free cores are those of the machine (nproc) less
its load (uptime). --check prints the cores, the load and that N; give
it as --jobs N. Without --jobs the batch takes the JOBS of run_batch.py
({jobs} as packed, about {cores} cores), which was set where the batch was
packed and knows nothing of this machine.

How long: a star takes longer with more nights and with more datasets.
Measured on a 24-core server in October 2026, each star with its archives:

    GJ 1214    154 nights in 3 datasets       6 minutes
    GJ 581     515 nights in 6 datasets      31 minutes
    GJ 699    1075 nights in 6 datasets      89 minutes

A star with many datasets and a signal of hundreds of m/s can take far
longer (GJ 876, 9 datasets: half of its first pass after 7 hours). The
others do not wait for it: each star is on its own.

With REPORT = True ({report} as packed), a star that has a candidate has
its detailed report too (step 6). A report takes far longer than a quick
look, tens of minutes to hours, and about four cores instead of three:
count the cores with that. --check says "its detailed report" of a star
that is at it.

To stop it:

    python run_batch.py --root "$PWD" --stop

It ends the processes of this batch on this machine and no other. The
stars done are kept. To go on, start the batch again with the same line:
the stars done are not computed again, those that were being computed
start over, the others follow.


4. ON A CLUSTER, WITH SLURM
---------------------------

submit.sh is the same batch as a job array: a star a task, {jobs} at once.
Set in it:

    #SBATCH --account=def-CHANGE_ME     your allocation
    #SBATCH --time=03:00:00             more for stars with many datasets
    module load python scipy-stack      what gives Python where you are

then, from this folder (the tasks run where they are submitted from, and
find the folder by that):

    sbatch submit.sh
    squeue -u $USER                     the tasks waiting and running
    python run_batch.py --root "$PWD" --check
    tail logs/slurm_*_0.log             the log of task 0

A task ended by its time limit loses the star it was computing and no
other: submit again with a longer --time, and only the stars not done are
computed (a star whose quick look was done and whose detailed report was
cut short keeps its quick look: only its report is taken again). With
REPORT = True, give the tasks more time, and --cpus-per-task=4 is right.
When every task has ended,

    python run_batch.py --root "$PWD"

computes nothing more and writes the table of them all and
results/summary.pdf. To stop the tasks: scancel (--stop ends only what
runs on the machine it is typed on).

The nodes of most clusters have no network: GATHER stays False, as packed.
The archives are those of archives/, gathered before packing.


5. BRING THE RESULTS BACK
-------------------------

Only results/ is needed, into the folder of the batch on the machine that
made it. From that machine:

    rsync -av me@server:{root}/results/ \\
          {name}/results/

It can be done at any time: the stars done so far come back, the others
the next time.


6. LOOK AT THEM
---------------

In koloa's page (koloanui), Survey tab, "The results of a batch folder":
give the folder {name} and "Open them". The batch opens as a table in
the Batch FIP tab; Open, on a line, puts the star in the Analysis tab with
its FIP as it was computed, to fold, tick and report like any other.

Without the page:

    results/summary.pdf            the batch in one PDF: the table of the
                                   stars, the candidates first, then the
                                   page of each
    results/table.csv              a line for each star done
    results/<star>/quicklook.pdf   its page of summary first, then its
                                   series, its FIP, its folds, its numbers
    results/<star>/report/         its detailed report, when it has a
                                   candidate and REPORT is True
    results/<star>/result.json     its line of the table, how its peaks
                                   read, and what was done with each of
                                   its datasets and why
    results/<star>/quick.json      the FIP itself, as the page draws it

The page of a star says how its quick look reads. Its verdict:

    candidate   a peak with a FIP below {report_fip} that is neither a
                known planet nor a drift
    known       its peaks below that FIP are known planets (within 1 % of
                the period of one of the NASA Exoplanet Archive's)
    drift       they are as long as the series (more than half of it)
    nothing     no peak below that FIP

then each numbered peak with what it is (a known planet, at the rotation
of the star or one of its harmonics, at a year, a drift, or none of
these), its FIP with the known planets and the rotation marked, its two
best folds, its velocities, its datasets and what the rules left out, its
acceleration, the transit looked for, its published rotation periods, and
what its detailed report found. It is a reading of a quick look, which
has no GP: it sorts the stars, it does not decide on a planet.

    python run_batch.py --root "$PWD" --summary

makes results/summary.pdf again from what is done, at any time, and
computes nothing.

The detailed report of a star (REPORT = True) is koloa's own (koloa
--detailed), on the series the quick look used, whatever it came from:
files, archives, or both. That series is beside it
(report/velocities.csv), with its log (run.log), its figures, and
<star>_report.pdf, or <star>_report.txt and .tex where there is no
pdflatex. The GP of its FIP is an SHO at the rotation period of the star
when one is published (the P rot of section 9; else the NASA Exoplanet
Archive's), a local GP by period band otherwise; result.json says which.
With GATHER = False it asks nothing of the network, so it has no TESS
light curve to check its signals against.

The columns of table.csv:

    name, sptype, distance   the star, its spectral type, its distance [pc]
    status                   done, or failed (error says why)
    files                    how many files of velocities it came with
    datasets                 used/all: 6/17 is 6 datasets used of the 17
                             the star has. The others are releases of the
                             same spectra that another dataset has more
                             precisely, or constrain nothing; result.json
                             names each and says why
    nights                   the nightly means used
    baseline                 the days they span
    period                   the best peak of the quick FIP [days]
    fip, fip_alone           its false inclusion probability: of the
                             period or any of its aliases, and of the
                             period alone (0 reads: below 1e-300)
    K, K_err                 the semi-amplitude of its sinusoid [m/s]
    rms                      the scatter of the nights about it [m/s]
    accel, accel_err         the acceleration of the star, fitted with it
                             [m/s/yr]
    accel_sigma              how many sigma that acceleration is from zero
    transit, transit_snr     the transit looked for at that period in the
                             light curve of the star: plausible, or none,
                             with its signal to noise. "no light curve":
                             the batch was packed without the light curves
    verdict                  how its quick look reads: candidate, known,
                             drift, nothing (above)
    report                   its detailed report: done, failed, none (no
                             candidate, or not asked for), running
    elapsed                  how long the star took [s]
    error                    why a star failed

The table is a list of questions, not of planets. A best peak about as
long as the baseline is a slow drift between datasets, not an orbit. A
peak near 1 day, or near the rotation period of the star, wants a second
look. Each line is to be opened and looked at.


7. THE SETTINGS, AT THE TOP OF run_batch.py
-------------------------------------------

    ROOT     where this folder is (step 2)
    JOBS     the stars at once ({jobs} as packed), about three cores each;
             --jobs N on the command line takes its place (step 3)
    RULES    True: of the datasets of a star, the best release of the same
             spectra, and not those that constrain nothing. False: all
    GATHER   False: the archives as packed. True asks the archives for
             what a star lacks: it needs the network
    ONLY     [] for every star, or a few of them: ['GJ 581', 'GJ 876']
    AGAIN    True computes again the stars already done
    TREND    1: an acceleration fitted with the signals. 2: its change too
    REPORT   True: the detailed report of each star that has a candidate
             (step 6). Hours, for a batch with many candidates
    REPORT_FIP   the FIP below which a peak counts ({report_fip} as packed)

and on the command line: --root PATH, --jobs N (the stars at once, in
place of JOBS), --check, --stop, --summary, and --part I/N (the stars I,
I+N, I+2N...: a task of a job array).

To compute one star again: ONLY = ['GJ 581'] and AGAIN = True, then the
batch as before; or delete results/GJ_581/ and start the batch again.


8. WHAT IS IN THIS FOLDER
-------------------------

    README.txt      this file
    CLAUDE.md       its section for a Claude session, read by Claude Code
                    when it is started in this folder
    run_batch.py    the script; its first setting is ROOT
    submit.sh       the same as a SLURM job array
    targets.json    the stars, and the files of each
    files/          the files of velocities of each star
    archives/       the archives of each star,
                    {when}
    cache/          what koloa fetched once: the NASA Exoplanet Archive,
                    the lists of the surveys, APERO's names
    koloa_src/      koloa as it was when the batch was packed
    results/        written by the batch: a folder for each star, and
                    table.csv, table.json, summary.pdf
    logs/           written by the batch: the log of each part, and of
                    each task of a job array


9. THE STARS
------------

{stars}


10. WHEN SOMETHING GOES WRONG
-----------------------------

"ROOT is ...: there is no targets.json there"
    The script was told a folder that is not this one: step 2.

"This Python (...) has no numpy" (or scipy, matplotlib), or
ModuleNotFoundError
    This Python is not the one that has them: step 3.

A detailed report is "failed"
    results/<star>/report/run.log has the whole, and result.json its
    last line. The quick look of the star is kept. To take it again:
    ONLY = ['<star>'] and AGAIN = True.

A star is "failed"
    Its line of results/table.csv, and results/<star>/result.json, say
    why; the log of its part (logs/part_*.log, or logs/slurm_*.log) has
    the whole. The other stars are not affected. A star with no file and
    no archive has no velocity to look at, and fails saying so.

Nothing seems to happen
    A star says nothing between the tenths of its FIP, which can be an
    hour apart for a large one. --check shows the last line of each star
    being computed, with the time left.

The machine was restarted, or the batch was killed
    Start it again with the same line. --check first tells which stars
    were being computed: they start over.

URLError, "could not be asked", a timeout
    GATHER is True on a machine with no network: set it back to False.

The batch starts itself again and again
    The lines of main() were copied out of it. koloa computes in processes
    of its own, which read run_batch.py again as they start: everything
    that starts the batch must stay under main(), as packed.


NOT PUBLIC BY DEFAULT
---------------------

files/ holds your own files of velocities, and archives/ what the archives
gave you, which with a DACE key can be more than the public data. This
folder and its tar are no more public than they are.

koloa: https://github.com/eartigau/koloa
'''


def _commit() -> Optional[str]:
    """the commit of this koloa, when it can be told: a clone of its
    repository, or an install from it (pip keeps where it came from)"""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if os.path.isdir(os.path.join(here, '.git')):
        try:
            out = subprocess.run(['git', '-C', here, 'rev-parse', '--short',
                                  'HEAD'], capture_output=True, text=True,
                                 timeout=10)
            if out.returncode == 0 and out.stdout.strip():
                return out.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    try:
        from importlib import metadata
        text = metadata.distribution('koloa').read_text('direct_url.json')
        return str(json.loads(text)['vcs_info']['commit_id'])[:7]
    except Exception:  # not said: the version alone
        return None


#: the section of the README of a batch written for an agent, as it starts
#: and as it ends (CLAUDE.md is that section)
AGENT_FROM = 'FOR A CLAUDE SESSION (OR ANOTHER AGENT) ASKED TO RUN THIS BATCH'
AGENT_TO = '1. COPY IT TO THE MACHINE THAT RUNS IT'


def claude_md(readme: str, name: str, nstar: int, made: str) -> str:
    """
    CLAUDE.md of a batch folder: the section of its README written for a
    Claude session, in the file Claude Code reads on its own when it is
    started in the folder

    :param readme: str, the README of the batch, filled
    :param name: str, the name of the batch
    :param nstar: int, its stars
    :param made: str, when it was made

    :return: str, the text of CLAUDE.md
    """
    start = readme.index(AGENT_FROM)
    start = readme.index('\n', readme.index('\n', start) + 1) + 1
    section = readme[start:readme.index(AGENT_TO)].strip('\n')
    return (f'# koloa: the batch {name}\n\n'
            f'This folder is a batch of koloa, made {made}: {nstar} stars '
            f'to compute with `run_batch.py`. `README.txt`, beside this '
            f'file, is the whole guide. This file is its section for you, '
            f'a Claude session (or another agent) asked to run the batch on '
            f'this machine: read it first, then follow it in order. The '
            f'steps it names by number (step 3, step 5...) are the numbered '
            f'sections of `README.txt`.\n\n'
            f'## What you are asked to do\n\n{section}\n')


def _packed(folder: str, star: Dict[str, Any]) -> Tuple[str, bool]:
    """a star of a batch folder, as a line of its README: its files, and
    what its archives hold as packed (told whether they were gathered)"""
    here = os.path.join(folder, 'archives', _folder(star['name']))
    gathered = False
    told = 'none'
    try:
        with open(os.path.join(here, 'manifest.json')) as handle:
            held = json.load(handle).get('archives') or {}
        parts = [f'{label} {held[key]["npoints"]}' for key, label in (
            ('dace', 'DACE'), ('carmenes', 'CARMENES DR1'),
            ('published', 'VizieR')) if (held.get(key) or {}).get('npoints')]
        gathered = True
        told = ', '.join(parts) or 'gathered: no velocity'
        if (held.get('tess') or {}).get('status') == 'ok':
            told += ', TESS'
    except (OSError, ValueError):
        if os.path.isdir(here):
            told = 'what the cross-match kept of DACE only (not gathered)'
    dist = star.get('distance')
    name, sptype = str(star['name'])[:24], str(star.get('sptype') or '')[:9]
    # its published rotation period, the first: the one its detailed
    #   report would take for its GP
    spins = [one for one in star.get('rotation') or [] if one.get('period')]
    spin = f'{float(spins[0]["period"]):.4g}' if spins else ''
    return (f'    {name:24s} {sptype:9s} '
            f'{f"{dist:7.2f}" if dist else "       "} '
            f'{len(star.get("files") or []):5d} {spin:>9s}  {told}'), gathered


def pack(stars: Sequence[Dict[str, Any]], name: str, out: str = '.',
         root: str = 'archives', server_root: Optional[str] = None,
         gather: bool = True, tess: bool = False, jobs: int = 6,
         rules: bool = True, trend: int = 1, refresh: bool = False,
         api_key: Any = None, tar: bool = True,
         progress: Optional[Callable[[int, int, str], None]] = None,
         workers: int = 3, report: bool = False,
         report_fip: float = REPORT_FIP) -> Dict[str, Any]:
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
    :param report: bool, the REPORT of the script: the detailed report of
                   each star with a candidate
    :param report_fip: float, its REPORT_FIP

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
            # its manifest says where it was gathered: a path of this
            #   machine, which the batch does not read and need not carry
            kept = os.path.join(folder, 'archives', tag, 'manifest.json')
            if os.path.exists(kept):
                with open(kept) as handle:
                    manifest = json.load(handle)
                manifest['folder'] = os.path.join('archives', tag)
                with open(kept, 'w') as handle:
                    json.dump(manifest, handle, indent=1)
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
            archives=star.get('archives'),
            # its published rotation periods: its page says whether a peak
            #   is at one, and its detailed report takes the first
            rotation=star.get('rotation') or []))
    # who each star is (its names, its type), kept with what koloa fetched
    #   once: asked here, where there is the network, for the stars that
    #   were never asked (the batch asks nothing where it runs)
    def named(star):
        try:
            from koloa.archive import resolve
            resolve(str(star['name']))
        except Exception:  # the batch goes on without (said in its log)
            pass
    with ThreadPoolExecutor(max_workers=6) as pool:
        list(pool.map(named, stars))
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
    commit = _commit()
    lines = [_packed(folder, star) for star in targets]
    ngath = sum(gathered for _, gathered in lines)
    when = ('gathered when the batch was packed' if gather else
            'as they were on the machine that packed it')
    words = dict(name=name, made=made, nstar=len(targets),
                 root=server_root or folder, jobs=int(jobs),
                 rules=bool(rules), gather=False, trend=int(trend),
                 report=bool(report), report_fip=float(report_fip),
                 last=max(len(targets) - 1, 0), when=when,
                 version=__version__, needs=NEEDS,
                 commit=f', commit {commit}' if commit else '',
                 cores=CORES * int(jobs), rule='=' * (17 + len(name)),
                 intro=textwrap.fill(
                     f'Made {made} with koloa {__version__}'
                     + (f' (commit {commit})' if commit else '')
                     + f'. {len(targets)} '
                     f'stars: the quick FIP of each, with its files of '
                     f'velocities and its archives, and a table of their '
                     f'best peaks. With them: {nfile} files of velocities, '
                     + (f'and the archives of {ngath} of the stars ({when}: '
                        f'DACE, CARMENES DR1, the surveys and papers on '
                        f'VizieR' + (', the light curves of TESS' if tess
                                     else '') + ').' if ngath else
                        'and no archives (none was gathered where the '
                        'batch was packed).'), 76),
                 stars='\n'.join(
                     [f'    {"star":24s} {"type":9s} {"d [pc]":>7s} '
                      f'{"files":>5s} {"P rot [d]":>9s}  its archives']
                     + [line for line, _ in lines]))
    with open(os.path.join(folder, TARGETS), 'w') as handle:
        json.dump(dict(name=name, made=made, koloa=__version__,
                       koloa_commit=commit,
                       options=dict(rules=bool(rules), trend=int(trend),
                                    report=bool(report),
                                    report_fip=float(report_fip)),
                       targets=targets), handle, indent=1, default=_plain)
    for file, text in (('run_batch.py', SCRIPT), ('submit.sh', SUBMIT),
                       ('README.txt', README)):
        with open(os.path.join(folder, file), 'w') as handle:
            handle.write(text.format(**words))
    # the section of the README for a Claude session, where Claude Code
    #   reads it on its own
    with open(os.path.join(folder, 'CLAUDE.md'), 'w') as handle:
        handle.write(claude_md(README.format(**words), name, len(targets),
                               made))
    made_tar = None
    if tar:
        made_tar = folder + '.tar.gz'
        with tarfile.open(made_tar, 'w:gz') as handle:
            handle.add(folder, arcname=name, filter=lambda info: (
                None if os.path.basename(info.name) in ('results', 'logs')
                and info.isdir() else info))
            # what the batch writes in, empty: there once unpacked (a
            #   scheduler does not make the folder of its logs)
            for empty in ('results', 'logs'):
                info = tarfile.TarInfo(f'{name}/{empty}')
                info.type, info.mode = tarfile.DIRTYPE, 0o755
                info.mtime = time.time()
                handle.addfile(info)
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
    parser.add_argument('--report', action='store_true',
                        help='the detailed report of each star with a '
                             'candidate')
    parser.add_argument('--report-fip', type=float, default=REPORT_FIP,
                        help='the FIP below which a peak counts')
    args = parser.parse_args(argv)
    part = None
    if args.part:
        one, of = args.part.split('/')
        part = (int(one), int(of))
    try:
        run(args.root, jobs=args.jobs, rules=not args.no_rules,
            gather=args.gather, only=args.only, again=args.again, part=part,
            trend=args.trend, report=args.report,
            report_fip=args.report_fip)
    except KeyboardInterrupt:  # stopped: said by run(), and all
        raise SystemExit(130)


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

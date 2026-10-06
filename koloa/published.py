#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The published velocities of a star on VizieR, one star at a time, beyond
DACE and CARMENES: the large surveys that give their velocities star by
star, and the tables of the papers that SIMBAD lists for the star

The surveys (precise velocities, a few m/s; the instruments of the United
States that DACE does not have among them):

- the California Legacy Survey (Rosenthal et al. 2021, ApJS 255, 8;
  J/ApJS/255/8): Keck HIRES (before and after its 2004 upgrade), the APF
  and the Lick Hamilton, 719 stars;
- Keck HIRES carried to 2023 and corrected for its nightly zero points
  (Teklu et al. 2025, A&A 702, A68; J/A+A/702/A68), and the release before
  it (Tal-Or et al. 2019, MNRAS 484, L8; J/MNRAS/484/L8);
- the Lick Hamilton of the Lick planet search (Fischer et al. 2014, ApJS
  210, 5; J/ApJS/210/5), one instrument per dewar;
- HARPS by SERVAL with its nightly zero points (Trifonov et al. 2020, A&A
  636, A74; J/A+A/636/A74), before and after the 2015 fibre change.

The list of the stars of each survey (its names and positions) is kept on
this machine (CACHE), asked once: a star not in a survey's list (none
within RADIUS) is not asked of it. A survey's star has its velocities
asked by the name the survey gives it; their columns are asked for by name
(VizieR leaves some errors out of a request for every column), in m/s
whatever unit a ReadMe declares.

The papers: the bibliography of the star in SIMBAD, each paper's VizieR
catalogue (J/<journal>/<volume>/<page>), and in it the tables of
velocities of the star (koloa.literature.vizier_velocities), each paper's
instruments named after it ('HARPS (Trifonov+ 2018)') to keep their own
offsets.

The papers whose catalogue the CDS does not have, or has without a table
of velocities, are kept (CACHE/catalogues.json) and not asked again, for
any star. CARMENES DR1 (Ribas et al. 2023) is not looked at as a paper:
koloa gathers it as an archive of its own.

Only an instrument whose median error is below MAX_ERROR is kept (no km/s
velocities of a binary survey). Each source is kept whole: the same
spectra are often published more than once (HIRES in three surveys, HARPS
on DACE and by SERVAL), and which release of them is used is chosen when
the series is put together (koloa.datasets: the most precise, the others
left on disk to be asked back). Two velocities are the same spectrum
within SAME_ANY, or within SAME_FAMILY for one spectrograph (the start of
an exposure, its middle); a note says how many of a source's spectra the
sources before it have.

Created on 2026-10-03

@author: artigau
"""
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ElementTree
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from koloa import literature
from koloa.archive import CACHE as ARCHIVE_CACHE
from koloa.data import RVData, merge
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: SIMBAD's TAP service (the bibliography of a star)
SIMBAD_TAP = 'https://simbad.cds.unistra.fr/simbad/sim-tap/sync'
#: the largest distance to a survey's star [arcsec]
RADIUS = 30.0
#: the largest median error of an instrument kept [m/s]
MAX_ERROR = 20.0
#: the earliest year of a paper whose tables are looked at
MIN_YEAR = 1995
#: where the lists of the surveys' stars and what is known of each
#: catalogue are kept
CACHE = os.path.join(os.path.dirname(ARCHIVE_CACHE), 'published')
#: a spectrum already kept: within this of a velocity of any instrument,
#: or of one of the same spectrograph [days]
SAME_ANY = 1.0 / 1440
SAME_FAMILY = 10.0 / 1440
#: the fewest velocities of an instrument added to a series (each has its
#: own offset: one or two say nothing)
MIN_POINTS = 3
#: catalogues not looked at as papers: archives of koloa's own
NOT_PAPERS = {'J/A+A/670/A139'}
#: spectrographs whose velocities, published twice, are the same spectra
#: (the longer names first)
FAMILIES = ('HARPS-N', 'HARPS', 'HIRES', 'CARMENES', 'ESPRESSO', 'SPIRou',
            'NIRPS', 'SOPHIE', 'CORALIE', 'NEID', 'EXPRES', 'KPF', 'HPF',
            'PFS', 'APF', 'Lick', 'UVES', 'MAROON-X', 'IRD')
#: the surveys, in their order of preference (a spectrum in two is kept
#: from the first): the table of stars (found by position) and its column
#: of names, the table of velocities and its columns, the time offset to
#: BJD - 2400000 [days], and the instruments
SURVEYS: List[Dict[str, Any]] = [
    dict(key='teklu25', reference='Teklu et al. 2025', catalogue='J/A+A/702/A68',
         stars='J/A+A/702/A68/table1', name='Name',
         rvs='J/A+A/702/A68/tablea1', link='Name', time='BJD', rv='RVcor',
         err='e_RVcor', offset=2400000.0, label='HIRES (Teklu+ 2025)'),
    dict(key='cls21', reference='Rosenthal et al. 2021 (California Legacy '
         'Survey)', catalogue='J/ApJS/255/8', stars='J/ApJS/255/8/table2',
         name='CPS', rvs='J/ApJS/255/8/table6', link='CPS', time='BJD',
         rv='RVel', err='e_RVel', offset=2400000.0, inst='Inst',
         names={'k': 'HIRES-k (CLS)', 'j': 'HIRES (CLS)', 'apf': 'APF (CLS)',
                'lick': 'Lick (CLS)'}),
    dict(key='talor19', reference='Tal-Or et al. 2019', catalogue='J/MNRAS/484/L8',
         stars='J/MNRAS/484/L8/stars', name='Name',
         rvs='J/MNRAS/484/L8/table1', link='Name', time='BJD', rv='RVcor',
         err='e_RVcor', offset=2400000.0, label='HIRES (Tal-Or+ 2019)'),
    dict(key='fischer14', reference='Fischer et al. 2014 (Lick)',
         catalogue='J/ApJS/210/5', stars='J/ApJS/210/5/table1', name='Name',
         rvs='J/ApJS/210/5/table2', link='Name', time='JD', rv='RVel',
         err='e_RVel', offset=0.0, inst='CCD',
         pattern='Lick-d{} (Fischer+ 2014)'),
    dict(key='rvbank20', reference='Trifonov et al. 2020 (HARPS-RVBank)',
         catalogue='J/A+A/636/A74', stars='J/A+A/636/A74/list', name='Name',
         rvs='J/A+A/636/A74/rvbank', link='Name', time='BJD',
         rv='DRVmlcnzp', err='e_DRVmlcnzp', offset=2400000.0,
         eras=(57174.5, 'HARPS03 (RVBank)', 'HARPS15 (RVBank)')),
]
#: the sources of published velocities, each fetched or not on its own: the
#: surveys (their keys) and the tables of the star's papers ('papers')
SOURCES: List[str] = [survey['key'] for survey in SURVEYS] + ['papers']


# =============================================================================
# VizieR and SIMBAD
# =============================================================================
def _rows(source: str, constraints: Dict[str, str], timeout: float
          ) -> List[Dict[str, str]]:
    """the rows of a VizieR table under some constraints, as dicts"""
    params = {'-source': source, '-out.max': 'unlimited'}
    params.update(constraints)
    url = literature.VIZIER + '?' + urllib.parse.urlencode(params, safe='/*')
    with urllib.request.urlopen(url, timeout=timeout) as response:
        root = ElementTree.fromstring(response.read())
    for elem in root.iter():
        if '}' in elem.tag:
            elem.tag = elem.tag.split('}', 1)[1]
    table = root.find('.//TABLE')
    if table is None:
        return []
    names = [field.get('name') for field in table.findall('FIELD')]
    return [dict(zip(names, [(td.text or '').strip()
                             for td in tr.findall('TD')]))
            for tr in table.iter('TR')]


def _number(values: Sequence[str]) -> np.ndarray:
    """numbers from the cells of a table, nan for an empty one"""
    out = []
    for val in values:
        try:
            out.append(float(val))
        except (TypeError, ValueError):
            out.append(np.nan)
    return np.asarray(out, dtype=float)


def survey_stars(survey: Dict[str, Any], refresh: bool = False,
                 timeout: float = 120.0) -> List[Tuple[str, float, float]]:
    """
    The stars of a survey (their names in it and their J2000 positions),
    asked of VizieR once and kept in CACHE

    :param survey: dict, one of SURVEYS
    :param refresh: bool, ask VizieR again
    :param timeout: float [s]

    :return: list of (name, ra, dec) [deg]
    """
    path = os.path.join(CACHE, f'{survey["key"]}_stars.json')
    if os.path.exists(path) and not refresh:
        with open(path) as handle:
            return [tuple(row) for row in json.load(handle)]
    rows = _rows(survey['stars'], {'-out': survey['name'],
                                   '-out.add': '_RAJ,_DEJ',
                                   '-oc.form': 'dec'}, timeout)
    out = []
    for row in rows:
        try:
            out.append((row[survey['name']], float(row['_RAJ2000']),
                        float(row['_DEJ2000'])))
        except (KeyError, TypeError, ValueError):
            continue
    os.makedirs(CACHE, exist_ok=True)
    with open(path + '.part', 'w') as handle:
        json.dump(out, handle)
    os.replace(path + '.part', path)
    return out


def refresh_lists(timeout: float = 120.0) -> None:
    """the lists of the surveys' stars and what is known of the
    catalogues, asked again (--refresh-archive)"""
    for survey in SURVEYS:
        try:
            survey_stars(survey, refresh=True, timeout=timeout)
        except (OSError, ValueError, ElementTree.ParseError) as err:
            log(f'published: the stars of {survey["reference"]} could not '
                f'be asked ({err})', 'warn')
    path = os.path.join(CACHE, 'catalogues.json')
    if os.path.exists(path):
        os.remove(path)


def nearest_star(stars: Sequence[Tuple[str, float, float]], ra: float,
                 dec: float, radius: float = RADIUS) -> Optional[str]:
    """the name of the star of a list nearest a position, within a
    radius [arcsec], or None"""
    if not stars:
        return None
    pos = np.radians(np.array([[row[1], row[2]] for row in stars]))
    ra0, dec0 = np.radians(ra), np.radians(dec)
    cosd = (np.sin(dec0) * np.sin(pos[:, 1]) + np.cos(dec0)
            * np.cos(pos[:, 1]) * np.cos(pos[:, 0] - ra0))
    sep = np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0))) * 3600.0
    best = int(np.argmin(sep))
    return stars[best][0] if sep[best] <= radius else None


def survey_velocities(survey: Dict[str, Any], ra: float, dec: float,
                      timeout: float = 60.0
                      ) -> Tuple[Optional[RVData], str]:
    """
    The velocities of a star in a survey: its star in the survey's list
    (kept here), then its velocities by the survey's name of it

    :param survey: dict, one of SURVEYS
    :param ra: float, the right ascension of the star (J2000) [deg]
    :param dec: float, its declination [deg]
    :param timeout: float [s]

    :return: tuple, the series (None when the survey does not have the
             star) and what was found, in words
    """
    name = nearest_star(survey_stars(survey, timeout=max(timeout, 120.0)),
                        ra, dec)
    if name is None:
        return None, 'not in the survey (its list of stars)'
    cols = [survey['link'], survey['time'], survey['rv'], survey['err']]
    if survey.get('inst'):
        cols.append(survey['inst'])
    rows = _rows(survey['rvs'], {survey['link']: name,
                                 '-out': ','.join(cols)}, timeout)
    if not rows:
        return None, f'{name}: no velocity'
    time = _number([row.get(survey['time']) for row in rows]) - \
        survey['offset']
    rv = _number([row.get(survey['rv']) for row in rows])
    err = _number([row.get(survey['err']) for row in rows])
    if survey.get('inst'):
        raw = [row.get(survey['inst']) or '?' for row in rows]
        if 'names' in survey:
            inst = [survey['names'].get(val, f'{val} ({survey["key"]})')
                    for val in raw]
        else:
            inst = [survey['pattern'].format(val) for val in raw]
    elif survey.get('eras'):
        split, before, after = survey['eras']
        inst = [before if tt < split else after for tt in time]
    else:
        inst = [survey['label']] * len(rows)
    good = np.isfinite(time) & np.isfinite(rv) & np.isfinite(err) & (err > 0)
    if not np.any(good):
        return None, f'{name}: no usable velocity'
    data = RVData(time=time[good], rv=rv[good], err=err[good],
                  inst=np.array(inst)[good], name=survey['key'])
    return data, f'{name}: ' + ', '.join(
        f'{inst} {int(np.sum(data.inst == inst))}' for inst in data.instruments)


def star_bibcodes(main: str, timeout: float = 60.0) -> List[Tuple[str, int]]:
    """
    The papers SIMBAD lists for a star (by its main identifier), the latest
    first

    :return: list of (bibcode, year)
    """
    query = ('SELECT r.bibcode, r."year" FROM ref AS r JOIN has_ref AS h '
             'ON r.oidbib = h.oidbibref JOIN basic AS b ON h.oidref = b.oid '
             f"WHERE b.main_id = '{main.replace(chr(39), chr(39) * 2)}'")
    data = urllib.parse.urlencode(dict(request='doQuery', lang='adql',
                                       format='json', query=query)).encode()
    with urllib.request.urlopen(SIMBAD_TAP, data=data,
                                timeout=timeout) as response:
        answer = json.loads(response.read())
    out = [(str(row[0]), int(row[1] or 0)) for row in answer.get('data', [])]
    return sorted(set(out), key=lambda item: -item[1])


def catalogue_title(catalogue: str, timeout: float = 30.0) -> Optional[str]:
    """
    The authors and year of a VizieR catalogue ('Trifonov+ 2020', from the
    first line of its ReadMe), or None when the CDS does not have it
    """
    try:
        with urllib.request.urlopen(literature.CDSARC.format(catalogue),
                                    timeout=timeout) as response:
            first = response.readline().decode('utf-8', 'replace')
    except urllib.error.HTTPError as err:
        if err.code == 404:
            return None
        raise
    found = re.search(r'\(([^()]+?)\s*,\s*((?:19|20)\d\d)\)\s*$', first.strip())
    if found:
        return f'{found.group(1).strip()} {found.group(2)}'
    return catalogue


# =============================================================================
# One star: the surveys, then its papers
# =============================================================================
def _keep_precise(data: RVData) -> Optional[RVData]:
    """the instruments of a series whose median error is below MAX_ERROR"""
    keep = np.zeros(data.n, dtype=bool)
    for inst in data.instruments:
        sel = data.inst == inst
        if np.median(data.err[sel]) <= MAX_ERROR:
            keep |= sel
    if not np.any(keep):
        return None
    return data if np.all(keep) else data.select(keep)


def enough(data: Optional[RVData]) -> Optional[RVData]:
    """the instruments of a series with MIN_POINTS velocities or more"""
    if data is None:
        return None
    keep = np.zeros(data.n, dtype=bool)
    for inst in data.instruments:
        sel = data.inst == inst
        if np.sum(sel) >= MIN_POINTS:
            keep |= sel
    if not np.any(keep):
        return None
    return data if np.all(keep) else data.select(keep)


def family(inst: Any) -> str:
    """the spectrograph of an instrument's name ('HIRES-k (CLS)' HIRES,
    'HARPS03 (RVBank)' HARPS), '' when none is named"""
    text = str(inst).upper()
    return next((name for name in FAMILIES if re.search(
        rf'(?<![A-Z]){re.escape(name.upper())}', text)), '')


def _new(data: RVData, kept: List[Tuple[np.ndarray, np.ndarray]]
         ) -> Optional[RVData]:
    """the velocities of a series that are not spectra already kept: none
    within SAME_ANY of one, or within SAME_FAMILY of one of the same
    spectrograph"""
    if not kept:
        return data
    times = np.concatenate([part[0] for part in kept])
    fams = np.concatenate([part[1] for part in kept])
    mine = np.array([family(inst) for inst in data.inst])
    dist = np.abs(data.time[:, None] - times[None, :])
    same = (dist < SAME_ANY) | ((dist < SAME_FAMILY)
                                & (mine[:, None] == fams[None, :])
                                & (mine[:, None] != ''))
    fresh = ~np.any(same, axis=1)
    if not np.any(fresh):
        return None
    return data if np.all(fresh) else data.select(fresh)


def new_spectra(data: RVData, kept: Sequence[RVData]) -> Optional[RVData]:
    """
    The velocities of a series that are not spectra of other series (a
    spectrum published twice, or one that a file or an archive has): none
    within SAME_ANY of a velocity of any of them, or within SAME_FAMILY of
    one of the same spectrograph

    :param data: RVData, the series
    :param kept: list of RVData, the series it is set against

    :return: RVData or None (when every velocity is another's)
    """
    return _new(data, [(part.time, np.array([family(val)
                                             for val in part.inst]))
                       for part in kept])


def _tidy(data: RVData, title: str) -> RVData:
    """a paper's instruments: its own (their own offsets), and one when
    its 'instrument' column is not one (numbers, or a value per row)"""
    raw = [str(val) for val in data.inst]
    distinct = sorted(set(raw))
    numeric = sum(bool(re.fullmatch(r'[-+]?[\d.]+', val)) for val in distinct)
    if len(distinct) > 6 or numeric > len(distinct) / 2:
        named = family(' '.join(distinct)) or ''
        raw = [named or title] * data.n
    inst = np.array([val if title in val else f'{val} ({title})'
                     for val in raw])
    zero = np.array([data.zero_point.get(str(val), 0.0) for val in data.inst])
    return RVData(time=data.time, rv=data.rv + zero, err=data.err, inst=inst,
                  name=data.name)


def _write(data: RVData, path: str) -> None:
    """a series as a CSV (rjd, vrad, svrad, inst), its zero points back"""
    zero = np.array([data.zero_point.get(str(inst), 0.0)
                     for inst in data.inst])
    with open(path, 'w') as handle:
        handle.write('rjd,vrad,svrad,inst\n')
        for row in zip(data.time, data.rv + zero, data.err, data.inst):
            handle.write(f'{row[0]:.6f},{row[1]:.6f},{row[2]:.6f},{row[3]}\n')


def fetch(ident: Dict[str, Any], folder: str, refresh: bool = False,
          papers: bool = True, timeout: float = 60.0, workers: int = 8,
          sources: Optional[Sequence[str]] = None
          ) -> List[Dict[str, Any]]:
    """
    The published velocities of a star: the surveys, then the tables of
    its papers; each source kept as a CSV in a folder, with what was found
    (published.json), and read back from there unless refresh

    :param ident: dict, the star (koloa.archive.resolve: main, aliases, ra,
                  dec)
    :param folder: str, where the files are kept
    :param refresh: bool, ask VizieR and SIMBAD again
    :param papers: bool, look at the papers too (not the surveys alone)
    :param timeout: float [s]
    :param workers: int, the questions to the CDS at once
    :param sources: list of str or None, the sources asked (keys of
                    SOURCES: the surveys, and 'papers'), None for all

    :return: list of dict, per source: kind (survey or paper), key,
             reference, catalogue, n, instruments, file and note
    """
    if sources is None:
        sources = list(SOURCES)
    unknown = [key for key in sources if key not in SOURCES]
    if unknown:
        raise ValueError(f'no source {unknown} of published velocities: '
                         f'{", ".join(SOURCES)}')
    papers = papers and 'papers' in sources
    asked = [key for key in SOURCES if key in sources
             and (key != 'papers' or papers)]
    os.makedirs(folder, exist_ok=True)
    index = os.path.join(folder, 'published.json')
    # the sources asked last time (all of them for a folder of before
    #   they could be chosen): read back when they are the same
    before = os.path.join(folder, 'sources.json')
    if os.path.exists(index) and not refresh:
        last = list(SOURCES)
        if os.path.exists(before):
            with open(before) as handle:
                last = json.load(handle)
        if set(last) == set(asked):
            with open(index) as handle:
                return json.load(handle)
    names = [ident.get('main') or ''] + list(ident.get('aliases') or [])
    notes, kept = [], []

    def keep(data, note):
        """a source's velocities, its imprecise instruments left out: kept
        whole, with how many of its spectra the sources before it have
        (the release used is chosen later: koloa.datasets)"""
        if data is not None:
            data = _keep_precise(data)
            if data is None:
                note['note'] += f'; no instrument below {MAX_ERROR:.0f} m/s'
        if data is not None:
            fresh = _new(data, kept)
            note['shared'] = int(data.n - (fresh.n if fresh is not None
                                           else 0))
            if note['shared']:
                note['note'] += (f'; {note["shared"]} of its spectra in '
                                 f'the sources before it')
            kept.append((data.time, np.array([family(val)
                                              for val in data.inst])))
            # every release of a spectrum is on disk (a folder gathered
            #   before 2026-10-06 has the first one only)
            note['whole'] = True
            note['file'] = f'{note["key"]}.csv'
            note['n'] = int(data.n)
            note['instruments'] = {str(inst): int(np.sum(data.inst == inst))
                                   for inst in data.instruments}
            _write(data, os.path.join(folder, note['file']))
        notes.append(note)
        log(f'published, {note["reference"]}: {note["note"]}'
            + (f' ({note["n"]} kept)' if note['n'] else ''),
            'value' if note['n'] else 'info')

    # the surveys, by the position of the star
    for survey in SURVEYS:
        if survey['key'] not in asked:
            continue
        note = dict(kind='survey', key=survey['key'],
                    reference=survey['reference'],
                    catalogue=survey['catalogue'], n=0, instruments={},
                    file=None, note='')
        if ident.get('ra') is None:
            note['note'] = 'no position of the star'
            notes.append(note)
            continue
        try:
            data, note['note'] = survey_velocities(survey, ident['ra'],
                                                   ident['dec'], timeout)
        except (OSError, ValueError, ElementTree.ParseError) as err:
            data, note['note'] = None, f'VizieR could not be asked ({err})'
        keep(data, note)
    # the papers of the star: their VizieR catalogues, their tables of it
    if papers and ident.get('main'):
        try:
            codes = star_bibcodes(ident['main'], timeout)
        except (OSError, ValueError) as err:
            codes = []
            log(f'published: SIMBAD\'s bibliography could not be asked '
                f'({err})', 'warn')
        skip = {survey['catalogue'] for survey in SURVEYS} | NOT_PAPERS
        cats = []
        for code, year in codes:
            cat = literature.vizier_catalogue(code)
            if (cat and year >= MIN_YEAR and cat not in skip
                    and cat not in [item[1] for item in cats]):
                cats.append((code, cat))
        # what is known of each catalogue (kept for every star): the CDS
        #   does not have it, or has it without a table of velocities
        known = _catalogues()
        ask = [item for item in cats if item[1] not in known]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            answers = list(pool.map(
                lambda item: _safe(_describe, item[1], timeout), ask))
        for (code, cat), answer in zip(ask, answers):
            if answer is not None:
                known[cat] = answer
        _catalogues(known)
        found = [(code, cat, known[cat]['title']) for code, cat in cats
                 if known.get(cat, {}).get('rv')]
        log(f'published: {len(codes)} papers of {ident["main"]} in SIMBAD, '
            f'{len(cats)} with a catalogue VizieR could have, {len(found)} '
            f'with tables of velocities', 'info')

        def look(item):
            code, cat, title = item
            try:
                return literature.vizier_velocities(cat, names, title,
                                                    timeout, strict=True)
            except (OSError, ValueError, ElementTree.ParseError) as err:
                return None, f'VizieR could not be asked ({err})'
            except (KeyError, IndexError, TypeError) as err:
                # a table koloa does not read (a date that is not of an
                #   observation, say): left out, not the others
                return None, f'a table not read ({type(err).__name__}: {err})'
        with ThreadPoolExecutor(max_workers=max(1, workers // 2)) as pool:
            answers = list(pool.map(look, found))
        for (code, cat, title), (data, words) in zip(found, answers):
            if data is None:
                continue
            # each paper's instruments its own (their own offsets)
            data = _tidy(data, title)
            key = re.sub(r'[^A-Za-z0-9]+', '_', cat).strip('_')
            keep(data, dict(kind='paper', key=key, reference=title,
                            bibcode=code, catalogue=cat, n=0,
                            instruments={}, file=None, note=words))
    with open(index, 'w') as handle:
        json.dump(notes, handle, indent=1)
    with open(before, 'w') as handle:
        json.dump(asked, handle)
    return notes


def _catalogues(update: Optional[Dict[str, Any]] = None
                ) -> Dict[str, Any]:
    """what is known of the catalogues of papers (CACHE/catalogues.json):
    read, or written when given"""
    path = os.path.join(CACHE, 'catalogues.json')
    if update is not None:
        os.makedirs(CACHE, exist_ok=True)
        # several stars may be gathered at once: a part of its own for
        #   each writer, put in place whole
        import threading
        part = f'{path}.{os.getpid()}.{threading.get_ident()}.part'
        with open(part, 'w') as handle:
            json.dump(update, handle)
        os.replace(part, path)
        return update
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except ValueError:  # asked again, rather than stopped by a bad file
        return {}


def _describe(catalogue: str, timeout: float = 60.0) -> Dict[str, Any]:
    """a catalogue: its authors and year (None when the CDS does not
    have it), and whether it has a table of velocities"""
    title = catalogue_title(catalogue, timeout)
    if title is None:
        return dict(title=None, rv=False)
    root = literature._votable({'-source': catalogue, '-meta.all': ''},
                               timeout)
    rv = any(literature._columns(literature._fields(table)) is not None
             for table in root.iter('TABLE')
             if table.get('name', '').startswith(catalogue + '/'))
    return dict(title=title, rv=bool(rv))


def _safe(func, *args):
    """a question to the CDS that may fail: None then"""
    try:
        return func(*args)
    except (OSError, ValueError, ElementTree.ParseError):
        return None


def origins(folder: str) -> Dict[str, str]:
    """
    Where each instrument of the published velocities kept in a folder
    came from: its survey or its paper (two datasets of one source are not
    the same spectra; two sources can publish the same ones)

    :return: dict, {instrument: reference}
    """
    index = os.path.join(folder, 'published.json')
    if not os.path.exists(index):
        return {}
    with open(index) as handle:
        notes = json.load(handle)
    return {str(inst): str(note['reference']) for note in notes
            if note.get('file') for inst in note.get('instruments') or {}}


def whole(folder: str) -> bool:
    """whether a folder of published velocities has every source whole
    (gathered since each release of a spectrum is kept: before, a spectrum
    was kept from the first source that had it)"""
    index = os.path.join(folder, 'published.json')
    if not os.path.exists(index):
        return True
    with open(index) as handle:
        notes = json.load(handle)
    return all(note.get('whole') for note in notes if note.get('file'))


def load(folder: str) -> Optional[RVData]:
    """
    The published velocities kept in a folder (fetch), as one series

    :return: RVData, or None when there are none
    """
    index = os.path.join(folder, 'published.json')
    if not os.path.exists(index):
        return None
    with open(index) as handle:
        notes = json.load(handle)
    parts = [RVData.from_csv(os.path.join(folder, note['file']), inst='inst',
                             name=note['reference'])
             for note in notes if note.get('file')
             and os.path.exists(os.path.join(folder, note['file']))]
    if not parts:
        return None
    return parts[0] if len(parts) == 1 else merge(parts, name='published')


# =============================================================================
# End of code
# =============================================================================

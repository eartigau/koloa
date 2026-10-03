#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Velocities from the literature: given, or found on VizieR

- read(source): a series from a file of published velocities (a VizieR
  .dat without a header: time, velocity, error and an optional instrument;
  or a csv or .rdb with named columns, as RVData.from_csv reads them), or
  an RVData; times in BJD, BJD - 2400000, BJD - 2450000 or MJD become
  koloa's BJD - 2400000, velocities in km/s become m/s;
- vizier_catalogue(bibcode): the VizieR catalogue of the tables of a paper
  (J/A+A/680/A28 for 2023A&A...680A..28G: VizieR names the tables of
  journals after their volume and page);
- vizier_velocities(catalogue, names): the tables of velocities of a
  catalogue that are those of the star (by the name of the table, its
  description, or a column of names), as a series;
- fetch(known, names, folder): the velocities published with every
  solution of the known planets of a star (koloa.archive), kept as CSV
  files in a folder and read back.

A published instrument that shares its name with one of the series keeps
its own offset (renamed 'SPIRou (Author et al. 2025)'), and a published
velocity within a minute of an exposure of the series is the same spectrum
and is left out.

Created on 2026-09-30

@author: artigau
"""
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ElementTree
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from koloa.data import ERR_NAMES, RV_NAMES, TIME_NAMES, RVData, merge
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
VIZIER = 'https://vizier.cds.unistra.fr/viz-bin/votable'
#: the ReadMe of every catalogue (asked first: VizieR takes a catalogue it
#: does not have for a pattern, slowly, and answers with others)
CDSARC = 'https://cdsarc.cds.unistra.fr/ftp/{}/ReadMe'
#: VizieR's names of the journals of the ADS bibcodes
JOURNALS = {'A&A': 'A+A', 'MNRAS': 'MNRAS', 'AJ': 'AJ', 'ApJ': 'ApJ',
            'ApJS': 'ApJS', 'ApJL': 'ApJ', 'PASP': 'PASP', 'A&AS': 'A+AS',
            'AN': 'AN', 'PASJ': 'PASJ', 'RNAAS': 'RNAAS'}
#: names of an instrument column
INST_NAMES = ['inst', 'instrument', 'instr', 'ins', 'spectrograph',
              'spectro', 'tel', 'telescope', 'dataset']
#: names of a velocity column, in order of preference
VIZIER_RV = ['rv', 'rvel', 'vr', 'vrad', 'hrv', 'rvc', 'rvcorr', 'vel',
             'rvs', 'rvobs']
#: names of a column of star names
VIZIER_STAR = ['name', 'star', 'target', 'object', 'id', 'simbad', 'hd',
               'gj', 'tic', 'karmn']
#: velocity columns that are not velocities of the star
NOT_RV = ('bis', 'fwhm', 'ccf', 'contr', 'crx', 'dlw', 'span', 'drift')
#: a published velocity within this of an exposure of the series is the same
#: spectrum [days]
SAME = 1.0 / 1440
#: spectrographs a table without an instrument column may name in its
#: description (the longer names first: HARPS-N before HARPS)
SPECTROGRAPHS = ('HARPS-N', 'HARPS', 'CARMENES', 'ESPRESSO', 'SPIRou',
                 'NIRPS', 'HIRES', 'IRD', 'MAROON-X', 'SOPHIE', 'CORALIE',
                 'NEID', 'EXPRES', 'KPF', 'HPF', 'PFS', 'APF', 'UVES',
                 'iSHELL', 'ELODIE', 'FEROS')


# =============================================================================
# Times and units
# =============================================================================
def to_rjd(time: np.ndarray, hint: str = '') -> np.ndarray:
    """
    Times in BJD - 2400000, from BJD (or JD), from BJD minus another offset
    (read from the hint: 'BJD-2450000'), or from MJD

    :param time: np.ndarray, the times
    :param hint: str, the name, unit and description of the column

    :return: np.ndarray
    """
    time = np.asarray(time, dtype=float)
    mid = np.nanmedian(time)
    hint = hint or ''
    if mid > 2.3e6:
        return time - 2400000.0
    offset = re.search(r'24(\d{5})', hint.replace(' ', ''))
    if offset:
        return time + float('24' + offset.group(1)) - 2400000.0
    if 3.0e4 < mid < 1.0e5:
        # MJD = JD - 2400000.5
        return time + 0.5 if 'mjd' in hint.lower() else time
    # a short time with no offset given: BJD - 2450000, the usual one
    log(f'literature: times near {mid:.0f} with no offset given ({hint}); '
        f'taken as BJD - 2450000', 'warn')
    return time + 50000.0


def _velocity_scale(unit: str) -> float:
    """m/s per unit of a velocity column"""
    return 1000.0 if (unit or '').lower().startswith('km') else 1.0


def _name_pattern(names: Sequence[str]):
    """a pattern that finds any of the names of a star in a text, whatever
    the spaces and dashes ('GJ 3988' in 'gj3988rv' or 'RVs of GJ 3988')"""
    alts = []
    for name in names:
        tokens = [tok for tok in re.split(r'[\s_\-]+', str(name).lower())
                  if tok]
        if not tokens or len(''.join(tokens)) < 4:
            continue
        alts.append(r'[\s_\-]*'.join(re.escape(tok) for tok in tokens))
    if not alts:
        return None
    return re.compile(r'(?<![a-z0-9])(?:' + '|'.join(alts) + r')(?![0-9])')


# =============================================================================
# Files
# =============================================================================
def read(source: Union[str, RVData], inst: Optional[str] = None,
         name: Optional[str] = None) -> RVData:
    """
    A series of published velocities

    :param source: str or RVData, a file (a VizieR .dat: time, velocity,
                   error and an optional instrument, without a header; or a
                   csv or .rdb with named columns), or a series
    :param inst: str or None, the instrument when the file has no column
                 of it (the name of the file when None)
    :param name: str or None, the name of the series

    :return: RVData, the times in BJD - 2400000
    """
    if isinstance(source, RVData):
        return source
    stem = os.path.splitext(os.path.basename(source))[0]
    with open(source) as handle:
        lines = [line for line in handle
                 if line.strip() and not line.lstrip().startswith('#')]
    if not lines:
        raise ValueError(f'{source} holds no velocity')
    first = re.split(r'[\s,;|]+', lines[0].strip().lower())
    known = set(TIME_NAMES + RV_NAMES + ERR_NAMES)
    if any(tok in known for tok in first):
        # named columns: RVData.from_csv, then the times to BJD - 2400000
        header = [tok.strip() for tok in re.split(
            '\t' if '\t' in lines[0] else ',', lines[0].strip())]
        icol = next((col for col in header if col.lower() in INST_NAMES),
                    None)
        data = RVData.from_csv(source, inst=icol, name=name or stem)
        tcol = next((col for col in header if col.lower() in TIME_NAMES), '')
        data.time = to_rjd(data.time, tcol)
    else:
        table = np.genfromtxt(source, dtype=None, encoding=None,
                              comments='#')
        cols = [table[field] for field in table.dtype.names]
        if len(cols) < 3:
            raise ValueError(f'{source}: needs time, velocity and error')
        texts = [col for col in cols[3:] if col.dtype.kind in 'US']
        data = RVData(time=to_rjd(cols[0]), rv=np.asarray(cols[1], float),
                      err=np.asarray(cols[2], float),
                      inst=(np.char.strip(texts[0].astype(str)) if texts
                            else None),
                      name=name or stem)
    if data.instruments == ['inst']:
        label = inst or stem
        data.inst = np.array([label] * data.n)
        data.zero_point = {label: data.zero_point.get('inst', 0.0)}
    return data


# =============================================================================
# VizieR
# =============================================================================
def vizier_catalogue(code: str) -> Optional[str]:
    """
    The VizieR catalogue of the tables of a paper, from its bibcode (VizieR
    names the catalogues of journals J/<journal>/<volume>/<page>)

    :param code: str, an ADS bibcode (19 characters)

    :return: str or None (a journal VizieR does not name this way)
    """
    if not code or len(code) != 19:
        return None
    # YYYY JJJJJ VVVV Q PPPP A: the qualifier Q (A&A's article letter, L
    #   for a letter) is part of VizieR's page
    journal = code[4:9].strip('.')
    volume = code[9:13].strip('.')
    page = code[13].strip('.') + code[14:18].strip('.')
    if journal not in JOURNALS or not volume or not page:
        return None
    return f'J/{JOURNALS[journal]}/{volume}/{page}'


def has_catalogue(catalogue: str, timeout: float = 30.0) -> bool:
    """whether the CDS has a catalogue (its ReadMe answers)"""
    request = urllib.request.Request(CDSARC.format(catalogue), method='HEAD')
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status == 200
    except urllib.error.HTTPError as err:
        if err.code == 404:
            return False
        raise


def _votable(params: Dict[str, str], timeout: float) -> ElementTree.Element:
    """a VizieR VOTable, its namespaces taken out"""
    url = VIZIER + '?' + urllib.parse.urlencode(params, safe='/*')
    with urllib.request.urlopen(url, timeout=timeout) as response:
        text = response.read()
    root = ElementTree.fromstring(text)
    for elem in root.iter():
        if '}' in elem.tag:
            elem.tag = elem.tag.split('}', 1)[1]
    return root


def _fields(table: ElementTree.Element) -> List[Dict[str, str]]:
    """the columns of a VOTable table"""
    out = []
    for field in table.findall('FIELD'):
        desc = field.find('DESCRIPTION')
        out.append(dict(name=field.get('name', ''), ucd=field.get('ucd', ''),
                        unit=field.get('unit', ''),
                        description=(desc.text or '') if desc is not None
                        else ''))
    return out


def _columns(fields: List[Dict[str, str]]) -> Optional[Dict[str, Any]]:
    """the time, velocity, error, instrument and star columns of a table,
    or None when it has no velocities"""
    def low(field):
        return field['name'].lower()

    times = [fd for fd in fields if 'time.epoch' in fd['ucd']
             or low(fd) in ('bjd', 'jd', 'hjd', 'mjd', 'rjd', 'time')]
    vels = [fd for fd in fields
            if ('phys.veloc' in fd['ucd'] or 'spect.dopplerveloc'
                in fd['ucd'].lower() or low(fd) in VIZIER_RV)
            and 'stat.error' not in fd['ucd'] and not low(fd).startswith('e_')
            and not any(word in low(fd) for word in NOT_RV)]
    if not times or not vels:
        return None
    order = {name: rank for rank, name in enumerate(VIZIER_RV)}
    vel = min(vels, key=lambda fd: order.get(low(fd), len(order)))
    names = [fd['name'] for fd in fields]
    err = next((fd for fd in fields if fd['name'] == f'e_{vel["name"]}'),
               None)
    if err is None:
        # the error column right after the velocity
        ivel = names.index(vel['name'])
        err = next((fd for fd in fields[ivel + 1:ivel + 3]
                    if 'stat.error' in fd['ucd']), None)
    if err is None:
        return None
    tpref = ('bjd', 'jd', 'hjd', 'time', 'rjd', 'mjd')
    tcol = min(times, key=lambda fd: next(
        (rank for rank, word in enumerate(tpref) if word in low(fd)),
        len(tpref)))
    inst = next((fd for fd in fields if 'instr' in fd['ucd']
                 or low(fd) in INST_NAMES), None)
    star = next((fd for fd in fields if 'meta.id' in fd['ucd']
                 and 'instr' not in fd['ucd'] and low(fd) in VIZIER_STAR),
                None)
    return dict(time=tcol, rv=vel, err=err, inst=inst, star=star)


def vizier_velocities(catalogue: str, names: Sequence[str],
                      reference: str = '', timeout: float = 60.0,
                      strict: bool = False
                      ) -> Tuple[Optional[RVData], str]:
    """
    The published velocities of a star in a VizieR catalogue

    A table is the star's when its name or description names the star, or
    when a column of names does (its rows are then selected), or when it is
    the only table of velocities of the catalogue (strict: and the
    catalogue's description names the star; a paper that only mentions the
    star has the velocities of another).

    :param catalogue: str, e.g. J/A+A/680/A28
    :param names: list of str, the names of the star (SIMBAD's aliases)
    :param reference: str, the paper (for the instruments' names, when one
                      is not given)
    :param timeout: float [s]
    :param strict: bool, the only table of velocities of a catalogue is the
                   star's only when the catalogue names it (a paper found in
                   the star's bibliography, not one of its planets)

    :return: tuple, the series (None when there is none) and what was found,
             in words
    """
    if not has_catalogue(catalogue, timeout):
        return None, f'VizieR has no catalogue {catalogue}'
    root = _votable({'-source': catalogue, '-meta.all': ''}, timeout)
    # the description of a table is that of its RESOURCE (and its own)
    tables = [(res, table) for res in root.iter('RESOURCE')
              for table in res.findall('TABLE')
              if table.get('name', '').startswith(catalogue + '/')]
    if not tables:
        return None, f'{catalogue} has no table'
    pattern = _name_pattern(names)
    candidates = []
    for res, table in tables:
        cols = _columns(_fields(table))
        if cols is None:
            continue
        own = table.find('DESCRIPTION')
        shared = res.find('DESCRIPTION')
        own = (own.text or '') if own is not None else ''
        words = ' '.join([(shared.text or '') if shared is not None
                          else '', own])
        # the star in the name of the table or its own description (a
        #   description shared by the tables of several stars names them all)
        text = (table.get('name', '').split('/')[-1] + ' ' + own).lower()
        named = bool(pattern and pattern.search(text))
        # no instrument column: the one spectrograph its description names
        spectro = [name for name in SPECTROGRAPHS
                   if re.search(rf'(?<![\w-]){re.escape(name)}(?![\w-])',
                                words, re.IGNORECASE)]
        cols['spectrograph'] = spectro[0] if len(spectro) == 1 else None
        cols['catalogue_named'] = bool(pattern and pattern.search(
            words.lower()))
        candidates.append((table.get('name'), cols, named))
    if not candidates:
        return None, f'{catalogue} has no table of velocities'
    chosen = [cand for cand in candidates if cand[2] or cand[1]['star']]
    if not chosen and len(candidates) == 1 and (
            not strict or candidates[0][1]['catalogue_named']):
        chosen = candidates
    if not chosen:
        return None, (f'{catalogue}: {len(candidates)} tables of velocities, '
                      f'none of this star')
    series, words = [], []
    for tname, cols, named in chosen:
        table = _votable({'-source': tname, '-out.all': '',
                          '-out.max': 'unlimited'}, timeout).find('.//TABLE')
        fields = [fd['name'] for fd in _fields(table)]
        rows = [[(td.text or '').strip() for td in tr.findall('TD')]
                for tr in table.iter('TR')]
        if not rows:
            continue
        cells = {name: np.array([row[it] if it < len(row) else ''
                                 for row in rows])
                 for it, name in enumerate(fields)}
        keep = np.ones(len(rows), dtype=bool)
        if cols['star'] is not None and not named:
            keep = np.array([bool(pattern and pattern.search(val.lower()))
                             for val in cells[cols['star']['name']]])
        if not np.any(keep):
            continue

        def number(field):
            vals = np.array([float(val) if val else np.nan
                             for val in cells[field['name']]])
            return vals[keep]
        time = to_rjd(number(cols['time']), ' '.join(
            cols['time'][key] for key in ('name', 'unit', 'description')))
        rv = number(cols['rv']) * _velocity_scale(cols['rv']['unit'])
        err = number(cols['err']) * _velocity_scale(
            cols['err']['unit'] or cols['rv']['unit'])
        label = cols['spectrograph'] or reference or catalogue
        inst = (cells[cols['inst']['name']][keep] if cols['inst'] is not None
                else np.array([label] * int(np.sum(keep))))
        good = np.isfinite(time) & np.isfinite(rv) & np.isfinite(err) \
            & (err > 0)
        if not np.any(good):
            continue
        series.append(RVData(time=time[good], rv=rv[good], err=err[good],
                             inst=np.array([str(val) or reference
                                            for val in inst[good]]),
                             name=tname))
        words.append(f'{tname.split("/")[-1]} ({int(np.sum(good))})')
    if not series:
        return None, f'{catalogue}: no velocity of this star'
    data = series[0] if len(series) == 1 else merge(series, name=catalogue)
    return data, f'{catalogue}: ' + ', '.join(words)


# =============================================================================
# The velocities published with the known planets
# =============================================================================
def _papers(known: Dict[str, Any]) -> List[Dict[str, str]]:
    """every paper of a solution of the known planets, once"""
    from koloa.archive import bibcode
    out, seen = [], set()
    for pl in known.get('planets', []):
        refs = list(pl.get('solutions') or [])
        refs.append(dict(reference=pl.get('reference'),
                         reference_url=pl.get('reference_url')))
        for sol in refs:
            code = bibcode(sol.get('reference_url') or '')
            if code and code not in seen:
                seen.add(code)
                out.append(dict(bibcode=code,
                                reference=sol.get('reference') or code))
    return out


def fetch(known: Dict[str, Any], names: Sequence[str], folder: str,
          refresh: bool = False, timeout: float = 60.0
          ) -> Tuple[List[RVData], List[Dict[str, Any]]]:
    """
    The velocities published with every solution of the known planets of a
    star, from VizieR, kept in a folder (a CSV per paper, and what was
    found in literature.json) and read back from it

    :param known: dict, from koloa.archive.known_planets
    :param names: list of str, the names of the star
    :param folder: str, where the files are kept
    :param refresh: bool, ask VizieR even when the files are there
    :param timeout: float [s]

    :return: tuple, a series per paper that has velocities of the star, and
             for every paper: reference, bibcode, catalogue, n, instruments,
             file and note
    """
    os.makedirs(folder, exist_ok=True)
    index = os.path.join(folder, 'literature.json')
    notes: List[Dict[str, Any]] = []
    if os.path.exists(index) and not refresh:
        with open(index) as handle:
            notes = json.load(handle)
    done = {note['bibcode'] for note in notes}
    for paper in _papers(known):
        if paper['bibcode'] in done:
            continue
        note = dict(paper, catalogue=vizier_catalogue(paper['bibcode']),
                    n=0, instruments={}, file=None)
        if note['catalogue'] is None:
            note['note'] = 'not a journal whose tables VizieR names'
        else:
            try:
                data, note['note'] = vizier_velocities(
                    note['catalogue'], names, paper['reference'], timeout)
            except (OSError, ValueError, ElementTree.ParseError) as err:
                data, note['note'] = None, f'VizieR could not be asked ({err})'
            if data is not None:
                safe = note['catalogue'].replace('/', '_').replace('+', '')
                note['file'] = f'{safe}.csv'
                note['n'] = int(data.n)
                note['instruments'] = {str(inst): int(np.sum(data.inst == inst))
                                       for inst in data.instruments}
                zero = np.array([data.zero_point.get(str(inst), 0.0)
                                 for inst in data.inst])
                with open(os.path.join(folder, note['file']), 'w') as handle:
                    handle.write('rjd,vrad,svrad,inst\n')
                    for row in zip(data.time, data.rv + zero, data.err,
                                   data.inst):
                        handle.write(f'{row[0]:.6f},{row[1]:.6f},'
                                     f'{row[2]:.6f},{row[3]}\n')
        log(f'literature, {paper["reference"]}: {note["note"]}',
            'value' if note['n'] else 'warn')
        notes.append(note)
    with open(index, 'w') as handle:
        json.dump(notes, handle, indent=1)
    series = []
    for note in notes:
        if note.get('file'):
            data = RVData.from_csv(os.path.join(folder, note['file']),
                                   inst='inst', name=note['reference'])
            series.append(data)
    return series, notes


def add(data: RVData, others: Sequence[RVData], labels: Sequence[str]
        ) -> Tuple[RVData, List[Dict[str, Any]]]:
    """
    A series with published ones merged in: an exposure within a minute of
    one of the series is left out (the same spectrum), and an instrument
    already in the series is renamed after its paper (its own offset)

    :param data: RVData, the series
    :param others: list of RVData, the published series
    :param labels: list of str, their papers (or files)

    :return: tuple, the merged series, and per published series: label, n
             (kept), same (left out), renamed ({old: new})
    """
    out, info = data, []
    for other, label in zip(others, labels):
        near = np.min(np.abs(other.time[:, None] - out.time[None, :]),
                      axis=1)
        same = near < SAME
        kept = other.select(~same) if np.any(same) else other
        renamed = {}
        if kept.n:
            insts = np.asarray(kept.inst).astype(object)
            for inst in kept.instruments:
                if inst in out.instruments:
                    new = f'{inst} ({label})'
                    renamed[str(inst)] = new
                    insts[insts == inst] = new
            if renamed:
                zero = np.array([kept.zero_point.get(str(inst), 0.0)
                                 for inst in kept.inst])
                kept = RVData(time=kept.time, rv=kept.rv + zero, err=kept.err,
                              inst=insts.astype(str), name=kept.name)
            out = merge([out, kept], name=data.name)
        info.append(dict(label=label, n=int(kept.n), same=int(np.sum(same)),
                         renamed=renamed))
    return out, info


# =============================================================================
# End of code
# =============================================================================

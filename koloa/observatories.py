#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The velocities of the pipelines of other spectrographs, from the archives
of their observatories, for koloa (HARPS-N at the TNG is koloa.tng).

Each archive is asked with no account, for what is public, and gives the
velocity its pipeline made of each exposure, in a table or in the header
of a small file:

- eso(): the archive of ESO. HARPS from its catalogue of velocities
  (Barbieri 2023: the pipeline's, 2003 to 2023, one question), ESPRESSO,
  NIRPS and HARPS since 2026 from the header of the CCF file of each
  exposure (120 to 310 kB, fetched whole: the velocity is at the end of
  a long header). The files of ESPRESSO are of the pipeline of their
  day: those of espdr 3 are 24 m/s from those of espdr 2 for Barnard's
  star, and a dataset of their own;
- neid(): NEID at the WIYN, a table of the archive of NExScI;
- sophie(): SOPHIE at the 1.93 m of the OHP, a table of its archive
  (before and after its octagonal fibres of June 2011, SOPHIE+; its two
  resolutions apart);
- elodie(): ELODIE, its predecessor (1994 to 2006), a table of its
  archive, the exposures with the thorium lamp beside the star. The table
  has no error: it is that of Perrier et al. (2003, A&A 410, 1039), the
  instrumental error of ELODIE, 6 m/s, and the photon noise added in
  quadrature, here 311 m/s over the signal-to-noise the table gives (so
  that 70 Vir has their 2.8 m/s); the mean errors of four stars of Naef
  et al. (2003) come back within 1.3 m/s;
- cfht(): SPIRou at the CFHT, from the headers of the CCF files of APERO
  at the CADC.

A star is found by where the telescope was pointed, as in koloa.tng: the
object names are the observers'.

Each dataset has the indicators of its CCF where the archive gives them,
as koloa keeps those of DACE (koloa.tng.ccf_indicators): its width, its
contrast and its bisector span (the width alone for ELODIE; none for
NEID, whose table has none).

Two of them are asked only when they are named (koloa.published
.OPTIONAL). ESO has what DACE has of HARPS, ESPRESSO and NIRPS: the same
pipeline and the same spectra (for GJ 436, 169, 23 and 99 velocities
from each), so it is of use where DACE cannot be reached, or to see the
velocities of a pipeline as its archive has them. SPIRou at the CFHT:
the CADC answers a header in 2.5 s, half an hour for the 1800 exposures
of a star.

Two things differ from one archive to the next, and each dataset says
them (the dict that goes with it):

- perspective: whether its velocities still hold the perspective
  acceleration of the star (koloa.secular). Measured on Barnard's star
  (4.54 m/s/yr): the pipelines of ESO do (HARPS 4.0 to 4.2 m/s/yr,
  NIRPS 4.0 +- 0.2), NEID's does not (-0.16 +- 0.11). APERO is told how
  the star moves and does not either: its CCF velocities of SPIRou drift
  by 2.51 +- 0.06 m/s/yr over 1761 exposures, between the 1.2 of its
  line-by-line velocities and far from the 5.7 they would show with it.
  SOPHIE and ELODIE are taken to hold it, as the pipelines of their
  family do (SOPHIE+ drifts by 2.4 m/s/yr for GJ 411, where 1.35 is
  expected: its zero point moves too, and does not settle it);
- defer: the velocity of a CCF in the near infrared is a poor one (SPIRou
  here, like NIRPS on DACE): a file given of the same spectrograph (LBL)
  is believed first, and the dataset is left out altogether when there is
  one.

These are velocities to look at a star with, not a release: the archive
of the OHP asks that what it serves be not handed on.

Created on 2026-10-13

@author: artigau
"""
import csv
import io
import json
import os
import re
import time as _time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from koloa import tng
from koloa.data import RVData
from koloa.log import log
from koloa.paths import cache

# =============================================================================
# Define variables
# =============================================================================
#: the services
ESO_CAT = 'https://archive.eso.org/tap_cat/sync'
ESO_OBS = 'https://archive.eso.org/tap_obs/sync'
ESO_FILE = 'https://dataportal.eso.org/dataPortal/file/'
NEID_TAP = 'https://neid.ipac.caltech.edu/TAP/sync'
SOPHIE_URL = 'http://atlas.obs-hp.fr/sophie/sophie.cgi'
ELODIE_URL = 'http://atlas.obs-hp.fr/elodie/E.cgi'
CADC_TAP = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync'
CADC_FILE = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/raven/files/'
#: how far from where a star was the telescope may have been pointed
#: [arcsec] (the coordinates of a night, as the telescope had them), and
#: from its place in a catalogue of positions
POINTED = 15.0
LISTED = 5.0
#: the days the instruments changed [BJD - 2400000]: the fibres of HARPS
#: (2015), the fibre link of ESPRESSO (2019), the octagonal fibres of
#: SOPHIE (2011)
HARPS_2015, ESPRESSO_2019, SOPHIE_2011 = 57174.5, 58660.0, 55730.0
#: the error of an ELODIE velocity [m/s]: its instrumental error, and its
#: photon noise at a signal-to-noise of one (Perrier et al. 2003)
ELODIE_FLOOR, ELODIE_PHOTON = 6.0, 311.0
#: the files fetched at once
WORKERS = 4
#: where the list of the pointings of NEID is kept, and for how long
#: [days]
NEID_POINTINGS = cache('observatories', 'neid_pointings.json')
MAX_AGE = 30.0


# =============================================================================
# Define functions
# =============================================================================
def _http(url: str, data: Optional[Dict[str, Any]] = None,
          where: str = 'archive', timeout: float = 60.0,
          tries: int = 2) -> bytes:
    """a page of an archive, asked again once when it stumbles (stays
    silent for timeout seconds, or the connection fails)"""
    body = urllib.parse.urlencode(data).encode() if data else None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(
                    url, data=body, headers={'User-Agent': 'koloa'}),
                    timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as err:
            if err.code in (401, 403):
                raise RuntimeError(f'{where}: not public ({err.code})'
                                   ) from err
            last: Exception = err
        except (urllib.error.URLError, TimeoutError, OSError) as err:
            last = err
        if attempt + 1 < tries:
            _time.sleep(3.0)
    raise RuntimeError(f'the {where} could not be reached ({last})')


def _tap(url: str, query: str, where: str, fmt: str = 'csv'
         ) -> List[Dict[str, Any]]:
    """a question to a TAP service (ADQL): its rows"""
    raw = _http(url, dict(REQUEST='doQuery', LANG='ADQL', FORMAT=fmt,
                          MAXREC=200000, QUERY=query), where).decode(
        'utf8', 'replace')
    if fmt == 'json':
        try:
            held = json.loads(raw)
        except ValueError as err:
            raise RuntimeError(f'the {where} could not be reached (it '
                               f'answered an error)') from err
        names = [col['name'] for col in held['metadata']]
        return [dict(zip(names, row)) for row in held['data']]
    if raw.lstrip().startswith('<'):
        raise RuntimeError(f'the {where} could not be reached (it answered '
                           f'an error)')
    return list(csv.DictReader(io.StringIO(raw)))


def _table(url: str, where: str) -> List[Dict[str, str]]:
    """a table of the archives of the OHP, asked as CSV: the lines that
    begin with # are about the question"""
    raw = _http(url, None, where).decode('utf8', 'replace')
    lines = [line for line in raw.splitlines()
             if line.strip() and not line.startswith('#')]
    return [{key.strip(): str(val).strip() for key, val in row.items()
             if key} for row in csv.DictReader(lines)]


def _num(value: Any) -> float:
    """a float, nan when it is not one"""
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _year(rjd: Any) -> np.ndarray:
    """the year of a date [BJD - 2400000]"""
    return 2000.0 + (np.asarray(rjd, float) - 51544.5) / 365.25


def _reach(star: Dict[str, Any], margin: float = POINTED) -> float:
    """how far a star may be from its place of 2000 by now, and a margin
    [deg]"""
    return (margin + (tng._now() - tng.EPOCH) / 1000.0 * float(np.hypot(
        float(star.get('pmra') or 0.0), float(star.get('pmdec') or 0.0)))
            ) / 3600.0


def _pointed(ra: Any, dec: Any, rjd: Any, star: Dict[str, Any],
             radius: float = POINTED) -> np.ndarray:
    """whether the telescope was pointed at a star: within radius of
    where it was that night, or at one of the years of a catalogue"""
    ra, dec = np.atleast_1d(np.asarray(ra, float)), np.atleast_1d(
        np.asarray(dec, float))
    years = np.column_stack([np.full(len(ra), year) for year in tng.EPOCHS]
                            + [np.atleast_1d(_year(rjd))])
    return tng.away(ra, dec, star, years) < radius


def _part(label: str, key: str, reference: str, time: Any, rv: Any,
          err: Any, names: Optional[Sequence[str]] = None,
          name: Optional[str] = None, perspective: bool = True,
          defer: bool = False, said: str = '', fwhm: Any = None,
          contrast: Any = None, bis: Any = None
          ) -> Optional[Tuple[RVData, Dict[str, Any]]]:
    """
    A dataset of an archive: its exposures with a velocity and an error,
    not those under a name whose velocities are of another star, nor one
    far from the others (as koloa.tng.series)

    :param time: [BJD - 2400000]; rv and err [m/s]
    :param names: the name the observer gave each, or None
    :param fwhm: the widths of the CCFs [km/s], contrast [%] and bis (their
                 bisector spans [km/s]), or None: the indicators of the
                 dataset (koloa.tng.ccf_indicators)

    :return: (RVData, dict: key, label, reference, n, left, perspective,
             defer), or None when nothing is left
    """
    time, rv, err = (np.asarray(val, float) for val in (time, rv, err))
    good = np.isfinite(time + rv + err) & (err > 0)
    left = [said] if said else []
    if np.sum(~good):
        left.append(f'{int(np.sum(~good))} with no velocity or no error')
    if names is not None and np.any(good):
        names = np.asarray(names, dtype=str)
        kinds, counts = np.unique(names[good], return_counts=True)
        main = np.median(rv[good & (names == kinds[np.argmax(counts)])])
        for kind in kinds:
            sel = good & (names == kind)
            if abs(np.median(rv[sel]) - main) > tng.OTHER * 1000.0:
                left.append(f'{int(np.sum(sel))} as {kind}, of another '
                            f'star by their velocity')
                good &= ~sel
    if np.any(good):
        far = good & (np.abs(rv - np.median(rv[good])) > tng.FAR * 1000.0)
        if np.any(far):
            left.append(f'{int(np.sum(far))} more than {tng.FAR:g} km/s '
                        f'from the others')
        good &= ~far
    if not np.any(good):
        return None
    order = np.argsort(time[good])
    more = tng.ccf_indicators(err[good][order], *(
        None if val is None else np.asarray(val, float)[good][order]
        for val in (fwhm, contrast, bis)))
    data = RVData(time[good][order], rv[good][order], err[good][order],
                  inst=np.full(int(np.sum(good)), label), name=name,
                  indicators=more)
    return data, dict(key=key, label=label, reference=reference,
                      n=int(data.n), left='; '.join(left),
                      perspective=bool(perspective), defer=bool(defer))


def _most(values: Sequence[Any]) -> Any:
    """the value most of them have"""
    kinds, counts = np.unique(np.asarray(values, dtype=str),
                              return_counts=True)
    return str(kinds[np.argmax(counts)])


def _told(parts: Sequence[Tuple[RVData, Dict[str, Any]]], none: str) -> str:
    """what was found, in words"""
    return ('; '.join(f'{about["label"]}: {about["n"]}'
                      for _, about in parts) or none)


def _headers(listed: Sequence[Dict[str, Any]], path: str,
             read: Callable[[Dict[str, Any]], Dict[str, Any]], words: str,
             refresh: bool = False, key: Optional[str] = None
             ) -> List[Dict[str, Any]]:
    """
    What the headers of the files of a star say, kept in a table (a CSV)
    so that a file is fetched once: only those that are public and not in
    it yet are (four at a time, with a bar)

    :param listed: list of dict: file (its name), public (bool), and what
                   read needs
    :param read: what reads a file, to a dict of its values (error when
                 it could not be)
    :param key: str or None, the archive in the table of koloa's GUI

    :return: list of dict, a row per public file read
    """
    known: Dict[str, Dict[str, Any]] = {}
    if os.path.exists(path) and not refresh:
        with open(path) as handle:
            known = {row['file']: row for row in csv.DictReader(handle)
                     if not row.get('error')}
    ask = [one for one in listed
           if one['public'] and one['file'] not in known]
    if ask:
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            for one, row in zip(ask, tng._progress(
                    pool.map(read, ask), len(ask), words, key)):
                known[one['file']] = dict(row, file=one['file'])
    rows = [known[one['file']] for one in listed if one['file'] in known]
    if rows:
        cols = list(dict.fromkeys(key for row in rows for key in row))
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=cols)
            writer.writeheader()
            writer.writerows(rows)
    failed = [row for row in rows if row.get('error')]
    if failed:
        log(f'{words}: {len(failed)} files could not be read '
            f'({failed[0]["error"]})', 'warn')
    return [row for row in rows if not row.get('error')]


# -----------------------------------------------------------------------------
# NEID
# -----------------------------------------------------------------------------
def neid(star: Dict[str, Any], folder: str, refresh: bool = False,
         name: Optional[str] = None
         ) -> Tuple[List[Tuple[RVData, Dict[str, Any]]], str]:
    """
    The velocities of NEID (WIYN, Kitt Peak) of a star: the table of the
    archive of NExScI, its public rows, not those its pipeline rejected or
    flagged (two of Barnard's star, 105 m/s from the others); its two
    resolutions apart. Its pipeline is told how the star moves: no
    perspective acceleration

    :param star: dict, ra, dec [deg, J2000], pmra, pmdec [mas/yr], and
                 have (the spectrographs files are given of, for eso())

    :return: tuple, the datasets and what was found, in words
    """
    rows = _tap(NEID_TAP, (
        "select object, ccfjdsum, ccfrvmod*1000.0 as rv, dvrms*1000.0 as "
        "err, obsmode, swversion, rejected, flagged, qrad, qdecd from neidl2 "
        "where "
        f"contains(point('icrs', qrad, qdecd), circle('icrs', "
        f"{float(star['ra']):.6f}, {float(star['dec']):.6f}, "
        f"{_reach(star):.5f}))=1"), 'NEID archive')
    rjd = np.array([_num(row['ccfjdsum']) for row in rows]) - 2400000.0
    keep = (_pointed([_num(row['qrad']) for row in rows],
                     [_num(row['qdecd']) for row in rows], rjd, star)
            & np.array([str(row.get('rejected', '0')).strip() in ('0', '')
                        and str(row.get('flagged', '0')).strip()
                        in ('0', '') for row in rows])
            ) if rows else np.zeros(0, bool)
    mode = np.array([str(row.get('obsmode') or '').strip().lower()
                     for row in rows])
    parts = []
    for one, label in (('hr', 'NEID (NExScI)'), ('he', 'NEID-HE (NExScI)')):
        sel = keep & (mode == one)
        if np.any(sel):
            got = _part(label, f'neid_{one}', 'NEID archive (NExScI)',
                        rjd[sel], [_num(row['rv']) for row, ok
                                   in zip(rows, sel) if ok],
                        [_num(row['err']) for row, ok in zip(rows, sel)
                         if ok], [row['object'] for row, ok
                                  in zip(rows, sel) if ok], name,
                        perspective=False)
            if got:
                parts.append(got)
    return parts, _told(parts, 'NEID has no public velocity of the star')


# -----------------------------------------------------------------------------
# The OHP: SOPHIE and ELODIE
# -----------------------------------------------------------------------------
def _ohp(base: str, table: str, star: Dict[str, Any], cols: str,
         where: str) -> List[Dict[str, str]]:
    """the rows of a table of an archive of the OHP in a box about the
    place of a star in 2000 (wide enough for where it is now)"""
    ra, dec = float(star['ra']) / 15.0, float(star['dec'])
    sign, dec = ('-' if dec < 0 else '+'), abs(dec)
    place = (f'J{int(ra):02d} {int(ra * 60) % 60:02d} '
             f'{ra * 3600 % 60:04.1f} {sign}{int(dec):02d} '
             f'{int(dec * 60) % 60:02d} {dec * 3600 % 60:02.0f}')
    size = max(2.0, 2.0 * 60.0 * _reach(star, 30.0))
    return _table(f'{base}?n={table}&c=o&p={urllib.parse.quote(place)}'
                  f'&f={size:.1f}&a=csv[,]&d={cols}', where)


def sophie(star: Dict[str, Any], folder: str, refresh: bool = False,
           name: Optional[str] = None
           ) -> Tuple[List[Tuple[RVData, Dict[str, Any]]], str]:
    """
    The velocities of SOPHIE (the 1.93 m of the OHP) of a star: the table
    of the CCFs of its archive, what is past its year of embargo (the
    rows of a programme still protected have no date: left out). SOPHIE
    and SOPHIE+ apart (the fibres of June 2011), and its two resolutions

    :return: tuple, the datasets and what was found, in words
    """
    rows = _ohp(SOPHIE_URL, 'sophiecc', star,
                'seq,objname,bjd,rv,err,mask,fwhm,contrast,span',
                'SOPHIE archive')
    modes = {row.get('seq'): str(row.get('mode') or '').upper()
             for row in _ohp(SOPHIE_URL, 'sophie', star, 'seq,mode',
                             'SOPHIE archive')} if rows else {}
    rows = [row for row in rows if _num(row.get('bjd')) > 2400000.0
            and 0 < _num(row.get('err')) < 0.09]
    parts = []
    if rows:
        mask = _most([row['mask'] for row in rows])
        other = sum(row['mask'] != mask for row in rows)
        # a spectrum cross-correlated twice with a mask: once
        seen, kept = set(), []
        for row in rows:
            if row['mask'] == mask and row['seq'] not in seen:
                seen.add(row['seq'])
                kept.append(row)
        rjd = np.array([_num(row['bjd']) for row in kept]) - 2400000.0
        mode = np.array([modes.get(row['seq'], '') or 'HR' for row in kept])
        for res in ('HR', 'HE'):
            for era, late in (('SOPHIE', False), ('SOPHIE+', True)):
                sel = (mode == res) & ((rjd >= SOPHIE_2011) == late)
                if not np.any(sel):
                    continue
                label = f'{era}{"" if res == "HR" else " HE"} (OHP)'
                got = _part(
                    label, re.sub(r'[^a-z+]+', '_', label.lower()).replace(
                        '+', 'p').strip('_'), 'SOPHIE archive (OHP)',
                    rjd[sel], [_num(row['rv']) * 1000.0 for row, ok
                               in zip(kept, sel) if ok],
                    [_num(row['err']) * 1000.0 for row, ok
                     in zip(kept, sel) if ok],
                    [row['objname'] for row, ok in zip(kept, sel) if ok],
                    name, said=(f'{other} with another mask than {mask}'
                                if other else ''),
                    **{key: [_num(row.get(col)) for row, ok
                             in zip(kept, sel) if ok]
                       for key, col in (('fwhm', 'fwhm'), ('bis', 'span'),
                                        ('contrast', 'contrast'))})
                if got:
                    parts.append(got)
    return parts, _told(parts, 'SOPHIE has no public velocity of the star')


def elodie(star: Dict[str, Any], folder: str, refresh: bool = False,
           name: Optional[str] = None
           ) -> Tuple[List[Tuple[RVData, Dict[str, Any]]], str]:
    """
    The velocities of ELODIE (the 1.93 m of the OHP, 1994 to 2006) of a
    star: the table of its archive, the exposures with the thorium lamp
    beside the star. Their errors are not in the table: the instrumental
    error and the photon noise of Perrier et al. (2003), see the module

    :return: tuple, the datasets and what was found, in words
    """
    rows = [row for row in _ohp(
        ELODIE_URL, 'e501', star,
        'objname,jdb,imatyp,masque,vfit,sn,sigfit', 'ELODIE archive')
        if _num(row.get('sn')) > 0]
    lamp = [row for row in rows if row.get('imatyp', '').startswith('OBTH')]
    parts = []
    if lamp:
        mask = _most([row['masque'] for row in lamp])
        kept = [row for row in lamp if row['masque'] == mask]
        said = [f'{len(rows) - len(lamp)} with no lamp beside the star'
                if len(rows) > len(lamp) else '',
                f'{len(lamp) - len(kept)} with another mask than {mask}'
                if len(lamp) > len(kept) else '']
        got = _part(
            'ELODIE (OHP)', 'elodie', 'ELODIE archive (OHP)',
            [_num(row['jdb']) for row in kept],
            [_num(row['vfit']) * 1000.0 for row in kept],
            [float(np.hypot(ELODIE_FLOOR, ELODIE_PHOTON / _num(row['sn'])))
             for row in kept], [row['objname'] for row in kept], name,
            said='; '.join(word for word in said if word),
            # the width of its CCF, from the sigma of the gaussian fitted
            fwhm=[2.0 * np.sqrt(2.0 * np.log(2.0)) * _num(row.get('sigfit'))
                  for row in kept])
        if got:
            parts.append(got)
    return parts, _told(parts, 'ELODIE has no velocity of the star with its '
                        'lamp')


# -----------------------------------------------------------------------------
# ESO: HARPS, ESPRESSO, NIRPS
# -----------------------------------------------------------------------------
def _eso_read(one: Dict[str, Any]) -> Dict[str, Any]:
    """what the header of a CCF file of ESO says"""
    try:
        head = tng.header(_http(ESO_FILE + one['id'], None, 'ESO archive',
                                120.0))
    except RuntimeError as err:
        return dict(error=str(err))
    return dict(
        error='', ins=one['ins'], target=one['target'],
        bjd=head.get('ESO QC BJD', ''), rv=head.get('ESO QC CCF RV', ''),
        err=head.get('ESO QC CCF RV ERROR', ''),
        mask=head.get('ESO QC CCF MASK', ''),
        mode=head.get('ESO INS MODE', ''),
        version=head.get('ESO PRO REC1 PIPE ID', ''),
        qc=head.get('ESO QC SCIRED CHECK', ''),
        fwhm=head.get('ESO QC CCF FWHM', ''),
        contrast=head.get('ESO QC CCF CONTRAST', ''),
        bis=head.get('ESO QC CCF BIS SPAN', ''))


def eso(star: Dict[str, Any], folder: str, refresh: bool = False,
        name: Optional[str] = None
        ) -> Tuple[List[Tuple[RVData, Dict[str, Any]]], str]:
    """
    The velocities of the pipelines of ESO of a star: HARPS from the
    catalogue of its velocities (2003 to 2023; before and after its
    fibres of 2015, its high-accuracy mode), ESPRESSO (before and after
    2019), NIRPS and HARPS by the pipeline of ESPRESSO from the headers
    of their CCF files, kept in a table of the folder (eso_exposures.csv)

    :return: tuple, the datasets and what was found, in words
    """
    parts = []
    ra, dec = float(star['ra']), float(star['dec'])
    rows = _tap(ESO_CAT, (
        "SELECT tel_object, ins_mode, drs_bjd, drs_ccf_rvc, drs_ccf_noise, "
        "drs_mask, ra_simbad, dec_simbad, drs_ccf_fwhm, drs_ccf_contrast, "
        "drs_bis_span FROM safcat.HARPS_RVCAT_V1 WHERE "
        "CONTAINS(POINT('ICRS', ra_simbad, dec_simbad), CIRCLE('ICRS', "
        f"{ra:.6f}, {dec:.6f}, 0.01))=1"), 'ESO archive', 'json')
    rows = [row for row in rows if str(row.get('ins_mode')).strip() == 'HAM'
            and tng.away(_num(row['ra_simbad']), _num(row['dec_simbad']),
                         star, tng.EPOCH) < LISTED]
    if rows:
        rjd = np.array([_num(row['drs_bjd']) for row in rows]) - 2400000.0
        for label, late in (('HARPS03 (ESO)', False),
                            ('HARPS15 (ESO)', True)):
            sel = (rjd >= HARPS_2015) == late
            if not np.any(sel):
                continue
            mine = [row for row, ok in zip(rows, sel) if ok]
            mask = _most([row['drs_mask'] for row in mine])
            same = np.array([str(row['drs_mask']) == mask for row in mine])
            got = _part(
                label, label.split()[0].lower() + '_eso',
                'ESO archive (HARPS, the catalogue of its velocities)',
                rjd[sel][same],
                [_num(row['drs_ccf_rvc']) * 1000.0 for row, ok
                 in zip(mine, same) if ok],
                [_num(row['drs_ccf_noise']) * 1000.0 for row, ok
                 in zip(mine, same) if ok],
                [str(row['tel_object']) for row, ok in zip(mine, same)
                 if ok], name,
                said=(f'{int(np.sum(~same))} with another mask than {mask}'
                      if np.any(~same) else ''),
                **{key: [_num(row.get(col)) for row, ok in zip(mine, same)
                         if ok]
                   for key, col in (('fwhm', 'drs_ccf_fwhm'),
                                    ('contrast', 'drs_ccf_contrast'),
                                    ('bis', 'drs_bis_span'))})
            if got:
                parts.append(got)
    # the CCF files of the other pipelines
    now = _time.strftime('%Y-%m-%dT%H:%M:%S', _time.gmtime())
    found = _tap(ESO_OBS, (
        "SELECT o.instrument_name AS ins, o.target_name AS target, o.t_min "
        "AS t_min, o.obs_release_date AS release, o.s_ra AS ra, o.s_dec AS "
        "dec, f.archive_id AS id, f.original_filename AS orig FROM "
        "ivoa.ObsCore AS o JOIN phase3v2.product_files AS f ON "
        "f.product_id = o.dp_id WHERE CONTAINS(POINT('ICRS', o.s_ra, "
        f"o.s_dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {_reach(star):.5f}"
        "))=1 AND o.instrument_name IN ('HARPS', 'ESPRESSO', 'NIRPS') AND "
        "f.eso_category = 'ANCILLARY.CCF'"), 'ESO archive', 'json')
    listed = []
    for row in found:
        orig = str(row.get('orig') or '')
        # the star's fibre; a NIRPS exposure has its CCF twice, with and
        #   without the tellurics taken out: the first
        if re.search(r'(_B\.fits|SCCB|SKYSUB)', orig) or (
                row['ins'] == 'NIRPS' and 'TELL_CORR' not in orig) or (
                row['ins'] != 'NIRPS' and 'TELL_CORR' in orig):
            continue
        if not _pointed(_num(row['ra']), _num(row['dec']),
                        _num(row['t_min']), star)[0]:
            continue
        # NIRPS defers to a file of NIRPS: not fetched when there is one
        if row['ins'] == 'NIRPS' and 'NIRPS' in (star.get('have') or ()):
            continue
        listed.append(dict(file=str(row['id']), id=str(row['id']),
                           ins=str(row['ins']), target=str(row['target']),
                           public=str(row.get('release') or '9') <= now))
    private = sum(not one['public'] for one in listed)
    heads = _headers(listed, os.path.join(folder, 'eso_exposures.csv'),
                     _eso_read, 'ESO archive', refresh, 'eso')
    heads = [row for row in heads if str(row.get('qc', '1')) != '0']
    for ins in ('ESPRESSO', 'NIRPS', 'HARPS'):
        mine = [row for row in heads if row['ins'] == ins
                and _num(row['rv']) == _num(row['rv'])]
        if not mine:
            continue
        rjd = np.array([_num(row['bjd']) for row in mine]) - 2400000.0
        mode = np.array([str(row['mode']) for row in mine])
        main = _most(mode)
        # the files are of the pipeline of their day: a new one (espdr 3
        #   after espdr 2) is another zero point, a dataset of its own
        major = np.array([(re.findall(r'/(\d+)', str(row['version']))
                           or [''])[0] for row in mine])
        first = min(major) if ins == 'ESPRESSO' else ''
        for kind, era, late, ver in [
                (kind, era, late, ver) for kind in np.unique(mode)
                for era, late in ((('18', False), ('19', True))
                                  if ins == 'ESPRESSO' else (('', False),))
                for ver in (np.unique(major) if ins == 'ESPRESSO'
                            else ('',))]:
            sel = mode == kind
            if ins == 'ESPRESSO':
                sel &= ((rjd >= ESPRESSO_2019) == late) & (major == ver)
            if not np.any(sel):
                continue
            these = [row for row, ok in zip(mine, sel) if ok]
            mask = _most([row['mask'] for row in these])
            same = np.array([str(row['mask']) == mask
                             for row in these])
            tail = ('' if kind == main or (ins == 'NIRPS'
                                           and kind == 'HE')
                    else f'_{kind}')
            label = (f'{"HARPS-espdr" if ins == "HARPS" else ins}{era}'
                     f'{tail} (ESO{"" if ver == first else " v" + ver})')
            got = _part(
                label, re.sub(r'[^a-z0-9]+', '_', label.lower()
                              ).strip('_'),
                f'ESO archive ({ins}, its CCF files)', rjd[sel][same],
                [_num(row['rv']) * 1000.0 for row, ok
                 in zip(these, same) if ok],
                [_num(row['err']) * 1000.0 for row, ok
                 in zip(these, same) if ok],
                [str(row['target']) for row, ok in zip(these, same)
                 if ok], name, defer=ins == 'NIRPS',
                said=(f'{int(np.sum(~same))} with another mask than '
                      f'{mask}' if np.any(~same) else ''),
                **{key: [_num(row.get(key)) for row, ok
                         in zip(these, same) if ok]
                   for key in ('fwhm', 'contrast', 'bis')})
            if got:
                parts.append(got)
    words = _told(parts, 'ESO has no public velocity of the star')
    return parts, words + (f'; {private} files not public yet' if private
                           else '')


# -----------------------------------------------------------------------------
# The CFHT: SPIRou
# -----------------------------------------------------------------------------
def _cadc_read(one: Dict[str, Any]) -> Dict[str, Any]:
    """what the headers of a CCF file of APERO say (asked alone, as
    text: the header of each extension, one card a line)"""
    try:
        text = _http(CADC_FILE + one['uri'] + '?META=true', None,
                     'CFHT archive (CADC)', 120.0).decode('ascii', 'replace')
    except RuntimeError as err:
        return dict(error=str(err))
    head: Dict[str, str] = {}
    for line in text.splitlines():
        if line[8:10] == '= ':
            head[line[:8].strip()] = line[10:].split(' /')[0].strip(" '")
    return dict(error='', target=one['target'], bjd=head.get('BJD', ''),
                rv=head.get('RV_CORR') or head.get('RV_OBJ', ''),
                err=head.get('DVRMS_CC', ''), mask=head.get('CCFMASK', ''),
                version=head.get('VERSION', ''), qc=head.get('QCC_ALL', ''),
                dpr=head.get('DPRTYPE', ''), fwhm=head.get('CCFMFWHM', ''),
                contrast=head.get('CCFMCONT', ''),
                bis=head.get('CCFMBISV', ''))


def cfht(star: Dict[str, Any], folder: str, refresh: bool = False,
         name: Optional[str] = None
         ) -> Tuple[List[Tuple[RVData, Dict[str, Any]]], str]:
    """
    The velocities of SPIRou (CFHT) of a star: those of the CCF of APERO,
    from the headers of its public files at the CADC, kept in a table of
    the folder (cfht_exposures.csv). APERO is told how the star moves: no
    perspective acceleration. A velocity of a CCF in the near infrared
    is a poor one: the dataset defers to a file of SPIRou (see the
    module)

    :return: tuple, the datasets and what was found, in words
    """
    rows = _tap(CADC_TAP, (
        "SELECT o.target_name AS target, o.targetPosition_coordinates_cval1"
        " AS ra, o.targetPosition_coordinates_cval2 AS dec, "
        "p.time_bounds_lower AS mjd, p.dataRelease AS release, a.uri AS uri"
        " FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID = "
        "p.obsID JOIN caom2.Artifact AS a ON a.planeID = p.planeID WHERE "
        "o.collection = 'CFHT' AND o.instrument_name = 'SPIRou' AND "
        "p.productID LIKE '%v' AND a.contentType = 'application/fits' AND "
        f"INTERSECTS(CIRCLE('ICRS', {float(star['ra']):.6f}, "
        f"{float(star['dec']):.6f}, {_reach(star):.5f}), "
        "p.position_bounds) = 1"), 'CFHT archive (CADC)')
    now = _time.strftime('%Y-%m-%dT%H:%M:%S', _time.gmtime())
    listed = [dict(file=str(row['uri']).rsplit('/', 1)[-1],
                   uri=str(row['uri']), target=str(row['target']),
                   public=str(row.get('release') or '9') <= now)
              for row in rows
              if _pointed(_num(row['ra']), _num(row['dec']),
                          _num(row['mjd']), star)[0]]
    private = sum(not one['public'] for one in listed)
    heads = [row for row in _headers(
        listed, os.path.join(folder, 'cfht_exposures.csv'), _cadc_read,
        'SPIRou, CFHT archive', refresh, 'cfht')
        if str(row.get('qc', 'T')).strip().upper() not in ('F', 'FALSE')
        and str(row.get('dpr', 'OBJ')).upper().startswith(('OBJ', 'POLAR'))]
    parts = []
    if heads:
        mask = _most([row['mask'] for row in heads])
        kept = [row for row in heads if row['mask'] == mask]
        got = _part(
            'SPIRou (CFHT)', 'spirou_cfht', 'CFHT archive (SPIRou, the CCF '
            'of APERO, CADC)',
            [_num(row['bjd']) - 2400000.0 for row in kept],
            [_num(row['rv']) * 1000.0 for row in kept],
            [_num(row['err']) for row in kept],
            [row['target'] for row in kept], name, perspective=False,
            defer=True, said=(f'{len(heads) - len(kept)} with another mask '
                              f'than {mask}' if len(heads) > len(kept)
                              else ''),
            fwhm=[_num(row.get('fwhm')) for row in kept],
            # APERO gives the contrast as a fraction
            contrast=[100.0 * _num(row.get('contrast')) for row in kept],
            bis=[_num(row.get('bis')) for row in kept])
        if got:
            parts.append(got)
    words = _told(parts, 'SPIRou has no public velocity of the star')
    return parts, words + (f'; {private} files not public yet' if private
                           else '')


#: the archives, by their key in koloa.published.ARCHIVES
SOURCES: Dict[str, Callable] = dict(eso=eso, neid=neid, sophie=sophie,
                                    elodie=elodie, cfht=cfht)


def neid_pointings(refresh: bool = False) -> Dict[str, np.ndarray]:
    """
    Every place NEID has velocities of (a thousand; one question, 10 s,
    kept MAX_AGE days), so that whether it has a star is known at once
    for a whole sample

    :return: dict, ra and dec [deg], and n (their exposures)
    """
    kept = None
    if os.path.exists(NEID_POINTINGS) and not refresh:
        with open(NEID_POINTINGS) as handle:
            kept = json.load(handle)
        if _time.time() - kept.get('asked', 0) > MAX_AGE * 86400:
            kept = None
    if kept is None:
        rows = _tap(NEID_TAP, 'select qrad, qdecd, count(*) as n from '
                    'neidl2 group by qrad, qdecd', 'NEID archive')
        kept = dict(asked=_time.time(),
                    ra=[_num(row['qrad']) for row in rows],
                    dec=[_num(row['qdecd']) for row in rows],
                    n=[int(_num(row['n'])) for row in rows])
        os.makedirs(os.path.dirname(NEID_POINTINGS), exist_ok=True)
        with open(NEID_POINTINGS, 'w') as handle:
            json.dump(kept, handle)
        log(f'NEID archive: {len(rows)} pointings', 'value')
    return dict(ra=np.array(kept['ra'], float),
                dec=np.array(kept['dec'], float),
                n=np.array(kept['n'], int))


def has(key: str, star: Dict[str, Any]) -> bool:
    """
    Whether an archive has public velocities of a star, asked lightly
    (for the stars of a sample, before anything is gathered): HARPS-N and
    NEID from the lists of where they were pointed, SOPHIE and ELODIE
    from their tables (a second each). The archives asked only when
    named (ESO, the CFHT) are not looked at: False

    :param key: str, the archive (koloa.published.ARCHIVES)
    :param star: dict, ra, dec [deg, J2000], pmra, pmdec [mas/yr]

    :return: bool
    """
    if key == 'tng':
        return bool(tng.near(star)['free'])
    if key == 'neid':
        there = neid_pointings()
        if not len(there['ra']):
            return False
        return bool(np.any(
            (tng.away(there['ra'], there['dec'], star) < POINTED)
            | (tng.offpath(there['ra'], there['dec'], star, start=2015.0)
               < POINTED)))
    if key in ('sophie', 'elodie'):
        return bool(SOURCES[key](star, '')[0])
    return False

# =============================================================================
# End of code
# =============================================================================

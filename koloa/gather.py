#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Everything public about a star, from its SIMBAD name: the velocities of
DACE and of CARMENES, and the TESS light curves, kept in one folder.

    from koloa.gather import gather, load
    manifest = gather('GJ 436', 'archives')   # archives/GJ_436/...
    star = load('archives/GJ_436')            # star['rv'], an RVData

The folder of a star (the name as given, spaces as underscores):

    target.json          SIMBAD: main name, every identifier, position, TIC,
                         and the periods of variability it lists (ROT for a
                         rotation), with the rotation period of CARMENES
    manifest.json        what each archive gave: files, instruments, number
                         of points, time span, or why nothing
    rv/all_rv.csv        every velocity, one instrument per instrument era,
                         with the indicators of each (nan for the others):
                         RVData.from_csv(path, inst='inst')
    rv/dace/raw_<name>.csv   DACE's answer as it came (koloa.dace.fetch)
    rv/dace/<INST>.csv   one instrument era of DACE (HARPS03, HARPS15,
                         ESPRESSO19, NIRPS...): the velocities that pass the
                         pipeline's quality control, with the indicators
    rv/carmenes/raw.csv  the CARMENES DR1 rows of the star (GAVO's TAP
                         service, Ribas et al. 2023), every scalar column
    rv/carmenes/CARMENES.csv   the velocities corrected for the nightly zero
                         points (SERVAL's a_rv), with the indicators
    rv/published/        the velocities published on VizieR (koloa.published:
                         the surveys, Keck HIRES, the APF, the Lick
                         Hamilton, HARPS by SERVAL; the tables of the star's
                         papers), a CSV per source and published.json
    phot/tess/fits/      the TESS light curves of MAST, one per sector (the
                         best of SPOC 2-minute, TESS-SPOC, QLP)
    phot/tess/s<sector>_<pipeline>.csv   each sector: rjd, flux, sflux [ppt]
    phot/tess.csv        every sector together, with its sector and pipeline
    phot/kepler.csv, phot/k2.csv, phot/corot.csv   the same from Kepler, K2
                         and CoRoT when the star lay in their fields (and
                         phot/<mission>/<stretch>.csv, each quarter,
                         campaign or run)

Times are BJD - 2400000 (rjd), velocities m/s. A file that is there is read
back rather than asked for again, unless refresh=True.

DACE: with no API key, the public data only. A key is looked for in the
environment variable DACE_API_KEY, then in ~/.dacerc (koloa.dace.find_key);
with one, DACE adds what its account may see, data that may not be public,
and the manifest says which key was used (where it came from, never the
key). Some networks filter DACE: it is then reported as unreachable and the
rest goes on. CARMENES DR1 (the
GTO of 2016 to 2020) covers about 360 M dwarfs of the north; its velocities
(a_rv) are per exposure, corrected for the nightly zero points, and leave
out the exposures of nights the survey did not keep (transit sequences,
mostly). GAVO's metadata give them in km/s, but they are m/s (errors of
about 1 m/s).

Created on 2026-10-01

@author: artigau
"""
import csv
import json
import os
import re
import time as _time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from koloa.data import RVData, merge
from koloa.log import log, outcome, step, row
from koloa.paths import cache

# =============================================================================
# Define variables
# =============================================================================
#: the TAP service that serves CARMENES DR1
GAVO_TAP = 'https://dc.g-vo.org/tap/sync'
#: SIMBAD's TAP service
SIMBAD_TAP = 'https://simbad.cds.unistra.fr/simbad/sim-tap/sync'
#: where the list of the stars of CARMENES DR1 is kept, once fetched
CARMENES_CACHE = cache('carmenes_objects.json')
#: how far from SIMBAD's position a CARMENES star may be [degrees]
CARMENES_RADIUS = 0.005
#: the scalar columns of carmenes.rvs (its per-order columns are left out)
CARMENES_COLUMNS = [
    'bjd', 'rv', 'err_rv', 'a_rv', 'err_a_rv', 'dlw', 'err_dlw', 'ccf_rv',
    'err_ccf_rv', 'ccf_fwhm', 'err_ccf_fwhm', 'ccf_contrast',
    'err_ccf_contrast', 'ccf_bis', 'err_ccf_bis', 'halpha', 'err_halpha',
    'cai', 'err_cai', 'cai_8498', 'err_cai_8498', 'cai_8542', 'err_cai_8542',
    'cai_8698', 'err_cai_8698', 'nad1', 'err_nad1', 'nad2', 'err_nad2',
    'spec_flag', 'drift', 'err_drift', 'berv', 'sadrift', 'rvmean',
    'err_rvmean', 'rvmed', 'err_rvmed', 'spec_accref']
#: the CARMENES indicators kept with the velocities
CARMENES_INDICATORS = ['dlw', 'ccf_fwhm', 'ccf_contrast', 'ccf_bis',
                       'halpha', 'cai', 'cai_8498', 'cai_8542', 'cai_8698',
                       'nad1', 'nad2']


# =============================================================================
# Define functions
# =============================================================================
def folder_name(target: str) -> str:
    """the folder of a star: its name as given, spaces as underscores"""
    return re.sub(r'[^A-Za-z0-9+\-.]+', '_', target.strip()).strip('_')


def _tap(url: str, query: str, timeout: float = 120.0) -> List[Dict[str, str]]:
    """the rows of an ADQL query to a TAP service, as text"""
    body = urllib.parse.urlencode(dict(REQUEST='doQuery', LANG='ADQL',
                                       FORMAT='csv', QUERY=query)).encode()
    request = urllib.request.Request(url, data=body,
                                     headers={'User-Agent': 'koloa'})
    with urllib.request.urlopen(request, timeout=timeout) as resp:
        text = resp.read().decode('utf-8', 'replace')
    return list(csv.DictReader(text.splitlines()))


def _float(value: Any) -> float:
    """a number, nan when there is none"""
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _cell(value: Any) -> str:
    """a number written in full, an empty cell for nan"""
    if isinstance(value, (float, np.floating)):
        return repr(float(value)) if np.isfinite(value) else ''
    return str(value)


def write_rv(data: RVData, path: str, source: Optional[str] = None) -> str:
    """
    A series as a CSV file that RVData.from_csv(path, inst='inst') reads
    back: rjd, vrad (as given, the zero points put back), svrad, inst, the
    source, and each indicator with its error (s<name>)

    :return: str, the file
    """
    vrad = data.rv + np.array([data.zero_point[str(name)]
                               for name in data.inst])
    header = ['rjd', 'vrad', 'svrad', 'inst'] + (['source'] if source else [])
    for key in data.indicators:
        header += [key, 's' + key]
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for ip in range(data.n):
            row = [_cell(data.time[ip]), _cell(vrad[ip]), _cell(data.err[ip]),
                   str(data.inst[ip])] + ([source] if source else [])
            for val, err in data.indicators.values():
                row += [_cell(val[ip]), _cell(err[ip])]
            writer.writerow(row)
    return path


def _summary(data: RVData, files: List[str], root: str) -> Dict[str, Any]:
    """what a series holds, for the manifest"""
    return dict(status='ok', files=[os.path.relpath(fl, root) for fl in files],
                instruments={name: int(np.sum(data.inst == name))
                             for name in data.instruments},
                npoints=data.n, rjd=[float(data.time.min()),
                                     float(data.time.max())],
                indicators=sorted(data.indicators))


# -----------------------------------------------------------------------------
# SIMBAD
# -----------------------------------------------------------------------------
def variability(main: str, timeout: float = 10.0) -> List[Dict[str, Any]]:
    """
    The periods of variability SIMBAD lists for a star (its mesVar table):
    type (ROT for a rotation), period [days] and reference

    SIMBAD's TAP service is at times very slow (a minute and more for the
    simplest question): beyond the timeout, an error, which a caller takes
    as no answer for now.

    :param main: str, SIMBAD's main name of the star
    :param timeout: float, the longest wait for an answer [s]
    """
    name = main.replace("'", "''")
    rows = _tap(SIMBAD_TAP, (
        'SELECT v.vartyp, v.period, v.bibcode FROM mesVar AS v JOIN basic '
        f"AS b ON v.oidref = b.oid WHERE b.main_id = '{name}'"),
                timeout=timeout)
    return [dict(type=row['vartyp'], period=_float(row['period']),
                 bibcode=row['bibcode'])
            for row in rows if np.isfinite(_float(row['period']))]


# -----------------------------------------------------------------------------
# The archives
# -----------------------------------------------------------------------------
def dace_rv(ident: Dict[str, Any], target: str, folder: str,
            api_key: Any = None, refresh: bool = False
            ) -> Optional[RVData]:
    """
    The velocities of a star on DACE, trying its names in turn (HD69830,
    GJ687...): the raw answer and one file per instrument era

    :param ident: dict, from koloa.archive.resolve
    :param target: str, the name as given
    :param folder: str, where the files go
    :param api_key: str, None or False: a DACE API key, None to look for
                    one, False for none (koloa.dace.find_key)
    :param refresh: bool, ask DACE even when a file is there

    :return: RVData or None (DACE does not know the star, or unreachable:
             RuntimeError then)
    """
    from koloa.dace import fetch, names, rvdata
    for name in names(ident, target):
        # a name as a file's: nothing a disk or a shell could take for
        #   something else (a * of SIMBAD's, say)
        path = os.path.join(folder, f'raw_{folder_name(name)}.csv')
        try:
            fetch(name, path, api_key=api_key, refresh=refresh)
        except RuntimeError as err:
            if 'could not be reached' in str(err):
                raise
            log(f'DACE, {name}: {str(err).split(";")[0]}', 'info')
            continue
        data = rvdata(path, name=target)
        for inst in data.instruments:
            write_rv(data.select(data.inst == inst),
                     os.path.join(folder, f'{inst}.csv'))
        return data
    return None


def carmenes_objects(refresh: bool = False) -> List[Dict[str, Any]]:
    """
    The stars of CARMENES DR1 (361, GAVO's carmenes.objects), fetched once
    and kept in CARMENES_CACHE: a star is then found without the network

    :param refresh: bool, fetch the list again

    :return: list of dict, carmenes_id, name, nobs, nspect, mass, radius,
             teff, p_rot, p_rot_source, survey, remarks, ra, dec
    """
    if os.path.exists(CARMENES_CACHE) and not refresh:
        with open(CARMENES_CACHE) as handle:
            return json.load(handle)
    log('fetching the list of the stars of CARMENES DR1 (once; kept in '
        f'{CARMENES_CACHE})')
    rows = _tap(GAVO_TAP, (
        'SELECT carmenes_id, name, nobs, nspect, mass, radius, teff, p_rot, '
        'p_rot_source, survey, remarks, ra, dec FROM carmenes.objects'))
    os.makedirs(os.path.dirname(CARMENES_CACHE), exist_ok=True)
    with open(CARMENES_CACHE, 'w') as handle:
        json.dump(rows, handle)
    return rows


def carmenes_star(ra: float, dec: float) -> Optional[Dict[str, Any]]:
    """
    The CARMENES DR1 star at a position (its Karmn name, number of
    velocities, rotation period...), None if DR1 has none there (the list
    of its stars kept here, carmenes_objects)

    :param ra: float, J2000 [degrees]
    :param dec: float, J2000 [degrees]
    """
    best, near = None, CARMENES_RADIUS
    cosd = np.cos(np.radians(dec))
    for row in carmenes_objects():
        dra = (_float(row['ra']) - ra + 180.0) % 360.0 - 180.0
        dist = float(np.hypot(dra * cosd, _float(row['dec']) - dec))
        if dist <= near:
            best, near = dict(row), dist
    return best


def carmenes_rv(star: Dict[str, Any], target: str, folder: str,
                refresh: bool = False) -> Optional[RVData]:
    """
    The CARMENES DR1 velocities of a star (see the module): the rows as
    served, and the velocities corrected for the nightly zero points with
    their indicators

    :param star: dict, from carmenes_star
    :param target: str, the name as given
    :param folder: str, where the files go
    :param refresh: bool, ask GAVO even when the file is there

    :return: RVData or None (no corrected velocity)
    """
    raw = os.path.join(folder, 'raw.csv')
    if refresh or not os.path.exists(raw):
        karmn = star['carmenes_id'].replace("'", "''")
        rows = _tap(GAVO_TAP, f'SELECT {", ".join(CARMENES_COLUMNS)} FROM '
                              f"carmenes.rvs WHERE object = '{karmn}'")
        os.makedirs(folder, exist_ok=True)
        with open(raw, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=CARMENES_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
    with open(raw) as handle:
        rows = list(csv.DictReader(handle))
    col = {key: np.array([_float(row.get(key)) for row in rows])
           for key in CARMENES_COLUMNS if key != 'spec_accref'}
    good = np.isfinite(col['a_rv']) & np.isfinite(col['err_a_rv'])
    if not np.any(good):
        return None
    indicators = {key: (col[key][good], col['err_' + key][good])
                  for key in CARMENES_INDICATORS
                  if np.sum(np.isfinite(col[key][good])
                            & np.isfinite(col['err_' + key][good])) >= 5}
    data = RVData(col['bjd'][good] - 2400000.0, col['a_rv'][good],
                  col['err_a_rv'][good],
                  inst=np.full(int(np.sum(good)), 'CARMENES'),
                  indicators=indicators, name=target)
    write_rv(data, os.path.join(folder, 'CARMENES.csv'))
    return data


def tess_photometry(ident: Dict[str, Any], target: str, folder: str,
                    refresh: bool = False) -> Dict[str, Any]:
    """
    The TESS light curves of a star (koloa.tess.light_curves, every sector):
    the FITS files, one CSV per sector and one with every sector

    :return: dict, from koloa.tess.light_curves
    """
    from koloa.tess import light_curves
    name = ident.get('tic') or target
    lcs = light_curves(name, folder=os.path.join(folder, 'tess', 'fits'),
                       max_sectors=1000, refresh=refresh)
    allrows = []
    for lc in lcs['sectors']:
        path = os.path.join(folder, 'tess', f's{lc["sector"]:04d}_'
                            f'{lc["provenance"]}.csv')
        with open(path, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['rjd', 'flux', 'sflux'])
            for row in zip(lc['time'], lc['flux'], lc['err']):
                writer.writerow([_cell(val) for val in row])
        lc['csv'] = path
        allrows += [[_cell(tt), _cell(ff), _cell(ee), str(lc['sector']),
                     lc['provenance']]
                    for tt, ff, ee in zip(lc['time'], lc['flux'], lc['err'])]
    if allrows:
        with open(os.path.join(folder, 'tess.csv'), 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['rjd', 'flux', 'sflux', 'sector', 'pipeline'])
            writer.writerows(allrows)
    return lcs


def space_photometry(ident: Dict[str, Any], target: str, folder: str,
                     refresh: bool = False) -> Dict[str, Dict[str, Any]]:
    """
    The light curves of a star from Kepler, K2 and CoRoT, when it lay in
    their fields (koloa.photometry): one CSV per quarter, campaign or run,
    and one per mission with all of them (phot/kepler.csv, phot/k2.csv,
    phot/corot.csv: rjd, flux, sflux [ppt], sector, pipeline; the sector is
    the number of the stretch, koloa.photometry.label names it)

    :return: dict, kepler, k2, corot (koloa.photometry.light_curves)
    """
    from koloa.photometry import light_curves
    lcs = light_curves(ident, refresh=refresh)
    for key, val in lcs.items():
        allrows = []
        for lc in val['sectors']:
            path = os.path.join(folder, key,
                                lc['label'].replace(' ', '_') + '.csv')
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, 'w', newline='') as handle:
                writer = csv.writer(handle)
                writer.writerow(['rjd', 'flux', 'sflux'])
                for row in zip(lc['time'], lc['flux'], lc['err']):
                    writer.writerow([_cell(val_) for val_ in row])
            lc['csv'] = path
            allrows += [[_cell(tt), _cell(ff), _cell(ee), str(lc['sector']),
                         lc['provenance']]
                        for tt, ff, ee in zip(lc['time'], lc['flux'],
                                              lc['err'])]
        if allrows:
            with open(os.path.join(folder, f'{key}.csv'), 'w',
                      newline='') as handle:
                writer = csv.writer(handle)
                writer.writerow(['rjd', 'flux', 'sflux', 'sector',
                                 'pipeline'])
                writer.writerows(allrows)
    return lcs


# -----------------------------------------------------------------------------
# All of them
# -----------------------------------------------------------------------------
def gather(target: str, root: str = '.', dace: bool = True,
           carmenes: bool = True, tess: bool = True,
           api_key: Any = None, refresh: bool = False,
           vizier: bool = True,
           vizier_sources: Optional[Sequence[str]] = None,
           have: Sequence[str] = ()) -> Dict[str, Any]:
    """
    Everything public about a star, in root/<target> (see the module)

    The archives are asked all at once; one that fails is reported in the
    manifest and the others go on.

    :param target: str, a name SIMBAD knows (GJ 436, Ross 905, TOI-700...)
    :param root: str, the folder of the folders of the stars
    :param dace: bool, the velocities of DACE
    :param carmenes: bool, the velocities of CARMENES DR1
    :param tess: bool, the light curves of TESS, and of Kepler, K2 and
                 CoRoT when the star lay in their fields
    :param vizier: bool, the velocities published on VizieR (the surveys:
                   Keck HIRES, the APF, the Lick Hamilton, HARPS by SERVAL;
                   the tables of the star's papers), koloa.published
    :param vizier_sources: list of str or None, the sources of VizieR
                           asked (koloa.published.SOURCES: teklu25, cls21,
                           talor19, fischer14, rvbank20, papers), None for
                           all
    :param api_key: str, None or False: a DACE API key, None to look for
                    one (DACE_API_KEY, ~/.dacerc), False for the public
                    data only (the key is never written anywhere)
    :param refresh: bool, ask again for what is already there
    :param have: the spectrographs you have files of: SPIRou at the CFHT
                 and NIRPS at ESO, which defer to them, are not asked

    :return: dict, the manifest (also in manifest.json): folder, target,
             and per archive its status (ok, none, unreachable, error),
             files, instruments, number of points, time span or message
    """
    from koloa.archive import resolve
    folder = os.path.join(root, folder_name(target))
    rvdir, photdir = os.path.join(folder, 'rv'), os.path.join(folder, 'phot')
    os.makedirs(folder, exist_ok=True)
    step('SIMBAD')
    ident = resolve(target)
    log(f'{target}: SIMBAD {ident["main"]}, TIC '
        f'{(ident.get("tic") or "none").replace("TIC ", "")}', 'info')
    info = dict(ident)
    # the spectral type of the star and a rough mass from it
    from koloa.stars import stellar
    info['star'] = stellar(ident)
    outcome(f'{ident["main"]}'
            + (f', {ident["tic"]}' if ident.get('tic') else '')
            + (f', {ident["sptype"]}' if ident.get('sptype') else '')
            + (f', {info["star"]["mass"]:.2f} Msun (rough)'
               if info['star'].get('mass') else ''))
    # how the star moves (Gaia DR3): an archive of an observatory is
    #   searched by where it was, and the velocities of a pipeline hold
    #   its perspective acceleration, taken out when the series is read
    from koloa.secular import star_astrometry, star_perspective
    try:
        info['astrometry'] = star_astrometry(ident)
        info['perspective'] = star_perspective(info)
    except Exception as err:  # a help too
        info['astrometry'] = info['perspective'] = None
        log(f'astrometry: {type(err).__name__}: {err}', 'warn')
    if info.get('perspective'):
        log(f'perspective acceleration: {info["perspective"]["value"]:.3f} '
            f'm/s/yr ({info["perspective"]["source"]})', 'value')
    ident = dict(ident, astrometry=info['astrometry'])
    manifest: Dict[str, Any] = dict(
        target=target, simbad=ident['main'], folder=os.path.abspath(folder),
        created=_time.strftime('%Y-%m-%d %H:%M:%S'), archives={})

    def attempt(key, func):
        """one archive, its failure reported and not raised"""
        try:
            return func()
        except Exception as err:
            text = str(err).split(';')[0]
            status = 'unreachable' if 'reached' in text else 'error'
            manifest['archives'][key] = dict(status=status, message=text)
            log(f'{key}: {text}', 'warn')
            return None

    def ask_variability():
        """the variability of the star (SIMBAD): a help, not a need"""
        try:
            info['variability'] = variability(ident['main'])
        except Exception as err:
            info['variability'] = []
            log(f'SIMBAD variability: {err}', 'warn')

    def ask_dace():
        """the velocities of DACE"""
        row('dace', 'asking')
        ddir = os.path.join(rvdir, 'dace')
        data = attempt('dace', lambda: dace_rv(ident, target, ddir,
                                               api_key=api_key,
                                               refresh=refresh))
        if data is not None:
            files = [os.path.join(ddir, fl) for fl in sorted(os.listdir(ddir))]
            from koloa.dace import find_key
            origin = find_key(api_key)[1]
            manifest['archives']['dace'] = dict(
                _summary(data, files, folder),
                key=(f'the key of {origin}: may hold data that are not '
                     f'public' if origin else 'none: public data only'))
        elif 'dace' not in manifest['archives']:
            manifest['archives']['dace'] = dict(
                status='none', message='DACE has no public velocities under '
                'the names of the star')
        log(f'DACE: {_told(manifest["archives"]["dace"])}', 'value')
        row('dace', *state_of(manifest['archives']['dace']))
        return data

    def ask_carmenes():
        """the velocities of CARMENES DR1"""
        row('carmenes', 'asking')
        cdir = os.path.join(rvdir, 'carmenes')
        star = (attempt('carmenes', lambda: carmenes_star(ident['ra'],
                                                          ident['dec']))
                if ident.get('ra') is not None else None)
        data = None
        if star is not None:
            info['carmenes'] = star
            data = attempt('carmenes', lambda: carmenes_rv(star, target, cdir,
                                                           refresh=refresh))
            if data is not None:
                manifest['archives']['carmenes'] = dict(
                    _summary(data, [os.path.join(cdir, 'raw.csv'),
                                    os.path.join(cdir, 'CARMENES.csv')],
                             folder), karmn=star['carmenes_id'])
        if 'carmenes' not in manifest['archives']:
            manifest['archives']['carmenes'] = dict(
                status='none', message='not in CARMENES DR1, or no '
                'velocity corrected for the nightly zero points')
        log(f'CARMENES DR1: {_told(manifest["archives"]["carmenes"])}',
            'value')
        row('carmenes', *state_of(manifest['archives']['carmenes']))
        return data

    def ask_published():
        """the velocities of the surveys and papers on VizieR, and of the
        archives of the observatories (koloa.published)"""
        from koloa.published import fetch as published_fetch
        pdir = os.path.join(rvdir, 'published')
        notes = attempt('published', lambda: published_fetch(
            ident, pdir, refresh=refresh, sources=vizier_sources,
            have=have))
        if notes is not None:
            manifest['archives']['published'] = _published(notes)
            kept = [note for note in notes if note.get('file')]
            nsurv = sum(note['kind'] == 'survey' for note in kept)
            narch = sum(note['kind'] == 'archive' for note in kept)
            npap = len(kept) - nsurv - narch
            log('VizieR and the archives of the observatories: '
                + _told(manifest['archives']['published'])
                + (f'; from {nsurv} survey{"s" * (nsurv != 1)}, '
                   f'{npap} paper{"s" * (npap != 1)} and {narch} '
                   f'reduction{"s" * (narch != 1)} of the archive of '
                   f'an observatory' if kept else ''), 'value')

    def ask_tess():
        """the light curves of TESS"""
        row('tess', 'asking')
        lcs = attempt('tess', lambda: tess_photometry(ident, target, photdir,
                                                      refresh=refresh))
        if lcs is not None:
            sectors = lcs['sectors']
            manifest['archives']['tess'] = (dict(
                status='ok', tic=lcs['tic'],
                sectors=[dict(sector=lc['sector'], pipeline=lc['provenance'],
                              exposure=lc['exposure'], npoints=len(lc['time']),
                              file=os.path.relpath(lc['csv'], folder))
                         for lc in sectors],
                files=['phot/tess.csv']) if sectors else dict(
                status='none', tic=lcs['tic'],
                message='no light curve at MAST'))
        if 'tess' in manifest['archives']:
            log(f'TESS: {_told(manifest["archives"]["tess"])}', 'value')
            row('tess', *state_of(manifest['archives']['tess']))

    def ask_space():
        """the light curves of the other space missions, where the star
        lay in their fields"""
        row('space', 'asking')
        space = attempt('space', lambda: space_photometry(
            ident, target, photdir, refresh=refresh))
        told = []
        for key, val in (space or {}).items():
            sectors = val['sectors']
            manifest['archives'][key] = (dict(
                status='ok', target=val['target'],
                sectors=[dict(sector=lc['sector'], label=lc['label'],
                              pipeline=lc['provenance'],
                              exposure=lc['exposure'], npoints=len(lc['time']),
                              file=os.path.relpath(lc['csv'], folder))
                         for lc in sectors],
                files=[f'phot/{key}.csv']) if sectors else dict(
                status='none', target=val.get('target'),
                message=val.get('error')
                or f'not in the fields of {val["mission"]}'))
            if sectors:
                told.append(f'{val["mission"]}: ' + ', '.join(
                    lc['label'].split()[-1] for lc in sectors)
                    + f' ({sum(len(lc["time"]) for lc in sectors)} points)')
        if space is not None:
            log('Kepler, K2, CoRoT: ' + ('; '.join(told) or 'the star is in '
                                         'none of their fields'), 'value')
        row('space', *space_state(manifest['archives']))

    # the archives all at once, not one after the other (DACE can take a
    #   minute, the papers too, the archive of the TNG half of one): each
    #   writes its own files and its own entry of the manifest
    asks = dict(variability=ask_variability)
    for key, wanted, func in (('dace', dace, ask_dace),
                              ('carmenes', carmenes, ask_carmenes),
                              ('published', vizier, ask_published),
                              ('tess', tess, ask_tess),
                              ('space', tess, ask_space)):
        if wanted:
            asks[key] = func
    names = dict(dace='DACE', carmenes='CARMENES DR1',
                 published='VizieR and the archives of the observatories',
                 tess='TESS', space='Kepler, K2, CoRoT')
    step(', '.join(names[key] for key in asks if key in names)
         + ', side by side')
    with ThreadPoolExecutor(max_workers=len(asks)) as pool:
        jobs = {key: pool.submit(func) for key, func in asks.items()}
        found = {key: job.result() for key, job in jobs.items()}
    # the velocities in the same order as ever: DACE, then CARMENES
    series = [found[key] for key in ('dace', 'carmenes')
              if found.get(key) is not None]
    # what each gave, in one line (the table of the GUI says more)
    told = []
    for key in ('dace', 'carmenes', 'published', 'tess', 'kepler', 'k2',
                'corot'):
        entry = manifest['archives'].get(key)
        if not entry or key in ('kepler', 'k2', 'corot') \
                and entry.get('status') != 'ok':
            continue
        sectors = entry.get('sectors') or []
        if entry.get('status') != 'ok':
            what = 'nothing'
        elif 'sectors' in entry:
            what = f'{len(sectors)} sector{"s" * (len(sectors) != 1)}'
        else:
            what = f'{int(entry.get("npoints") or 0)} points'
        told.append(dict(names, published='VizieR and the observatories',
                         kepler='Kepler', k2='K2', corot='CoRoT')[key]
                    + f': {what}')
    outcome('; '.join(told))
    if series:
        both = merge(series, name=target)
        write_rv(both, os.path.join(rvdir, 'all_rv.csv'))
        manifest['rv'] = dict(file='rv/all_rv.csv', npoints=both.n,
                              instruments=both.instruments)
    with open(os.path.join(folder, 'target.json'), 'w') as handle:
        json.dump(info, handle, indent=1, default=str)
    with open(os.path.join(folder, 'manifest.json'), 'w') as handle:
        json.dump(manifest, handle, indent=1, default=str)
    for key, val in manifest['archives'].items():
        if val['status'] == 'ok':
            what = (', '.join(f'{name} {num}' for name, num
                              in val['instruments'].items())
                    if 'instruments' in val else
                    f'{len(val["sectors"])} sectors')
            log(f'{key}: {what}', 'value')
        elif key in ('kepler', 'k2', 'corot', 'space') \
                and val['status'] == 'none':
            log(f'{key}: {val["message"]}')
        else:
            log(f'{key}: {val["status"]}, {val["message"]}', 'warn')
    return manifest


def state_of(entry: Optional[Dict[str, Any]]) -> Tuple[str, str]:
    """
    Where an archive is, from what the manifest says of it: ok and what
    it gave, none and why, failed and why (it could not be reached)

    :return: tuple, the state (ok, none, failed) and words
    """
    if not entry:
        return 'none', ''
    if entry.get('status') == 'ok':
        if entry.get('instruments'):
            return 'ok', (f'{int(entry.get("npoints") or 0)} velocities: '
                          + ', '.join(f'{name} {num}' for name, num
                                      in entry['instruments'].items()))
        sectors = entry.get('sectors') or []
        return 'ok', (f'{len(sectors)} sector{"s" * (len(sectors) != 1)}, '
                      f'{sum(int(one.get("npoints") or 0) for one in sectors)}'
                      f' points')
    return (('none' if entry.get('status') == 'none' else 'failed'),
            str(entry.get('message') or ''))


def space_state(archives: Dict[str, Any]) -> Tuple[str, str]:
    """where Kepler, K2 and CoRoT are, in one line: the missions that
    have the star, or none"""
    told, failed = [], ''
    for key, name in (('kepler', 'Kepler'), ('k2', 'K2'),
                      ('corot', 'CoRoT')):
        state, words = state_of(archives.get(key))
        if state == 'ok':
            told.append(f'{name}: {words}')
        elif state == 'failed':
            failed = words
    if 'space' in archives and not told:
        failed = failed or str(archives['space'].get('message') or '')
    if told:
        return 'ok', '; '.join(told)
    return ('failed', failed) if failed else (
        'none', 'the star is in none of their fields')


def _published(notes: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """what the manifest says of the published velocities of a star
    (koloa.published.fetch): the surveys and the papers on VizieR, the
    archives of the observatories"""
    kept = [note for note in notes if note.get('file')]
    insts: Dict[str, int] = {}
    for note in kept:
        for name, num in note['instruments'].items():
            insts[name] = insts.get(name, 0) + num
    return dict(
        status='ok' if kept else 'none',
        npoints=int(sum(insts.values())), instruments=insts,
        sources=[dict(reference=note['reference'],
                      catalogue=note['catalogue'], n=note['n'],
                      kind=note['kind']) for note in kept],
        files=[os.path.join('rv', 'published', note['file'])
               for note in kept],
        message=('' if kept else 'no published velocity of the star in '
                 'the surveys or its papers on VizieR, nor in the archives '
                 'of the observatories'))


def _wanted(sources: Optional[Sequence[str]] = None) -> List[str]:
    """the archives of the observatories a folder should have read: those
    of the sources asked, or all but those asked only when named"""
    from koloa.published import ARCHIVES, OPTIONAL
    return [key for key in ARCHIVES if (key not in OPTIONAL
                                        if sources is None
                                        else key in sources)]


def stale(folder: str, sources: Optional[Sequence[str]] = None) -> bool:
    """whether the folder of a gathered star is of before koloa kept the
    astrometry of the star with it, or did not read the archives of the
    observatories asked (koloa.published.ARCHIVES; sources: those of the
    sources asked, None for all but ESO and the CFHT): update() brings it
    up to date"""
    target = os.path.join(folder, 'target.json')
    if not os.path.exists(target):
        return False
    with open(target) as handle:
        if not json.load(handle).get('astrometry'):
            return True
    asked = os.path.join(folder, 'rv', 'published', 'sources.json')
    if not os.path.exists(asked):
        return False
    with open(asked) as handle:
        return not set(_wanted(sources)) <= set(json.load(handle))


def update(folder: str, have: Sequence[str] = (),
           sources: Optional[Sequence[str]] = None) -> bool:
    """
    Bring the folder of a star gathered earlier up to date, with nothing
    asked again of what it has: the astrometry of the star (target.json),
    and the archives of the observatories that koloa did not read then,
    beside its published velocities

    :param folder: str, the folder of the star
    :param have: the spectrographs you have files of (see gather)
    :param sources: the sources asked (see gather's vizier_sources): the
                    archives of the observatories among them, None for all
                    but ESO and the CFHT

    :return: bool, whether anything was added
    """
    from koloa.published import fetch as published_fetch
    from koloa.secular import star_astrometry, star_perspective
    usual = _wanted(sources)
    if not stale(folder, sources):
        return False
    target, mfile = (os.path.join(folder, name)
                     for name in ('target.json', 'manifest.json'))
    with open(target) as handle:
        info = json.load(handle)
    if not info.get('astrometry'):
        info['astrometry'] = star_astrometry(info)
        info['perspective'] = star_perspective(info)
        with open(target, 'w') as handle:
            json.dump(info, handle, indent=1, default=str)
    pdir = os.path.join(folder, 'rv', 'published')
    asked = os.path.join(pdir, 'sources.json')
    if os.path.exists(asked):
        with open(asked) as handle:
            last = json.load(handle)
        if not set(usual) <= set(last):
            # the sources it had, and the archives: VizieR is not asked
            notes = published_fetch(info, pdir, sources=list(
                dict.fromkeys(list(last) + usual)), have=have)
            if os.path.exists(mfile):
                with open(mfile) as handle:
                    manifest = json.load(handle)
                manifest.setdefault('archives', {})['published'] = \
                    _published(notes)
                with open(mfile, 'w') as handle:
                    json.dump(manifest, handle, indent=1, default=str)
    return True


def archive_berv(folder: str, times: np.ndarray, tol: float = 1e-4
                 ) -> np.ndarray:
    """
    The barycentric Earth radial velocity (BERV) of each velocity of a
    star's archives, from what they gave as it came: DACE's cal_berv
    (rv/dace/raw_*.csv, by rjd) and CARMENES' berv (rv/carmenes/raw.csv, by
    bjd); matched by time

    :param folder: str, the folder of a star
    :param times: np.ndarray, the times of the velocities [rjd]
    :param tol: float, the largest difference of time to match [days]

    :return: np.ndarray, the BERV [km/s] (nan where there is none)
    """
    import glob
    ref_t, ref_v = [], []
    sources = [(path, 'rjd', 'cal_berv', 0.0) for path in
               sorted(glob.glob(os.path.join(folder, 'rv', 'dace',
                                             'raw_*.csv')))]
    sources.append((os.path.join(folder, 'rv', 'carmenes', 'raw.csv'),
                    'bjd', 'berv', 2400000.0))
    for path, tcol, bcol, offset in sources:
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            for row in csv.DictReader(handle):
                tt, bb = _float(row.get(tcol)), _float(row.get(bcol))
                if np.isfinite(tt) and np.isfinite(bb):
                    ref_t.append(tt - offset)
                    ref_v.append(bb)
    times = np.asarray(times, dtype=float)
    out = np.full(len(times), np.nan)
    if not ref_t or not len(times):
        return out
    order = np.argsort(ref_t)
    ref_t, ref_v = np.asarray(ref_t)[order], np.asarray(ref_v)[order]
    idx = np.clip(np.searchsorted(ref_t, times), 1, max(len(ref_t) - 1, 1))
    idx = np.minimum(idx, len(ref_t) - 1)
    left = np.maximum(idx - 1, 0)
    pick = np.where(np.abs(times - ref_t[left]) <= np.abs(ref_t[idx] - times),
                    left, idx)
    near = np.abs(ref_t[pick] - times) <= tol
    out[near] = ref_v[pick][near]
    return out


def _told(entry: Dict[str, Any]) -> str:
    """what an archive gave, in words: its points, by instrument, or its
    sectors; or why nothing"""
    if entry.get('status') != 'ok':
        return f'nothing: {entry.get("message") or entry.get("status")}'
    if 'sectors' in entry:
        sectors = entry['sectors']
        return (f'{len(sectors)} sector{"s" if len(sectors) != 1 else ""}, '
                f'{sum(item["npoints"] for item in sectors)} points ('
                + ', '.join(f's{item["sector"]} {item["pipeline"]}'
                            for item in sectors) + ')')
    insts = entry.get('instruments') or {}
    return (f'{entry.get("npoints", sum(insts.values()))} points'
            + (': ' + ', '.join(f'{name} {num}' for name, num
                                in insts.items()) if insts else ''))


def load(folder: str, photometry: bool = True) -> Dict[str, Any]:
    """
    What gather() put in a folder

    :param folder: str, the folder of a star
    :param photometry: bool, read the TESS sectors too (False: the
                       velocities only, faster); a sector the manifest lists
                       but the folder no longer has is left out

    :return: dict, target (target.json), manifest, rv (RVData of every
             velocity, None if there is none), tess (list of dict: sector,
             pipeline, time, flux, err)
    """
    def read_json(name):
        path = os.path.join(folder, name)
        if not os.path.exists(path):
            return {}
        with open(path) as handle:
            return json.load(handle)
    manifest = read_json('manifest.json')
    target = read_json('target.json')
    rvfile = os.path.join(folder, 'rv', 'all_rv.csv')
    rv = (RVData.from_csv(rvfile, inst='inst',
                          name=manifest.get('target'))
          if os.path.exists(rvfile) else None)
    sectors = []
    for item in (manifest.get('archives', {}).get('tess', {})
                 .get('sectors', []) if photometry else []):
        path = os.path.join(folder, item['file'])
        if not os.path.exists(path):
            log(f'TESS sector {item["sector"]}: {item["file"]} is not in '
                f'{folder}, left out', 'warn')
            continue
        with open(path) as handle:
            rows = list(csv.DictReader(handle))
        sectors.append(dict(sector=item['sector'], pipeline=item['pipeline'],
                            time=np.array([_float(r['rjd']) for r in rows]),
                            flux=np.array([_float(r['flux']) for r in rows]),
                            err=np.array([_float(r['sflux']) for r in rows])))
    return dict(target=target, manifest=manifest, rv=rv, tess=sectors)


# =============================================================================
# End of code
# =============================================================================

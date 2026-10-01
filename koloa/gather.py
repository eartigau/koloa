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
    phot/tess/fits/      the TESS light curves of MAST, one per sector (the
                         best of SPOC 2-minute, TESS-SPOC, QLP)
    phot/tess/s<sector>_<pipeline>.csv   each sector: rjd, flux, sflux [ppt]
    phot/tess.csv        every sector together, with its sector and pipeline

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
from typing import Any, Dict, List, Optional

import numpy as np

from koloa.data import RVData, merge
from koloa.log import log, step

# =============================================================================
# Define variables
# =============================================================================
#: the TAP service that serves CARMENES DR1
GAVO_TAP = 'https://dc.g-vo.org/tap/sync'
#: SIMBAD's TAP service
SIMBAD_TAP = 'https://simbad.cds.unistra.fr/simbad/sim-tap/sync'
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
def variability(main: str) -> List[Dict[str, Any]]:
    """
    The periods of variability SIMBAD lists for a star (its mesVar table):
    type (ROT for a rotation), period [days] and reference

    :param main: str, SIMBAD's main name of the star
    """
    name = main.replace("'", "''")
    rows = _tap(SIMBAD_TAP, (
        'SELECT v.vartyp, v.period, v.bibcode FROM mesVar AS v JOIN ident '
        f"AS i ON v.oidref = i.oidref WHERE i.id = '{name}'"))
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
        path = os.path.join(folder, f'raw_{name}.csv')
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


def carmenes_star(ra: float, dec: float) -> Optional[Dict[str, Any]]:
    """
    The CARMENES DR1 star at a position (its Karmn name, number of
    velocities, rotation period...), None if DR1 has none there

    :param ra: float, J2000 [degrees]
    :param dec: float, J2000 [degrees]
    """
    rows = _tap(GAVO_TAP, (
        'SELECT carmenes_id, name, nobs, nspect, mass, radius, teff, p_rot, '
        'p_rot_source, survey, remarks, ra, dec FROM carmenes.objects WHERE '
        f"1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, "
        f'{dec:.6f}, {CARMENES_RADIUS}))'))
    return dict(rows[0]) if rows else None


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


# -----------------------------------------------------------------------------
# All of them
# -----------------------------------------------------------------------------
def gather(target: str, root: str = '.', dace: bool = True,
           carmenes: bool = True, tess: bool = True,
           api_key: Any = None, refresh: bool = False
           ) -> Dict[str, Any]:
    """
    Everything public about a star, in root/<target> (see the module)

    Each archive is asked in turn; one that fails is reported in the
    manifest and the others go on.

    :param target: str, a name SIMBAD knows (GJ 436, Ross 905, TOI-700...)
    :param root: str, the folder of the folders of the stars
    :param dace: bool, the velocities of DACE
    :param carmenes: bool, the velocities of CARMENES DR1
    :param tess: bool, the light curves of TESS
    :param api_key: str, None or False: a DACE API key, None to look for
                    one (DACE_API_KEY, ~/.dacerc), False for the public
                    data only (the key is never written anywhere)
    :param refresh: bool, ask again for what is already there

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
    try:
        info['variability'] = variability(ident['main'])
    except Exception as err:  # SIMBAD's TAP is a help, not a need
        info['variability'] = []
        log(f'SIMBAD variability: {err}', 'warn')
    manifest: Dict[str, Any] = dict(
        target=target, simbad=ident['main'], folder=os.path.abspath(folder),
        created=_time.strftime('%Y-%m-%d %H:%M:%S'), archives={})
    series = []

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
    if dace:
        step('DACE')
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
            series.append(data)
        elif 'dace' not in manifest['archives']:
            manifest['archives']['dace'] = dict(
                status='none', message='DACE has no public velocities under '
                'the names of the star')
    if carmenes:
        step('CARMENES DR1')
        cdir = os.path.join(rvdir, 'carmenes')
        star = (attempt('carmenes', lambda: carmenes_star(ident['ra'],
                                                          ident['dec']))
                if ident.get('ra') is not None else None)
        if star is not None:
            info['carmenes'] = star
            data = attempt('carmenes', lambda: carmenes_rv(star, target, cdir,
                                                           refresh=refresh))
            if data is not None:
                manifest['archives']['carmenes'] = dict(
                    _summary(data, [os.path.join(cdir, 'raw.csv'),
                                    os.path.join(cdir, 'CARMENES.csv')],
                             folder), karmn=star['carmenes_id'])
                series.append(data)
        if 'carmenes' not in manifest['archives']:
            manifest['archives']['carmenes'] = dict(
                status='none', message='not in CARMENES DR1, or no '
                'velocity corrected for the nightly zero points')
    if tess:
        step('TESS')
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
        else:
            log(f'{key}: {val["status"]}, {val["message"]}', 'warn')
    return manifest


def load(folder: str) -> Dict[str, Any]:
    """
    What gather() put in a folder

    :param folder: str, the folder of a star

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
    for item in manifest.get('archives', {}).get('tess', {}).get('sectors',
                                                                 []):
        with open(os.path.join(folder, item['file'])) as handle:
            rows = list(csv.DictReader(handle))
        sectors.append(dict(sector=item['sector'], pipeline=item['pipeline'],
                            time=np.array([_float(r['rjd']) for r in rows]),
                            flux=np.array([_float(r['flux']) for r in rows]),
                            err=np.array([_float(r['sflux']) for r in rows])))
    return dict(target=target, manifest=manifest, rv=rv, tess=sectors)


# =============================================================================
# End of code
# =============================================================================

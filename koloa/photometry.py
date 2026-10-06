#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The light curves of a star from Kepler, K2 and CoRoT, beside those of TESS
(koloa.tess): a transit is looked for in them the same way (koloa.transit).

A star is in them only when it lay in their fields: Kepler's one field in
Cygnus and Lyra (2009 to 2013, quarters of three months), the campaigns of
K2 along the ecliptic (2014 to 2018, 80 days each), and the two eyes of
CoRoT towards the centre and the anticentre of the Galaxy (2007 to 2012,
runs of 20 to 150 days). They are asked by position.

- kepler_light_curves(ident): the long-cadence (30 min) light curves of
  Kepler and of K2 at MAST, the nearest target within a few arcseconds:
  the PDCSAP flux of each quarter of Kepler (SAP where there is none); for
  each campaign of K2, the K2SFF light curve when there is one (Vanderburg
  & Johnson 2014: corrected for the roll of the spacecraft, which the
  mission's PDCSAP flux of the early campaigns keeps), else the mission's.
- corot_light_curves(ident): the CoRoT N2 light curves (their version 4.4,
  VizieR B/corot, the files at the CDS): for a star of the exoplanet
  channel (faint stars) the flux corrected for systematics (512 s), for a
  star of the seismology channel (bright stars) the regular flux (32 s,
  binned to 512 s); the points flagged (the South Atlantic Anomaly, jumps)
  left out.
- light_curves(ident): all of them.

Each stretch (a quarter, a campaign, a run) is normalised on its own, its
flux in ppt about its median, its times in BJD - 2400000, like a sector of
TESS. A stretch has a number as a sector of TESS has (koloa.transit splits
a light curve by it): 1000 + the quarter of Kepler, 2000 + the campaign of
K2, 3000 + the rank of the run of CoRoT; label() names it.

The FITS files of Kepler and K2 are kept (a few hundred kB each); those of
CoRoT are large (tens of MB, three colours at 32 s): what is read of them
is kept as CSV, not the files.

Created on 2026-10-06

@author: artigau
"""
import csv
import os
import re
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from koloa.log import log
from koloa.paths import cache

# =============================================================================
# Define variables
# =============================================================================
#: where the light curves are kept (they are public, and shared by runs)
CACHE = cache()
#: the missions, in the order they are tried after TESS
MISSIONS = ('kepler', 'k2', 'corot')
#: their names
NAMES = dict(tess='TESS', kepler='Kepler', k2='K2', corot='CoRoT')
#: the number of a stretch: this plus the quarter, the campaign, the rank
#: of the run (a sector of TESS is below the first)
BASE = dict(kepler=1000, k2=2000, corot=3000)
#: a target of Kepler or K2 is the star when within this [arcsec] (their
#: positions are of 2000, a star may have moved since by a few arcseconds)
MATCH_ARCSEC = 8.0
#: the same for CoRoT (its positions are those of its input catalogue)
COROT_ARCSEC = 10.0
#: the CoRoT observation log at VizieR, and its files at the CDS
VIZIER = 'https://vizier.cds.unistra.fr/viz-bin/asu-tsv'
COROT_FILES = 'https://cdsarc.cds.unistra.fr/ftp/B/corot/files/'
#: the tables of the log: the exoplanet channel, the seismology channel
COROT_TABLES = ('B/corot/Faint_star', 'B/corot/Bright_star')
#: the bright stars of CoRoT (32 s) binned to this [s]
COROT_BIN = 512.0


# =============================================================================
# Define functions
# =============================================================================
def mission_of(sector: float) -> str:
    """the mission of a stretch, from its number"""
    sector = int(sector)
    for name in ('corot', 'k2', 'kepler'):
        if sector >= BASE[name]:
            return name
    return 'tess'


def label(sector: float, run: Optional[str] = None) -> str:
    """
    A stretch by its name: 'sector 22' (TESS), 'Kepler Q3', 'K2 C5',
    'CoRoT LRa01' (its run, when given: its number does not say it)
    """
    sector = int(sector)
    mission = mission_of(sector)
    if mission == 'tess':
        return f'sector {sector}'
    if mission == 'kepler':
        return f'Kepler Q{sector - BASE["kepler"]}'
    if mission == 'k2':
        return f'K2 C{sector - BASE["k2"]}'
    return f'CoRoT {run or sector - BASE["corot"]}'


def _normalised(time: np.ndarray, flux: np.ndarray,
                err: Optional[np.ndarray] = None) -> Optional[Dict[str, Any]]:
    """a stretch in ppt about its median (its error the scatter of its
    point-to-point differences when it has none), or None with too few
    points"""
    good = np.isfinite(time) & np.isfinite(flux)
    if good.sum() < 100:
        return None
    time, flux = time[good], flux[good]
    med = float(np.median(flux))
    if not med > 0:
        return None
    rel = 1e3 * (flux / med - 1.0)
    if err is not None:
        sig = 1e3 * np.asarray(err, float)[good] / med
        bad = ~np.isfinite(sig) | (sig <= 0)
    else:
        sig, bad = np.zeros(len(rel)), np.ones(len(rel), bool)
    if bad.any():
        diff = np.diff(rel)
        sig = np.where(bad, 1.4826 * np.median(np.abs(
            diff - np.median(diff))) / np.sqrt(2.0), sig)
    order = np.argsort(time)
    return dict(time=time[order], flux=rel[order], err=sig[order])


def _download(url: str, path: str, timeout: float = 300.0) -> None:
    """a file fetched and kept"""
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        content = resp.read()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + '.part', 'wb') as handle:
        handle.write(content)
    os.replace(path + '.part', path)


# -----------------------------------------------------------------------------
# Kepler and K2
# -----------------------------------------------------------------------------
def _read_kepler(path: str) -> Optional[Dict[str, Any]]:
    """a long-cadence light curve of Kepler or K2: its quarter or campaign,
    its PDCSAP flux (SAP where it has none), the cadences flagged left out;
    its times BJD - 2400000"""
    from astropy.io import fits
    with fits.open(path) as hdul:
        head, tbl = hdul[0].header, hdul[1].data
        cols = tbl.columns.names
        ref = (hdul[1].header.get('BJDREFI', 2454833)
               + hdul[1].header.get('BJDREFF', 0.0))
        time_ = np.asarray(tbl['TIME'], float) + ref - 2400000.0
        qual = (np.asarray(tbl['SAP_QUALITY']) if 'SAP_QUALITY' in cols
                else np.zeros(len(time_), int))
        column = 'PDCSAP_FLUX'
        flux = np.asarray(tbl[column], float)
        if np.isfinite(flux).sum() < 100 and 'SAP_FLUX' in cols:
            column = 'SAP_FLUX'
            flux = np.asarray(tbl[column], float)
        err = (np.asarray(tbl[column + '_ERR'], float)
               if column + '_ERR' in cols else None)
        quarter, campaign = head.get('QUARTER'), head.get('CAMPAIGN')
        exposure = 60.0 * 29.4244 if head.get('OBSMODE', 'long cadence') \
            .startswith('long') else 58.85
    keep = qual == 0
    out = _normalised(np.where(keep, time_, np.nan), flux, err)
    if out is None:
        return None
    out.update(column=column, exposure=exposure,
               quarter=None if quarter is None else int(quarter),
               campaign=None if campaign is None else int(campaign))
    return out


def _read_k2sff(path: str) -> Optional[Dict[str, Any]]:
    """a K2SFF light curve (Vanderburg & Johnson 2014: the flux of K2
    corrected for the roll of the spacecraft, which the mission's own
    PDCSAP flux of the early campaigns keeps): its best aperture, the
    cadences where the spacecraft moved left out; times BJD - 2400000"""
    from astropy.io import fits
    with fits.open(path) as hdul:
        tbl = hdul['BESTAPER'].data
        ref = (hdul['BESTAPER'].header.get('BJDREFI', 2454833)
               + hdul['BESTAPER'].header.get('BJDREFF', 0.0))
        time_ = np.asarray(tbl['T'], float) + ref - 2400000.0
        flux = np.asarray(tbl['FCOR'], float)
        still = np.asarray(tbl['MOVING']) == 0
    out = _normalised(np.where(still, time_, np.nan), flux)
    if out is None:
        return None
    out.update(column='FCOR (K2SFF)', exposure=60.0 * 29.4244, quarter=None,
               campaign=None)
    return out


#: what MAST has at a position, once asked: (ra, dec) -> kepler, k2
_PRODUCTS: Dict[Any, Dict[str, Dict[str, Any]]] = {}


def kepler_products(ident: Dict[str, Any], timeout: float = 120.0
                    ) -> Dict[str, Dict[str, Any]]:
    """
    What MAST has of a star from Kepler and K2, without fetching it: the
    long-cadence light curve of each quarter and campaign of the nearest
    target within MATCH_ARCSEC (for a campaign of K2, K2SFF's before the
    mission's); asked once per position

    :return: dict, kepler and k2: each dict(target, files: {tag: product
             of MAST (productFilename, dataURI), with sff for K2SFF's})
    """
    from koloa.tess import _mast
    out = {key: dict(target=None, files={}) for key in ('kepler', 'k2')}
    if ident.get('ra') is None or ident.get('dec') is None:
        return out
    where = (round(float(ident['ra']), 5), round(float(ident['dec']), 5))
    if where in _PRODUCTS:
        return _PRODUCTS[where]
    log(f'Kepler, K2: asking MAST for the light curves at '
        f'{ident["ra"]:.5f}, {ident["dec"]:+.5f}')
    query = dict(
        columns='obsid,obs_collection,provenance_name,target_name,'
                'sequence_number,t_exptime,obs_id,distance',
        filters=[dict(paramName='dataproduct_type', values=['timeseries']),
                 dict(paramName='obs_collection',
                      values=['Kepler', 'K2', 'HLSP'])],
        position=f'{ident["ra"]}, {ident["dec"]}, '
                 f'{MATCH_ARCSEC / 3600.0}')
    try:
        rows = _mast('Mast.Caom.Filtered.Position', query, timeout=timeout)
    except OSError:
        # MAST drops a connection now and then: asked once more
        rows = _mast('Mast.Caom.Filtered.Position', query, timeout=timeout)
    for key, collection in (('kepler', 'Kepler'), ('k2', 'K2')):
        # the mission's own long cadence (and, for K2, the K2SFF light
        #   curves made of it), of the nearest target
        mine = [row for row in rows
                if (row.get('obs_collection') == collection
                    and row.get('provenance_name') == collection
                    and abs(float(row.get('t_exptime') or 0) - 1800) < 100)
                or (key == 'k2' and row.get('provenance_name') == 'K2SFF')]
        if not mine:
            continue
        near = min(float(row.get('distance') or 0) for row in mine)
        mine = [row for row in mine
                if float(row.get('distance') or 0) <= near + 1.0]
        out[key]['target'] = next(
            (row['target_name'] for row in mine
             if row.get('provenance_name') == collection),
            mine[0]['target_name'])
        files = _mast('Mast.Caom.Products', dict(
            obsid=','.join(str(row['obsid']) for row in mine)),
            timeout=timeout)
        # one file per quarter or campaign: K2SFF's before the mission's
        chosen = out[key]['files']
        for item in sorted(files, key=lambda one: one['productFilename']):
            name = item['productFilename']
            if not name.endswith('_llc.fits'):
                continue
            found = re.search(r'-c(\d+)_', name)
            if name.startswith('hlsp_k2sff') and found:
                chosen[f'c{int(found.group(1))}'] = dict(item, sff=True)
            elif name.startswith('ktwo') and found:
                chosen.setdefault(f'c{int(found.group(1))}', item)
            elif name.startswith('kplr'):
                chosen[name] = item
    _PRODUCTS[where] = out
    return out


def kepler_light_curves(ident: Dict[str, Any], folder: Optional[str] = None,
                        refresh: bool = False, timeout: float = 120.0,
                        missions: Sequence[str] = ('kepler', 'k2')
                        ) -> Dict[str, Dict[str, Any]]:
    """
    The long-cadence light curves of a star from Kepler and from K2 (MAST,
    by position: the nearest target within MATCH_ARCSEC), one per quarter
    or campaign

    :param ident: dict, the star (koloa.archive.resolve: ra, dec [deg])
    :param folder: str or None, where the FITS files are kept
                   (CACHE/<mission>/<target> when None)
    :param refresh: bool, download the files again
    :param missions: which of kepler and k2 are fetched

    :return: dict, kepler and k2: each dict(mission, target, sectors (list
             of dict: sector, label, provenance, exposure, column, time,
             flux, err, file)); no sector when the star is not in it
    """
    from koloa.tess import DOWNLOAD
    out = {key: dict(mission=NAMES[key], target=None, sectors=[])
           for key in ('kepler', 'k2')}
    products = kepler_products(ident, timeout)
    for key, collection in (('kepler', 'Kepler'), ('k2', 'K2')):
        target = products[key]['target']
        out[key]['target'] = target
        if key not in missions or not products[key]['files']:
            continue
        where = folder or os.path.join(CACHE, key, str(target))
        for tag, item in products[key]['files'].items():
            name = item['productFilename']
            path = os.path.join(where, name)
            try:
                if refresh or not os.path.exists(path):
                    _download(DOWNLOAD + urllib.parse.quote(
                        item['dataURI'], safe=':/'), path, timeout)
                lc = _read_k2sff(path) if item.get('sff') \
                    else _read_kepler(path)
            except Exception as err:  # a broken file: the next one
                log(f'{collection}: {name} unreadable ({err})', 'warn')
                continue
            if lc is None:
                continue
            number = lc['quarter'] if key == 'kepler' else int(tag[1:])
            if number is None:
                continue
            sector = BASE[key] + int(number)
            lc.update(sector=sector, label=label(sector),
                      provenance='K2SFF' if item.get('sff') else collection,
                      file=path)
            out[key]['sectors'].append(lc)
        out[key]['sectors'].sort(key=lambda lc: lc['time'][0])
        if out[key]['sectors']:
            log(f'{collection}: {len(out[key]["sectors"])} '
                f'{"quarters" if key == "kepler" else "campaigns"} of '
                f'{target} ('
                + ', '.join(lc['label'].split()[-1]
                            for lc in out[key]['sectors']) + ')', 'value')
    return out


# -----------------------------------------------------------------------------
# CoRoT
# -----------------------------------------------------------------------------
#: the runs of CoRoT at a position, once asked
_COROT_LOG: Dict[Any, List[Dict[str, str]]] = {}


def _corot_log(ident: Dict[str, Any], timeout: float = 60.0
               ) -> List[Dict[str, str]]:
    """the runs of CoRoT that hold a position (VizieR B/corot, its two
    tables): CoRoT (its number of the star), Run, FileName, date1"""
    where = (round(float(ident['ra']), 5), round(float(ident['dec']), 5))
    if where in _COROT_LOG:
        return _COROT_LOG[where]
    found = []
    for table in COROT_TABLES:
        url = VIZIER + '?' + urllib.parse.urlencode({
            '-source': table, '-c': f'{ident["ra"]} {ident["dec"]:+}',
            '-c.rs': f'{COROT_ARCSEC}', '-out': 'CoRoT,Run,FileName,date1',
            '-out.max': '50'})
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            lines = [line for line in resp.read().decode('utf-8', 'replace')
                     .splitlines() if line and not line.startswith('#')]
        if len(lines) < 4:
            continue
        names = [name.strip() for name in lines[0].split('\t')]
        for line in lines[3:]:
            row = dict(zip(names, (cell.strip() for cell in line.split('\t'))))
            if row.get('FileName'):
                found.append(dict(row, table=table))
    _COROT_LOG[where] = sorted(found, key=lambda row: row.get('date1', ''))
    return _COROT_LOG[where]


def _read_corot(path: str) -> Optional[Dict[str, Any]]:
    """
    A CoRoT N2 file: the flux corrected for systematics of a faint star
    (SYSTEMATIC, 512 s; else its filled or raw white flux), the regular
    flux of a bright one (BARREG, binned to COROT_BIN), its valid points
    (status 0); its times are BJD - 2400000 already
    """
    from astropy.io import fits
    tries = (('SYSTEMATIC', 'DATEBARTT', 'WHITEFLUXSYS', 'STATUSSYS'),
             ('BARREG', 'DATEBARREGTT', 'FLUXBARREG', 'STATUSBARREG'),
             ('BARFILL', 'DATEBARTT', 'WHITEFLUXFIL', 'STATUSFIL'),
             ('BAR', 'DATEBARTT', 'WHITEFLUX', 'STATUS'),
             ('BAR', 'DATEBARTT', 'FLUXBAR', 'STATUSBAR'))
    with fits.open(path) as hdul:
        have = {hdu.name: hdu for hdu in hdul[1:]}
        for ext, tcol, fcol, scol in tries:
            if ext not in have or fcol not in have[ext].columns.names:
                continue
            tbl = have[ext].data
            time_ = np.asarray(tbl[tcol], float)
            flux = np.asarray(tbl[fcol], float)
            keep = np.asarray(tbl[scol]) == 0
            if keep.sum() < 100:
                continue
            time_, flux = time_[keep], flux[keep]
            step = float(np.median(np.diff(np.sort(time_)))) * 86400.0
            if step < 0.5 * COROT_BIN:
                # a bright star at 32 s: means over 512 s
                idx = np.floor((time_ - time_.min()) * 86400.0
                               / COROT_BIN).astype(int)
                num = np.bincount(idx)
                full = num >= 0.5 * COROT_BIN / max(step, 1.0)
                time_ = (np.bincount(idx, time_)[full] / num[full])
                flux = (np.bincount(idx, flux)[full] / num[full])
                step = COROT_BIN
            out = _normalised(time_, flux)
            if out is not None:
                out.update(column=f'{ext}.{fcol}', exposure=step)
                return out
    return None


def corot_light_curves(ident: Dict[str, Any], folder: Optional[str] = None,
                       refresh: bool = False, timeout: float = 600.0
                       ) -> Dict[str, Any]:
    """
    The CoRoT light curves of a star, one per run (VizieR's observation
    log by position, the N2 files of the CDS; what is read of each file is
    kept as CSV, the file itself is not: tens of MB)

    :param ident: dict, the star (ra, dec [deg])
    :param folder: str or None, where the light curves are kept
                   (CACHE/corot/<its CoRoT number> when None)

    :return: dict, mission, target (its CoRoT number), sectors (list of
             dict: sector, label, run, provenance, exposure, column, time,
             flux, err, file)
    """
    out: Dict[str, Any] = dict(mission=NAMES['corot'], target=None,
                               sectors=[])
    if ident.get('ra') is None or ident.get('dec') is None:
        return out
    log(f'CoRoT: asking VizieR for its runs at {ident["ra"]:.5f}, '
        f'{ident["dec"]:+.5f}')
    rows = _corot_log(ident)
    if not rows:
        return out
    out['target'] = rows[0].get('CoRoT')
    where = folder or os.path.join(CACHE, 'corot', str(out['target']))
    for rank, row in enumerate(rows):
        run = row.get('Run') or f'run{rank + 1}'
        path = os.path.join(where, f'{run}_{row.get("CoRoT")}.csv')
        try:
            if refresh or not os.path.exists(path):
                log(f'CoRoT: fetching the run {run} (tens of MB)')
                fits_path = os.path.join(where, os.path.basename(
                    row['FileName']))
                _download(COROT_FILES + row['FileName'], fits_path, timeout)
                lc = _read_corot(fits_path)
                os.remove(fits_path)
                if lc is None:
                    continue
                with open(path, 'w', newline='') as handle:
                    writer = csv.writer(handle)
                    writer.writerow(['rjd', 'flux', 'sflux', 'column',
                                     'exposure'])
                    for k, (tt, ff, ee) in enumerate(zip(
                            lc['time'], lc['flux'], lc['err'])):
                        writer.writerow([f'{tt:.6f}', f'{ff:.5f}',
                                         f'{ee:.5f}']
                                        + ([lc['column'], lc['exposure']]
                                           if k == 0 else ['', '']))
            else:
                with open(path, newline='') as handle:
                    read = list(csv.DictReader(handle))
                lc = dict(time=np.array([float(r['rjd']) for r in read]),
                          flux=np.array([float(r['flux']) for r in read]),
                          err=np.array([float(r['sflux']) for r in read]),
                          column=read[0]['column'],
                          exposure=float(read[0]['exposure']))
        except Exception as err:  # a run that cannot be had: the next one
            log(f'CoRoT: run {run} unreadable ({err})', 'warn')
            continue
        sector = BASE['corot'] + rank + 1
        lc.update(sector=sector, label=label(sector, run), run=run,
                  provenance=run, file=path)
        out['sectors'].append(lc)
    if out['sectors']:
        log(f'CoRoT: {len(out["sectors"])} runs of {out["target"]} ('
            + ', '.join(lc['run'] for lc in out['sectors']) + ')', 'value')
    return out


def available(ident: Dict[str, Any]) -> Dict[str, Any]:
    """
    Which of Kepler, K2 and CoRoT have a star, without fetching their
    light curves (MAST and VizieR asked by position, once)

    :return: dict, kepler, k2, corot: how many quarters, campaigns, runs
             (None for a mission that could not be asked)
    """
    out: Dict[str, Any] = {}
    try:
        products = kepler_products(ident)
        out.update({key: len(products[key]['files'])
                    for key in ('kepler', 'k2')})
    except Exception as err:
        log(f'Kepler, K2: MAST could not be asked ({err})', 'warn')
        out.update(kepler=None, k2=None)
    try:
        out['corot'] = (len(_corot_log(ident))
                        if ident.get('ra') is not None else 0)
    except Exception as err:
        log(f'CoRoT: VizieR could not be asked ({err})', 'warn')
        out['corot'] = None
    return out


def light_curves(ident: Dict[str, Any], refresh: bool = False,
                 folders: Optional[Dict[str, str]] = None,
                 missions: Sequence[str] = MISSIONS
                 ) -> Dict[str, Dict[str, Any]]:
    """
    The light curves of a star from Kepler, K2 and CoRoT (those it is in);
    a mission that cannot be asked is left without any, with its error

    :param folders: dict or None, where each mission's files are kept
                    (kepler, k2, corot)
    :param missions: which of them are fetched

    :return: dict, kepler, k2, corot: each dict(mission, target, sectors)
    """
    folders = folders or {}
    out = {key: dict(mission=NAMES[key], target=None, sectors=[])
           for key in MISSIONS}
    for key in ('kepler', 'k2'):
        if key not in missions:
            continue
        try:
            out[key] = kepler_light_curves(ident, folders.get(key), refresh,
                                           missions=(key,))[key]
        except Exception as err:
            log(f'{NAMES[key]}: MAST could not be asked ({err})', 'warn')
            out[key]['error'] = f'{type(err).__name__}: {err}'
    if 'corot' in missions:
        try:
            out['corot'] = corot_light_curves(ident, folders.get('corot'),
                                              refresh=refresh)
        except Exception as err:
            log(f'CoRoT: VizieR could not be asked ({err})', 'warn')
            out['corot']['error'] = f'{type(err).__name__}: {err}'
    return out


# =============================================================================
# End of code
# =============================================================================

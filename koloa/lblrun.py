#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LBL velocities of a star from the public spectra of the archives: the
spectra fetched, and lbl-rapide (the branch of LBL made faster; Artigau
et al. 2022) run on them.

koloa never imports LBL: LBL needs versions of its own (Python 3.12,
numpy 2.3, numba, astropy 7), in an environment of its own (lbl-rapide by
default). koloa writes the parameters of the run in a short script, as one
writes a wrap file, and runs it with that Python, in a process of its own.
The rdb LBL makes is kept in the folder of the star, under its name, one
file per dataset koloa knows (archives/GJ_436/lbl/GJ_436_HARPSN_lbl.rdb;
SOPHIE and SOPHIE+ apart), to be given to the analysis as a file; the
names of their datasets are in labels.json beside them.

    koloa "GJ 436" --lbl HARPSN SOPHIE --archives archives --lbl-dir lbl_data

An instrument is one or more groups, a run of LBL each (its template, its
mask): the spectra of one setting (SOPHIE's HR and HE apart). Each class
of LBL has a folder of its own in lbl-dir (lbl-dir/HARPSN_ESO,
lbl-dir/SOPHIE_None), with the models of LBL shared between them
(lbl-dir/models, a few hundred MB, fetched by LBL the first time).

The instruments:

HARPSN: the espdr S2D_A files of the archive of the TNG (blaze corrected,
on a barycentric wavelength grid: LBL's HARPSN ESO source, no calibration
file), 10 MB each, of the exposures koloa keeps of the star from their
CCFs (koloa.tng.series: those on the star, with its mask, the companion of
a binary left out by its velocity); the telluric lines cleaned by LBL
first.

SOPHIE: the e2ds of the archive of the OHP (1.3 MB each once kept), of the exposures
koloa keeps from the table of its CCFs (koloa.observatories.sophie: past
their year of embargo, a programme still protected left out), HR and HE
apart; the archive puts the blaze beside the spectrum in one file: koloa
writes it as the calibration file LBL reads (calib/<seq>_blaze_A.fits,
named in the header). SOPHIE and SOPHIE+ (the fibres of June 2011) are
two datasets of the rdb.

SPIROU: the t.fits of APERO at the CADC (the telluric lines taken out by
APERO; 18 MB each), those public, pointed at where the star was; LBL's
SPIROU CADC source reads their blaze in them.

NIRPS: the 2D spectra of ESO (in the ancillary tar of each product, its
files fetched alone by Range requests): the S2D with its blaze divided by
ESO's model of the telluric lines (S2D_TELL_SPECTRUM_A: the S2D of ESO are
not corrected for them, and LBL does not clean them for NIRPS), the lines
deeper than TELLURIC_FLOOR left out, and the blaze of the night as the
calibration file the header names; HE and HA apart (LBL's NIRPS_HE and
NIRPS_HA ESO).

ESPRESSO: the S2D_A of ESO (blaze corrected, 65 MB each), a run per mode
and binning (HR 2x1, UHR...), the telluric lines cleaned by LBL first;
ESPRESSO18 and ESPRESSO19 (the intervention of 2019) are two datasets.

Created on 2026-10-14

@author: artigau
"""
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import time as _time
import urllib.error
import uuid
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from koloa.log import log, outcome, step
from koloa.paths import config

# =============================================================================
# Define variables
# =============================================================================
#: your choices: the Python of LBL, the folder of its data, its processes
SETTINGS = config('lbl.json')
#: where an environment called lbl-rapide may be, when none is chosen
PYTHONS = [os.path.join(os.path.expanduser(base), 'envs', 'lbl-rapide',
                        'bin', 'python')
           for base in ('/opt/anaconda3', '/opt/miniconda3', '~/anaconda3',
                        '~/miniconda3', '~/miniforge3', '~/mambaforge',
                        '/opt/homebrew/Caskroom/miniforge/base')]
#: the files of the archive of the TNG, and the e2ds of SOPHIE's
TNG_FILES = 'http://archives.ia2.inaf.it/files/tng/'
SOPHIE_E2DS = ('http://atlas.obs-hp.fr/sophie/sophie.cgi?c=i&a=mime:'
               'application/fits&o=sophie:%5Be2ds,{seq}%5D')
#: the instruments: what they are called, where their spectra are from,
#: about how large one is [MB]
INSTRUMENTS: Dict[str, Dict[str, Any]] = {
    'HARPSN': dict(label='HARPS-N', size=9.8, classes=['HARPSN_ESO'],
                   keys=[''], archive='the archive of the TNG (espdr S2D_A)'),
    'SOPHIE': dict(label='SOPHIE', size=4.2, classes=['SOPHIE_None'],
                   keys=['', 'HE'], archive='the archive of the OHP (e2ds)'),
    'SPIROU': dict(label='SPIRou', size=18.4, classes=['SPIROU_CADC'],
                   keys=[''], archive='the CADC (APERO t.fits)'),
    'NIRPS': dict(label='NIRPS', size=24.0,
                  classes=['NIRPS_HE_ESO', 'NIRPS_HA_ESO'], keys=['', 'HA'],
                  archive='the archive of ESO (S2D, telluric model)'),
    'ESPRESSO': dict(label='ESPRESSO', size=65.0, classes=['ESPRESSO_None'],
                     keys=['', 'HR11', 'HR21', 'HR42', 'UHR11', 'MR42', 'MR84'],
                     archive='the archive of ESO (S2D_A)'),
}
#: the deepest telluric line kept in a spectrum of NIRPS from ESO (its
#: transmission in ESO's model); deeper, the pixel is left out
TELLURIC_FLOOR = 0.5
#: the bytes of a tar of ESO read first: the header of its first file
HEAD_BYTES = 512 + 2880 * 24
#: in the near infrared (SPIRou, NIRPS), a star whose spectra span less
#: BERV than this [km/s] is not its own template: the residuals of the
#: telluric lines stay where they are in each spectrum, and its velocities
#: follow its BERV (TOI-1452 on SPIRou, 3 km/s: 119 m/s per km/s, against
#: -1.7 with the template of Gl 699). The template of TEMPLATE_STAR is
#: used instead, TEMPLATE_N of its spectra spread over its seasons
BERV_MIN = 10.0
TEMPLATE_STAR = 'GJ 699'
TEMPLATE_N = 60
#: the header keys of the BERV [km/s], per class of LBL (and the extension)
BERV_KEYS = {'SPIROU': ('FluxAB', 'BERV'), 'NIRPS_HE': (0, 'ESO QC BERV'),
             'NIRPS_HA': (0, 'ESO QC BERV')}
#: the temperatures of LBL's gradients of temperature [K]
GRADIENTS = (3000, 3500, 4000, 4500, 5000, 5500, 6000)
#: the files fetched at once
WORKERS = 4
#: the names of the datasets of the files kept, beside them
LABELS = 'labels.json'
#: how long a spectrum left out (not public yet, a reduction LBL does not
#: read) is not asked for again [days]
LEFT_AGE = 30.0
ALWAYS = (-1e9, 1e9)


# =============================================================================
# Define functions
# =============================================================================
def settings() -> Dict[str, Any]:
    """your choices (SETTINGS): python, data (the folder of LBL's data),
    nproc (its processes)"""
    out: Dict[str, Any] = dict(python='', data='', nproc=0)
    if os.path.exists(SETTINGS):
        try:
            with open(SETTINGS) as handle:
                out.update(json.load(handle))
        except (OSError, ValueError):
            pass
    return out


def save_settings(**kwargs: Any) -> Dict[str, Any]:
    """your choices, kept (those given replace those kept)"""
    out = settings()
    out.update({key: val for key, val in kwargs.items()
                if key in ('python', 'data', 'nproc')})
    os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
    with open(SETTINGS, 'w') as handle:
        json.dump(out, handle, indent=1)
    return out


def find_python(given: Optional[str] = None) -> str:
    """the Python of LBL: the one given, else the one chosen, else that of
    an environment called lbl-rapide where conda keeps them ('' when none)"""
    for one in [given, settings().get('python')] + PYTHONS:
        if one and os.path.isfile(os.path.expanduser(one)):
            return os.path.expanduser(one)
    return ''


def check_python(python: str, timeout: float = 120.0) -> Dict[str, Any]:
    """
    Whether a Python runs LBL: its version of LBL, and where it is

    :return: dict, ok, version, where, error
    """
    if not python or not os.path.isfile(python):
        return dict(ok=False, version='', where='',
                    error='no Python of LBL (an environment with lbl-rapide)')
    code = ('import os, lbl; from lbl.recipes import lbl_wrap; '
            'here = os.path.dirname(lbl.__file__); '
            'print(open(os.path.join(here, "version.txt")).read().split()[0]);'
            ' print(here)')
    try:
        done = subprocess.run([python, '-c', code], capture_output=True,
                              text=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError) as err:
        return dict(ok=False, version='', where='', error=str(err))
    lines = done.stdout.split()
    if done.returncode != 0 or len(lines) < 2:
        why = (done.stderr.strip().splitlines() or ['?'])[-1]
        return dict(ok=False, version='', where='', error=why)
    return dict(ok=True, version=lines[0], where=lines[1], error='')


def object_name(target: str) -> str:
    """the name of the star for LBL (letters and digits: GJ436)"""
    return re.sub(r'[^A-Za-z0-9]+', '', str(target)) or 'star'


def data_dir(lbl_dir: str, lbl: str, source: str) -> str:
    """the folder of a class of LBL (lbl-dir/<INSTRUMENT>_<SOURCE>), the
    models shared"""
    here = os.path.join(os.path.abspath(os.path.expanduser(lbl_dir)),
                        f'{lbl}_{source}')
    models = os.path.join(os.path.dirname(here), 'models')
    os.makedirs(models, exist_ok=True)
    os.makedirs(here, exist_ok=True)
    link = os.path.join(here, 'models')
    if not os.path.lexists(link):
        os.symlink(models, link)
    return here


def _ident(target: str) -> Dict[str, Any]:
    """the star as SIMBAD knows it, with how it moves"""
    from koloa.archive import resolve
    from koloa.secular import star_astrometry
    ident = resolve(target)
    try:
        ident['astrometry'] = star_astrometry(ident)
    except Exception as err:  # its place then is that of 2000
        log(f'astrometry: {type(err).__name__}: {err}', 'warn')
        ident['astrometry'] = None
    return ident


def _star(ident: Dict[str, Any]) -> Dict[str, Any]:
    """the star as the archives are asked about it"""
    astro = ident.get('astrometry') or {}
    return dict(ra=ident['ra'], dec=ident['dec'], pmra=astro.get('pmra'),
                pmdec=astro.get('pmdec'),
                names=[ident.get('main') or ''] + list(
                    ident.get('aliases') or []))


# -----------------------------------------------------------------------------
# The spectra of each instrument, and how a file of them is fetched
# -----------------------------------------------------------------------------
def harpsn_spectra(target: str, root: str = 'archives',
                   ident: Optional[Dict[str, Any]] = None) -> List[str]:
    """
    The espdr S2D_A files of HARPS-N of a star: those of the exposures
    koloa keeps from their CCFs (koloa.tng.series), read from the folder of
    the star when its archives were gathered, else from what koloa's users
    shared or from the archive

    :return: list of str, the names of the files (gzipped, as the archive
             has them), by date
    """
    from koloa import tng
    from koloa.gather import folder_name
    folder = os.path.join(root, folder_name(target), 'rv', 'published')
    table = os.path.join(folder, 'tng_exposures.csv')
    if os.path.exists(table):
        import csv
        with open(table) as handle:
            rows = list(csv.DictReader(handle))
    else:
        os.makedirs(folder, exist_ok=True)
        rows = tng.exposures(_star(ident or _ident(target)), folder)
    kept = [about for _, about in tng.series(rows, target)
            if about['reduction'] == 'espdr']
    if not kept:
        return []
    return [name.replace('_CCF_A.fits', '_S2D_A.fits')
            for name in kept[0]['files']]


def _download(url: str, path: str, timeout: float = 300.0,
              gunzip: bool = False) -> str:
    """a file of an archive, gunzipped if asked, written whole or not at
    all; '' or why not"""
    part = path + '.part'
    try:
        with urllib.request.urlopen(urllib.request.Request(
                url, headers={'User-Agent': 'koloa'}),
                timeout=timeout) as response, open(part, 'wb') as handle:
            source = gzip.GzipFile(fileobj=response) if gunzip else response
            shutil.copyfileobj(source, handle, 1 << 20)
        os.replace(part, path)
        return ''
    except (urllib.error.URLError, TimeoutError, OSError, EOFError) as err:
        if os.path.exists(part):
            os.remove(part)
        return f'{type(err).__name__}: {err}'


def _range(url: str, start: int, size: int, timeout: float = 300.0
           ) -> bytes:
    """bytes of a file of an archive (an HTTP Range request)"""
    with urllib.request.urlopen(urllib.request.Request(
            url, headers={'User-Agent': 'koloa', 'Range':
                          f'bytes={start}-{start + size - 1}'}),
            timeout=timeout) as response:
        data = response.read()
    if len(data) != size:
        raise OSError(f'{len(data)} bytes of {size} ({url})')
    return data


def tar_members(url: str) -> List[Tuple[str, int, int]]:
    """
    The files of a tar of an archive, read from their headers alone (a
    Range request of 512 bytes each): name, where its data start, size

    :return: list of (name, offset, size)
    """
    out, where = [], 0
    for _ in range(1000):
        head = _range(url, where, 512)
        if not head.strip(b'\0'):
            break
        size = int(head[124:136].rstrip(b'\0 ').decode() or '0', 8)
        name = head[:100].rstrip(b'\0').decode('utf8', 'replace')
        kind = head[156:157]
        if kind in (b'L', b'x'):
            # a long name (GNU) or the attributes of the next (pax)
            extra = _range(url, where + 512, size).rstrip(b'\0').decode(
                'utf8', 'replace')
            where += 512 + -(-size // 512) * 512
            head = _range(url, where, 512)
            size = int(head[124:136].rstrip(b'\0 ').decode() or '0', 8)
            name = head[:100].rstrip(b'\0').decode('utf8', 'replace')
            if kind == b'L':
                name = extra
            for line in extra.splitlines():
                key, _, val = line.partition(' ')[2].partition('=')
                if key == 'path':
                    name = val
                elif key == 'size':
                    size = int(val)
        else:
            prefix = head[345:500].rstrip(b'\0').decode('utf8', 'replace')
            name = f'{prefix}/{name}' if prefix else name
        out.append((name.lstrip('./'), where + 512, size))
        where += 512 + -(-size // 512) * 512
    return out


def _get_tng(one: Dict[str, Any], science: str, calib: str) -> str:
    """an S2D_A of the TNG, gunzipped"""
    return _download(TNG_FILES + one['id'], os.path.join(science, one['name']),
                     gunzip=True)


def _get_sophie(one: Dict[str, Any], science: str, calib: str) -> str:
    """an e2ds of SOPHIE: its spectrum (fibre A) as the science file, its
    blaze as the calibration file its header names; a spectrum of a
    programme still protected (no BERV, its date rounded) left out"""
    try:
        from astropy.io import fits
    except ImportError:
        return 'astropy is needed to read the spectra of SOPHIE'
    try:
        with urllib.request.urlopen(urllib.request.Request(
                SOPHIE_E2DS.format(seq=one['id']),
                headers={'User-Agent': 'koloa'}), timeout=300) as response:
            raw = response.read()
    except (urllib.error.URLError, TimeoutError, OSError) as err:
        return f'{type(err).__name__}: {err}'
    if not raw:
        return 'not public yet'
    try:
        with fits.open(io.BytesIO(raw)) as hdul:
            head = hdul[0].header.copy()
            bjd = head.get('HIERARCH OHP DRS BJD')
            if head.get('HIERARCH OHP DRS BERV') is None or bjd is None \
                    or float(bjd) == int(float(bjd)):
                return 'of a programme still protected (no BERV)'
            blaze = f'{one["id"]}_blaze_A.fits'
            head['HIERARCH OHP DRS BLAZE FILE'] = blaze
            fits.PrimaryHDU(hdul['BLAZE_A'].data).writeto(
                os.path.join(calib, blaze), overwrite=True)
            path = os.path.join(science, one['name'])
            fits.PrimaryHDU(hdul[0].data, head).writeto(path + '.part',
                                                        overwrite=True)
            os.replace(path + '.part', path)
    except (OSError, KeyError, ValueError) as err:
        return f'{type(err).__name__}: {err}'
    return ''


def _designation(name: str) -> Optional[Tuple[str, str]]:
    """a name as a catalogue and a number (GJ 699, Gl699, GJ699-HE: GJ and
    699; TOI-1452: TOI and 1452; GJ 15 B: GJ and 15B), None when it is not
    one (K4_02, BARNARD)"""
    found = re.match(r'^\s*([A-Za-z]+)[\s_-]*(\d+)(?:[\s_.-]*([A-Ca-c])\b)?',
                     str(name))
    if not found:
        return None
    family = found[1].upper()
    family = {'GL': 'GJ', 'GLIESE': 'GJ'}.get(family, family)
    return family, found[2].lstrip('0') + (found[3] or '').upper()


def _other(name: str, mine: Dict[str, set]) -> bool:
    """whether the name of a target is that of another star: of a
    catalogue the star is in, under another number (TOI-1760 beside
    TOI-1452, GJ 15 B beside GJ 15 A); a name of no catalogue of the star,
    or none, is not (an observer's own: K4_02)"""
    found = _designation(name)
    if found is None or found[0] not in mine:
        return False
    number = found[1]
    for one in mine[found[0]]:
        # a component given on one side only is the same star
        if number == one or number.rstrip('ABC') == one or \
                one.rstrip('ABC') == number:
            return False
    return True


def _named(found: List[Dict[str, Any]], star: Dict[str, Any],
           what: str) -> List[Dict[str, Any]]:
    """
    The exposures pointed at a star that are not of another star by their
    name: one of a catalogue the star is in, under another number (the
    companion of a binary, TOI-1760 beside TOI-1452), would spoil the
    template of LBL; the observers' own names are kept (GJ699-HE, Gl699)

    :param found: the exposures, each with its target
    """
    mine: Dict[str, set] = {}
    for name in star.get('names') or []:
        got = _designation(name)
        if got:
            mine.setdefault(got[0], set()).add(got[1])
    kept = [one for one in found if not _other(one['target'], mine)]
    gone = sorted({one['target'] for one in found} -
                  {one['target'] for one in kept})
    if gone:
        log(f'{what}: {len(found) - len(kept)} exposures of another star '
            f'pointed nearby left out ({", ".join(gone)})', 'warn')
    return kept


def spirou_spectra(star: Dict[str, Any]) -> List[str]:
    """
    The t.fits of SPIRou of a star at the CADC (APERO's spectra, the
    telluric lines taken out): those public, pointed at where the star was
    (koloa.observatories, as its CCFs are)

    :return: list of str, their URIs (cadc:CFHT/<odometer>t.fits), by date
    """
    from koloa import observatories as obs
    rows = obs._tap(obs.CADC_TAP, (
        "SELECT o.target_name AS target, "
        "o.targetPosition_coordinates_cval1 AS ra, "
        "o.targetPosition_coordinates_cval2 AS dec, p.time_bounds_lower AS "
        "mjd, p.dataRelease AS release, a.uri AS uri FROM caom2.Observation "
        "AS o JOIN caom2.Plane AS p ON o.obsID = p.obsID JOIN "
        "caom2.Artifact AS a ON a.planeID = p.planeID WHERE o.collection = "
        "'CFHT' AND o.instrument_name = 'SPIRou' AND p.productID LIKE '%t' "
        "AND a.contentType = 'application/fits' AND "
        f"INTERSECTS(CIRCLE('ICRS', {float(star['ra']):.6f}, "
        f"{float(star['dec']):.6f}, {obs._reach(star):.5f}), "
        "p.position_bounds) = 1"), 'CFHT archive (CADC)')
    now = _time.strftime('%Y-%m-%dT%H:%M:%S', _time.gmtime())
    found = [dict(mjd=obs._num(row['mjd']), uri=str(row['uri']),
                  target=str(row.get('target') or '')) for row in rows
             if str(row.get('release') or '9') <= now
             and str(row['uri']).endswith('t.fits')
             and obs._pointed(obs._num(row['ra']), obs._num(row['dec']),
                              obs._num(row['mjd']), star)[0]]
    kept = _named(found, star, 'SPIRou')
    return list(dict.fromkeys(one['uri'] for one in sorted(
        kept, key=lambda one: one['mjd'])))


def _get_cadc(one: Dict[str, Any], science: str, calib: str) -> str:
    """a t.fits of SPIRou at the CADC (a signed address it redirects to);
    one of an older reduction, whose FluxAB LBL cannot read (no MJDMID),
    left out"""
    from koloa.observatories import CADC_FILE
    path = os.path.join(science, one['name'])
    why = _download(CADC_FILE + one['id'], path)
    if why:
        return why
    try:
        from astropy.io import fits
        with fits.open(path) as hdul:
            head = hdul['FluxAB'].header
            if 'MJDMID' not in head:
                version = hdul[0].header.get('VERSION', '?')
                os.remove(path)
                return f'APERO {version}: an older reduction LBL does not read'
    except ImportError:
        pass
    except (OSError, KeyError) as err:
        os.remove(path)
        return f'{type(err).__name__}: {err}'
    return ''


def _safe(name: str) -> str:
    """a name of a file with only plain characters (LBL asks for them)"""
    return re.sub(r'[^A-Za-z0-9._-]+', '-', name)


def eso_spectra(star: Dict[str, Any], ins: str) -> List[Dict[str, Any]]:
    """
    The exposures of an instrument of ESO (NIRPS, ESPRESSO) of a star with
    their 2D spectra (the ancillary tar of their product), public, pointed
    at where the star was; with what the header of the first spectrum of
    each tar says (its mode, its binning), one Range request a tar

    :return: list of dict: id (the tar), target, mjd, mode, bin, object, by
             date
    """
    from koloa import observatories as obs
    from koloa.tng import header
    ra, dec = float(star['ra']), float(star['dec'])
    rows = obs._tap(obs.ESO_OBS, (
        "SELECT o.target_name AS target, o.t_min AS mjd, o.obs_release_date "
        "AS release, o.s_ra AS ra, o.s_dec AS dec, f.archive_id AS tar FROM "
        "ivoa.ObsCore AS o JOIN phase3v2.product_files AS f ON f.product_id"
        " = o.dp_id WHERE CONTAINS(POINT('ICRS', o.s_ra, o.s_dec), "
        f"CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {obs._reach(star):.5f}))=1 AND "
        f"o.instrument_name = '{ins}' AND f.eso_category = "
        "'ANCILLARY.2DECHELLE.TAR'"), 'ESO archive', 'json')
    now = _time.strftime('%Y-%m-%dT%H:%M:%S', _time.gmtime())
    found = {}
    for row in rows:
        if str(row.get('release') or '9') > now or not obs._pointed(
                obs._num(row['ra']), obs._num(row['dec']),
                obs._num(row['mjd']), star)[0]:
            continue
        found[str(row['tar'])] = dict(id=str(row['tar']),
                                      target=str(row['target']),
                                      mjd=obs._num(row['mjd']))

    def head(one):
        try:
            raw = _range(obs.ESO_FILE + one['id'], 0, HEAD_BYTES, 120.0)
        except (urllib.error.URLError, TimeoutError, OSError):
            return None
        cards = header(raw[512:])
        return dict(one, mode=str(cards.get('ESO INS MODE') or '').strip(),
                    bin=f'{cards.get("ESO DET BINX", 1)}'
                        f'{cards.get("ESO DET BINY", 1)}',
                    object=str(cards.get('OBJECT') or ''))
    with ThreadPoolExecutor(max_workers=8) as pool:
        out = [one for one in pool.map(head, _named(list(found.values()),
                                                     star, ins)) if one]
    return sorted(out, key=lambda one: one['mjd'])


def _member(url: str, members: List[Tuple[str, int, int]], *ends: str
            ) -> Optional[Tuple[str, int, int]]:
    """the file of a tar whose name ends one of the ways given"""
    return next((one for end in ends for one in members
                 if one[0].endswith(end)), None)


def _get_nirps(one: Dict[str, Any], science: str, calib: str) -> str:
    """
    A spectrum of NIRPS from ESO: its S2D with the blaze (S2D_BLAZE_A)
    divided by ESO's model of the telluric lines (S2D_TELL_SPECTRUM_A),
    the lines deeper than TELLURIC_FLOOR left out, as the science file;
    the blaze of its night, as the calibration file its header names
    """
    try:
        from astropy.io import fits
        import numpy as np
    except ImportError:
        return 'astropy is needed to read the spectra of ESO'
    from koloa.observatories import ESO_FILE
    url = ESO_FILE + one['id']
    try:
        members = tar_members(url)
        flux = _member(url, members, '_S2D_BLAZE_A.fits')
        tell = _member(url, members, '_S2D_TELL_SPECTRUM_A.fits')
        # the blaze of the night (not the S2D with its blaze)
        blaze = next((part for part in members
                      if part[0].endswith('_BLAZE_A.fits')
                      and '_S2D_' not in part[0]), None)
        if not flux or not tell or not blaze:
            return 'no S2D_BLAZE_A, telluric model or blaze in its tar'
        bname = os.path.basename(blaze[0])
        if not os.path.exists(os.path.join(calib, bname)):
            # the blaze of a night, which several spectra share: each its
            #   own part, the last one written wins
            part = os.path.join(calib, f'{bname}.{uuid.uuid4().hex}.part')
            with open(part, 'wb') as handle:
                handle.write(_range(url, blaze[1], blaze[2]))
            os.replace(part, os.path.join(calib, bname))
        with fits.open(io.BytesIO(_range(url, flux[1], flux[2]))) as hdul, \
                fits.open(io.BytesIO(_range(url, tell[1], tell[2]))) as model:
            trans = np.asarray(model['SCIDATA'].data, float)
            bad = ~(trans > TELLURIC_FLOOR)
            with np.errstate(divide='ignore', invalid='ignore'):
                for ext in ('SCIDATA', 'ERRDATA'):
                    data = np.asarray(hdul[ext].data, float) / trans
                    data[bad] = np.nan
                    hdul[ext].data = data.astype(np.float32)
            hdul[0].header['HIERARCH KOLOA TELL FLOOR'] = (
                TELLURIC_FLOOR, 'koloa: divided by S2D_TELL_SPECTRUM_A')
            path = os.path.join(science, one['name'])
            hdul.writeto(path + '.part', overwrite=True)
            os.replace(path + '.part', path)
    except (urllib.error.URLError, TimeoutError, OSError, KeyError,
            ValueError) as err:
        return f'{type(err).__name__}: {err}'
    return ''


def _get_espresso(one: Dict[str, Any], science: str, calib: str) -> str:
    """a spectrum of ESPRESSO from ESO: the S2D_A of its tar (blaze
    corrected), as it is"""
    from koloa.observatories import ESO_FILE
    url = ESO_FILE + one['id']
    try:
        members = tar_members(url)
        flux = _member(url, members, '_S2D_A.fits') or next(
            (part for part in members
             if os.path.basename(part[0]).startswith('ES_S2DA')), None)
        if not flux:
            return 'no S2D_A in its tar'
        path = os.path.join(science, one['name'])
        with open(path + '.part', 'wb') as handle:
            handle.write(_range(url, flux[1], flux[2], 900.0))
        os.replace(path + '.part', path)
    except (urllib.error.URLError, TimeoutError, OSError) as err:
        return f'{type(err).__name__}: {err}'
    return ''


def groups(target: str, inst: str, root: str = 'archives',
           ident: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """
    The runs of LBL of an instrument for a star: for each, its spectra
    (id in the archive, name here), the class of LBL that reads them, its
    parameters, and the datasets its rdb makes (eras: suffix of the file,
    name of the dataset, from and to [BJD - 2400000])

    :return: list of dict: key, obj, lbl, source, input, tellu, params,
             spectra, eras, get
    """
    obj = object_name(target)
    if inst == 'HARPSN':
        names = harpsn_spectra(target, root, ident)
        return [dict(key='', obj=obj, lbl='HARPSN', source='ESO',
                     input='*S2D_A.fits', tellu=True, params={},
                     spectra=[dict(id=name, name=name[:-3]) for name in names],
                     eras=[('', 'HARPS-N (LBL)') + ALWAYS], get=_get_tng, via='tng')]
    if inst == 'SOPHIE':
        from koloa import observatories
        parts, _ = observatories.sophie(_star(ident or _ident(target)), '',
                                        name=target)
        out = []
        for mode, key in (('HR', ''), ('HE', 'HE')):
            seqs = [seq for _, about in parts
                    if (' HE ' in about['label']) == (mode == 'HE')
                    for seq in about.get('files', [])]
            if not seqs:
                continue
            tail = ' HE' if mode == 'HE' else ''
            cut = observatories.SOPHIE_2011
            out.append(dict(
                key=key, obj=obj + key, lbl='SOPHIE', source='None',
                input='*_e2ds_A.fits', tellu=True, params={},
                spectra=[dict(id=seq, name=f'{seq}_e2ds_A.fits')
                         for seq in seqs],
                eras=[(key, f'SOPHIE{tail} (LBL)', ALWAYS[0], cut),
                      (key + 'p', f'SOPHIE+{tail} (LBL)', cut, ALWAYS[1])],
                get=_get_sophie, via='sophie'))
        return out
    if inst == 'SPIROU':
        uris = spirou_spectra(_star(ident or _ident(target)))
        return [dict(key='', obj=obj, lbl='SPIROU', source='CADC',
                     input='*t.fits', tellu=False, params={},
                     spectra=[dict(id=uri, name=uri.rsplit('/', 1)[-1])
                              for uri in uris],
                     eras=[('', 'SPIRou (LBL)') + ALWAYS], get=_get_cadc, via='cadc')]
    if inst == 'NIRPS':
        found = eso_spectra(_star(ident or _ident(target)), 'NIRPS')
        out = []
        for mode, key in (('HE', ''), ('HA', 'HA')):
            mine = [one for one in found if one['mode'] == mode]
            if mine:
                out.append(dict(
                    key=key, obj=obj + key, lbl=f'NIRPS_{mode}', source='ESO',
                    input='*TCORR_A.fits', tellu=False,
                    params=dict(BLAZE_CORRECTED=False),
                    spectra=[dict(id=one['id'], name=_safe(
                        f'r.NIRPS.{one["id"]}_S2D_BLAZE_TCORR_A.fits'))
                        for one in mine],
                    eras=[(key, 'NIRPS' + ('' if mode == 'HE' else '_HA')
                           + ' (LBL)') + ALWAYS], get=_get_nirps, via='nirps'))
        return out
    if inst == 'ESPRESSO':
        from koloa.observatories import ESPRESSO_2019
        found = eso_spectra(_star(ident or _ident(target)), 'ESPRESSO')
        kinds = [one['mode'].replace('SINGLE', '').replace('MULTI', '')
                 + one['bin'] for one in found]
        out = []
        if kinds:
            main = max(set(kinds), key=kinds.count)
            for kind in sorted(set(kinds), key=lambda one: (one != main, one)):
                mine = [one for one, it in zip(found, kinds) if it == kind]
                key = '' if kind == main else kind
                tail = '' if kind == main else f'_{kind}'
                out.append(dict(
                    key=key, obj=obj + key, lbl='ESPRESSO', source='None',
                    input='*S2D_A.fits', tellu=True, params={},
                    spectra=[dict(id=one['id'], name=_safe(
                        f'ESPRESSO.{one["id"]}_S2D_A.fits'))
                        for one in mine],
                    eras=[(key + '18', f'ESPRESSO18{tail} (LBL)', ALWAYS[0],
                           ESPRESSO_2019),
                          (key + '19', f'ESPRESSO19{tail} (LBL)',
                           ESPRESSO_2019, ALWAYS[1])],
                    get=_get_espresso, via='espresso'))
        return out
    raise ValueError(f'no instrument {inst} for LBL: {", ".join(INSTRUMENTS)}')


#: how the spectra of a run are fetched, by its name in a plan
GETTERS = dict(tng=_get_tng, sophie=_get_sophie, cadc=_get_cadc,
               nirps=_get_nirps, espresso=_get_espresso)


def plain_group(group: Dict[str, Any]) -> Dict[str, Any]:
    """a run as a plan keeps it (JSON: the way its spectra are fetched by
    its name)"""
    return {key: val for key, val in group.items() if key != 'get'}


def live_group(group: Dict[str, Any]) -> Dict[str, Any]:
    """a run of a plan, its way of fetching its spectra back"""
    return dict(group, get=GETTERS[group['via']],
                eras=[tuple(era) for era in group['eras']])


def folders(group: Dict[str, Any], lbl_dir: str) -> Tuple[str, str, str]:
    """the folders of a run: that of its class of LBL, its spectra, its
    calibrations"""
    here = data_dir(lbl_dir, group['lbl'], group['source'])
    science = os.path.join(here, 'science', group['obj'])
    calib = os.path.join(here, 'calib')
    os.makedirs(science, exist_ok=True)
    os.makedirs(calib, exist_ok=True)
    return here, science, calib


def _transient(why: str) -> bool:
    """whether a spectrum not fetched may be fetched when asked again (the
    archive silent, a file not written), not left out for good (not public
    yet, of a protected programme, a reduction LBL does not read)"""
    return bool(re.match(r'[A-Za-z]*(Error|Exception|Timeout)\b', why)
                or 'timed out' in why or 'bytes of' in why or 'astropy' in why)


def fetch(target: str, inst: str, group: Dict[str, Any], lbl_dir: str
          ) -> Dict[str, Any]:
    """
    The spectra of a run of LBL in its science folder
    (lbl-dir/<INSTRUMENT>_<SOURCE>/science/<star>/): only those not there
    yet are fetched, a few at a time

    :return: dict, folder, n (the spectra there), fetched, failed (and why,
             the first)
    """
    here, science, calib = folders(group, lbl_dir)
    # the spectra left out before, and why (asked again after LEFT_AGE)
    book = os.path.join(here, f'koloa_left_{group["obj"]}.json')
    left: Dict[str, Any] = {}
    if os.path.exists(book):
        with open(book) as handle:
            left = {key: val for key, val in json.load(handle).items()
                    if _time.time() - val[0] < LEFT_AGE * 86400}
    # a spectrum there that is no longer one of the star's (another
    #   target pointed nearby, a reduction left out): moved aside, as LBL
    #   reads every file of the folder
    wanted = {one['name'] for one in group['spectra']}
    aside = [name for name in os.listdir(science)
             if name.endswith('.fits') and name not in wanted]
    if aside:
        moved = os.path.join(here, 'left', group['obj'])
        os.makedirs(moved, exist_ok=True)
        for name in aside:
            os.replace(os.path.join(science, name),
                       os.path.join(moved, name))
        log(f'{len(aside)} spectra no longer of the star moved to {moved}',
            'warn')
    todo = [one for one in group['spectra']
            if not os.path.exists(os.path.join(science, one['name']))
            and one['id'] not in left]
    label = INSTRUMENTS[inst]['label']
    log(f'{label}: {len(group["spectra"])} spectra of {target} in '
        f'{INSTRUMENTS[inst]["archive"]}, {len(todo)} to fetch (about '
        f'{len(todo) * INSTRUMENTS[inst]["size"]:.0f} MB)', 'info')
    errors: List[str] = []
    if todo:
        from koloa.tng import _steps
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            whys = list(_steps(pool.map(lambda one: group['get'](
                one, science, calib), todo), len(todo),
                f'{label} spectra', f'lbl_{inst}'))
        # those that failed (the archive silent, a file not written):
        #   asked once more, one at a time
        for one, why in zip(todo, whys):
            if why and _transient(why):
                why = group['get'](one, science, calib)
            if why:
                errors.append(why)
                # not public, or not read by LBL: not asked again soon
                if not _transient(why):
                    left[one['id']] = (_time.time(), why)
    with open(book, 'w') as handle:
        json.dump(left, handle)
    have = len([name for name in os.listdir(science)
                if name.endswith('.fits')])
    if errors:
        common = max(set(errors), key=errors.count)
        log(f'{len(errors)} spectra not fetched ({common})', 'warn')
    return dict(folder=science, n=have, fetched=len(todo) - len(errors),
                failed=len(errors),
                why=max(set(errors), key=errors.count) if errors else '')


# -----------------------------------------------------------------------------
# LBL
# -----------------------------------------------------------------------------
def parameters(target: str, group: Dict[str, Any], lbl_dir: str,
               teff: float, nproc: int = 4, tellu: Optional[bool] = None
               ) -> Dict[str, Any]:
    """the parameters of LBL for a run (as in a wrap file): its template
    and masks made from the star's own spectra, the telluric lines cleaned
    first (optical spectra), the gradient of temperature of LBL nearest
    its temperature projected out (DTEMP)"""
    near = min(GRADIENTS, key=lambda one: abs(one - float(teff)))
    tellu = group['tellu'] if tellu is None else bool(tellu) and group[
        'tellu']
    return dict(
        INSTRUMENT=group['lbl'], DATA_SOURCE=group['source'],
        DATA_DIR=folders(group, lbl_dir)[0], INPUT_FILE=group['input'],
        DATA_TYPES=['SCIENCE'], OBJECT_SCIENCE=[group['obj']],
        OBJECT_COMPARISON=[group['obj']],
        OBJECT_TEFF=[int(round(min(max(float(teff), 3000.0), 6000.0)))],
        RUN_LBL_RESET=False, RUN_LBL_TELLUCLEAN=tellu,
        RUN_LBL_TEMPLATE=True, RUN_LBL_MASK=True, RUN_LBL_COMPUTE=True,
        RUN_LBL_COMPILE=True, SKIP_LBL_TEMPLATE=True, SKIP_LBL_MASK=True,
        SKIP_LBL_COMPUTE=True, SKIP_LBL_COMPILE=False,
        RESPROJ_TABLES={f'DTEMP{near}': f'temperature_gradient_{near}.fits'},
        COMPUTE_NPROC=max(1, int(nproc)), PLOT=False, **group['params'])


def wrap(params: Dict[str, Any]) -> str:
    """the script that runs LBL with its parameters (a wrap file)"""
    return ('#!/usr/bin/env python\n'
            '"""LBL, run by koloa (koloa.lblrun): its parameters as koloa '
            'wrote them"""\n'
            'import json\n\n'
            'from lbl import lbl_wrap\n\n'
            f'rparams = json.loads({json.dumps(json.dumps(params, indent=1))})'
            '\n\n'
            'if __name__ == "__main__":\n'
            '    lbl_wrap(rparams)\n')


def rdb_of(params: Dict[str, Any]) -> str:
    """where LBL writes the rdb of a run (with the telluric cleaning, the
    names end in _tc; the template's is its star's, another one's when
    the template is of another star)"""
    tail = '_tc' if params.get('RUN_LBL_TELLUCLEAN') else ''
    obj = params['OBJECT_SCIENCE'][0] + tail
    tpl = params['OBJECT_COMPARISON'][0] + tail
    return os.path.join(params['DATA_DIR'], 'lblrdb', f'lbl_{obj}_{tpl}.rdb')


def berv_span(group: Dict[str, Any], lbl_dir: str) -> Optional[float]:
    """the span of the BERV of the spectra of a run of the near infrared
    [km/s] (None for an optical one, or when not read)"""
    if group['lbl'] not in BERV_KEYS:
        return None
    try:
        from astropy.io import fits
    except ImportError:
        return None
    ext, key = BERV_KEYS[group['lbl']]
    science = folders(group, lbl_dir)[1]
    values = []
    for name in os.listdir(science):
        if not name.endswith('.fits'):
            continue
        try:
            value = fits.getheader(os.path.join(science, name), ext).get(
                f'HIERARCH {key}' if ' ' in key else key)
        except (OSError, KeyError):
            continue
        if value is not None:
            values.append(float(value))
    return float(max(values) - min(values)) if values else None


def template_group(group: Dict[str, Any], inst: str, star: str,
                   root: str = 'archives') -> Optional[Dict[str, Any]]:
    """the run of another star whose template a star with too little BERV
    uses: its spectra of the same setting, TEMPLATE_N spread over its
    seasons (None when it has none)"""
    for one in groups(star, inst, root):
        if one['lbl'] == group['lbl'] and one['spectra']:
            every = max(1, len(one['spectra']) // TEMPLATE_N)
            return dict(one, spectra=one['spectra'][::every])
    return None


def run(params: Dict[str, Any], python: str,
        echo: Optional[Callable[[str], None]] = None) -> str:
    """
    LBL, with its Python, in a process of its own: its lines said as they
    come (echo, else printed), the script kept beside its data

    :return: str, the rdb it made
    """
    script = os.path.join(params['DATA_DIR'],
                          f'koloa_wrap_{params["OBJECT_SCIENCE"][0]}.py')
    with open(script, 'w') as handle:
        handle.write(wrap(params))
    echo = echo or (lambda line: print(line, flush=True))
    colour = re.compile(r'\x1b\[[0-9;]*m')
    last: List[str] = []
    env = dict(os.environ, PYTHONUNBUFFERED='1')
    with subprocess.Popen([python, script], cwd=params['DATA_DIR'],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, bufsize=1, env=env) as proc:
        for line in proc.stdout:
            # its words, without the colours of a terminal
            line = colour.sub('', line).rstrip()
            if line:
                echo(line)
                last = (last + [line])[-5:]
        code = proc.wait()
    if code != 0:
        raise RuntimeError(f'LBL stopped (exit {code}): '
                           + ' / '.join(last[-2:]))
    rdb = rdb_of(params)
    if not os.path.exists(rdb):
        raise RuntimeError(f'LBL made no rdb ({rdb})')
    return rdb


def keep(rdb: str, target: str, inst: str, group: Dict[str, Any],
         root: str = 'archives') -> List[Tuple[str, str, int]]:
    """
    The rdb of LBL in the folder of the star, under its name, a file per
    dataset (era) with velocities (archives/GJ_436/lbl/
    GJ_436_SOPHIEp_lbl.rdb), their names in labels.json beside them

    :return: list of (path, the name of its dataset, its velocities)
    """
    from koloa.gather import folder_name
    name = folder_name(target)
    here = os.path.join(root, name, 'lbl')
    os.makedirs(here, exist_ok=True)
    with open(rdb) as handle:
        lines = handle.read().splitlines()
    head, dashes, rows = lines[0], lines[1], [line for line in lines[2:]
                                              if line.strip()]
    col = head.split('\t').index('rjd')
    labels = {}
    if os.path.exists(os.path.join(here, LABELS)):
        with open(os.path.join(here, LABELS)) as handle:
            labels = json.load(handle)
    out = []
    for suffix, label, low, high in group['eras']:
        mine = [row for row in rows
                if low <= float(row.split('\t')[col]) < high]
        path = os.path.join(here, f'{name}_{inst}{suffix}_lbl.rdb')
        if not mine:
            if os.path.exists(path):
                os.remove(path)
            labels.pop(os.path.basename(path), None)
            continue
        with open(path, 'w') as handle:
            handle.write('\n'.join([head, dashes] + mine) + '\n')
        labels[os.path.basename(path)] = label
        out.append((path, label, len(mine)))
    with open(os.path.join(here, LABELS), 'w') as handle:
        json.dump(labels, handle, indent=1)
    return out


def kept(target: str, root: str = 'archives') -> List[Dict[str, Any]]:
    """the rdb of LBL kept for a star: path, inst, label (its dataset),
    npoints, made"""
    from koloa.gather import folder_name
    name = folder_name(target)
    here = os.path.join(root, name, 'lbl')
    if not os.path.isdir(here):
        return []
    labels = {}
    if os.path.exists(os.path.join(here, LABELS)):
        with open(os.path.join(here, LABELS)) as handle:
            labels = json.load(handle)
    out = []
    for file in sorted(os.listdir(here)):
        found = re.match(re.escape(name) + r'_([A-Z]+?)([A-Za-z0-9]*)_lbl\.rdb$',
                         file)
        if not found:
            continue
        inst = next((key for key in INSTRUMENTS
                     if (found[1] + found[2]).startswith(key)), found[1])
        path = os.path.join(here, file)
        with open(path) as handle:
            npoints = max(0, sum(1 for _ in handle) - 2)
        out.append(dict(path=os.path.abspath(path), inst=inst,
                        label=labels.get(file, INSTRUMENTS.get(
                            inst, {}).get('label', inst) + ' (LBL)'),
                        npoints=npoints, made=_time.strftime(
                            '%Y-%m-%d %H:%M', _time.localtime(
                                os.path.getmtime(path)))))
    return out


def run_group(target: str, inst: str, group: Dict[str, Any], lbl_dir: str,
              teff: float, python: str, nproc: int = 4,
              tellu: Optional[bool] = None, template: Optional[str] = None,
              root: str = 'archives', ident: Optional[Dict[str, Any]] = None,
              push: Optional[bool] = None) -> List[Tuple[str, str, int]]:
    """
    One run of LBL, its spectra there: its parameters, the template of
    another star in the near infrared when the BERV of its spectra spans
    less than BERV_MIN (or the one named), LBL, its rdb kept in the folder
    of the star and added to the team's repository (koloa.rvstore: when you
    keep yours there, or push)

    :return: list of (path, the name of its dataset, its velocities)
    """
    name = INSTRUMENTS[inst]['label'] + (f' {group["key"]}' if group['key']
                                         else '')
    params = parameters(target, group, lbl_dir, teff, nproc, tellu)
    span = berv_span(group, lbl_dir)
    other = template
    if other is None and span is not None and span < BERV_MIN and \
            object_name(TEMPLATE_STAR) != group['obj']:
        other = TEMPLATE_STAR
        log(f'{name}: the BERV of its spectra spans {span:.1f} km/s, less '
            f'than {BERV_MIN:g}: the template of {other}', 'warn')
    if other:
        step(f'{name}: the template of {other}')
        mine = template_group(group, inst, other, root)
        if mine is None:
            raise RuntimeError(f'{other} has no spectra of {name} for a '
                               f'template')
        got = fetch(other, inst, mine, lbl_dir)
        outcome(f'{got["n"]} spectra of {other}')
        params['OBJECT_COMPARISON'] = [mine['obj']]
    step(f'{name}: LBL')
    start = _time.time()
    rdb = run(params, python)
    files = keep(rdb, target, inst, group, root)
    try:
        from koloa import rvstore
        from koloa.gather import folder_name
        rvstore.push(ident or dict(main=target), os.path.join(
            root, folder_name(target)), [target], always=push)
    except Exception as err:  # a help, not a need
        log(f'koloa-rvs: {type(err).__name__}: {err}', 'warn')
    outcome('; '.join(f'{label} {n}: {path}' for path, label, n in files)
            + f' ({_time.time() - start:.0f} s)')
    return files


def lbl_star(target: str, instruments: Sequence[str], lbl_dir: str,
             root: str = 'archives', python: Optional[str] = None,
             nproc: Optional[int] = None, teff: Optional[float] = None,
             tellu: Optional[bool] = None, get: bool = True,
             compute: bool = True,
             template: Optional[str] = None) -> Dict[str, Any]:
    """
    The LBL velocities of a star from the public spectra of the archives:
    for each instrument, each of its runs, its spectra fetched (get), LBL
    run on them (compute), its rdb kept in the folder of the star

    :param target: str, a name SIMBAD knows
    :param instruments: the instruments (keys of INSTRUMENTS)
    :param lbl_dir: str, the folder of LBL's data (spectra, templates...)
    :param root: str, the folder of the folders of the stars
    :param python: str or None, the Python of LBL (find_python)
    :param nproc: int or None, the processes of LBL's compute (half the
                  cores of the machine when not said)
    :param teff: float or None, the temperature of the star [K]
                 (koloa.stars.stellar_teff when not said)
    :param tellu: bool or None, clean the telluric lines (the instrument's
                  default when None; never where LBL does not)
    :param template: str or None, the star whose template LBL uses (its
                     spectra fetched too); None: the star's own, but in the
                     near infrared when its BERV spans less than BERV_MIN,
                     that of TEMPLATE_STAR

    :return: dict, per instrument, per run: n (spectra), files (those
             kept), or error
    """
    unknown = [inst for inst in instruments if inst not in INSTRUMENTS]
    if unknown:
        raise ValueError(f'no instrument {unknown} for LBL: '
                         f'{", ".join(INSTRUMENTS)}')
    python = find_python(python)
    if compute:
        state = check_python(python)
        if not state['ok']:
            raise RuntimeError(f'LBL cannot be run: {state["error"]}')
        log(f'LBL {state["version"]} ({state["where"]}), with {python}',
            'info')
    ident = _ident(target)
    if teff is None and compute:
        from koloa.stars import stellar_teff
        guess = stellar_teff(ident)
        teff = guess.get('teff') or 4000.0
        log(f'{target}: Teff {teff:.0f} K ({guess.get("source") or "a guess"})',
            'value')
    nproc = int(nproc or settings().get('nproc') or
                max(1, (os.cpu_count() or 2) // 2))
    out: Dict[str, Any] = {}
    for inst in instruments:
        label = INSTRUMENTS[inst]['label']
        out[inst] = {}
        try:
            step(f'{label}: what the archive has')
            runs = groups(target, inst, root, ident)
            outcome(', '.join(f'{len(group["spectra"])} spectra'
                              + (f' ({group["key"]})' if group['key'] else '')
                              for group in runs) or 'no spectrum')
        except (RuntimeError, ValueError, OSError) as err:
            log(f'{label}: {err}', 'warn')
            out[inst] = dict(error=str(err))
            continue
        for group in runs:
            name = label + (f' {group["key"]}' if group['key'] else '')
            try:
                if get:
                    step(f'{name}: the spectra')
                    got = fetch(target, inst, group, lbl_dir)
                    outcome(f'{got["n"]} spectra ({got["fetched"]} fetched'
                            + (f', {got["failed"]} not: {got["why"]}'
                               if got['failed'] else '') + ')')
                    if not got['n']:
                        out[inst][group['key']] = dict(n=0,
                                                       error='no spectrum')
                        continue
                if compute:
                    files = run_group(target, inst, group, lbl_dir, teff,
                                      python, nproc=nproc, tellu=tellu,
                                      template=template, root=root,
                                      ident=ident)
                    out[inst][group['key']] = dict(files=files)
            except (RuntimeError, ValueError, OSError) as err:
                log(f'{name}: {err}', 'warn')
                out[inst][group['key']] = dict(error=str(err))
    return out

# =============================================================================
# End of code
# =============================================================================

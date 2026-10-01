#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Who the star is, and what is known of its planets.

- resolve(name): the names of a star, from CDS's Sesame (SIMBAD): its main
  identifier and aliases (GJ, HD, HIP, TIC, Gaia DR3...), whatever the
  spelling in a file header ('Gl687' is GJ 687);
- known_planets(name): its planets in the NASA Exoplanet Archive, found by
  their Gaia DR3 or TIC identifier (then by host name): the archive's
  default solution of each planet (pscomppars) and every published one
  that gives K (ps), with their references;
- compare(orbit, known): a fitted orbit against the planet of its period,
  its K against every published solution, the most recent first (the
  archive's default is often the discovery paper, from far fewer velocities
  than a star has today);
- conjunction(planet): the published time of conjunction of a planet, to
  carry its ephemeris to a new series; bibcode(url): the paper of a
  solution (koloa.literature finds its velocities on VizieR).

known_planets also keeps what the archive knows of the star (rotation
period, spectral type, Teff, radius, metallicity, v sin i, magnitudes).

The archive is fetched once and kept: its two tables (pscomppars, the
default solution of every planet, and the solutions of ps that give K, a
few MB in all) in ~/.cache/koloa/archive, where every star is then looked
up without the network. tables(refresh=True) (koloa --refresh-archive, or
--refresh with --detailed) fetches them again. known_planets also keeps
what it found for a star in a JSON file when given a path.

Created on 2026-09-29

@author: artigau
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request
from html import unescape
from typing import Any, Dict, List, Optional

import numpy as np

from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
SESAME = 'https://cds.unistra.fr/cgi-bin/nph-sesame/-oxI/S'
ARCHIVE = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync'
#: where the archive's tables are kept, once fetched
CACHE = os.path.join(os.path.expanduser('~'), '.cache', 'koloa', 'archive')
#: the columns of pscomppars kept (the ones that find a star, and the ones
#: known_planets gives)
PSCOMP_COLUMNS = ['hostname', 'gaia_dr3_id', 'tic_id', 'pl_name',
                  'pl_orbper', 'pl_orbpererr1', 'pl_orbpererr2', 'pl_rvamp',
                  'pl_rvamperr1', 'pl_rvamperr2', 'pl_orbeccen', 'pl_bmasse',
                  'pl_rvamp_reflink', 'st_rotperr1', 'st_rotperr2',
                  'pl_tranmiderr1', 'pl_tranmiderr2', 'pl_orbtpererr1',
                  'pl_orbtpererr2']
#: the columns of ps kept (its solutions that give K)
PS_COLUMNS = ['hostname', 'pl_name', 'pl_orbper', 'pl_rvamp', 'pl_rvamperr1',
              'pl_rvamperr2', 'pl_refname', 'pl_pubdate']
#: the columns of the TESS Objects of Interest kept (the toi table)
TOI_COLUMNS = ['toi', 'tid', 'tfopwg_disp', 'pl_orbper', 'pl_orbpererr1',
               'pl_orbpererr2', 'pl_tranmid', 'pl_tranmiderr1',
               'pl_tranmiderr2', 'pl_trandurh', 'pl_trandep', 'pl_rade']
#: the dispositions of a TOI that is not a planet (false positive, false
#: alarm), left out unless asked for by its number
TOI_NOT_PLANETS = ('FP', 'FA')
#: the tables, once read
_TABLES: Optional[Dict[str, Any]] = None
#: a fitted signal is a known planet when their periods are within this
#: fraction (the periods of old solutions can be off by a few per cent)
MATCH = 0.05
#: what a file of known_planets holds: 2 adds the star (rotation, type...)
#: and the ephemeris and discovery of each planet
SCHEMA = 2
#: the columns of the star kept from pscomppars, and their names in koloa
STAR_COLUMNS = dict(st_spectype='spectral_type', st_teff='teff',
                    st_mass='mass', st_rad='radius', st_met='metallicity',
                    st_logg='logg', st_vsin='vsini', st_rotp='rotation',
                    st_age='age', sy_dist='distance', sy_vmag='V',
                    sy_jmag='J', sy_kmag='K')
#: the columns of each planet kept from pscomppars (beyond P, K, e, mass)
PLANET_COLUMNS = dict(pl_tranmid='tc', pl_orbtper='tp', pl_orblper='omega',
                      pl_orbsmax='a', pl_eqt='teq', pl_insol='insolation',
                      discoverymethod='discovery', disc_year='disc_year',
                      disc_facility='disc_facility')


# =============================================================================
# Define functions
# =============================================================================
def resolve(name: str, timeout: float = 30.0, refresh: bool = False
            ) -> Dict[str, Any]:
    """
    The names of a star (CDS Sesame, SIMBAD's resolver, which takes most
    spellings: Gl687, GJ 687, HD69830, Proxima...), kept on disk once asked
    (CACHE/sesame), since they do not change

    :param name: str, a name
    :param timeout: float [s]
    :param refresh: bool, ask Sesame again

    :return: dict, name (as given), main (SIMBAD's main identifier), aliases
             (every identifier), and the ones koloa uses: gaia_dr3, tic,
             hip, hd, gj (None when SIMBAD has none), and ra, dec (J2000,
             degrees)
    """
    kept = os.path.join(CACHE, 'sesame', re.sub(r'[^A-Za-z0-9+\-.]+', '_',
                                                name.strip()) + '.json')
    if os.path.exists(kept) and not refresh:
        with open(kept) as handle:
            return dict(json.load(handle), name=name)
    url = f'{SESAME}?{urllib.parse.quote(name)}'
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        text = resp.read().decode('utf-8', 'replace')
    main = re.search(r'<oname>(.*?)</oname>', text)
    if main is None:
        raise ValueError(f'SIMBAD does not know {name}')
    aliases = [' '.join(alias.split())
               for alias in re.findall(r'<alias>(.*?)</alias>', text)]

    def first(prefix):
        return next((alias for alias in aliases
                     if alias.startswith(prefix + ' ')), None)
    # the position (J2000, degrees), for planning observations
    radeg = re.search(r'<jradeg>(.*?)</jradeg>', text)
    dedeg = re.search(r'<jdedeg>(.*?)</jdedeg>', text)
    out = dict(name=name, main=' '.join(main.group(1).split()),
               aliases=aliases, gaia_dr3=first('Gaia DR3'),
               tic=first('TIC'), hip=first('HIP'), hd=first('HD'),
               gj=first('GJ'),
               ra=float(radeg.group(1)) if radeg else None,
               dec=float(dedeg.group(1)) if dedeg else None)
    os.makedirs(os.path.dirname(kept), exist_ok=True)
    with open(kept, 'w') as handle:
        json.dump(out, handle)
    return out


def _query(query: str, timeout: float = 60.0) -> List[Dict[str, Any]]:
    """the rows of a query of the NASA Exoplanet Archive (TAP, JSON)"""
    url = ARCHIVE + '?' + urllib.parse.urlencode(dict(query=query,
                                                      format='json'))
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return json.loads(response.read())


def _key(name: Any) -> str:
    """a name compared without its case and spaces"""
    return re.sub(r'\s+', '', str(name)).lower()


def tables(refresh: bool = False) -> Dict[str, Any]:
    """
    The NASA Exoplanet Archive, kept: pscomppars and the solutions of ps
    that give K, fetched the first time (or with refresh) into CACHE and
    read from there after

    :param refresh: bool, fetch the tables again

    :return: dict, fetched (date), pscomppars and ps (lists of rows), and
             the indexes of the hosts (by Gaia DR3, TIC and name)
    """
    global _TABLES
    path = os.path.join(CACHE, 'tables.json')
    if _TABLES is not None and not refresh:
        return _TABLES
    if os.path.exists(path) and not refresh:
        with open(path) as handle:
            out = json.load(handle)
        if 'toi' not in out:
            # a copy from before the TOIs were kept: they are added (and
            #   without the network, the copy goes on without them)
            log('fetching the TESS Objects of Interest of the archive')
            try:
                out['toi'] = _query(f'select {",".join(TOI_COLUMNS)} from '
                                    f'toi', timeout=600)
                with open(path + '.part', 'w') as handle:
                    json.dump(out, handle)
                os.replace(path + '.part', path)
            except Exception as err:
                log(f'the TOIs could not be fetched ({err}): none for now',
                    'warn')
                out['toi'] = []
    else:
        log(f'fetching the NASA Exoplanet Archive (a few MB, once; kept in '
            f'{CACHE})')
        cols = ','.join(PSCOMP_COLUMNS + [col for col in
                                          list(STAR_COLUMNS)
                                          + list(PLANET_COLUMNS)
                                          if col not in PSCOMP_COLUMNS])
        out = dict(fetched=time.strftime('%Y-%m-%d %H:%M'),
                   pscomppars=_query(f'select {cols} from pscomppars',
                                     timeout=600),
                   ps=_query(f'select {",".join(PS_COLUMNS)} from ps where '
                             f'pl_rvamp is not null', timeout=600),
                   toi=_query(f'select {",".join(TOI_COLUMNS)} from toi',
                              timeout=600))
        os.makedirs(CACHE, exist_ok=True)
        with open(path + '.part', 'w') as handle:
            json.dump(out, handle)
        os.replace(path + '.part', path)
        with open(os.path.join(CACHE, 'fetched.txt'), 'w') as handle:
            handle.write(out['fetched'] + '\n')
        log(f'NASA Exoplanet Archive: {len(out["pscomppars"])} planets, '
            f'{len(out["ps"])} published solutions with K and '
            f'{len(out["toi"])} TOIs, kept', 'value')
    index: Dict[str, Dict[str, str]] = dict(gaia={}, tic={}, name={})
    for row in out['pscomppars']:
        for key, col in (('gaia', 'gaia_dr3_id'), ('tic', 'tic_id'),
                         ('name', 'hostname')):
            if row.get(col):
                index[key][_key(row[col])] = row['hostname']
    out['index'] = index
    _TABLES = out
    return out


def tois(tic: Any) -> List[Dict[str, Any]]:
    """
    The TESS Objects of Interest of a star (the toi table of the archive,
    its copy kept here): the ephemeris of each from the transits

    :param tic: str or int, its TIC number ('TIC 307210830' will do)

    :return: list of dict: toi ('175.01'), disposition (CP, KP, PC, APC,
             FP, FA...), P and P_err [d], tc and tc_err (the time of a
             transit, BJD - 2400000), duration [h], depth [ppm],
             radius [Earth radii]; by TOI number
    """
    number = int(str(tic).replace('TIC', '').strip())

    def err(row, col):
        vals = [abs(row[key]) for key in (col + '1', col + '2')
                if row.get(key) is not None]
        return float(np.mean(vals)) if vals else None
    out = []
    for row in tables()['toi']:
        if row.get('tid') != number or row.get('pl_orbper') is None:
            continue
        tc = row.get('pl_tranmid')
        out.append(dict(
            toi=f'{float(row["toi"]):.2f}', disposition=row.get('tfopwg_disp'),
            P=float(row['pl_orbper']), P_err=err(row, 'pl_orbpererr'),
            tc=None if tc is None else float(tc) - 2400000.0,
            tc_err=err(row, 'pl_tranmiderr'), duration=row.get('pl_trandurh'),
            depth=row.get('pl_trandep'), radius=row.get('pl_rade')))
    return sorted(out, key=lambda item: float(item['toi']))


def fetched() -> Optional[str]:
    """when the kept archive was fetched (None: never)"""
    path = os.path.join(CACHE, 'fetched.txt')
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return handle.read().strip()


def _reference(link: str):
    """the text and the ADS address of the archive's <a href=...> Author et
    al. year </a>"""
    ref = re.search(r'>\s*(.*?)\s*</a>', link or '')
    url = re.search(r'href=(\S+)', link or '')
    # the archive writes accents as HTML entities (L&oacute;pez-Morales)
    return (unescape(ref.group(1)) if ref else ''), (url.group(1) if url else '')


def _error(row: Dict[str, Any], prefix: str) -> Optional[float]:
    """the mean of the archive's two errors of a column (<prefix>1 and
    <prefix>2), or None"""
    errs = [abs(row[key]) for key in (f'{prefix}1', f'{prefix}2')
            if row.get(key) is not None]
    return float(np.mean(errs)) if errs else None


def bibcode(url: str) -> Optional[str]:
    """
    The bibcode of an ADS address (the archive's links to the papers)

    :param url: str, e.g. https://ui.adsabs.harvard.edu/abs/2023A&A...680A..28G/abstract

    :return: str or None
    """
    found = re.search(r'/abs/([^/?#\s]+)', url or '')
    if found is None:
        return None
    code = urllib.parse.unquote(found.group(1))
    return code if len(code) == 19 else None


def conjunction(planet: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    The time of conjunction of a known planet, for its ephemeris: the
    transit when the archive has it, or the time of periastron of a circular
    orbit (for which it is the conjunction, omega = 90 deg: the convention
    of most papers, not of all)

    :param planet: dict, one of known_planets()['planets']

    :return: dict or None, tc and tc_err [BJD - 2400000], P and P_err
             [days], source (transit or periastron)
    """
    if not planet.get('P'):
        return None
    base = dict(P=float(planet['P']), P_err=planet.get('P_err'))
    if planet.get('tc') is not None:
        return dict(base, tc=planet['tc'], tc_err=planet.get('tc_err'),
                    source='transit')
    if planet.get('tp') is not None and not planet.get('e'):
        return dict(base, tc=planet['tp'], tc_err=planet.get('tp_err'),
                    source='periastron of a circular orbit')
    return None


def host_name(ident: Dict[str, Any]) -> Optional[str]:
    """
    The archive's name of a star: by its Gaia DR3 identifier, its TIC, then
    its aliases as host names

    :param ident: dict, from resolve()

    :return: str or None (the archive has no planet of it)
    """
    index = tables()['index']
    tries = [('gaia', ident.get('gaia_dr3')), ('tic', ident.get('tic'))]
    tries += [('name', alias) for alias in
              [ident.get('main'), ident.get('name')] + ident.get('aliases',
                                                                 [])]
    for kind, value in tries:
        if value and _key(value) in index[kind]:
            return index[kind][_key(value)]
    return None


def solutions(host: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Every published solution of the planets of a host that gives K (the ps
    table of the archive)

    :return: dict, planet name -> list of (reference, reference_url, date,
             P, K, K_err), oldest first
    """
    rows = [row for row in tables()['ps'] if row['hostname'] == host]
    out: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        if row.get('pl_rvamp') is None:
            continue
        kerr = [abs(row[key]) for key in ('pl_rvamperr1', 'pl_rvamperr2')
                if row.get(key) is not None]
        ref, url = _reference(row.get('pl_refname'))
        out.setdefault(row['pl_name'], []).append(dict(
            reference=ref, reference_url=url, date=row.get('pl_pubdate') or '',
            P=row['pl_orbper'], K=row['pl_rvamp'],
            K_err=float(np.mean(kerr)) if kerr else None))
    for sols in out.values():
        sols.sort(key=lambda sol: sol['date'])
    return out


def known_planets(name: Optional[str] = None, host: Optional[str] = None,
                  path: Optional[str] = None,
                  refresh: bool = False) -> Dict[str, Any]:
    """
    The planets of a star in the NASA Exoplanet Archive: the default
    solution of each (pscomppars), every published solution that gives K
    (ps), and the mass and distance of the star

    :param name: str or None, any name of the star (resolved by Sesame)
    :param host: str or None, the archive's host name, when known
    :param path: str or None, a JSON file that keeps the answer (read back
                 when it exists, unless refresh)
    :param refresh: bool, fetch the archive again (tables(refresh=True)),
                    even when the file exists

    :return: dict, host, fetched (date), star (mass, distance, and the
             STAR_COLUMNS the archive has: rotation, spectral type, Teff...),
             planets (name, P, P_err, K, K_err, e, mass_earth, reference,
             reference_url, solutions, and the PLANET_COLUMNS: tc and tp
             [BJD - 2400000] with their errors, a, teq, discovery...); no
             planet when the archive has none
    """
    if path and os.path.exists(path) and not refresh:
        with open(path) as handle:
            out = json.load(handle)
        if out.get('schema', 1) >= SCHEMA or not out.get('host'):
            return out
        # a file of an older koloa: ask again, for what it lacks
        host = out['host']
    kept = tables(refresh=refresh)
    if host is None:
        if name is None:
            raise ValueError('known_planets needs a name or a host')
        host = host_name(resolve(name))
    out = dict(host=host, fetched=kept['fetched'], star={}, planets=[],
               schema=SCHEMA)
    if host is not None:
        rows = sorted((row for row in kept['pscomppars']
                       if row['hostname'] == host),
                      key=lambda row: (row.get('pl_orbper') is None,
                                       row.get('pl_orbper') or 0))
        sols = solutions(host)
        if rows:
            out['star'] = {name: rows[0].get(col)
                           for col, name in STAR_COLUMNS.items()}
            out['star']['rotation_err'] = _error(rows[0], 'st_rotperr')
        for row in rows:
            kerr = [abs(row[key]) for key in ('pl_rvamperr1', 'pl_rvamperr2')
                    if row.get(key) is not None]
            perr = [abs(row[key]) for key in ('pl_orbpererr1',
                                              'pl_orbpererr2')
                    if row.get(key) is not None]
            ref, url = _reference(row.get('pl_rvamp_reflink'))
            planet = dict(
                name=row['pl_name'], P=row['pl_orbper'],
                P_err=float(np.mean(perr)) if perr else None,
                K=row['pl_rvamp'],
                K_err=float(np.mean(kerr)) if kerr else None,
                e=row['pl_orbeccen'], mass_earth=row['pl_bmasse'],
                reference=ref, reference_url=url,
                solutions=sols.get(row['pl_name'], []))
            planet.update({name: row.get(col)
                           for col, name in PLANET_COLUMNS.items()})
            # the times in koloa's BJD - 2400000
            for key, col in (('tc', 'pl_tranmiderr'), ('tp', 'pl_orbtpererr')):
                if planet[key] is not None:
                    planet[key] = float(planet[key]) - 2400000.0
                planet[f'{key}_err'] = _error(row, col)
            out['planets'].append(planet)
    if path:
        with open(path, 'w') as handle:
            json.dump(out, handle, indent=1)
    log(f'{len(out["planets"])} planets of {host or name} in the NASA '
        f'Exoplanet Archive (kept, of {kept["fetched"]})', 'value')
    return out


def compare(orb: Dict[str, Any], known: Dict[str, Any],
            match: float = MATCH) -> Dict[str, Any]:
    """
    A fitted orbit against the known planet of its period: its K against
    every published solution, the most recent first

    :param orb: dict, the orbit (P, K, ... as (median, minus, plus))
    :param known: dict, from known_planets()
    :param match: float, the largest relative difference of the periods

    :return: dict, the orbit with known (the planet), comparisons (per
             solution: reference, date, K, K_err, z) and K_z (against the
             most recent); unchanged when no planet has its period
    """
    found = [pl for pl in known.get('planets', [])
             if pl['P'] and abs(pl['P'] / orb['P'][0] - 1) < match]
    if not found:
        return orb
    orb['known'] = pl = found[0]
    width = 0.5 * (orb['K'][1] + orb['K'][2])
    sols = pl.get('solutions') or [dict(reference=pl['reference'],
                                        reference_url=pl['reference_url'],
                                        date='', P=pl['P'], K=pl['K'],
                                        K_err=pl['K_err'])]
    orb['comparisons'] = [
        dict(sol, z=float((orb['K'][0] - sol['K'])
                          / np.hypot(width, sol['K_err'])))
        for sol in reversed(sols) if sol['K'] is not None and sol['K_err']]
    if orb['comparisons']:
        orb['K_z'] = orb['comparisons'][0]['z']
    return orb


# =============================================================================
# End of code
# =============================================================================

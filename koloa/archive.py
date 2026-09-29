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
  than a star has today).

Everything needs the network; known_planets keeps what it found in a JSON
file when given a path, and reads it back.

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
#: a fitted signal is a known planet when their periods are within this
#: fraction (the periods of old solutions can be off by a few per cent)
MATCH = 0.05


# =============================================================================
# Define functions
# =============================================================================
def resolve(name: str, timeout: float = 30.0) -> Dict[str, Any]:
    """
    The names of a star (CDS Sesame, SIMBAD's resolver, which takes most
    spellings: Gl687, GJ 687, HD69830, Proxima...)

    :param name: str, a name
    :param timeout: float [s]

    :return: dict, name (as given), main (SIMBAD's main identifier), aliases
             (every identifier), and the ones koloa uses: gaia_dr3, tic,
             hip, hd, gj (None when SIMBAD has none)
    """
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
    return dict(name=name, main=' '.join(main.group(1).split()),
                aliases=aliases, gaia_dr3=first('Gaia DR3'),
                tic=first('TIC'), hip=first('HIP'), hd=first('HD'),
                gj=first('GJ'))


def _query(query: str, timeout: float = 60.0) -> List[Dict[str, Any]]:
    """the rows of a query of the NASA Exoplanet Archive (TAP, JSON)"""
    url = ARCHIVE + '?' + urllib.parse.urlencode(dict(query=query,
                                                      format='json'))
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return json.loads(response.read())


def _reference(link: str):
    """the text and the ADS address of the archive's <a href=...> Author et
    al. year </a>"""
    ref = re.search(r'>\s*(.*?)\s*</a>', link or '')
    url = re.search(r'href=(\S+)', link or '')
    # the archive writes accents as HTML entities (L&oacute;pez-Morales)
    return (unescape(ref.group(1)) if ref else ''), (url.group(1) if url else '')


def host_name(ident: Dict[str, Any]) -> Optional[str]:
    """
    The archive's name of a star: by its Gaia DR3 identifier, its TIC, then
    its aliases as host names

    :param ident: dict, from resolve()

    :return: str or None (the archive has no planet of it)
    """
    tries = [('gaia_dr3_id', ident.get('gaia_dr3')),
             ('tic_id', ident.get('tic'))]
    tries += [('hostname', alias) for alias in
              [ident.get('main'), ident.get('name')] + ident.get('aliases',
                                                                 [])]
    for column, value in tries:
        if not value:
            continue
        safe = str(value).replace("'", "''")
        rows = _query(f"select hostname from pscomppars where {column} = "
                      f"'{safe}'")
        if rows:
            return rows[0]['hostname']
    return None


def solutions(host: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Every published solution of the planets of a host that gives K (the ps
    table of the archive)

    :return: dict, planet name -> list of (reference, reference_url, date,
             P, K, K_err), oldest first
    """
    rows = _query(
        f"select pl_name,pl_orbper,pl_rvamp,pl_rvamperr1,pl_rvamperr2,"
        f"pl_refname,pl_pubdate from ps where hostname = '{host}'")
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
    :param refresh: bool, ask the archive even when the file exists

    :return: dict, host, fetched (date), star (mass, distance), planets
             (name, P, P_err, K, K_err, e, mass_earth, reference,
             reference_url, solutions); no planet when the archive has none
    """
    if path and os.path.exists(path) and not refresh:
        with open(path) as handle:
            out = json.load(handle)
        # a file written before the solutions were kept: add them
        if out.get('host') and not all('solutions' in pl
                                       for pl in out['planets']):
            sols = solutions(out['host'])
            for pl in out['planets']:
                pl['solutions'] = sols.get(pl['name'], [])
            with open(path, 'w') as handle:
                json.dump(out, handle, indent=1)
        return out
    if host is None:
        if name is None:
            raise ValueError('known_planets needs a name or a host')
        host = host_name(resolve(name))
    out = dict(host=host, fetched=time.strftime('%Y-%m-%d'), star={},
               planets=[])
    if host is not None:
        cols = ('pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_rvamp,'
                'pl_rvamperr1,pl_rvamperr2,pl_orbeccen,pl_bmasse,'
                'pl_rvamp_reflink,st_mass,sy_dist')
        rows = _query(f"select {cols} from pscomppars where hostname = "
                      f"'{host}' order by pl_orbper")
        sols = solutions(host)
        if rows:
            out['star'] = dict(mass=rows[0]['st_mass'],
                               distance=rows[0]['sy_dist'])
        for row in rows:
            kerr = [abs(row[key]) for key in ('pl_rvamperr1', 'pl_rvamperr2')
                    if row.get(key) is not None]
            perr = [abs(row[key]) for key in ('pl_orbpererr1',
                                              'pl_orbpererr2')
                    if row.get(key) is not None]
            ref, url = _reference(row.get('pl_rvamp_reflink'))
            out['planets'].append(dict(
                name=row['pl_name'], P=row['pl_orbper'],
                P_err=float(np.mean(perr)) if perr else None,
                K=row['pl_rvamp'],
                K_err=float(np.mean(kerr)) if kerr else None,
                e=row['pl_orbeccen'], mass_earth=row['pl_bmasse'],
                reference=ref, reference_url=url,
                solutions=sols.get(row['pl_name'], [])))
    if path:
        with open(path, 'w') as handle:
            json.dump(out, handle, indent=1)
    log(f'{len(out["planets"])} planets of {host or name} in the NASA '
        f'Exoplanet Archive', 'value')
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

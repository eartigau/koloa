#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Public radial velocities from DACE, for koloa.

DACE (the Data and Analysis Center for Exoplanets, University of Geneva)
serves the velocities of the pipelines of HARPS, HARPS-N, ESPRESSO, CORALIE
and others. fetch() asks it for the time series of a target, as the
dace-query client does (Spectroscopy.get_timeseries), anonymously (public
data only) unless an API key is given, and writes it to a CSV file that
later runs read instead. rvdata() turns that file into an RVData with one
instrument per DACE instrument era: HARPS03 and HARPS15 (before and after
the 2015 fibre upgrade), ESPRESSO18 and ESPRESSO19 (before and after the
2019 intervention), CORALIE98, CORALIE07 and CORALIE14, HARPN, and so on.
Each era has its own zero point, so each is an instrument for koloa.

DACE filters some networks (the TLS handshake then stalls): fetch from
another network, or keep the CSV file and work from it.

    from koloa.dace import fetch, rvdata
    path = fetch('HD69830', 'hd69830_dace.csv')     # once, with the network
    data = rvdata(path, name='HD 69830')            # an RVData, per era

Created on 2026-09-29

@author: artigau
"""
import csv
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from koloa.data import EMPTY_CELLS, RVData, find_sequences, split_sequences
from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: the spectroscopy service of DACE (the one dace-query uses)
SPECTROSCOPY_API = 'https://spectroscopy-webapp-pg.obsuksprd2.unige.ch'
#: the velocity products: the pipeline (DRS) velocities by default
SOURCES = {'standard': 'POSTDRS_A', 'telluric': 'POSTDRS_TELL_CORR_A'}
#: the indicators kept, DACE column: koloa name
INDICATORS = {'ccf_fwhm': 'fwhm', 'ccf_bispan': 'bis', 'ccf_contrast': 'contrast',
              'spectro_rhk': 'rhk', 'spectro_smw': 'smw',
              'spectro_halpha': 'halpha'}


# =============================================================================
# Define functions
# =============================================================================
def _columns(payload) -> Dict[str, list]:
    """
    Decode DACE's parameter lists into columns (its values are run-length
    encoded: values[i] repeats occurrences[i] times)

    :param payload: dict, the JSON DACE returns

    :return: dict, name: values
    """
    columns = {}
    if not isinstance(payload, dict):
        return columns
    for par in payload.get('parameters') or []:
        name = par.get('variableName')
        if not name:
            continue
        values = next((par[key] for key in ('doubleValues', 'floatValues',
                                             'intValues', 'stringValues',
                                             'boolValues')
                       if par.get(key) is not None), [])
        occ = par.get('occurrences')
        if occ:
            values = [val for val, num in zip(values, occ) for _ in range(num)]
        columns.setdefault(name, []).extend(values)
        errors = par.get('minErrorValues')
        if errors is not None:
            if occ:
                errors = [val for val, num in zip(errors, occ)
                          for _ in range(num)]
            columns.setdefault(f'{name}_err', []).extend(errors)
    return columns


def fetch(target: str, path: str, api_key: Optional[str] = None,
          source: str = 'standard', limit: int = 20000,
          timeout: float = 120.0, refresh: bool = False) -> str:
    """
    The radial velocities of a target from DACE, written to a CSV file

    Anonymous requests see the public data only; an API key (DACE, 'My
    account') adds what that account may see.

    :param target: str, the name DACE knows the target by (e.g. HD69830)
    :param path: str, the CSV file written (and read back, if it exists and
                 refresh is False)
    :param api_key: str or None, a DACE API key (anonymous if None)
    :param source: str, standard (the pipeline velocities, POSTDRS_A) or
                   telluric (the pipeline's telluric-corrected ones)
    :param limit: int, the most rows asked for
    :param timeout: float, the longest wait for DACE [s]
    :param refresh: bool, ask DACE even when the file exists

    :return: str, the CSV file
    """
    if os.path.exists(path) and not refresh:
        log(f'DACE velocities of {target} from {path}', 'info')
        return path
    query = urllib.parse.urlencode({
        'limit': limit,
        'filters': json.dumps({'source_product_file_ext':
                               {'equals': [SOURCES[source]]}}),
        'sort': json.dumps({})})
    url = (f'{SPECTROSCOPY_API}/target/'
           f'{urllib.parse.quote(target, safe="")}/timeseries/'
           f'radial-velocities?{query}')
    request = urllib.request.Request(url, headers={
        'Accept': 'application/json', 'User-Agent': 'koloa'})
    if api_key:
        request.add_header('Authorization', f'apiKey:{api_key}')
    log(f'asking DACE for the velocities of {target}', 'info')
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read())
    except urllib.error.HTTPError as err:
        # DACE answers: it does not know the target under this name (404),
        #   which is not a network that filters it
        if err.code == 404:
            raise RuntimeError(f'DACE does not know {target}') from err
        raise RuntimeError(f'DACE could not be reached ({err}); some '
                           f'networks are filtered: fetch from another '
                           f'one, or give the CSV file') from err
    except (urllib.error.URLError, TimeoutError) as err:
        raise RuntimeError(f'DACE could not be reached ({err}); some '
                           f'networks are filtered: fetch from another '
                           f'one, or give the CSV file') from err
    columns = _columns(payload)
    nrow = max((len(val) for val in columns.values()), default=0)
    if nrow == 0:
        raise RuntimeError(f'DACE has no velocities of {target} that this '
                           f'account may see')
    names = list(columns)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for irow in range(nrow):
            writer.writerow([columns[name][irow]
                             if irow < len(columns[name]) else ''
                             for name in names])
    log(f'{nrow} velocities of {target} written to {path}', 'value')
    return path


def names(ident: Dict[str, Any], target: str) -> List[str]:
    """the names DACE may know a star by: its catalogue names without
    spaces (HD69830, GJ687, HIP86162...), then as SIMBAD writes them"""
    out = []
    for key in ('hd', 'gj', 'hip', 'main'):
        if ident.get(key):
            out.append(ident[key].replace('NAME ', '').replace(' ', ''))
    out.append(target.replace(' ', ''))
    for alias in ident.get('aliases', []):
        if alias.split()[0] in ('GJ', 'Gl', 'HD', 'HIP', 'TOI', 'K2',
                                'Kepler', 'LHS', 'Wolf', 'Ross'):
            out.append(alias.replace(' ', ''))
    return list(dict.fromkeys(name for name in out if name))


def _float(value) -> float:
    """A number, or nan"""
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def rvdata(path: str, name: Optional[str] = None, min_points: int = 5,
           instruments: Optional[Sequence[str]] = None,
           exclude: Sequence[str] = (), max_visit: Optional[int] = 10
           ) -> RVData:
    """
    An RVData from the DACE velocities of fetch(), one instrument per DACE
    instrument era

    Only the velocities that pass the quality control of the pipeline and
    have an error are kept, one per spectrum. An instrument era observed in
    several modes (HARPS and its EGGS mode, the ESPRESSO resolutions) keeps
    each mode of at least min_points velocities as an instrument of its
    own: the main mode keeps the name of the era, the others are named
    era_mode; smaller groups are dropped.

    :param path: str, the CSV file written by fetch()
    :param name: str or None, the name of the target
    :param min_points: int, the fewest velocities an instrument keeps
    :param instruments: list of str or None, the instruments kept (all when
                        None)
    :param exclude: list of str, the instruments left out
    :param max_visit: int or None, visits of more exposures are split into
                      consecutive ones of at most this many
                      (koloa.data.split_sequences), so that every exposure
                      and every visit can be an outlier; None keeps them

    :return: RVData, the series (m/s, BJD - 2400000), with the indicators of
             INDICATORS and every other column of DACE in its meta (S/N per
             order, exposure time, drift noise, ...)
    """
    with open(path) as handle:
        rows = list(csv.DictReader(handle))
    # the quality control of the pipeline, a velocity and its error
    good = []
    seen = set()
    for row in rows:
        if str(row.get('drs_qc', 'True')).lower() not in ('true', '1'):
            continue
        tt, vv, ee = (_float(row.get(key)) for key in ('rjd', 'rv', 'rv_err'))
        if not (np.isfinite(tt) and np.isfinite(vv) and np.isfinite(ee)
                and ee > 0):
            continue
        # one velocity per spectrum
        key = row.get('file_rootname') or (row.get('instrument_name'),
                                           round(tt, 6))
        if key in seen:
            continue
        seen.add(key)
        good.append(row)
    # the instrument of each velocity: its era, and its mode when the era
    #   has several
    era = np.array([row.get('instrument_name', 'DACE') for row in good])
    mode = np.array([row.get('ins_mode', '') or '' for row in good])
    # objects, not fixed-width strings: era_mode is longer than era
    label = np.array(era, dtype=object)
    for ename in np.unique(era):
        sel = era == ename
        modes, counts = np.unique(mode[sel], return_counts=True)
        # the main mode keeps the name of the era, the others are suffixed
        main = modes[np.argmax(counts)]
        for mname in modes:
            if mname != main:
                label[sel & (mode == mname)] = f'{ename}_{mname}'
    names, counts = np.unique(label.astype(str), return_counts=True)
    keep = {nm for nm, num in zip(names, counts) if num >= min_points}
    if instruments is not None:
        keep &= set(instruments)
    keep -= set(exclude)
    dropped = sorted(set(names) - keep)
    if dropped:
        log(f'left out: {", ".join(dropped)} (fewer than {min_points} '
            f'velocities, or excluded)', 'warn')
    label = label.astype(str)
    sel = np.isin(label, sorted(keep))
    rows_kept = [row for row, ok in zip(good, sel) if ok]
    time = np.array([_float(row['rjd']) for row in rows_kept])
    rv = np.array([_float(row['rv']) for row in rows_kept])
    err = np.array([_float(row['rv_err']) for row in rows_kept])
    indicators = {}
    for col, key in INDICATORS.items():
        val = np.array([_float(row.get(col)) for row in rows_kept])
        vale = np.array([_float(row.get(f'{col}_err')) for row in rows_kept])
        if np.sum(np.isfinite(val) & np.isfinite(vale)) >= min_points:
            indicators[key] = (val, vale)
    # every other column, as numbers where most of the filled cells are
    meta = {}
    for col in (rows_kept[0].keys() if rows_kept else []):
        if col in ('rjd', 'rv', 'rv_err'):
            continue
        cells = ['' if row.get(col) is None else str(row.get(col)).strip()
                 for row in rows_kept]
        values = np.array([_float(cell) for cell in cells])
        filled = sum(cell.lower() not in EMPTY_CELLS for cell in cells)
        meta[col] = (values if np.sum(np.isfinite(values)) >= 0.5 * filled
                     else np.array(cells))
    # the visits, in time order, long ones split
    order = np.argsort(time, kind='stable')
    time, rv, err = time[order], rv[order], err[order]
    inst = label[sel][order]
    indicators = {key: (val[0][order], val[1][order])
                  for key, val in indicators.items()}
    meta = {key: val[order] for key, val in meta.items()}
    seq = find_sequences(time, inst)
    if max_visit:
        longest = int(np.max(np.bincount(seq))) if len(seq) else 0
        if longest > max_visit:
            seq = split_sequences(seq, max_visit)
            log(f'visits of more than {max_visit} exposures (the longest '
                f'has {longest}) split into consecutive ones', 'warn')
    data = RVData(time, rv, err, inst=inst, seq=seq, indicators=indicators,
                  name=name or os.path.splitext(os.path.basename(path))[0],
                  meta=meta)
    log(f'{data.n} velocities in {data.nseq} visits: ' + ', '.join(
        f'{inst} {np.sum(data.inst == inst)}' for inst in data.instruments),
        'value')
    return data


def instruments_table(data: RVData) -> List[dict]:
    """
    The points, visits, years and median error of each instrument

    :param data: RVData, the series

    :return: list of dict, one per instrument
    """
    out = []
    for inst in data.instruments:
        sel = data.inst == inst
        years = 2000.0 + (data.time[sel] + 2400000.0 - 2451544.5) / 365.25
        out.append(dict(inst=inst, n=int(np.sum(sel)),
                        nvisits=int(len(np.unique(data.seq[sel]))),
                        years=[float(years.min()), float(years.max())],
                        median_err=float(np.median(data.err[sel])),
                        rms=float(np.std(data.rv[sel]))))
    return out

# =============================================================================
# End of code
# =============================================================================

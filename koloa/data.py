#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A radial velocity time series, as koloa sees it.

Four things are carried with the velocities because every outlier-aware
tool downstream needs them:

- the instrument of each point, because each instrument has its own offset
  and its own jitter;
- the observing sequence of each point (the exposures taken back to back
  in one visit), because a bad night moves all its exposures together and
  an outlier is then a sequence, not a point;
- the activity indicators, because a signal that is also in them is not a
  planet;
- the offset taken out of each instrument, so that the numbers stay small
  and the priors on the offsets can be centred on zero.

Created on 2026-09-27

@author: artigau
"""
import csv
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# =============================================================================
# Define variables
# =============================================================================
#: points closer than this in time belong to the same sequence [days]
SEQUENCE_GAP = 0.3
#: the largest gap inside a night [days]: two visits a few hours apart in
#:  the same night are one nightly mean
NIGHT_GAP = 0.5
#: koloa's analyses (oafip, duck_test, analyze, detailed_analysis) run on
#:  the nightly means unless told otherwise: nightly=False in a call,
#:  koloa.data.NIGHTLY = False, or the environment variable KOLOA_NIGHTLY=0
NIGHTLY = os.environ.get('KOLOA_NIGHTLY', '1').strip().lower() not in (
    '0', 'false', 'no', 'off')
#: the cells of a column that count as empty
EMPTY_CELLS = ('', 'nan', 'none', 'null', 'na', 'n/a')

#: names tried, in order, for the time, velocity and error columns
TIME_NAMES = ['rjd', 'bjd', 'mjd', 'jd', 'time', 't']
RV_NAMES = ['vrad', 'rv', 'vel', 'velocity', 'y']
ERR_NAMES = ['svrad', 'sig_vrad', 'erv', 'rv_err', 'e_rv', 'err', 'error',
             'sigma']
#: the names of a column that gives the instrument of each velocity
INST_NAMES = ('inst', 'instrument', 'ins_name', 'instrument_name')


# =============================================================================
# Define functions
# =============================================================================
def robust_std(value: np.ndarray) -> float:
    """
    1.4826 times the median absolute deviation, the sigma of a gaussian

    :param value: np.ndarray, the values

    :return: float, the robust standard deviation
    """
    value = np.asarray(value, dtype=float)
    value = value[np.isfinite(value)]
    if len(value) < 2:
        return np.nan
    return float(1.4826 * np.median(np.abs(value - np.median(value))))


#: the instruments by the first word of the names of their files (an LBL
#: file keeps them: r.HARPS.2018-..., r.ESPRE.2019-...)
FILE_INSTRUMENTS = {'HARPS': 'HARPS', 'HARPN': 'HARPN', 'ESPRE': 'ESPRESSO',
                    'ESPRESSO': 'ESPRESSO', 'NIRPS': 'NIRPS',
                    'CORALIE': 'CORALIE', 'SPIROU': 'SPIRou'}
#: the eras of an instrument, each with its own zero point, named as DACE
#: names them: (the first night of the next era [BJD - 2400000], names)
ERAS = {'HARPS': ([57174.5], ['HARPS03', 'HARPS15']),
        'ESPRESSO': ([58662.0], ['ESPRESSO18', 'ESPRESSO19'])}


def instrument_names(data: 'RVData') -> Optional[np.ndarray]:
    """
    The instrument of every point of a file of one instrument (LBL), from
    its own columns: the names of its files (r.HARPS..., r.ESPRE...,
    r.NIRPS...), else the S/N keys of APERO (EXTSN035: SPIRou, EXTSN060:
    NIRPS); an instrument whose zero point changed is split in its eras
    (HARPS03 and HARPS15, ESPRESSO18 and ESPRESSO19), as on DACE

    :return: np.ndarray of str, or None (nothing says which)
    """
    base = None
    for col in ('FILENAME', 'local_file_name', 'ARCFILE', 'filename'):
        vals = data.meta.get(col)
        if vals is None or np.asarray(vals).dtype.kind not in 'USO':
            continue
        words = [str(val).strip() for val in vals if str(val).strip()]
        if not words:
            continue
        word = words[0]
        word = word[2:] if word.startswith('r.') else word
        head = word.split('.')[0].split('_')[0].upper()
        base = FILE_INSTRUMENTS.get(head)
        if base:
            break
    if base is None:
        if 'EXTSN035' in data.meta:
            base = 'SPIRou'
        elif 'EXTSN060' in data.meta:
            base = 'NIRPS'
        else:
            return None
    names = np.full(data.n, base, dtype=object)
    if base in ERAS:
        edges, labels = ERAS[base]
        names = np.array(labels, dtype=object)[np.searchsorted(edges,
                                                                data.time)]
    return names.astype(str)


def night_index(data: 'RVData', gap: float = NIGHT_GAP) -> np.ndarray:
    """
    The night of each exposure (per instrument; see RVData.nightly), in the
    order of the nightly means

    :param data: RVData, the series
    :param gap: float, the largest gap inside a night [days]

    :return: np.ndarray, an integer per exposure
    """
    return find_sequences(data.time, data.inst, gap)


def find_sequences(time: np.ndarray, inst: Optional[np.ndarray] = None,
                   gap: float = SEQUENCE_GAP) -> np.ndarray:
    """
    Group points into observing sequences

    Two points are in the same sequence when they come from the same
    instrument and nothing separates them by more than `gap` in time. With
    the default of 0.3 day this is one visit per night, which is the scale
    on which a bad calibration, a passing cloud or moonlight moves every
    exposure at once.

    :param time: np.ndarray, the time [days]
    :param inst: np.ndarray or None, the instrument of each point
    :param gap: float, the largest gap inside a sequence [days]

    :return: np.ndarray, an integer sequence number per point (0, 1, ...,
             in time order)
    """
    time = np.asarray(time, dtype=float)
    if inst is None:
        inst = np.zeros(len(time), dtype=int)
    inst = np.asarray(inst)
    seq = np.full(len(time), -1, dtype=int)
    count = 0
    # one instrument at a time, so that two instruments on the same night
    #   are two sequences
    starts = []
    for name in np.unique(inst):
        idx = np.where(inst == name)[0]
        idx = idx[np.argsort(time[idx])]
        new = np.concatenate([[True], np.diff(time[idx]) > gap])
        labels = np.cumsum(new) - 1 + count
        seq[idx] = labels
        for lab in np.unique(labels):
            starts.append((float(np.min(time[idx][labels == lab])), lab))
        count = int(np.max(labels)) + 1
    # renumber in time order so that sequence 0 is the first visit
    order = {old: new for new, (_, old) in enumerate(sorted(starts))}
    return np.array([order[s] for s in seq], dtype=int)


def _pick(names: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    """
    The first of the candidate names that is a column

    :param names: list of str, the column names
    :param candidates: list of str, what to look for (any case)

    :return: str or None, the column found
    """
    lower = {name.lower(): name for name in names}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    return None


def _error_column(name: str, names: Sequence[str]) -> Optional[str]:
    """
    The error column that goes with a value column

    LBL writes sX, sig_X or X_err; any of them will do.

    :param name: str, the value column
    :param names: list of str, all the columns

    :return: str or None, the error column
    """
    for cand in ['s' + name, 'sig_' + name, 'sig' + name, name + '_err',
                 'e_' + name, 'err_' + name, 's_' + name]:
        if cand in names:
            return cand
    return None


def split_sequences(seq: np.ndarray, max_size: int = 10) -> np.ndarray:
    """
    Sequences longer than max_size split into consecutive chunks of nearly
    equal size

    koloa's mixture sums over every way the exposures of a visit can be
    outliers, which is affordable up to about ten exposures; a longer visit
    (a run at a cadence of minutes) becomes several shorter ones, each with
    the visit jitter, and no velocity is lost.

    :param seq: np.ndarray, the sequence of each point (0, 1, ..., in time
                order, as find_sequences gives them; the points in time
                order)
    :param max_size: int, the most points in a sequence

    :return: np.ndarray, the new sequence of each point (0, 1, ..., in time
             order)
    """
    seq = np.asarray(seq, dtype=int)
    out = np.empty_like(seq)
    count = 0
    for sval in np.unique(seq):
        idx = np.where(seq == sval)[0]
        for chunk in np.array_split(idx, int(np.ceil(len(idx) / max_size))):
            out[chunk] = count
            count += 1
    return out


def merge(series: Sequence['RVData'], name: Optional[str] = None
          ) -> 'RVData':
    """
    Several series in one: each keeps its instruments (their names must
    differ), its visits, its indicators and its other columns; an indicator
    or a column that a series lacks is nan (or empty text) for its
    exposures

    :param series: list of RVData
    :param name: str or None, the name (that of the first series when None)

    :return: RVData, in time order, visits numbered in time order
    """
    parts = [part for part in series if part is not None and part.n]
    if not parts:
        raise ValueError('merge: no series with data')
    seen = set()
    for part in parts:
        both = seen & set(part.instruments)
        if both:
            raise ValueError(f'merge: instruments in two series: '
                             f'{sorted(both)}')
        seen |= set(part.instruments)
    time = np.concatenate([part.time for part in parts])
    # the velocities as given: each series' zero points put back
    rv = np.concatenate([part.rv + np.array([part.zero_point[str(inst)]
                                             for inst in part.inst])
                         for part in parts])
    err = np.concatenate([part.err for part in parts])
    inst = np.concatenate([part.inst.astype(object) for part in parts])
    seq, offset = [], 0
    for part in parts:
        seq.append(part.seq + offset)
        offset += part.nseq
    seq = np.concatenate(seq)

    def column(get, text):
        out = []
        for part in parts:
            val = get(part)
            if val is None:
                val = (np.full(part.n, '', dtype=object) if text
                       else np.full(part.n, np.nan))
            out.append(np.asarray(val, dtype=object) if text
                       else np.asarray(val, dtype=float))
        return np.concatenate(out)
    indicators = {}
    for key in dict.fromkeys(key for part in parts for key in part.indicators):
        indicators[key] = (
            column(lambda part: part.indicators[key][0]
                   if key in part.indicators else None, False),
            column(lambda part: part.indicators[key][1]
                   if key in part.indicators else None, False))
    meta = {}
    for key in dict.fromkeys(key for part in parts for key in part.meta):
        text = any(part.meta[key].dtype.kind not in 'fiu'
                   for part in parts if key in part.meta)
        meta[key] = column(lambda part: part.meta.get(key), text)
        if text:
            meta[key] = np.array([str(val) for val in meta[key]])
    # in time order, the visits numbered in the order they start
    order = np.argsort(time, kind='stable')
    _, first = np.unique(seq[order], return_index=True)
    rank = np.empty(offset, dtype=int)
    rank[np.unique(seq[order])[np.argsort(first)]] = np.arange(len(first))
    # each part's zero points, of its own instruments only: a part cut out
    #   of a larger series (select) may still hold another's under the
    #   same name, which would shift that instrument here
    zero = {}
    for part in parts:
        zero.update({str(inst): part.zero_point[str(inst)]
                     for inst in part.instruments})
    return RVData(time[order], rv[order], err[order],
                  inst=inst[order].astype(str), seq=rank[seq[order]],
                  indicators={key: (val[0][order], val[1][order])
                              for key, val in indicators.items()},
                  meta={key: val[order] for key, val in meta.items()},
                  name=name or parts[0].name, zero_point=zero,
                  sequence_gap=parts[0].sequence_gap)


# =============================================================================
# Define classes
# =============================================================================
@dataclass
class RVData:
    """
    Radial velocities with their instruments, sequences and indicators

    The velocities are stored with the weighted median of each instrument
    taken out (kept in `zero_point`), so every number a fit sees is of the
    order of the signal and not of the systemic velocity.
    """
    time: np.ndarray
    rv: np.ndarray
    err: np.ndarray
    inst: np.ndarray = None
    seq: np.ndarray = None
    indicators: Dict[str, Tuple[np.ndarray, np.ndarray]] = field(
        default_factory=dict)
    name: str = 'target'
    zero_point: Dict[str, float] = field(default_factory=dict)
    sequence_gap: float = SEQUENCE_GAP
    #: the other columns of the file, one value per point (numbers, or text
    #: where a column is not numeric): LBL copies header keywords there
    #: (EXPTIME, AIRMASS, EXTSN035...), and koloa.outliers reads them
    meta: Dict[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self):
        self.time = np.asarray(self.time, dtype=float)
        self.rv = np.asarray(self.rv, dtype=float)
        self.err = np.asarray(self.err, dtype=float)
        if self.inst is None:
            self.inst = np.array(['inst'] * len(self.time))
        self.inst = np.asarray(self.inst).astype(str)
        # only the points that have everything
        good = (np.isfinite(self.time) & np.isfinite(self.rv)
                & np.isfinite(self.err) & (self.err > 0))
        order = np.argsort(self.time[good], kind='stable')
        keep = np.where(good)[0][order]
        self.time, self.rv = self.time[keep], self.rv[keep]
        self.err, self.inst = self.err[keep], self.inst[keep]
        if self.seq is not None:
            self.seq = np.asarray(self.seq, dtype=int)[keep]
        self.indicators = {key: (np.asarray(val[0], dtype=float)[keep],
                                 np.asarray(val[1], dtype=float)[keep])
                           for key, val in self.indicators.items()}
        self.meta = {key: np.asarray(val)[keep]
                     for key, val in (self.meta or {}).items()}
        # the sequences, found from the time if they are not given
        if self.seq is None:
            self.seq = find_sequences(self.time, self.inst, self.sequence_gap)
        # the zero point of each instrument
        for name in self.instruments:
            mask = self.inst == name
            if name not in self.zero_point:
                self.zero_point[name] = float(np.median(self.rv[mask]))
            self.rv[mask] = self.rv[mask] - self.zero_point[name]

    # -------------------------------------------------------------------------
    @property
    def n(self) -> int:
        """The number of points"""
        return len(self.time)

    @property
    def instruments(self) -> List[str]:
        """The instruments, in order of their first point"""
        _, first = np.unique(self.inst, return_index=True)
        return [str(self.inst[i]) for i in sorted(first)]

    @property
    def inst_index(self) -> np.ndarray:
        """The instrument of each point, as an index into `instruments`"""
        lookup = {name: it for it, name in enumerate(self.instruments)}
        return np.array([lookup[name] for name in self.inst], dtype=int)

    @property
    def nseq(self) -> int:
        """The number of sequences"""
        return int(np.max(self.seq)) + 1 if self.n else 0

    @property
    def baseline(self) -> float:
        """The time span of the observations [days]"""
        return float(np.max(self.time) - np.min(self.time))

    @property
    def tref(self) -> float:
        """The reference time of trends and phases [days]"""
        return float(np.round(np.mean(self.time), 1))

    def summary(self) -> Dict[str, float]:
        """
        The numbers that say what kind of series this is

        :return: dict, the number of points, sequences, the baseline, the
                 median error, the rms and its robust version
        """
        return dict(npoints=self.n, nseq=self.nseq, baseline=self.baseline,
                    ninst=len(self.instruments),
                    median_err=float(np.median(self.err)),
                    rms=float(np.std(self.rv)),
                    robust_rms=robust_std(self.rv))

    # -------------------------------------------------------------------------
    def select(self, mask: np.ndarray) -> 'RVData':
        """
        A copy with only some of the points

        :param mask: np.ndarray, boolean or index array

        :return: RVData, the selection (sequences renumbered)
        """
        mask = np.asarray(mask)
        if mask.dtype == bool:
            mask = np.where(mask)[0]
        indicators = {key: (val[0][mask], val[1][mask])
                      for key, val in self.indicators.items()}
        # the zero point is already out of the velocities: keep it as is,
        #   for the instruments kept only
        kept = set(np.asarray(self.inst)[mask].astype(str))
        zero = {key: val for key, val in self.zero_point.items()
                if key in kept}
        out = RVData.__new__(RVData)
        out.time, out.rv = self.time[mask].copy(), self.rv[mask].copy()
        out.err, out.inst = self.err[mask].copy(), self.inst[mask].copy()
        _, out.seq = np.unique(self.seq[mask], return_inverse=True)
        out.indicators, out.name, out.zero_point = indicators, self.name, zero
        out.sequence_gap = self.sequence_gap
        out.meta = {key: val[mask].copy()
                    for key, val in getattr(self, 'meta', {}).items()}
        return out

    def with_values(self, rv: np.ndarray, err: Optional[np.ndarray] = None
                    ) -> 'RVData':
        """
        A copy with other velocities (an indicator, a residual, a simulation)

        :param rv: np.ndarray, the new values
        :param err: np.ndarray or None, their errors (the same when None)

        :return: RVData, the copy
        """
        out = self.select(np.arange(self.n))
        out.rv = np.asarray(rv, dtype=float).copy()
        if err is not None:
            out.err = np.asarray(err, dtype=float).copy()
        return out

    def indicator(self, name: str) -> 'RVData':
        """
        An activity indicator as a series of its own, ready for any koloa tool

        :param name: str, the indicator

        :return: RVData, the indicator in place of the velocity (median out)
        """
        value, error = self.indicators[name]
        good = np.isfinite(value) & np.isfinite(error) & (error > 0)
        out = self.select(good)
        value, error = value[good], error[good]
        for inst in out.instruments:
            mask = out.inst == inst
            value = value.copy()
            value[mask] = value[mask] - np.median(value[mask])
        out.rv, out.err = value, error
        out.name = f'{self.name} {name}'
        return out

    def binned(self) -> 'RVData':
        """
        One point per sequence: the weighted mean of its exposures

        The error of the mean is the formal one; the scatter inside a
        sequence is not added, because a sequence of a few exposures says
        little about it and the jitter of the fit is there to absorb it.

        :return: RVData, the binned series
        """
        weight = 1.0 / self.err ** 2
        nseq = self.nseq
        wsum = np.bincount(self.seq, weights=weight, minlength=nseq)
        time = np.bincount(self.seq, weights=weight * self.time,
                           minlength=nseq) / wsum
        rv = np.bincount(self.seq, weights=weight * self.rv,
                         minlength=nseq) / wsum
        err = 1.0 / np.sqrt(wsum)
        inst = np.empty(nseq, dtype=object)
        inst[self.seq] = self.inst
        indicators = {}
        for key, (val, sval) in self.indicators.items():
            good = np.isfinite(val) & np.isfinite(sval) & (sval > 0)
            wind = np.where(good, 1.0 / np.where(good, sval, 1) ** 2, 0.0)
            wsum_i = np.bincount(self.seq, weights=wind, minlength=nseq)
            vmean = np.bincount(self.seq, weights=wind * np.where(good, val, 0),
                                minlength=nseq)
            with np.errstate(invalid='ignore', divide='ignore'):
                indicators[key] = (vmean / wsum_i, 1.0 / np.sqrt(wsum_i))
        # the other columns: the mean of each sequence for numbers, the
        #   first exposure's for text
        meta = {}
        _, first = np.unique(self.seq, return_index=True)
        for key, val in getattr(self, 'meta', {}).items():
            if val.dtype.kind in 'fiu':
                good = np.isfinite(val)
                total = np.bincount(self.seq, weights=np.where(good, val, 0.0),
                                    minlength=nseq)
                count = np.bincount(self.seq, weights=good.astype(float),
                                    minlength=nseq)
                with np.errstate(invalid='ignore', divide='ignore'):
                    meta[key] = total / count
            else:
                meta[key] = val[first]
        out = RVData.__new__(RVData)
        out.time, out.rv, out.err = time, rv, err
        out.inst = inst.astype(str)
        out.seq = np.arange(nseq)
        out.indicators, out.name = indicators, self.name + ' (binned)'
        out.zero_point = dict(self.zero_point)
        out.sequence_gap = self.sequence_gap
        out.meta = meta
        return out

    def nightly(self, gap: float = NIGHT_GAP) -> 'RVData':
        """
        One point per night and instrument: the weighted mean of its
        exposures (koloa's analyses run on these by default; nightly=False
        keeps the exposures)

        A night is the exposures of one instrument that no gap longer than
        `gap` separates, so two visits a few hours apart are one point. The
        error of a mean is the formal one, as in binned(): the jitter of the
        fits absorbs the rest. The series keeps its name.

        :param gap: float, the largest gap inside a night [days]

        :return: RVData, one point per night (itself when no night holds two
                 exposures)
        """
        night = night_index(self, gap)
        if len(np.unique(night)) == self.n:
            return self
        tmp = RVData.__new__(RVData)
        tmp.__dict__.update(self.__dict__)
        tmp.seq = night
        out = tmp.binned()
        out.name = self.name
        return out

    # -------------------------------------------------------------------------
    @classmethod
    def from_csv(cls, filename: str, time: Optional[str] = None,
                 rv: Optional[str] = None, err: Optional[str] = None,
                 inst: Optional[str] = None, name: Optional[str] = None,
                 indicators: Optional[List[str]] = None,
                 delimiter: Optional[str] = None,
                 sequence_gap: float = SEQUENCE_GAP) -> 'RVData':
        """
        Read a series from a csv (or LBL .rdb) file

        The columns are found by name when they are not given: rjd, vrad and
        svrad for LBL, and the usual alternatives otherwise. Every other
        column that has an error column next to it becomes an indicator, and
        every column but the time, velocity, error and instrument is kept in
        `meta` (numbers, or text where most of a column is not numeric).

        :param filename: str, the file
        :param time: str or None, the time column
        :param rv: str or None, the velocity column
        :param err: str or None, the error column
        :param inst: str or None, the instrument column; None: a column
                     named inst, instrument or ins_name when it holds
                     names, else one instrument
        :param name: str or None, the name of the target
        :param indicators: list of str or None, which indicators to keep
                           (all of them when None)
        :param delimiter: str or None, the delimiter (guessed when None)
        :param sequence_gap: float, the gap that separates sequences [days]

        :return: RVData, the series
        """
        filename = os.fspath(filename)
        with open(filename, 'r') as handle:
            lines = [line for line in handle if line.strip()
                     and not line.startswith('#')]
        if delimiter is None:
            delimiter = '\t' if '\t' in lines[0] else ','
        reader = csv.reader(lines, delimiter=delimiter)
        rows = list(reader)
        header = [col.strip() for col in rows[0]]
        body = rows[1:]
        # an rdb file has a line of dashes under its header
        if body and all(set(cell.strip()) <= {'-'} for cell in body[0]):
            body = body[1:]
        columns = {}
        for it, col in enumerate(header):
            values = []
            for row in body:
                try:
                    values.append(float(row[it]))
                except (ValueError, IndexError):
                    values.append(np.nan)
            columns[col] = np.array(values)
        texts = {col: np.array([row[it].strip() if it < len(row) else ''
                                for row in body])
                 for it, col in enumerate(header)}
        time = time or _pick(header, TIME_NAMES)
        rv = rv or _pick(header, RV_NAMES)
        err = err or _pick(header, ERR_NAMES)
        if time is None or rv is None or err is None:
            raise ValueError(f'Could not find the time, velocity and error '
                             f'columns in {filename}: {header}')
        # a column named for it gives the instrument of each velocity, when
        #   it holds names (HARPS03, HARPS15: each its own zero point)
        if inst is None:
            inst = next((col for col in header
                         if col.lower() in INST_NAMES and col not in
                         (time, rv, err)
                         and not np.any(np.isfinite(columns[col]))
                         and np.all(texts[col] != '')), None)
        # every column with an error column of its own is an indicator
        used = {time, rv, err}
        found = {}
        for col in header:
            if col in used or col == inst:
                continue
            ecol = _error_column(col, header)
            if ecol is None or ecol in used:
                continue
            if indicators is not None and col not in indicators:
                continue
            found[col] = (columns[col], columns[ecol])
        instv = texts[inst] if inst is not None else None
        # the other columns, as numbers unless most of the filled cells are
        #   not ('nan' and the like count as empty)
        meta = {}
        for col in header:
            if col in used or col == inst:
                continue
            cells = texts[col]
            filled = ~np.isin(np.char.lower(cells), list(EMPTY_CELLS))
            numeric = np.isfinite(columns[col])
            meta[col] = (columns[col] if np.sum(numeric) >= 0.5 * np.sum(filled)
                         else cells)
        if name is None:
            name = filename.split('/')[-1].split('.')[0]
        return cls(time=columns[time], rv=columns[rv], err=columns[err],
                   inst=instv, indicators=found, name=name,
                   sequence_gap=sequence_gap, meta=meta)


# =============================================================================
# End of code
# =============================================================================

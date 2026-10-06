#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The survey tab of koloa's GUI: a sample asked of SIMBAD by its
constraints, what the archives have of each of its stars, the files one
has of them, and the batch of those ticked: run here (the batch tab), or
packed for another machine (koloa.survey.pack: a folder and its tar, with
a script whose first setting is where the folder is there). The results
of a batch brought back are opened as a batch of the page.

Each sample is kept by the server (its stars are many: the page asks for
them once, then for what changed).

Created on 2026-10-06

@author: artigau
"""
import json
import os
import threading
import time
import uuid
from typing import Any, Dict, List, Optional

from koloa import survey

# =============================================================================
# Define variables
# =============================================================================
#: the samples of this session: id -> its stars, its check, its packing
SURVEYS: Dict[str, Dict[str, Any]] = {}
#: what the page shows of a star
SHOWN = ('name', 'main', 'sptype', 'spnum', 'distance', 'ra', 'dec', 'V',
         'G', 'J', 'K', 'near', 'archives')


# =============================================================================
# Define functions
# =============================================================================
def _number(value: Any, default: Optional[float] = None) -> Optional[float]:
    """a number of a field of the page, the default when it is empty"""
    try:
        return float(str(value).replace(',', '.'))
    except (TypeError, ValueError):
        return default


def _star(star: Dict[str, Any]) -> Dict[str, Any]:
    """a star as the page shows it: with the names of its files"""
    out = {key: star.get(key) for key in SHOWN}
    out['files'] = [os.path.basename(path) for path in star.get('files')
                    or []]
    return out


def _picked(job: Dict[str, Any], names: Any) -> List[Dict[str, Any]]:
    """the stars of a sample that are named (all of them for None)"""
    if names is None:
        return list(job['stars'])
    keep = {survey.name_key(name) for name in names}
    return [star for star in job['stars']
            if survey.name_key(star['name']) in keep]


def state(sid: str, stars: bool = True) -> Dict[str, Any]:
    """
    What the page shows of a sample: its stars (unless not asked: they are
    many), how far the check of the archives and the packing are

    :param sid: str, the sample
    :param stars: bool, with its stars
    """
    job = SURVEYS[sid]
    out = dict(id=sid, n=len(job['stars']), asked=job['asked'],
               check=job['check'], pack=job['pack'],
               unmatched=job.get('unmatched') or [])
    if stars:
        out['stars'] = [_star(star) for star in job['stars']]
    return out


def sample(body: Dict[str, Any]) -> Dict[str, Any]:
    """a sample asked of SIMBAD (koloa.survey.sample) from the fields of
    the page: sp_from, sp_to, dmax, dec_min, dec_max, vmax, dwarfs"""
    asked = dict(
        sptype=(str(body.get('sp_from') or 'M0').strip(),
                str(body.get('sp_to') or 'M9').strip()),
        dmax=_number(body.get('dmax'), 15.0),
        dec=(_number(body.get('dec_min'), -90.0),
             _number(body.get('dec_max'), 90.0)),
        vmax=_number(body.get('vmax')),
        dwarfs=body.get('dwarfs', True) is not False)
    stars = survey.sample(**asked)
    sid = uuid.uuid4().hex[:8]
    SURVEYS[sid] = dict(id=sid, stars=stars, asked=asked, check=None,
                        pack=None, made=time.time())
    return state(sid)


def check(body: Dict[str, Any]) -> Dict[str, Any]:
    """what the archives have of the stars of a sample (those named, or
    all), in a thread: DACE, CARMENES DR1, the surveys on VizieR"""
    job = SURVEYS[body['id']]
    if (job['check'] or {}).get('status') == 'running':
        raise ValueError('the archives of this sample are being checked')
    stars = _picked(job, body.get('names'))
    root = str(body.get('root') or 'archives')
    job['check'] = dict(status='running', done=0, total=len(stars),
                        error=None, start=time.time())

    def work():
        try:
            survey.check(stars, root, dace=body.get('dace', True)
                         is not False,
                         progress=lambda done, total: job['check'].update(
                             done=done))
            job['check']['status'] = 'done'
        except Exception as err:  # said on the page
            job['check'].update(status='failed',
                                error=f'{type(err).__name__}: {err}')
        job['check']['end'] = time.time()
    threading.Thread(target=work, daemon=True).start()
    return state(body['id'], stars=False)


def files(body: Dict[str, Any]) -> Dict[str, Any]:
    """the files of folders put with the stars of a sample
    (koloa.survey.match_files): folders (separated by commas, or a list),
    pattern"""
    job = SURVEYS[body['id']]
    folders = body.get('folders') or []
    if isinstance(folders, str):
        folders = [part.strip() for part in folders.replace(';', ',')
                   .split(',') if part.strip()]
    found = survey.match_files(job['stars'], folders,
                               str(body.get('pattern') or '*.rdb'))
    job['unmatched'] = [dict(name=os.path.basename(row['path']),
                             raw=row.get('raw'), target=row.get('target'))
                        for row in found['unmatched']]
    out = state(body['id'])
    out['matched'] = sum(len(val) for val in found['matched'].values())
    return out


def _targets(job: Dict[str, Any], names: Any) -> List[Dict[str, Any]]:
    """the stars ticked, as a batch takes them"""
    stars = _picked(job, names)
    if not stars:
        raise ValueError('no star ticked')
    return stars


def run(body: Dict[str, Any]) -> Dict[str, Any]:
    """the batch of the stars ticked, here: each with its files and every
    archive of it (the batch tab shows it)"""
    from koloa import gui
    stars = _targets(SURVEYS[body['id']], body.get('names'))
    return gui.batch_fip([], body.get('options') or {}, True,
                         bool(body.get('regather')),
                         str(body.get('root') or 'archives'),
                         body.get('rules', True) is not False,
                         targets=[dict(name=star['name'],
                                       files=star.get('files') or [],
                                       sptype=star.get('sptype'))
                                  for star in stars])


def pack(body: Dict[str, Any]) -> Dict[str, Any]:
    """the batch of the stars ticked, packed for another machine
    (koloa.survey.pack), in a thread: name, out, root, server_root, jobs,
    rules, gather, tess"""
    from koloa import gui
    job = SURVEYS[body['id']]
    if (job['pack'] or {}).get('status') == 'running':
        raise ValueError('this sample is being packed')
    stars = _targets(job, body.get('names'))
    name = str(body.get('name') or '').strip()
    if not name:
        raise ValueError('a name for the batch')
    job['pack'] = dict(status='running', done=0, total=len(stars), star='',
                       error=None, result=None, start=time.time())

    def told(done, total, star):
        job['pack'].update(done=done, star=star)

    def work():
        try:
            made = survey.pack(
                stars, name, out=str(body.get('out') or '.'),
                root=str(body.get('root') or 'archives'),
                server_root=str(body.get('server_root') or '').strip()
                or None, gather=body.get('gather', True) is not False,
                tess=bool(body.get('tess')),
                jobs=int(_number(body.get('jobs'), 6)),
                rules=body.get('rules', True) is not False,
                trend=gui.trend_order(body.get('options') or {}),
                progress=told)
            # as the page shows them: from the folder koloa runs in (where
            #   its terminals start too)
            made.update(folder_shown=gui._path_shown(made['folder']),
                        tar_shown=(gui._path_shown(made['tar'])
                                   if made['tar'] else None))
            job['pack'].update(status='done', result=made)
        except Exception as err:  # said on the page
            job['pack'].update(status='failed',
                               error=f'{type(err).__name__}: {err}')
        job['pack']['end'] = time.time()
    threading.Thread(target=work, daemon=True).start()
    return state(body['id'], stars=False)


def results(body: Dict[str, Any]) -> Dict[str, Any]:
    """
    The results of a batch folder (brought back from the machine that ran
    it, or run here by its script) as a batch of the page: its table, each
    star to be opened in the Analysis tab with its FIP as it was computed

    :param body: dict, root: the batch folder
    """
    from koloa import gui
    root = os.path.abspath(os.path.expanduser(str(body.get('root') or '')))
    held = survey.targets_of(root)
    items = []
    for star, res in zip(held['targets'], survey.results(root)):
        if res is None:
            continue
        paths = [os.path.join(root, path) for path in star.get('files') or []]
        qid = None
        quick = os.path.join(os.path.dirname(survey.result_path(
            root, star['name'])), 'quick.json')
        if os.path.exists(quick):
            with open(quick) as handle:
                kept = json.load(handle)
            qid = uuid.uuid4().hex[:8]
            # as a batch's: kept when the page starts a new target
            gui.QUICKS[qid] = dict(kept, id=qid, batch='results',
                                   status=kept.get('status') or 'done',
                                   start=0.0, end=float(kept.get('elapsed')
                                                        or 0.0))
        items.append(dict(
            path=paths[0] if paths else '', files=paths, name=star['name'],
            given=dict(name=star['name'], sptype=star.get('sptype')),
            status=res['status'], qid=qid, error=res.get('error'),
            summary=res.get('summary'), star=res.get('star'), stage=None,
            note=res.get('note'), datasets=res.get('datasets'),
            left_out=[]))
    if not items:
        raise ValueError(f'no star done yet in {root}')
    bid = uuid.uuid4().hex[:8]
    options = held.get('options') or {}
    elapsed = sum(float(res.get('elapsed') or 0.0)
                  for res in survey.results(root) if res)
    gui.BATCHES[bid] = dict(
        id=bid, status='done', start=0.0, end=elapsed, cancel=False,
        trend=int(options.get('trend', 1)), archives=True, regather=False,
        gather=False, rules=bool(options.get('rules', True)),
        root=os.path.join(root, 'archives'), items=items)
    out = gui.batch_state(bid)
    out['name'] = held.get('name')
    out['todo'] = len(held['targets']) - len(items)
    return out


def route(path: str, body: Dict[str, Any]) -> Dict[str, Any]:
    """the questions of the survey tab (/api/survey...)"""
    if path == '/api/survey':
        return state(body['id'], stars=body.get('stars') not in (False, '0',
                                                               0))
    try:
        func = dict(sample=sample, check=check, files=files, run=run,
                    pack=pack, results=results)[path.rsplit('/', 1)[-1]]
    except KeyError:
        raise ValueError(f'no such question: {path}') from None
    return func(body)


# =============================================================================
# End of code
# =============================================================================

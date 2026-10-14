#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LBL for many stars, run elsewhere (a server): a package that fetches the
spectra of each star where it runs, runs LBL on them, keeps the rdb as the
files of the star, adds them to the team's repository of velocities
(koloa.rvstore), and then runs koloa's batch on them (koloa.survey: the
quick FIP of the joint sequences, the detailed report of a candidate).

    koloa "GJ 436" "Gl 408" --lbl-pack HARPSN SOPHIE --name lbl_m --out .

Here, nothing heavy: who each star is, which of its spectra the archives
have (their list and their size: lbl_plan.json, fetch_list.csv), the
archives of its velocities (gathered here, small), koloa and what it
fetched once (the batch folder of koloa.survey.pack), and the script,
run_lbl.py. There: the spectra (GB a star), LBL, koloa's batch.

With fetch (the LBL tab's Fetch the spectra and prepare the package), the
spectra are fetched here and carried in the package (lbl_data/): it runs
on any machine with an LBL, nothing asked of the archives there (but the
models of LBL, the first time).

Created on 2026-10-14

@author: artigau
"""
import csv
import json
import os
import shutil
import tarfile
import time as _time
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from koloa import lblrun
from koloa.log import log, outcome, step

# =============================================================================
# Define variables
# =============================================================================
#: the plan of a package, its list of spectra, what each star gave
PLAN = 'lbl_plan.json'
LIST = 'fetch_list.csv'
DONE = 'results_lbl'
#: what the Python of the batch needs (astropy: the spectra of SOPHIE and
#: of ESO are rewritten for LBL)
NEEDS = ('numpy', 'scipy', 'matplotlib', 'astropy')


# =============================================================================
# Define functions
# =============================================================================
def plan_star(name: str, instruments: Sequence[str],
              root: str = 'archives') -> Dict[str, Any]:
    """
    What a star of a package is: its names, its temperature, and for each
    instrument the runs of LBL and their spectra in the archives (asked,
    not fetched); whether it is among the names never to share

    :return: dict: name, main, names, teff, never, groups (as
             koloa.lblrun.groups makes them, the way their spectra are
             fetched by its name), errors
    """
    from koloa import shared
    from koloa.stars import stellar_teff
    ident = lblrun._ident(name)
    teff = stellar_teff(ident).get('teff') or 4000.0
    names = [ident.get('main') or name] + list(ident.get('aliases') or []) \
        + [name]
    out = dict(name=name, main=ident.get('main') or name, names=names[:12],
               teff=float(teff), never=shared.never(names), groups=[],
               errors={})
    for inst in instruments:
        try:
            for group in lblrun.groups(name, inst, root, ident):
                if group['spectra']:
                    out['groups'].append(dict(lblrun.plain_group(group),
                                              inst=inst))
        except (RuntimeError, ValueError, OSError) as err:
            out['errors'][inst] = str(err)
    return out


def plan(names: Sequence[str], instruments: Sequence[str],
         root: str = 'archives', workers: int = 4) -> List[Dict[str, Any]]:
    """the plans of the stars of a package (plan_star), a few at a time"""
    def one(name):
        try:
            return plan_star(name, instruments, root)
        except Exception as err:  # the star without its runs, and why
            return dict(name=name, main=name, names=[name], teff=4000.0,
                        never=False, groups=[],
                        errors=dict(star=f'{type(err).__name__}: {err}'))
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
        return list(pool.map(one, names))


def size(star: Dict[str, Any]) -> float:
    """about how much a star's spectra weigh [MB]"""
    return sum(len(group['spectra']) * lblrun.INSTRUMENTS[group['inst']][
        'size'] for group in star['groups'])


def fetch_list(stars: Sequence[Dict[str, Any]], path: str):
    """the list of what is fetched: a line per run of a star"""
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['star', 'instrument', 'run', 'spectra', 'MB',
                         'from'])
        for star in stars:
            for group in star['groups']:
                about = lblrun.INSTRUMENTS[group['inst']]
                writer.writerow([star['name'], about['label'],
                                 group['key'] or '-', len(group['spectra']),
                                 round(len(group['spectra']) * about['size']),
                                 about['archive']])


def pack(names: Sequence[str], instruments: Sequence[str], name: str,
         out: str = '.', root: str = 'archives',
         server_root: Optional[str] = None, lbl_dir: str = 'lbl_data',
         lbl_python: str = 'python', nproc: int = 8, push: bool = True,
         analyse: bool = True, fetch: bool = False, jobs: int = 4,
         rules: bool = True, trend: int = 1, report: bool = False,
         tar: bool = True,
         progress: Optional[Callable[[int, int, str], None]] = None
         ) -> Dict[str, Any]:
    """
    A package of LBL for some stars, to run elsewhere (see the module)

    :param names: the stars (names SIMBAD knows)
    :param instruments: the instruments (keys of koloa.lblrun.INSTRUMENTS)
    :param name: str, the name of the package: its folder, its tar
    :param out: str, where the folder is made
    :param root: str, the folder of the archives here
    :param server_root: str or None, where the folder will be where it runs
    :param lbl_dir: str, where the spectra go where it runs (relative: in
                    the folder of the package)
    :param lbl_python: str, the Python that runs LBL where it runs
    :param nproc: int, the processes of LBL's compute there
    :param push: bool, add the rdb to the team's repository there (gh
                 logged in; never a star among the names never to share
                 here)
    :param analyse: bool, run koloa's batch once LBL is done
    :param fetch: bool, fetch the spectra here and carry them (lbl_data/)
    :param jobs, rules, trend, report: those of koloa's batch
                                       (koloa.survey.pack)

    :return: dict: folder, tar, nstar, nspectra, mb, errors (by star)
    """
    from koloa import __version__, survey
    unknown = [inst for inst in instruments
               if inst not in lblrun.INSTRUMENTS]
    if unknown:
        raise ValueError(f'no instrument {unknown} for LBL: '
                         f'{", ".join(lblrun.INSTRUMENTS)}')
    step('the spectra the archives have')
    stars = plan(names, instruments, root)
    nspec = sum(len(group['spectra']) for star in stars
                for group in star['groups'])
    total = sum(size(star) for star in stars)
    outcome(f'{len(stars)} stars, {nspec} spectra, about {total / 1e3:.1f} '
            f'GB')
    step('the batch of koloa (archives, koloa, its cache)')
    made = survey.pack([dict(name=star['name'], main=star['main'])
                        for star in stars], name, out=out, root=root,
                       server_root=server_root, gather=True, tess=False,
                       jobs=jobs, rules=rules, trend=trend, report=report,
                       tar=False, progress=progress)
    folder = made['folder']
    if fetch:
        step('the spectra, fetched here')
        here = os.path.join(folder, 'lbl_data')
        for star in stars:
            for group in star['groups']:
                live = lblrun.live_group(group)
                lblrun.fetch(star['name'], group['inst'], live, here)
                # the template of another star where LBL will want it
                span = lblrun.berv_span(live, here)
                if span is not None and span < lblrun.BERV_MIN:
                    other = lblrun.template_group(live, group['inst'],
                                                  lblrun.TEMPLATE_STAR, root)
                    if other is not None:
                        lblrun.fetch(lblrun.TEMPLATE_STAR, group['inst'],
                                     other, here)
        outcome(f'{nspec} spectra in lbl_data/')
        lbl_dir = 'lbl_data'
    with open(os.path.join(folder, PLAN), 'w') as handle:
        json.dump(dict(name=name, made=_time.strftime('%Y-%m-%d %H:%M'),
                       koloa=__version__, instruments=list(instruments),
                       fetched=bool(fetch), stars=stars), handle, indent=1)
    fetch_list(stars, os.path.join(folder, LIST))
    words = dict(name=name, made=_time.strftime('%Y-%m-%d %H:%M'),
                 nstar=len(stars), nspec=nspec, gb=total / 1e3,
                 root=server_root or folder, lbl_dir=lbl_dir,
                 lbl_python=lbl_python, nproc=int(nproc), push=bool(push),
                 analyse=bool(analyse), needs=NEEDS,
                 fetched=('fetched when the package was made, in lbl_data/'
                          if fetch else 'fetched where it runs: the network'),
                 never=sum(star['never'] for star in stars))
    with open(os.path.join(folder, 'run_lbl.py'), 'w') as handle:
        handle.write(SCRIPT.format(**words))
    section = README.format(**words)
    for file in ('README.txt', 'CLAUDE.md'):
        path = os.path.join(folder, file)
        with open(path) as handle:
            text = handle.read()
        with open(path, 'w') as handle:
            handle.write(section + '\n' + text)
    made_tar = None
    if tar:
        made_tar = folder + '.tar.gz'
        step('the tar')
        with tarfile.open(made_tar, 'w:gz') as handle:
            handle.add(folder, arcname=name)
        outcome(made_tar)
    errors = {star['name']: star['errors'] for star in stars
              if star['errors']}
    return dict(folder=folder, tar=made_tar, nstar=len(stars),
                nspectra=nspec, mb=round(total), errors=errors,
                stars=[dict(name=star['name'], spectra=sum(
                    len(group['spectra']) for group in star['groups']),
                    mb=round(size(star)), never=star['never'],
                    runs=[dict(inst=group['inst'], key=group['key'],
                               n=len(group['spectra']))
                          for group in star['groups']],
                    errors=star['errors']) for star in stars])


def _plan(root: str) -> Dict[str, Any]:
    with open(os.path.join(root, PLAN)) as handle:
        return json.load(handle)


def _mine(stars: List[Dict[str, Any]], part: Optional[Tuple[int, int]],
          only: Sequence[str]) -> List[Dict[str, Any]]:
    """the stars of a part of the batch (the 3rd of 20: 3, 23, 43...)"""
    if only:
        stars = [star for star in stars if star['name'] in only]
    if part is None:
        return stars
    one, of = part
    return [star for it, star in enumerate(stars) if it % of == one - 1]


def check(root: str, lbl_dir: str, python: str) -> bool:
    """
    Whether the package can run here: its Python of LBL, the room for the
    spectra not there yet, gh and the team's repository; said, line by
    line

    :return: bool, whether LBL can run (the room and its Python)
    """
    from koloa import rvstore
    stars = _plan(root)['stars']
    left = 0.0
    for star in stars:
        for group in star['groups']:
            about = lblrun.INSTRUMENTS[group['inst']]
            here = os.path.join(lbl_dir, f'{group["lbl"]}_{group["source"]}',
                                'science', group['obj'])
            there = set(os.listdir(here)) if os.path.isdir(here) else set()
            left += sum(one['name'] not in there
                        for one in group['spectra']) * about['size']
    probe = os.path.abspath(lbl_dir)
    while not os.path.exists(probe):
        probe = os.path.dirname(probe)
    free = shutil.disk_usage(probe).free / 1e6
    state = lblrun.check_python(python)
    gh = rvstore.access()
    done = sum(os.path.exists(os.path.join(root, DONE, f'{_key(star)}.json'))
               for star in stars)
    log(f'{len(stars)} stars, {done} done; the spectra not here yet: about '
        f'{left / 1e3:.1f} GB, {free / 1e3:.0f} GB free in {probe}',
        'value' if free > left * 1.2 else 'warn')
    log(f'LBL: {state["version"]} ({state["where"]}) with {python}'
        if state['ok'] else f'LBL cannot run with {python}: '
                            f'{state["error"]}',
        'value' if state['ok'] else 'warn')
    log(f'the team\'s repository ({rvstore.REPO}): '
        + ('read and written' if gh.get('write') else 'read only'
           if gh.get('read') else 'no access (gh auth login, a token '
                                  'that may write it)'),
        'value' if gh.get('write') else 'warn')
    return bool(state['ok'] and free > left * 1.2)


def _key(star: Dict[str, Any]) -> str:
    from koloa.gather import folder_name
    return folder_name(star['name'])


def run(root: str, lbl_dir: str, python: str, nproc: int = 8,
        push: bool = True, part: Optional[Tuple[int, int]] = None,
        only: Sequence[str] = (), again: bool = False) -> List[str]:
    """
    The stars of a package (a part of them), each in turn: its spectra
    fetched (those not there), LBL on each of its runs, its rdb kept as
    files of the star (archives/<star>/lbl/) and added to the team's
    repository (push; never a star among the names never to share where
    the package was made); what each gave in results_lbl/<star>.json. Then,
    for the whole batch, merge()

    :return: list of str, the stars done here
    """
    stars = _mine(_plan(root)['stars'], part, only)
    os.makedirs(os.path.join(root, DONE), exist_ok=True)
    archives = os.path.join(root, 'archives')
    done = []
    for star in stars:
        mark = os.path.join(root, DONE, f'{_key(star)}.json')
        if os.path.exists(mark) and not again:
            continue
        ident = dict(main=star['main'], aliases=star['names'][1:])
        target = os.path.join(archives, _key(star), 'target.json')
        if os.path.exists(target):
            with open(target) as handle:
                ident = json.load(handle)
        files, errors = [], {}
        for group in star['groups']:
            live = lblrun.live_group(group)
            label = lblrun.INSTRUMENTS[group['inst']]['label'] + (
                f' {group["key"]}' if group['key'] else '')
            try:
                step(f'{star["name"]}, {label}: the spectra')
                got = lblrun.fetch(star['name'], group['inst'], live, lbl_dir)
                outcome(f'{got["n"]} spectra ({got["fetched"]} fetched)')
                if not got['n']:
                    errors[label] = 'no spectrum'
                    continue
                files += lblrun.run_group(
                    star['name'], group['inst'], live, lbl_dir, star['teff'],
                    python, nproc=nproc, root=archives, ident=ident,
                    push=bool(push) and not star['never'])
            except (RuntimeError, ValueError, OSError) as err:
                log(f'{star["name"]}, {label}: {err}', 'warn')
                errors[label] = str(err)
        with open(mark, 'w') as handle:
            json.dump(dict(name=star['name'], errors=errors,
                           made=_time.strftime('%Y-%m-%d %H:%M:%S'),
                           files=[dict(path=os.path.relpath(path, root),
                                       label=label, n=n)
                                  for path, label, n in files]), handle,
                      indent=1)
        done.append(star['name'])
    if part is None:
        merge(root)
    return done


def merge(root: str) -> int:
    """
    What the runs of LBL gave (results_lbl/), as the files of the stars of
    koloa's batch (targets.json: files and labels), for run_batch.py

    :return: int, the stars with files
    """
    from koloa import survey
    path = os.path.join(root, survey.TARGETS)
    with open(path) as handle:
        batch = json.load(handle)
    count = 0
    for star in batch['targets']:
        mark = os.path.join(root, DONE, f'{_key(star)}.json')
        if not os.path.exists(mark):
            continue
        with open(mark) as handle:
            got = json.load(handle)
        star['files'] = [one['path'] for one in got['files']]
        star['labels'] = [one['label'] for one in got['files']]
        count += bool(star['files'])
    with open(path, 'w') as handle:
        json.dump(batch, handle, indent=1)
    return count


SCRIPT = '''#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa: LBL for the batch {name}, made {made}

{nstar} stars: for each, its spectra ({nspec} in all, about {gb:.1f} GB;
{fetched}), LBL run on them with LBL_PYTHON, its rdb kept as files of the
star (archives/<star>/lbl/), added to the team's repository of velocities
when PUSH (gh logged in, a token that may write it), and then koloa's
batch on them when ANALYSE (run_batch.py: the quick FIP of the joint
sequences). fetch_list.csv lists the spectra; README.txt says the rest.

    python run_lbl.py --check       # the Python of LBL, the room for the
                                    #   spectra, gh: can it run here?
    python run_lbl.py               # every star not done yet, then
                                    #   run_batch.py
    python run_lbl.py --part 3/20   # the stars 3, 23, 43... (a job
                                    #   array); once every part is done:
    python run_lbl.py --merge       #   their rdb as the files of the
                                    #   batch, then python run_batch.py
    python run_lbl.py --root /elsewhere/{name} ...
"""
import os
import sys

# =============================================================================
# The settings
# =============================================================================
# where this folder is ON THE MACHINE THAT RUNS IT (or give --root)
ROOT = {root!r}

# where the spectra go: a disk with room ({gb:.1f} GB); relative, in ROOT
LBL_DIR = {lbl_dir!r}

# the Python that runs LBL (an environment with lbl-rapide, or LBL)
LBL_PYTHON = {lbl_python!r}

# the processes of LBL's compute (one star at a time)
NPROC = {nproc}

# add the rdb to the team's repository (koloa-rvs; gh logged in there)
PUSH = {push}

# run koloa's batch on the rdb once LBL is done (run_batch.py)
ANALYSE = {analyse}

# these stars only ([]: all of them), and those already done again
ONLY = []
AGAIN = False


def main():
    root, part = ROOT, None
    args = sys.argv[1:]
    if '--root' in args:
        root = args[args.index('--root') + 1]
    if '--part' in args:
        one, of = args[args.index('--part') + 1].split('/')
        part = (int(one), int(of))
    root = os.path.abspath(os.path.expanduser(root))
    if not os.path.exists(os.path.join(root, 'lbl_plan.json')):
        raise SystemExit('ROOT is ' + root + ': there is no lbl_plan.json '
                         'there. Set ROOT at the top of this script, or '
                         'give --root PATH.')
    lacking = []
    for module in {needs!r}:
        try:
            __import__(module)
        except ImportError:
            lacking.append(module)
    if lacking:
        raise SystemExit('This Python (' + sys.executable + ') has no '
                         + ', '.join(lacking) + ': koloa needs it here '
                         '(LBL runs with LBL_PYTHON).')
    os.environ['KOLOA_CACHE'] = os.path.join(root, 'cache')
    src = os.path.join(root, 'koloa_src')
    if os.path.isdir(src):
        sys.path.insert(0, src)
        os.environ['PYTHONPATH'] = src + os.pathsep + os.environ.get(
            'PYTHONPATH', '')
    from koloa import lblbatch
    lbl_dir = os.path.expanduser(LBL_DIR)
    if not os.path.isabs(lbl_dir):
        lbl_dir = os.path.join(root, lbl_dir)
    if '--check' in args:
        raise SystemExit(0 if lblbatch.check(root, lbl_dir, LBL_PYTHON)
                         else 1)
    if '--merge' in args:
        print(lblbatch.merge(root), 'stars with the rdb of LBL as files')
        return
    lblbatch.run(root, lbl_dir, LBL_PYTHON, NPROC, PUSH, part=part,
                 only=ONLY, again=AGAIN)
    if ANALYSE and part is None:
        import subprocess
        subprocess.run([sys.executable, os.path.join(root, 'run_batch.py'),
                        '--root', root], check=False)


if __name__ == '__main__':
    main()
'''

README = '''LBL first: run_lbl.py
=====================

This batch is a batch of LBL: {nstar} stars, {nspec} spectra (about
{gb:.1f} GB; {fetched}). run_lbl.py fetches the spectra of each star into
LBL_DIR, runs LBL on them (LBL_PYTHON: an environment with lbl-rapide,
or LBL; NPROC processes), keeps the rdb as the files of the star
(archives/<star>/lbl/), adds them to the team's private repository of
velocities when PUSH (github.com/eartigau/koloa-rvs: gh logged in with a
token that may write it; {never} stars never written), and then runs
koloa's batch on them (run_batch.py, below) when ANALYSE.

  1. Set ROOT, LBL_DIR and LBL_PYTHON at the top of run_lbl.py.
  2. python run_lbl.py --check   (the Python of LBL, the room, gh)
  3. python run_lbl.py           (or --part i/N in a job array, then
                                  python run_lbl.py --merge and
                                  python run_batch.py)

fetch_list.csv lists what is fetched, star by star and run by run;
lbl_plan.json holds it all (the archives' names of the spectra). A
spectrum already in LBL_DIR is not fetched again; a star done (in
results_lbl/) is not run again unless AGAIN.

'''

# =============================================================================
# End of code
# =============================================================================

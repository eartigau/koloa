#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run every example of the web page and keep what it printed.

Each docs/examples/ex_<module>.py runs in its own process, from the root of
the repository and on one thread. What it prints (stdout and stderr, in
order, without colour codes, with the paths made relative to the
repository) and how long it takes (wall-clock and CPU seconds) go to
docs/examples/outputs.json, which docs/make_page.py turns into the tabs
of the web page. An example that fails is recorded as failed, and the
others still run.

The title and the blurb of a tab are the first two paragraphs of the
example's docstring, and the caption of its figure the paragraph that
starts with 'Figure:'. The figure is docs/figures/examples/<module>.svg,
kept only when the example wrote it during the run.

    python docs/examples/run_examples.py             # every example
    python docs/examples/run_examples.py fit fip     # only these two

The exit status is 1 when an example failed.

Created on 2026-09-27

@author: artigau
"""
import ast
import glob
import json
import os
import re
import resource
import subprocess
import sys
import time
from typing import Any, Dict, List

from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.dirname(HERE)
ROOT = os.path.dirname(DOCS)
OUTPUT = os.path.join(HERE, 'outputs.json')
FIGDIR = os.path.join(DOCS, 'figures', 'examples')
#: the order of the tabs; an example not listed comes after, alphabetically
ORDER = ['data', 'simulate', 'periodogram', 'fip', 'fit', 'instruments',
         'diagnostics', 'outliers',
         'plotting', 'analyze', 'completeness', 'secular', 'gp', 'kepler',
         'weights', 'radvel_bridge', 'doppler', 'log', 'noise', 'linear']
#: an example still running after this long is stopped [s]
TIMEOUT = 900
#: colour codes of a terminal
ANSI = re.compile(r'\x1b\[[0-9;?]*[A-Za-z]')
#: one thread for the linear algebra, as the page says
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')


# =============================================================================
# Define functions
# =============================================================================
def module_of(path: str) -> str:
    """
    The module an example is about, from its file name

    :param path: str, docs/examples/ex_<module>.py

    :return: str, the module
    """
    return os.path.basename(path)[len('ex_'):-len('.py')]


def describe(path: str) -> Dict[str, Any]:
    """
    The title, blurb, figure caption and code of an example

    The code is the file without its header (the first two comment lines
    and the docstring), which the page shows as the title and the blurb.

    :param path: str, the example

    :return: dict, title, blurb, caption (or None) and code
    """
    with open(path) as handle:
        source = handle.read()
    tree = ast.parse(source)
    doc = ast.get_docstring(tree) or ''
    paragraphs = [' '.join(par.split()) for par in doc.split('\n\n')
                  if par.strip()]
    title = paragraphs[0].rstrip('.') if paragraphs else module_of(path)
    blurb = paragraphs[1] if len(paragraphs) > 1 else ''
    caption = next((par[len('Figure:'):].strip() for par in paragraphs
                    if par.startswith('Figure:')), None)
    if caption:
        caption = caption[0].upper() + caption[1:]
    lines = source.split('\n')
    start = 0
    if doc and tree.body and isinstance(tree.body[0], ast.Expr):
        start = tree.body[0].end_lineno
    code = '\n'.join(lines[start:]).strip('\n') + '\n'
    return dict(title=title, blurb=blurb, caption=caption, code=code)


def clean(text: str) -> str:
    """
    The output as a reader should see it: no colour codes, no trailing
    spaces, and paths relative to the repository

    :param text: str, what the example printed

    :return: str, the cleaned output
    """
    text = ANSI.sub('', text).replace('\r\n', '\n')
    text = text.replace(ROOT + os.sep, '')
    return '\n'.join(line.rstrip() for line in text.split('\n')).strip('\n')


def run(path: str) -> Dict[str, Any]:
    """
    Run one example in its own process and record what it did

    :param path: str, the example

    :return: dict, the entry of outputs.json
    """
    module = module_of(path)
    entry = dict(module=module, script=os.path.relpath(path, DOCS))
    entry.update(describe(path))
    # the examples were written on the exposures: koloa's default of the
    #   nightly means is off for them (KOLOA_NIGHTLY=0)
    env = dict(os.environ, PYTHONUNBUFFERED='1', MPLBACKEND='Agg',
               KOLOA_NIGHTLY='0')
    env.update({name: '1' for name in THREADS})
    figure = os.path.join(FIGDIR, f'{module}.svg')
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    try:
        proc = subprocess.run([sys.executable, path], cwd=ROOT, env=env,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True,
                              timeout=TIMEOUT)
        output, ok = proc.stdout, proc.returncode == 0
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ''
        if isinstance(output, bytes):
            output = output.decode(errors='replace')
        output += f'\n[stopped after {TIMEOUT} s]'
        ok = False
    entry['runtime'] = round(time.time() - start, 2)
    # the CPU time barely depends on what else the machine is doing
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    entry['cpu'] = round(after.ru_utime + after.ru_stime - usage.ru_utime
                         - usage.ru_stime, 2)
    entry['output'] = clean(output)
    entry['ok'] = bool(ok)
    # the figure only if this run wrote it
    fresh = os.path.exists(figure) and os.path.getmtime(figure) >= start
    entry['figure'] = (os.path.relpath(figure, DOCS).replace(os.sep, '/')
                       if ok and fresh else None)
    return entry


def ordered(entries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    The entries in the order of the tabs

    :param entries: list of dict, the entries

    :return: list of dict, sorted
    """
    def key(entry):
        module = entry['module']
        return (ORDER.index(module) if module in ORDER else len(ORDER),
                module)
    return sorted(entries, key=key)


def main(argv: List[str]) -> int:
    """
    Run the examples asked for (all of them by default)

    :param argv: list of str, the modules to run (all when empty)

    :return: int, the exit status (1 when an example failed)
    """
    paths = sorted(glob.glob(os.path.join(HERE, 'ex_*.py')))
    if argv:
        missing = [name for name in argv
                   if not os.path.exists(os.path.join(HERE, f'ex_{name}.py'))]
        if missing:
            log(f'No example for: {", ".join(missing)}', 'error')
            return 1
        paths = [path for path in paths if module_of(path) in argv]
    os.makedirs(FIGDIR, exist_ok=True)
    # a partial run keeps the entries of the examples it does not run
    kept = {}
    if argv and os.path.exists(OUTPUT):
        with open(OUTPUT) as handle:
            kept = {entry['module']: entry for entry in json.load(handle)}
    log(f'Running {len(paths)} example(s) from {ROOT}')
    for path in paths:
        entry = run(path)
        kept[entry['module']] = entry
        if entry['ok']:
            log(f'  {os.path.basename(path):24s} ok in '
                f'{entry["runtime"]:5.1f} s'
                + ('' if entry['figure'] is None else
                   f', figure {entry["figure"]}'), 'value')
        else:
            last = entry['output'].split('\n')[-1] if entry['output'] else ''
            log(f'  {os.path.basename(path):24s} FAILED after '
                f'{entry["runtime"]:.1f} s: {last}', 'warn')
    # entries whose example is gone are dropped
    present = {module_of(path) for path in
               glob.glob(os.path.join(HERE, 'ex_*.py'))}
    entries = ordered([entry for name, entry in kept.items()
                       if name in present])
    with open(OUTPUT, 'w') as handle:
        json.dump(entries, handle, indent=1)
        handle.write('\n')
    failed = [entry['module'] for entry in entries if not entry['ok']]
    log(f'{len(entries) - len(failed)} of {len(entries)} examples ran; '
        f'written to {os.path.relpath(OUTPUT, ROOT)}',
        'warn' if failed else 'info')
    if failed:
        log(f'Failed: {", ".join(failed)}', 'warn')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

# =============================================================================
# End of code
# =============================================================================

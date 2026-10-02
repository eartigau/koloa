#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The summary of a batch of detailed reports.

    koloa --batch-summary BATCH          # or koloa.batch.summary('BATCH')

A batch (the script koloa's GUI writes on its batch page, or any other) puts
the detailed report of each target in BATCH/<target>/, its log in
BATCH/logs/<target>.log, and a copy of its PDF in BATCH/pdf/. summary()
reads every <star>_summary.json, tells which targets failed (no summary:
the end of their log says why), and writes BATCH/batch_summary.pdf (LaTeX,
or a table drawn by matplotlib where there is no pdflatex), .txt and .csv:
for each target its data, its signals (period, K, the FIP of the period or
any of its aliases, the verdict of the duck test, the known planet it is),
the known planets of the star, the acceleration, and its report.

Created on 2026-10-02

@author: artigau
"""
import csv
import glob
import json
import os
import re
import time
from typing import Any, Dict, List, Optional

import numpy as np

from koloa.log import log

# =============================================================================
# Define variables
# =============================================================================
#: the folders of a batch that are not targets
NOT_TARGETS = ('pdf', 'logs')


# =============================================================================
# Define functions
# =============================================================================
def _num(val: Any) -> Optional[float]:
    """a number, None when there is none"""
    try:
        out = float(val)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _row(target: str, summ: Dict[str, Any], folder: str) -> Dict[str, Any]:
    """one target, from its summary"""
    signals = []
    for orb in summ.get('orbits') or []:
        per = _num(orb['P'][0])
        kval = orb.get('K') or [None, None, None]
        signals.append(dict(
            P=per, K=_num(kval[0]),
            K_err=_num(0.5 * (kval[1] + kval[2])) if kval[1] is not None
            else None,
            fip=_num(orb.get('family_fip')), origin=orb.get('origin', ''),
            verdict=(summ.get('duck') or {}).get(f'{per:.4f}', ''),
            known=(orb.get('known') or {}).get('name', '')))
    acc = (summ.get('acceleration') or {}).get('accel')
    pdf = summ.get('report_pdf')
    copy = os.path.join(folder, 'pdf', os.path.basename(pdf)) if pdf else None
    return dict(
        target=target, star=summ.get('star', target), status='done',
        n=summ.get('n'), nvisits=summ.get('nvisits'),
        instruments=summ.get('instruments') or {},
        baseline=_num(summ.get('baseline')), signals=signals,
        known=[pl.get('name') for pl in
               (summ.get('known') or {}).get('planets', [])],
        accel=None if not acc else (_num(acc[0]),
                                    _num(0.5 * (acc[1] + acc[2]))),
        report=(os.path.relpath(copy, folder) if copy and os.path.exists(copy)
                else (os.path.relpath(pdf, folder) if pdf
                      and os.path.exists(pdf) else None)))


def collect(folder: str) -> List[Dict[str, Any]]:
    """
    Every target of a batch: its summary, or why it failed

    :param folder: str, the batch folder

    :return: list of dict, one per target (target, status, and what its
             summary says)
    """
    rows = []
    names = sorted(name for name in os.listdir(folder)
                   if os.path.isdir(os.path.join(folder, name))
                   and name not in NOT_TARGETS and not name.startswith('.'))
    logs = {os.path.splitext(os.path.basename(path))[0]: path for path in
            glob.glob(os.path.join(folder, 'logs', '*.log'))}
    for name in sorted(set(names) | set(logs)):
        found = sorted(glob.glob(os.path.join(folder, name,
                                              '*_summary.json')))
        if found:
            with open(found[0]) as handle:
                rows.append(_row(name, json.load(handle), folder))
            continue
        why = 'no report (not run yet, or still running)'
        if name in logs:
            with open(logs[name], errors='replace') as handle:
                tail = [line.strip() for line in handle if line.strip()]
            # the error itself when there is one, else the end of the log
            errors = [line for line in tail
                      if re.search(r'(Error|Exception)\b.*:', line)]
            why = (errors[-1] if errors else ' / '.join(tail[-2:]))[-400:] \
                or 'an empty log'
        rows.append(dict(target=name, star=name, status='failed', why=why))
    return rows


def _signals_text(row: Dict[str, Any]) -> str:
    """the signals of a target in words"""
    words = []
    for sig in row.get('signals', []):
        text = f'{sig["P"]:.4f} d'
        if sig['K'] is not None:
            text += f' K {sig["K"]:.2f}'
            if sig['K_err'] is not None:
                text += f'+-{sig["K_err"]:.2f}'
        if sig['fip'] is not None:
            text += f' FIP {sig["fip"]:.1e}'
        for key in ('verdict', 'known'):
            if sig.get(key):
                text += f' [{sig[key]}]'
        words.append(text)
    return '; '.join(words) or 'none'


def _write_text(rows: List[Dict[str, Any]], folder: str) -> str:
    """the summary as text"""
    path = os.path.join(folder, 'batch_summary.txt')
    lines = [f'koloa batch summary, {folder}, {time.strftime("%Y-%m-%d %H:%M")}',
             f'{len(rows)} targets, '
             f'{sum(row["status"] == "done" for row in rows)} done', '']
    for row in rows:
        if row['status'] != 'done':
            lines.append(f'{row["target"]}: FAILED: {row["why"]}')
            continue
        insts = ', '.join(f'{key} {val}' for key, val in
                          row['instruments'].items())
        acc = (f'; dv/dt {row["accel"][0]:+.2f}+-{row["accel"][1]:.2f} '
               f'm/s/yr' if row['accel'] and row['accel'][0] is not None
               else '')
        lines.append(f'{row["star"]}: {row["n"]} points ({insts}) over '
                     f'{row["baseline"]:.0f} d{acc}')
        lines.append(f'  signals: {_signals_text(row)}')
        if row['known']:
            lines.append(f'  known planets: {", ".join(row["known"])}')
        if row['report']:
            lines.append(f'  report: {row["report"]}')
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    return path


def _write_csv(rows: List[Dict[str, Any]], folder: str) -> str:
    """the summary as a table, one line per signal (or per target)"""
    path = os.path.join(folder, 'batch_summary.csv')
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['target', 'star', 'status', 'n', 'instruments',
                         'baseline_d', 'accel_mps_yr', 'accel_err',
                         'P_d', 'K_mps', 'K_err', 'fip_period_or_alias',
                         'origin', 'verdict', 'known', 'report', 'why'])
        for row in rows:
            base = [row['target'], row.get('star'), row['status'],
                    row.get('n'),
                    ' '.join(f'{key}:{val}' for key, val in
                             (row.get('instruments') or {}).items()),
                    row.get('baseline'),
                    *(row.get('accel') or (None, None))]
            sigs = row.get('signals') or [None]
            for sig in sigs:
                sig = sig or {}
                writer.writerow(base + [sig.get('P'), sig.get('K'),
                                        sig.get('K_err'), sig.get('fip'),
                                        sig.get('origin'), sig.get('verdict'),
                                        sig.get('known'), row.get('report'),
                                        row.get('why', '')])
    return path


def _write_latex(rows: List[Dict[str, Any]], folder: str) -> Optional[str]:
    """the summary as a LaTeX table, compiled (None without pdflatex)"""
    from koloa.latex import _path, _preamble, compile_pdf, escape
    pre = _preamble('batch').replace('detailed analysis of batch',
                                     'batch summary')
    done = sum(row['status'] == 'done' for row in rows)
    out = [pre, '\\section*{koloa: batch summary}',
           f'{_path(os.path.abspath(folder))}, '
           f'{time.strftime("%Y-%m-%d %H:%M")}: {len(rows)} targets, '
           f'{done} done' + (f', \\status{{flag}}{{{len(rows) - done} '
                             f'failed}}' if done < len(rows) else '') + '.',
           '',
           'For each target: its data (points, instruments, span), its '
           'signals (period, K in m/s, the FIP of the period or any of its '
           'aliases, the verdict of the duck test, the known planet it is), '
           'the acceleration of the star, and its detailed report.\n',
           '\\small',
           '\\begin{longtable}{@{}p{2.6cm}p{3.6cm}p{6.6cm}p{2.0cm}p{1.0cm}@{}}',
           '\\toprule Target & Data & Signals & dv/dt [m/s/yr] & Report '
           '\\\\\\midrule\\endhead']
    for row in rows:
        name = escape(row['star'] if row['status'] == 'done'
                      else row['target'])
        if row['status'] != 'done':
            out.append(f'{name} & \\multicolumn{{4}}{{p{{13.2cm}}}}'
                       f'{{\\status{{flag}}{{failed:}} '
                       f'{escape(row["why"])}}} \\\\[2pt]')
            continue
        insts = ', '.join(f'{escape(key)} {val}' for key, val in
                          row['instruments'].items())
        data = f'{row["n"]} points over {row["baseline"]:.0f}\\,d: {insts}'
        sigs = []
        for sig in row['signals']:
            text = f'{sig["P"]:.4f}\\,d'
            if sig['K'] is not None:
                text += f', K {sig["K"]:.2f}'
                if sig['K_err'] is not None:
                    text += f'$\\pm${sig["K_err"]:.2f}'
            if sig['fip'] is not None:
                fiptext = f'{sig["fip"]:.1e}'
                text += f', FIP {escape(fiptext)}'
            if sig.get('verdict'):
                text += f', \\textbf{{{escape(sig["verdict"])}}}'
            if sig.get('known'):
                text += f' ({escape(sig["known"])})'
            sigs.append(text)
        if row['known']:
            sigs.append('\\textcolor{muted}{known: '
                        + escape(', '.join(row['known'])) + '}')
        acc = ('--' if not row['accel'] or row['accel'][0] is None else
               f'${row["accel"][0]:+.2f}\\pm{row["accel"][1]:.2f}$')
        link = (f'\\href{{{row["report"]}}}{{PDF}}' if row['report']
                else '--')
        out.append(f'{name} & {data} & ' + ('\\newline '.join(sigs) or
                                            'none')
                   + f' & {acc} & {link} \\\\[2pt]')
    out += ['\\bottomrule', '\\end{longtable}', '\\end{document}']
    tex = os.path.join(folder, 'batch_summary.tex')
    with open(tex, 'w') as handle:
        handle.write('\n'.join(out) + '\n')
    return compile_pdf(tex)


def _write_mpl(rows: List[Dict[str, Any]], folder: str) -> str:
    """the summary drawn by matplotlib (no pdflatex): the text, by page"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    path = os.path.join(folder, 'batch_summary.pdf')
    with open(_write_text(rows, folder)) as handle:
        lines = handle.read().splitlines()
    per_page = 60
    with PdfPages(path) as pdf:
        for start in range(0, max(len(lines), 1), per_page):
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.05, 0.97, '\n'.join(lines[start:start + per_page]),
                     va='top', ha='left', family='monospace', fontsize=6.5)
            pdf.savefig(fig)
            plt.close(fig)
    return path


def summary(folder: str) -> Dict[str, Optional[str]]:
    """
    The summary of a batch, in BATCH/batch_summary.pdf, .txt and .csv

    :param folder: str, the batch folder

    :return: dict, the paths written (pdf, txt, csv)
    """
    rows = collect(folder)
    out = dict(txt=_write_text(rows, folder), csv=_write_csv(rows, folder))
    pdf = None
    try:
        pdf = _write_latex(rows, folder)
    except Exception as err:  # the table is a help, not a stop
        log(f'batch summary in LaTeX: {err}', 'warn')
    out['pdf'] = pdf or _write_mpl(rows, folder)
    done = sum(row['status'] == 'done' for row in rows)
    log(f'batch summary: {len(rows)} targets, {done} done, in '
        f'{out["pdf"]}', 'value')
    return out


def script(targets: List[Dict[str, Any]], options: Dict[str, Any],
           koloa: str = 'koloa', batch: str = '', workdir: str = '',
           jobs: int = 1, bashrc: bool = True, where: str = 'this machine'
           ) -> Dict[str, Any]:
    """
    The bash script of a batch: the detailed report of every target, its
    log, a copy of its PDF, then the summary (koloa --batch-summary)

    Every target's command line is the one koloa's GUI would run for it
    (koloa.gui.command), with its own folder in the batch.

    :param targets: list of dict, target (a SIMBAD name, may be empty) and
                    files (list of dict, path and label; may be empty)
    :param options: dict, the options of the detailed report, for all
    :param koloa: str, how koloa is called where the script runs (koloa,
                  /path/python -m koloa.cli, PYTHONPATH=... before it)
    :param batch: str, the batch folder (batch_<date> when empty)
    :param workdir: str, the folder the script runs in (where it is)
    :param jobs: int, the targets at a time
    :param bashrc: bool, source ~/.bashrc first (DACE_API_KEY, conda...)
    :param where: str, the machine, for the header

    :return: dict, script (the text), name (its file name), targets (the
             folder of each), warnings
    """
    import shlex
    from koloa.gather import folder_name
    from koloa.gui import command
    batch = batch.strip() or time.strftime('batch_%Y%m%d_%H%M')
    words = shlex.split(koloa.strip() or 'koloa')
    envs = [word for word in words if '=' in word.split('/')[0]
            and not word.startswith('-')]
    words = words[len(envs):] if words[:len(envs)] == envs else words
    lines, folders, warnings = [], [], []
    for item in targets:
        name = str(item.get('target') or '').strip()
        files = [f for f in item.get('files') or []
                 if str(f.get('path') or '').strip()]
        if not name and not files:
            continue
        stem = name or os.path.splitext(os.path.basename(
            files[0]['path']))[0]
        folder = folder_name(stem) or 'target'
        while folder in folders:
            folder += '_2'
        folders.append(folder)
        args = command('detailed', dict(options, target=name, files=files,
                                        outdir='OUTDIR'))
        cut = args.index('--outdir')
        args = args[:cut] + args[cut + 2:]
        label = name or os.path.basename(files[0]['path'])
        lines.append(f'slot; one {shlex.quote(label)} {shlex.quote(folder)} '
                     + ' '.join(shlex.quote(arg) for arg in args) + ' &')
    if not lines:
        raise ValueError('no target: give each one a SIMBAD name, a file, or '
                         'both')
    name = f'koloa_{batch}.sh'
    run_line = (f'nohup bash {name} > {os.path.splitext(name)[0]}.out 2>&1 &')
    out = ['#!/usr/bin/env bash',
           f'# koloa batch: {len(lines)} targets, written by koloa\'s GUI on '
           f'{time.strftime("%Y-%m-%d %H:%M")}, to run on {where}:',
           f'#   cd {workdir or "<this folder>"} && {run_line}',
           '# Each target: its detailed report in BATCH/<target>/, its log in',
           '# BATCH/logs/<target>.log, its PDF copied to BATCH/pdf/; at the '
           'end, the',
           '# summary of the batch, BATCH/batch_summary.pdf (.txt, .csv).', '']
    if bashrc:
        out += ['# the environment of your shell (a script does not read it '
                'otherwise):',
                '#   DACE_API_KEY, conda, PATH...',
                '[ -f "$HOME/.bashrc" ] && source "$HOME/.bashrc"', '']
    out += ['set -u']
    if workdir.strip():
        out += [f'cd {shlex.quote(workdir.strip())} || exit 1']
    out += [f'export {env}' for env in envs]
    out += [f'KOLOA=({" ".join(shlex.quote(word) for word in words)})'
            '    # how koloa is called here',
            f'BATCH={shlex.quote(batch)}',
            f'JOBS={max(int(jobs or 1), 1)}    # targets at a time',
            'mkdir -p "$BATCH/logs" "$BATCH/pdf"', '',
            '# a line as koloa writes them (koloa\'s GUI reads its steps)',
            'say() { echo "$(date \'+%y%m%d %H:%M:%S\').00 | $*"; }', '',
            '# one target: its report, its log, a copy of its PDF',
            'one() {',
            '  local name=$1 folder=$2',
            '  shift 2',
            '  say "step: $name"',
            '  if "${KOLOA[@]}" "$@" --outdir "$BATCH/$folder" '
            '> "$BATCH/logs/$folder.log" 2>&1; then',
            '    cp "$BATCH/$folder"/*_report.pdf "$BATCH/pdf/" 2>/dev/null',
            '    say "$name: done"',
            '  else',
            '    say "$name: FAILED (its log: $BATCH/logs/$folder.log)"',
            '  fi',
            '}', '',
            '# wait for a free slot',
            'slot() { while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do '
            'sleep 5; done; }', '',
            f'say "koloa batch: {len(lines)} targets in $BATCH, $JOBS at a '
            f'time"']
    out += lines
    out += ['wait', 'say "step: the summary of the batch"',
            '"${KOLOA[@]}" --batch-summary "$BATCH"',
            'say "done: $BATCH/batch_summary.pdf"', '']
    return dict(script='\n'.join(out), name=name, batch=batch,
                targets=folders, warnings=warnings, start=run_line)


# =============================================================================
# End of code
# =============================================================================

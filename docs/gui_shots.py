#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The screenshots of koloa's GUI for the site, and the numbers they show.

Public data only: the server runs without a DACE key (DACE_API_KEY empty:
DACE answers anonymously), and its remembered targets are kept in the work
folder, not the user's. GJ 436 is gathered from its name (DACE, CARMENES
DR1, VizieR, TESS); five stars with public HARPS velocities (GJ 876, GJ 581,
HD 69830, GJ 667 C, HD 40307) become files of velocities, each with its
OBJECT, for the batch and for the star found from a file.

    python gui_shots.py [--work docs/gui_demo] [--port 8790] [--skip-prep]
                        [--langs en,fr] [--parts name,file,runs,batch]
                        [--server http://127.0.0.1:8790]

writes docs/figures/gui/*.png (English), docs/figures/gui/fr/*.png (the
GUI in French) and docs/figures/gui/numbers.json (what the walk-throughs
quote; make_page.py fills the page with them). Needs playwright and Chrome
(python -m playwright install chrome), and the network.

Created on 2026-10-04

@author: artigau
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'figures', 'gui')
#: the stars whose public HARPS velocities become files
FILE_STARS = ('GJ 876', 'GJ 581', 'HD 69830', 'GJ 667 C', 'HD 40307')
#: the width of the browser [px]
WIDTH = 1360
#: how long a quick FIP may take [s] (GJ 436 with all its archives, about
#:   ten instruments: some 15 minutes)
FIP_LIMIT = 3600.0


# =============================================================================
# The data: public only
# =============================================================================
def prepare(work: str) -> None:
    """GJ 436's archives, and the files of the batch (public HARPS)"""
    os.environ['DACE_API_KEY'] = ''
    from koloa.gather import folder_name, gather, load
    cwd = os.getcwd()
    os.chdir(work)
    try:
        gather('GJ 436', 'archives', api_key=False)
        os.makedirs('files', exist_ok=True)
        for star in FILE_STARS:
            gather(star, 'archives_files', dace=True, carmenes=False,
                   tess=False, vizier=False, api_key=False)
            folder = os.path.join('archives_files', folder_name(star))
            if not os.path.exists(os.path.join(folder, 'rv', 'all_rv.csv')):
                continue
            rv = load(folder, photometry=False)['rv']
            sel = np.char.startswith(rv.inst.astype(str), 'HARPS')
            if sel.sum() < 20:
                continue
            zero = np.array([rv.zero_point[str(inst)]
                             for inst in rv.inst[sel]])
            name = folder_name(star).lower() + '_harps.csv'
            with open(os.path.join('files', name), 'w') as handle:
                handle.write('rjd,vrad,svrad,inst,OBJECT\n')
                for row in zip(rv.time[sel], rv.rv[sel] + zero, rv.err[sel],
                               rv.inst[sel]):
                    handle.write(f'{row[0]:.6f},{row[1]:.3f},{row[2]:.3f},'
                                 f'{row[3]},{star}\n')
    finally:
        os.chdir(cwd)


# =============================================================================
# The server and the browser
# =============================================================================
def start_server(work: str, port: int, remembered: str, koloa: str = ROOT):
    """koloa's GUI in the work folder: no DACE key, its own remembered
    targets"""
    # a home of its own (no DACE key in it, a plain prompt in the terminal
    #   of the page, its own routes), what koloa fetched once as it is here
    home = os.path.join(work, 'home')
    os.makedirs(home, exist_ok=True)
    with open(os.path.join(home, '.zshrc'), 'w') as handle:
        handle.write("PS1='%1~ %# '\n")
    env = dict(os.environ, DACE_API_KEY='', PYTHONPATH=koloa, HOME=home,
               ZDOTDIR=home, SHELL='/bin/zsh',
               KOLOA_CACHE=os.path.join(os.path.expanduser('~'), '.cache',
                                        'koloa'))
    # the key of its terminal, known here (the page is opened with it),
    #   and its routes in the work folder (not the user's)
    import secrets
    key = secrets.token_urlsafe(12)
    code = ('from koloa import gui, terminal; '
            f'gui.REMEMBERED = {remembered!r}; terminal.KEY = {key!r}; '
            f'terminal.ROUTES = {os.path.join(work, "routes_demo.json")!r}; '
            f'gui.serve(port={port}, browser=False)')
    proc = subprocess.Popen([sys.executable, '-c', code], cwd=work, env=env,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    proc.koloa_key = key
    time.sleep(8)
    return proc


class Shots:
    """the screenshots of one language"""

    def __init__(self, page, lang: str):
        self.page = page
        self.lang = lang
        self.folder = OUT if lang == 'en' else os.path.join(OUT, lang)
        os.makedirs(self.folder, exist_ok=True)

    def path(self, name: str) -> str:
        return os.path.join(self.folder, name + '.png')

    def plain(self) -> None:
        """the number fields as text in English: Chrome on macOS writes
        their decimals as the system does (0,470 on a French one)"""
        if self.lang == 'en':
            self.page.evaluate(
                "() => document.querySelectorAll('input[type=number]')"
                ".forEach((box) => { box.type = 'text'; })")

    def element(self, selector: str, name: str) -> None:
        self.plain()
        loc = self.page.locator(selector).first
        loc.scroll_into_view_if_needed()
        time.sleep(0.6)
        loc.screenshot(path=self.path(name))
        print('  ', name)

    def region(self, top: str, bottom: str, name: str, pad: int = 8,
               inside: str = '#fipcard') -> None:
        """from the top of one element to the bottom of another, as wide as
        the card they are in"""
        self.plain()
        page = self.page
        page.locator(top).first.scroll_into_view_if_needed()
        time.sleep(0.6)
        box_top = page.locator(top).first.bounding_box()
        box_bot = page.locator(bottom).first.bounding_box()
        card = page.locator(inside).first.bounding_box()
        scroll = page.evaluate('window.scrollY')
        clip = dict(x=card['x'], y=box_top['y'] + scroll - pad,
                    width=card['width'],
                    height=box_bot['y'] + box_bot['height'] - box_top['y']
                    + 2 * pad)
        page.screenshot(path=self.path(name), clip=clip, full_page=True)
        print('  ', name)


def wait(page, condition: str, limit: float = 900.0, retry: str = '',
         every: float = 120.0) -> None:
    """until a condition of the page holds, for so many seconds awake (a
    laptop asleep, its lid closed, does not use them up: the seconds are
    counted a step at a time, not read from the clock)

    :param retry: str, JavaScript run again every so many seconds while
                  nothing came: what asks for it once more (a request sent
                  as the laptop fell asleep, or as the tunnel dropped, is
                  lost for good)
    """
    spent = 0.0
    while spent < limit:
        if page.evaluate(f'() => !!({condition})'):
            return
        time.sleep(1.0)
        spent += 1.0
        if retry and spent % every == 0:
            try:
                page.evaluate(f'() => {{ {retry} }}')
            except Exception as err:  # asked again the next time
                print('   (asked again:', str(err).splitlines()[0][:80], ')')
    raise TimeoutError(condition)


def settle(page, selector: str = '#tsnote', limit: float = 600.0) -> None:
    """until a note has no spinner (seconds awake, as wait counts them)"""
    spent = 0.0
    while spent < limit:
        if 'spin' not in page.inner_html(selector):
            return
        time.sleep(1.0)
        spent += 1.0


# =============================================================================
# The walk-throughs
# =============================================================================
#: the Keplerian orbit of the fold shown, asked again
KEPLER = 'if (foldShown !== null) showFold(foldShown);'


def plotted(page, shots: Shots) -> None:
    """GJ 436 from its name alone, as far as its velocities plotted: the
    star, its archives, the series and the table of its datasets (no
    computation: taken again after a change of the page)"""
    page.fill('#root', 'archives')
    page.fill('#target', 'GJ 436')
    page.press('#target', 'Enter')
    wait(page, "document.querySelector('#ident .stats-grid')",
         retry='resolveStar();')
    time.sleep(2)
    shots.element('#starcard', 'star')
    shots.element('section:has(#run-gather)', 'gather')
    for key in ('dace', 'carmenes', 'vizier'):
        page.check(f'[data-mirror="{key}"]')
    page.click('#plot')
    wait(page, "document.querySelector('#rvplot.on') && lastRV",
         retry="document.getElementById('plot').click();")
    time.sleep(2)
    # the pointer off the plot: its tools show under it
    page.mouse.move(5, 5)
    time.sleep(0.5)
    shots.element('section:has(#filelist)', 'velocities')
    shots.element('#plotcard', 'series')


def by_name(page, shots: Shots, numbers: dict) -> None:
    """GJ 436 from its name alone: the star, its archives, the quick FIP,
    the folds, the transit in TESS, the report, a result remembered"""
    plotted(page, shots)
    page.click('.startfip')
    wait(page, "quick && quick.status === 'done'", FIP_LIMIT,
         retry="if (!quick) startQuick(); else if (quick.status === "
               "'running') { quick.misses = 0; pollQuick(quick.id); }")
    settle(page)
    time.sleep(3)
    shots.region('#fipstatus', '#fipplot', 'fip')
    page.click('[data-fipview="each"]')
    time.sleep(1.5)
    shots.region('#fipstatus', '#fipeach', 'fip_each')
    page.click('[data-fipview="joint"]')
    shots.region('.foldhead', '#foldplot', 'fold')
    page.click('[data-fcol="date"]')
    time.sleep(1.5)
    shots.region('.foldhead', '#foldplot', 'fold_date')
    page.click('[data-fcol="inst"]')
    page.click('[data-fmodel="kepler"]')
    wait(page, "currentFold() && currentFold().model && "
               "currentFold().model.kind === 'kepler'", 1200, retry=KEPLER)
    time.sleep(2)
    shots.region('.foldhead', '#foldplot', 'fold_kepler')
    page.check('#foldoverlay')
    time.sleep(4)
    shots.element('#plotcard', 'series_solution')
    numbers['kepler'] = page.evaluate(
        '() => { const f = currentFold(); return { P: f.period, '
        'P_err: f.P_err, K: f.K, K_err: f.K_err, e: f.e, e_err: f.e_err }; }')
    # zoomed in on the fortnight with the most nights: the orbit itself
    span = page.evaluate(
        '() => { const xs = [...new Set(lastRV.flatMap((inst) => '
        'inst.time.map(Math.floor)))].sort((a, b) => a - b); '
        'let best = [0, 0]; '
        'for (let i = 0, j = 0; i < xs.length; i++) { '
        'while (xs[i] - xs[j] > 14) j++; '
        'if (i - j > best[1] - best[0]) best = [j, i]; } '
        'return [xs[best[0]] - 1, xs[best[1]] + 1]; }')
    page.evaluate('([a, b]) => Plotly.relayout("rvplot", '
                  '{"xaxis.range[0]": a, "xaxis.range[1]": b})', span)
    time.sleep(5)
    shots.element('#plotcard', 'series_zoom')
    page.evaluate('() => Plotly.relayout("rvplot", {"xaxis.autorange": true})')
    time.sleep(2)
    page.uncheck('#foldoverlay')
    page.click('[data-fmodel="sine"]')
    time.sleep(2)
    numbers['quick'] = page.evaluate(
        '() => { const r = quick.result; const pk = r.peak_list[0]; '
        'const f = r.folds.find((x) => x.id === pk.id); return { '
        'n: r.n, nexp: r.nexp, instruments: r.instruments, '
        'period: pk.period, fip: pk.family, K: f.K, K_err: f.K_err, '
        'acceleration: r.acceleration, known: r.known.map((p) => p.name), '
        'transits: (r.transits || []).map((t) => t.name), '
        'msini: minimumMass(f), star: starNow, peaks: r.peak_list.map('
        '(p) => ({ id: p.id, period: p.period, fip: p.family })) }; }')
    if page.locator('[data-known]').count():
        page.locator('[data-known]').first.click()
        wait(page, "currentFold() && currentFold().known", 300, every=60,
             retry="const b = document.querySelector('[data-known]'); "
                   "if (b) b.click();")
        time.sleep(4)
        shots.region('.foldhead', '#foldplot', 'fold_known')
    if page.locator('#foldbuttons [data-transit]').count():
        page.locator('#foldbuttons [data-transit]').first.click()
        wait(page, "currentFold() && currentFold().transit", 300, every=60,
             retry="const b = document.querySelector("
                   "'#foldbuttons [data-transit]'); if (b) b.click();")
        time.sleep(4)
        shots.region('.foldhead', '#foldplot', 'fold_transit')
    # the transit in TESS, at peak #1 (about the conjunction of its fold)
    page.locator('[data-tsearch]').first.click()
    time.sleep(1)
    settle(page)
    wait(page, "tsLast && tsLast.res && !document.querySelector("
               "'#tsnote .bad') && !document.querySelector('#tsnote .spin')",
         900, retry="showTransit(transitCandidates()[0].key, true);")
    time.sleep(2)
    shots.region('#tsbuttons', '#tsplot', 'tess_fold')
    numbers['tess'] = page.evaluate(
        '() => { const r = tsLast.res; return { plausible: r.plausible, '
        'why: r.why, fit: r.fit, best: r.best, sectors: r.sectors, '
        'fold_period: r.fold_period, scan: r.scan || null, '
        'window: r.window, star: r.star, duration: r.duration }; }')
    page.click('[data-tsview="series"]')
    time.sleep(3)
    shots.region('#tsbuttons', '#tsplot', 'tess_series')
    page.click('[data-tsview="fold"]')
    time.sleep(2)
    # the report and its options, the analysis script
    shots.element('section:has(#run-detailed)', 'report')
    # a help bubble
    page.locator('[data-help="transit_search"]').hover()
    time.sleep(0.8)
    box = page.locator('#tip').bounding_box()
    scroll = page.evaluate('window.scrollY')
    page.screenshot(path=shots.path('help'), full_page=True, clip=dict(
        x=max(0, box['x'] - 30), y=box['y'] + scroll - 40,
        width=min(WIDTH - box['x'] + 30, box['width'] + 60),
        height=box['height'] + 60))
    page.mouse.move(5, 5)
    # the Keplerian orbit of peak #1 subtracted (a sinusoid would leave
    #   the harmonic of an eccentric orbit), the FIP of what is left
    page.locator('[data-fold]').first.click()
    time.sleep(2)
    page.click('[data-fmodel="kepler"]')
    wait(page, "currentFold() && currentFold().model && "
               "currentFold().model.kind === 'kepler'", 1200, retry=KEPLER)
    time.sleep(2)
    # ticked: its orbit taken out, the FIP of what is left drawn over the
    #   FIP of the series, its solution on the series
    page.locator('[data-tick]').first.check()
    wait(page, "stageLists().length === 1 && lastStage() && "
               "lastStage().job.status === 'done'", FIP_LIMIT,
         retry="if (!ticked.length) document.querySelector('[data-tick]')"
               ".click(); else if (!tickBusy) { for (const [key, st] of "
               "stages) { if (st.error) stages.delete(key); } "
               "updateStages(); }")
    time.sleep(3)
    shots.region('#fipstatus', '#fipplot', 'fip_residuals')
    numbers['residuals'] = page.evaluate(
        '() => { const pk = lastStage().job.result.peaks[0]; '
        'return { period: pk.period, fip: pk.family }; }')
    time.sleep(1)
    # remembered
    page.fill('#remnote', 'GJ 436 b: the hot Neptune that transits'
              if shots.lang == 'en' else
              'GJ 436 b : le Neptune chaud qui transite')
    page.click('#remember')
    wait(page, "document.querySelector('#remstate').textContent"
               ".includes('\u2713')", 600,
         retry="if (!document.querySelector('#remstate .spin')) "
               "rememberResult();")
    page.click('[data-tab="remembered"]')
    time.sleep(1.5)
    shots.element('#tab-remembered', 'remembered')
    page.click('[data-tab="analysis"]')


def recalled(page, shots: Shots) -> None:
    """the result the walk-through remembered, recalled: what needs no
    computation taken again (after a change of the page): the Keplerian
    fold of peak #1, unticked"""
    page.click('[data-tab="remembered"]')
    wait(page, "document.querySelector('[data-recall]')", 120)
    page.locator('[data-recall]').first.click()
    wait(page, "quick && quick.result && "
               "document.querySelector('#fipplot.on')", 300)
    page.click('[data-tab="analysis"]')
    time.sleep(3)
    if page.locator('[data-tick]:checked').count():
        # by its own click: the buttons are drawn again as it is unticked
        page.evaluate("() => document.querySelector("
                      "'[data-tick]:checked').click()")
        wait(page, "!ticked.length && !tickBusy", 120)
        time.sleep(2)
    page.locator('[data-fold]').first.click()
    time.sleep(2)
    page.click('[data-fmodel="kepler"]')
    wait(page, "currentFold() && currentFold().model && "
               "currentFold().model.kind === 'kepler'", 1200, retry=KEPLER)
    time.sleep(3)
    shots.region('.foldhead', '#foldplot', 'fold_kepler')
    # the card that gathers the archives (its folder as one types it, not
    #   the path the remembered result holds), the transit panel, its help
    page.fill('#root', 'archives')
    time.sleep(2)
    shots.element('section:has(#run-gather)', 'gather')
    page.locator('[data-tsearch]').first.click()
    time.sleep(1)
    settle(page)
    wait(page, "tsLast && tsLast.res && !document.querySelector("
               "'#tsnote .bad') && !document.querySelector('#tsnote .spin')",
         900, retry="showTransit(transitCandidates()[0].key, true);")
    time.sleep(2)
    shots.region('#tsbuttons', '#tsplot', 'tess_fold')
    page.click('[data-tsview="series"]')
    time.sleep(3)
    shots.region('#tsbuttons', '#tsplot', 'tess_series')
    page.click('[data-tsview="fold"]')
    time.sleep(2)
    page.locator('[data-help="transit_search"]').hover()
    time.sleep(0.8)
    box = page.locator('#tip').bounding_box()
    scroll = page.evaluate('window.scrollY')
    page.screenshot(path=shots.path('help'), full_page=True, clip=dict(
        x=max(0, box['x'] - 30), y=box['y'] + scroll - 40,
        width=min(WIDTH - box['x'] + 30, box['width'] + 60),
        height=box['height'] + 60))
    page.mouse.move(5, 5)


#: the star whose datasets show the rules of koloa.datasets (many
#: releases of the same spectra, old velocities)
RULES_STAR = 'GJ 876'


def datasets(page, shots: Shots, numbers: dict) -> None:
    """which datasets are used: a star from its archives alone, all of
    them ticked; the table says what the rules leave out, and why"""
    page.fill('#root', 'archives')
    page.fill('#target', RULES_STAR)
    page.press('#target', 'Enter')
    wait(page, "document.querySelector('#ident .stats-grid')",
         retry='resolveStar();')
    time.sleep(2)
    for key in ('dace', 'carmenes', 'vizier'):
        page.check(f'[data-mirror="{key}"]')
    page.click('#plot')
    wait(page, "document.querySelector('#rvplot.on') && lastRV",
         retry="document.getElementById('plot').click();")
    time.sleep(2)
    shots.element('#rvtable', 'datasets')
    numbers['datasets'] = dict(star=RULES_STAR, rows=page.evaluate(
        '() => lastRV.map((one) => ({ name: one.name, n: one.n, '
        'total: one.total, source: one.source, status: one.status, '
        'rule: one.rule }))'))


#: the survey of the screenshots: the M dwarfs this near [pc]
SURVEY_DMAX = 7


def survey_tab(page, shots: Shots, numbers: dict) -> None:
    """a survey: the M dwarfs within a few parsecs asked of SIMBAD, their
    archives checked (public data), the files of the batch put with them,
    a route to a server (not a real one), the batch packed, and the copy
    of its tar typed in the terminal (not run)"""
    page.click('[data-tab="survey"]')
    page.fill('#sv-dmax', str(SURVEY_DMAX))
    page.click('#sv-ask')
    wait(page, 'survey && survey.stars', 300,
         retry="document.getElementById('sv-ask').click();")
    page.click('#sv-check')
    wait(page, "survey.check && survey.check.status !== 'running'", 900)
    time.sleep(3)
    page.fill('#sv-folders', 'files')
    page.fill('#sv-pattern', '*.csv')
    page.click('#sv-match')
    settle(page, '#sv-matchstatus')
    page.click('#sv-tickfiles')
    time.sleep(1)
    shots.region('#tab-survey section', '#sv-card', 'survey',
                 inside='#sv-card')
    numbers['survey'] = page.evaluate(
        '() => ({ dmax: +document.getElementById("sv-dmax").value, '
        'n: survey.stars.length, '
        'data: survey.stars.filter(svHas).length, '
        'dace: survey.stars.filter((s) => s.archives && s.archives.dace)'
        '.length, carmenes: survey.stars.filter((s) => s.archives && '
        's.archives.carmenes).length, surveys: survey.stars.filter((s) => '
        's.archives && (s.archives.surveys || []).length).length, '
        'files: survey.stars.filter((s) => (s.files || []).length).length, '
        'rotation: survey.stars.filter((s) => (s.rotation || []).some('
        '(one) => !one.source.startsWith("CARMENES"))).length, '
        'unmatched: (survey.unmatched || []).length })')
    # a terminal, and how to get to a server (a made-up one)
    page.click('#term-new')
    wait(page, 'terms.length === 1', 60)
    time.sleep(2)
    page.click('#term-remember')
    time.sleep(1)
    page.fill('#route-name', 'server')
    page.fill('#route-host', 'me@server')
    page.fill('#route-folder', '/scratch/me/koloa_batches')
    page.fill('#route-lines', 'ssh me@server\ncd /scratch/me/koloa_batches'
                              '\nmodule load python scipy-stack')
    time.sleep(0.5)
    shots.element('#route-dialog', 'survey_route')
    page.click('#route-save')
    time.sleep(1.5)
    # the batch of the stars ticked, packed (their archives as checked)
    page.fill('#sv-name', 'm_dwarfs_demo')
    page.dispatch_event('#sv-name', 'input')
    page.uncheck('#sv-gather')
    page.click('#sv-pack')
    wait(page, "survey.pack && survey.pack.status !== 'running'", 600)
    time.sleep(1)
    shots.element('#sv-batchcard', 'survey_batch')
    page.click('#term-copy')
    time.sleep(3)
    shots.element('#term-card', 'survey_terminal')


def runs(page, shots: Shots) -> None:
    """a run of the command line: the archives of GJ 436 gathered again
    (public), its steps and its log"""
    page.fill('#root', 'archives')
    page.fill('#target', 'GJ 436')
    page.press('#target', 'Enter')
    wait(page, "document.querySelector('#ident .stats-grid')")
    time.sleep(2)
    page.click('#run-gather')
    wait(page, "[...jobs.values()].length && [...jobs.values()].every("
               "(job) => job.status !== 'running')", 900)
    time.sleep(2)
    shots.element('section:has(#jobs)', 'runs')


def from_file(page, shots: Shots, numbers: dict, path: str) -> None:
    """a file of velocities: its star found by its APERO name"""
    page.fill('input.fpath[data-row="0"]', path)
    page.dispatch_event('input.fpath[data-row="0"]', 'change')
    wait(page, "document.querySelector('#ident .stats-grid')", 300)
    time.sleep(2)
    shots.element('#starcard', 'star_from_file')
    numbers['from_file'] = dict(source=page.inner_text('#targetsrc'),
                                target=page.input_value('#target'))


def batch_run(page, shots: Shots, numbers: dict, files,
              archives: bool = False) -> str:
    """the quick FIP of every file: each on its own, or with every archive
    of its star (many datasets a star: hours for a few stars)"""
    page.click('[data-tab="batch"]')
    page.evaluate(f'() => addBatchFiles({json.dumps(files)})')
    if archives:
        page.check('#batcharchives')
    page.click('#batchrun')
    wait(page, "batch && batch.status !== 'running'", 7200,
         retry="if (!batch) document.getElementById('batchrun').click(); "
               "else if (batch.status === 'running') { batchMisses = 0; "
               "pollBatch(batch.id); }")
    time.sleep(2)
    shots.element('#tab-batch', 'batch')
    numbers['batch'] = page.evaluate(
        '() => batch.items.map((it) => ({ name: it.name, star: it.star, '
        'summary: it.summary }))')
    # how long it took, to show it again without running it
    numbers['batch_elapsed'] = page.evaluate('() => batch.elapsed')
    return page.evaluate('() => batch.id')


def batch_again(page, shots: Shots, bid, files, kept=None) -> None:
    """a batch already run, shown again (in another language): the one
    the server still has (its id), else the one kept in numbers.json (its
    lines and how long it took: no server has it any more)"""
    page.click('[data-tab="batch"]')
    page.evaluate(f'() => addBatchFiles({json.dumps(files)})')
    if bid:
        page.evaluate("async () => { batch = await api('/api/batch?id="
                      f"{bid}'); renderBatch(); }}")
    else:
        state = dict(id='kept', status='done', archives=False,
                     elapsed=kept.get('batch_elapsed'),
                     items=[dict(path=os.path.join('files', item['name']),
                                 name=item['name'], status='done', qid=None,
                                 error=None, summary=item['summary'],
                                 star=item['star'], stage=None, note=None)
                            for item in kept['batch']])
        page.evaluate(f'() => {{ batch = {json.dumps(state)}; '
                      'renderBatch(); }')
    time.sleep(1.5)
    shots.element('#tab-batch', 'batch')


# =============================================================================
# Start of code
# =============================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('--work', default=os.path.join(HERE, 'gui_demo'))
    parser.add_argument('--port', type=int, default=8790)
    parser.add_argument('--skip-prep', action='store_true')
    parser.add_argument('--langs', default='en,fr',
                        help='the languages, of en and fr')
    parser.add_argument('--parts', default='name,file,runs,batch',
                        help='the walk-throughs, of name, file, runs and '
                             'batch; recall takes again, from the result '
                             'remembered, what needs no computation, '
                             'batchshow the batch kept in numbers.json, '
                             'datasets the table of the rules of a star, '
                             'survey the survey tab and its terminal')
    parser.add_argument('--server', default='',
                        help='a koloa GUI already running (its address: '
                             'http://127.0.0.1:8790), in the same folder '
                             'of data, in place of one started here: the '
                             'FIPs computed on another machine')
    parser.add_argument('--batch-archives', action='store_true',
                        help='the batch with every archive of each star')
    parser.add_argument('--koloa', default=ROOT,
                        help='the folder of the koloa package served (a '
                             'copy that does not change while this runs)')
    args = parser.parse_args()
    langs, parts = args.langs.split(','), args.parts.split(',')
    work = os.path.abspath(args.work)
    os.makedirs(work, exist_ok=True)
    if not args.skip_prep:
        prepare(work)
    files = sorted(os.path.join('files', name) for name in
                   os.listdir(os.path.join(work, 'files')))
    from playwright.sync_api import sync_playwright
    # the remembered targets of this run: a folder of its own
    mem = os.path.join(work, 'remembered_' + '_'.join(langs))
    proc = None if args.server else start_server(
        work, args.port, mem, os.path.abspath(args.koloa))
    address = args.server or f'http://127.0.0.1:{args.port}'
    # the key of the server started here, for the terminal of the page
    key = getattr(proc, 'koloa_key', '') if proc is not None else ''
    numbers: dict = dict(files=[os.path.basename(f) for f in files])
    # the batch of this server (run in the first language, shown again in
    #   the next)
    bid = None
    try:
        with sync_playwright() as play:
            browser = play.chromium.launch(channel='chrome')
            for part in parts:
                for lang in langs:
                    print(part, lang)
                    shots = Shots(None, lang)
                    ctx = browser.new_context(
                        viewport=dict(width=WIDTH, height=1000),
                        locale='en-GB' if lang == 'en' else 'fr-CA')
                    ctx.add_init_script(
                        "try { localStorage.clear(); localStorage.setItem("
                        f"'koloa-lang', '{lang}'); }} catch (e) {{}}")
                    page = ctx.new_page()
                    page.goto(address.rstrip('/') + '/'
                              + (f'?key={key}' if key else ''))
                    shots.page = page
                    keep = numbers if lang == 'en' else {}
                    if part == 'name':
                        # the remembered targets of this language only: an
                        #   earlier list set aside
                        if os.path.isdir(mem) and not args.server:
                            os.rename(mem, f'{mem}_{time.time():.0f}')
                        by_name(page, shots, keep)
                    elif part == 'file':
                        # the file of GJ 581
                        from_file(page, shots, keep, next(
                            (one for one in files if '581' in one), files[0]))
                    elif part == 'runs':
                        runs(page, shots)
                    elif part == 'recall':
                        recalled(page, shots)
                    elif part == 'datasets':
                        datasets(page, shots, keep)
                    elif part == 'plot':
                        plotted(page, shots)
                    elif part == 'survey':
                        survey_tab(page, shots, keep)
                    elif part == 'batchshow':
                        # the batch kept in numbers.json, not run again
                        with open(os.path.join(OUT, 'numbers.json')) as handle:
                            batch_again(page, shots, None, files,
                                        json.load(handle))
                    elif bid is None:
                        bid = batch_run(page, shots, keep, files,
                                        args.batch_archives)
                    else:
                        batch_again(page, shots, bid, files)
                    ctx.close()
            browser.close()
    finally:
        if proc is not None:
            proc.terminate()
    # the numbers of the parts not run now, as they were
    path = os.path.join(OUT, 'numbers.json')
    kept = json.load(open(path)) if os.path.exists(path) else {}
    kept.update(numbers)
    with open(path, 'w') as handle:
        json.dump(kept, handle, indent=1, default=float)
    print('numbers.json written')


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

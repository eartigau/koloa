#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fill the numbers of the web page from the computations.

The blocks of index.html between <!--AUTO:name--> and <!--/AUTO--> are
rewritten from the campaign summary (paper/campaign_summary.json) and the
outputs of the demos, so the page never shows a number typed by hand.
The examples card is written from docs/examples/outputs.json (run
docs/examples/run_examples.py first).

    python make_page.py

Created on 2026-09-27

@author: artigau
"""
import json
import os
import re
from html import escape

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
#: the figures of the examples are shown at most this many times their size
EXAMPLE_SCALE = 1.3
#: the timestamp that starts every line of koloa's log
TIMESTAMP = re.compile(r'(\d{6} \d{2}:\d{2}:\d{2}\.\d{2} \| )(.*)')
METHODS = [('gauss', 'one signal, fixed jitter'),
           ('soft', 'one signal, soft clip'),
           ('hard', 'one signal, hard clip'),
           ('white', 'sampled, white jitter'),
           ('binned', 'sampled, visit means'),
           ('visit', 'sampled, visit jitter'),
           ('koloa', 'koloa')]


def fill(html, name, content):
    """Replace one AUTO block"""
    pattern = re.compile(r'(<!--AUTO:' + name + r'-->)(.*?)(<!--/AUTO-->)',
                         re.S)
    return pattern.sub(lambda mt: mt.group(1) + content + mt.group(3), html)


def campaign(html):
    """The rates of the campaign"""
    path = os.path.join(ROOT, 'paper', 'campaign_summary.json')
    if not os.path.exists(path):
        return html
    summ = json.load(open(path))
    stats = summ['stats']
    fpv = stats['visits_null']
    html = fill(html, 'fp_gauss', f'{100 * fpv["gauss"]["fp"][0]:.0f}%')
    html = fill(html, 'fp_koloa', f'{100 * fpv["koloa"]["fp"][0]:.0f}%')
    flags = summ['flags'].get('visits_null')
    if flags:
        html = fill(html, 'recall', f'{100 * flags["recall"]:.0f}%')
        html = fill(html, 'falseflag',
                    f'false flags {100 * flags["false_rate"]:.1f}%')
    cal = summ.get('calibration', {})
    for key in ('koloa', 'gauss'):
        if key in cal:
            high = cal[key][2]
            if high['n']:
                html = fill(html, f'cal_{key}',
                            f'{100 * high["true"] / high["n"]:.0f}% '
                            f'({high["true"]} of {high["n"]})')
    bycase = summ.get('calibration_by_case', {})
    for key, case in (('gauss', 'clean'), ('gauss', 'visits'),
                      ('white', 'clean')):
        if case in bycase:
            high = bycase[case][key][2]
            html = fill(html, f'cal_{key}_{case}',
                        f'{high["true"]} of {high["n"]}')
    for key in ('gauss', 'koloa'):
        det = summ['stats'].get('clean_planet', {}).get(key)
        if det:
            html = fill(html, f'det_{key}_clean', f'{100 * det["det"][0]:.0f}%')
    cases = [(case, label) for case, label in
             (('clean', 'clean'), ('visits', 'bad visits'),
              ('spikes', 'spikes'), ('border', 'borderline'))
             if f'{case}_null' in stats]
    ncase = len(cases)
    heads = ''.join(f'<th>{label}</th>' for _, label in cases)
    rows = ['<table class="koloa">',
            f'<tr><th>method</th><th colspan="{ncase}">false positives '
            f'(no planet)</th><th colspan="{ncase}">detections (K = 3 m/s)'
            f'</th></tr>', f'<tr><th></th>{heads}{heads}</tr>']
    for key, label in METHODS:
        cells = []
        for case, _ in cases:
            val = stats.get(f'{case}_null', {}).get(key)
            cells.append(f'{100 * val["fp"][0]:.0f}%' if val else '')
        for case, _ in cases:
            val = stats.get(f'{case}_planet', {}).get(key)
            cells.append(f'{100 * val["det"][0]:.0f}%' if val else '')
        rows.append(f'<tr><td>{label}</td>' + ''.join(
            f'<td class="num">{cell}</td>' for cell in cells) + '</tr>')
    rows.append('</table>')
    html = fill(html, 'campaign_table', '\n'.join(rows))
    matched = summ.get('matched')
    if matched:
        mrows = ['<table class="koloa">',
                 f'<tr><th>detections at a matched false-positive rate '
                 f'(none of 100 empty series)</th>{heads}</tr>']
        for key, label in METHODS:
            cells = [f'{matched[case][key]["det_0.00"][0]:.2f}'
                     if key in matched.get(case, {}) else ''
                     for case, _ in cases]
            mrows.append(f'<tr><td>{label}</td>' + ''.join(
                f'<td class="num">{cell}</td>' for cell in cells) + '</tr>')
        mrows.append('</table>')
        html = fill(html, 'matched_table', '\n'.join(mrows))
    return html


def toi2120(html):
    """What koloa finds in TOI-2120, as observed and with added outliers"""
    import numpy as np
    base = os.path.join(ROOT, 'demos', 'output', 'toi2120')
    path = os.path.join(base, 'clean_series.json')
    if not os.path.exists(path):
        return html
    dat = json.load(open(path))
    vis = ', '.join(f'{2400000 + vv["time"]:.1f} ({vv["offset"]:+.0f}&nbsp;m/s)'
                    for vv in dat['visits'])
    text = (f'<p>In the series <em>as observed</em>, koloa flags '
            f'{len(dat["visits"])} visits, at BJD {vis}. With them set aside, '
            f'the GP amplitude falls from {dat["gaussian"]["gp_sigma"]:.1f} to '
            f'{dat["koloa"]["gp_sigma"]:.1f}&nbsp;m/s and the visit jitter from '
            f'{dat["gaussian"]["seq_jitter"]:.1f} to '
            f'{dat["koloa"]["seq_jitter"]:.1f}&nbsp;m/s, and '
            f'K = {dat["koloa"]["K"][0]:.1f} &plusmn; '
            f'{0.5 * (dat["koloa"]["K"][1] + dat["koloa"]["K"][2]):.1f}&nbsp;m/s '
            f'(Gaussian: {dat["gaussian"]["K"][0]:.1f} &plusmn; '
            f'{0.5 * (dat["gaussian"]["K"][1] + dat["gaussian"]["K"][2]):.1f}).')
    res = [rr for rr in np.load(os.path.join(base, 'results.npy'),
                                allow_pickle=True) if 'scan' in rr]
    sub = [rr for rr in res if rr['scan'] == 'fraction'
           and np.isclose(rr['level'], 0.2)]
    if sub:
        gsc = np.std([rr['gaussian'][0] for rr in sub])
        ksc = np.std([rr['koloa'][0] for rr in sub])
        text += (f' With 20% of the visits moved by 60 to 90&nbsp;m/s, the '
                 f'Gaussian K scatters by {gsc:.1f}&nbsp;m/s from one '
                 f'realisation to the next, and koloa\'s by {ksc:.1f}. At 30%, '
                 f'the bad visits become a population comparable to the good '
                 f'one, and koloa sometimes falls back on the Gaussian '
                 f'answer.</p>')
    return fill(html, 'toi_summary', text)


def _svg_size(path):
    """
    The size of a matplotlib SVG in CSS pixels, from its size in points

    :param path: str, the SVG file

    :return: tuple, width and height [px], or (None, None)
    """
    if not os.path.exists(path):
        return None, None
    with open(path) as handle:
        head = handle.read(4000)
    tag = re.search(r'<svg\b[^>]*>', head)
    if tag is None:
        return None, None
    width = re.search(r'\bwidth="([\d.]+)pt"', tag.group(0))
    height = re.search(r'\bheight="([\d.]+)pt"', tag.group(0))
    if width is None or height is None:
        return None, None
    return (round(float(width.group(1)) * 4 / 3),
            round(float(height.group(1)) * 4 / 3))


def _output_html(text):
    """
    What an example printed, one block per line so that a long line wraps
    under its message, with the timestamps of the log set apart

    :param text: str, the output

    :return: str, the HTML (escaped)
    """
    rows = []
    for line in text.split('\n'):
        match = TIMESTAMP.match(line)
        if match:
            rows.append(f'<span class="ex-line ts"><span class="ex-ts">'
                        f'{escape(match.group(1))}</span>'
                        f'{escape(match.group(2))}</span>')
        else:
            rows.append(f'<span class="ex-line">{escape(line) or " "}</span>')
    return ''.join(rows)


def _example_view(entry):
    """
    The page of one module: its code, what it printed and its figure

    :param entry: dict, one entry of docs/examples/outputs.json

    :return: str, the page (HTML)
    """
    mod = escape(entry['module'])
    status = '' if entry['ok'] else 'this example failed'
    parts = [f'<section class="view ex-panel" id="ex-panel-{mod}" '
             f'data-group="examples" data-title="{mod}">',
             '<p class="crumb"><a href="#examples">Every module, by example</a></p>',
             f'<h2><code>koloa.{mod}</code> {escape(entry["title"])}</h2>']
    if entry.get('blurb'):
        parts.append(f'<p>{escape(entry["blurb"])}</p>')
    script = escape(entry['script'], quote=True)
    parts += [f'<div class="ex-label"><span>code</span>'
              f'<a href="{script}">docs/{escape(entry["script"])}</a></div>',
              f'<pre class="codeblock ex-code"><code>'
              f'{escape(entry["code"].rstrip())}</code></pre>',
              f'<div class="ex-label"><span>output</span>'
              f'<span class="ex-meta">{status}</span></div>',
              f'<pre class="ex-output{"" if entry["ok"] else " ex-failed"}">'
              f'{_output_html(entry["output"])}</pre>']
    if entry.get('figure'):
        width, height = _svg_size(os.path.join(HERE, entry['figure']))
        size = ''
        if width:
            size = (f' width="{width}" height="{height}" style="max-width: '
                    f'{round(EXAMPLE_SCALE * width)}px"')
        caption = entry.get('caption') or ''
        alt = escape(caption or entry['title'], quote=True)
        parts.append(f'<figure class="lecture-figure ex-figure">'
                     f'<img src="{escape(entry["figure"], quote=True)}" '
                     f'alt="{alt}" loading="lazy"{size}>'
                     f'<figcaption>{escape(caption)}</figcaption></figure>')
    parts.append('</section>')
    return '\n'.join(parts)


def examples(html):
    """
    The pages of the modules, one each, and the list of them on the page
    that introduces them, from docs/examples/outputs.json (written by
    docs/examples/run_examples.py)

    :param html: str, the page

    :return: str, the page with the examples_views and examples_index
             blocks filled
    """
    path = os.path.join(HERE, 'examples', 'outputs.json')
    if not os.path.exists(path):
        return html
    with open(path) as handle:
        entries = json.load(handle)
    if not entries:
        return html
    views = [_example_view(entry) for entry in entries]
    index = ['<ul class="page-index">'] + [
        f'<li><a href="#ex-panel-{escape(entry["module"])}">'
        f'<code>koloa.{escape(entry["module"])}</code></a>'
        f'<span>{escape(entry["title"])}</span></li>' for entry in entries]
    index.append('</ul>')
    html = fill(html, 'examples_views', '\n'.join(views))
    return fill(html, 'examples_index', '\n'.join(index))


def _comment_spans(code):
    """
    The (line, start column) of every comment of a snippet: Python's own
    tokenizer, so a # inside a string is not one; a snippet that does not
    tokenize (a shell command) has its comments at a # that starts a line
    or follows a space

    :param code: str, the snippet (not escaped)

    :return: dict, line index: start column
    """
    import io
    import tokenize
    out = {}
    try:
        for tok in tokenize.generate_tokens(io.StringIO(code).readline):
            if tok.type == tokenize.COMMENT:
                out[tok.start[0] - 1] = tok.start[1]
        return out
    except (tokenize.TokenError, IndentationError, SyntaxError):
        out = {}
    for it, line in enumerate(code.split('\n')):
        mt = re.search(r'(^|\s)#', line)
        if mt:
            out[it] = mt.start() + len(mt.group(1))
    return out


def colour_comments(html):
    """
    The comments of every code snippet of the page in grey (span.c), found
    by _comment_spans; any span.c already there is replaced
    """
    import html as htmlmod

    def recolour(mt):
        head, body, tail = mt.group(1), mt.group(2), mt.group(3)
        code = htmlmod.unescape(re.sub(r'</?span[^>]*>', '', body))
        spans = _comment_spans(code)
        lines = []
        for it, line in enumerate(code.split('\n')):
            col = spans.get(it)
            if col is None:
                lines.append(escape(line))
            else:
                lines.append(escape(line[:col]) + '<span class="c">'
                             + escape(line[col:]) + '</span>')
        return head + '\n'.join(lines) + tail
    html = re.sub(r'(<div class="codeblock">)(.*?)(</div>)', recolour, html,
                  flags=re.S)
    return re.sub(r'(<pre class="codeblock ex-code"><code>)(.*?)(</code></pre>)',
                  recolour, html, flags=re.S)


def main():
    path = os.path.join(HERE, 'index.html')
    html = open(path).read()
    html = campaign(html)
    html = toi2120(html)
    html = examples(html)
    # the cards of the other demos, and the mathematics
    from page_cards import fill_cards
    html = fill_cards(html)
    html = colour_comments(html)
    open(path, 'w').write(html)
    print('index.html filled')
    # the French page, from the catalogue of translations (docs/i18n/fr.json)
    from i18n import french
    page, missing = french(html)
    open(os.path.join(HERE, 'index_fr.html'), 'w').write(page)
    print(f'index_fr.html written; {len(missing)} blocks not translated yet '
          f'(docs/i18n/missing_fr.json)')


if __name__ == '__main__':
    main()

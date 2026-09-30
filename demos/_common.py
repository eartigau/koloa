#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
What every demo shares: where the outputs go, and the two figure styles.

Each demo draws its figures twice: in the 'paper' style as PDF, into
demos/output/<demo>/, and in the 'web' style as SVG, into docs/figures/,
where the web page picks them up (both styles are light).

Created on 2026-09-27

@author: artigau
"""
import os

from koloa import data as kdata
from koloa import plotting as kplot
from koloa.data import RVData

# the demos are the numbers of the site and the paper, made on the
#   exposures (each exposure or visit an outlier); koloa's default became
#   the nightly means (2026-09-30), so they keep the exposures, here and in
#   their worker processes
os.environ['KOLOA_NIGHTLY'] = '0'
kdata.NIGHTLY = False

# =============================================================================
# Define variables
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
#: a real NIRPS sampling (times and error bars, no velocities), the
#: template of every simulation
DATA = os.path.join(ROOT, 'data', 'nirps_template.csv')
WEBDIR = os.path.join(ROOT, 'docs', 'figures')


# =============================================================================
# Define functions
# =============================================================================
def outdir(name: str) -> str:
    """The output folder of a demo (created)"""
    path = os.path.join(HERE, 'output', name)
    os.makedirs(path, exist_ok=True)
    os.makedirs(WEBDIR, exist_ok=True)
    return path


def template() -> RVData:
    """A real NIRPS sampling, whose times, visits and error bars the
    simulations borrow"""
    return RVData.from_csv(DATA, name='NIRPS template')


def save_both(make, name: str, demo: str):
    """
    Draw a figure in both styles and save it: PDF for the paper, SVG for
    the web

    :param make: callable, returns a matplotlib figure (called once per
                 style, after the style is set)
    :param name: str, the file name without extension
    :param demo: str, the demo
    """
    kplot.set_style('paper')
    kplot.savefig(make(), os.path.join(outdir(demo), name + '.pdf'))
    kplot.set_style('web')
    kplot.savefig(make(), os.path.join(WEBDIR, f'{demo}_{name}.svg'))
    kplot.set_style('paper')

# =============================================================================
# End of code
# =============================================================================


def web_fips(fips: dict) -> dict:
    """
    The FIPs a figure shows: all of them for the paper, and on the site only
    the fixed-jitter FIP and koloa's (the clips were found confusing there;
    the tables of the site still compare every method)

    :param fips: dict, method key: FIPResult

    :return: dict, the FIPs to draw in the current style
    """
    if kplot._STYLE[0] != 'web':
        return fips
    return {key: val for key, val in fips.items()
            if key not in ('soft', 'hard')}


def flag_summary(prob, labels, threshold: float = 0.5) -> dict:
    """
    What koloa flagged against what was injected, the clear and the
    borderline outliers counted apart

    :param prob: np.ndarray, the outlier probability of every exposure
    :param labels: np.ndarray, the recipe of every exposure ('' when good;
                   simulate's outlier_label)
    :param threshold: float, flagged above this probability

    :return: dict, per kind (clear, borderline, good): n, flagged and the
             mean outlier probability
    """
    import numpy as np
    prob = np.asarray(prob, dtype=float)
    labels = np.asarray(labels)
    flagged = prob > threshold
    out = {}
    for kind, sel in (('clear', labels == 'clear'),
                      ('borderline', labels == 'borderline'),
                      ('good', labels == '')):
        out[kind] = dict(n=int(np.sum(sel)),
                         flagged=int(np.sum(flagged & sel)),
                         mean_prob=float(np.mean(prob[sel]))
                         if np.any(sel) else float('nan'))
    return out


def log_flags(summ: dict, who: str = 'koloa'):
    """One log line from flag_summary"""
    from koloa.log import log
    clear, border, good = summ['clear'], summ['borderline'], summ['good']
    log(f'  {who} flags: clear outliers {clear["flagged"]}/{clear["n"]}, '
        f'borderline {border["flagged"]}/{border["n"]} (mean P(outlier) '
        f'{border["mean_prob"]:.2f}), good exposures {good["flagged"]}/'
        f'{good["n"]}', 'value')


def save_summary(demo: str, summary: dict, name: str = 'summary.json'):
    """The numbers of a demo, for the web page and the paper"""
    import json
    with open(os.path.join(outdir(demo), name), 'w') as handle:
        json.dump(summary, handle, indent=1, default=float)

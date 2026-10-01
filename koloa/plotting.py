#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figures where every point wears how much it is trusted.

The colour of a point is its probability of being an outlier, on a
diverging scale: blue is reliable, grey is undecided (0.5), red is an
outlier. Colour is never alone: a point more likely an outlier than not is
also drawn hollow.

Every method has one colour everywhere, in every figure: koloa's
outlier-aware results are blue, the plain gaussian orange, the soft clip
aqua, the hard clip yellow, and the sampled gaussian models magenta (white
jitter) and green (visit jitter). Two styles exist, 'paper' (for PDF) and
'web' (for the site, SVG); both are light, with the same steps, validated
for colour-blind separation on a white background. (The web style was dark
until 2026-09-29, when the site turned light; STYLES_DARK keeps it.)

The figures:
- timeseries: velocities in time (per point, or per sequence);
- phase: a phase fold with binned weighted means and the fitted curve
  carrying its one sigma envelope;
- sequences: every visit, its exposures and its outlier probability;
- periodograms: GLS and OAP on the same axis, the window, and the FIP
  periodograms side by side;
- jackknife: the periodogram envelope with each unit left out, and which
  units hold a peak up;
- coherence: the amplitude and phase of a signal chunk by chunk;
- corner: several posteriors on the same corner plot, each in its colour;
- recovery_map, detection_limits: what a series could have found;
- outlier_keys: why the outliers are outliers (koloa.outliers).

Created on 2026-09-27

@author: artigau
"""
import os
from typing import Any, Dict, Optional, Sequence, Tuple

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from koloa.data import RVData  # noqa: E402

# =============================================================================
# Define variables
# =============================================================================
STYLES = {
    'paper': dict(
        surface='#ffffff', text='#0b0b0b', muted='#52514e', grid='#e6e5e1',
        koloa='#2a78d6', gaussian='#eb6834', soft='#1baf7a', hard='#eda100',
        white='#e87ba4', visit='#008300', reliable='#2a78d6', neutral='#9a9893', outlier='#e34948',
        model='#0b0b0b', envelope='#2a78d6', binned='#0b0b0b',
        errbar='#b9b8b3', threshold='#52514e'),
}
#: the site is light: its figures take the steps of the paper
STYLES['web'] = dict(STYLES['paper'])
#: the dark style the site used until 2026-09-29 (to recolour its figures)
STYLES_DARK = dict(
    surface='#0e182a', text='#e8eef8', muted='#a8b4ca', grid='#26334a',
    koloa='#3987e5', gaussian='#d95926', soft='#199e70', hard='#c98500',
    white='#d55181', visit='#008300', reliable='#3987e5', neutral='#8a93a3',
    outlier='#e66767', model='#e8eef8', envelope='#3987e5', binned='#e8eef8',
    errbar='#3a4760', threshold='#a8b4ca')
#: the steps of the sequential colour map, light (both styles) and dark
SEQUENTIAL_STOPS = ['#f3f7fd', '#a9c8ef', '#2a78d6', '#0f3d78']
SEQUENTIAL_STOPS_DARK = ['#13233d', '#1d4f8f', '#3987e5', '#a9cdf7']

#: the colours of the current style
C = dict(STYLES['paper'])
_STYLE = ['paper']

METHOD_LABELS = {'koloa': 'koloa (outlier-aware)', 'gaussian': 'Gaussian',
                 'soft': 'soft clip', 'hard': 'hard clip'}


# =============================================================================
# Style
# =============================================================================
def set_style(name: str = 'paper'):
    """
    Switch between the paper style (PDF figures) and the web style (SVG
    figures for the site); both are light

    :param name: str, paper or web
    """
    C.clear()
    C.update(STYLES[name])
    _STYLE[0] = name
    plt.rcParams.update({
        'figure.facecolor': C['surface'], 'axes.facecolor': C['surface'],
        'savefig.facecolor': C['surface'], 'axes.edgecolor': C['muted'],
        'axes.labelcolor': C['text'], 'text.color': C['text'],
        'xtick.color': C['muted'], 'ytick.color': C['muted'],
        'xtick.labelcolor': C['text'], 'ytick.labelcolor': C['text'],
        'axes.grid': True, 'grid.color': C['grid'], 'grid.linewidth': 0.6,
        'axes.axisbelow': True, 'axes.spines.top': False,
        'axes.spines.right': False, 'font.size': 9,
        'axes.titlesize': 9.5, 'axes.labelsize': 9, 'legend.fontsize': 8,
        'legend.frameon': False, 'lines.linewidth': 1.4,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'figure.dpi': 110, 'savefig.bbox': 'tight',
        'pdf.fonttype': 42, 'svg.fonttype': 'path'})


def reliability_cmap() -> LinearSegmentedColormap:
    """Blue (reliable) to grey (undecided) to red (outlier)"""
    return LinearSegmentedColormap.from_list(
        'koloa_outlier', [C['reliable'], C['neutral'], C['outlier']])


def savefig(fig: Any, path: str, close: bool = True) -> str:
    """
    Save a figure (PDF for the paper style, SVG for the web style unless
    the path says otherwise)

    :param fig: matplotlib figure
    :param path: str, the file (its extension is kept if given)
    :param close: bool, close the figure after

    :return: str, the file written
    """
    path = os.fspath(path)
    if '.' not in os.path.basename(path):
        path += '.svg' if _STYLE[0] == 'web' else '.pdf'
    fig.savefig(path)
    if close:
        plt.close(fig)
    return path


set_style('paper')


# =============================================================================
# Helpers
# =============================================================================
def _scatter_reliability(ax, time, value, err, prob, size=16, zorder=3,
                         errorbars=True):
    """Points coloured by outlier probability, hollow when more likely out"""
    prob = np.zeros(len(time)) if prob is None else np.asarray(prob)
    cmap = reliability_cmap()
    if errorbars:
        ax.errorbar(time, value, err, fmt='none', ecolor=C['errbar'],
                    elinewidth=0.8, zorder=zorder - 1)
    good = prob <= 0.5
    ax.scatter(time[good], value[good], c=prob[good], cmap=cmap, vmin=0,
               vmax=1, s=size, edgecolor=C['surface'], linewidth=0.5,
               zorder=zorder)
    if np.any(~good):
        colors = cmap(prob[~good])
        ax.scatter(time[~good], value[~good], facecolor='none',
                   edgecolor=colors, s=size * 1.8, linewidth=1.4,
                   zorder=zorder + 1)
    return cmap


def _colorbar(fig, ax, cmap, label='P(outlier)'):
    """The outlier probability scale"""
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cbar = fig.colorbar(sm, ax=ax, pad=0.01, fraction=0.035, aspect=30)
    cbar.set_label(label)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(colors=C['muted'], labelcolor=C['text'])
    return cbar


def sequence_means(data: RVData, value: np.ndarray, prob: Optional[np.ndarray]
                   = None, jitter: float = 0.0) -> Dict[str, np.ndarray]:
    """
    Weighted means of each sequence, with the sequence's outlier probability

    :param data: RVData, the series
    :param value: np.ndarray, the values to average (n)
    :param prob: np.ndarray or None, the outlier probability per point
    :param jitter: float, a jitter added to the errors [m/s]

    :return: dict, time, value, err, prob (one per sequence)
    """
    weight = 1.0 / (data.err ** 2 + jitter ** 2)
    nseq = data.nseq
    wsum = np.bincount(data.seq, weights=weight, minlength=nseq)
    out = dict(time=np.bincount(data.seq, weights=weight * data.time,
                                minlength=nseq) / wsum,
               value=np.bincount(data.seq, weights=weight * value,
                                 minlength=nseq) / wsum,
               err=1 / np.sqrt(wsum))
    if prob is None:
        out['prob'] = np.zeros(nseq)
    else:
        out['prob'] = np.bincount(data.seq, weights=prob, minlength=nseq) / \
            np.bincount(data.seq, minlength=nseq)
    return out


# =============================================================================
# Time series
# =============================================================================
#: the instruments, in order: eight hues and a marker each
INST_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4',
                '#008300', '#4a3aa7', '#e34948']
INST_MARKERS = ['o', 's', 'D', '^', 'v', 'P', '*', 'h']


def trend_curve(model: Any, theta: np.ndarray, time: np.ndarray
                ) -> np.ndarray:
    """the polynomial trend of a fit (and a perspective acceleration fitted
    beside it) at any time"""
    tnorm = (np.asarray(time) - model.tref) / max(model.data.baseline, 1e-9)
    out = np.zeros(len(tnorm))
    for deg in range(1, model.trend + 1):
        out += theta[model.index[f'trend_{deg}']] * tnorm ** deg
    if 'secacc' in model.index:
        out += theta[model.index['secacc']] * (np.asarray(time)
                                                - model.tref) / 365.25
    return out


def model_series(res: Any, accel: Optional[Dict[str, Any]] = None,
                 title: Optional[str] = None):
    """
    The velocities with the best model of a fit, in three panels sharing
    time: every instrument with its offset (and decorrelation) taken out,
    with the Keplerians and the trend; the planets taken out, with the
    trend (the acceleration of the star and its change, with their errors
    when given); the residuals. An exposure more likely an outlier than not
    is hollow.

    :param res: FitResult, the fit (its maximum a posteriori is drawn)
    :param accel: dict or None, from koloa.secular.acceleration
    :param title: str or None, the title

    :return: matplotlib figure
    """
    model, theta = res.model, res.theta
    data = model.data
    npl = len(model.planets)
    # what is taken out of every panel: the offsets and decorrelation
    trend_at = trend_curve(model, theta, data.time)
    shift = model.systematics(theta) - trend_at
    planets_at = sum((model.planet_rv(theta, ip) for ip in range(npl)),
                     np.zeros(data.n))
    shown = data.rv - shift
    resid = data.rv - model.mean_model(theta)
    # the curves, fine enough for the shortest period
    pers = [model.orbit(theta, ip)[0] for ip in range(npl)]
    nfine = int(np.clip(20 * data.baseline / min(pers + [data.baseline]),
                        2000, 200000))
    fine = np.linspace(data.time.min(), data.time.max(), nfine)
    trend_fine = trend_curve(model, theta, fine)
    total_fine = trend_fine + sum((model.planet_rv(theta, ip, fine)
                                   for ip in range(npl)), np.zeros(nfine))
    prob = res.outlier_prob if res.outlier_prob is not None else \
        np.zeros(data.n)
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.6), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.4, 1.0, 0.8]))
    ax1, ax2, ax3 = axes
    for it, inst in enumerate(data.instruments):
        sel = data.inst == inst
        colour = INST_COLOURS[it % len(INST_COLOURS)]
        marker = INST_MARKERS[it % len(INST_MARKERS)]
        for ax, yval in ((ax1, shown), (ax2, shown - planets_at),
                         (ax3, resid)):
            for bad in (False, True):
                part = sel & ((prob > 0.5) == bad)
                if not np.any(part):
                    continue
                ax.errorbar(data.time[part], yval[part], data.err[part],
                            fmt=marker, ms=3.2, lw=0.6, elinewidth=0.6,
                            color=colour,
                            mfc='none' if bad else colour,
                            label=(f'{inst} ({int(np.sum(sel))})'
                                   if ax is ax1 and not bad else None),
                            zorder=2)
    ax1.plot(fine, total_fine, color=C['model'], lw=0.6, alpha=0.8,
             zorder=3, label='Keplerians + trend')
    ax1.plot(fine, trend_fine, color=C['model'], lw=1.0, ls='--', zorder=4,
             label='trend')
    words = []
    if accel:
        for key, name, unit in (('accel', 'dv/dt', 'm/s/yr'),
                                ('jerk', 'd$^2$v/dt$^2$', 'm/s/yr$^2$')):
            if key in accel:
                val, low, high = accel[key]
                words.append(f'{name} = {val:+.3g} $\\pm$ '
                             f'{0.5 * (low + high):.2g} {unit}')
    ax2.plot(fine, trend_fine, color=C['model'], lw=1.2, ls='--', zorder=4,
             label='trend: ' + ', '.join(words) if words else 'trend')
    ax3.axhline(0, color=C['model'], lw=0.8, ls=':')
    ax1.set_ylabel('RV [m s$^{-1}$]')
    ax2.set_ylabel('planets out')
    ax3.set_ylabel('residuals')
    ax3.set_xlabel('BJD - 2400000')
    ax1.legend(fontsize=6.5, ncol=4, loc='upper left', frameon=False)
    ax2.legend(fontsize=7, loc='upper left', frameon=False)
    if title:
        ax1.set_title(title, fontsize=9)
    fig.tight_layout(h_pad=0.4)
    return fig


def timeseries(data: RVData, prob: Optional[np.ndarray] = None,
               fit: Any = None, level: str = 'point', ax=None,
               title: Optional[str] = None, ylabel: str = 'RV [m s$^{-1}$]'):
    """
    The velocities in time, coloured by their outlier probability

    :param data: RVData, the series
    :param prob: np.ndarray or None, the outlier probability per point
    :param fit: FitResult or None, a model to draw over the points
    :param level: str, point or sequence (one weighted mean per visit)
    :param ax: matplotlib axis or None
    :param title: str or None, the title
    :param ylabel: str, the label of the y axis

    :return: the figure
    """
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.2, 2.8))
    else:
        fig = ax.figure
    value = data.rv.copy()
    if fit is not None:
        value = value - fit.model.systematics(fit.theta)
    if level == 'sequence':
        seqm = sequence_means(data, value, prob)
        cmap = _scatter_reliability(ax, seqm['time'], seqm['value'],
                                    seqm['err'], seqm['prob'], size=22)
    else:
        cmap = _scatter_reliability(ax, data.time, value, data.err, prob)
    if fit is not None:
        tgrid = np.linspace(data.time.min(), data.time.max(), 4000)
        curve = np.zeros(len(tgrid))
        for ip in range(len(fit.model.planets)):
            curve += fit.model.planet_rv(fit.theta, ip, tgrid)
        if fit.model.gp is not None and fit.model.gp.get('kernel') == 'multi':
            # a GP per group of instruments: the model of each group over
            #   its own time span
            ax.plot(tgrid, curve, color=C['model'], lw=0.8, zorder=2,
                    label='orbit')
            for group in fit.model.gp_groups:
                tsub = data.time[group['idx']]
                inside = (tgrid >= tsub.min()) & (tgrid <= tsub.max())
                gpm, _ = fit.gp_prediction(tgrid[inside],
                                           inst=group['instruments'][0])
                ax.plot(tgrid[inside], curve[inside] + gpm, color=C['model'],
                        lw=0.6, alpha=0.6, zorder=2)
        else:
            if fit.model.gp is not None:
                gpm, _ = fit.gp_prediction(tgrid)
                curve += gpm
            ax.plot(tgrid, curve, color=C['model'], lw=0.8, zorder=2,
                    label='model')
    ax.set_xlabel('time [BJD - 2400000]')
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc='left')
    _colorbar(fig, ax, cmap)
    return fig


# =============================================================================
# Phase folds
# =============================================================================
def _phase_bins(phase, value, weight, nbins):
    """Weighted means in phase bins"""
    edges = np.linspace(0, 1, nbins + 1)
    idx = np.clip(np.digitize(phase, edges) - 1, 0, nbins - 1)
    wsum = np.bincount(idx, weights=weight, minlength=nbins)
    mean = np.bincount(idx, weights=weight * value, minlength=nbins)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = mean / wsum
        err = 1 / np.sqrt(wsum)
    centre = 0.5 * (edges[1:] + edges[:-1])
    good = wsum > 0
    return centre[good], mean[good], err[good]


def phase(fit: Any, planet: int = 0, level: str = 'point', nbins: int = 10,
          ax=None, title: Optional[str] = None, nsample: int = 300,
          select: Optional[np.ndarray] = None,
          ylim: Optional[Tuple[float, float]] = None, legend: bool = True,
          colorbar: bool = True):
    """
    A phase fold of one orbit, with binned means and the fitted envelope

    Everything but this orbit (offsets, trend, the other orbits, the GP) is
    taken out first. The binned means weight each point by its inverse
    variance times its reliability, so an outlier does not drag a bin. The
    curve is the best orbit; the band is the 16-84 % range of the orbit
    over parameter draws (the chain, or the Laplace covariance of a fit).

    :param fit: FitResult, the fit
    :param planet: int, the orbit
    :param level: str, point or sequence
    :param nbins: int, the number of phase bins
    :param ax: matplotlib axis or None
    :param title: str or None
    :param nsample: int, parameter draws for the envelope
    :param select: np.ndarray or None, the points shown and binned (a
                   boolean mask, e.g. one instrument of a joint fit); the
                   orbit and its envelope are the fit's
    :param ylim: tuple or None, the range of velocities (set by the points
                 when None); points beyond it are drawn at its edge
    :param legend: bool, draw the legend
    :param colorbar: bool, draw the scale of the outlier probability that
                     colours the points

    :return: the figure
    """
    model, data = fit.model, fit.model.data
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.6, 3.2))
    else:
        fig = ax.figure
    period, tperi, ecc, omega, amp = model.orbit(fit.theta, planet)
    from koloa.kepler import tp_to_tc
    tconj = tp_to_tc(tperi, period, ecc, omega)
    value = fit.residuals(keep_planet=planet)
    prob = fit.outlier_prob
    # the white jitter of each point's instrument
    jit = np.zeros(data.n)
    for it, inst in enumerate(data.instruments):
        name = f'log_jit_{inst}'
        if name in model.index:
            jit[data.inst_index == it] = float(np.exp(
                fit.theta[model.index[name]]))
    keep = (np.ones(data.n, dtype=bool) if select is None
            else np.asarray(select, dtype=bool))
    if level == 'sequence':
        seqm = sequence_means(data, value, prob, jit)
        # the sequences that hold a point shown
        shown = np.bincount(data.seq, weights=keep.astype(float),
                            minlength=data.nseq) > 0
        tt, vv, ee, pp = (seqm[key][shown] for key in
                          ('time', 'value', 'err', 'prob'))
    else:
        tt, vv, ee, pp = data.time[keep], value[keep], \
            np.sqrt(data.err ** 2 + jit ** 2)[keep], prob[keep]
    ph = ((tt - tconj) / period) % 1.0
    if ylim is None:
        spread = max(np.nanpercentile(np.abs(vv), 98), 1.2 * amp)
        # room above the data for the legend, laid flat, so it hides nothing
        ylo, yhi = -1.3 * spread, 1.95 * spread
    else:
        ylo, yhi = ylim
    # the highest point drawn: below the band of the legend when there is
    #   one; the points beyond are drawn at the edges as arrows, so that no
    #   flagged unit disappears from the figure and none hides the legend
    ceiling = ylo + (0.80 if legend else 0.96) * (yhi - ylo)
    shown = (vv >= ylo) & (vv <= ceiling)
    cmap = _scatter_reliability(ax, ph[shown], vv[shown], ee[shown],
                                pp[shown],
                                size=22 if level == 'sequence' else 14)
    # the model and its envelope
    grid = np.linspace(0, 1, 300)
    tgrid = tconj + grid * period
    draws = fit.samples(nsample)
    curves = np.array([model.planet_rv(th, planet, tconj + grid *
                                       model.orbit(th, planet)[0])
                       for th in draws])
    if len(curves) > 5:
        low, mid, high = np.percentile(curves, [16, 50, 84], axis=0)
        ax.fill_between(grid, low, high, color=C['model'], alpha=0.18,
                        lw=0, zorder=4, label='1$\\sigma$ envelope')
        # the pointwise median of the posterior curves: the most probable
        #   orbit of an eccentric fit can have a sharp periastron that the
        #   posterior as a whole does not have
        ax.plot(grid, mid, color=C['model'], lw=1.2, zorder=4,
                label='posterior median')
    else:
        ax.plot(grid, model.planet_rv(fit.theta, planet, tgrid),
                color=C['model'], lw=1.2, zorder=4, label='best orbit')
    weight = (1 - pp) / ee ** 2
    bc, bm, be = _phase_bins(ph, vv, weight, nbins)
    ax.errorbar(bc, bm, be, fmt='s', color=C['binned'], ms=4.5, mfc=C['binned'],
                mec=C['surface'], mew=0.6, elinewidth=1.2, capsize=0,
                zorder=5, label='binned mean')
    ax.set_xlim(0, 1)
    ax.set_ylim(ylo, yhi)
    low, high = vv < ylo, vv > ceiling
    for sel, ypos, marker in ((low, ylo + 0.04 * (yhi - ylo), 'v'),
                              (high, ceiling, '^')):
        if np.any(sel):
            ax.scatter(ph[sel], np.full(np.sum(sel), ypos), marker=marker,
                       s=28, color=C['outlier'], zorder=6,
                       label='beyond the range' if marker == 'v' or
                       not np.any(low) else None)
    ax.set_xlabel('phase (0 = conjunction)')
    ax.set_ylabel('RV [m s$^{-1}$]')
    if title is None:
        title = f'P = {period:.4f} d, K = {amp:.2f} m s$^{{-1}}$'
    ax.set_title(title, loc='left')
    if legend:
        ax.legend(loc='upper center', ncol=4, handlelength=1.2, fontsize=7,
                  columnspacing=1.0)
    if colorbar:
        _colorbar(fig, ax, cmap)
    return fig


# =============================================================================
# Sequences
# =============================================================================
def sequences(data: RVData, prob: Optional[np.ndarray] = None,
              fit: Any = None, title: Optional[str] = None,
              highlight: Optional[float] = 0.5):
    """
    Every visit: its exposures, their mean, and its outlier probability

    The top panel is the residual of the model (or the velocity) of each
    exposure, visit by visit, with the weighted mean of the visit; the
    bottom panel is the probability that the visit is an outlier. A visit
    whose exposures agree with each other but not with the rest is what a
    sequence-level outlier looks like, and why koloa can make the visit
    the unit of an outlier.

    :param data: RVData, the series
    :param prob: np.ndarray or None, the outlier probability per point
    :param fit: FitResult or None, the model whose residual is shown
    :param title: str or None
    :param highlight: float or None, label visits above this probability

    :return: the figure
    """
    value = data.rv if fit is None else fit.residuals()
    prob = np.zeros(data.n) if prob is None else prob
    seqm = sequence_means(data, value, prob)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 3.8), sharex=True,
                                   gridspec_kw=dict(height_ratios=[2.2, 1]))
    # the exposures of a visit side by side, left of the visit mean
    rank = np.zeros(data.n)
    for ss in range(data.nseq):
        idx = np.where(data.seq == ss)[0]
        idx = idx[np.argsort(data.time[idx])]
        rank[idx] = np.linspace(-0.35, 0.05, len(idx)) if len(idx) > 1 \
            else -0.15
    cmap = _scatter_reliability(ax1, data.seq + rank, value, data.err, prob,
                                size=9)
    ax1.errorbar(np.arange(data.nseq) + 0.25, seqm['value'], seqm['err'],
                 fmt='_', color=C['binned'], ms=7, mew=1.4, elinewidth=1.0,
                 zorder=5, label='visit mean')
    ax1.set_ylabel('residual [m s$^{-1}$]' if fit is not None
                   else 'RV [m s$^{-1}$]')
    ax1.legend(loc='lower left')
    if title:
        ax1.set_title(title, loc='left')
    colors = cmap(seqm['prob'])
    ax2.bar(np.arange(data.nseq), seqm['prob'], width=0.8, color=colors,
            edgecolor=C['surface'], linewidth=1.0)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel('P(outlier)')
    ax2.set_xlabel('visit (in time order)')
    if highlight is not None:
        flagged = np.where(seqm['prob'] > highlight)[0]
        for ss in flagged:
            # upright inside the bar: neighbours cannot collide, and the
            #   labels keep the order of the visits
            ax2.annotate(f'{seqm["time"][ss]:.1f}', (ss, seqm['prob'][ss]),
                         xytext=(0, -3), textcoords='offset points',
                         rotation=90, ha='center', va='top', fontsize=6,
                         color=C['text'], annotation_clip=True)
    ax2.set_xlim(-1, data.nseq)
    _colorbar(fig, [ax1, ax2], cmap)
    return fig


# =============================================================================
# Periodograms
# =============================================================================
def _decimate(xval: np.ndarray, yval: np.ndarray, npix: int = 1500):
    """
    A periodogram thinned for drawing, every peak kept

    The curve is cut into about npix bins along x, and each bin keeps its
    lowest and highest points, in order: the drawn envelope is the same,
    the file is a fraction of the size.

    :param xval: np.ndarray, x (monotonic)
    :param yval: np.ndarray, y
    :param npix: int, the number of bins

    :return: tuple, the thinned x and y
    """
    if len(xval) <= 2 * npix:
        return xval, yval
    edges = np.linspace(0, len(xval), npix + 1).astype(int)
    keep = []
    for start, stop in zip(edges[:-1], edges[1:]):
        seg = yval[start:stop]
        if len(seg) == 0:
            continue
        pair = sorted({start + int(np.argmin(seg)), start + int(np.argmax(seg))})
        keep.extend(pair)
    keep = np.array(keep)
    return xval[keep], yval[keep]


def _mark_period(ax, period, label=None, color=None):
    """A thin vertical line at a period"""
    ax.axvline(period, color=color or C['muted'], lw=0.8, ls=(0, (2, 2)),
               zorder=0)
    if label:
        ax.annotate(label, (period, 1), xycoords=('data', 'axes fraction'),
                    xytext=(2, -2), textcoords='offset points', va='top',
                    fontsize=7, color=C['text'])


def plain_log_ticks(ax, axis: str = 'x'):
    """
    Plain numbers on a log axis (1, 2, 5, 10, 20...) instead of 10^0, 10^1

    :param ax: matplotlib axes
    :param axis: str, 'x', 'y' or 'both'
    """
    from matplotlib.ticker import FuncFormatter, LogLocator
    major = FuncFormatter(lambda val, _: f'{val:g}')

    def minor(val, _):
        text = f'{val:g}'
        return text if text.lstrip('0.').startswith(('2', '5')) else ''
    for name in (('x', 'y') if axis == 'both' else (axis,)):
        sub = getattr(ax, f'{name}axis')
        sub.set_major_locator(LogLocator(base=10))
        sub.set_major_formatter(major)
        sub.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
        sub.set_minor_formatter(FuncFormatter(minor))


def fip_family(res: Any, marks: Sequence[float] = (),
               threshold: float = 0.01, title: Optional[str] = None):
    """
    The FIP of the period OR any of its aliases at every period (thick, what
    decides on a planet), and of the period alone (thin grey: which alias
    holds the probability), with the threshold and the periods marked

    :param res: FIPResult (with family)
    :param marks: list of float, periods marked (dashed), e.g. known planets
    :param threshold: float, the FIP of a detection (dotted)
    :param title: str or None

    :return: the figure
    """
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    period = 1.0 / np.asarray(res.freq)
    floor = 1e-12
    alone = -np.log10(np.clip(res.fip, floor, 1.0))
    ax.plot(*_decimate(period, alone), color=C['muted'], lw=0.6,
            label='the period alone')
    if getattr(res, 'family', None) is not None:
        fam = -np.log10(np.clip(res.family, floor, 1.0))
        ax.plot(*_decimate(period, fam), color=C['koloa'], lw=1.2,
                label='the period or any of its aliases')
    ax.axhline(-np.log10(threshold), color=C['text'], lw=0.8,
               ls=(0, (1, 2)), label=f'FIP = {threshold:g}')
    for per in marks:
        ax.axvline(per, color=C['outlier'], lw=0.8, ls=(0, (4, 2)))
    ax.set_xscale('log')
    plain_log_ticks(ax, 'x')
    ax.set_xlabel('period [d]')
    ax.set_ylabel('$-\\log_{10}$ FIP')
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', fontsize=7)
    if title:
        ax.set_title(title, loc='left')
    fig.tight_layout()
    return fig


def periodograms(freq: np.ndarray, gls_power: Optional[np.ndarray] = None,
                 oap_gauss: Optional[np.ndarray] = None,
                 oap_mix: Optional[np.ndarray] = None,
                 window_power: Optional[np.ndarray] = None,
                 fips: Optional[Dict[str, Any]] = None,
                 mark: Sequence[float] = (), title: Optional[str] = None,
                 xlim: Optional[Sequence[float]] = None):
    """
    Periodograms stacked on a shared period axis

    - the profile likelihood gain of a sinusoid, gaussian and outlier-aware
      (the same axis, Delta ln L, so they read against each other);
    - the window of the sampling;
    - one panel of -log10 FIP per method, with the 1 % and 0.1 % lines.

    :param freq: np.ndarray, the frequencies [1/day]
    :param gls_power: np.ndarray or None, the GLS power (drawn if no OAP)
    :param oap_gauss: np.ndarray or None, Delta ln L, gaussian
    :param oap_mix: np.ndarray or None, Delta ln L, outlier-aware
    :param window_power: np.ndarray or None, the window
    :param fips: dict or None, method key (koloa, gaussian, soft, hard):
                 FIPResult
    :param mark: list of float, periods to mark [days]
    :param title: str or None
    :param xlim: tuple or None, the period range [days]

    :return: the figure
    """
    period = 1 / freq
    fips = fips or {}
    panels = []
    if oap_gauss is not None or oap_mix is not None or gls_power is not None:
        panels.append('power')
    if window_power is not None:
        panels.append('window')
    order = [key for key in ('gaussian', 'soft', 'hard', 'koloa')
             if key in fips]
    panels += [f'fip:{key}' for key in order]
    heights = [1.6 if pn == 'power' else (0.7 if pn == 'window' else 1.0)
               for pn in panels]
    fig, axes = plt.subplots(len(panels), 1, figsize=(7.2, 1.25 *
                                                      sum(heights) + 0.6),
                             sharex=True,
                             gridspec_kw=dict(height_ratios=heights))
    axes = np.atleast_1d(axes)
    for ax, pn in zip(axes, panels):
        if pn == 'power':
            if oap_gauss is not None or oap_mix is not None:
                if oap_mix is not None:
                    ax.plot(*_decimate(period, oap_mix), color=C['koloa'],
                            lw=0.9, label='outlier-aware (OAP)')
                if oap_gauss is not None:
                    ax.plot(*_decimate(period, oap_gauss),
                            color=C['gaussian'], lw=0.8, alpha=0.9,
                            label='Gaussian')
                ax.set_ylabel('$\\Delta\\ln L$')
                peak = max(float(np.nanmax(arr)) for arr in (oap_mix, oap_gauss)
                           if arr is not None)
            else:
                ax.plot(*_decimate(period, gls_power), color=C['gaussian'],
                        lw=0.9, label='GLS')
                ax.set_ylabel('GLS power')
                peak = float(np.nanmax(gls_power))
            # the legend sits in room kept above the highest peak
            ax.set_ylim(0, 1.35 * max(peak, 1e-6))
            ax.legend(loc='upper left', ncol=2)
        elif pn == 'window':
            ax.plot(*_decimate(period, window_power), color=C['muted'],
                    lw=0.8)
            ax.set_ylabel('window')
            ax.set_ylim(0, max(0.05, 1.1 * np.max(window_power[period < (
                0.5 * np.max(period))])))
        else:
            key = pn.split(':')[1]
            res = fips[key]
            per, score = _decimate(res.period, -res.log10fip)
            ax.fill_between(per, 0, score, color=C[key], alpha=0.18, lw=0)
            ax.plot(per, score, color=C[key], lw=0.9)
            for thr in (2, 3):
                ax.axhline(thr, color=C['threshold'], lw=0.6,
                           ls=(0, (1, 2)), zorder=0)
            # the label of the method sits above the 0.1 % line
            top = max(4.6, float(np.nanmax(score)) * 1.35)
            ax.set_ylim(0, min(top, 18))
            ax.set_ylabel('$-\\log_{10}$FIP')
            best = res.best()
            label = METHOD_LABELS.get(key, key)
            if best:
                label += (f'   best: {best["period"]:.3f} d, FIP = '
                          f'{best["fip"]:.1e}')
            ax.text(0.01, 0.95, label, transform=ax.transAxes, va='top',
                    fontsize=7.5, color=C['text'],
                    bbox=dict(boxstyle='round,pad=0.2', fc=C['surface'],
                              ec='none', alpha=0.8))
    for ax in axes:
        ax.set_xscale('log')
        for per in mark:
            _mark_period(ax, per)
    if xlim is not None:
        axes[-1].set_xlim(*xlim)
    else:
        axes[-1].set_xlim(period.min(), period.max())
    axes[-1].set_xlabel('period [d]')
    if title:
        axes[0].set_title(title, loc='left')
    fig.align_ylabels(axes)
    return fig


def jackknife(freq: np.ndarray, jack: Dict[str, np.ndarray], period: float,
              data: RVData, unit: str = 'sequence',
              title: Optional[str] = None):
    """
    The GLS with every unit left out in turn, and who holds the peak up

    :param freq: np.ndarray, the frequencies [1/day]
    :param jack: dict, the output of koloa.periodogram.jackknife
    :param period: float, the peak to examine [days]
    :param data: RVData, the series (for the times of the units)
    :param unit: str, sequence or point
    :param title: str or None

    :return: the figure
    """
    from koloa.periodogram import influence
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 4.0),
                                   gridspec_kw=dict(height_ratios=[1.4, 1]))
    per = 1 / freq
    plo, low = _decimate(per, jack['low'])
    phi, high = _decimate(per, jack['high'])
    ax1.fill_between(phi, np.interp(phi, plo[::-1], low[::-1]), high,
                     color=C['gaussian'], alpha=0.3, lw=0,
                     label='leave-one-out range')
    ax1.plot(*_decimate(per, jack['power']), color=C['gaussian'], lw=0.8,
             label='all units')
    ax1.set_xscale('log')
    ax1.set_xlim(per.min(), per.max())
    ax1.set_ylabel('GLS power')
    ax1.set_xlabel('period [d]')
    _mark_period(ax1, period, f'{period:.3f} d')
    ax1.legend(loc='upper left')
    infl = influence(jack, freq, period)
    if unit == 'sequence':
        tunit = np.array([np.mean(data.time[data.seq == ss])
                          for ss in range(len(infl))])
    else:
        tunit = data.time
    colors = np.where(infl >= 0, C['koloa'], C['outlier'])
    ax2.bar(np.arange(len(infl)), infl, color=colors, width=0.8,
            edgecolor=C['surface'], linewidth=0.8)
    ax2.axhline(0, color=C['muted'], lw=0.6)
    ax2.set_ylabel(f'power held up\nat {period:.2f} d')
    ax2.set_xlabel(f'{unit} (in time order)')
    worst = np.argsort(-np.abs(infl))[:3]
    for ww in worst:
        ax2.annotate(f'{tunit[ww]:.1f}', (ww, infl[ww]), xytext=(0, 3 if
                     infl[ww] >= 0 else -9), textcoords='offset points',
                     ha='center', fontsize=6.5, color=C['text'])
    if title:
        ax1.set_title(title, loc='left')
    return fig


def coherence(coh: Dict[str, Any], title: Optional[str] = None):
    """
    Amplitude and phase of a fixed-period sinusoid, chunk by chunk

    A planet keeps both; rotation-modulated activity, whose spots come and
    go, does not.

    :param coh: dict, the output of koloa.diagnostics.coherence
    :param title: str or None

    :return: the figure
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.6))
    chunks = coh['chunks']
    tmid = np.array([ch['tmid'] for ch in chunks])
    amp = np.array([ch['K'] for ch in chunks])
    samp = np.array([ch['sK'] for ch in chunks])
    phs = np.degrees(np.array([ch['phase'] for ch in chunks]))
    sphs = np.degrees(np.array([ch['sphase'] for ch in chunks]))
    ax1.errorbar(tmid, amp, samp, fmt='o', color=C['koloa'], ms=5,
                 elinewidth=1.2)
    ax1.axhline(coh['K_all'], color=C['muted'], lw=0.8, ls=(0, (3, 2)),
                label=f'all data: {coh["K_all"]:.2f} m s$^{{-1}}$')
    ax1.set_ylabel('K [m s$^{-1}$]')
    ax1.set_xlabel('time [BJD - 2400000]')
    ax1.set_ylim(bottom=0)
    ax1.legend(loc='best')
    ax2.errorbar(tmid, phs, sphs, fmt='o', color=C['koloa'], ms=5,
                 elinewidth=1.2)
    ax2.set_ylabel('phase [deg]')
    ax2.set_xlabel('time [BJD - 2400000]')
    txt = (f'P = {coh["period"]:.4f} d\n'
           f'amplitude $\\chi^2$ p = {coh["p_amplitude"]:.1e}\n'
           f'vector $\\chi^2$ p = {coh["p_vector"]:.1e}')
    ax2.text(0.02, 0.96, txt, transform=ax2.transAxes, va='top', fontsize=7.5)
    if title:
        ax1.set_title(title, loc='left')
    fig.tight_layout()
    return fig


#: the posterior mass inside the one and two sigma contours of a 2-D
#: gaussian: the contour levels of the corner plots
CORNER_LEVELS = (1 - np.exp(-0.5), 1 - np.exp(-2.0))


def _density_levels(hist: np.ndarray, masses: Sequence[float]) -> np.ndarray:
    """The heights of a 2-D histogram above which given masses lie"""
    flat = np.sort(hist.ravel())[::-1]
    cum = np.cumsum(flat)
    if cum[-1] <= 0:
        return np.zeros(len(masses))
    cum /= cum[-1]
    return np.array([flat[min(int(np.searchsorted(cum, mass)),
                               len(flat) - 1)] for mass in masses])


def _post_color(post: Dict[str, Any]) -> str:
    """The colour of a posterior: a method of the style, or a colour"""
    color = post.get('color', 'koloa')
    return C.get(color, color)


def corner(posteriors: Sequence[Dict[str, Any]], names: Sequence[str],
           labels: Optional[Dict[str, str]] = None,
           truths: Optional[Dict[str, float]] = None,
           ranges: Optional[Dict[str, Sequence[float]]] = None,
           bins: int = 40, smooth: float = 1.2, size: float = 1.0,
           title: Optional[str] = None):
    """
    Several posteriors on one corner plot, each in its own colour

    The diagonal shows the marginal density of every posterior; the other
    panels show the contours that hold 39 % and 86 % of each posterior (the
    one and two sigma of a 2-D gaussian), filled unless the posterior is a
    reference drawn dashed. True values are thin lines.

    :param posteriors: list of dict, one per posterior: samples (dict,
                       name: draws, as FitResult.posterior gives), label,
                       color (a method of the style, such as koloa,
                       gaussian or neutral, or any colour), and optionally
                       dashed (bool, contour lines only, dashed); a
                       posterior may hold only some of the names, and is
                       drawn in the panels of those
    :param names: list of str, the quantities to show, in order
    :param labels: dict or None, the axis label of each name
    :param truths: dict or None, the true value of each name
    :param ranges: dict or None, the limits of each name (default: the
                   central 99 % of every posterior, and the truth)
    :param bins: int, the bins per axis
    :param smooth: float, the gaussian smoothing of the histograms [bins]
    :param size: float, the size of one panel [inches]
    :param title: str or None, the title of the legend

    :return: the figure
    """
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    from scipy.ndimage import gaussian_filter, gaussian_filter1d
    labels = labels or {}
    truths = truths or {}
    ranges = dict(ranges or {})
    nn = len(names)
    for name in names:
        if name in ranges:
            continue
        lows, highs = [], []
        for post in posteriors:
            if name not in post['samples']:
                continue
            val = np.asarray(post['samples'][name], dtype=float)
            low, high = np.percentile(val[np.isfinite(val)], [0.5, 99.5])
            lows.append(low)
            highs.append(high)
        if name in truths:
            lows.append(truths[name])
            highs.append(truths[name])
        low, high = min(lows), max(highs)
        pad = 0.05 * (high - low) if high > low else 1e-3 * max(abs(low), 1)
        ranges[name] = (low - pad, high + pad)
    fig, axes = plt.subplots(nn, nn, figsize=(size * nn + 0.5,
                                              size * nn + 0.5), squeeze=False)
    for row in range(nn):
        for col in range(nn):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue
            xname = names[col]
            ax.set_xlim(ranges[xname])
            if row == col:
                # the marginal densities
                top = 0.0
                for post in posteriors:
                    if xname not in post['samples']:
                        continue
                    color = _post_color(post)
                    val = np.asarray(post['samples'][xname], dtype=float)
                    hist, edges = np.histogram(val, bins,
                                               range=ranges[xname],
                                               density=True)
                    hist = gaussian_filter1d(hist, smooth)
                    mid = 0.5 * (edges[1:] + edges[:-1])
                    dashed = post.get('dashed', False)
                    ax.plot(mid, hist, color=color, lw=1.4,
                            ls=(0, (4, 2)) if dashed else '-')
                    if not dashed:
                        ax.fill_between(mid, 0, hist, color=color, alpha=0.16,
                                        lw=0)
                    top = max(top, float(np.max(hist)))
                ax.set_ylim(0, 1.1 * top if top > 0 else 1)
                ax.set_yticks([])
                ax.grid(False)
                if xname in truths:
                    ax.axvline(truths[xname], color=C['text'], lw=0.8,
                               alpha=0.8)
            else:
                yname = names[row]
                ax.set_ylim(ranges[yname])
                for post in posteriors:
                    if (xname not in post['samples']
                            or yname not in post['samples']):
                        continue
                    color = _post_color(post)
                    xval = np.asarray(post['samples'][xname], dtype=float)
                    yval = np.asarray(post['samples'][yname], dtype=float)
                    hist, xedge, yedge = np.histogram2d(
                        xval, yval, bins, range=[ranges[xname],
                                                 ranges[yname]])
                    hist = gaussian_filter(hist, smooth)
                    inner, outer = _density_levels(hist, CORNER_LEVELS)
                    xmid = 0.5 * (xedge[1:] + xedge[:-1])
                    ymid = 0.5 * (yedge[1:] + yedge[:-1])
                    if not outer < inner < hist.max():
                        continue
                    dashed = post.get('dashed', False)
                    if not dashed:
                        ax.contourf(xmid, ymid, hist.T,
                                    levels=[inner, hist.max() * 1.01],
                                    colors=[color], alpha=0.28)
                        ax.contourf(xmid, ymid, hist.T,
                                    levels=[outer, inner], colors=[color],
                                    alpha=0.10)
                    ax.contour(xmid, ymid, hist.T, levels=[outer, inner],
                               colors=[color], linewidths=[0.9, 1.3],
                               linestyles='--' if dashed else '-')
                if xname in truths:
                    ax.axvline(truths[xname], color=C['text'], lw=0.8,
                               alpha=0.8)
                if yname in truths:
                    ax.axhline(truths[yname], color=C['text'], lw=0.8,
                               alpha=0.8)
                if xname in truths and yname in truths:
                    ax.plot(truths[xname], truths[yname], 's',
                            color=C['text'], ms=3)
                ax.yaxis.set_major_locator(MaxNLocator(3))
            ax.xaxis.set_major_locator(MaxNLocator(3))
            ax.tick_params(labelsize=7)
            try:
                ax.ticklabel_format(useOffset=False)
            except AttributeError:
                pass
            if row < nn - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(labels.get(xname, xname))
                for tick in ax.get_xticklabels():
                    tick.set_rotation(40)
                    tick.set_ha('right')
            if col > 0 or row == 0:
                if row != col:
                    ax.set_yticklabels([])
            else:
                ax.set_ylabel(labels.get(names[row], names[row]))
    handles = [Line2D([], [], color=_post_color(post), lw=2,
                      ls=(0, (4, 2)) if post.get('dashed', False) else '-')
               for post in posteriors]
    texts = [post['label'] for post in posteriors]
    if truths:
        handles.append(Line2D([], [], color=C['text'], lw=0.8, marker='s',
                              ms=3))
        texts.append('true value')
    fig.legend(handles, texts, loc='upper right',
               bbox_to_anchor=(0.98, 0.98), title=title, fontsize=8.5,
               title_fontsize=9, alignment='left')
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.98,
                        wspace=0.07, hspace=0.07)
    fig.align_ylabels(axes[:, 0])
    fig.align_xlabels(axes[-1, :])
    return fig


def _sequential_cmap() -> LinearSegmentedColormap:
    """One hue, from the background to strong koloa blue (magnitude)"""
    return LinearSegmentedColormap.from_list('koloa_seq', SEQUENTIAL_STOPS)


def recovery_map(rmap: Any, title: Optional[str] = None, ax: Any = None,
                 colorbar: bool = True):
    """
    The recovered fraction on the grid of a koloa.completeness.RecoveryMap,
    with its 50 and 90 % contours

    :param rmap: RecoveryMap
    :param title: str or None
    :param ax: matplotlib axes or None (a new figure)
    :param colorbar: bool, draw the colour bar

    :return: the figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.2))
    else:
        fig = ax.figure
    cmap = _sequential_cmap()
    mesh = ax.pcolormesh(rmap.period_edges, rmap.amp_edges, rmap.rate.T,
                         cmap=cmap, vmin=0, vmax=1, shading='flat')
    xmid, ymid = rmap.periods, rmap.amplitudes
    if np.sum(np.isfinite(rmap.rate)) > 4:
        cont = ax.contour(xmid, ymid, rmap.rate.T, levels=[0.5, 0.9],
                          colors=[C['text']], linewidths=[1.0, 1.4],
                          linestyles=['--', '-'])
        ax.clabel(cont, fmt={0.5: '50 %', 0.9: '90 %'}, fontsize=7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('period [d]')
    ax.set_ylabel('$K$ [m s$^{-1}$]')
    ax.grid(False)
    if title:
        ax.set_title(title, loc='left')
    if np.isfinite(rmap.false_alarm):
        ax.text(0.02, 0.97, f'false alarms at $K = 0$: '
                f'{100 * rmap.false_alarm:.1f} %', transform=ax.transAxes,
                va='top', fontsize=7,
                bbox=dict(boxstyle='round,pad=0.25', fc=C['surface'],
                          ec='none', alpha=0.85))
    if colorbar:
        cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
        cbar.set_label('recovered fraction')
        cbar.outline.set_visible(False)
    return fig


def detection_limits(maps: Sequence[Any], labels: Sequence[str],
                     colors: Sequence[str], level: float = 0.5,
                     title: Optional[str] = None):
    """
    The semi-amplitude recovered with a given probability, against period,
    for several recovery maps

    :param maps: list of RecoveryMap
    :param labels: list of str
    :param colors: list of str, methods of the style (koloa, gaussian, ...)
                   or colours
    :param level: float, the recovered fraction
    :param title: str or None

    :return: the figure
    """
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    for rmap, label, color in zip(maps, labels, colors):
        dashed = label.endswith('(no outliers)')
        ax.plot(rmap.periods, rmap.k_at(level), marker='o', ms=3.5,
                color=C.get(color, color), lw=1.5,
                ls='--' if dashed else '-', label=label)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # a labelled tick below the lowest curve
    lowest = np.nanmin([np.nanmin(rmap.k_at(level)) for rmap in maps])
    if np.isfinite(lowest) and lowest > 1:
        ax.set_ylim(bottom=1.0)
    # plain numbers on the log axes, not powers of ten
    from matplotlib.ticker import FuncFormatter
    plain = FuncFormatter(lambda val, _: f'{val:g}')
    ax.yaxis.set_major_formatter(plain)
    ax.yaxis.set_minor_formatter(plain)
    ax.xaxis.set_major_formatter(plain)
    ax.tick_params(axis='y', which='minor', labelsize=7)
    ax.set_xlabel('period [d]')
    ax.set_ylabel(f'$K$ recovered {100 * level:.0f} % of the time '
                  '[m s$^{-1}$]')
    # below the axes, where it hides no curve
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2,
              fontsize=7, frameon=False)
    if title:
        ax.set_title(title, loc='left')
    fig.tight_layout()
    return fig



# =============================================================================
# Why the outliers are outliers (koloa.outliers)
# =============================================================================
#: the diverging steps of a deviation: low (teal), none (grey), high
#: (orange); blue and red stay for good data and outliers
DIVERGING_STOPS = ['#008a9e', '#f1f0ec', '#d0661c']


def outlier_keys(report: Any, inst: Optional[str] = None, maxkeys: int = 12,
                 zmax: float = 6.0, title: Optional[str] = None) -> Any:
    """
    Why the outliers are outliers (koloa.outliers.explain): each outlier (a
    row: a whole visit, or a single exposure) against each key (a column),
    coloured by how far it is from the good data (z: low in teal, high in
    orange, grey near zero), its z written in the cell beyond 2 sigma and
    framed when it is significantly off (a text key: 'rare' when its value
    is). The last row is every outlier together (the mean z, clipped at 3,
    of the association test), starred when the outliers share the key.

    :param report: OutlierReport
    :param inst: str or None, the instrument (the one with the most outliers
                 when None)
    :param maxkeys: int, the most keys shown: those off somewhere first, then
                    the most deviant
    :param zmax: float, the z of the strongest colour
    :param title: str or None

    :return: matplotlib figure
    """
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle
    if inst is None:
        counts = {name: sum(unit['inst'] == name for unit in report.units)
                  for name in report.counts}
        inst = max(counts, key=counts.get) if counts else None
    units = [unit for unit in report.units if unit['inst'] == inst]
    if not units:
        fig, ax = plt.subplots(figsize=(6.0, 1.2))
        ax.axis('off')
        ax.text(0.5, 0.5, f'{report.name}: no outlier', ha='center',
                va='center', color=C['muted'])
        return fig
    # the keys: off somewhere first, then the most deviant anywhere
    score = {}
    for unit in units:
        for key in unit['keys']:
            if key.get('text'):
                val = 1.0 - key['frequency']
            else:
                val = abs(key['z']) if key['z'] is not None and \
                    not np.isnan(key['z']) else 0.0
            old = score.get(key['key'], (False, 0.0, key['label']))
            score[key['key']] = (old[0] or key['significant'],
                                 max(old[1], val), key['label'])
    shared = {row['key']: row for row in report.association
              if row['inst'] == inst}
    names = sorted(score, key=lambda name: (not score[name][0],
                                            not shared.get(name, {}).get(
                                                'significant', False),
                                            -score[name][1]))[:maxkeys]
    # the row of every outlier together, when there was a test of them
    together = bool(shared)
    nrow, ncol = len(units) + int(together), len(names)
    zmat = np.full((nrow, ncol), np.nan)
    cmap = LinearSegmentedColormap.from_list('koloa_diverging',
                                             DIVERGING_STOPS)
    norm = TwoSlopeNorm(vmin=-zmax, vcenter=0.0, vmax=zmax)
    fig, ax = plt.subplots(figsize=(max(6.0, 0.62 * ncol + 2.6),
                                    0.36 * nrow + 1.9))
    notes = []
    for irow, unit in enumerate(units):
        keys = {key['key']: key for key in unit['keys']}
        for icol, name in enumerate(names):
            key = keys.get(name)
            if key is None:
                continue
            if key.get('text'):
                notes.append((irow, icol, 'rare' if key['significant']
                              else '', key['significant']))
                continue
            zval = key['z']
            if zval is None or np.isnan(zval):
                continue
            zmat[irow, icol] = np.clip(zval, -zmax, zmax)
            label = ('+inf' if zval == np.inf else '-inf' if zval == -np.inf
                     else f'{zval:+.0f}' if abs(zval) >= 9.5 else
                     f'{zval:+.1f}')
            notes.append((irow, icol, label if abs(zval) >= 2 else '',
                          key['significant']))
    # every outlier together
    for icol, name in enumerate(names if together else []):
        row = shared.get(name)
        if row is None:
            continue
        if not row.get('text'):
            zmat[-1, icol] = np.clip(row['mean_z'], -zmax, zmax)
            notes.append((nrow - 1, icol, f'{row["mean_z"]:+.1f}' +
                          ('*' if row['significant'] else ''), False))
        elif row['significant']:
            notes.append((nrow - 1, icol, 'shared*', False))
    ax.imshow(zmat, cmap=cmap, norm=norm, aspect='auto',
              interpolation='nearest')
    for irow, icol, label, framed in notes:
        zval = zmat[irow, icol]
        strong = np.isfinite(zval) and abs(zval) > 0.6 * zmax
        if label:
            ax.text(icol, irow, label, ha='center', va='center', fontsize=7,
                    color=C['surface'] if strong else C['text'])
        if framed:
            ax.add_patch(Rectangle((icol - 0.5, irow - 0.5), 1, 1,
                                   fill=False, edgecolor=C['text'], lw=1.4))
    # a gap between the outliers and all of them together
    if together:
        ax.axhline(nrow - 1.5, color=C['surface'], lw=3)
    rows = [f'{unit["date"]}  ' + (f'visit ({unit["n"]})'
                                   if unit['kind'] == 'visit' else 'exposure')
            for unit in units] + (['all outliers (mean z)'] if together
                                  else [])
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(rows, fontsize=7.5)
    ax.set_xticks(range(ncol))
    ax.set_xticklabels([_short_key(name) for name in names], rotation=40,
                       ha='right', fontsize=7.5)
    ax.tick_params(length=0)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label('deviation from the good data, z', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    from koloa.outliers import NAMES
    shown = inst if inst != 'inst' else NAMES.get(report.lists.get(inst),
                                                  'the series')
    ax.set_title(title or f'{report.name} ({shown}): why the outliers are '
                          f'outliers (framed: significantly off)',
                 loc='left', fontsize=9)
    fig.tight_layout()
    return fig


def _short_key(name: str) -> str:
    """a key name short enough for an axis (ESO's long keywords cut down)"""
    if name.startswith('HIERARCH ESO '):
        name = name[len('HIERARCH ESO '):]
        if name.endswith(' VAL'):
            name = name[:-4]
        if name.endswith(' START'):
            name = name[:-6]
    return name


# =============================================================================
# End of code
# =============================================================================

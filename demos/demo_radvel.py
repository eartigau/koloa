#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo 7: koloa's outliers inside a radvel fit.

A radvel user keeps their model, their priors and their sampler, and swaps
one line: RVLikelihood becomes koloa.radvel_bridge.OutlierRVLikelihood.
The fit then knows that visits can be bad. Shown here on a simulated
planet with bad visits: the plain radvel fit and the outlier-aware one,
from the same starting point, with the same optimiser.

    python demo_radvel.py

Created on 2026-09-27

@author: artigau
"""
import numpy as np
from scipy.optimize import minimize

from _common import template
from koloa.log import log
from koloa.simulate import REALISTIC, simulate


def start_phase(data, period):
    """
    The time of conjunction and amplitude of the best sinusoid at a period,
    so that both fits start from the same sensible point

    For a circular radvel orbit, v = -K sin(2 pi (t - tc) / P).
    """
    tref = float(np.mean(data.time))
    phase = 2 * np.pi * (data.time - tref) / period
    design = np.array([np.cos(phase), np.sin(phase),
                       np.ones(data.n)]).T / data.err[:, None]
    aval, bval, _ = np.linalg.lstsq(design, data.rv / data.err,
                                    rcond=None)[0]
    lam = np.arctan2(-bval, aval)
    return tref + (np.pi / 2 - lam) * period / (2 * np.pi), np.hypot(aval,
                                                                     bval)


def start_period(data):
    """The best peak of the outlier-aware periodogram between 11 and 13.5 d
    (the same start for both fits)"""
    from koloa.periodogram import oap
    freq = np.arange(1 / 13.5, 1 / 11.0, 1 / (20 * data.baseline))
    res = oap(data, freq, unit='sequence')
    return float(1 / freq[np.argmax(res['dlnl'])])


def radvel_setup(data, likelihood_class, period, **kwargs):
    """A one-planet radvel posterior on a series"""
    import radvel
    tconj, amp = start_phase(data, period)
    params = radvel.Parameters(1, basis='per tc secosw sesinw k')
    params['per1'] = radvel.Parameter(value=period)
    params['tc1'] = radvel.Parameter(value=float(tconj))
    params['secosw1'] = radvel.Parameter(value=0.0, vary=False)
    params['sesinw1'] = radvel.Parameter(value=0.0, vary=False)
    params['k1'] = radvel.Parameter(value=float(amp))
    params['dvdt'] = radvel.Parameter(value=0.0, vary=False)
    params['curv'] = radvel.Parameter(value=0.0, vary=False)
    model = radvel.RVModel(params)
    like = likelihood_class(model, data.time, data.rv, data.err, **kwargs)
    like.params['gamma'] = radvel.Parameter(value=0.0, vary=True,
                                            linear=False)
    like.params['jit'] = radvel.Parameter(value=2.0)
    like.vector.dict_to_vector()
    post = radvel.posterior.Posterior(like)
    post.priors += [radvel.prior.HardBounds('jit', 0.0, 50.0),
                    radvel.prior.HardBounds('k1', 0.0, 50.0)]
    return post


def main():
    try:
        import radvel
    except ImportError:
        log('radvel is not installed (pip install radvel)', 'error')
        return
    from koloa.radvel_bridge import OutlierRVLikelihood, outlier_priors
    tpl = template()
    sim = simulate(planets=[dict(P=12.3, K=4.0, tp=tpl.time[0] + 5)],
                   template=tpl, outliers=REALISTIC, visit_jitter=1.5,
                   seed=31)
    data = sim['data']
    log(f'Simulated: P = 12.3 d, K = 4.0 m/s, '
        f'{int(np.sum(sim["outlier_mask"]))} outlying exposures, clear and '
        f'borderline')
    period = start_period(data)
    log(f'  both fits start at P = {period:.3f} d (outlier-aware '
        f'periodogram)')
    for label, cls, kwargs in (
            ('radvel, gaussian', radvel.likelihood.RVLikelihood, {}),
            ('radvel + koloa outliers', OutlierRVLikelihood,
             dict(seq=data.seq, unit='sequence'))):
        post = radvel_setup(data, cls, period, **kwargs)
        if cls is OutlierRVLikelihood:
            post.priors += outlier_priors(post.likelihood)
        res = minimize(post.neglogprob_array, post.get_vary_params(),
                       method='Powell', options=dict(maxiter=20000))
        res = minimize(post.neglogprob_array, res.x, method='Nelder-Mead',
                       options=dict(maxiter=20000, adaptive=True))
        post.set_vary_params(res.x)
        per = post.params['per1'].value
        amp = post.params['k1'].value
        jit = post.params['jit'].value
        line = (f'  {label:26s} P = {per:.3f} d, K = {amp:.2f} m/s, '
                f'jitter = {jit:.2f} m/s')
        if cls is OutlierRVLikelihood:
            prob = post.likelihood.outlier_probability()
            caught = np.sum((prob > 0.5) & sim['outlier_mask'])
            line += (f', {caught} of {np.sum(sim["outlier_mask"])} bad '
                     f'exposures flagged')
        log(line, 'value')


if __name__ == '__main__':
    main()

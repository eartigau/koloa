#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The radvel likelihood with koloa's mixture gives koloa's numbers.

Created on 2026-09-27

@author: artigau
"""
import numpy as np
import pytest

from koloa import kepler
from koloa import log as klog
from koloa.fit import RVModel
from koloa.simulate import simulate

klog.VERBOSE = False
radvel = pytest.importorskip('radvel')


@pytest.mark.parametrize('unit', ['point', 'sequence'])
def test_radvel_mixture_equals_koloa(unit):
    from koloa.radvel_bridge import OutlierRVLikelihood, outlier_priors
    sim = simulate(planets=[dict(P=9.1, K=4.0, e=0.2, omega=0.7, tp=2.0)],
                   outliers=[dict(kind='visit', frac=0.1, amplitude=6.0)],
                   err=1.5, seed=9, nvisits=40, per_visit=2, baseline=400)
    data = sim['data']
    model = RVModel(data, [dict(period=9.1, eccentric=True)],
                    likelihood='mixture', unit=unit, trend=0,
                    seq_jitter=False)
    theta = model.init.copy()
    theta[model.index['xe_0']] = 0.3
    theta[model.index['ye_0']] = 0.2
    theta[model.index['log_jit_inst']] = np.log(0.8)
    period, tperi, ecc, omega, amp = model.orbit(theta, 0)
    params = radvel.Parameters(1, basis='per tp e w k')
    params['per1'] = radvel.Parameter(value=period)
    params['tp1'] = radvel.Parameter(value=tperi)
    params['e1'] = radvel.Parameter(value=ecc)
    params['w1'] = radvel.Parameter(value=omega)
    params['k1'] = radvel.Parameter(value=amp)
    params['dvdt'] = radvel.Parameter(value=0.0, vary=False)
    params['curv'] = radvel.Parameter(value=0.0, vary=False)
    rvmod = radvel.RVModel(params)
    like = OutlierRVLikelihood(rvmod, data.time, data.rv, data.err,
                               seq=data.seq, unit=unit)
    like.params['gamma'] = radvel.Parameter(
        value=theta[model.index['offset_inst']], vary=True, linear=False)
    like.params['jit'] = radvel.Parameter(value=0.8)
    frac, width = model.outlier_params(theta)
    like.params['logit_fout'] = radvel.Parameter(value=np.log(frac / (1 - frac)))
    like.params['log_wout'] = radvel.Parameter(value=np.log(width))
    like.vector.dict_to_vector()
    assert np.isclose(like.logprob(), model.log_likelihood(theta))
    assert np.allclose(like.outlier_probability(),
                       model.outlier_probability(theta))
    post = radvel.posterior.Posterior(like)
    post.priors += outlier_priors(like)
    assert np.isfinite(post.logprob())

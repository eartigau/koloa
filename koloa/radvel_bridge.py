#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa's outlier mixture, inside radvel.

radvel (Fulton et al. 2018) is what much of the field fits orbits with; a
setup file, its priors and its MCMC should not have to be rewritten to
become outlier-aware. `OutlierRVLikelihood` is a drop-in replacement for
radvel.likelihood.RVLikelihood: same model, same residuals, same jitter,
two more parameters (the logit of the outlier fraction and the log of the
outlier width), and a log likelihood in which every point, or every visit,
is summed over being good or an outlier.

    like = OutlierRVLikelihood(model, t, vel, errvel, seq=visit_ids,
                               unit='sequence')
    post = radvel.posterior.Posterior(like)
    post.priors += outlier_priors(like)

radvel is imported only here, and only when this module is used.

Created on 2026-09-27

@author: artigau
"""
from typing import List, Optional

import numpy as np
from scipy.special import expit, logit

from koloa.noise import mixture_loglike

try:
    import radvel
    from radvel.likelihood import RVLikelihood
except ImportError as _exc:  # pragma: no cover
    raise ImportError('koloa.radvel_bridge needs radvel (pip install radvel)'
                      ) from _exc


# =============================================================================
# Define classes
# =============================================================================
class OutlierRVLikelihood(RVLikelihood):
    """
    radvel's RV likelihood with koloa's outlier mixture

    Extra parameters (with the likelihood's suffix):
    - 'logit_fout': logit of the outlier fraction f;
    - 'log_wout': log of the outlier width W [m/s].
    """

    def __init__(self, model, t, vel, errvel, suffix: str = '',
                 seq: Optional[np.ndarray] = None, unit: str = 'point',
                 frac: float = 0.03, width: Optional[float] = None,
                 **kwargs):
        """
        :param model: radvel.model.RVModel
        :param t: np.ndarray, the time [days]
        :param vel: np.ndarray, the velocity [m/s]
        :param errvel: np.ndarray, its error [m/s]
        :param suffix: str, the suffix of this instrument's parameters
        :param seq: np.ndarray or None, the visit of each point (needed for
                    unit='sequence')
        :param unit: str, point or sequence
        :param frac: float, the starting outlier fraction
        :param width: float or None, the starting outlier width (8 times
                      the rms of the velocities when None)
        :param kwargs: passed to RVLikelihood
        """
        super().__init__(model, t, vel, errvel, suffix=suffix, **kwargs)
        self.unit = unit
        if seq is None:
            seq = np.arange(len(t))
        _, self.seq = np.unique(np.asarray(seq), return_inverse=True)
        self.nseq = int(np.max(self.seq)) + 1
        if width is None:
            width = 8.0 * float(np.std(vel))
        self.fout_param = 'logit_fout' + suffix
        self.wout_param = 'log_wout' + suffix
        nvec = self.vector.vector.shape[0]
        for key, val in ((self.fout_param, float(logit(frac))),
                         (self.wout_param, float(np.log(width)))):
            if key not in self.params:
                self.params[key] = radvel.Parameter(value=val)
            if key not in self.vector.indices:
                self.vector.indices.update({key: nvec})
                nvec += 1
        self.vector.dict_to_vector()
        self.vector.vector_names()
        self.fout_index = self.vector.indices[self.fout_param]
        self.wout_index = self.vector.indices[self.wout_param]

    def _outlier_values(self):
        """The fraction and the width"""
        frac = float(expit(self.vector.vector[self.fout_index][0]))
        width = float(np.exp(self.vector.vector[self.wout_index][0]))
        return frac, width

    def logprob(self) -> float:
        """
        The log likelihood, every unit summed over good and outlier

        :return: float, the log likelihood
        """
        jit = self.vector.vector[self.jit_index][0]
        resid = self.residuals()
        diag = self.yerr ** 2 + jit ** 2
        frac, width = self._outlier_values()
        if not (0 < frac < 1) or width <= 0:
            return -np.inf
        per_block, _, _ = mixture_loglike(resid, diag, self.seq, self.nseq,
                                          0.0, frac, width, self.unit)
        return float(np.sum(per_block))

    def outlier_probability(self) -> np.ndarray:
        """
        The probability that each point belongs to an outlier unit, at the
        current parameters

        :return: np.ndarray, (n)
        """
        jit = self.vector.vector[self.jit_index][0]
        diag = self.yerr ** 2 + jit ** 2
        frac, width = self._outlier_values()
        _, logp, _ = mixture_loglike(self.residuals(), diag, self.seq,
                                     self.nseq, 0.0, frac, width, self.unit)
        prob = np.exp(logp)
        return prob[self.seq] if self.unit == 'sequence' else prob


# =============================================================================
# Define functions
# =============================================================================
def outlier_priors(like: OutlierRVLikelihood, frac_prior=(1.0, 20.0),
                   width_range=None) -> List:
    """
    koloa's default priors on the two outlier parameters, as radvel priors

    A beta prior on the fraction (in logit space, the jacobian included)
    and hard bounds on the log width.

    :param like: OutlierRVLikelihood
    :param frac_prior: tuple, the beta prior (a, b) of the fraction
    :param width_range: tuple or None, the bounds of W [m/s] (3 to 300 times
                        the robust rms of the velocities when None)

    :return: list, radvel priors
    """
    apri, bpri = frac_prior
    if width_range is None:
        rstd = 1.4826 * np.median(np.abs(like.y - np.median(like.y)))
        scale = max(rstd, float(np.median(like.yerr)))
        width_range = (3 * scale, 300 * scale)

    def _beta(values):
        uval = values[0]
        return float(-apri * np.logaddexp(0, -uval)
                     - bpri * np.logaddexp(0, uval))

    return [radvel.prior.UserDefinedPrior([like.fout_param], _beta,
                                          'Beta prior on f'),
            radvel.prior.HardBounds(like.wout_param, np.log(width_range[0]),
                                    np.log(width_range[1]))]


# =============================================================================
# End of code
# =============================================================================

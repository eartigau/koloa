#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The noise model: visit blocks and the outlier mixture.

The covariance is a diagonal plus one constant block per visit, solved
block by block and never built. On top of it, each exposure and each visit
is good or an outlier whose variance grows by W^2, with every
configuration summed. Three exposures moved together cost one visit
outlier; a lone spike is an outlier exposure.

Created on 2026-09-27

@author: artigau
"""
# numpy, and koloa's log (timestamped lines, 'value' for a number) and
#   noise model
import numpy as np

from koloa.log import log
from koloa.noise import BlockCov, mixture_loglike_both

rng = np.random.default_rng(4)
# six visits of three exposures: 2 m/s photon noise, a 3 m/s visit jitter
# seq: the visit of each exposure; one jitter draw per visit, indexed by
#   seq so that the three exposures of a visit share it
seq = np.repeat(np.arange(6), 3)
err = np.full(len(seq), 2.0)
resid = 2.0 * rng.normal(size=len(seq)) + 3.0 * rng.normal(size=6)[seq]
resid[seq == 4] += 25.0     # a bad visit: its exposures move together
resid[4] += 20.0            # a spike: one exposure alone, in visit 1

# V = diag(err^2) plus (3 m/s)^2 across each visit: block, the visit of
#   each point, and blockval, the constant of each block [(m/s)^2]
cov = BlockCov(err ** 2, block=seq, blockval=np.full(6, 3.0 ** 2))
# the same matrix built in full, for comparison
dense = np.diag(err ** 2) + 3.0 ** 2 * (seq[:, None] == seq[None, :])
# the log determinant, block by block and dense
log(f'log|V| = {cov.logdet():.6f} (dense: {np.linalg.slogdet(dense)[1]:.6f})',
    'value')
# V^-1 r, block by block (Sherman-Morrison) and dense
diff = np.max(np.abs(cov.solve(resid) - np.linalg.solve(dense, resid)))
log(f'V^-1 r: largest difference from the dense solve {diff:.1e}', 'value')

# every outlier configuration summed: the residual, the noise variance of
#   each exposure [(m/s)^2], the visits and their number, the visit
#   jitter squared, then the fraction and the width W [m/s] of the
#   exposure outliers (_pt) and of the visit outliers (_seq)
# it returns the log likelihood of each visit, and the log probability
#   that each exposure is bad (itself, or in a bad visit) and that each
#   visit is an outlier
per_visit, logp_bad, logp_visit = mixture_loglike_both(
    resid, err ** 2, seq, 6, 3.0 ** 2, frac_pt=0.05, width_pt=30.0,
    frac_seq=0.05, width_seq=30.0)
for vis in range(6):
    sel = seq == vis
    log(f'visit {vis}: {" ".join(f"{val:+6.1f}" for val in resid[sel])} m/s'
        f'  P(visit out) {np.exp(logp_visit[vis]):.2f}  P(exposure bad) '
        f'{" ".join(f"{val:.2f}" for val in np.exp(logp_bad[sel]))}', 'value')
# the visits are independent: the sum is the log likelihood of the series
log(f'log likelihood, every configuration summed: {per_visit.sum():.2f}',
    'value')

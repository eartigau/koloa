#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa: Keplerian Orbits Lifted from Outlier-Affliction.

The koloa is the Hawaiian duck. The package asks of a radial velocity
signal the question of the duck test (does it look, swim and quack like a
planet?) and asks it of data in which some points, or whole visits, are
not to be trusted.

Every tool is outlier-aware: the noise is a mixture in which each point,
or each visit, is good or an outlier, and that indicator is summed over
(or sampled) rather than decided by a clip.

- koloa.data          RVData: velocities, instruments, visits, indicators
- koloa.fip           the outlier-aware false inclusion probability (OAFIP)
- koloa.periodogram   GLS, the outlier-aware periodogram (OAP), jackknife
- koloa.fit           keplerians + GP + outliers, fitted jointly;
                      mcmc_orbits: orbital elements by MCMC
- koloa.diagnostics   coherence, indicators, and the duck test
- koloa.outliers      why an outlier is one: the header keywords, S/N,
                      error bars and indicators that are off for it
- koloa.dace          public velocities of every instrument from DACE
- koloa.tess          the TESS light curves of a star, and whether its
                      brightness varies at a velocity period
- koloa.plotting      figures coloured by reliability
- koloa.simulate      series with known planets, activity and outliers
- koloa.secular       the acceleration of a star with its errors, what it
                      says of a companion, and the perspective acceleration
                      from Gaia
- koloa.analyze       everything, in one call
- koloa.detailed      everything koloa can say about a star from one file:
                      its known planets, more data from DACE, the FIP in
                      two passes, the activity, why each outlier is one

Created on 2026-09-27

@author: artigau
"""
__version__ = '0.1.0'

from koloa.data import RVData, find_sequences, robust_std  # noqa: F401
from koloa.fip import FIPResult, fip_comparison, fip_single, oafip  # noqa
from koloa.fip import inflate_to_fit  # noqa: F401
from koloa.fit import RVModel, FitResult, fit_planets  # noqa: F401
from koloa.fit import mcmc_orbits  # noqa: F401
from koloa.fit import period_prior_from_indicator, sequential_fit  # noqa
from koloa.periodogram import gls, oap, jackknife, window  # noqa: F401
from koloa.diagnostics import coherence, duck_test  # noqa: F401
from koloa.simulate import simulate  # noqa: F401
from koloa.outliers import explain as explain_outliers  # noqa: F401
from koloa.detailed import detailed_analysis  # noqa: F401,E402
from koloa.secular import acceleration, companion_min_mass  # noqa: F401
from koloa.secular import perspective_acceleration, gaia_astrometry  # noqa
from koloa.secular import perspective_from_astrometry  # noqa: F401
from koloa.secular import secular_acceleration  # noqa: F401 (old name)
from koloa.secular import secular_from_astrometry  # noqa: F401 (old name)
# koloa.log is the module (koloa.log.VERBOSE = False silences it); the
#   function is koloa.log.log
from koloa.log import KoloaWarning  # noqa: F401
import koloa.log  # noqa: F401,E402


def analyze(*args, **kwargs):
    """See koloa.analyze.analyze"""
    from koloa.analyze import analyze as _analyze
    return _analyze(*args, **kwargs)

# =============================================================================
# End of code
# =============================================================================

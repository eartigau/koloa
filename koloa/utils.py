#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Small things every module needs.

Created on 2026-09-27

@author: artigau
"""
import contextlib

# =============================================================================
# Define functions
# =============================================================================
@contextlib.contextmanager
def blas_threads(nthreads: int = 1):
    """
    Limit the threads of the linear algebra library inside a block

    koloa's samplers make many products of small matrices. A threaded BLAS
    spends more time waking its threads than computing them: on a 10-core
    machine with OpenBLAS the FIP sampler runs 30 times faster on one
    thread. threadpoolctl does the limiting when it is installed; without
    it, nothing is changed.

    :param nthreads: int, the number of threads inside the block
    """
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        yield
        return
    with threadpool_limits(limits=nthreads, user_api='blas'):
        yield


# =============================================================================
# End of code
# =============================================================================

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The log.

Every line koloa prints reads 'YYMMDD HH:MM:SS.SS | message', and on a
terminal its colour says what it is: green for what is happening, blue for
a number, orange for something skipped or doubtful, red for why a run
stops. The output below was captured without a terminal, so without
colours.

Created on 2026-09-27

@author: artigau
"""
# importlib and warnings from the standard library; koloa's log function,
#   its warning class, and a warning that is also logged
import importlib
import warnings

from koloa.log import KoloaWarning, log, loud_warning

# the module itself, whose switch VERBOSE silences koloa (the log imported
#   above is the function)
klog = importlib.import_module('koloa.log')

# the second argument is the role of the line: 'info' (the default),
#   'value', 'warn' or 'error'
log('reading the series')                            # green
log('median error 3.19 m/s', 'value')                # blue
log('indicator dW has no errors: skipped', 'warn')   # orange
log('no velocity column in the file', 'error')       # red

# silence koloa: warnings and errors still get through
klog.VERBOSE = False
log('this line is not printed')
log('but a warning is', 'warn')
# and the log back on
klog.VERBOSE = True

# a warning for the human (the log) and for the code (a KoloaWarning)
# catch_warnings(record=True) collects the warnings in a list instead of
#   printing them; 'always' keeps every KoloaWarning, repeats included
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always', KoloaWarning)
    # logged in orange, then raised as a KoloaWarning
    loud_warning('a sequential fit biases K low')
# each warning caught holds its category (the class) and its message
log(f'caught {len(caught)} {caught[0].category.__name__}: '
    f'{caught[0].message}', 'value')

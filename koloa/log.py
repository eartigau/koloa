#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The one way koloa talks to the console.

Every line reads 'YYMMDD HH:MM:SS.SS | message' and its colour says what it
is: green for what is happening, blue for a number being reported, orange
for something skipped or doubtful that does not stop the run, red for why
the run stops.

Created on 2026-09-27

@author: artigau
"""
import sys
import warnings
from datetime import datetime

# =============================================================================
# Define variables
# =============================================================================
_LOG_COLORS = {
    'info': '\033[92m',        # green, what is happening
    'value': '\033[94m',       # blue, a number being reported
    'warn': '\033[38;5;208m',  # orange, skipped or doubtful, not fatal
    'error': '\033[91m',       # red, why the run stops
}
_LOG_RESET = '\033[0m'

#: set to False to silence koloa (the tests do)
VERBOSE = True


class KoloaWarning(UserWarning):
    """A warning that koloa wants the user to actually read"""


# =============================================================================
# Define functions
# =============================================================================
def _timestamp() -> str:
    """
    The current instant as 'YYMMDD HH:MM:SS.SS'

    :return: str, the timestamp
    """
    now = datetime.now()
    return now.strftime('%y%m%d %H:%M:%S') + f'.{now.microsecond // 10000:02d}'


def log(message: str, level: str = 'info'):
    """
    Print one status line, timestamped and colour-coded by its role

    :param message: str, what to say
    :param level: str, info, value, warn or error
    """
    if not VERBOSE and level not in ('warn', 'error'):
        return
    color = _LOG_COLORS.get(level, _LOG_COLORS['info'])
    # colours only when a human is reading
    if sys.stdout.isatty():
        print(f'{color}{_timestamp()} | {message}{_LOG_RESET}', flush=True)
    else:
        print(f'{_timestamp()} | {message}', flush=True)


def step(name: str):
    """
    The start of a step of an analysis: one line, 'step: <name>', that
    koloa's GUI reads to show where a run is (and the terminal shows too)

    :param name: str, the step
    """
    log(f'step: {name}', 'info')


def outcome(message: str):
    """
    What a step of an analysis came to (the points an archive gave, say):
    one line, 'result: <message>', that koloa's GUI shows under the step,
    once it is over too

    :param message: str, what it came to
    """
    log(f'result: {message}', 'value')


def loud_warning(message: str):
    """
    A warning that is both logged in orange and raised as a KoloaWarning

    Used for the things a user can get wrong without noticing (a sequential
    GP then keplerian fit, a clip that removed too much weight): the log
    line is for the human, the warning for the code that calls koloa.

    :param message: str, what went wrong and why it matters
    """
    log(message, 'warn')
    warnings.warn(message, KoloaWarning, stacklevel=3)


# =============================================================================
# End of code
# =============================================================================

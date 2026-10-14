#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Where koloa keeps what it fetched once: ~/.cache/koloa, or the folder the
environment variable KOLOA_CACHE names. A batch carried to another machine
takes its cache with it (the tables of the NASA Exoplanet Archive, the
lists of the surveys, APERO's names): its script sets KOLOA_CACHE before
it imports koloa, and nothing is asked of the network there.

And where it keeps your choices (whether to share what it fetches, say):
~/.config/koloa, or the folder KOLOA_CONFIG names.

Created on 2026-10-06

@author: artigau
"""
import os


def cache(*parts: str) -> str:
    """a path in koloa's cache: cache('archive') is ~/.cache/koloa/archive,
    or <KOLOA_CACHE>/archive"""
    base = os.environ.get('KOLOA_CACHE') or os.path.join(
        os.path.expanduser('~'), '.cache', 'koloa')
    return os.path.join(os.path.abspath(os.path.expanduser(base)), *parts)


def config(*parts: str) -> str:
    """a path among koloa's settings: config('share.json') is
    ~/.config/koloa/share.json, or <KOLOA_CONFIG>/share.json"""
    base = os.environ.get('KOLOA_CONFIG') or os.path.join(
        os.path.expanduser('~'), '.config', 'koloa')
    return os.path.join(os.path.abspath(os.path.expanduser(base)), *parts)

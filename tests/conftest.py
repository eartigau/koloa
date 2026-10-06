#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
What every test shares: no copy of APERO's database of names is fetched
(koloa.apero_names reads an empty folder, unless a test gives its own)

Created on 2026-10-03

@author: artigau
"""
import pytest


@pytest.fixture(autouse=True)
def no_apero_names(tmp_path_factory, monkeypatch):
    """APERO's names read from an empty folder: nothing fetched"""
    monkeypatch.setenv('KOLOA_APERO_ASTROMETRICS',
                       str(tmp_path_factory.mktemp('apero_none')))


@pytest.fixture(autouse=True)
def no_encyclopaedia(monkeypatch):
    """exoplanet.eu's catalogue is not fetched (minutes, and the network):
    an empty one, unless a test gives its own"""
    from koloa import archive
    monkeypatch.setattr(archive, '_EU', dict(fetched='a test', planets=[]))


@pytest.fixture(autouse=True)
def no_space_missions(monkeypatch):
    """Kepler, K2 and CoRoT are not asked (MAST, VizieR): a star is in
    none of their fields, unless a test says otherwise"""
    from koloa import gui
    monkeypatch.setattr(gui, 'space_have', lambda target, ident=None:
                        dict(kepler=0, k2=0, corot=0))
    monkeypatch.setattr(gui, 'SPACE_LC', {})

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
What every test shares: no copy of APERO's database of names is fetched
(koloa.apero_names reads an empty folder, unless a test gives its own)

Created on 2026-10-03

@author: artigau
"""
import os

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


@pytest.fixture(autouse=True)
def no_observatory_archives(monkeypatch):
    """the archive of the TNG is not asked: HARPS-N was pointed nowhere,
    unless a test gives its own list of pointings"""
    import numpy as np
    from koloa import tng
    monkeypatch.setattr(tng, 'pointings', lambda refresh=False: dict(
        asked=0.0, object=[], ra=np.zeros(0), dec=np.zeros(0),
        free=np.zeros(0, int), private=np.zeros(0, int)))
    # nor are those of ESO, NExScI, the OHP and the CFHT
    from koloa import observatories
    monkeypatch.setattr(observatories, 'neid_pointings',
                        lambda refresh=False: dict(
                            ra=np.zeros(0), dec=np.zeros(0),
                            n=np.zeros(0, int)))
    monkeypatch.setattr(observatories, 'SOURCES', {
        key: (lambda star, folder, refresh=False, name=None:
              ([], 'not asked in a test'))
        for key in observatories.SOURCES})


@pytest.fixture(autouse=True)
def no_shared_velocities(tmp_path_factory, monkeypatch):
    """what koloa's users shared is not read, nothing is shared, and your
    settings are not those of the tests"""
    from koloa import shared
    monkeypatch.setattr(shared, 'index', lambda refresh=False: [])
    monkeypatch.setattr(shared, 'can_share', lambda: False)
    monkeypatch.setenv('KOLOA_CONFIG', str(tmp_path_factory.mktemp('config')))
    monkeypatch.delenv('KOLOA_SHARE', raising=False)
    monkeypatch.setattr(shared, 'SETTINGS', os.path.join(
        os.environ['KOLOA_CONFIG'], 'share.json'))
    # nor those of LBL, and its Python is asked again in each test
    from koloa import gui, lblrun
    monkeypatch.setattr(lblrun, 'SETTINGS', os.path.join(
        os.environ['KOLOA_CONFIG'], 'lbl.json'))
    monkeypatch.setattr(gui, '_LBL_PYTHON', {})
    # nor is the team's repository of velocities read or written
    from koloa import rvstore
    monkeypatch.setattr(rvstore, 'SETTINGS', os.path.join(
        os.environ['KOLOA_CONFIG'], 'rvs.json'))
    monkeypatch.setattr(rvstore, '_ACCESS', dict(read=False, write=False))


@pytest.fixture(autouse=True)
def no_astrometry(monkeypatch):
    """the archive of Gaia and SIMBAD are not asked how a star moves: no
    astrometry, so no perspective acceleration, unless a test gives its
    own (or the star carries it)"""
    from koloa import secular
    monkeypatch.setattr(secular, 'star_astrometry',
                        lambda ident, timeout=60.0, refresh=False: None)

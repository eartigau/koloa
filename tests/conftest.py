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

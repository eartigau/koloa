#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Everything in one call, and the koloa command.

analyze runs the periodograms, the FIP with every treatment of the
outliers, a fit of the best signal, the duck test and every figure, and
writes a text report and a JSON summary; here on a simulated planet at
11.2 d with bad visits. The koloa command does the same from a shell, on a
file. Small settings here; the defaults take a few minutes.

Created on 2026-09-27

@author: artigau
"""
# os and tempfile for a scratch folder; koloa's analyze, its series
#   (RVData), its log (timestamped lines, 'value' for a number) and its
#   simulation
import os
import tempfile
from pathlib import Path

from koloa.analyze import analyze
from koloa.data import RVData
from koloa.log import log
from koloa.simulate import simulate

# the repository root: this file sits two levels down, in docs/examples
ROOT = Path(__file__).resolve().parents[2]

# the same analysis from a shell takes a csv or .rdb file (pip install -e .
#   puts koloa on the path), for instance the public Kepler-21 series:
#   koloa data/kepler21_harpsn_dace_drs3.3.12.csv --outdir example_analysis \
#         --kmax 1 --unit sequence --nsweep 300 --nburn 100 --nchains 1 \
#         --no-gp
#   (add --name Kepler-21 to name the series and its files, --style web
#   for SVG)

# a circular planet (P = 11.2 d, K = 5 m/s, tp its time of periastron
#   [days]) on a real NIRPS sampling, 7 % of the visits moved as a block
#   by 6 median errors or more, and a 2 m/s jitter shared by each visit
#   (ex_simulate.py shows the recipes)
tpl = RVData.from_csv(ROOT / 'data' / 'nirps_template.csv',
                      name='NIRPS template')
sim = simulate(planets=[dict(P=11.2, K=5.0, tp=60100.0)], template=tpl,
               outliers=[dict(kind='visit', frac=0.07, amplitude=6.0)],
               visit_jitter=2.0, seed=1)
# the series (named 'simulation', which starts every file name)
data = sim['data']
# in a scratch folder, so that the example leaves nothing behind
with tempfile.TemporaryDirectory() as scratch:
    os.chdir(scratch)
    # outdir: the folder for the figures and the reports; kmax: the most
    #   signals the FIP sampler considers; unit='sequence': an outlier is
    #   a whole visit (the default 'both' also allows a lone exposure);
    #   nsweep, nburn, nchains: the FIP chains, short here (defaults 2000,
    #   400 and 2); gp=False skips the slow GP test of the duck test;
    #   style='web' writes SVG figures ('paper', the default: PDF)
    result = analyze(data, outdir='example_analysis', kmax=1,
                     unit='sequence', nsweep=300, nburn=100, nchains=1,
                     gp=False, style='web')
    log(f'written: {", ".join(sorted(os.listdir("example_analysis")))}')
    # out of the scratch folder before it is deleted
    os.chdir(ROOT)
# a dict of every result: the numbers of the series, the periodograms,
#   every FIP, the period examined, its fit, the coherence test and the
#   duck test, and the paths of the figures
log(f'the result holds: {", ".join(result)}')

#!/bin/sh
# Run every demo, in the order the web page and the paper use them.
# Each writes PDF figures to demos/output/<demo>/ and SVG ones to
# docs/figures/, and its log to demos/output/<demo>.log (the paper reads
# some numbers from the logs). The whole set takes about two hours on a
# ten-core laptop (the two posterior ensembles, 40 realisations each, take
# half of it).
set -e
# one thread per process: the demos run in parallel pools, and a BLAS pool
#   per process would oversubscribe the machine many times over
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
mkdir -p output
# a demo that fails (those of GL 725, TOI-2120 and GJ 687 need series of
#   the SPIRou and NIRPS teams, not in this repository) is reported, and
#   the others still run
run() {
    name=$1
    shift
    echo "running $name"
    python3 "$@" > "output/$name.log" 2>&1 \
        || echo "  $name failed: see demos/output/$name.log"
}
run false_alarm demo_false_alarm.py
run recovery demo_recovery.py
run two_planets demo_two_planets.py
run gp_joint demo_gp_joint.py
# the bias and coverage of the joint fit, 30 more realisations (added to
#   the summary of the run above)
run gp_joint_ensemble demo_gp_joint.py --ensemble 30 --workers 9
run completeness demo_completeness.py --workers 9
run contamination demo_contamination.py --workers 9
run secular demo_secular.py
run posterior demo_posterior.py --kamp 10 --nreal 40 --workers 9
run posterior_k6 demo_posterior.py --kamp 6 --nreal 24 --workers 9
run radvel demo_radvel.py
run toi2120 demo_toi2120.py --workers 4 --nreal 6
run gl725b demo_gl725b.py
run kepler21 demo_kepler21.py
run rotation demo_rotation.py --nsteps 16000 --nreal 6 --ncpu 9
run outliers demo_outliers.py
# DACE answers from some networks only: the fetch is kept in data/dace/
run dace demo_dace.py
echo "all demos done"

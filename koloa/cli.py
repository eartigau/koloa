#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
koloa from the command line.

    koloa series.csv --outdir out --kmax 3 --unit sequence
    koloa star.rdb --detailed --outdir star      # the known planets from the
                                                 # archive, more data from
                                                 # DACE and VizieR, the FIP in
                                                 # two passes, why each
                                                 # outlier, a PDF report
    koloa star.rdb --detailed --literature paper_rvs.dat
    koloa star.rdb --detailed --target "GJ 436"  # the SIMBAD name: all the
                                                 # instruments of the archives
    koloa --detailed --target "GJ 436"           # no file: the archives only
    koloa --gui                                  # the same in the browser,
                                                 # with the command lines
    koloanui                                     # koloa --gui, in one word
    koloa "GJ 436" --gather archives             # DACE, CARMENES DR1 and
                                                 # TESS in archives/GJ_436

Created on 2026-09-27

@author: artigau
"""
import argparse

from koloa.analyze import analyze
from koloa.data import RVData
from koloa.detailed import detailed_analysis


# =============================================================================
# Define functions
# =============================================================================
def main(argv=None):
    """The command-line entry point"""
    parser = argparse.ArgumentParser(
        prog='koloa', description='Outlier-aware radial velocity analysis: '
        'periodograms, FIPs, fits and the duck test.')
    parser.add_argument('filename', nargs='*', default=[],
                        help='csv or rdb file (rjd, vrad, svrad and '
                             'indicators with their errors), or with --gather '
                             'the name of a star; with --detailed, several '
                             '(SPIRou and NIRPS, HARPS through LBL...), or '
                             'none for the archives only (--target)')
    parser.add_argument('--instruments', nargs='+', default=None,
                        metavar='NAME',
                        help='with --detailed: the instrument of each file, '
                             'in order (auto: the name the file gives, from '
                             'its columns)')
    parser.add_argument('--outdir', default='koloa_output')
    parser.add_argument('--name', default=None, help='the target name')
    parser.add_argument('--kmax', type=int, default=3,
                        help='the largest number of signals')
    parser.add_argument('--unit', default='both',
                        choices=['both', 'sequence', 'point'],
                        help='what an outlier is')
    parser.add_argument('--pmin', type=float, default=1.1)
    parser.add_argument('--pmax', type=float, default=None)
    parser.add_argument('--nsweep', type=int, default=2000)
    parser.add_argument('--nburn', type=int, default=400)
    parser.add_argument('--nchains', type=int, default=2)
    parser.add_argument('--period', type=float, default=None,
                        help='examine this period instead of the best peak')
    parser.add_argument('--no-duck', action='store_true',
                        help='skip the duck test')
    parser.add_argument('--no-gp', action='store_true',
                        help='skip the GP part of the duck test')
    parser.add_argument('--style', default='paper', choices=['paper', 'web'])
    parser.add_argument('--inst', default=None,
                        help='the instrument column, if several')
    parser.add_argument('--gap', type=float, default=0.3,
                        help='the gap that separates visits [days]')
    parser.add_argument('--detailed', action='store_true',
                        help='the detailed analysis: who the star is and its '
                             'known planets (SIMBAD, NASA Exoplanet Archive), '
                             'more velocities from DACE and VizieR, the FIP '
                             'in two passes, the signals against the known '
                             'planets and the rotation, the activity '
                             'indicators, why each outlier is one, and a '
                             'LaTeX/PDF report')
    parser.add_argument('--target', default=None,
                        help='the SIMBAD name of the star, for the archive, '
                             'DACE, CARMENES and TESS (a guess from the '
                             "file's OBJECT column otherwise); with "
                             '--detailed, give a file, a name, or both')
    parser.add_argument('--dace', action='store_true',
                        help='with --detailed: add the velocities DACE has '
                             'of the star (the archives are used only when '
                             'asked)')
    parser.add_argument('--carmenes', action='store_true',
                        help='with --detailed: add the velocities of '
                             'CARMENES DR1')
    parser.add_argument('--vizier', action='store_true',
                        help='with --detailed: add the velocities published '
                             'on VizieR (the surveys: Keck HIRES, the APF, '
                             'the Lick Hamilton, HARPS by SERVAL; the tables '
                             'of the star\'s papers and of its known '
                             'planets\'), koloa.published')
    parser.add_argument('--no-dace', action='store_true',
                        help='with --gather: not DACE')
    parser.add_argument('--no-archive', action='store_true',
                        help='do not ask the NASA Exoplanet Archive')
    parser.add_argument('--mcmc', action='store_true',
                        help='sample the orbits by MCMC (detailed analysis)')
    parser.add_argument('--literature', nargs='+', default=None,
                        metavar='FILE',
                        help='published velocities to add (detailed '
                             'analysis): VizieR .dat files (time velocity '
                             'error instrument) or csv/rdb files')
    parser.add_argument('--no-vizier', action='store_true',
                        help='with --detailed: no VizieR (the default); with '
                             '--gather: not the published velocities')
    parser.add_argument('--no-latex', action='store_true',
                        help='no LaTeX/PDF report (detailed analysis)')
    parser.add_argument('--site', default=None,
                        help='with --detailed: the observatory of the plans '
                             'that lift an alias (CFHT, La Silla, Paranal, '
                             '...; from the instruments by default)')
    parser.add_argument('--exposures', action='store_true',
                        help='analyse every exposure (koloa runs on the '
                             'nightly means by default)')
    parser.add_argument('--detection-map', action='store_true',
                        help='with --detailed: a detection map, which planets '
                             'the series could have found (injections looked '
                             'for by the FIP with the GP; hours on a long '
                             'series)')
    parser.add_argument('--search-map', action='store_true',
                        help='with --detailed: the quicker detection map of a '
                             'blind periodogram search (no GP; an alias '
                             'counts as missed)')
    parser.add_argument('--no-trend', action='store_true',
                        help='with --detailed: no trend in time fitted with '
                             'the planets (by default the acceleration of the '
                             'star, dv/dt, in the likelihood, reported in '
                             'm/s/yr with its errors)')
    parser.add_argument('--curvature', action='store_true',
                        help='with --detailed: fit the change of the '
                             'acceleration too, d2v/dt2 [m/s/yr^2]')
    parser.add_argument('--toi', nargs='*', default=None, metavar='TOI',
                        help='with --detailed: fit the TESS Objects of '
                             'Interest of the star with the ephemerides of '
                             'TESS (P and a transit as priors): all of them '
                             '(but the false positives) with no number, or '
                             'those given (175.01 175.02)')
    parser.add_argument('--rotation', type=float, default=None,
                        metavar='P',
                        help='with --detailed: a rotation period that can be '
                             'trusted [d]: the GP of the FIP is then an SHO at '
                             'it and at its half, over every period (no bands)')
    parser.add_argument('--fip-gp', default=None,
                        choices=['banded', 'sho', 'none'],
                        help='with --detailed: the GP of the FIP (banded by '
                             'default, sho with --rotation)')
    parser.add_argument('--no-fip-gp', action='store_true',
                        help='with --detailed: no GP of the activity inside '
                             'the FIP')
    parser.add_argument('--no-tess', action='store_true',
                        help='with --detailed or --gather: do not fetch '
                             'the TESS light curves of the star (nor, with '
                             '--gather, those of Kepler, K2 and CoRoT)')
    parser.add_argument('--gather', nargs='?', const='.', default=None,
                        metavar='ROOT',
                        help='gather what DACE, CARMENES DR1, VizieR (the '
                             'published velocities: HIRES, APF, Lick, HARPS '
                             'by SERVAL, the star\'s papers) and TESS have '
                             'of the star named instead of a file, in '
                             'ROOT/<star> (koloa.gather); --no-dace, '
                             '--no-carmenes, --no-vizier, --no-tess leave one '
                             'out, --vizier-sources picks the sources of '
                             'VizieR')
    parser.add_argument('--vizier-sources', nargs='+', default=None,
                        metavar='SOURCE',
                        help='with --gather: the sources of the published '
                             'velocities (all by default): teklu25 (Keck '
                             'HIRES, Teklu et al. 2025), cls21 (HIRES, APF '
                             'and Lick of the California Legacy Survey, '
                             'Rosenthal et al. 2021), talor19 (Keck HIRES, '
                             'Tal-Or et al. 2019), fischer14 (the Lick '
                             'Hamilton, Fischer et al. 2014), rvbank20 '
                             '(HARPS by SERVAL, Trifonov et al. 2020), papers '
                             '(the tables of the star\'s papers)')
    parser.add_argument('--gui', action='store_true',
                        help="koloa's GUI in the browser (koloa.gui): the "
                             'SIMBAD resolver, the velocities by instrument, '
                             'and the command line of every run')
    parser.add_argument('--port', type=int, default=8765,
                        help='with --gui: the port of the page (the next free '
                             'one when taken)')
    parser.add_argument('--exclude', nargs='+', default=None,
                        metavar='INST',
                        help='with --detailed: instruments left out of the '
                             'analysis (NIRPS, HARPS03...), once the file, '
                             'DACE, CARMENES and VizieR are put together')
    parser.add_argument('--refresh-apero', action='store_true',
                        help="a new copy of APERO's database of names "
                             '(koloa.apero_names: the astrometrics of its '
                             'assets, a few MB), by which the star of a file '
                             'is found when --target is not given (its '
                             'OBJECT column, or lbl_<OBJECT>_<TEMPLATE>.rdb)')
    parser.add_argument('--refresh-archive', action='store_true',
                        help='fetch the NASA Exoplanet Archive again (koloa '
                             'keeps it in ~/.cache/koloa/archive)')
    parser.add_argument('--refresh', action='store_true',
                        help='with --gather or --detailed: ask the archives '
                             'again, whatever is already on disk')
    parser.add_argument('--no-carmenes', action='store_true',
                        help='with --gather: not CARMENES DR1')
    parser.add_argument('--periods', nargs='+', type=float, default=None,
                        help='more periods to test (detailed analysis: '
                             'candidates the archive does not list) [d]')
    args = parser.parse_args(argv)
    if args.gui:
        from koloa.gui import serve
        serve(args.port)
        return
    if args.refresh_apero:
        from koloa.apero_names import refresh
        refresh()
        if not args.filename and not args.refresh_archive:
            return
    if args.refresh_archive:
        from koloa.archive import encyclopaedia, tables
        from koloa.gather import carmenes_objects
        from koloa.published import refresh_lists
        tables(refresh=True)
        encyclopaedia(refresh=True)
        carmenes_objects(refresh=True)
        refresh_lists()
        if not args.filename:
            return
    if not args.filename and not (args.detailed and args.target):
        parser.error('give a file (or with --gather a star); --detailed '
                     'takes --target alone too')
    if args.gather is not None:
        from koloa.gather import gather
        gather(args.filename[0], args.gather, dace=not args.no_dace,
               carmenes=not args.no_carmenes, tess=not args.no_tess,
               vizier=not args.no_vizier, refresh=args.refresh,
               vizier_sources=args.vizier_sources)
        return
    if args.detailed:
        detailed_analysis(args.filename or None, outdir=args.outdir,
                          name=args.name, instruments=args.instruments,
                          target=args.target, archive=not args.no_archive,
                          dace=args.dace and not args.no_dace,
                          carmenes=args.carmenes and not args.no_carmenes,
                          exclude=args.exclude, refresh=args.refresh,
                          literature=args.literature,
                          vizier=args.vizier and not args.no_vizier,
                          periods=args.periods,
                          gp=not args.no_gp, kmax=args.kmax,
                          nsweep=args.nsweep, nburn=args.nburn,
                          pmin=args.pmin, pmax=args.pmax, mcmc=args.mcmc,
                          duck=not args.no_duck, style=args.style,
                          latex=not args.no_latex,
                          tess=not args.no_tess, site=args.site,
                          fip_gp=(None if args.no_fip_gp
                                  or args.fip_gp == 'none' else
                                  args.fip_gp or ('sho' if args.rotation
                                                  else 'banded')),
                          toi=(None if args.toi is None else
                               (args.toi or True)),
                          rotation=args.rotation, trend=not args.no_trend,
                          curvature=args.curvature,
                          nightly=False if args.exposures else None,
                          detection_map=('search' if args.search_map else
                                         'fip' if args.detection_map
                                         else False))
        return
    data = RVData.from_csv(args.filename[0], name=args.name, inst=args.inst,
                           sequence_gap=args.gap)
    analyze(data, outdir=args.outdir, kmax=args.kmax, unit=args.unit,
            pmin=args.pmin, pmax=args.pmax, nsweep=args.nsweep,
            nburn=args.nburn, nchains=args.nchains, duck=not args.no_duck,
            gp=not args.no_gp, style=args.style, period=args.period,
            nightly=False if args.exposures else None)


def gui(argv=None):
    """koloanui: koloa's GUI in the browser, as koloa --gui does (with its
    options: koloanui --port 8800)"""
    import sys
    main(['--gui'] + list(sys.argv[1:] if argv is None else argv))


if __name__ == '__main__':
    main()

# =============================================================================
# End of code
# =============================================================================

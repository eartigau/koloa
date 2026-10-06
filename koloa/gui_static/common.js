// koloa's GUI: what its pages share (texts, help, options, the file
// browser, the runs); koloa.gui serves it
'use strict';

const TEXT = {
  en: {
    tagline: 'Outlier-aware radial velocities · on this machine', docs: 'Docs',
    star: 'Star', star_hint: 'Its SIMBAD name: the archives, DACE, CARMENES and TESS find it by that name.',
    resolve: 'Resolve', velocities: 'Velocities',
    velocities_hint: 'Your files of velocities (LBL .rdb, csv, DACE csv), and the archives gathered for the star when DACE or CARMENES DR1 is ticked in the report: what the report will use. Untick an instrument (or click it in the legend) to leave it out of the report.',
    use: 'Used', exclude: 'Instruments left out', source: 'Source', src_file: 'input file',
    include: 'Datasets asked back', ds_why: 'By default', ds_back: 'asked back', ds_sets: 'datasets',
    ds_release: 'left out: the same spectra as {better}{how}',
    ds_more: ', the more precise (this one scatters {extra} m/s more, in quadrature, on the spectra they share)',
    ds_tie: ' (as precise, with more spectra)', ds_err: ' (median errors of {a} against {b} m/s)',
    ds_weak: 'left out: its {nights} nights constrain neither the mean ({mean} % of its error) nor the slope ({slope} %)',
    ds_weak_one: 'left out: one night constrains neither the mean nor the slope', ds_asked: ', asked back',
    ds_weak_part: 'left out: {same} of its points are spectra of {better}{how}; the {left} others constrain neither the mean nor the slope',
    ds_part: '{used} of {n} points: the others are spectra of {better}',
    batch_rules: 'of them, the best release of the same spectra, and not the datasets that constrain nothing',
    file: 'File (optional)', root: 'Archives folder', plot: 'Plot', browse: 'Browse...', picking: 'choosing...',
    var_pending: 'SIMBAD’s periods of variability are on their way (its service is slow today)',
    p_valid: 'P(valid)', p_valid_low: 'P(valid) < 50 %',
    fip_prompt: 'Several instruments: untick those to leave out (fewer, faster), then compute the quick FIP; it can take minutes on a long series.',
    start_fip: 'Compute the quick FIP', stop_fip: 'Stop', fip_stopped: 'stopped',
    envelope: '1σ envelope', published: 'published', known_fold: 'a fold at the published period of this known planet, with its published orbit (dashed)',
    tab_batch: 'Batch FIP', batch_pick: 'Add files...', batch_or: 'or the files of a folder:', batch_add: 'Add them', batch_clear: 'Clear',
    batch_run: 'Run the FIPs', batch_none: 'No file yet.', batch_files: 'files', batch_nomatch: 'No file of that folder matches the pattern.',
    batch_hint: 'The quick FIP of many files, one after the other: each file on its own, or with every archive of its star (found by its APERO name); a table of their objects, best peaks and accelerations; open one in the Analysis tab to look closer.',
    col_file: 'File', col_object: 'Object', col_insts: 'Nights: instruments', batch_gather: 'gathering the archives',
    batch_archives: 'every archive of each star (DACE, CARMENES DR1, VizieR), found by its APERO name', batch_regather: 'gather them again',
    apero_refresh: 'Refresh APERO names', apero_names: 'APERO names', objects: 'objects', apero_none: 'No copy of APERO’s names yet: made at the first file (a few MB).', apero_fetching: 'fetching APERO’s names...',
    col_span: 'Span [d]', col_fip: 'FIP (P or alias)', col_fip_alone: 'FIP (P alone)', open: 'Open', from_batch: 'from the batch FIP', batch_waiting: 'waiting',
    vizier_pub: 'VizieR (published)',
    zero_by: 'Zero', zero_fit: 'fitted offsets', zero_median: 'median', fitted_trend: 'the fitted trend', fitted_model: 'the fitted model (trend and signals)',
    zero_fit_note: 'one solution for every instrument (the trend and the signal), each instrument about it (the median of its velocities minus it, zero), the residuals under it', resid: 'residual', resid_axis: 'residuals [m/s]', rv_zero_axis: 'RV - zero [m/s]', trend_only: 'its trend: zoom in for its signal',
    mstar: 'Mass of the star [M☉]', sptype: 'Spectral type, mass', from_word: 'from', typed: 'typed',
    no_mstar: 'm sin i: give the mass of the star (the Star card)', mass_transit: 'mass (it transits: sin i ≈ 1)',
    fipview_joint: 'the period alone, and with its aliases', fipview_each: 'each instrument too',
    transit_word: 'transit', transit_fold: 'a fold on this transit ephemeris: phase 0 at the transit, its uncertainty at the velocities',
    on_transit: 'on its transit ephemeris', cycles_from: 'periods from', phase_word: 'phase', phase_fixed: 'with its phase fixed by the transits',
    conj_alone: 'the velocities alone put the conjunction',
    fold_at: 'Fold at P =', fold_go: 'Fold', or_click: 'or click the periodogram', asked: 'a period asked', folding: 'folding...',
    fmodel_sine: 'sinusoid', fmodel_kepler: 'Keplerian (e free)', fmodel_kepler_short: 'Keplerian',
    fitting_kepler: 'fitting the Keplerian orbit (outliers, jitters, trend; a few seconds)...',
    fold_note_kep: 'nightly means; koloa’s fit, the eccentricity and period free within the peak; phase 0 at the conjunction',
    tick_tip: 'take this signal out: its sum with the other ticked ones on the series, and the FIP of what is left', tick_hint: 'tick a signal to take it out',
    tick_fitting: 'fitting each ticked signal without the others...', without: 'without', stage_running: 'the FIP of what is left...',
    band_note: 'its range: zoom in for its cycles', fold_minus: 'the other ticked signals taken out:',
    bic_vs_none: 'against no planet', bic_vs_circ: 'against a circular orbit', bic_sign: '(positive: this orbit is the better model)',
    bic_very: 'its eccentricity very strongly favoured', bic_strong: 'its eccentricity strongly favoured', bic_positive: 'its eccentricity favoured',
    bic_none: 'its eccentricity not called for', bic_against: 'a circular orbit favoured',
    subtract: 'Subtract this solution, FIP of the residuals', residuals_of: 'FIP of the residuals, without:', back_series: 'Back to the series',
    colour_by: 'Colour by', fcol_inst: 'instrument', fcol_date: 'date', overlay: 'the solution on the time series', offscale: 'beyond the range',
    accel_star: 'acceleration of the star', no_trend_fit: 'no trend fitted (the box “acceleration (trend)” of the report is unticked)',
    accel_note: 'the fit of the quick look, perspective acceleration included',
    tab_survey: 'Survey', sv_title: 'A survey', sv_from: 'Spectral type from', sv_to: 'to', sv_dmax: 'within [pc]', sv_dec: 'declination from',
    sv_hint: 'The stars that meet a few constraints, asked of SIMBAD; which of them have velocities in the archives or in your files; then the quick FIP of those ticked, here or on another machine.',
    sv_vmax: 'V brighter than', sv_dwarfs: 'dwarfs only', sv_ask: 'Ask SIMBAD', sv_check: 'Check the archives',
    sv_check_hint: 'of the stars ticked (of all when none is): DACE, CARMENES DR1, the surveys on VizieR',
    sv_folders: 'Your files (LBL...), in the folders:', sv_match: 'Match the files', sv_tick: 'Tick:', sv_tick_data: 'those with data',
    sv_tick_files: 'those with files', sv_tick_all: 'all', sv_tick_none: 'none', sv_col_name: 'Star', sv_col_type: 'Type',
    sv_col_surveys: 'Surveys (VizieR)', sv_col_files: 'files', sv_stars: 'stars', sv_with_data: 'with data', sv_ticked: 'ticked', sv_shown: 'shown',
    sv_checked: 'stars checked', sv_matched: 'files put with their stars', sv_unmatched: 'files of stars not in the sample',
    sv_near: 'within 5 arcsec of', sv_none_ticked: 'No star ticked.', sv_missing: 'nothing for',
    sv_batch: 'The batch of the stars ticked', sv_run: 'Run it here', sv_run_hint: 'each star with its files and every archive of it, in the Batch FIP tab',
    sv_name: 'Or pack it as', sv_out: 'made in', sv_root: 'where it will be there (its ROOT)', sv_jobs: 'stars at once',
    sv_gather: 'gather the archives first', sv_tess: 'with the light curves', sv_pack: 'Pack it', sv_results: 'The results of a batch folder',
    sv_open: 'Open them', sv_todo: 'stars not done yet',
    term_title: 'A terminal, to the machine that runs the batch', term_new: 'New terminal', term_go: 'Go there',
    term_remember: 'Remember how I got here...', term_forget_btn: 'forget', term_forget: 'Forget the route', term_tab: 'terminal', term_here: 'here',
    term_types: 'Typed for you, to read then run with Enter:', term_copy: 'copy the tar there', term_launch: 'unpack and launch',
    term_follow: 'follow its log', term_fetch: 'bring the results back', term_no_route: 'no route kept yet',
    term_none: 'No terminal on this system.', term_key: 'The terminal answers only to the page opened from the address koloa printed when it started (it ends with ?key=...).',
    term_no_lib: 'The terminal of the page could not be loaded (no network?).', term_need_route: 'A route with its host and its folder first (Remember how I got here...).',
    term_need_batch: 'Pack a batch first.', term_need_server: 'In a terminal that is on the server (Go there).',
    route_title: 'How to get there', route_hint: 'The lines typed in this terminal that it showed (a password is never shown: never kept). Correct them: one line each, as you would type them.',
    route_name: 'Name', route_host: 'Host, as ssh and rsync name it', route_folder: 'Folder of the batches there', route_save: 'Remember', cancel: 'Cancel',
    tab_analysis: 'Analysis', tab_remembered: 'Remembered targets', remember: 'Remember the result', note: 'Note',
    remembered_ok: 'remembered: see the Remembered targets tab', recalled_from: 'recalled: remembered on', recall: 'Recall', forget: 'Forget',
    forget_confirm: 'Forget this result and its copies?', rem_title: 'Remembered targets',
    rem_hint: 'The quick looks remembered, the last first: recall one to find the page as it was (its velocities, its FIP, its folds).',
    no_remembered: 'Nothing remembered yet: once a quick FIP has ended, Remember the result keeps it here.',
    col_target: 'Target', col_when: 'Remembered', col_nights: 'Nights', col_peaks: 'Numbered peaks: P [d] (FIP)',
    col_known: 'Known planets', col_each: 'Each instrument: its best peak', col_data: 'Data', remembering: 'remembering...',
    fold_title: 'The fold at a peak', fold_note: 'nightly means; a sinusoid with an offset per instrument and the trend; phase 0 at the conjunction',
    phase_axis: 'phase (0 = conjunction)',
    quick_title: 'Quick FIP (no GP): a look before the report', quick_noise: 'fitting the noise of each instrument',
    quick_fip1: 'the FIP, first pass', quick_planets: 'fitting the signals found and the known planets', quick_fip2: 'the FIP, second pass',
    quick_inst: 'the FIP of each instrument on its own', inst_alone: 'alone', few_nights: 'fewer than 10 nights: no FIP of its own',
    best_peaks: 'its best peaks: P [d] (FIP, P or alias)', none_found: 'none', each_title: 'Each instrument on its own',
    arch_plot: 'Archives gathered for the star:', arch_auto: 'no file: the archives gathered for the star are shown',
    nothing_gathered: 'nothing gathered yet', points_word: 'points',
    refip: 'Recompute the FIP', stale: 'made before the last change of what is shown: recompute it', alone: 'the period alone',
    family: 'the period or any of its aliases', strongest: 'strongest', no_gp: 'no GP', signals_word: 'signals', sweeps_word: 'sweeps',
    period_axis: 'period [d]', fip_axis: '-log10 FIP (1 % dotted)', nights: 'nights', passes: 'pass(es)', waiting: 'waiting for the previous one',
    need_archive: 'No file: tick DACE, CARMENES DR1 or VizieR (in the options below) for the velocities of the archives.',
    clip: 'twice the 3 to 97 percentile range', full_range: 'Full range', sliders: 'the sliders',
    use_sho: 'use for the SHO', use_sho_head: 'pick one for the SHO', no_sho: 'none (no SHO)',
    new_target: 'New target', known: 'Known planets (NASA Exoplanet Archive, exoplanet.eu)', none_known: 'none in the archive nor in exoplanet.eu',
    eu_only: 'exoplanet.eu only', eu_pending: 'exoplanet.eu: its catalogue is being fetched (a few minutes, once)',
    toi_title: 'TESS Objects of Interest', toi_pick: 'click a TOI to fit it with the ephemeris of TESS', toi_none: 'none',
    toi_on: 'TESS ephemerides (TOIs)', tois: 'TOIs (empty: all)', transit: 'transit',
    from_disk: 'from the disk', disk_archives: 'the archives folder', disk_kept: 'the copy kept', asked_now: 'asked the network now', refresh_star: 'Refresh',
    solution: 'solution', solutions: 'solutions',
    files_in: 'Files (optional): one per instrument or reduction', add_file: '+ Add a file', inst_auto: 'instrument (auto)', remove: 'remove',
    gather: 'Gather the archives',
    gather_hint: 'What DACE (with your key when there is one), CARMENES DR1, VizieR (the velocities published: Keck HIRES, the APF and the Lick Hamilton of the US surveys, HARPS by SERVAL, the star’s papers) and TESS (and Kepler, K2 and CoRoT, when the star lay in their fields) have of the star, kept in the archives folder, one folder per star.',
    gather_phot: 'TESS, Kepler, K2, CoRoT', ts_others_none: 'It lay in no field of Kepler, K2 or CoRoT either.',
    ts_mission_tip: 'the light curve searched: TESS’s, or that of Kepler, K2 or CoRoT when the star lay in their fields',
    vz_papers: 'the star’s papers',
    transit_title: 'A transit in TESS', ts_scan: 'the period scanned:', ts_periods: 'periods', ts_folded: 'folded at', ts_none: 'No peak below a FIP of 10 %, and no transiting planet or TOI: no transit to look for.',
    ts_fetching: 'looking for a transit in TESS (the light curve fetched from MAST when the archives do not have it: up to a minute)...',
    ts_searching: 'looking for a transit in TESS...', ts_missing: 'no TESS light curve of the star', ts_fetch: 'Fetch it from MAST',
    ts_unobserved: 'TESS has not observed this star yet: no sector at MAST holds its position (its sectors leave gaps), so there is no light curve to look for a transit in.',
    ts_planned: 'The mission’s pointing table (tess-point) holds it in sectors',
    ts_frames: 'no TESS light curve of this star at MAST: it is in the full frames only, of sectors',
    colon: ':', src_archives: 'the archives of the star', src_cds: 'the CDS', src_mass: 'its mass (R ~ M^0.9)', src_not_gathered: 'not gathered',
    ts_why_none: 'no box with data where the transit would be',
    ts_why_weak: 'the deepest box is {snr} sigma deep (a transit: {need} or more)',
    ts_why_single: '{snr} sigma, but in a single transit',
    ts_why_one: '{snr} sigma, but from one transit: {drop} sigma without it',
    ts_why_deep: '{snr} sigma, but too deep for a planet (an eclipsing binary?)',
    ts_why_dips: '{snr} sigma in {n} transits, but as strong a box at another period ({drop} sigma without its deepest transit, {top} at one of {trials}): this light curve makes dips of its own',
    ts_why_chance: '{snr} sigma in {n} transits, but within reach of chance (the null of {trials} other periods reaches {drop} sigma one time in {odds})',
    ts_why_plausible: '{snr} sigma in {n} transits ({drop} without the deepest; at most {top} at {trials} other periods, chance {chance})',
    ts_plausible: 'a plausible transit', ts_not: 'no plausible transit', ts_box: 'the box:', ts_phase: 'phase', dace_copy_off: 'DACE’s copy of the file’s spectra (its pipeline, less precise than LBL) left out by default: tick it to add it', ts_view_fold: 'fold', ts_view_series: 'time series, each transit', ts_bins_hour: 'medians in bins of an hour', ts_sector: 'sector', teq_sub: 'eq', ts_fit: 'the box fitted to the medians:', ts_fitname: 'the box fitted to the medians', ts_expected: 'a central transit',
    ts_window: 'searched within', ts_whole: 'searched over the whole phase', ts_points: 'each point', ts_bins: 'medians in bins',
    ts_boxname: 'the box found', ts_axis: 'hours from the expected transit', ts_axis_found: 'hours from the transit found',
    ts_flux: 'flux, high-passed [ppt]', col_transit: 'Transit (TESS)', batch_tess: 'a transit in TESS', ts_no_tess: 'no TESS', ts_not_observed: 'not observed by TESS yet', inst_head: 'instrument', src_spt: 'its spectral type (Pecaut & Mamajek 2013)', src_archive: 'the NASA Exoplanet Archive', from_file: 'From the file', file_name: 'file name', for_simbad: 'for SIMBAD', not_apero: 'not in APERO’s names: as it is',
    copy: 'Copy', copied: 'Copied', run_gather: 'Gather', detailed: 'Detailed report',
    detailed_hint: 'Everything koloa can say about the star, as a PDF report: from the file, from the SIMBAD name (the archives only), or (best) from both.',
    outdir: 'Report folder', kmax: 'Signals (kmax)', nsweep: 'Sweeps', nburn: 'Burn-in',
    pmin: 'Shortest period [d]', pmax: 'Longest period [d]', periods: 'More periods to test [d]',
    dmap: 'Detection map', dmap_none: 'none', dmap_fip: 'by the FIP (hours)', dmap_search: 'blind search (minutes)',
    fip_gp: 'GP in the FIP, by band', exposures: 'every exposure (not nightly means)',
    archive: 'Exoplanet Archive', gpcheck: 'signals against a GP', duck: 'duck test', latex: 'PDF report',
    run_detailed: 'Make the report', analysis_script: 'Analysis script', runs: 'Runs',
    runs_hint: 'Each run is the command line above, in a process of its own; closing the page does not stop it, stopping the server does.',
    no_runs: 'No run yet.', cwd: 'Runs from', rel: '(relative paths start there)',
    stop: 'Stop', log: 'Log', report: 'Report (PDF)', files: 'Files',
    running: 'running', done: 'done', failed: 'failed', stopped: 'stopped',
    resolving: 'Asking SIMBAD...', loading: 'Reading the velocities...',
    main: 'SIMBAD', tic: 'TIC', gaia: 'Gaia DR3', pos: 'RA, Dec (J2000)', names: 'Other names',
    rotation: 'Periods of variability', carmenes: 'CARMENES DR1', none: 'none',
    not_in: 'not in DR1', points: 'points', inst: 'Instrument', n: 'N', rms: 'rms [m/s]',
    no_rv: 'No velocity: give a file, or gather the archives of the star first.',
    starting: 'starting', need: 'Give a file, a SIMBAD name, or (best) both.',
    left_fip: 'left for this FIP', pass_left: 'the whole pass: at most about', of_up_to: 'of at most',
    trend: 'acceleration (trend)', curvature: 'its change (curvature)',
    fip_gp_sel: 'GP in the FIP', gp_banded: 'by period band', gp_sho: 'SHO at the rotation', gp_none: 'none',
    rotation_p: 'Rotation period [d]', refresh: 'Refresh', on_disk: 'already on disk', busy: 'being gathered...',
    pick: 'click a period to use it (an SHO at it in the FIP)',
    archive_kept: 'NASA Exoplanet Archive: kept, of', archive_none: 'NASA Exoplanet Archive: not kept yet (fetched at the first report)',
    refresh_archive: 'Refresh the archive',
  },
  fr: {
    tagline: 'Vitesses radiales robustes aux valeurs aberrantes · sur cette machine', docs: 'Docs',
    star: 'Étoile', star_hint: 'Son nom SIMBAD : les archives, DACE, CARMENES et TESS la trouvent par ce nom.',
    resolve: 'Résoudre', velocities: 'Vitesses',
    velocities_hint: 'Vos fichiers de vitesses (LBL .rdb, csv, csv de DACE), et les archives récupérées pour l’étoile quand DACE ou CARMENES DR1 est cochée dans le rapport : ce que le rapport utilisera. Décochez un instrument (ou cliquez-le dans la légende) pour l’écarter du rapport.',
    use: 'Utilisé', exclude: 'Instruments écartés', source: 'Source', src_file: 'fichier d’entrée',
    include: 'Jeux de données rappelés', ds_why: 'Par défaut', ds_back: 'rappelé', ds_sets: 'jeux',
    ds_release: 'écarté : les mêmes spectres que {better}{how}',
    ds_more: ', le plus précis (celui-ci disperse {extra} m/s de plus, en quadrature, sur les spectres communs)',
    ds_tie: ' (aussi précis, avec plus de spectres)', ds_err: ' (erreurs médianes de {a} contre {b} m/s)',
    ds_weak: 'écarté : ses {nights} nuits ne contraignent ni la moyenne ({mean} % de son erreur) ni la pente ({slope} %)',
    ds_weak_one: 'écarté : une seule nuit ne contraint ni la moyenne ni la pente', ds_asked: ', rappelé',
    ds_weak_part: 'écarté : {same} de ses points sont des spectres de {better}{how} ; les {left} autres ne contraignent ni la moyenne ni la pente',
    ds_part: '{used} points sur {n} : les autres sont des spectres de {better}',
    batch_rules: 'parmi elles, la meilleure publication des mêmes spectres, et pas les jeux qui ne contraignent rien',
    file: 'Fichier (facultatif)', root: 'Dossier des archives', plot: 'Tracer', browse: 'Parcourir...', picking: 'choix en cours...',
    var_pending: 'les périodes de variabilité de SIMBAD arrivent (son service est lent aujourd’hui)',
    p_valid: 'P(valide)', p_valid_low: 'P(valide) < 50 %',
    fip_prompt: 'Plusieurs instruments : décochez ceux à écarter (moins, plus vite), puis calculez le FIP rapide ; il peut prendre des minutes sur une longue série.',
    start_fip: 'Calculer le FIP rapide', stop_fip: 'Arrêter', fip_stopped: 'arrêté',
    envelope: 'enveloppe à 1σ', published: 'publiée', known_fold: 'un repliement à la période publiée de cette planète connue, avec son orbite publiée (en tirets)',
    tab_batch: 'FIP en lot', batch_pick: 'Ajouter des fichiers...', batch_or: 'ou les fichiers d’un dossier :', batch_add: 'Les ajouter', batch_clear: 'Vider',
    batch_run: 'Lancer les FIP', batch_none: 'Aucun fichier encore.', batch_files: 'fichiers', batch_nomatch: 'Aucun fichier de ce dossier ne correspond au motif.',
    batch_hint: 'Le FIP rapide de nombreux fichiers, l’un après l’autre : chaque fichier seul, ou avec toutes les archives de son étoile (trouvée par son nom APERO) ; un tableau de leurs objets, de leurs meilleurs pics et de leurs accélérations ; ouvrez-en un dans l’onglet Analyse pour regarder de plus près.',
    col_file: 'Fichier', col_object: 'Objet', col_insts: 'Nuits : instruments', batch_gather: 'récupération des archives',
    batch_archives: 'toutes les archives de chaque étoile (DACE, CARMENES DR1, VizieR), trouvée par son nom APERO', batch_regather: 'les récupérer à nouveau',
    apero_refresh: 'Rafraîchir les noms APERO', apero_names: 'Noms APERO', objects: 'objets', apero_none: 'Pas encore de copie des noms d’APERO : faite au premier fichier (quelques Mo).', apero_fetching: 'récupération des noms d’APERO...',
    col_span: 'Durée [j]', col_fip: 'FIP (P ou alias)', col_fip_alone: 'FIP (P seule)', open: 'Ouvrir', from_batch: 'du FIP en lot', batch_waiting: 'en attente',
    vizier_pub: 'VizieR (publiées)',
    zero_by: 'Zéro', zero_fit: 'offsets ajustés', zero_median: 'médiane', fitted_trend: 'la tendance ajustée', fitted_model: 'le modèle ajusté (tendance et signaux)',
    zero_fit_note: 'une solution pour tous les instruments (la tendance et le signal), chaque instrument autour d’elle (la médiane de ses vitesses moins elle, nulle), les résidus dessous', resid: 'résidu', resid_axis: 'résidus [m/s]', rv_zero_axis: 'RV - zéro [m/s]', trend_only: 'sa tendance : zoomez pour son signal',
    mstar: 'Masse de l’étoile [M☉]', sptype: 'Type spectral, masse', from_word: 'de', typed: 'tapée',
    no_mstar: 'm sin i : donnez la masse de l’étoile (la carte Étoile)', mass_transit: 'masse (elle transite : sin i ≈ 1)',
    fipview_joint: 'la période seule, et avec ses alias', fipview_each: 'chaque instrument aussi',
    transit_word: 'transit', transit_fold: 'un repliement sur cette éphéméride de transit : phase 0 au transit, son incertitude aux vitesses',
    on_transit: 'sur son éphéméride de transit', cycles_from: 'périodes depuis', phase_word: 'phase', phase_fixed: 'sa phase fixée par les transits',
    conj_alone: 'les vitesses seules placent la conjonction',
    fold_at: 'Replier à P =', fold_go: 'Replier', or_click: 'ou cliquez le périodogramme', asked: 'une période demandée', folding: 'repliement...',
    fmodel_sine: 'sinusoïde', fmodel_kepler: 'képlérienne (e libre)', fmodel_kepler_short: 'Képlérienne',
    fitting_kepler: 'ajustement de l’orbite képlérienne (valeurs aberrantes, gigues, tendance ; quelques secondes)...',
    fold_note_kep: 'moyennes par nuit ; l’ajustement de koloa, l’excentricité et la période libres dans le pic ; phase 0 à la conjonction',
    tick_tip: 'retirer ce signal : sa somme avec les autres signaux cochés sur la série, et le FIP de ce qui reste', tick_hint: 'cochez un signal pour le retirer',
    tick_fitting: 'ajustement de chaque signal coché sans les autres...', without: 'sans', stage_running: 'le FIP de ce qui reste...',
    band_note: 'son étendue : zoomez pour voir ses cycles', fold_minus: 'les autres signaux cochés retirés :',
    bic_vs_none: 'face à aucune planète', bic_vs_circ: 'face à une orbite circulaire', bic_sign: '(positif : cette orbite est le meilleur modèle)',
    bic_very: 'son excentricité très fortement favorisée', bic_strong: 'son excentricité fortement favorisée', bic_positive: 'son excentricité favorisée',
    bic_none: 'son excentricité n’est pas requise', bic_against: 'une orbite circulaire est favorisée',
    subtract: 'Soustraire cette solution, FIP des résidus', residuals_of: 'FIP des résidus, sans :', back_series: 'Revenir à la série',
    colour_by: 'Couleur selon', fcol_inst: 'instrument', fcol_date: 'date', overlay: 'la solution sur la série temporelle', offscale: 'hors de la plage',
    accel_star: 'accélération de l’étoile', no_trend_fit: 'pas de tendance ajustée (la case « accélération (tendance) » du rapport est décochée)',
    accel_note: 'l’ajustement du coup d’œil, accélération de perspective incluse',
    tab_survey: 'Relevé', sv_title: 'Un relevé', sv_from: 'Type spectral de', sv_to: 'à', sv_dmax: 'à moins de [pc]', sv_dec: 'déclinaison de',
    sv_hint: 'Les étoiles qui remplissent quelques contraintes, demandées à SIMBAD ; lesquelles ont des vitesses dans les archives ou dans vos fichiers ; puis le FIP rapide de celles cochées, ici ou sur une autre machine.',
    sv_vmax: 'V plus brillant que', sv_dwarfs: 'naines seulement', sv_ask: 'Demander à SIMBAD', sv_check: 'Vérifier les archives',
    sv_check_hint: 'des étoiles cochées (de toutes si aucune ne l’est) : DACE, CARMENES DR1, les relevés sur VizieR',
    sv_folders: 'Vos fichiers (LBL...), dans les dossiers :', sv_match: 'Apparier les fichiers', sv_tick: 'Cocher :', sv_tick_data: 'celles avec données',
    sv_tick_files: 'celles avec fichiers', sv_tick_all: 'toutes', sv_tick_none: 'aucune', sv_col_name: 'Étoile', sv_col_type: 'Type',
    sv_col_surveys: 'Relevés (VizieR)', sv_col_files: 'fichiers', sv_stars: 'étoiles', sv_with_data: 'avec données', sv_ticked: 'cochées', sv_shown: 'montrées',
    sv_checked: 'étoiles vérifiées', sv_matched: 'fichiers appariés à leur étoile', sv_unmatched: 'fichiers d’étoiles hors de l’échantillon',
    sv_near: 'à moins de 5 arcsec de', sv_none_ticked: 'Aucune étoile cochée.', sv_missing: 'rien pour',
    sv_batch: 'Le lot des étoiles cochées', sv_run: 'Le lancer ici', sv_run_hint: 'chaque étoile avec ses fichiers et toutes ses archives, dans l’onglet FIP en lot',
    sv_name: 'Ou l’empaqueter sous le nom', sv_out: 'fait dans', sv_root: 'où il sera là-bas (son ROOT)', sv_jobs: 'étoiles à la fois',
    sv_gather: 'récupérer d’abord les archives', sv_tess: 'avec les courbes de lumière', sv_pack: 'Empaqueter', sv_results: 'Les résultats d’un dossier de lot',
    sv_open: 'Les ouvrir', sv_todo: 'étoiles pas encore faites',
    term_title: 'Un terminal, vers la machine qui exécute le lot', term_new: 'Nouveau terminal', term_go: 'Y aller',
    term_remember: 'Retenir comment je suis arrivé ici...', term_forget_btn: 'oublier', term_forget: 'Oublier la route', term_tab: 'terminal', term_here: 'ici',
    term_types: 'Tapé pour vous, à lire puis à lancer avec Entrée :', term_copy: 'copier le tar là-bas', term_launch: 'dépaqueter et lancer',
    term_follow: 'suivre son journal', term_fetch: 'rapporter les résultats', term_no_route: 'aucune route retenue',
    term_none: 'Pas de terminal sur ce système.', term_key: 'Le terminal ne répond qu’à la page ouverte depuis l’adresse que koloa a affichée au démarrage (elle finit par ?key=...).',
    term_no_lib: 'Le terminal de la page n’a pas pu être chargé (pas de réseau ?).', term_need_route: 'D’abord une route avec son hôte et son dossier (Retenir comment je suis arrivé ici...).',
    term_need_batch: 'Empaquetez d’abord un lot.', term_need_server: 'Dans un terminal qui est sur le serveur (Y aller).',
    route_title: 'Comment y aller', route_hint: 'Les lignes tapées dans ce terminal qu’il a montrées (un mot de passe n’est jamais montré : jamais retenu). Corrigez-les : une ligne chacune, comme vous les taperiez.',
    route_name: 'Nom', route_host: 'Hôte, comme ssh et rsync le nomment', route_folder: 'Dossier des lots là-bas', route_save: 'Retenir', cancel: 'Annuler',
    tab_analysis: 'Analyse', tab_remembered: 'Cibles retenues', remember: 'Retenir le résultat', note: 'Note',
    remembered_ok: 'retenu : voyez l’onglet Cibles retenues', recalled_from: 'rappelé : retenu le', recall: 'Rappeler', forget: 'Oublier',
    forget_confirm: 'Oublier ce résultat et ses copies ?', rem_title: 'Cibles retenues',
    rem_hint: 'Les coups d’œil retenus, le dernier en premier : rappelez-en un pour retrouver la page telle qu’elle était (ses vitesses, son FIP, ses repliements).',
    no_remembered: 'Rien de retenu encore : une fois un FIP rapide terminé, Retenir le résultat le garde ici.',
    col_target: 'Cible', col_when: 'Retenu le', col_nights: 'Nuits', col_peaks: 'Pics numérotés : P [j] (FIP)',
    col_known: 'Planètes connues', col_each: 'Chaque instrument : son meilleur pic', col_data: 'Données', remembering: 'en cours...',
    fold_title: 'Le repliement à un pic', fold_note: 'moyennes par nuit ; une sinusoïde avec un offset par instrument et la tendance ; phase 0 à la conjonction',
    phase_axis: 'phase (0 = conjonction)',
    quick_title: 'FIP rapide (sans GP) : un coup d’œil avant le rapport', quick_noise: 'ajustement du bruit de chaque instrument',
    quick_fip1: 'le FIP, premier passage', quick_planets: 'ajustement des signaux trouvés et des planètes connues', quick_fip2: 'le FIP, second passage',
    quick_inst: 'le FIP de chaque instrument seul', inst_alone: 'seul', few_nights: 'moins de 10 nuits : pas de FIP à lui',
    best_peaks: 'ses meilleurs pics : P [j] (FIP, P ou alias)', none_found: 'aucun', each_title: 'Chaque instrument seul',
    arch_plot: 'Archives récupérées pour l’étoile :', arch_auto: 'pas de fichier : les archives récupérées pour l’étoile sont montrées',
    nothing_gathered: 'rien de récupéré encore', points_word: 'points',
    refip: 'Refaire le FIP', stale: 'fait avant le dernier changement de ce qui est montré : refaites-le', alone: 'la période seule',
    family: 'la période ou un de ses alias', strongest: 'le plus fort', no_gp: 'sans GP', signals_word: 'signaux', sweeps_word: 'itérations',
    period_axis: 'période [j]', fip_axis: '-log10 FIP (1 % en pointillé)', nights: 'nuits', passes: 'passage(s)', waiting: 'en attente du précédent',
    need_archive: 'Pas de fichier : cochez DACE, CARMENES DR1 ou VizieR (dans les options ci-dessous) pour les vitesses des archives.',
    clip: 'deux fois l’écart des centiles 3 à 97', full_range: 'Tout voir', sliders: 'les curseurs',
    use_sho: 'utiliser pour le SHO', use_sho_head: 'choisissez-en une pour le SHO', no_sho: 'aucune (pas de SHO)',
    new_target: 'Nouvelle cible', known: 'Planètes connues (NASA Exoplanet Archive, exoplanet.eu)', none_known: 'aucune dans l’archive ni dans exoplanet.eu',
    eu_only: 'exoplanet.eu seulement', eu_pending: 'exoplanet.eu : son catalogue est en cours de récupération (quelques minutes, une fois)',
    toi_title: 'TESS Objects of Interest', toi_pick: 'cliquez un TOI pour l’ajuster avec l’éphéméride de TESS', toi_none: 'aucun',
    toi_on: 'éphémérides TESS (TOI)', tois: 'TOI (vide : tous)', transit: 'transit',
    from_disk: 'lu sur le disque', disk_archives: 'le dossier des archives', disk_kept: 'la copie gardée', asked_now: 'demandé au réseau à l’instant', refresh_star: 'Rafraîchir',
    solution: 'solution', solutions: 'solutions',
    files_in: 'Fichiers (facultatifs) : un par instrument ou réduction', add_file: '+ Ajouter un fichier', inst_auto: 'instrument (auto)', remove: 'retirer',
    gather: 'Récupérer les archives',
    gather_hint: 'Ce que DACE (avec votre clé s’il y en a une), CARMENES DR1, VizieR (les vitesses publiées : Keck HIRES, l’APF et le Lick Hamilton des relevés américains, HARPS par SERVAL, les articles de l’étoile) et TESS (et Kepler, K2 et CoRoT, si l’étoile était dans leurs champs) ont de l’étoile, rangé dans le dossier des archives, un dossier par étoile.',
    gather_phot: 'TESS, Kepler, K2, CoRoT', ts_others_none: 'Elle n’était dans aucun champ de Kepler, K2 ou CoRoT non plus.',
    ts_mission_tip: 'la courbe de lumière cherchée : celle de TESS, ou celle de Kepler, K2 ou CoRoT si l’étoile était dans leurs champs',
    vz_papers: 'les articles de l’étoile',
    transit_title: 'Un transit dans TESS', ts_scan: 'la période balayée :', ts_periods: 'périodes', ts_folded: 'repliée à', ts_none: 'Aucun pic sous un FIP de 10 %, et aucune planète qui transite ni TOI : pas de transit à chercher.',
    ts_fetching: 'recherche d’un transit dans TESS (la courbe de lumière récupérée à MAST quand les archives ne l’ont pas : jusqu’à une minute)...',
    ts_searching: 'recherche d’un transit dans TESS...', ts_missing: 'pas de courbe de lumière TESS de l’étoile', ts_fetch: 'La récupérer à MAST',
    ts_unobserved: 'TESS n’a pas encore observé cette étoile : aucun secteur à MAST ne contient sa position (ses secteurs laissent des trous), il n’y a donc pas de courbe de lumière où chercher un transit.',
    ts_planned: 'Le plan de pointage de la mission (tess-point) la place dans les secteurs',
    ts_frames: 'pas de courbe de lumière TESS de cette étoile à MAST : elle n’est que dans les images plein champ, des secteurs',
    colon: ' :', src_archives: 'les archives de l’étoile', src_cds: 'le CDS', src_mass: 'sa masse (R ~ M^0.9)', src_not_gathered: 'pas récupérée',
    ts_why_none: 'aucune boîte avec des données là où serait le transit',
    ts_why_weak: 'la boîte la plus profonde est à {snr} sigma (un transit : {need} ou plus)',
    ts_why_single: '{snr} sigma, mais en un seul transit',
    ts_why_one: '{snr} sigma, mais d’un seul transit : {drop} sigma sans lui',
    ts_why_deep: '{snr} sigma, mais trop profond pour une planète (une binaire à éclipses ?)',
    ts_why_dips: '{snr} sigma en {n} transits, mais une boîte aussi forte à une autre période ({drop} sigma sans son transit le plus profond, {top} à l’une de {trials}) : cette courbe de lumière fait ses propres creux',
    ts_why_chance: '{snr} sigma en {n} transits, mais à la portée du hasard (le nul de {trials} autres périodes atteint {drop} sigma une fois sur {odds})',
    ts_why_plausible: '{snr} sigma en {n} transits ({drop} sans le plus profond ; au plus {top} à {trials} autres périodes, chance {chance})',
    ts_plausible: 'un transit plausible', ts_not: 'pas de transit plausible', ts_box: 'la boîte :', ts_phase: 'phase', dace_copy_off: 'la copie DACE des spectres du fichier (son pipeline, moins précis que LBL) écartée par défaut : cochez-la pour l’ajouter', ts_view_fold: 'repliement', ts_view_series: 'série temporelle, chaque transit', ts_bins_hour: 'médianes par heure', ts_sector: 'secteur', teq_sub: 'éq', ts_fit: 'la boîte ajustée aux médianes :', ts_fitname: 'la boîte ajustée aux médianes', ts_expected: 'un transit central',
    ts_window: 'cherché à', ts_whole: 'cherché sur toute la phase', ts_points: 'chaque point', ts_bins: 'médianes par intervalle',
    ts_boxname: 'la boîte trouvée', ts_axis: 'heures depuis le transit attendu', ts_axis_found: 'heures depuis le transit trouvé',
    ts_flux: 'flux, passe-haut [ppt]', col_transit: 'Transit (TESS)', batch_tess: 'un transit dans TESS', ts_no_tess: 'pas de TESS', ts_not_observed: 'pas encore observée par TESS', inst_head: 'instrument', src_spt: 'son type spectral (Pecaut & Mamajek 2013)', src_archive: 'la NASA Exoplanet Archive', from_file: 'Du fichier', file_name: 'nom de fichier', for_simbad: 'pour SIMBAD', not_apero: 'pas dans les noms d’APERO : tel quel',
    copy: 'Copier', copied: 'Copié', run_gather: 'Récupérer', detailed: 'Rapport détaillé',
    detailed_hint: 'Tout ce que koloa peut dire de l’étoile, en un rapport PDF : à partir du fichier, du nom SIMBAD (les archives seules), ou (le mieux) des deux.',
    outdir: 'Dossier du rapport', kmax: 'Signaux (kmax)', nsweep: 'Itérations', nburn: 'Rodage',
    pmin: 'Période minimale [j]', pmax: 'Période maximale [j]', periods: 'Autres périodes à tester [j]',
    dmap: 'Carte de détection', dmap_none: 'aucune', dmap_fip: 'par le FIP (heures)', dmap_search: 'recherche aveugle (minutes)',
    fip_gp: 'GP dans le FIP, par bande', exposures: 'chaque pose (pas les moyennes par nuit)',
    archive: 'Exoplanet Archive', gpcheck: 'signaux face à un GP', duck: 'duck test', latex: 'rapport PDF',
    run_detailed: 'Faire le rapport', analysis_script: 'Script d’analyse', runs: 'Exécutions',
    runs_hint: 'Chaque exécution est la ligne de commande ci-dessus, dans son propre processus ; fermer la page ne l’arrête pas, arrêter le serveur oui.',
    no_runs: 'Aucune exécution pour l’instant.', cwd: 'Lancé depuis', rel: '(les chemins relatifs partent de là)',
    stop: 'Arrêter', log: 'Journal', report: 'Rapport (PDF)', files: 'Fichiers',
    running: 'en cours', done: 'terminé', failed: 'échec', stopped: 'arrêté',
    resolving: 'Question à SIMBAD...', loading: 'Lecture des vitesses...',
    main: 'SIMBAD', tic: 'TIC', gaia: 'Gaia DR3', pos: 'AD, Déc (J2000)', names: 'Autres noms',
    rotation: 'Périodes de variabilité', carmenes: 'CARMENES DR1', none: 'aucune',
    not_in: 'pas dans DR1', points: 'points', inst: 'Instrument', n: 'N', rms: 'rms [m/s]',
    no_rv: 'Aucune vitesse : donnez un fichier, ou récupérez d’abord les archives de l’étoile.',
    starting: 'démarrage', need: 'Donnez un fichier, un nom SIMBAD, ou (le mieux) les deux.',
    left_fip: 'restantes pour ce FIP', pass_left: 'toute la passe : au plus environ', of_up_to: 'sur au plus',
    trend: 'accélération (tendance)', curvature: 'sa variation (courbure)',
    fip_gp_sel: 'GP dans le FIP', gp_banded: 'par bande de période', gp_sho: 'SHO à la rotation', gp_none: 'aucun',
    rotation_p: 'Période de rotation [j]', refresh: 'Rafraîchir', on_disk: 'déjà sur le disque', busy: 'récupération en cours...',
    pick: 'cliquez une période pour l’utiliser (un SHO à cette période dans le FIP)',
    archive_kept: 'NASA Exoplanet Archive : copie locale du', archive_none: 'NASA Exoplanet Archive : pas encore de copie locale (faite au premier rapport)',
    refresh_archive: 'Rafraîchir l’archive',
  },
};
// the instruments: eight hues checked for colour-blind separation on the
// dark surface, each with its own marker
const COLOURS = ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767'];
const SYMBOLS = ['circle', 'square', 'diamond', 'triangle-up', 'triangle-down', 'cross', 'star', 'hexagon'];

let lang = 'en';
try { lang = localStorage.getItem('koloa-lang') || ((navigator.language || '').startsWith('fr') ? 'fr' : 'en'); } catch (e) { /* no storage */ }
const t = (key) => (TEXT[lang] && TEXT[lang][key]) || TEXT.en[key] || key;
// a sentence with its {names} filled
const fill = (text, vals) => String(text).replace(/\{(\w+)\}/g, (all, name) => (name in vals ? vals[name] : all));
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
// a FIP in powers of ten (HTML): one too small for a float, below 1e-300
const fipExp = (v) => (+v > 0 ? (+v).toExponential(1) : '&lt; 1e-300');


// each page says what to redraw when the language changes
const langHooks = [];
function applyLang() {
  document.documentElement.lang = lang;
  document.querySelectorAll('[data-i18n]').forEach((el) => { el.textContent = t(el.dataset.i18n); });
  if ($('lang')) $('lang').textContent = lang === 'fr' ? 'EN' : 'FR';
  langHooks.forEach((hook) => hook());
}

async function api(path, body) {
  const resp = await fetch(path, body === undefined ? {} : {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  });
  const data = await resp.json();
  if (!resp.ok || data.error) throw new Error(data.error || resp.statusText);
  return data;
}

// -----------------------------------------------------------------------------
// the runs
// -----------------------------------------------------------------------------
const jobs = new Map(); // id -> state, with all its lines
// the logs closed by hand (a log is open by default)
const closedLogs = new Set();

function clock(sec) {
  sec = Math.max(0, Math.round(sec));
  const h = Math.floor(sec / 3600), m = Math.floor((sec % 3600) / 60), s = sec % 60;
  return h ? `${h} h ${String(m).padStart(2, '0')} min` : m ? `${m} min ${String(s).padStart(2, '0')} s` : `${s} s`;
}

// the sweeps of the FIP that runs: a bar, the time left for it, and for the
// whole pass when it is one FIP of several (the descent of the banded FIP)
function progressHtml(p) {
  const frac = Math.min(1, p.done / Math.max(p.total, 1));
  const pct = Math.round(100 * frac);
  const perFip = frac > 0.02 ? p.seconds / frac : null;
  const m = p.label.match(/FIP (\d+) of up to (\d+) \((.*)\)$/);
  let name = p.label;
  let pass = '';
  if (m) {
    const k = +m[1], n = +m[2];
    name = `FIP ${k} ${t('of_up_to')} ${n} (${m[3]})`;
    if (perFip && n > k) pass = ` \u00b7 ${t('pass_left')} ${clock((1 - frac) * perFip + (n - k) * perFip)}`;
  }
  const left = perFip ? ` \u00b7 ~${clock((1 - frac) * perFip)} ${t('left_fip')}` : '';
  return `<span class="pbar"><span style="width:${pct}%"></span></span>`
    + `<span class="ptext">${esc(name)}: ${pct} %${left}${pass}</span>`;
}

function renderJobs() {
  const box = $('jobs');
  if (!jobs.size) { box.innerHTML = `<p class="hint">${esc(t('no_runs'))}</p>`; return; }
  box.innerHTML = [...jobs.values()].reverse().map((job) => {
    const running = job.status === 'running';
    const steps = job.steps.map((st, i) => {
      const last = i === job.steps.length - 1;
      const icon = last && running ? '<span class="hourglass">⏳</span>'
        : last && job.status !== 'done' ? '<span class="bad">✗</span>' : '<span class="ok">✓</span>';
      const detail = (last && running && st.progress ? progressHtml(st.progress) : '')
        + (last && running && st.detail && !(st.progress && st.detail.includes('% of the sweeps'))
          ? `<span class="detail">${esc(st.detail)}</span>` : '');
      const result = st.result ? `<span class="result">${esc(st.result)}</span>` : '';
      return `<li><span class="icon">${icon}</span><span>${esc(st.name)}</span><span class="time">${clock(st.elapsed)}</span>${result}${detail}</li>`;
    }).join('') || (running ? `<li><span class="icon"><span class="hourglass">⏳</span></span><span>${esc(t('starting'))}</span><span></span></li>` : '');
    const base = `/api/output?id=${job.id}&name=`;
    const files = (job.files || []).map((f) => `<a href="${base}${encodeURIComponent(f)}" target="_blank">${esc(f)}</a>`).join('');
    const report = job.report ? `<a href="${base}${encodeURIComponent(job.report)}" target="_blank"><button type="button">${esc(t('report'))}</button></a>` : '';
    const what = { gather: t('gather'), archive: t('refresh_archive') }[job.action] || t('detailed');
    return `<div class="job" id="job-${job.id}">
      <div class="job-head"><span class="what">${esc(what)}</span>
        <span class="status ${job.status}">${esc(t(job.status))}</span>
        <span class="when">${clock(job.elapsed)}</span><span class="spacer"></span>
        ${report}${running ? `<button class="stop" type="button" data-stop="${job.id}">${esc(t('stop'))}</button>` : ''}</div>
      <pre class="codeblock">${esc(job.line)}</pre>
      <ol class="steps">${steps}</ol>
      ${files ? `<div class="files"><span class="hint">${esc(t('files'))}:</span>${files}</div>` : ''}
      <details data-log="${job.id}"${closedLogs.has(job.id) ? '' : ' open'}><summary>${esc(t('log'))} (${job.lines.length})</summary>
        <pre class="codeblock log">${esc(job.lines.slice(-400).join('\n'))}</pre></details>
    </div>`;
  }).join('');
  box.querySelectorAll('details[data-log]').forEach((d) => {
    d.addEventListener('toggle', () => { if (d.open) closedLogs.delete(d.dataset.log); else closedLogs.add(d.dataset.log); });
    const pre = d.querySelector('pre');
    pre.scrollTop = pre.scrollHeight;
  });
}

function keep(state) {
  const old = jobs.get(state.id);
  const lines = old ? old.lines.concat(state.lines) : state.lines;
  jobs.set(state.id, { ...state, lines });
}

async function poll() {
  for (const job of jobs.values()) {
    if (job.status !== 'running') continue;
    try {
      keep(await api(`/api/job?id=${job.id}&since=${job.lines.length}`));
      if (jobs.get(job.id).status !== 'running') jobDoneHooks.forEach((hook) => hook(jobs.get(job.id)));
    } catch (err) { /* the server may be gone */ }
  }
  renderJobs();
}

// what a page does when one of its runs ends, or starts
const jobDoneHooks = [];
const jobStartHooks = [];
async function run(action, opts, path) {
  try {
    keep(await api(path || '/api/run', { action, options: opts }));
    renderJobs();
    jobStartHooks.forEach((hook) => hook(action));
    $('jobs').scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (err) {
    alert(err.message);
  }
}


// -----------------------------------------------------------------------------
// the explanations: an (i) next to a field shows its own on hover or focus
// -----------------------------------------------------------------------------
const HELP = {
  target: {
    en: 'The name of the star as SIMBAD knows it: GJ 436, Ross 905, HD 69830, TOI-700 or TIC 307210830 all work. koloa resolves it (CDS Sesame) into its main identifier, its Gaia DR3 and TIC numbers and its aliases, and every archive (DACE, CARMENES DR1, the NASA Exoplanet Archive, the TOIs, TESS) is then searched with those. Its known planets come from the NASA Exoplanet Archive and from exoplanet.eu: the archive’s values, exoplanet.eu’s beside them (they mostly agree, not always), and the planets or times of transit only one of them has. Empty when a file is given: koloa fills it from the file, its OBJECT column (LBL keeps the header key of each spectrum), else its name (lbl_<OBJECT>_<TEMPLATE>.rdb), found in APERO’s database of names (the uniform names of SPIRou and NIRPS: AN_SEX, GL699, PROXIMA) and given as its SIMBAD name; the line under it says how. A name that database does not have is used as it is, which SIMBAD may not know. Once asked, the answer is kept on disk.',
    fr: 'Le nom de l’étoile tel que SIMBAD le connaît : GJ 436, Ross 905, HD 69830, TOI-700 ou TIC 307210830 fonctionnent tous. koloa le résout (CDS Sesame) en son identifiant principal, ses numéros Gaia DR3 et TIC et ses alias, et chaque archive (DACE, CARMENES DR1, la NASA Exoplanet Archive, les TOI, TESS) est ensuite cherchée avec eux. Ses planètes connues viennent de la NASA Exoplanet Archive et d’exoplanet.eu : les valeurs de l’archive, celles d’exoplanet.eu à côté (elles concordent le plus souvent, pas toujours), et les planètes ou temps de transit que seule l’une des deux a. Vide quand un fichier est donné : koloa le remplit d’après le fichier, sa colonne OBJECT (LBL garde la clé d’en-tête de chaque spectre), sinon son nom (lbl_<OBJET>_<GABARIT>.rdb), trouvé dans la base de noms d’APERO (les noms uniformes de SPIRou et NIRPS : AN_SEX, GL699, PROXIMA) et donné par son nom SIMBAD ; la ligne dessous dit comment. Un nom que cette base n’a pas est utilisé tel quel, et SIMBAD ne le connaît pas toujours. Une fois demandée, la réponse est gardée sur le disque.',
  },
  files: {
    en: 'Your own velocities, one row per file: an LBL .rdb (SPIRou, NIRPS, HARPS, ESPRESSO, HARPS-N reduced with LBL), a csv with rjd, vrad and svrad columns (and indicators with their errors, as name and sname), or a csv written by DACE. Several files go into one analysis, each instrument with its own offset and jitter. All empty: the report comes from the archives alone (DACE, CARMENES DR1, VizieR). A spectrum both in a file and on DACE is counted once. Browse... opens a dialog of this machine.',
    fr: 'Vos propres vitesses, une ligne par fichier : un .rdb de LBL (SPIRou, NIRPS, HARPS, ESPRESSO, HARPS-N réduits avec LBL), un csv avec les colonnes rjd, vrad et svrad (et des indicateurs avec leurs erreurs, nom et snom), ou un csv écrit par DACE. Plusieurs fichiers vont dans une même analyse, chaque instrument avec son offset et son jitter. Tout vide : le rapport vient des archives seules (DACE, CARMENES DR1, VizieR). Un spectre à la fois dans un fichier et sur DACE est compté une fois. Parcourir... ouvre une boîte de dialogue de cette machine.',
  },
  file_label: {
    en: 'The name of this file’s instrument in the analysis. Empty: koloa reads it from the file’s own columns (the names of its reduced files, r.HARPS..., r.ESPRE..., r.NIRPS..., r.HARPN..., or APERO’s keys for SPIRou) and splits HARPS and ESPRESSO in their eras as DACE does: HARPS03 and HARPS15 at the 2015 fibre change, ESPRESSO18 and ESPRESSO19 at the 2019 intervention, each with its own zero point. Give a name (NIRPS_LBL2, SPIRou_APERO07) to set two reductions of one instrument apart, or to correct a wrong guess.',
    fr: 'Le nom de l’instrument de ce fichier dans l’analyse. Vide : koloa le lit dans les colonnes du fichier (les noms de ses fichiers réduits, r.HARPS..., r.ESPRE..., r.NIRPS..., r.HARPN..., ou les clés d’APERO pour SPIRou) et sépare HARPS et ESPRESSO en leurs époques comme DACE : HARPS03 et HARPS15 au changement de fibres de 2015, ESPRESSO18 et ESPRESSO19 à l’intervention de 2019, chacune avec son point zéro. Donnez un nom (NIRPS_LBL2, SPIRou_APERO07) pour distinguer deux réductions d’un même instrument, ou pour corriger une mauvaise devinette.',
  },
  root: {
    en: 'The folder where the archives gathered for each star are kept, one folder per star (archives/GJ_436/): the DACE and CARMENES velocities, the TESS light curves and SIMBAD’s answer. A relative path starts from the folder koloa runs from (shown under the report). The plot reads what is there, and the resolver looks there before asking the network. Default: archives.',
    fr: 'Le dossier où sont gardées les archives récupérées pour chaque étoile, un dossier par étoile (archives/GJ_436/) : les vitesses de DACE et de CARMENES, les courbes de lumière TESS et la réponse de SIMBAD. Un chemin relatif part du dossier d’où koloa est lancé (indiqué sous le rapport). Le graphique lit ce qui s’y trouve, et le résolveur y regarde avant de demander au réseau. Par défaut : archives.',
  },
  gather_sources: {
    en: 'The archives to gather. DACE: HARPS, ESPRESSO, NIRPS, CORALIE... (with your key, what your account may see). CARMENES DR1: about 360 M dwarfs of the north, 2016 to 2020. TESS: the light curve of every sector at MAST.\n\nVizieR: the velocities published for the star, each source its own box: Keck HIRES (Teklu et al. 2025; Tal-Or et al. 2019), the California Legacy Survey (HIRES, APF, Lick), the Lick Hamilton (Fischer et al. 2014), HARPS by SERVAL (Trifonov et al. 2020), and the tables of the star’s papers. A spectrum two sources have is kept once. A minute or a few for a well-studied star.',
    fr: 'Les archives à récupérer. DACE : HARPS, ESPRESSO, NIRPS, CORALIE... (avec votre clé, ce que votre compte peut voir). CARMENES DR1 : environ 360 naines M du nord, 2016 à 2020. TESS : la courbe de lumière de chaque secteur à MAST.\n\nVizieR : les vitesses publiées pour l’étoile, chaque source sa case : Keck HIRES (Teklu et al. 2025 ; Tal-Or et al. 2019), le California Legacy Survey (HIRES, APF, Lick), le Lick Hamilton (Fischer et al. 2014), HARPS par SERVAL (Trifonov et al. 2020), et les tables des articles de l’étoile. Un spectre que deux sources ont est gardé une fois. Une minute ou quelques-unes pour une étoile très étudiée.',
  },
  analysis_script: {
    en: 'A kit to run the analysis yourself, and change it: a .tar.gz with analysis.py, a Python script that calls koloa’s routines step by step (the velocities read and put together, a fit of their noise, the FIP in two passes, the Keplerian orbits of the signals with their errors and minimum masses, the acceleration of the star, the folds, the FIP of the residuals), every step explained in its comments (in the language of this page), its figures saved as PDF; the files of velocities given; the archives gathered for the star (DACE, CARMENES DR1, VizieR, the TESS light curve); and star.yaml, what SIMBAD, the NASA Exoplanet Archive and TESS say of the star (star.json, the same, without PyYAML). Its settings are those of this card (the archives ticked, the instruments left out, the trend, the FIP), and it runs offline: python analysis.py in its folder. It may hold data that are not public (your files, DACE with a key).',
    fr: 'Un paquet pour faire l’analyse vous-même, et la modifier : un .tar.gz avec analysis.py, un script Python qui appelle les routines de koloa étape par étape (les vitesses lues et mises ensemble, un ajustement de leur bruit, le FIP en deux passages, les orbites képlériennes des signaux avec leurs erreurs et leurs masses minimales, l’accélération de l’étoile, les repliements, le FIP des résidus), chaque étape expliquée dans ses commentaires (dans la langue de cette page), ses figures sauvegardées en PDF ; les fichiers de vitesses donnés ; les archives récupérées pour l’étoile (DACE, CARMENES DR1, VizieR, la courbe de lumière TESS) ; et star.yaml, ce que SIMBAD, la NASA Exoplanet Archive et TESS disent de l’étoile (star.json, le même, sans PyYAML). Ses réglages sont ceux de cette carte (les archives cochées, les instruments écartés, la tendance, le FIP), et il tourne sans réseau : python analysis.py dans son dossier. Il peut contenir des données non publiques (vos fichiers, DACE avec une clé).',
  },
  outdir: {
    en: 'Where the report goes: the PDF, the text report, the summary (JSON), every figure as its own PDF, and what was fetched for it (DACE, CARMENES, VizieR). Empty: reports/<star>, as the command line shows. A relative path starts from the folder koloa runs from.',
    fr: 'Où va le rapport : le PDF, le rapport texte, le résumé (JSON), chaque figure en PDF, et ce qui a été récupéré pour lui (DACE, CARMENES, VizieR). Vide : reports/<étoile>, comme le montre la ligne de commande. Un chemin relatif part du dossier d’où koloa est lancé.',
  },
  kmax: {
    en: 'The largest number of signals the FIP looks for at once. It samples models with 0 to kmax sinusoids together with the noise, so that a strong signal does not hide a weaker one. Default 3. More costs more time per sweep (about in proportion), and rarely changes what is found unless the star has many planets.',
    fr: 'Le plus grand nombre de signaux que le FIP cherche à la fois. Il échantillonne des modèles de 0 à kmax sinusoïdes avec le bruit, pour qu’un signal fort n’en cache pas un plus faible. Par défaut 3. Plus coûte plus de temps par itération (à peu près en proportion), et change rarement ce qui est trouvé, sauf pour une étoile à nombreuses planètes.',
  },
  nsweep: {
    en: 'The sweeps kept in each of the two Markov chains of the FIP. The FIP of a period is the fraction of sweeps with no signal there, so the smallest it can tell is about 1 / (2 x sweeps): 2000, the default, tells FIPs down to a few 1e-4. A sweep takes about half a second for 200 nights.',
    fr: 'Les itérations gardées dans chacune des deux chaînes de Markov du FIP. Le FIP d’une période est la fraction des itérations sans signal à cet endroit, donc le plus petit qu’il puisse dire est environ 1 / (2 x itérations) : 2000, la valeur par défaut, distingue des FIP jusqu’à quelques 1e-4. Une itération prend environ une demi-seconde pour 200 nuits.',
  },
  nburn: {
    en: 'The first sweeps of each chain, thrown away while it settles (the jitters, the outliers and the signals finding their place). Default 400. Raise it when the log shows the number of signals still drifting at the end of the burn-in.',
    fr: 'Les premières itérations de chaque chaîne, jetées le temps qu’elle se stabilise (les jitters, les valeurs aberrantes et les signaux trouvent leur place). Par défaut 400. Augmentez-le si le journal montre le nombre de signaux qui dérive encore à la fin du rodage.',
  },
  pmin: {
    en: 'The shortest period searched [days]. Default 1.1, clear of the 1-day alias. Lower it for ultra-short periods (0.3 for instance): the grid of frequencies grows as 1 / pmin, and so does the time.',
    fr: 'La période la plus courte cherchée [jours]. Par défaut 1,1, à l’écart de l’alias d’un jour. Baissez-la pour les périodes ultra-courtes (0,3 par exemple) : la grille de fréquences croît comme 1 / pmin, et le temps aussi.',
  },
  pmax: {
    en: 'The longest period searched [days], 5000 for instance. Empty: twice the baseline of the data. Beyond the baseline, an orbit is mostly a curvature, which the trend and its change take (acceleration, curvature).',
    fr: 'La période la plus longue cherchée [jours], 5000 par exemple. Vide : deux fois la durée couverte par les données. Au-delà, une orbite est surtout une courbure, que prennent la tendance et sa variation (accélération, courbure).',
  },
  periods: {
    en: 'Periods to test whatever the FIP says [days], separated by spaces: 113.46 27.1 for instance, a candidate from a paper the archive does not list, or a period from photometry. Each is fitted with the other signals, free within 2 %, and goes through the duck test. Empty: none; the known planets of the archive are tested anyway.',
    fr: 'Des périodes à tester quoi qu’en dise le FIP [jours], séparées par des espaces : 113.46 27.1 par exemple, un candidat d’un article que l’archive ne liste pas, ou une période de la photométrie. Chacune est ajustée avec les autres signaux, libre à 2 % près, et passe le duck test. Vide : aucune ; les planètes connues de l’archive sont testées de toute façon.',
  },
  fip_gp: {
    en: 'The model of the stellar activity inside the FIP, fitted together with the signals (its hyperparameters sampled with the jitters, one GP per instrument). By period band (default): from the longest periods down, a local GP that cannot reach the band it decides, made more flexible only when the data ask for it. SHO at the rotation: an SHO at a rotation period you trust and one at its half (celerite’s RotationTerm), over every period at once, with no bands; it needs the rotation period below. None: no GP, the activity left to the jitters and to the duck test.',
    fr: 'Le modèle de l’activité stellaire dans le FIP, ajusté avec les signaux (ses hyperparamètres échantillonnés avec les jitters, un GP par instrument). Par bande de période (par défaut) : des périodes les plus longues vers les plus courtes, un GP local qui ne peut pas atteindre la bande qu’il décide, rendu plus souple seulement si les données le demandent. SHO à la rotation : un SHO à une période de rotation fiable et un à sa moitié (le RotationTerm de celerite), sur toutes les périodes à la fois, sans bandes ; il demande la période de rotation ci-dessous. Aucun : pas de GP, l’activité laissée aux jitters et au duck test.',
  },
  rotation: {
    en: 'A rotation period you trust [days], a published one: 2.704 for Wolf 359 for instance. It becomes the star’s rotation for every check (a signal at P_rot, P_rot/2 or 2 P_rot is flagged) and the prior of the GP check, and with SHO in the GP menu, the period of the GP inside the FIP (held). Tick a period of the literature in the resolver (use for the SHO) to fill it in and choose the SHO. Empty: the archive’s rotation, when it has one.',
    fr: 'Une période de rotation fiable [jours], publiée : 2,704 pour Wolf 359 par exemple. Elle devient la rotation de l’étoile pour chaque vérification (un signal à P_rot, P_rot/2 ou 2 P_rot est signalé) et l’a priori du test GP, et avec SHO dans le menu GP, la période du GP dans le FIP (fixée). Cochez une période de la littérature dans le résolveur (utiliser pour le SHO) pour la remplir et choisir le SHO. Vide : la rotation de l’archive, si elle en a une.',
  },
  toi_on: {
    en: 'Fit the TESS Objects of Interest of the star with the ephemerides of TESS: the period and the time of a transit, with their errors from the TOI list, as gaussian priors; K free and positive, the phase held by the transit. Each TOI is fitted whether the FIP finds it or not, before any period found near it, and its fold has phase 0 at the transit.',
    fr: 'Ajuster les TESS Objects of Interest de l’étoile avec les éphémérides de TESS : la période et le temps d’un transit, avec leurs erreurs tirées de la liste des TOI, comme a priori gaussiens ; K libre et positif, la phase tenue par le transit. Chaque TOI est ajusté que le FIP le trouve ou non, avant toute période trouvée près de lui, et son repliement a la phase 0 au transit.',
  },
  tois: {
    en: 'The TOIs to fit, by number, separated by spaces: 175.01 175.02 for instance. Empty with the box ticked: every TOI of the star but the false positives (FP, FA). Click a TOI in the resolver to add it here.',
    fr: 'Les TOI à ajuster, par numéro, séparés par des espaces : 175.01 175.02 par exemple. Vide avec la case cochée : tous les TOI de l’étoile sauf les faux positifs (FP, FA). Cliquez un TOI dans le résolveur pour l’ajouter ici.',
  },
  exclude: {
    en: 'Instruments left out of the report, by their names in the plot, separated by commas: NIRPS_DACE, HARPS03, HIRES (CLS) for instance. Those the rules leave out by default are put here when the velocities are plotted (the table says why). Untick an instrument in the table under the plot, or click it in the legend, to add it here. They are left out once the files, DACE, CARMENES and VizieR are put together. Empty: none.',
    fr: 'Les instruments écartés du rapport, par leurs noms dans le graphique, séparés par des virgules : NIRPS_DACE, HARPS03, HIRES (CLS) par exemple. Ceux que les règles écartent par défaut y sont mis quand les vitesses sont tracées (le tableau dit pourquoi). Décochez un instrument dans le tableau sous le graphique, ou cliquez-le dans la légende, pour l’ajouter ici. Ils sont écartés une fois les fichiers, DACE, CARMENES et VizieR rassemblés. Vide : aucun.',
  },
  survey: {
    en: 'A survey: every star of SIMBAD within a distance (from its parallax) and in a range of spectral types (M0 to M9: every M; M3 to M5.5; K5 to M2), dwarfs unless unticked, between two declinations, brighter than a V if given. A system and its components are objects of their own in SIMBAD (\u29c9: another object of the sample within 5 arcsec).\n\nCheck the archives: which of the stars have velocities, without gathering them: CARMENES DR1 and the surveys on VizieR (Keck HIRES, the APF, Lick, HARPS by SERVAL) by the star\u2019s position in their lists, and DACE, asked star by star, six at a time (a few seconds a star; what it answers is kept in the archives folder, where the gathering reads it back). The papers of a star and its light curves come when its archives are gathered.\n\nYour files: folders of LBL files (.rdb) of SPIRou, of NIRPS...; the star of each file is its OBJECT column, else its name, through APERO\u2019s names, and it is put with the star of the sample that has that name. The files of stars that are not in the sample are counted (their names under the cursor).\n\nTick the stars to run: those with data (the others have nothing to run), those with files, all, none; click a column to sort by it.',
    fr: 'Un relevé : toutes les étoiles de SIMBAD à moins d\u2019une distance (d\u2019après la parallaxe) et dans un intervalle de types spectraux (M0 à M9 : toutes les M ; M3 à M5.5 ; K5 à M2), des naines sauf si décoché, entre deux déclinaisons, plus brillantes qu\u2019un V s\u2019il est donné. Un système et ses composantes sont des objets distincts dans SIMBAD (\u29c9 : un autre objet de l\u2019échantillon à moins de 5 arcsec).\n\nVérifier les archives : lesquelles de ces étoiles ont des vitesses, sans les récupérer : CARMENES DR1 et les relevés sur VizieR (Keck HIRES, l\u2019APF, Lick, HARPS par SERVAL) par la position de l\u2019étoile dans leurs listes, et DACE, demandé étoile par étoile, six à la fois (quelques secondes par étoile ; sa réponse est gardée dans le dossier des archives, où la collecte la relit). Les articles d\u2019une étoile et ses courbes de lumière viennent quand ses archives sont récupérées.\n\nVos fichiers : des dossiers de fichiers LBL (.rdb) de SPIRou, de NIRPS... ; l\u2019étoile de chaque fichier est sa colonne OBJECT, sinon son nom, par les noms d\u2019APERO, et il est mis avec l\u2019étoile de l\u2019échantillon qui porte ce nom. Les fichiers d\u2019étoiles hors de l\u2019échantillon sont comptés (leurs noms sous le curseur).\n\nCochez les étoiles à lancer : celles avec données (les autres n\u2019ont rien à lancer), celles avec fichiers, toutes, aucune ; cliquez une colonne pour trier.',
  },
  survey_batch: {
    en: 'The batch of the stars ticked: for each, its files and every archive of it (gathered when not on disk), its datasets chosen by the rules, its quick FIP, a transit looked for at its best peak.\n\nRun it here: in the Batch FIP tab, one star after the other.\n\nOr pack it for another machine (a server, a cluster): one folder, and its .tar.gz, with everything the batch needs and no path of this machine in it: the files of each star, its archives (gathered here first, so that nothing is asked of the network there), what koloa fetched once (the NASA Exoplanet Archive, the lists of the surveys), koloa itself, and run_batch.py. The first setting of that script is ROOT, where the folder is on the machine that runs it: given here (the folder of the route chosen below, and the name of the batch), to be changed in the script when the folder goes elsewhere. There: tar xzf, then python run_batch.py (several stars at once; submit.sh does the same as a SLURM job array). A batch stopped goes on where it was.\n\nThe results of a batch folder (its results/, brought back) open as a batch of the page: its table, each star opened in the Analysis tab with its FIP as it was computed.',
    fr: 'Le lot des étoiles cochées : pour chacune, ses fichiers et toutes ses archives (récupérées si elles ne sont pas sur le disque), ses jeux de données choisis par les règles, son FIP rapide, un transit cherché à son meilleur pic.\n\nLe lancer ici : dans l\u2019onglet FIP en lot, une étoile après l\u2019autre.\n\nOu l\u2019empaqueter pour une autre machine (un serveur, une grappe) : un dossier, et son .tar.gz, avec tout ce dont le lot a besoin et aucun chemin de cette machine : les fichiers de chaque étoile, ses archives (récupérées ici d\u2019abord, pour que rien ne soit demandé au réseau là-bas), ce que koloa a téléchargé une fois (la NASA Exoplanet Archive, les listes des relevés), koloa lui-même, et run_batch.py. Le premier réglage de ce script est ROOT, où se trouve le dossier sur la machine qui l\u2019exécute : donné ici (le dossier de la route choisie plus bas, et le nom du lot), à changer dans le script si le dossier va ailleurs. Là-bas : tar xzf, puis python run_batch.py (plusieurs étoiles à la fois ; submit.sh fait la même chose en tableau de tâches SLURM). Un lot arrêté reprend où il en était.\n\nLes résultats d\u2019un dossier de lot (son results/, rapporté) s\u2019ouvrent comme un lot de la page : son tableau, chaque étoile ouverte dans l\u2019onglet Analyse avec son FIP tel qu\u2019il a été calculé.',
  },
  terminal: {
    en: 'A terminal in the page: a shell of this machine, to go to the machine that will run the batch (ssh; its password or second factor is typed as in any terminal), and several of them (tabs).\n\nRemember how I got here: the lines you typed in this terminal that it showed (ssh maestria, cd /data/me/batches, what loads Python) are proposed, to be corrected, with the host as ssh and rsync name it and the folder of the batches there. A password is typed with the echo off: it is never shown, so never kept. Go there types the lines of a route again, each once the terminal is quiet (type a password when it is asked: the lines go on after it).\n\nThe four buttons type a command and do not run it: read it, change it if need be, press Enter. Copy the tar there: rsync from this machine (in a terminal of here) to the folder of the route. Unpack and launch: in the terminal that is on the server, tar xzf, then run_batch.py under nohup (it goes on when you leave). Follow its log: tail -f. Bring the results back: rsync of its results/ to the batch folder here, which Open them then reads.\n\nA shell is more than this page is otherwise trusted with: the terminal answers only to the page opened from the address koloa printed when it started (it ends with ?key=...).',
    fr: 'Un terminal dans la page : un shell de cette machine, pour aller sur la machine qui exécutera le lot (ssh ; son mot de passe ou son second facteur se tape comme dans tout terminal), et plusieurs à la fois (des onglets).\n\nRetenir comment je suis arrivé ici : les lignes que vous avez tapées dans ce terminal et qu\u2019il a montrées (ssh maestria, cd /data/moi/lots, ce qui charge Python) sont proposées, à corriger, avec l\u2019hôte tel que ssh et rsync le nomment et le dossier des lots là-bas. Un mot de passe se tape sans écho : il n\u2019est jamais montré, donc jamais retenu. Y aller retape les lignes d\u2019une route, chacune une fois le terminal au calme (tapez un mot de passe quand il est demandé : les lignes continuent après).\n\nLes quatre boutons tapent une commande sans la lancer : lisez-la, changez-la au besoin, appuyez sur Entrée. Copier le tar là-bas : rsync depuis cette machine (dans un terminal d\u2019ici) vers le dossier de la route. Dépaqueter et lancer : dans le terminal qui est sur le serveur, tar xzf, puis run_batch.py sous nohup (il continue quand vous partez). Suivre son journal : tail -f. Rapporter les résultats : rsync de son results/ vers le dossier du lot ici, que Les ouvrir lit ensuite.\n\nUn shell, c\u2019est plus que ce que cette page a autrement le droit de faire : le terminal ne répond qu\u2019à la page ouverte depuis l\u2019adresse que koloa a affichée au démarrage (elle finit par ?key=...).',
  },
  include: {
    en: 'Datasets used though the rules would leave them out (koloa.datasets: another release of the same spectra is more precise, or the dataset constrains nothing), by their names in the plot, separated by commas. Tick such a dataset in the table under the plot to add it here: it is then preferred to the other releases of its spectra. Empty: the rules as they are.',
    fr: 'Les jeux de données utilisés bien que les règles les écartent (koloa.datasets : une autre publication des mêmes spectres est plus précise, ou le jeu ne contraint rien), par leurs noms dans le graphique, séparés par des virgules. Cochez un tel jeu dans le tableau sous le graphique pour l’ajouter ici : il est alors préféré aux autres publications de ses spectres. Vide : les règles telles quelles.',
  },
  datasets: {
    en: 'Which datasets are used. The archives often hold the same spectra more than once (HARPS on DACE and by SERVAL, HIRES in three surveys), and old velocities too imprecise to constrain anything. Everything is plotted and can be ticked; by default:\n\n1. A spectrum several datasets have is taken from the most precise of them. On the spectra two releases share the star does the same in both, so their difference is noise alone, and it tells the noise of each: no model of the star is needed. When the two cannot be told apart: the one that names its spectrograph, the latest, the one with the more spectra. A release left with nothing of its own is unticked; one that has spectra the better one lacks keeps those (40/200 in the table).\n\n2. A dataset that constrains nothing is unticked. A line is fitted to the nightly means (an offset per dataset, one slope): a dataset is left out when taking it away makes the error of the mean and that of the slope grow by less than 1 %. The error of a night holds what the star does about a line, the same for every instrument, so a few nights among hundreds weigh little whatever their precision.\n\nTick a dataset to ask it back: a release asked back is preferred to the others for the spectra they share (untick the other to see it whole). The datasets of your own files are never left out. The report makes the same choice (--exclude, --include; --no-rules for no rule), and so does the batch.',
    fr: 'Quels jeux de données sont utilisés. Les archives ont souvent les mêmes spectres plus d’une fois (HARPS sur DACE et par SERVAL, HIRES dans trois relevés), et de vieilles vitesses trop peu précises pour contraindre quoi que ce soit. Tout est tracé et peut être coché ; par défaut :\n\n1. Un spectre que plusieurs jeux ont est pris du plus précis. Sur les spectres que deux publications partagent, l’étoile fait la même chose dans les deux : leur différence n’est que du bruit, et elle donne le bruit de chacune, sans modèle de l’étoile. Quand on ne peut pas les départager : celle qui nomme son spectrographe, la plus récente, celle qui a le plus de spectres. Une publication à qui il ne reste rien en propre est décochée ; celle qui a des spectres que la meilleure n’a pas garde ceux-là (40/200 dans le tableau).\n\n2. Un jeu qui ne contraint rien est décoché. Une droite est ajustée aux moyennes par nuit (un offset par jeu, une pente) : un jeu est écarté quand le retirer fait croître l’erreur de la moyenne et celle de la pente de moins de 1 %. L’erreur d’une nuit contient ce que fait l’étoile autour d’une droite, pareil pour chaque instrument : quelques nuits parmi des centaines pèsent peu, quelle que soit leur précision.\n\nCochez un jeu pour le rappeler : une publication rappelée est préférée aux autres pour les spectres communs (décochez l’autre pour la voir entière). Les jeux de vos propres fichiers ne sont jamais écartés. Le rapport fait le même choix (--exclude, --include ; --no-rules pour aucune règle), et le lot aussi.',
  },
  detection_map: {
    en: 'Which planets the series could have found, as a map of period and K. By the FIP: planets injected into the data (the signals found taken out) and looked for by the same FIP, with the same GP, as the decision itself; adaptive rounds by period band give the K found at 50 % and 90 %; hours on a long series. Blind search: a periodogram search instead, in minutes; it misses a planet that lands on an alias and has no GP, so it is pessimistic. Default: none.',
    fr: 'Les planètes que la série aurait pu trouver, en carte de période et de K. Par le FIP : des planètes injectées dans les données (les signaux trouvés retirés) et cherchées par le même FIP, avec le même GP, que la décision elle-même ; des tours adaptatifs par bande de période donnent le K trouvé à 50 % et 90 % ; des heures sur une longue série. Recherche aveugle : une recherche par périodogramme, en minutes ; elle rate une planète qui tombe sur un alias et n’a pas de GP, elle est donc pessimiste. Par défaut : aucune.',
  },
  exposures: {
    en: 'Analyse every exposure instead of the nightly means. By default the exposures of one night and one instrument are averaged into one point: an outlier is then a whole night, and the noise within a night weighs nothing. Ticked, the visits are kept, and an exposure or a whole visit can be an outlier.',
    fr: 'Analyser chaque pose au lieu des moyennes par nuit. Par défaut, les poses d’une nuit et d’un instrument sont moyennées en un point : une valeur aberrante est alors une nuit entière, et le bruit dans une nuit ne pèse pas. Cochée, les visites sont gardées, et une pose ou une visite entière peut être aberrante.',
  },
  dace: {
    en: 'Add the velocities DACE has of the star, from every instrument it holds. Off by default: the report uses your files only, and fetches nothing from DACE unless this is ticked. An instrument your files have too is set apart as <inst>_DACE (another pipeline, its own offset), and the spectra your files already have are not counted twice. Your DACE key, when there is one, adds what your account may see. Ticked, the plot shows the DACE velocities gathered for the star.',
    fr: 'Ajouter les vitesses que DACE a de l’étoile, de tous les instruments qu’il a. Désactivée par défaut : le rapport n’utilise que vos fichiers, et ne demande rien à DACE si elle n’est pas cochée. Un instrument que vos fichiers ont aussi est mis à part en <inst>_DACE (un autre pipeline, son propre offset), et les spectres que vos fichiers ont déjà ne sont pas comptés deux fois. Votre clé DACE, s’il y en a une, ajoute ce que votre compte peut voir. Cochée, le graphique montre les vitesses DACE récupérées pour l’étoile.',
  },
  carmenes: {
    en: 'Add the velocities of CARMENES DR1 (Ribas et al. 2023: the GTO of 2016 to 2020, about 360 M dwarfs of the north): SERVAL’s velocities corrected for the nightly zero points, with their activity indicators. Off by default: nothing is fetched unless this is ticked. Nothing is added for a star that is not in DR1. Ticked, the plot shows them.',
    fr: 'Ajouter les vitesses de CARMENES DR1 (Ribas et al. 2023 : le GTO de 2016 à 2020, environ 360 naines M du nord) : les vitesses de SERVAL corrigées des points zéro nocturnes, avec leurs indicateurs d’activité. Désactivée par défaut : rien n’est demandé si elle n’est pas cochée. Rien n’est ajouté pour une étoile absente de DR1. Cochée, le graphique les montre.',
  },
  tess_opt: {
    en: 'Fetch the TESS light curves of the star (MAST) and look in each sector for a photometric peak at the period of every signal, at its half, third or double: a spot rotating with the star moves both its light and its velocities. Beyond about 13 days, a period is out of reach of a 27-day sector.',
    fr: 'Récupérer les courbes de lumière TESS de l’étoile (MAST) et chercher dans chaque secteur un pic photométrique à la période de chaque signal, à sa moitié, son tiers ou son double : une tache qui tourne avec l’étoile change à la fois sa lumière et ses vitesses. Au-delà d’environ 13 jours, une période est hors de portée d’un secteur de 27 jours.',
  },
  vizier: {
    en: 'Add the velocities published on VizieR. The surveys, the star found by its position: Keck HIRES (Teklu et al. 2025; Tal-Or et al. 2019), the California Legacy Survey (HIRES, APF, Lick; Rosenthal et al. 2021), the Lick Hamilton (Fischer et al. 2014), HARPS by SERVAL (Trifonov et al. 2020). Then the tables of the star’s papers (those that are the star’s).\n\nEach source keeps its own offset; an instrument above 20 m/s of median error is left out; a spectrum published twice, or that your files or the archives have, is kept once. Off by default.',
    fr: 'Ajouter les vitesses publiées sur VizieR. Les relevés, l’étoile trouvée par sa position : Keck HIRES (Teklu et al. 2025 ; Tal-Or et al. 2019), le California Legacy Survey (HIRES, APF, Lick ; Rosenthal et al. 2021), le Lick Hamilton (Fischer et al. 2014), HARPS par SERVAL (Trifonov et al. 2020). Puis les tables des articles de l’étoile (celles qui sont les siennes).\n\nChaque source garde son propre offset ; un instrument à plus de 20 m/s d’erreur médiane est écarté ; un spectre publié deux fois, ou que vos fichiers ou les archives ont, est gardé une fois. Désactivé par défaut.',
  },

  archive_opt: {
    en: 'Use the NASA Exoplanet Archive (its copy kept on this machine) for the known planets of the star: each is tested at its period, and every signal is set against their published solutions (K, ephemeris). It also gives the star’s mass, for m sin i, and its rotation when known.',
    fr: 'Utiliser la NASA Exoplanet Archive (sa copie gardée sur cette machine) pour les planètes connues de l’étoile : chacune est testée à sa période, et chaque signal est comparé à leurs solutions publiées (K, éphéméride). Elle donne aussi la masse de l’étoile, pour m sin i, et sa rotation si elle est connue.',
  },
  gpcheck: {
    en: 'Fit the signals again with a GP of the activity (quasi-periodic, the rotation as its prior): how much likelihood each signal adds over the GP alone, periodograms whitened by the GP with false-alarm levels from simulations of its noise, and the whole series with the GP drawn over it. A planet holds against the GP; activity melts into it.',
    fr: 'Réajuster les signaux avec un GP de l’activité (quasi périodique, la rotation comme a priori) : la vraisemblance que chaque signal ajoute au GP seul, des périodogrammes blanchis par le GP avec des niveaux de fausse alarme tirés de simulations de son bruit, et toute la série avec le GP superposé. Une planète tient face au GP ; l’activité s’y fond.',
  },
  duck: {
    en: 'The duck test of every signal: does it walk like a planet? Its amplitude and phase stay put through the campaign (by halves, and visit by visit), no activity indicator peaks at its period, its eccentricity is not at the edge of its prior, a GP does not swallow it, TESS does not see it, and which of its aliases is the true period. Each check comes with its figures.',
    fr: 'Le duck test de chaque signal : marche-t-il comme une planète ? Son amplitude et sa phase restent stables au long de la campagne (par moitiés, et visite par visite), aucun indicateur d’activité n’a de pic à sa période, son excentricité n’est pas au bord de son a priori, un GP ne l’avale pas, TESS ne le voit pas, et lequel de ses alias est la vraie période. Chaque vérification vient avec ses figures.',
  },
  trend: {
    en: 'Fit a trend in time with the planets, in the likelihood: the acceleration of the star, dv/dt, one for all instruments beside their offsets (a distant companion, or the perspective acceleration of a nearby star when the pipeline leaves it in). It is reported in m/s/yr with its errors, and drawn with the model. On by default.',
    fr: 'Ajuster une tendance en temps avec les planètes, dans la vraisemblance : l’accélération de l’étoile, dv/dt, une pour tous les instruments à côté de leurs offsets (un compagnon lointain, ou l’accélération de perspective d’une étoile proche quand le pipeline la laisse). Elle est donnée en m/s/an avec ses erreurs, et tracée avec le modèle. Activée par défaut.',
  },
  curvature: {
    en: 'Fit the change of the acceleration too, d2v/dt2 [m/s/yr^2]: a second-order trend, for a companion whose orbit curves over the baseline. It takes some of the power of the longest periods. Off by default.',
    fr: 'Ajuster aussi la variation de l’accélération, d2v/dt2 [m/s/an^2] : une tendance du second ordre, pour un compagnon dont l’orbite se courbe sur la durée couverte. Elle prend une part de la puissance des plus longues périodes. Désactivée par défaut.',
  },
  mcmc: {
    en: 'Sample the orbits by MCMC (emcee) instead of the maximum a posteriori with Laplace errors: the errors of P, K, e and of the acceleration then come from the posterior itself, honest when it is not gaussian (a weak signal, a K near zero). Several minutes more per report.',
    fr: 'Échantillonner les orbites par MCMC (emcee) au lieu du maximum a posteriori avec les erreurs de Laplace : les erreurs de P, K, e et de l’accélération viennent alors de la loi a posteriori elle-même, justes quand elle n’est pas gaussienne (un signal faible, un K proche de zéro). Plusieurs minutes de plus par rapport.',
  },
  clip: {
    en: 'Show the velocity axis over twice the spread between the 3rd and the 97th percentiles of the velocities kept, about its middle: the bulk of the points fills the plot, and a few huge outliers no longer squash it. Points beyond are off the plot, not out of the analysis: a red triangle at the top or bottom edge points to each (hover: its velocity). Unticked: every point kept in view. Unticking an instrument recomputes it.',
    fr: 'Montrer l’axe des vitesses sur deux fois l’écart entre les centiles 3 et 97 des vitesses gardées, autour de leur milieu : le gros des points remplit le graphique, et quelques énormes valeurs aberrantes ne l’écrasent plus. Les points au-delà sont hors du graphique, pas hors de l’analyse : un triangle rouge au bord du haut ou du bas pointe vers chacun (au survol : sa vitesse). Décochée : tous les points gardés sont visibles. Décocher un instrument la recalcule.',
  },
  sliders: {
    en: 'The sliders along the axes cut the ranges shown: the one under the plot for the time, the one beside its left edge for the velocities (drag either end). The velocity slider is log-smart (asinh about the scatter of the points): it moves finely around the bulk and in big steps toward the outliers, so the 3 m/s scatter and a point at 1000 m/s both stay within reach. Zooming with the mouse moves them too; a double-click on the plot, or Full range, shows everything again.',
    fr: 'Les curseurs le long des axes coupent les intervalles montrés : celui sous le graphique pour le temps, celui le long de son bord gauche pour les vitesses (tirez l’un ou l’autre bout). Le curseur des vitesses est logarithmique (asinh autour de la dispersion des points) : il bouge finement autour du gros des points et à grands pas vers les valeurs aberrantes, si bien que la dispersion de 3 m/s et un point à 1000 m/s restent tous deux à portée. Un zoom à la souris les déplace aussi ; un double-clic sur le graphique, ou Tout voir, montre tout à nouveau.',
  },
  quickfip: {
    en: 'A first look before the report: the FIP of exactly what the plot shows (your files, the archives ticked, the instruments kept), outlier-aware, in nightly means, with no GP. It starts on its own for one instrument; with several, untick those to leave out, then Compute the quick FIP (Stop ends it).\n\nTwo passes, as in the report: the errors inflated to the noise of a fit without planets, a FIP; then to that of a fit with the signals found and the known planets, the FIP again. About a minute.\n\nWhite: the period or any of its aliases (what decides); grey: the period alone; dotted: FIP = 1 % and the window; dashed: the known planets. With several instruments, each has its own FIP in its colour. Without a GP, a peak at the rotation is expected.\n\nThe trend follows the report’s boxes; under the FIP, the acceleration of the star (dv/dt, m/s/yr) with its error, as the velocities measure it (the perspective acceleration included).',
    fr: 'Un premier coup d’œil avant le rapport : le FIP de ce que montre le graphique (vos fichiers, les archives cochées, les instruments gardés), robuste aux valeurs aberrantes, en moyennes par nuit, sans GP. Il démarre seul pour un instrument ; avec plusieurs, décochez ceux à écarter, puis Calculer le FIP rapide (Arrêter le termine).\n\nDeux passages, comme dans le rapport : les erreurs gonflées au bruit d’un ajustement sans planète, un FIP ; puis à celui d’un ajustement avec les signaux trouvés et les planètes connues, le FIP de nouveau. Environ une minute.\n\nBlanc : la période ou un de ses alias (ce qui décide) ; gris : la période seule ; pointillé : FIP = 1 % et la fenêtre ; tirets : les planètes connues. Avec plusieurs instruments, chacun a son FIP dans sa couleur. Sans GP, un pic à la rotation est attendu.\n\nLa tendance suit les cases du rapport ; sous le FIP, l’accélération de l’étoile (dv/dt, m/s/an) avec son erreur, telle que la mesurent les vitesses (l’accélération de perspective comprise).',
  },
  batch: {
    en: 'The quick FIP of many files: add them (the dialog, or a folder and a pattern, or several patterns separated by spaces: GL*nightly.rdb GJ*nightly.rdb) and run. One file after the other, the same quick FIP as the Analysis tab; each on its own, or with every archive of its star when ticked.\n\nThe object of each file is its APERO name (its OBJECT column, else its name), with its SIMBAD name under the cursor, its known planets (★, their periods) and its candidate TOIs (◐); the best period is marked ★ or ◐ when it falls on one (within a peak’s width).\n\nThe table: the object, the nights (and each instrument’s, with the archives), the span, the best peak and its FIP (green below 1 %), the K and rms of a sinusoid there, the acceleration of the star; and, for a peak below 10 %, a transit sought in TESS (one found makes a weak signal strong): ⚑ when plausible, else its S/N (the reason under the cursor).\n\nClick a header to sort, CSV to save it, Open to take a file to the Analysis tab. Stop ends the batch.',
    fr: 'Le FIP rapide de nombreux fichiers : ajoutez-les (la boîte de dialogue, ou un dossier et un motif, ou plusieurs motifs séparés par des espaces : GL*nightly.rdb GJ*nightly.rdb) et lancez. Un fichier après l’autre, le même FIP rapide que l’onglet Analyse ; chacun seul, ou avec toutes les archives de son étoile si c’est coché.\n\nL’objet de chaque fichier est son nom APERO (sa colonne OBJECT, sinon son nom), avec son nom SIMBAD sous le curseur, ses planètes connues (★, leurs périodes) et ses TOI candidates (◐) ; la meilleure période est marquée ★ ou ◐ quand elle tombe sur l’une d’elles (à la largeur d’un pic près).\n\nLe tableau : l’objet, les nuits (et celles de chaque instrument, avec les archives), la durée, le meilleur pic et son FIP (en vert sous 1 %), le K et le rms d’une sinusoïde là, l’accélération de l’étoile ; et, pour un pic sous 10 %, un transit cherché dans TESS (en trouver un rend fort un signal faible) : ⚑ s’il est plausible, sinon son S/B (la raison sous le curseur).\n\nCliquez un en-tête pour trier, CSV pour l’enregistrer, Ouvrir pour passer un fichier à l’onglet Analyse. Arrêter termine le lot.',
  },
  ts_view: {
    en: 'The fold, or the high-passed light curve of each sector on its own: each putative transit shaded, the fitted box drawn on it, and its own depth under it (orange when it is half the fitted one or more). A transit seen in every epoch, at about the same depth, is a consistent one; one deep event and nothing at the others, a glitch.',
    fr: 'Le repliement, ou la courbe de lumière filtrée de chaque secteur à part : chaque transit présumé ombré, la boîte ajustée dessinée dessus, et sa propre profondeur dessous (en orange quand elle fait au moins la moitié de celle ajustée). Un transit vu à chaque époque, à peu près à la même profondeur, est cohérent ; un seul événement profond et rien aux autres, un artefact.',
  },
  transit_search: {
    en: 'A transit in the TESS light curve at the period of a signal: each peak with a FIP below 10 % (a transit found would make a weak signal strong; about the conjunction of its fold, where a transiting planet transits), and each transiting planet or TOI of the star (at its ephemeris). The light curve comes from the star’s archives, else from MAST; a star TESS has not observed (its sectors leave gaps) has none, and the panel says so. Kepler, K2 and CoRoT are searched too when the star lay in their fields: the first mission with a plausible transit is shown, a button for each.\n\nThe search: the light curve high-passed (a running median three times the expected transit) and folded; a box sought within 3 sigma of the expected time carried to TESS, or over the whole phase. When the period’s error would smear the transits across the sectors, the period is scanned within 3 sigma and the light curve folded at the best.\n\nPlausible: 7 sigma or more, in 2 transits or more, still 3.5 without the deepest (one event is not a transit), and beyond the best box at periods where there is nothing (chance below 1 %). The other known planets are masked first. Tested on LHS 1140, TOI-700 and GJ 436.\n\nThe plot: every point (grey), medians in bins (blue), the box found (orange), and the depths of 1 Earth and 1 Jupiter radius before the star (its best-guess radius).',
    fr: 'Un transit dans la courbe de lumière TESS à la période d’un signal : chaque pic au FIP sous 10 % (en trouver un rendrait fort un signal faible ; autour de la conjonction de son repliement, là où une planète transite), et chaque planète qui transite ou TOI de l’étoile (à son éphéméride). La courbe vient des archives de l’étoile, sinon de MAST ; une étoile que TESS n’a pas observée (ses secteurs laissent des trous) n’en a pas, et le panneau le dit. Kepler, K2 et CoRoT sont cherchés aussi si l’étoile était dans leurs champs : la première mission avec un transit plausible est montrée, un bouton pour chacune.\n\nLa recherche : la courbe filtrée passe-haut (une médiane glissante de trois fois le transit attendu) et repliée ; une boîte cherchée à 3 sigma du moment attendu reporté à TESS, ou sur toute la phase. Quand l’erreur de la période étalerait les transits entre les secteurs, la période est balayée à 3 sigma et la courbe repliée à la meilleure.\n\nPlausible : 7 sigma ou plus, en 2 transits ou plus, encore 3,5 sans le plus profond (un événement n’est pas un transit), et au-delà de la meilleure boîte à des périodes où il n’y a rien (chance sous 1 %). Les autres planètes connues sont masquées d’abord. Testé sur LHS 1140, TOI-700 et GJ 436.\n\nLe graphique : chaque point (gris), les médianes par intervalle (bleu), la boîte trouvée (orange), et les profondeurs de 1 rayon terrestre et de 1 rayon de Jupiter devant l’étoile (son rayon le mieux estimé).',
  },
  batch_archives: {
    en: 'Every archive of each file’s star with it: the star found by its APERO name, its archives gathered as the Gather card does (DACE, CARMENES DR1, VizieR; not TESS), kept in the archives folder and read from there the next time (gather them again to ask anew). They are set apart from the file as the report does, and the star’s known planets go into the second pass. A file whose star is not found runs alone. About a minute more per star the first time.\n\nAPERO’s names: a copy of APERO’s astrometric database (a few MB from its assets server), made at the first file, again with Refresh APERO names.',
    fr: 'Toutes les archives de l’étoile de chaque fichier avec lui : l’étoile trouvée par son nom APERO, ses archives récupérées comme le fait la carte Récupérer (DACE, CARMENES DR1, VizieR ; pas TESS), gardées dans le dossier des archives et relues de là la fois suivante (cochez les récupérer à nouveau pour les redemander). Elles sont mises à part du fichier comme le fait le rapport, et les planètes connues de l’étoile entrent dans le second passage. Un fichier dont l’étoile n’est pas trouvée tourne seul. Environ une minute de plus par étoile la première fois.\n\nLes noms APERO : une copie de la base astrométrique d’APERO (quelques Mo de son serveur d’assets), faite au premier fichier, de nouveau avec Rafraîchir les noms APERO.',
  },
  series_zero: {
    en: 'The zero of each instrument on the plot of the series. Fitted offsets: the series about one solution, the same for every instrument: the fold shown on the series (its sinusoid or Keplerian orbit and its trend, when its box is ticked), else the fit of the quick look (its trend and the signals it found, once it has run); drawn once with its 1-sigma envelope (grey: draws from the full covariance of its fit), its signal over the time shown once few enough of its cycles are shown to see them (zoom in; before, its trend alone, dashed). Each instrument sits about it: the median of its velocities minus the solution, over its own epochs, is zero, so that the instruments run on from one to the next on one trend and one orbit, not each on its own piece of a solution. Under the series, the residuals to the solution, each instrument in its colour. Median: each about its own median (before the quick FIP, or without a trend). Kept; the PDF follows it.',
    fr: 'Le zéro de chaque instrument sur le graphique de la série. Offsets ajustés : la série autour d’une solution, la même pour tous les instruments : le repliement montré sur la série (sa sinusoïde ou son orbite képlérienne et sa tendance, quand sa case est cochée), sinon l’ajustement du coup d’œil (sa tendance et les signaux trouvés, une fois fait) ; tracée une fois avec son enveloppe à 1 sigma (en gris : des tirages dans la covariance complète de son ajustement), son signal sur le temps montré une fois assez peu de ses cycles montrés pour les voir (zoomez ; avant, sa tendance seule, en tirets). Chaque instrument est autour d’elle : la médiane de ses vitesses moins la solution, sur ses propres époques, est nulle, pour que les instruments se suivent sur une seule tendance et une seule orbite, pas chacun sur son propre morceau de solution. Sous la série, les résidus à la solution, chaque instrument dans sa couleur. Médiane : chacun autour de sa médiane (avant le FIP rapide, ou sans tendance). Gardé ; le PDF le suit.',
  },
  mstar: {
    en: 'The mass of the star, for the minimum mass of a companion (m sin i) under a Keplerian fold, and its mass under a transit fold (sin i ≈ 1). Filled by the resolver: the NASA Exoplanet Archive’s mass when the star hosts known planets, else a rough one from SIMBAD’s spectral type (Pecaut & Mamajek 2013, a dwarf; 10 % taken as its error). Type your own (a better one from the literature, say): it is kept, and the masses follow at once. m sin i is solved exactly (the companion’s mass not neglected against the star’s), its error from draws of K, P, e and this mass; in Earth, Neptune and Jupiter masses (IAU 2015 nominal values; Neptune from JPL).',
    fr: 'La masse de l’étoile, pour la masse minimale d’un compagnon (m sin i) sous un repliement képlérien, et sa masse sous un repliement de transit (sin i ≈ 1). Remplie par le résolveur : la masse de la NASA Exoplanet Archive si l’étoile a des planètes connues, sinon une masse grossière tirée du type spectral de SIMBAD (Pecaut & Mamajek 2013, une naine ; 10 % pris comme erreur). Tapez la vôtre (une meilleure tirée de la littérature, par exemple) : elle est gardée, et les masses suivent aussitôt. m sin i est résolue exactement (la masse du compagnon non négligée devant celle de l’étoile), son erreur tirée de tirages de K, P, e et de cette masse ; en masses de la Terre, de Neptune et de Jupiter (valeurs nominales de l’IAU 2015 ; Neptune selon le JPL).',
  },
  fip_view: {
    en: 'What the plot of the FIP shows. The period alone, and with its aliases: the FIP of the series (grey, of the period alone; white, of the period or any of its aliases, what decides on a planet). Each instrument too: the same, with the FIP of each instrument on its own (its nightly means, its own jitter) in its colour, when there are several (a signal every instrument sees is the star’s or a planet’s; one only one instrument sees may be that instrument’s). Kept from one visit to the next; the PDF follows it.',
    fr: 'Ce que montre le graphique du FIP. La période seule, et avec ses alias : le FIP de la série (gris, de la période seule ; blanc, de la période ou d’un de ses alias, ce qui décide d’une planète). Chaque instrument aussi : le même, avec le FIP de chaque instrument seul (ses moyennes par nuit, sa propre gigue) dans sa couleur, quand il y en a plusieurs (un signal que voient tous les instruments est celui de l’étoile ou d’une planète ; un signal qu’un seul voit peut venir de cet instrument). Gardé d’une visite à l’autre ; le PDF le suit.',
  },
  transit_fold: {
    en: 'A transiting planet (a known planet with a time of transit in the archive, or a TOI of TESS) folds on its transit ephemeris (◐, green): its period, and phase 0 at the transit itself, not at a conjunction the velocities put. The transit nearest the middle of the velocities is T0 = tc + N P, its error the ephemeris carried there, sqrt(σ_tc² + N² σ_P²), in time and in phase. The sinusoid then has its phase fixed by the transits (K, the offsets and the trend fitted); the conjunction the velocities put on their own (the phase free) is set against T0, in hours and in σ: a large offset is a wrong ephemeris, an eccentric orbit, or a signal that is not the planet’s. With Keplerian, koloa’s fit takes the period and the time of transit as gaussian priors (the eccentricity free). Several ephemerides of a planet (the archive’s, TESS’s) can be compared.',
    fr: 'Une planète en transit (une planète connue avec un temps de transit dans l’archive, ou un TOI de TESS) se replie sur son éphéméride de transit (◐, vert) : sa période, et la phase 0 au transit lui-même, pas à une conjonction que placent les vitesses. Le transit le plus proche du milieu des vitesses est T0 = tc + N P, son erreur l’éphéméride portée jusque-là, sqrt(σ_tc² + N² σ_P²), en temps et en phase. La sinusoïde a alors sa phase fixée par les transits (K, les offsets et la tendance ajustés) ; la conjonction que placent les vitesses seules (la phase libre) est comparée à T0, en heures et en σ : un grand écart est une éphéméride fausse, une orbite excentrique, ou un signal qui n’est pas celui de la planète. En képlérienne, l’ajustement de koloa prend la période et le temps de transit comme priors gaussiens (l’excentricité libre). Plusieurs éphémérides d’une planète (celle de l’archive, celle de TESS) peuvent être comparées.',
  },
  force_fold: {
    en: 'Fold at a period of your own: click the periodogram (the fold is at the dip of the FIP of the period alone nearest the click, within half a peak’s width, so that a click near a peak folds at its top), or type a period and press Fold (that period exactly). It gets the next number, a dashed button, and everything a numbered fold has: its colours, its Keplerian orbit, its solution on the series, the subtraction, the PDF. A known planet (★, orange) folds at its published period whatever the FIP shows there, with its published orbit (the archive’s default solution: P, K, e, ω, its epoch) dashed in orange on the fold and, with the solution on the time series, on the series (the fold’s offsets and trend). Every fit has its 1σ envelope in light grey: from the covariance of the sinusoid, or from draws of the Keplerian orbit.',
    fr: 'Replier à une période de votre choix : cliquez le périodogramme (le repliement est au creux du FIP de la période seule le plus proche du clic, à moins d’une demi-largeur de pic, pour qu’un clic près d’un pic replie à son sommet), ou tapez une période et appuyez sur Replier (cette période exactement). Il prend le numéro suivant, un bouton en tirets, et tout ce qu’a un repliement numéroté : ses couleurs, son orbite képlérienne, sa solution sur la série, la soustraction, le PDF. Une planète connue (★, orange) se replie à sa période publiée quoi que montre le FIP, avec son orbite publiée (la solution par défaut de l’archive : P, K, e, ω, son époque) en tirets orange sur le repliement et, avec la solution sur la série temporelle, sur la série (avec les offsets et la tendance du repliement). Chaque ajustement a son enveloppe à 1σ en gris pâle : de la covariance de la sinusoïde, ou de tirages de l’orbite képlérienne.',
  },
  fold_model: {
    en: 'The model of the fold. Sinusoid: a circular orbit, a weighted least-squares fit (each night weighted by its probability of being valid), instant. Keplerian (e free): koloa’s own fit, as in the report: the eccentricity free, the period free within the peak (half its width each side), the outliers of each night and the jitter of each instrument fitted, the trend; P, K and e with their errors (draws of the Laplace covariance), ω and the rms; a few seconds, once per fold. Under it, ΔBIC: the BIC of no planet, and of a circular orbit, minus the orbit’s, on the same nights, noise and trend (positive: the orbit is the better model; against a circular orbit, 2 to 6 favours the eccentricity, 6 to 10 strongly, above 10 very strongly, the scale of Kass & Raftery). Its curve, its points (about its offsets and trend) and its probabilities replace the sinusoid’s; the solution on the series, the subtraction and the PDF follow the choice.',
    fr: 'Le modèle du repliement. Sinusoïde : une orbite circulaire, un ajustement par moindres carrés pondérés (chaque nuit pondérée par sa probabilité d’être valide), instantané. Képlérienne (e libre) : l’ajustement de koloa, comme dans le rapport : l’excentricité libre, la période libre dans le pic (une demi-largeur de chaque côté), les valeurs aberrantes de chaque nuit et la gigue de chaque instrument ajustées, la tendance ; P, K et e avec leurs erreurs (des tirages de la covariance de Laplace), ω et la dispersion ; quelques secondes, une fois par repliement. Dessous, ΔBIC : le BIC sans planète, et celui d’une orbite circulaire, moins celui de l’orbite, sur les mêmes nuits, le même bruit et la même tendance (positif : l’orbite est le meilleur modèle ; face à une orbite circulaire, 2 à 6 favorise l’excentricité, 6 à 10 fortement, plus de 10 très fortement, l’échelle de Kass et Raftery). Sa courbe, ses points (autour de ses offsets et de sa tendance) et ses probabilités remplacent ceux de la sinusoïde ; la solution sur la série, la soustraction et le PDF suivent le choix.',
  },
  subtract: {
    en: 'Tick a signal (the box before its button) to take it out of the velocities: its sinusoid, or its Keplerian orbit when Keplerian is chosen as it is ticked (not its offsets and trend). Several can be ticked.\n\nEach fold is then fitted on the velocities without the other ticked signals, and shown so. With the solution on the time series, the sum of the ticked signals is drawn there.\n\nThe FIP of what is left is computed after each subtraction, in the order ticked, and drawn over the FIP of the series in the colour of its tick: a second signal the first one hid shows up, or nothing does. Click a peak of the last one to fold there, and tick it to go on; untick to put a signal back.\n\nTake the Keplerian of an eccentric orbit: its sinusoid leaves its harmonic at half the period. A result remembered keeps its ticks and their FIPs.',
    fr: 'Cochez un signal (la case devant son bouton) pour le retirer des vitesses : sa sinusoïde, ou son orbite képlérienne si Képlérienne est choisi au moment de cocher (pas ses offsets ni sa tendance). Plusieurs peuvent être cochés.\n\nChaque repliement est alors ajusté sur les vitesses sans les autres signaux cochés, et montré ainsi. Avec la solution sur la série temporelle, la somme des signaux cochés y est tracée.\n\nLe FIP de ce qui reste est calculé après chaque soustraction, dans l’ordre des coches, et tracé sur le FIP de la série dans la couleur de sa coche : un second signal que le premier cachait apparaît, ou rien. Cliquez un pic du dernier pour y replier, et cochez-le pour continuer ; décochez pour remettre un signal.\n\nPrenez la képlérienne d’une orbite excentrique : sa sinusoïde laisse son harmonique à la moitié de la période. Un résultat retenu garde ses coches et leurs FIP.',
  },
  series_colour: {
    en: 'The colour of the points of the series: their instrument, or their BERV, the barycentric Earth radial velocity in km/s (red when positive, blue when negative, on one scale for every instrument). The BERV is that of your LBL files (their BERV column: SPIRou, NIRPS...) and of the archives gathered (DACE’s cal_berv, CARMENES’ berv); an instrument without it is grey, and the button is off when no point has one. A drift or jumps that follow the BERV (the season, the year) rather than time are the Earth’s, telluric residuals, not the star’s. Hover over a point for its date and BERV. Kept from one visit to the next; the PDF follows it.',
    fr: 'La couleur des points de la série : leur instrument, ou leur BERV, la vitesse radiale barycentrique de la Terre en km/s (rouge si elle est positive, bleu si négative, sur une même échelle pour tous les instruments). La BERV est celle de vos fichiers LBL (leur colonne BERV : SPIRou, NIRPS...) et des archives récupérées (cal_berv de DACE, berv de CARMENES) ; un instrument qui ne l’a pas est gris, et le bouton est désactivé si aucun point ne l’a. Une dérive ou des sauts qui suivent la BERV (la saison, l’année) plutôt que le temps sont ceux de la Terre, des résidus telluriques, pas ceux de l’étoile. Survolez un point pour voir sa date et sa BERV. Gardée d’une visite à l’autre ; le PDF la suit.',
  },
  fold_colour: {
    en: 'The colour of the points of the fold: their instrument (the colours of the plot of the series), their date (one scale, from the first night to the last, with its calendar), or their BERV, the barycentric Earth radial velocity in km/s (red when positive, blue when negative). The BERV is that of your LBL files (their BERV column: SPIRou, NIRPS...) and of the archives gathered (DACE’s cal_berv, CARMENES’ berv); an instrument without it is grey, and the button is off when no point has one. A signal that follows the BERV more than the phase is the Earth’s (telluric residuals), not the star’s. Kept from one fold and one visit to the next; the PDF follows it.',
    fr: 'La couleur des points du repliement : leur instrument (les couleurs du graphique de la série), leur date (une échelle, de la première nuit à la dernière, avec son calendrier) ou leur BERV, la vitesse radiale barycentrique de la Terre en km/s (rouge si elle est positive, bleu si négative). La BERV est celle de vos fichiers LBL (leur colonne BERV : SPIRou, NIRPS...) et des archives récupérées (cal_berv de DACE, berv de CARMENES) ; un instrument qui ne l’a pas est gris, et le bouton est désactivé si aucun point ne l’a. Un signal qui suit la BERV plus que la phase est celui de la Terre (des résidus telluriques), pas celui de l’étoile. Gardée d’un repliement et d’une visite à l’autre ; le PDF la suit.',
  },
  overlay: {
    en: 'Draw the solution of the fold shown on the plot of the series above: one curve for every instrument, its trend and its sinusoid (or Keplerian orbit) at that period, with its 1-sigma envelope; each instrument about it (the median of its velocities minus it, zero), the residuals to it under the series. Another fold, another solution. With signals ticked, it is the sum of the ticked ones. Where the series spans too many periods to draw each cycle, a band shows the range of the solution (what its ups and downs fill), its trend dashed; zoom on the time (its slider) to see the cycles. The PDF draws it too.',
    fr: 'Tracer la solution du repliement montré sur le graphique de la série au-dessus : une seule courbe pour tous les instruments, sa tendance et sa sinusoïde (ou son orbite képlérienne) à cette période, avec son enveloppe à 1 sigma ; chaque instrument autour d’elle (la médiane de ses vitesses moins elle, nulle), les résidus à elle sous la série. Un autre repliement, une autre solution. Avec des signaux cochés, c’est la somme des signaux cochés. Quand la série couvre trop de périodes pour tracer chaque cycle, une bande montre l’étendue de la solution (ce que remplissent ses montées et descentes), sa tendance en tirets ; zoomez sur le temps (son curseur) pour voir les cycles. Le PDF la trace aussi.',
  },
  remember: {
    en: 'Keep this quick look, to come back to it: the page as it is (the star, the files, the archives ticked, the instruments left out, the options of the report, the ranges of the plots), the quick FIP (joint and of each instrument, its numbered peaks and folds), the velocities plotted, and a copy of the files and of the velocities of the star’s archives (not the TESS photometry: the report fetches its own). It goes to the Remembered targets tab, with the note written beside the button (why it is worth keeping, say). Kept on this machine, in ~/.cache/koloa/remembered. The quick FIP must have ended, and be of what is shown (recompute it after a change).',
    fr: 'Garder ce coup d’œil, pour y revenir : la page telle qu’elle est (l’étoile, les fichiers, les archives cochées, les instruments écartés, les options du rapport, les plages des graphiques), le FIP rapide (commun et de chaque instrument, ses pics numérotés et ses repliements), les vitesses tracées, et une copie des fichiers et des vitesses des archives de l’étoile (pas la photométrie TESS : le rapport va chercher la sienne). Il va dans l’onglet Cibles retenues, avec la note écrite à côté du bouton (pourquoi il vaut d’être gardé, par exemple). Gardé sur cette machine, dans ~/.cache/koloa/remembered. Le FIP rapide doit être terminé, et être de ce qui est montré (refaites-le après un changement).',
  },
  remembered: {
    en: 'The quick looks remembered with Remember the result, the last first: the star and its note, when, the nights of each instrument, the numbered peaks of the joint FIP, the known planets, the best peak of each instrument alone, and the data. Recall takes you back to the Analysis tab with the page as it was, nothing recomputed: the velocities, the FIP, the folds, the command line of the report. A file is the original when it has not changed since, else its copy; the archives likewise. Forget removes the result and its copies.',
    fr: 'Les coups d’œil retenus avec Retenir le résultat, le dernier en premier : l’étoile et sa note, quand, les nuits de chaque instrument, les pics numérotés du FIP commun, les planètes connues, le meilleur pic de chaque instrument seul, et les données. Rappeler vous ramène à l’onglet Analyse avec la page telle qu’elle était, sans rien recalculer : les vitesses, le FIP, les repliements, la ligne de commande du rapport. Un fichier est l’original s’il n’a pas changé depuis, sinon sa copie ; de même pour les archives. Oublier efface le résultat et ses copies.',
  },
  arch_plot: {
    en: 'The velocities that the Gather card put on disk for the star (in the archives folder above), to plot and to look at: DACE (with your key, what your account may see) and CARMENES DR1, each with the number of its points. These are the same boxes as DACE and CARMENES DR1 of the detailed report: what is plotted is what the report will use. With no file, the archives gathered are ticked when you plot, as they are then all there is. The DACE copy of an instrument of your files (NIRPS_DACE: the same spectra through DACE’s pipeline, less precise than LBL’s) is left out by default, and said: tick it to add it.',
    fr: 'Les vitesses que la carte Récupérer a mises sur le disque pour l’étoile (dans le dossier des archives ci-dessus), à tracer et à regarder : DACE (avec votre clé, ce que votre compte peut voir) et CARMENES DR1, chacune avec son nombre de points. Ce sont les mêmes cases que DACE et CARMENES DR1 du rapport détaillé : ce qui est tracé est ce que le rapport utilisera. Sans fichier, les archives récupérées sont cochées quand vous tracez, puisqu’elles sont alors tout ce qu’il y a. La copie DACE d’un instrument de vos fichiers (NIRPS_DACE : les mêmes spectres par le pipeline de DACE, moins précis que LBL) est écartée par défaut, et c’est dit : cochez-la pour l’ajouter.',
  },
  fold: {
    en: 'The nightly means folded at a numbered peak of the quick FIP (#1 the strongest; one number per family of aliases). A sinusoid is fitted there with an offset per instrument and the trend, each night weighted by its probability of being valid; phase 0 at the conjunction (the velocity falling through zero). Click another number to fold there.\n\nThe y range is the series’ above. A white circle marks a night the FIP took for an outlier more often than not; a red triangle at an edge, a night beyond the range. Hover over a point for its probability, date and BERV.',
    fr: 'Les moyennes par nuit repliées à un pic numéroté du FIP rapide (#1 le plus fort ; un numéro par famille d’alias). Une sinusoïde y est ajustée avec un offset par instrument et la tendance, chaque nuit pondérée par sa probabilité d’être valide ; phase 0 à la conjonction (la vitesse qui descend à zéro). Cliquez un autre numéro pour replier là.\n\nLa plage en y est celle de la série au-dessus. Un cercle blanc marque une nuit que le FIP a prise pour aberrante plus souvent que non ; un triangle rouge au bord, une nuit hors de la plage. Survolez un point pour sa probabilité, sa date et son BERV.',
  },
  pdf: {
    en: 'A PDF of what the page shows, as a LaTeX document: the velocities with the ranges of the plot, the quick FIP with its period range, its numbered peaks, the window and the known planets, and the folds at the numbered peaks; then a page of numbers: the star, every instrument (its file and path, points, nights, dates, errors, scatter, the noise added for the FIP), what is shown, the settings of the FIP, the table of its peaks with their K, the known planets and the peak that finds each, and the command line of the detailed report. Downloaded as koloa_quicklook_<star>.pdf; without pdflatex, the figures alone.',
    fr: 'Un PDF de ce que montre la page, en document LaTeX : les vitesses avec les intervalles du graphique, le FIP rapide avec son intervalle de périodes, ses pics numérotés, la fenêtre et les planètes connues, et les repliements aux pics numérotés ; puis une page de chiffres : l’étoile, chaque instrument (son fichier et son chemin, ses points, ses nuits, ses dates, ses erreurs, sa dispersion, le bruit ajouté pour le FIP), ce qui est montré, les réglages du FIP, la table de ses pics avec leur K, les planètes connues et le pic qui trouve chacune, et la ligne de commande du rapport détaillé. Téléchargé sous le nom koloa_quicklook_<étoile>.pdf ; sans pdflatex, les figures seules.',
  },
  latex: {
    en: 'Write the report in LaTeX and compile it to PDF (pdflatex). Unticked: the text report, the summary (JSON) and every figure are written all the same.',
    fr: 'Écrire le rapport en LaTeX et le compiler en PDF (pdflatex). Décochée : le rapport texte, le résumé (JSON) et toutes les figures sont écrits quand même.',
  },
};

// the info button of a field
const info = (key) => `<button type="button" class="info" data-help="${key}" aria-label="?">i</button>`;

let tip = null;
function showHelp(btn) {
  const help = HELP[btn.dataset.help];
  if (!help) return;
  if (!tip) {
    tip = document.createElement('div');
    tip.id = 'tip';
    tip.setAttribute('role', 'tooltip');
    document.body.appendChild(tip);
  }
  const text = help[lang] || help.en;
  tip.textContent = text;
  tip.style.display = 'block';
  const box = btn.getBoundingClientRect();
  // wider for a longer text (420 to 760 px), within the window
  const width = Math.min(Math.max(420, Math.min(760, 300 + 0.3 * text.length)), window.innerWidth - 24);
  tip.style.width = `${width}px`;
  let left = Math.min(Math.max(12, box.left - 20), window.innerWidth - width - 12);
  let top = box.bottom + 8;
  if (top + tip.offsetHeight > window.innerHeight - 8) top = Math.max(8, box.top - tip.offsetHeight - 8);
  tip.style.left = `${left + window.scrollX}px`;
  tip.style.top = `${top + window.scrollY}px`;
}
function hideHelp() { if (tip) tip.style.display = 'none'; }
document.addEventListener('mouseover', (e) => { const b = e.target.closest('.info'); if (b) showHelp(b); });
document.addEventListener('mouseout', (e) => { if (e.target.closest('.info')) hideHelp(); });
document.addEventListener('focusin', (e) => { const b = e.target.closest('.info'); if (b) showHelp(b); });
document.addEventListener('focusout', (e) => { if (e.target.closest('.info')) hideHelp(); });
document.addEventListener('keydown', (e) => { if (e.key === 'Escape') hideHelp(); });
document.addEventListener('click', (e) => { const b = e.target.closest('.info'); if (b) { e.preventDefault(); showHelp(b); } });

// -----------------------------------------------------------------------------
// the options of a detailed report, the same on every page
// -----------------------------------------------------------------------------
const OPTIONS = {
  fields: [
    { key: 'kmax', label: 'kmax', kind: 'number', value: '3', attrs: 'min="1" max="8"' },
    { key: 'nsweep', label: 'nsweep', kind: 'number', value: '2000', attrs: 'min="100" step="100"' },
    { key: 'nburn', label: 'nburn', kind: 'number', value: '400', attrs: 'min="50" step="50"' },
    { key: 'pmin', label: 'pmin', kind: 'number', value: '1.1', attrs: 'step="0.1"' },
    { key: 'pmax', label: 'pmax', kind: 'number', value: '', attrs: 'step="1"' },
    { key: 'periods', label: 'periods', kind: 'text', value: '' },
    { key: 'fip_gp', label: 'fip_gp_sel', help: 'fip_gp', kind: 'select', id: 'fipgp',
      choices: [['banded', 'gp_banded'], ['sho', 'gp_sho'], ['none', 'gp_none']] },
    { key: 'rotation', label: 'rotation_p', help: 'rotation', kind: 'number', id: 'rotation', value: '', attrs: 'step="any" min="0"' },
    { key: 'tois', label: 'tois', kind: 'text', id: 'tois', value: '' },
    { key: 'exclude', label: 'exclude', kind: 'text', id: 'exclude', value: '' },
    { key: 'include', label: 'include', kind: 'text', id: 'include', value: '' },
    { key: 'detection_map', label: 'dmap', help: 'detection_map', kind: 'select',
      choices: [['none', 'dmap_none'], ['fip', 'dmap_fip'], ['search', 'dmap_search']] },
  ],
  checks: [
    { key: 'toi_on', label: 'toi_on', id: 'toi_on', checked: false },
    { key: 'trend', label: 'trend', checked: true },
    { key: 'curvature', label: 'curvature', checked: false },
    { key: 'exposures', label: 'exposures', checked: false },
    { key: 'dace', text: 'DACE', checked: false },
    { key: 'carmenes', text: 'CARMENES DR1', checked: false },
    { key: 'tess', text: 'TESS', help: 'tess_opt', checked: true },
    { key: 'vizier', text: 'VizieR', checked: false },
    { key: 'archive', label: 'archive', help: 'archive_opt', checked: true },
    { key: 'gpcheck', label: 'gpcheck', checked: true },
    { key: 'duck', label: 'duck', checked: true },
    { key: 'mcmc', text: 'MCMC', checked: false },
    { key: 'latex', label: 'latex', checked: true },
  ],
};

// the fields and the boxes, for the options of `forName` (detailed)
function renderOptions(forName) {
  const field = (o) => {
    const id = o.id ? ` id="${o.id}"` : '';
    const name = `<span><span data-i18n="${o.label}">${esc(t(o.label))}</span> ${info(o.help || o.key)}</span>`;
    if (o.kind === 'select') {
      return `<label>${name}<select${id} data-opt="${o.key}" data-for="${forName}">`
        + o.choices.map(([val, txt]) => `<option value="${val}" data-i18n="${txt}">${esc(t(txt))}</option>`).join('')
        + '</select></label>';
    }
    return `<label>${name}<input${id} data-opt="${o.key}" data-for="${forName}" type="${o.kind}" value="${o.value}" ${o.attrs || ''} autocomplete="off"></label>`;
  };
  const check = (o) => {
    const id = o.id ? ` id="${o.id}"` : '';
    const text = o.label ? `<span data-i18n="${o.label}">${esc(t(o.label))}</span>` : esc(o.text);
    return `<span class="check"><label><input type="checkbox"${id} data-opt="${o.key}" data-for="${forName}"${o.checked ? ' checked' : ''}> ${text}</label>${info(o.help || o.key)}</span>`;
  };
  return `<div class="fields">${OPTIONS.fields.map(field).join('')}</div>`
    + `<div class="checks">${OPTIONS.checks.map(check).join('')}</div>`;
}

// the options of `forName` as the server reads them
function readOptions(forName) {
  const opts = {};
  document.querySelectorAll(`[data-for="${forName}"]`).forEach((el) => {
    opts[el.dataset.opt] = el.type === 'checkbox' ? el.checked : el.value;
  });
  return opts;
}

// -----------------------------------------------------------------------------
// what every page does the same way
// -----------------------------------------------------------------------------
document.addEventListener('click', async (e) => {
  if (e.target.id === 'lang') {
    lang = lang === 'fr' ? 'en' : 'fr';
    try { localStorage.setItem('koloa-lang', lang); } catch (err) { /* no storage */ }
    applyLang();
  }
  const copy = e.target.closest('[data-copy]');
  if (copy) {
    const text = $(copy.dataset.copy).textContent;
    try { await navigator.clipboard.writeText(text); } catch (err) { /* old browsers */ }
    copy.textContent = t('copied');
    setTimeout(() => { copy.textContent = t('copy'); }, 1200);
  }
  const stop = e.target.closest('[data-stop]');
  if (stop) { await api('/api/stop', { id: stop.dataset.stop }); poll(); }
});

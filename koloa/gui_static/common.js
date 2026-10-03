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
    file: 'File (optional)', root: 'Archives folder', plot: 'Plot', browse: 'Browse...', picking: 'choosing...',
    var_pending: 'SIMBAD’s periods of variability are on their way (its service is slow today)',
    p_valid: 'P(valid)', p_valid_low: 'P(valid) < 50 %',
    fip_prompt: 'Several instruments: untick those to leave out (fewer, faster), then compute the quick FIP; it can take minutes on a long series.',
    start_fip: 'Compute the quick FIP', stop_fip: 'Stop', fip_stopped: 'stopped',
    envelope: '1σ envelope', published: 'published', known_fold: 'a fold at the published period of this known planet, with its published orbit (dashed)',
    fipview_joint: 'the period alone, and with its aliases', fipview_each: 'each instrument too',
    transit_word: 'transit', transit_fold: 'a fold on this transit ephemeris: phase 0 at the transit, its uncertainty at the velocities',
    on_transit: 'on its transit ephemeris', cycles_from: 'periods from', phase_word: 'phase', phase_fixed: 'with its phase fixed by the transits',
    conj_alone: 'the velocities alone put the conjunction',
    fold_at: 'Fold at P =', fold_go: 'Fold', or_click: 'or click the periodogram', asked: 'a period asked', folding: 'folding...',
    fmodel_sine: 'sinusoid', fmodel_kepler: 'Keplerian (e free)', fmodel_kepler_short: 'Keplerian',
    fitting_kepler: 'fitting the Keplerian orbit (outliers, jitters, trend; a few seconds)...',
    fold_note_kep: 'nightly means; koloa’s fit, the eccentricity and period free within the peak; phase 0 at the conjunction',
    subtract: 'Subtract this solution, FIP of the residuals', residuals_of: 'FIP of the residuals, without:', back_series: 'Back to the series',
    colour_by: 'Colour by', fcol_inst: 'instrument', fcol_date: 'date', overlay: 'the solution on the time series', offscale: 'beyond the range',
    accel_star: 'acceleration of the star', no_trend_fit: 'no trend fitted (the box “acceleration (trend)” of the report is unticked)',
    accel_note: 'the fit of the quick look, perspective acceleration included',
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
    use_sho: 'use for the SHO', no_sho: 'none (no SHO)',
    new_target: 'New target', known: 'Known planets (NASA Exoplanet Archive)', none_known: 'none in the archive',
    toi_title: 'TESS Objects of Interest', toi_pick: 'click a TOI to fit it with the ephemeris of TESS', toi_none: 'none',
    toi_on: 'TESS ephemerides (TOIs)', tois: 'TOIs (empty: all)', transit: 'transit',
    from_disk: 'from the disk', disk_archives: 'the archives folder', disk_kept: 'the copy kept', asked_now: 'asked the network now', refresh_star: 'Refresh',
    solution: 'solution', solutions: 'solutions',
    files_in: 'Files (optional): one per instrument or reduction', add_file: '+ Add a file', inst_auto: 'instrument (auto)', remove: 'remove',
    gather: 'Gather the archives',
    gather_hint: 'What DACE (with your key when there is one), CARMENES DR1 and TESS have of the star, kept in the archives folder, one folder per star.',
    copy: 'Copy', copied: 'Copied', run_gather: 'Gather', detailed: 'Detailed report',
    detailed_hint: 'Everything koloa can say about the star, as a PDF report: from the file, from the SIMBAD name (the archives only), or (best) from both.',
    outdir: 'Report folder', kmax: 'Signals (kmax)', nsweep: 'Sweeps', nburn: 'Burn-in',
    pmin: 'Shortest period [d]', pmax: 'Longest period [d]', periods: 'More periods to test [d]',
    dmap: 'Detection map', dmap_none: 'none', dmap_fip: 'by the FIP (hours)', dmap_search: 'blind search (minutes)',
    fip_gp: 'GP in the FIP, by band', exposures: 'every exposure (not nightly means)',
    archive: 'Exoplanet Archive', gpcheck: 'signals against a GP', duck: 'duck test', latex: 'PDF report',
    run_detailed: 'Make the report', runs: 'Runs',
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
    file: 'Fichier (facultatif)', root: 'Dossier des archives', plot: 'Tracer', browse: 'Parcourir...', picking: 'choix en cours...',
    var_pending: 'les périodes de variabilité de SIMBAD arrivent (son service est lent aujourd’hui)',
    p_valid: 'P(valide)', p_valid_low: 'P(valide) < 50 %',
    fip_prompt: 'Plusieurs instruments : décochez ceux à écarter (moins, plus vite), puis calculez le FIP rapide ; il peut prendre des minutes sur une longue série.',
    start_fip: 'Calculer le FIP rapide', stop_fip: 'Arrêter', fip_stopped: 'arrêté',
    envelope: 'enveloppe à 1σ', published: 'publiée', known_fold: 'un repliement à la période publiée de cette planète connue, avec son orbite publiée (en tirets)',
    fipview_joint: 'la période seule, et avec ses alias', fipview_each: 'chaque instrument aussi',
    transit_word: 'transit', transit_fold: 'un repliement sur cette éphéméride de transit : phase 0 au transit, son incertitude aux vitesses',
    on_transit: 'sur son éphéméride de transit', cycles_from: 'périodes depuis', phase_word: 'phase', phase_fixed: 'sa phase fixée par les transits',
    conj_alone: 'les vitesses seules placent la conjonction',
    fold_at: 'Replier à P =', fold_go: 'Replier', or_click: 'ou cliquez le périodogramme', asked: 'une période demandée', folding: 'repliement...',
    fmodel_sine: 'sinusoïde', fmodel_kepler: 'képlérienne (e libre)', fmodel_kepler_short: 'Képlérienne',
    fitting_kepler: 'ajustement de l’orbite képlérienne (valeurs aberrantes, gigues, tendance ; quelques secondes)...',
    fold_note_kep: 'moyennes par nuit ; l’ajustement de koloa, l’excentricité et la période libres dans le pic ; phase 0 à la conjonction',
    subtract: 'Soustraire cette solution, FIP des résidus', residuals_of: 'FIP des résidus, sans :', back_series: 'Revenir à la série',
    colour_by: 'Couleur selon', fcol_inst: 'instrument', fcol_date: 'date', overlay: 'la solution sur la série temporelle', offscale: 'hors de la plage',
    accel_star: 'accélération de l’étoile', no_trend_fit: 'pas de tendance ajustée (la case « accélération (tendance) » du rapport est décochée)',
    accel_note: 'l’ajustement du coup d’œil, accélération de perspective incluse',
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
    use_sho: 'utiliser pour le SHO', no_sho: 'aucune (pas de SHO)',
    new_target: 'Nouvelle cible', known: 'Planètes connues (NASA Exoplanet Archive)', none_known: 'aucune dans l’archive',
    toi_title: 'TESS Objects of Interest', toi_pick: 'cliquez un TOI pour l’ajuster avec l’éphéméride de TESS', toi_none: 'aucun',
    toi_on: 'éphémérides TESS (TOI)', tois: 'TOI (vide : tous)', transit: 'transit',
    from_disk: 'lu sur le disque', disk_archives: 'le dossier des archives', disk_kept: 'la copie gardée', asked_now: 'demandé au réseau à l’instant', refresh_star: 'Rafraîchir',
    solution: 'solution', solutions: 'solutions',
    files_in: 'Fichiers (facultatifs) : un par instrument ou réduction', add_file: '+ Ajouter un fichier', inst_auto: 'instrument (auto)', remove: 'retirer',
    gather: 'Récupérer les archives',
    gather_hint: 'Ce que DACE (avec votre clé s’il y en a une), CARMENES DR1 et TESS ont de l’étoile, rangé dans le dossier des archives, un dossier par étoile.',
    copy: 'Copier', copied: 'Copié', run_gather: 'Récupérer', detailed: 'Rapport détaillé',
    detailed_hint: 'Tout ce que koloa peut dire de l’étoile, en un rapport PDF : à partir du fichier, du nom SIMBAD (les archives seules), ou (le mieux) des deux.',
    outdir: 'Dossier du rapport', kmax: 'Signaux (kmax)', nsweep: 'Itérations', nburn: 'Rodage',
    pmin: 'Période minimale [j]', pmax: 'Période maximale [j]', periods: 'Autres périodes à tester [j]',
    dmap: 'Carte de détection', dmap_none: 'aucune', dmap_fip: 'par le FIP (heures)', dmap_search: 'recherche aveugle (minutes)',
    fip_gp: 'GP dans le FIP, par bande', exposures: 'chaque pose (pas les moyennes par nuit)',
    archive: 'Exoplanet Archive', gpcheck: 'signaux face à un GP', duck: 'duck test', latex: 'rapport PDF',
    run_detailed: 'Faire le rapport', runs: 'Exécutions',
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
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));


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
const openLogs = new Set();

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
      return `<li><span class="icon">${icon}</span><span>${esc(st.name)}</span><span class="time">${clock(st.elapsed)}</span>${detail}</li>`;
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
      <details data-log="${job.id}"${openLogs.has(job.id) ? ' open' : ''}><summary>${esc(t('log'))} (${job.lines.length})</summary>
        <pre class="codeblock log">${esc(job.lines.slice(-400).join('\n'))}</pre></details>
    </div>`;
  }).join('');
  box.querySelectorAll('details[data-log]').forEach((d) => {
    d.addEventListener('toggle', () => { if (d.open) openLogs.add(d.dataset.log); else openLogs.delete(d.dataset.log); });
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
    en: 'The name of the star as SIMBAD knows it: GJ 436, Ross 905, HD 69830, TOI-700 or TIC 307210830 all work. koloa resolves it (CDS Sesame) into its main identifier, its Gaia DR3 and TIC numbers and its aliases, and every archive (DACE, CARMENES DR1, the NASA Exoplanet Archive, the TOIs, TESS) is then searched with those. Empty: koloa guesses a name from the OBJECT column of the file, which SIMBAD may not know; give it whenever you can. Once asked, the answer is kept on disk.',
    fr: 'Le nom de l’étoile tel que SIMBAD le connaît : GJ 436, Ross 905, HD 69830, TOI-700 ou TIC 307210830 fonctionnent tous. koloa le résout (CDS Sesame) en son identifiant principal, ses numéros Gaia DR3 et TIC et ses alias, et chaque archive (DACE, CARMENES DR1, la NASA Exoplanet Archive, les TOI, TESS) est ensuite cherchée avec eux. Vide : koloa devine un nom à partir de la colonne OBJECT du fichier, que SIMBAD ne connaît pas toujours ; donnez-le dès que possible. Une fois demandée, la réponse est gardée sur le disque.',
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
    en: 'Which archives to gather. DACE: the velocities of HARPS, ESPRESSO, NIRPS, CORALIE..., the public ones, and with your key (DACE_API_KEY or ~/.dacerc) what your account may see. CARMENES DR1: the GTO velocities of 2016 to 2020, about 360 M dwarfs of the north, corrected for the nightly zero points. TESS: the light curve of every sector at MAST (SPOC 2-minute when there is one, else TESS-SPOC or QLP).',
    fr: 'Les archives à récupérer. DACE : les vitesses de HARPS, ESPRESSO, NIRPS, CORALIE..., les publiques, et avec votre clé (DACE_API_KEY ou ~/.dacerc) ce que votre compte peut voir. CARMENES DR1 : les vitesses du GTO de 2016 à 2020, environ 360 naines M du nord, corrigées des points zéro nocturnes. TESS : la courbe de lumière de chaque secteur à MAST (SPOC 2 minutes s’il y en a une, sinon TESS-SPOC ou QLP).',
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
    en: 'Instruments left out of the report, by their names in the plot, separated by spaces: NIRPS_DACE HARPS03 for instance. Untick an instrument in the table under the plot, or click it in the legend, to add it here. They are left out once the files, DACE, CARMENES and VizieR are put together. Empty: none.',
    fr: 'Les instruments écartés du rapport, par leurs noms dans le graphique, séparés par des espaces : NIRPS_DACE HARPS03 par exemple. Décochez un instrument dans le tableau sous le graphique, ou cliquez-le dans la légende, pour l’ajouter ici. Ils sont écartés une fois les fichiers, DACE, CARMENES et VizieR rassemblés. Vide : aucun.',
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
    en: 'Add the velocities published with the known planets of the star, found on VizieR from the bibcodes of their solutions in the archive. Off by default: nothing is fetched unless this is ticked. A published velocity within a minute of an exposure already in the series is the same spectrum, and is left out.',
    fr: 'Ajouter les vitesses publiées avec les planètes connues de l’étoile, trouvées sur VizieR à partir des bibcodes de leurs solutions dans l’archive. Désactivée par défaut : rien n’est demandé si elle n’est pas cochée. Une vitesse publiée à moins d’une minute d’une pose déjà dans la série est le même spectre, et elle est écartée.',
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
    en: 'A first look before the detailed report, made on its own when the plot holds one instrument; with several, once you have unticked those to leave out, with Compute the quick FIP (Stop ends it, while it runs or waits for its turn; a new one stops the one before). It is the FIP of exactly what the plot shows (your files, the archives ticked, the instruments not left out), outlier-aware, in nightly means, with no GP. As in the report, it runs twice: the errors of each instrument inflated to the noise of a fit without planets, a first FIP; then to the noise of a fit with the signals it found (FIP < 1 %) and the known planets, so that their variance is not taken for noise, and the FIP again. Two signals, 500 sweeps: a minute or so. White, the FIP of the period or any of its aliases (what decides on a planet; black in the PDF); grey, of the period alone; dotted, FIP = 1 % and the window (a day, a synodic month, a year); dashed, the known planets. The slider under it zooms on the periods. Without the GP of the activity, a peak at the rotation or its harmonics is expected: the report sorts that out. With several instruments, each one then has a FIP of its own (its nightly means, its own jitter sampled, at least 10 nights), drawn in its colour over the joint one, with its best peaks in the table below: a signal seen by every instrument is the star’s or a planet’s, one seen by a single instrument may be that instrument’s. Click an instrument in the legend to hide or show its curve. The trend follows the boxes of the report: an acceleration by default (acceleration (trend)), its change too with its change (curvature), none when the first is unticked; it is fitted in the noise, the FIP and the folds. Under the FIP, the acceleration of the star, dv/dt in m/s/yr (and d²v/dt² with the curvature), from the fit of the quick look (its signals and the known planets, the outliers and jitters of each instrument, no GP), with its Laplace errors and significance. It is what the velocities measure: the perspective acceleration (the proper motion times the distance, mu² d) is in it, not removed.',
    fr: 'Un premier coup d’œil avant le rapport détaillé, fait de lui-même quand le graphique a un seul instrument ; avec plusieurs, une fois décochés ceux à écarter, avec Calculer le FIP rapide (Arrêter l’interrompt, qu’il tourne ou attende son tour ; un nouveau arrête celui d’avant). C’est le FIP de ce que montre exactement le graphique (vos fichiers, les archives cochées, les instruments non écartés), robuste aux valeurs aberrantes, en moyennes par nuit, sans GP. Comme dans le rapport, il passe deux fois : les erreurs de chaque instrument gonflées au bruit d’un ajustement sans planète, un premier FIP ; puis au bruit d’un ajustement avec les signaux trouvés (FIP < 1 %) et les planètes connues, pour que leur variance ne soit pas prise pour du bruit, et le FIP à nouveau. Deux signaux, 500 itérations : une minute environ. En blanc, le FIP de la période ou d’un de ses alias (ce qui décide d’une planète ; en noir dans le PDF) ; en gris, celui de la période seule ; en pointillés, FIP = 1 % et la fenêtre (un jour, un mois synodique, un an) ; en tirets, les planètes connues. Le curseur dessous zoome sur les périodes. Sans le GP de l’activité, un pic à la rotation ou à ses harmoniques est attendu : le rapport fait le tri. Avec plusieurs instruments, chacun a ensuite son propre FIP (ses moyennes par nuit, sa propre gigue échantillonnée, au moins 10 nuits), tracé dans sa couleur par-dessus le FIP commun, avec ses meilleurs pics dans le tableau dessous : un signal que voient tous les instruments est celui de l’étoile ou d’une planète, un signal qu’un seul voit peut venir de cet instrument. Cliquez un instrument dans la légende pour cacher ou montrer sa courbe. La tendance suit les cases du rapport : une accélération par défaut (accélération (tendance)), sa variation aussi avec sa variation (courbure), aucune si la première est décochée ; elle est ajustée dans le bruit, le FIP et les repliements. Sous le FIP, l’accélération de l’étoile, dv/dt en m/s/an (et d²v/dt² avec la courbure), de l’ajustement du coup d’œil (ses signaux et les planètes connues, les valeurs aberrantes et gigues de chaque instrument, sans GP), avec ses erreurs de Laplace et sa signification. C’est ce que mesurent les vitesses : l’accélération de perspective (le mouvement propre fois la distance, mu² d) y est, non retirée.',
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
    en: 'The model of the fold. Sinusoid: a circular orbit, a weighted least-squares fit (each night weighted by its probability of being valid), instant. Keplerian (e free): koloa’s own fit, as in the report: the eccentricity free, the period free within the peak (half its width each side), the outliers of each night and the jitter of each instrument fitted, the trend; P, K and e with their errors (draws of the Laplace covariance), ω and the rms; a few seconds, once per fold. Its curve, its points (about its offsets and trend) and its probabilities replace the sinusoid’s; the solution on the series, the subtraction and the PDF follow the choice.',
    fr: 'Le modèle du repliement. Sinusoïde : une orbite circulaire, un ajustement par moindres carrés pondérés (chaque nuit pondérée par sa probabilité d’être valide), instantané. Képlérienne (e libre) : l’ajustement de koloa, comme dans le rapport : l’excentricité libre, la période libre dans le pic (une demi-largeur de chaque côté), les valeurs aberrantes de chaque nuit et la gigue de chaque instrument ajustées, la tendance ; P, K et e avec leurs erreurs (des tirages de la covariance de Laplace), ω et la dispersion ; quelques secondes, une fois par repliement. Sa courbe, ses points (autour de ses offsets et de sa tendance) et ses probabilités remplacent ceux de la sinusoïde ; la solution sur la série, la soustraction et le PDF suivent le choix.',
  },
  subtract: {
    en: 'Take the signal of the fold shown (its sinusoid or its Keplerian orbit, not its offsets and trend, which the FIP fits again) out of the velocities, and run the quick FIP on what is left: a second signal the first one hid shows up, or nothing does. Again on a fold of the residuals to take a second signal out, and so on. The FIP card says what was subtracted; Back to the series runs it on the velocities as they are. The plot of the series keeps the velocities themselves. A result remembered keeps its subtractions.',
    fr: 'Retirer des vitesses le signal du repliement montré (sa sinusoïde ou son orbite képlérienne, pas ses offsets et sa tendance, que le FIP ajuste de nouveau), et faire le FIP rapide sur ce qui reste : un second signal que le premier cachait apparaît, ou rien. De nouveau sur un repliement des résidus pour retirer un second signal, et ainsi de suite. La carte du FIP dit ce qui a été soustrait ; Revenir à la série le refait sur les vitesses telles quelles. Le graphique de la série garde les vitesses elles-mêmes. Un résultat retenu garde ses soustractions.',
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
    en: 'Draw the solution of the fold shown on the plot of the series above: for each instrument, its offset, the trend and the sinusoid at that period, about the instrument’s median as its points are, over the time it covers. Another fold, another solution. Zoom on the time (its slider) to see the sinusoid itself where the series spans many periods. The PDF draws it too.',
    fr: 'Tracer la solution du repliement montré sur le graphique de la série au-dessus : pour chaque instrument, son offset, la tendance et la sinusoïde à cette période, autour de la médiane de l’instrument comme le sont ses points, sur le temps qu’il couvre. Un autre repliement, une autre solution. Zoomez sur le temps (son curseur) pour voir la sinusoïde elle-même quand la série couvre beaucoup de périodes. Le PDF la trace aussi.',
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
    en: 'The velocities that the Gather card put on disk for the star (in the archives folder above), to plot and to look at: DACE (with your key, what your account may see) and CARMENES DR1, each with the number of its points. These are the same boxes as DACE and CARMENES DR1 of the detailed report: what is plotted is what the report will use. With no file, the archives gathered are ticked when you plot, as they are then all there is.',
    fr: 'Les vitesses que la carte Récupérer a mises sur le disque pour l’étoile (dans le dossier des archives ci-dessus), à tracer et à regarder : DACE (avec votre clé, ce que votre compte peut voir) et CARMENES DR1, chacune avec son nombre de points. Ce sont les mêmes cases que DACE et CARMENES DR1 du rapport détaillé : ce qui est tracé est ce que le rapport utilisera. Sans fichier, les archives récupérées sont cochées quand vous tracez, puisqu’elles sont alors tout ce qu’il y a.',
  },
  fold: {
    en: 'The nightly means folded at a numbered peak of the quick FIP (#1 is the strongest; the peaks below a FIP of 10 %, or the best three, are numbered, one per family of aliases, the true period of each family being the one whose own FIP is the lowest). A sinusoid is fitted at that period with an offset per instrument and the trend (weighted least squares, each night also weighted by its probability of being valid, so that an outlier hardly counts; the error of K scaled by the reduced chi-square when it is above one); the points are shown about their offsets, phase 0 at the conjunction (the velocity falling through zero). Click another number to fold there. The velocities are on the y range of the plot of the series above (its sliders and zoom move both). A white circle marks a night with less than a 50 % probability of being valid: the FIP (outlier-aware) took it for an outlier more often than not, and gave it little weight. Hover over a point for its probability, its date and its BERV. A red triangle at the top or bottom edge points to a night beyond the range shown (hover: its velocity).',
    fr: 'Les moyennes par nuit repliées à un pic numéroté du FIP rapide (#1 est le plus fort ; les pics sous un FIP de 10 %, ou les trois meilleurs, sont numérotés, un par famille d’alias, la vraie période de chaque famille étant celle dont le FIP propre est le plus bas). Une sinusoïde est ajustée à cette période avec un offset par instrument et la tendance (moindres carrés pondérés, l’erreur de K multipliée par le chi carré réduit s’il dépasse un ; chaque nuit pondérée aussi par sa probabilité d’être valide, pour qu’une valeur aberrante compte à peine) ; les points sont montrés autour de leurs offsets, la phase 0 à la conjonction (la vitesse descendant par zéro). Cliquez un autre numéro pour replier là. Les vitesses sont sur la plage en y du graphique de la série au-dessus (ses curseurs et son zoom déplacent les deux). Un cercle blanc marque une nuit dont la probabilité d’être valide est sous 50 % : le FIP (robuste aux valeurs aberrantes) l’a prise pour une valeur aberrante plus souvent qu’autrement, et lui a donné peu de poids. Survolez un point pour voir sa probabilité, sa date et sa BERV. Un triangle rouge au bord du haut ou du bas pointe vers une nuit hors de la plage montrée (au survol : sa vitesse).',
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
  tip.textContent = help[lang] || help.en;
  tip.style.display = 'block';
  const box = btn.getBoundingClientRect();
  const width = Math.min(420, window.innerWidth - 24);
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

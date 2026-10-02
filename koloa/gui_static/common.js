// koloa's GUI: what its pages share (texts, help, options, the file
// browser, the runs); koloa.gui serves it
'use strict';

const TEXT = {
  en: {
    tagline: 'Outlier-aware radial velocities · on this machine', docs: 'Docs',
    star: 'Star', star_hint: 'Its SIMBAD name: the archives, DACE, CARMENES and TESS find it by that name.',
    resolve: 'Resolve', velocities: 'Velocities',
    velocities_hint: 'A file of velocities (LBL .rdb, csv, DACE csv), the archives gathered for the star, or both. Untick an instrument (or click it in the legend) to leave it out of the report.',
    use: 'Used', exclude: 'Instruments left out', source: 'Source', src_file: 'input file',
    file: 'File (optional)', root: 'Archives folder', plot: 'Plot', browse: 'Browse...', picking: 'choosing...',
    tab_single: 'Single target', tab_batch: 'Batch', batch_tagline: 'Many targets at once: a script for this machine or a server',
    where: 'Where it runs', this_machine: 'this machine', a_server: 'a server (ssh)', host: 'Server', test: 'Test',
    workdir: 'Working folder', koloa_cmd: 'How koloa is called there', batchdir: 'Batch folder', jobs: 'Targets at a time',
    bashrc: 'load ~/.bashrc first', every_target: 'Options for every target',
    every_target_hint: 'The options of the detailed report, the same for every target of the batch; each report goes to its own folder in the batch.',
    targets: 'Targets', add_target: '+ Add a target', clear_targets: 'Clear the list', paste_title: 'Paste a list', add_lines: 'Add these lines',
    the_script: 'The script', make_script: 'Make the script', download: 'Download', save_here: 'Save it in the working folder', send: 'Send it to the server',
    to_start: 'To start it, and to follow it (koloa does not start it: you do):', saved: 'saved:', inst_head: 'instrument', choose_folder: 'Choose this folder', batch_run: 'Batch',
    n_targets: 'targets', name_col: 'SIMBAD name', files_col: 'files (0, 1 or more)', add_file_short: '+ file',
    sent: 'sent to', confirm_clear: 'Clear the whole list of targets?', testing: 'connecting...', need_host: 'Give the server first.',
    new_target: 'New target', known: 'Known planets (NASA Exoplanet Archive)', none_known: 'none in the archive',
    toi_title: 'TESS Objects of Interest', toi_pick: 'click a TOI to fit it with the ephemeris of TESS', toi_none: 'none',
    toi_on: 'TESS ephemerides (TOIs)', tois: 'TOIs (empty: all)', transit: 'transit',
    from_disk: 'from the disk', asked_now: 'asked the network now', refresh_star: 'Refresh',
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
    velocities_hint: 'Un fichier de vitesses (LBL .rdb, csv, csv de DACE), les archives récupérées pour l’étoile, ou les deux. Décochez un instrument (ou cliquez-le dans la légende) pour l’écarter du rapport.',
    use: 'Utilisé', exclude: 'Instruments écartés', source: 'Source', src_file: 'fichier d’entrée',
    file: 'Fichier (facultatif)', root: 'Dossier des archives', plot: 'Tracer', browse: 'Parcourir...', picking: 'choix en cours...',
    tab_single: 'Une cible', tab_batch: 'Lot', batch_tagline: 'Plusieurs cibles d’un coup : un script pour cette machine ou un serveur',
    where: 'Où ça tourne', this_machine: 'cette machine', a_server: 'un serveur (ssh)', host: 'Serveur', test: 'Tester',
    workdir: 'Dossier de travail', koloa_cmd: 'Comment koloa y est appelé', batchdir: 'Dossier du lot', jobs: 'Cibles à la fois',
    bashrc: 'charger ~/.bashrc d’abord', every_target: 'Options pour chaque cible',
    every_target_hint: 'Les options du rapport détaillé, les mêmes pour chaque cible du lot ; chaque rapport va dans son propre dossier du lot.',
    targets: 'Cibles', add_target: '+ Ajouter une cible', clear_targets: 'Vider la liste', paste_title: 'Coller une liste', add_lines: 'Ajouter ces lignes',
    the_script: 'Le script', make_script: 'Faire le script', download: 'Télécharger', save_here: 'L’enregistrer dans le dossier de travail', send: 'L’envoyer au serveur',
    to_start: 'Pour le lancer, et le suivre (koloa ne le lance pas : c’est vous) :', saved: 'enregistré :', inst_head: 'instrument', choose_folder: 'Choisir ce dossier', batch_run: 'Lot',
    n_targets: 'cibles', name_col: 'nom SIMBAD', files_col: 'fichiers (0, 1 ou plus)', add_file_short: '+ fichier',
    sent: 'envoyé à', confirm_clear: 'Vider toute la liste des cibles ?', testing: 'connexion...', need_host: 'Donnez d’abord le serveur.',
    new_target: 'Nouvelle cible', known: 'Planètes connues (NASA Exoplanet Archive)', none_known: 'aucune dans l’archive',
    toi_title: 'TESS Objects of Interest', toi_pick: 'cliquez un TOI pour l’ajuster avec l’éphéméride de TESS', toi_none: 'aucun',
    toi_on: 'éphémérides TESS (TOI)', tois: 'TOI (vide : tous)', transit: 'transit',
    from_disk: 'lu sur le disque', asked_now: 'demandé au réseau à l’instant', refresh_star: 'Rafraîchir',
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
    const what = { gather: t('gather'), archive: t('refresh_archive'), batch: t('batch_run') }[job.action] || t('detailed');
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
    en: 'A rotation period you trust [days], a published one: 2.704 for Wolf 359 for instance. It becomes the star’s rotation for every check (a signal at P_rot, P_rot/2 or 2 P_rot is flagged) and the prior of the GP check, and with SHO in the GP menu, the period of the GP inside the FIP (held). Click a period in the resolver to fill it in. Empty: the archive’s rotation, when it has one.',
    fr: 'Une période de rotation fiable [jours], publiée : 2,704 pour Wolf 359 par exemple. Elle devient la rotation de l’étoile pour chaque vérification (un signal à P_rot, P_rot/2 ou 2 P_rot est signalé) et l’a priori du test GP, et avec SHO dans le menu GP, la période du GP dans le FIP (fixée). Cliquez une période dans le résolveur pour la remplir. Vide : la rotation de l’archive, si elle en a une.',
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
    en: 'Ask DACE for the velocities of the star, from every instrument it holds. An instrument the files have too is set apart as <inst>_DACE (another pipeline, its own offset), and the spectra the files already have are not counted twice. Your DACE key, when there is one, adds what your account may see. DACE answers from some networks only; without it, the analysis goes on.',
    fr: 'Demander à DACE les vitesses de l’étoile, de tous les instruments qu’il a. Un instrument que les fichiers ont aussi est mis à part en <inst>_DACE (un autre pipeline, son propre offset), et les spectres que les fichiers ont déjà ne sont pas comptés deux fois. Votre clé DACE, s’il y en a une, ajoute ce que votre compte peut voir. DACE ne répond que depuis certains réseaux ; sans lui, l’analyse continue.',
  },
  carmenes: {
    en: 'Add the velocities of CARMENES DR1 (Ribas et al. 2023: the GTO of 2016 to 2020, about 360 M dwarfs of the north): SERVAL’s velocities corrected for the nightly zero points, with their activity indicators. Nothing is added for a star that is not in DR1.',
    fr: 'Ajouter les vitesses de CARMENES DR1 (Ribas et al. 2023 : le GTO de 2016 à 2020, environ 360 naines M du nord) : les vitesses de SERVAL corrigées des points zéro nocturnes, avec leurs indicateurs d’activité. Rien n’est ajouté pour une étoile absente de DR1.',
  },
  tess_opt: {
    en: 'Fetch the TESS light curves of the star (MAST) and look in each sector for a photometric peak at the period of every signal, at its half, third or double: a spot rotating with the star moves both its light and its velocities. Beyond about 13 days, a period is out of reach of a 27-day sector.',
    fr: 'Récupérer les courbes de lumière TESS de l’étoile (MAST) et chercher dans chaque secteur un pic photométrique à la période de chaque signal, à sa moitié, son tiers ou son double : une tache qui tourne avec l’étoile change à la fois sa lumière et ses vitesses. Au-delà d’environ 13 jours, une période est hors de portée d’un secteur de 27 jours.',
  },
  vizier: {
    en: 'Add the velocities published with the known planets of the star, found on VizieR from the bibcodes of their solutions in the archive. A published velocity within a minute of an exposure already in the series is the same spectrum, and is left out.',
    fr: 'Ajouter les vitesses publiées avec les planètes connues de l’étoile, trouvées sur VizieR à partir des bibcodes de leurs solutions dans l’archive. Une vitesse publiée à moins d’une minute d’une pose déjà dans la série est le même spectre, et elle est écartée.',
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
  latex: {
    en: 'Write the report in LaTeX and compile it to PDF (pdflatex). Unticked: the text report, the summary (JSON) and every figure are written all the same.',
    fr: 'Écrire le rapport en LaTeX et le compiler en PDF (pdflatex). Décochée : le rapport texte, le résumé (JSON) et toutes les figures sont écrits quand même.',
  },
  mode: {
    en: 'Where the script will run. The page writes the script and never runs it: you start it yourself, with the command it shows. This machine: the paths are local, and the script is downloaded or saved in the working folder. A server (ssh): the files are picked on the server through ssh, and the script is written for it, its paths and its koloa, then downloaded or sent there.',
    fr: 'Où le script tournera. La page écrit le script et ne le lance jamais : vous le lancez vous-même, avec la commande qu’elle montre. Cette machine : les chemins sont locaux, et le script est téléchargé ou enregistré dans le dossier de travail. Un serveur (ssh) : les fichiers sont choisis sur le serveur par ssh, et le script est écrit pour lui, ses chemins et son koloa, puis téléchargé ou envoyé là-bas.',
  },
  host: {
    en: 'The server, as ssh knows it: a name of your ~/.ssh/config (rali, for instance) or user@host. koloa uses your ssh keys and never asks for a password (ssh -o BatchMode=yes): if ssh asks you for one, set up a key first. Test checks the connection and koloa there.',
    fr: 'Le serveur, tel que ssh le connaît : un nom de votre ~/.ssh/config (rali, par exemple) ou utilisateur@hôte. koloa utilise vos clés ssh et ne demande jamais de mot de passe (ssh -o BatchMode=yes) : si ssh vous en demande un, installez d’abord une clé. Tester vérifie la connexion et koloa là-bas.',
  },
  workdir: {
    en: 'The folder of that machine where the script is saved and where it works when you start it: the batch folder is made inside it. An absolute path, /spirou2/batches for instance. Empty: the home folder of the server, or on this machine the folder koloa runs from.',
    fr: 'Le dossier de cette machine où le script est enregistré et où il travaille quand vous le lancez : le dossier du lot y est créé. Un chemin absolu, /spirou2/batches par exemple. Vide : le dossier personnel du serveur, ou sur cette machine le dossier d’où koloa est lancé.',
  },
  koloa_cmd: {
    en: 'How koloa is called on that machine: koloa when it is installed (pip install), or a Python and its module, /path/to/python -m koloa.cli, with PYTHONPATH=/path/to/koloa before it for a copy that is not installed. Empty: koloa. Test shows what answers there.',
    fr: 'Comment koloa est appelé sur cette machine : koloa s’il est installé (pip install), ou un Python et son module, /chemin/python -m koloa.cli, précédé de PYTHONPATH=/chemin/koloa pour une copie non installée. Vide : koloa. Tester montre ce qui répond là-bas.',
  },
  batchdir: {
    en: 'The name of the batch folder, made in the working folder: the report of each target (BATCH/<target>/), the logs (BATCH/logs/), a copy of every report PDF (BATCH/pdf/) and the summary of the batch (batch_summary.pdf, .txt, .csv). Empty: batch_ and today’s date.',
    fr: 'Le nom du dossier du lot, créé dans le dossier de travail : le rapport de chaque cible (LOT/<cible>/), les journaux (LOT/logs/), une copie de chaque PDF de rapport (LOT/pdf/) et le résumé du lot (batch_summary.pdf, .txt, .csv). Vide : batch_ et la date du jour.',
  },
  jobs: {
    en: 'How many targets run at a time. Each report runs the two chains of the FIP and, for the GP check, several processes: on a laptop keep 1; on a server with many cores, about a quarter of them is a fair start (10 on 40 cores). Default 1.',
    fr: 'Combien de cibles tournent à la fois. Chaque rapport fait tourner les deux chaînes du FIP et, pour le test GP, plusieurs processus : sur un portable, gardez 1 ; sur un serveur à nombreux cœurs, environ le quart est un bon début (10 sur 40 cœurs). Par défaut 1.',
  },
  bashrc: {
    en: 'Load the shell’s environment at the start of the script (source ~/.bashrc): a script run by bash does not read it otherwise, and your DACE key (DACE_API_KEY) and your conda paths may live there. Untick it if your .bashrc does something a script should not.',
    fr: 'Charger l’environnement du shell au début du script (source ~/.bashrc) : un script lancé par bash ne le lit pas autrement, et votre clé DACE (DACE_API_KEY) et vos chemins conda peuvent s’y trouver. Décochez-la si votre .bashrc fait quelque chose qu’un script ne devrait pas faire.',
  },
  paste: {
    en: 'Add many targets at once, one per line: the SIMBAD name, then a | and the files of that target separated by spaces, e.g. GJ 436 | /data/lbl_GJ436.rdb /data/lbl2_GJ436.rdb. A line without | is a name alone, its report from the archives only. A file may end with =NAME to give its instrument: /data/x.rdb=NIRPS_LBL2.',
    fr: 'Ajouter plusieurs cibles d’un coup, une par ligne : le nom SIMBAD, puis un | et les fichiers de cette cible séparés par des espaces, p. ex. GJ 436 | /data/lbl_GJ436.rdb /data/lbl2_GJ436.rdb. Une ligne sans | est un nom seul, son rapport tiré des archives seules. Un fichier peut finir par =NOM pour donner son instrument : /data/x.rdb=NIRPS_LBL2.',
  },
  targets: {
    en: 'Each line is a target: its SIMBAD name and 0, 1 or more files, each with its instrument (empty: read from the file). With no file, its report comes from the archives alone; with no name, the name is guessed from the first file. Every target takes the options above. The list is kept in this browser until you clear it.',
    fr: 'Chaque ligne est une cible : son nom SIMBAD et 0, 1 ou plusieurs fichiers, chacun avec son instrument (vide : lu dans le fichier). Sans fichier, son rapport vient des archives seules ; sans nom, le nom est deviné à partir du premier fichier. Chaque cible prend les options ci-dessus. La liste est gardée dans ce navigateur jusqu’à ce que vous la vidiez.',
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
    { key: 'dace', text: 'DACE', checked: true },
    { key: 'carmenes', text: 'CARMENES DR1', checked: true },
    { key: 'tess', text: 'TESS', help: 'tess_opt', checked: true },
    { key: 'vizier', text: 'VizieR', checked: true },
    { key: 'archive', label: 'archive', help: 'archive_opt', checked: true },
    { key: 'gpcheck', label: 'gpcheck', checked: true },
    { key: 'duck', label: 'duck', checked: true },
    { key: 'mcmc', text: 'MCMC', checked: false },
    { key: 'latex', label: 'latex', checked: true },
  ],
};

// the fields and the boxes, for the options of `forName` (detailed, batch)
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
// a file browser in the page: this machine or a server, through ssh
// -----------------------------------------------------------------------------
function browse({ host = '', start = '', kind = 'file' } = {}) {
  return new Promise((resolve) => {
    const shade = document.createElement('div');
    shade.className = 'browser-shade';
    shade.innerHTML = `<div class="browser" role="dialog">
      <div class="browser-head"><b>${esc(host ? `${host}:` : t('this_machine'))}</b> <span class="bpath"></span>
        <span class="spacer"></span>${kind === 'folder' ? `<button type="button" class="small bchoose">${esc(t('choose_folder'))}</button>` : ''}
        <button type="button" class="small bclose">&times;</button></div>
      <div class="blist"><p class="hint"><span class="spin"></span></p></div></div>`;
    document.body.appendChild(shade);
    let here = start;
    const close = (val) => { shade.remove(); resolve(val); };
    const open = async (path) => {
      const list = shade.querySelector('.blist');
      list.innerHTML = '<p class="hint"><span class="spin"></span></p>';
      try {
        const res = await api(`/api/ls?${new URLSearchParams({ host, path })}`);
        here = res.path;
        shade.querySelector('.bpath').textContent = res.path;
        const rows = [`<div class="bentry bdir" data-path="${esc(res.parent)}">..</div>`]
          .concat(res.entries.map((ent) => `<div class="bentry ${ent.dir ? 'bdir' : 'bfile'}${/\.(rdb|csv|dat|txt)$/i.test(ent.name) ? ' bdata' : ''}" data-path="${esc(res.path.replace(/\/$/, '') + '/' + ent.name)}">${esc(ent.name)}${ent.dir ? '/' : ''}</div>`));
        list.innerHTML = rows.join('');
      } catch (err) {
        list.innerHTML = `<p class="hint bad">${esc(err.message)}</p>`;
      }
    };
    shade.addEventListener('click', (e) => {
      if (e.target === shade || e.target.closest('.bclose')) close(null);
      if (e.target.closest('.bchoose')) close(here);
      const ent = e.target.closest('.bentry');
      if (!ent) return;
      if (ent.classList.contains('bdir')) open(ent.dataset.path);
      else if (kind === 'file') close(ent.dataset.path);
    });
    open(start);
  });
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

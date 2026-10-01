// koloa's GUI: the fields, the command lines, the runs (koloa.gui serves it)
'use strict';

const TEXT = {
  en: {
    tagline: 'Outlier-aware radial velocities · on this machine', docs: 'Docs',
    star: 'Star', star_hint: 'Its SIMBAD name: the archives, DACE, CARMENES and TESS find it by that name.',
    resolve: 'Resolve', velocities: 'Velocities',
    velocities_hint: 'A file of velocities (LBL .rdb, csv, DACE csv), the archives gathered for the star, or both. Untick an instrument (or click it in the legend) to leave it out of the report.',
    use: 'Used', exclude: 'Instruments left out',
    file: 'File (optional)', root: 'Archives folder', plot: 'Plot',
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
    use: 'Utilisé', exclude: 'Instruments écartés',
    file: 'Fichier (facultatif)', root: 'Dossier des archives', plot: 'Tracer',
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

function applyLang() {
  document.documentElement.lang = lang;
  document.querySelectorAll('[data-i18n]').forEach((el) => { el.textContent = t(el.dataset.i18n); });
  $('lang').textContent = lang === 'fr' ? 'EN' : 'FR';
  showCwd();
  renderJobs();
  checkArchives();
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
// the fields and the command lines
// -----------------------------------------------------------------------------
function folderName(text) {
  return text.trim().replace(/[^A-Za-z0-9+\-.]+/g, '_').replace(/^_+|_+$/g, '');
}

function defaultOutdir() {
  const target = $('target').value.trim();
  const file = $('file').value.trim();
  const stem = target || (file ? file.split('/').pop().replace(/\.[^.]*$/, '') : '');
  return stem ? `reports/${folderName(stem)}` : 'reports/star';
}

function options(action) {
  const opts = { target: $('target').value, file: $('file').value, root: $('root').value,
                 outdir: $('outdir').value.trim() || defaultOutdir() };
  document.querySelectorAll(`[data-for="${action}"]`).forEach((el) => {
    opts[el.dataset.opt] = el.type === 'checkbox' ? el.checked : el.value;
  });
  if (action === 'gather') opts.refresh = !!(onDisk && onDisk.exists);
  return opts;
}

let cmdTimer = null;
function updateCommands() {
  clearTimeout(cmdTimer);
  cmdTimer = setTimeout(updateCommandsNow, 150);
}

async function updateCommandsNow() {
  {
    $('outdir').placeholder = defaultOutdir();
    for (const action of ['gather', 'detailed']) {
      const box = $(`cmd-${action}`);
      try {
        const res = await api('/api/command', { action, options: options(action) });
        box.textContent = res.line;
        box.classList.remove('error');
        $(`run-${action}`).disabled = false;
      } catch (err) {
        box.textContent = action === 'detailed' ? t('need') : err.message;
        box.classList.add('error');
        $(`run-${action}`).disabled = true;
      }
    }
    if (onDisk && onDisk.busy) $('run-gather').disabled = true;
  }
}

let onDisk = null;
let diskTimer = null;
function checkArchives() {
  clearTimeout(diskTimer);
  diskTimer = setTimeout(async () => {
    const target = $('target').value.trim();
    const btn = $('run-gather');
    if (!target) { onDisk = null; $('ondisk').textContent = ''; btn.textContent = t('run_gather'); return; }
    try {
      onDisk = await api(`/api/archives?${new URLSearchParams({ target, root: $('root').value.trim() })}`);
    } catch (err) { onDisk = null; }
    if (onDisk && onDisk.busy) {
      btn.textContent = t('busy'); btn.disabled = true;
      $('ondisk').textContent = '';
    } else if (onDisk && onDisk.exists) {
      btn.textContent = t('refresh');
      const what = Object.entries(onDisk.archives || {}).map(([k, v]) => `${k} ${v}`).join(', ');
      $('ondisk').textContent = `${t('on_disk')} (${onDisk.created || ''}): ${what}`;
    } else {
      btn.textContent = t('run_gather');
      $('ondisk').textContent = '';
    }
    updateCommandsNow();
  }, 250);
}

let cwd = '';
let archiveDate = null;
function showCwd() {
  if (cwd) $('cwd').textContent = `${t('cwd')} ${cwd} ${t('rel')}`;
  $('archive-date').textContent = archiveDate ? `${t('archive_kept')} ${archiveDate}` : t('archive_none');
}

async function refreshInfo() {
  try {
    const info = await api('/api/info');
    cwd = info.cwd; archiveDate = info.archive; showCwd();
  } catch (err) { /* later */ }
}

// -----------------------------------------------------------------------------
// the star
// -----------------------------------------------------------------------------
function tile(label, value, wide) {
  return `<div class="stat${wide ? ' wide' : ''}"><div class="label">${esc(label)}</div><div class="value">${value}</div></div>`;
}

async function resolveStar() {
  const name = $('target').value.trim();
  if (!name) return;
  const box = $('ident');
  box.innerHTML = `<p class="hint"><span class="spin"></span> ${esc(t('resolving'))}</p>`;
  try {
    const id = await api(`/api/resolve?name=${encodeURIComponent(name)}`);
    const pos = id.ra != null ? `${id.ra.toFixed(5)}, ${id.dec.toFixed(5)}` : '-';
    const others = [id.gj, id.hd, id.hip].filter(Boolean).map(esc).join(', ') || '-';
    const chip = (per, text) => `<button type="button" class="chip" data-prot="${esc(per)}">${text}</button>`;
    const rot = (id.variability || []).map((v) => chip(v.period, `${esc(v.type)} ${v.period} d <span class="hint">(${esc(v.bibcode || '')})</span>`));
    if (id.carmenes && id.carmenes.p_rot) rot.push(chip(id.carmenes.p_rot, `ROT ${esc(id.carmenes.p_rot)} d <span class="hint">(CARMENES, ${esc(id.carmenes.p_rot_source || '')})</span>`));
    if (rot.length) rot.push(`<span class="hint">${esc(t('pick'))}</span>`);
    const carm = id.carmenes ? `${esc(id.carmenes.carmenes_id)}, ${esc(id.carmenes.nobs)} ${esc(t('points'))}` : esc(t('not_in'));
    box.innerHTML = '<div class="stats-grid">'
      + tile(t('main'), esc(id.main)) + tile(t('tic'), esc((id.tic || '-').replace('TIC ', '')))
      + tile(t('pos'), pos) + tile(t('gaia'), esc((id.gaia_dr3 || '-').replace('Gaia DR3 ', '')))
      + tile(t('names'), others, true) + tile(t('rotation'), rot.join('<br>') || esc(t('none')), true)
      + tile(t('carmenes'), carm, true) + '</div>';
  } catch (err) {
    box.innerHTML = `<p class="hint bad">${esc(err.message)}</p>`;
  }
}

// -----------------------------------------------------------------------------
// the velocities
// -----------------------------------------------------------------------------
async function plotVelocities() {
  const note = $('rvnote');
  note.innerHTML = `<span class="spin"></span> ${esc(t('loading'))}`;
  const q = new URLSearchParams({ file: $('file').value.trim(), target: $('target').value.trim(), root: $('root').value.trim() });
  try {
    const res = await api(`/api/rv?${q}`);
    note.textContent = (res.notes || []).join(' · ');
    const div = $('rvplot');
    if (!res.instruments.length) {
      div.classList.remove('on'); if (window.Plotly) Plotly.purge(div);
      $('rvtable').innerHTML = `<p class="hint">${esc(t('no_rv'))}</p>`;
      return;
    }
    const traces = res.instruments.map((inst, i) => ({
      x: inst.time, y: inst.rv, name: `${inst.name} (${inst.n})`, type: 'scatter', mode: 'markers',
      error_y: { type: 'data', array: inst.err, visible: true, thickness: 1, width: 0, color: COLOURS[i % 8] },
      marker: { color: COLOURS[i % 8], symbol: SYMBOLS[i % 8], size: 7, line: { color: '#08111f', width: 1 } },
      hovertemplate: `${inst.name}<br>rjd %{x:.4f}<br>%{y:.2f} m/s<extra></extra>`,
    }));
    div.classList.add('on');
    const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#a8b4ca' };
    if (window.Plotly) {
      Plotly.newPlot(div, traces, {
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)',
        font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
        margin: { l: 60, r: 10, t: 10, b: 50 }, legend: { orientation: 'h', y: -0.2 }, showlegend: traces.length > 1,
        xaxis: { ...axis, title: 'BJD - 2400000', tickformat: '.0f', exponentformat: 'none' }, yaxis: { ...axis, title: 'RV - median [m/s]' },
      }, { responsive: true, displaylogo: false });
    }
    lastRV = res.instruments;
    $('rvtable').innerHTML = `<table class="mini"><tr><th>${esc(t('use'))}</th><th>${esc(t('inst'))}</th><th>${esc(t('n'))}</th><th>${esc(t('rms'))}</th></tr>`
      + res.instruments.map((inst, i) => `<tr data-row="${esc(inst.name)}"><td><input type="checkbox" data-inst="${esc(inst.name)}" checked></td>`
        + `<td><span class="swatch" style="background:${COLOURS[i % 8]}"></span>${esc(inst.name)}</td>`
        + `<td class="num">${inst.n}</td><td class="num">${inst.rms.toFixed(2)}</td></tr>`).join('') + '</table>';
    if (window.Plotly && div.on) {
      div.removeAllListeners && div.removeAllListeners('plotly_legendclick');
      div.on('plotly_legendclick', (ev) => {
        const name = lastRV[ev.curveNumber].name;
        setExcluded(name, !excludedSet().has(name.toUpperCase()));
        return false;
      });
    }
    styleExcluded();
  } catch (err) {
    note.innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
}

// the instruments left out: the field of the report is what counts; the
// boxes of the table and the legend of the plot write into it
let lastRV = null;
function excludedSet() {
  return new Set($('exclude').value.split(/[\s,]+/).filter(Boolean).map((s) => s.toUpperCase()));
}

function setExcluded(name, off) {
  const kept = $('exclude').value.split(/[\s,]+/).filter((s) => s && s.toUpperCase() !== name.toUpperCase());
  if (off) kept.push(name);
  $('exclude').value = kept.join(' ');
  updateCommands();
  styleExcluded();
}

function styleExcluded() {
  if (!lastRV) return;
  const ex = excludedSet();
  const off = lastRV.map((inst) => ex.has(inst.name.toUpperCase()));
  if (window.Plotly && $('rvplot').data) Plotly.restyle('rvplot', { opacity: off.map((o) => (o ? 0.12 : 1)) });
  lastRV.forEach((inst, i) => {
    const row = document.querySelector(`tr[data-row="${CSS.escape(inst.name)}"]`);
    if (!row) return;
    row.classList.toggle('off', off[i]);
    row.querySelector('input').checked = !off[i];
  });
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
      if (job.action === 'gather' && jobs.get(job.id).status !== 'running') checkArchives();
      if (job.action === 'archive' && jobs.get(job.id).status !== 'running') refreshInfo();
    } catch (err) { /* the server may be gone */ }
  }
  renderJobs();
}

async function run(action) {
  try {
    keep(await api('/api/run', { action, options: options(action) }));
    renderJobs();
    if (action === 'gather') checkArchives();
    $('jobs').scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (err) {
    alert(err.message);
  }
}

// -----------------------------------------------------------------------------
// start
// -----------------------------------------------------------------------------
document.addEventListener('input', (e) => {
  updateCommands();
  if (e.target.id === 'target' || e.target.id === 'root') checkArchives();
  if (e.target.id === 'exclude') styleExcluded();
});
document.addEventListener('change', (e) => {
  if (e.target.dataset && e.target.dataset.inst) setExcluded(e.target.dataset.inst, !e.target.checked);
});
document.addEventListener('change', updateCommands);
$('resolve').addEventListener('click', resolveStar);
$('target').addEventListener('keydown', (e) => { if (e.key === 'Enter') resolveStar(); });
$('plot').addEventListener('click', plotVelocities);
$('run-gather').addEventListener('click', () => run('gather'));
$('run-detailed').addEventListener('click', () => run('detailed'));
$('run-archive').addEventListener('click', () => run('archive'));
$('lang').addEventListener('click', () => {
  lang = lang === 'fr' ? 'en' : 'fr';
  try { localStorage.setItem('koloa-lang', lang); } catch (e) { /* no storage */ }
  applyLang();
});
document.addEventListener('click', async (e) => {
  const copy = e.target.closest('[data-copy]');
  if (copy) {
    const text = $(copy.dataset.copy).textContent;
    try { await navigator.clipboard.writeText(text); } catch (err) { /* old browsers */ }
    copy.textContent = t('copied');
    setTimeout(() => { copy.textContent = t('copy'); }, 1200);
  }
  const prot = e.target.closest('[data-prot]');
  if (prot) {
    $('rotation').value = prot.dataset.prot;
    $('fipgp').value = 'sho';
    updateCommands();
  }
  const stop = e.target.closest('[data-stop]');
  if (stop) { await api('/api/stop', { id: stop.dataset.stop }); poll(); }
});

(async () => {
  applyLang();
  await refreshInfo();
  try { (await api('/api/jobs')).forEach(keep); } catch (err) { /* nothing yet */ }
  renderJobs();
  // the fields from the address (?target=GJ%20436&file=...&root=...&plot=1)
  const params = new URLSearchParams(location.search);
  for (const key of ['target', 'file', 'root', 'outdir', 'exclude']) {
    if (params.get(key)) $(key).value = params.get(key);
  }
  updateCommands();
  checkArchives();
  if (params.get('target')) resolveStar();
  if (params.get('plot')) plotVelocities();
  setInterval(poll, 1500);
})();

// koloa's GUI: the fields, the command lines, the runs (koloa.gui serves it)
'use strict';

// the shared parts (TEXT, t, $, esc, api, the runs) are in common.js

// -----------------------------------------------------------------------------
// the fields and the command lines
// -----------------------------------------------------------------------------
// the files: one row each (its path, and its instrument when the file
//   does not say), as many as wanted
let fileRows = [{ path: '', label: '' }];
function renderFiles() {
  $('filelist').innerHTML = fileRows.map((row, i) => `<div class="filerow">
      <input type="text" class="fpath" data-row="${i}" value="${esc(row.path)}" autocomplete="off">
      <button type="button" data-pick="file" data-row="${i}">${esc(t('browse'))}</button>
      <input type="text" class="flabel" data-row="${i}" value="${esc(row.label)}" autocomplete="off">
      <button type="button" class="small" data-del="${i}" title="${esc(t('remove'))}">&times;</button></div>`).join('');
}

function filesNow() {
  return fileRows.filter((row) => row.path.trim()).map((row) => ({ path: row.path.trim(), label: row.label.trim() }));
}

function folderName(text) {
  return text.trim().replace(/[^A-Za-z0-9+\-.]+/g, '_').replace(/^_+|_+$/g, '');
}

function defaultOutdir() {
  const target = $('target').value.trim();
  const file = (filesNow()[0] || {}).path || '';
  const stem = target || (file ? file.split('/').pop().replace(/\.[^.]*$/, '') : '');
  return stem ? `reports/${folderName(stem)}` : 'reports/star';
}

function options(action) {
  const opts = { target: $('target').value, files: filesNow(), root: $('root').value,
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
    for (const action of ['gather', 'detailed']) {
      const box = $(`cmd-${action}`);
      try {
        const res = await api('/api/command', { action, options: options(action) });
        box.textContent = res.line;
        box.classList.remove('error');
        $(`run-${action}`).disabled = false;
      } catch (err) {
        box.textContent = action !== 'detailed' ? err.message : (err.message.includes('tick') ? t('need_archive') : t('need'));
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
      showDiskPoints();
      const what = Object.entries(onDisk.archives || {}).map(([k, v]) => `${k} ${v}`).join(', ');
      $('ondisk').textContent = `${t('on_disk')} (${onDisk.created || ''}): ${what}`;
    } else {
      btn.textContent = t('run_gather');
      $('ondisk').textContent = '';
    }
    showDiskPoints();
    updateCommandsNow();
  }, 250);
}

// the archives gathered, by the velocities: the same choice as the report's
//   DACE and CARMENES DR1 boxes (what is plotted is what the report uses)
const ARCHIVES = ['dace', 'carmenes', 'vizier'];
// the archives of the manifest of a gather, by box
const ARCHIVE_KEY = { dace: 'dace', carmenes: 'carmenes', vizier: 'published' };
function reportBox(key) { return document.querySelector(`input[data-for="detailed"][data-opt="${key}"]`); }
function syncMirrors() {
  document.querySelectorAll('[data-mirror]').forEach((box) => { const twin = reportBox(box.dataset.mirror); if (twin) box.checked = twin.checked; });
}
function showDiskPoints() {
  const pts = (onDisk && onDisk.exists && onDisk.points) || {};
  for (const key of ARCHIVES) {
    const num = pts[ARCHIVE_KEY[key]];
    $(`pts-${key}`).textContent = num ? `(${num} ${t('points_word')})` : onDisk && onDisk.exists ? '(0)' : '';
  }
  if (!(onDisk && onDisk.exists)) $('pts-dace').textContent = $('target').value.trim() ? `(${t('nothing_gathered')})` : '';
}
async function diskState() {
  const target = $('target').value.trim();
  if (!target) return null;
  try { return await api(`/api/archives?${new URLSearchParams({ target, root: $('root').value.trim() })}`); } catch (err) { return null; }
}
// no file and no archive ticked: the archives gathered are all there is
function archivesByDefault(disk) {
  if (filesNow().length || !disk || !disk.exists) return false;
  if (ARCHIVES.some((key) => reportBox(key) && reportBox(key).checked)) return false;
  let ticked = false;
  for (const key of ARCHIVES) {
    if ((disk.points || {})[ARCHIVE_KEY[key]] && reportBox(key)) { reportBox(key).checked = true; ticked = true; }
  }
  syncMirrors();
  updateCommands();
  return ticked;
}

let cwd = '';
let archiveDate = null;
let euDate = null;          // when exoplanet.eu's catalogue was fetched
let euTarget = null, euTries = 0;
function showCwd() {
  if (cwd) $('cwd').textContent = `${t('cwd')} ${cwd} ${t('rel')}`;
  $('archive-date').textContent = (archiveDate ? `${t('archive_kept')} ${archiveDate}` : t('archive_none'))
    + (euDate ? ` \u00b7 exoplanet.eu: ${euDate}` : '');
}

async function refreshInfo() {
  try {
    const info = await api('/api/info');
    cwd = info.cwd_shown || info.cwd; archiveDate = info.archive; euDate = info.eu; showCwd();
  } catch (err) { /* later */ }
}

// -----------------------------------------------------------------------------
// the star
// -----------------------------------------------------------------------------
// the mass of the star: the resolver's (the archive's, or rough from its
//   spectral type), unless one was typed
function setMass(star) {
  if (!star || !star.mass || $('mstar').dataset.typed) return;
  $('mstar').value = (+star.mass).toFixed(3);
  $('mstar_err').value = (+(star.mass_err || 0)).toFixed(3);
  // the source in the page's language when it is one of koloa's own
  const src = { 'its spectral type (Pecaut & Mamajek 2013)': t('src_spt'), 'NASA Exoplanet Archive': t('src_archive') }[star.source] || star.source;
  $('mstarsrc').textContent = star.source ? `${t('from_word')} ${src}` : '';
}

// m sin i of a companion [Earth masses], exactly (m not << M), and its
//   error from draws of K, P, e and M (16th to 84th percentiles)
const GM_SUN = 1.32712440018e20, MEARTH_MSUN = 3.003489e-6, MNEP_MEARTH = 17.14775, MJUP_MEARTH = 317.8284;
function msini(K, P, e, M) {
  const f = P * 86400 * Math.abs(K) ** 3 * (1 - e * e) ** 1.5 / (2 * Math.PI * GM_SUN);
  let m = Math.cbrt(f * M * M);
  for (let k = 0; k < 40; k++) m = Math.cbrt(f * (M + m) ** 2);
  return m / MEARTH_MSUN;
}
function gauss() {
  let u = 0, v = 0;
  while (!u) u = Math.random();
  while (!v) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
function minimumMass(f) {
  const M = +$('mstar').value, Merr = +$('mstar_err').value || 0;
  if (!(M > 0)) return null;
  const e = f.kind === 'kepler' ? f.e || 0 : 0, eerr = f.kind === 'kepler' ? f.e_err || 0 : 0;
  const best = msini(f.K, f.period, e, M);
  const vals = [];
  for (let k = 0; k < 3000; k++) {
    vals.push(msini(f.K + (f.K_err || 0) * gauss(), f.period + (f.P_err || 0) * gauss(),
      Math.min(0.99, Math.abs(e + eerr * gauss())), Math.max(0.01, M + Merr * gauss())));
  }
  vals.sort((a, b) => a - b);
  return { best, lo: best - percentile(vals, 15.87), hi: percentile(vals, 84.13) - best, M, Merr };
}
// the DACE copies of the instruments of the files already seen (left out
//   once, by default; not when a result is recalled as it was)
const daceSeen = new Set();
let daceAuto = true;
// the star of the page (its radius and effective temperature, for the
//   equilibrium temperature of a planet)
let starNow = {};
// the equilibrium temperature of a planet: Teff sqrt(R* / 2a) (1 - A)^(1/4),
//   the heat spread over the planet, a from the period and the mass of the
//   star of the card; for a Bond albedo of 0 and of 0.3
function teqText(P, star) {
  const M = +$('mstar').value || (star && star.mass);
  const R = star && star.radius, T = star && star.teff;
  if (!(M > 0 && R > 0 && T > 0 && P > 0)) return '';
  const a = Math.cbrt(GM_SUN * M * (P * 86400) ** 2 / (4 * Math.PI ** 2));
  const teq = T * Math.sqrt(R * 6.957e8 / (2 * a));
  return ` <span class="massnote">T<sub>${esc(t('teq_sub'))}</sub> \u2248 <b>${teq.toFixed(0)} K</b> (A = 0) \u00b7 ${(teq * 0.7 ** 0.25).toFixed(0)} K (A = 0.3)`
    + ` <span class="hint">(a = ${(a / 1.495978707e11).toFixed(3)} au, T\u2605 = ${T.toFixed(0)} K, R\u2605 = ${R.toFixed(3)} R\u2609)</span></span>`;
}
// the BIC of a Keplerian orbit against no planet and against a circular
//   orbit (theirs minus its own: positive, the orbit is the better model);
//   in words, the scale of Kass & Raftery (1995): 2 to 6 positive, 6 to 10
//   strong, above 10 very strong
function bicText(f) {
  const b = f.bic;
  if (!b || !Number.isFinite(b.d_none)) return '';
  const fmt = (v) => `${v >= 0 ? '+' : '\u2212'}${Math.abs(v).toFixed(1)}`;
  const word = (v) => t(v > 10 ? 'bic_very' : v > 6 ? 'bic_strong' : v > 2 ? 'bic_positive' : v >= -2 ? 'bic_none' : 'bic_against');
  return `<span class="massnote">\u0394BIC = <b>${fmt(b.d_none)}</b> ${esc(t('bic_vs_none'))}`
    + (Number.isFinite(b.d_circular) ? ` \u00b7 <b>${fmt(b.d_circular)}</b> ${esc(t('bic_vs_circ'))} (${esc(word(b.d_circular))})` : '')
    + ` <span class="hint">${esc(t('bic_sign'))}</span></span>`;
}
// an asymmetric error, the + over the - after the value
const pmStack = (lo, hi) => `<span class="pm"><span>+${hi}</span><span>\u2212${lo}</span></span>`;
// an error after its value: one value after a +- when its two sides are
//   within 5 % of each other, else the + over the -
const pmErr = (lo, hi, fmt) => (Math.abs(lo - hi) <= 0.05 * Math.max(lo, hi)
  ? ` \u00b1 ${fmt(0.5 * (lo + hi))}` : pmStack(fmt(lo), fmt(hi)));
function massText(f) {
  const mm = minimumMass(f);
  if (!mm) return `<span class="massnote hint">${esc(t('no_mstar'))}</span>`;
  const fmt = (v) => (v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2) : v.toPrecision(2));
  const units = [['M\u2295', 1], ['M\u2646', MNEP_MEARTH], ['M\u2643', MJUP_MEARTH]];
  const what = f.transit ? t('mass_transit') : 'm sin i';
  return `<span class="massnote">${esc(what)} = ` + units.map(([u, s]) => `<b>${fmt(mm.best / s)}</b>${pmErr(mm.lo / s, mm.hi / s, fmt)} ${u}`).join(' = ')
    + ` <span class="hint">(M\u2605 = ${mm.M.toFixed(3)} \u00b1 ${mm.Merr.toFixed(3)} M\u2609)</span></span>`;
}

// a name as a trace's uid (Plotly makes CSS selectors of them: 'HIRES
//   (Teklu+ 2025)' would break one)
const safeId = (name) => String(name).replace(/[^A-Za-z0-9_-]/g, '_');

function tile(label, value, wide) {
  const span = wide === 'full' ? ' full' : wide ? ' wide' : '';
  return `<div class="stat${span}"><div class="label">${esc(label)}</div><div class="value">${value}</div></div>`;
}

// the star of a file, when the page has none: its OBJECT column, else its
//   name (lbl_<OBJECT>_<TEMPLATE>.rdb), in APERO's names, as SIMBAD knows it
async function starFromFile(path) {
  if (!path || !path.trim() || $('target').value.trim()) return;
  let res;
  try { res = await api(`/api/star_of_file?${new URLSearchParams({ path: path.trim() })}`); } catch (err) { return; }
  if (!res.target || $('target').value.trim()) return;
  $('target').value = res.target;
  const from = res.source === 'OBJECT column' ? `OBJECT ${res.raw}` : `${t('file_name')} ${res.raw}`;
  $('targetsrc').textContent = res.apero
    ? `${t('from_file')}: ${from} \u2192 ${res.apero} (APERO), ${res.target} ${t('for_simbad')}`
    : `${t('from_file')}: ${from} (${t('not_apero')})`;
  checkArchives();
  updateCommands();
  resolveStar();
}

let varTarget = '', varTries = 0;
async function resolveStar(refresh) {
  const name = $('target').value.trim();
  if (!name) return;
  const box = $('ident');
  box.innerHTML = `<p class="hint"><span class="spin"></span> ${esc(t('resolving'))}</p>`;
  try {
    const id = await api(`/api/resolve?${new URLSearchParams({ name, root: $('root').value.trim(), refresh: refresh === true ? '1' : '' })}`);
    const pos = id.ra != null ? `${id.ra.toFixed(5)}, ${id.dec.toFixed(5)}` : '-';
    const others = [id.gj, id.hd, id.hip].filter(Boolean).map(esc).join(', ') || '-';
    // each period of the literature: a tick that makes it the SHO's
    const now = $('rotation').value.trim();
    const chip = (per, text) => `<label class="prot"><input type="radio" name="prot" value="${esc(per)}"${now && +now === +per ? ' checked' : ''}> `
      + `${text}</label>`;
    // each period once (SIMBAD may list one twice)
    const seen = new Set();
    const rot = (id.variability || []).filter((v) => { const k = `${v.type} ${v.period} ${v.bibcode}`; const n = !seen.has(k); seen.add(k); return n; })
      .map((v) => chip(v.period, `${esc(v.type)} ${v.period} d <span class="hint">(${esc(v.bibcode || '')})</span>`));
    if (id.carmenes && id.carmenes.p_rot) rot.push(chip(id.carmenes.p_rot, `ROT ${esc(id.carmenes.p_rot)} d <span class="hint">(CARMENES, ${esc(id.carmenes.p_rot_source || '')})</span>`));
    if (id.archive_rotation) rot.push(chip(id.archive_rotation, `ROT ${esc(id.archive_rotation)} d <span class="hint">(NASA Exoplanet Archive)</span>`));
    const num = (v, d) => (v == null || v === '' ? '?' : (+v).toFixed(d));
    const planets = (id.planets || []).map((pl) => `<div class="planet"><b>${esc(pl.name)}</b> `
      + `<span class="pnum">P ${num(pl.P, 4)} d \u00b7 K ${num(pl.K, 2)} m/s`
      + (pl.mass_earth != null ? ` \u00b7 m sin i ${num(pl.mass_earth, 1)} M\u2295` : '') + `</span><br>`
      + `<span class="hint">${pl.reference_url ? `<a href="${esc(pl.reference_url)}" target="_blank">${esc(pl.reference || '')}</a>` : esc(pl.reference || '')}`
      + (pl.source === 'exoplanet.eu' ? ` \u00b7 ${esc(t('eu_only'))}` : ` \u00b7 ${pl.solutions} ${esc(t(pl.solutions === 1 ? 'solution' : 'solutions'))}`)
      + (pl.discovery ? ` \u00b7 ${esc(pl.discovery)}${pl.year ? ` ${esc(pl.year)}` : ''}` : '') + `</span>`
      // exoplanet.eu's answer, beside the archive's
      + (pl.eu && pl.source !== 'exoplanet.eu' ? `<br><span class="hint">exoplanet.eu: <span class="pnum">P ${num(pl.eu.P, 4)} d \u00b7 K ${num(pl.eu.K, 2)} m/s`
        + (pl.eu.e != null ? ` \u00b7 e ${num(pl.eu.e, 2)}` : '') + `</span></span>` : '')
      + `</div>`).join('')
      || `<span class="hint">${esc(id.planets_error || t('none_known'))}</span>`;
    // exoplanet.eu's catalogue on its way (its server is slow): the star
    //   read again when it is there
    const euNote = id.eu && id.eu.pending ? `<span class="hint"><span class="spin"></span> ${esc(t('eu_pending'))}</span>` : '';
    if (id.eu && id.eu.pending) {
      const asked = name;
      euTries = (euTarget === asked ? euTries : 0) + 1;
      euTarget = asked;
      if (euTries <= 30) setTimeout(() => { if ($('target').value.trim() === asked) resolveStar(); }, 30000);
    }
    const tois = (id.tois || []).map((ti) => `<button type="button" class="chip" data-toi="${esc(ti.toi)}"><b>TOI-${esc(ti.toi)}</b> ${esc(ti.disposition || '')} `
      + `<span class="pnum">P ${num(ti.P, 6)} d \u00b7 ${esc(t('transit'))} ${num(ti.tc, 4)}`
      + (ti.depth != null ? ` \u00b7 ${num(ti.depth, 0)} ppm` : '') + (ti.radius != null ? ` \u00b7 ${num(ti.radius, 2)} R\u2295` : '')
      + `</span></button>`).join('');
    const toiTile = tois ? tois + `<span class="hint">${esc(t('toi_pick'))}</span>`
      : `<span class="hint">${esc(id.tois_error || t('toi_none'))}</span>`;
    const disk = { 'the archives folder': t('disk_archives'), 'the copy kept': t('disk_kept') }[id.disk] || id.disk;
  const where = id.disk ? `${t('from_disk')} (${disk}, ${id.disk_date})` : t('asked_now');
    if (rot.length) rot.push(`<label class="prot"><input type="radio" name="prot" value=""${now ? '' : ' checked'}> ${esc(t('no_sho'))}</label>`);
    if (id.variability_pending) {
      rot.push(`<span class="hint"><span class="spin"></span> ${esc(t('var_pending'))}</span>`);
      // read the star again from the disk until SIMBAD's periods are there
      const asked = name;
      varTries = (varTarget === asked ? varTries : 0) + 1;
      varTarget = asked;
      if (varTries <= 40) setTimeout(() => { if ($('target').value.trim() === asked) resolveStar(); }, 5000);
    }
    const carm = id.carmenes ? `${esc(id.carmenes.carmenes_id)}, ${esc(id.carmenes.nobs)} ${esc(t('points'))}` : esc(t('not_in'));
    setMass(id.star);
    starNow = id.star || {};
    const star = id.star || {};
    const typeTile = `${esc(star.sptype || id.sptype || '-')}` + (star.mass ? ` \u00b7 ${star.mass.toFixed(2)} \u00b1 ${(star.mass_err || 0).toFixed(2)} M\u2609` : '');
    box.innerHTML = '<div class="stats-grid">'
      + tile(t('main'), esc(id.main)) + tile(t('tic'), esc((id.tic || '-').replace('TIC ', '')))
      + tile(t('sptype'), typeTile)
      + tile(t('pos'), pos) + tile(t('gaia'), esc((id.gaia_dr3 || '-').replace('Gaia DR3 ', '')))
      + tile(t('names'), others, true)
      + tile(rot.length > 1 ? `${t('rotation')} \u00b7 ${t('use_sho_head')}` : t('rotation'), rot.length ? `<div class="prots">${rot.join('')}</div>` : esc(t('none')), 'full')
      + tile(t('known'), planets + euNote, true)
      + tile(t('toi_title'), toiTile, (id.tois || []).length > 0)
      + tile(t('carmenes'), carm) + '</div>'
      + `<p class="hint">${esc(where)} <button type="button" class="small" id="refresh-star">${esc(t('refresh_star'))}</button></p>`;
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
  const auto = !filesNow().length && archivesByDefault(await diskState());
  // the archives the report will use, and only those
  const asked = readOptions('detailed');
  const q = new URLSearchParams({ files: JSON.stringify(filesNow()), target: $('target').value.trim(), root: $('root').value.trim(),
    dace: asked.dace ? '1' : '', carmenes: asked.carmenes ? '1' : '', vizier: asked.vizier ? '1' : '' });
  try {
    const res = await api(`/api/rv?${q}`);
    // the quick FIP on its own for one instrument; for several, once those
    //   to leave out are unticked (it can take minutes)
    if (await drawVelocities(res, auto ? [t('arch_auto')] : [])) {
      if (res.instruments.length === 1) startQuick(); else quickPrompt();
    }
  } catch (err) {
    note.innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
}

// the velocities by instrument (from /api/rv, or a result recalled)
async function drawVelocities(res, extra) {
  const note = $('rvnote');
  note.textContent = (extra || []).concat(res.notes || []).join(' · ');
  const div = $('rvplot');
  if (!res.instruments.length) {
    div.classList.remove('on'); if (window.Plotly) Plotly.purge(div);
    $('rvtable').innerHTML = `<p class="hint">${esc(t('no_rv'))}</p>`;
    return false;
  }
  lastRes = { res, extra };
  // each instrument about the solution shown (the same for all: the trend
  //   and the signals, no offsets), the median of its velocities minus it
  //   zero; or about its median
  const curve = seriesCurve();
  const aligned = !!curve && (zeroMode === 'fit' || $('foldoverlay').checked);
  res = { ...res, instruments: res.instruments.map((inst) => {
    if (!aligned) return { ...inst, zero: inst.median };
    const z = percentile(inst.time.map((tt, k) => inst.rv[k] - curve.at(tt)).sort((a, b) => a - b), 50);
    return { ...inst, zero: (inst.median || 0) + z, rv: inst.rv.map((v) => Math.round((v - z) * 1000) / 1000) };
  }) };
  if (aligned) note.textContent += ` \u00b7 ${t('zero_fit_note')}: ${curve.name}`;
  seriesAligned = aligned;
  document.querySelectorAll('[data-zero]').forEach((b) => b.classList.toggle('on', b.dataset.zero === zeroMode));
  // the points by instrument, or by their BERV (on one scale for all; an
  //   instrument without it grey)
  const hasBerv = res.instruments.some((inst) => inst.berv);
  $('scol-berv').disabled = !hasBerv;
  const bervMode = seriesColour === 'berv' && hasBerv;
  document.querySelectorAll('[data-scol]').forEach((b) => b.classList.toggle('on', b.dataset.scol === (bervMode ? 'berv' : 'inst')));
  const top = bervMode ? Math.max(1e-3, ...res.instruments.flatMap((inst) => (inst.berv || []).filter((v) => v !== null).map(Math.abs))) : 1;
  rvRight = bervMode ? 100 : PLOT_MARGIN.r;
  let barShown = false;
  const traces = res.instruments.map((inst, i) => {
    let marker = { color: COLOURS[i % 8], symbol: SYMBOLS[i % 8], size: 7, line: { color: '#08111f', width: 1 } };
    if (bervMode && inst.berv) {
      marker = { ...marker, color: inst.berv, colorscale: 'RdBu', reversescale: false, cmin: -top, cmax: top, showscale: !barShown,
        colorbar: { title: { text: 'BERV [km/s]' }, thickness: 12, len: 0.9, x: 1.01, outlinewidth: 0, tickfont: { size: 10 } } };
      barShown = true;
    } else if (bervMode) marker = { ...marker, color: '#5a6476' };
    return { x: inst.time, y: inst.rv, name: `${inst.name} (${inst.source}, ${inst.n})`, type: 'scatter', mode: 'markers',
      error_y: { type: 'data', array: inst.err, visible: true, thickness: 1, width: 0, color: bervMode ? 'rgba(200,220,255,0.35)' : COLOURS[i % 8] },
      marker,
      // the calendar date (and the BERV) of each point, under the cursor
      customdata: inst.time.map((tt, k) => [dateText(tt), inst.berv ? inst.berv[k] : null]),
      hovertemplate: `${inst.name}<br>%{customdata[0]}<br>rjd %{x:.4f}<br>%{y:.2f} m/s`
        + `${inst.berv ? '<br>BERV %{customdata[1]:.2f} km/s' : ''}<extra></extra>` };
  });
  const ninst = traces.length;
  // nothing to see on the axis of the dates: Plotly draws an axis a trace uses
  const times = res.instruments.flatMap((inst) => inst.time);
  traces.push({ x: [Math.min(...times), Math.max(...times)], y: [null, null], xaxis: 'x2', type: 'scatter', mode: 'markers',
    showlegend: false, hoverinfo: 'skip', marker: { opacity: 0 } });
  traces.push(offscaleTrace());
  // under the series, the residuals to the solution, each instrument in
  //   its colour
  let residRange = null;
  if (aligned) {
    // their range: twice the 3 to 97 percentile range, an outlier beyond
    const all = res.instruments.flatMap((inst) => inst.time.map((tt, k) => inst.rv[k] - curve.at(tt))).sort((a, b) => a - b);
    if (all.length > 3) {
      const p3 = percentile(all, 3), p97 = percentile(all, 97), mid = (p3 + p97) / 2, half = Math.max(p97 - p3, 1);
      residRange = [mid - half, mid + half];
    }
    res.instruments.forEach((inst, i) => {
      traces.push({ x: inst.time, y: inst.time.map((tt, k) => Math.round((inst.rv[k] - curve.at(tt)) * 1000) / 1000),
        xaxis: 'x', yaxis: 'y3', type: 'scatter', mode: 'markers', uid: `resid-${i}`, showlegend: false,
        error_y: { type: 'data', array: inst.err, visible: true, thickness: 1, width: 0, color: COLOURS[i % 8] },
        marker: { color: COLOURS[i % 8], symbol: SYMBOLS[i % 8], size: 5, line: { color: '#08111f', width: 0.5 } },
        hovertemplate: `${inst.name}<br>rjd %{x:.4f}<br>${t('resid')} %{y:.2f} m/s<extra></extra>` });
    });
  }
  $('plotcard').classList.add('on');
  div.classList.add('on');
  div.classList.toggle('withres', aligned);
  const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
  lastRV = res.instruments;
  if (window.Plotly) {
    await Plotly.newPlot(div, traces, {
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)',
      font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
      margin: { ...PLOT_MARGIN, r: rvRight }, showlegend: ninst > 1,
      // the legend at the top of the figure, the dates under it
      legend: { orientation: 'h', x: 0, y: 1, yref: 'container', yanchor: 'top' },
      xaxis: { ...axis, title: 'BJD - 2400000', tickformat: '.0f', exponentformat: 'none', ...(aligned ? { anchor: 'y3' } : {}) },
      yaxis: { ...axis, title: aligned ? t('rv_zero_axis') : 'RV - median [m/s]', ...(aligned ? { domain: [RES_TOP, 1] } : {}) },
      ...(aligned ? { yaxis3: { ...axis, domain: [0, RES_TOP - 0.07], anchor: 'x', title: t('resid_axis'), zeroline: true,
        ...(residRange ? { range: residRange } : {}),
        zerolinecolor: 'rgba(200,220,255,0.45)' } } : {}),
      // the calendar dates of the same times, on top
      xaxis2: { ...axis, overlaying: 'x', matches: 'x', side: 'top', showgrid: false, zeroline: false,
        tickmode: 'array', tickvals: [], ticktext: [], ticks: 'outside', ticklen: 4 },
    }, { responsive: true, displaylogo: false });
    div.removeAllListeners && div.removeAllListeners('plotly_relayout');
    // a legend of several rows (many instruments, a solution added): the
    //   dates of the top axis under it, after every redraw
    await roomForLegend(div);
    div.removeAllListeners && div.removeAllListeners('plotly_afterplot');
    div.on('plotly_afterplot', () => setTimeout(() => roomForLegend(div), 0));
    div.on('plotly_relayout', fromZoom);
  }
  view = null;
  $('rvtable').innerHTML = `<table class="mini"><tr><th>${esc(t('use'))}</th><th>${esc(t('inst'))}</th><th>${esc(t('source'))}</th><th>${esc(t('n'))}</th><th>${esc(t('rms'))}</th></tr>`
    + res.instruments.map((inst, i) => `<tr data-row="${esc(inst.name)}"><td><input type="checkbox" data-inst="${esc(inst.name)}" checked></td>`
      + `<td><span class="swatch" style="background:${COLOURS[i % 8]}"></span>${esc(inst.name)}</td>`
      + `<td class="src${(inst.source || '').startsWith('file') ? ' src-file' : ''}">${esc((inst.source || '').startsWith('file') ? inst.source.replace('file', t('src_file')) : inst.source)}</td>`
      + `<td class="num">${inst.n}</td><td class="num">${inst.rms.toFixed(2)}</td></tr>`).join('') + '</table>';
  if (window.Plotly && div.on) {
    div.removeAllListeners && div.removeAllListeners('plotly_legendclick');
    div.on('plotly_legendclick', (ev) => {
      if (ev.curveNumber >= lastRV.length) return true;   // the solution of a fold
      const name = lastRV[ev.curveNumber].name;
      setExcluded(name, !excludedSet().has(name.toUpperCase()));
      return false;
    });
  }
  // the DACE copy of an instrument of the files (NIRPS_DACE: the same
  //   spectra through DACE's pipeline, less precise than LBL's) left out
  //   the first time it is seen, and said: one tick adds it back
  const copies = res.instruments.map((inst) => inst.name).filter((name) => /_DACE$/i.test(name) && !daceSeen.has(name));
  copies.forEach((name) => daceSeen.add(name));
  if (copies.length && daceAuto) {
    copies.forEach((name) => { if (!excludedSet().has(name.toUpperCase())) setExcluded(name, true); });
    $('rvnote').textContent += ` \u00b7 ${copies.join(', ')}: ${t('dace_copy_off')}`;
  }
  styleExcluded();
  showModel();
  return true;
}

// -----------------------------------------------------------------------------
// the ranges of the plot: the sliders along its axes (the y one in asinh,
//   fine around the bulk of the points, coarse toward the outliers), twice
//   the 3 to 97 percentile range, and the instruments kept
// -----------------------------------------------------------------------------
const PLOT_MARGIN = { l: 62, r: 12, t: 60, b: 48 };
// the top margin of the series as tall as its legend, the dates under it
function roomForLegend(div) {
  const full = div._fullLayout;
  if (!full || !full.legend || !full.margin) return null;
  const want = Math.max(PLOT_MARGIN.t, Math.ceil(full.legend._height || 0) + 26);
  return Math.abs(full.margin.t - want) > 2 ? Plotly.relayout(div, { 'margin.t': want }) : null;
}
let rvRight = PLOT_MARGIN.r;   // the right margin of the series (its colour bar)
let lastRes = null;            // the velocities drawn, to draw them again
let seriesColour = 'inst';     // the series by instrument or by BERV
let zeroMode = 'fit';          // each instrument about its fitted offset or its median
let seriesAligned = false;     // the series drawn about the fit of the quick look
const MODEL_CYCLES = 60;       // the fitted model drawn when this few cycles are shown
const RES_TOP = 0.3;           // the residuals under the series: the bottom of the plot
try { zeroMode = localStorage.getItem('koloa-zero') || 'fit'; } catch (err) { /* no storage */ }

// the series drawn again (its zero points, its colours), its ranges kept
async function redrawSeries() {
  if (!lastRes) return;
  const keep = view ? { x: view.x.slice(), y: view.y.slice() } : null;
  const zeroBefore = (lastRV || []).map((inst) => inst.zero);
  await drawVelocities(lastRes.res, lastRes.extra);
  // the y range moved with the zero points: kept only when they did not
  const same = zeroBefore.length === (lastRV || []).length && (lastRV || []).every((inst, i) => inst.zero === zeroBefore[i]);
  if (keep && view) { view.x = keep.x; if (same) view.y = keep.y; applyView(); }
}
try { seriesColour = localStorage.getItem('koloa-seriescolour') || 'inst'; } catch (err) { /* no storage */ }
const STEPS = 1000;
let view = null;   // x and y shown, their domains, the scale of the y slider

function kept() {
  const ex = excludedSet();
  return (lastRV || []).filter((inst) => !ex.has(inst.name.toUpperCase()));
}

function percentile(sorted, q) {
  const pos = (sorted.length - 1) * q / 100;
  const lo = Math.floor(pos), hi = Math.ceil(pos);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
}

// the domains of the sliders: the times and the velocities of the
//   instruments kept (every one when none is)
function domains() {
  const all = kept().length ? kept() : (lastRV || []);
  const times = all.flatMap((inst) => inst.time);
  const vals = kept().flatMap((inst) => inst.rv).sort((a, b) => a - b);
  if (!times.length || !vals.length) return null;
  const tpad = 0.02 * ((Math.max(...times) - Math.min(...times)) || 1);
  const ypad = 0.05 * ((vals[vals.length - 1] - vals[0]) || 1);
  const med = percentile(vals, 50);
  const mad = percentile(vals.map((v) => Math.abs(v - med)).sort((a, b) => a - b), 50);
  // twice the 3-97 percentile range, about its middle
  const p3 = percentile(vals, 3), p97 = percentile(vals, 97);
  const mid = (p3 + p97) / 2, half = Math.max(p97 - p3, 1e-3);
  return {
    xdom: [Math.min(...times) - tpad, Math.max(...times) + tpad],
    ydom: [vals[0] - ypad, vals[vals.length - 1] + ypad],
    clip: [mid - half, mid + half],
    scale: Math.max(1.4826 * mad, 0.1),
  };
}

// a velocity to and from the y slider (asinh around the median scale)
const yToU = (y) => Math.asinh(y / view.scale);
const uToY = (u) => view.scale * Math.sinh(u);
function ySlider(y) {
  const [a, b] = view.ydom.map(yToU);
  return Math.round(STEPS * Math.min(1, Math.max(0, (yToU(y) - a) / (b - a))));
}
function yFromSlider(k) {
  const [a, b] = view.ydom.map(yToU);
  return uToY(a + (b - a) * k / STEPS);
}
const xSlider = (x) => Math.round(STEPS * Math.min(1, Math.max(0, (x - view.xdom[0]) / (view.xdom[1] - view.xdom[0]))));
const xFromSlider = (k) => view.xdom[0] + (view.xdom[1] - view.xdom[0]) * k / STEPS;

// the y range of the velocities kept: twice the 3-97 percentile range, or all
function refitY() {
  const dom = domains();
  if (!dom) return;
  view = { ...(view || {}), ...dom, x: view && view.x ? view.x : dom.xdom.slice() };
  view.y = $('clip').checked ? dom.clip.slice() : dom.ydom.slice();
  applyView();
}

function applyView() {
  if (!view || !window.Plotly || !$('rvplot').data) return;
  Plotly.relayout('rvplot', { 'xaxis.range': view.x.slice(), 'yaxis.range': view.y.slice() });
  syncSliders();
  foldY();
  dateAxis();
  rvOffscale();
  scheduleModel();
}

// the calendar dates on top of the time series: rjd (JD - 2400000) is
//   40587.5 at 1970-01-01 0h; years, months or days by the span shown
const RJD_UNIX = 40587.5;
const rjdDate = (r) => new Date((r - RJD_UNIX) * 86400000);
const dateRjd = (d) => d.getTime() / 86400000 + RJD_UNIX;
function dateTicks(x0, x1, width) {
  const span = x1 - x0, most = Math.max(2, Math.floor(width / 105));
  const locale = lang === 'fr' ? 'fr-CA' : 'en-GB';
  const vals = [], text = [];
  if (span > 75) {
    const step = [1, 2, 3, 6, 12, 24, 60, 120, 240].find((k) => span / (30.44 * k) <= most) || 240;
    const first = rjdDate(x0);
    for (let k = Math.ceil((first.getUTCFullYear() * 12 + first.getUTCMonth()) / step) * step; ; k += step) {
      const d = new Date(Date.UTC(Math.floor(k / 12), k % 12, 1));
      const r = dateRjd(d);
      if (r > x1) break;
      if (r < x0) continue;
      vals.push(r);
      text.push(step >= 12 ? String(d.getUTCFullYear())
        : d.toLocaleDateString(locale, { month: 'short', year: 'numeric', timeZone: 'UTC' }));
    }
  } else {
    const step = [1, 2, 5, 10, 15].find((k) => span / k <= most) || 30;
    for (let r = Math.ceil((x0 - RJD_UNIX) / step) * step + RJD_UNIX; r <= x1; r += step) {
      vals.push(r);
      text.push(rjdDate(r).toLocaleDateString(locale, { day: 'numeric', month: 'short', year: 'numeric', timeZone: 'UTC' }));
    }
  }
  return { vals, text };
}
function dateAxis() {
  const div = $('rvplot');
  if (!window.Plotly || !div.data || !div.layout || !div.layout.xaxis || !div.layout.xaxis.range) return;
  const [x0, x1] = div.layout.xaxis.range.map(Number);
  const tk = dateTicks(x0, x1, div.clientWidth - PLOT_MARGIN.l - rvRight);
  Plotly.relayout(div, { 'xaxis2.tickvals': tk.vals, 'xaxis2.ticktext': tk.text });
}

function syncSliders() {
  if (!view) return;
  const set = (id, val) => { $(id).value = val; };
  set('xlo', xSlider(view.x[0])); set('xhi', xSlider(view.x[1]));
  set('ylo', ySlider(view.y[0])); set('yhi', ySlider(view.y[1]));
  for (const ax of ['x', 'y']) {
    const lo = +$(`${ax}lo`).value / STEPS, hi = +$(`${ax}hi`).value / STEPS;
    $(`${ax}sel`).style.left = `${100 * lo}%`;
    $(`${ax}sel`).style.width = `${100 * (hi - lo)}%`;
  }
  $('rangetext').textContent = `x ${view.x[0].toFixed(0)} – ${view.x[1].toFixed(0)} · y ${view.y[0].toFixed(1)} – ${view.y[1].toFixed(1)} m/s`;
  sizeYSlider();
}

// the y slider as tall as the plotting area
function sizeYSlider() {
  const div = $('rvplot');
  const area = div.clientHeight - PLOT_MARGIN.t - PLOT_MARGIN.b;
  const h = Math.max(60, seriesAligned ? area * (1 - RES_TOP) : area);
  const box = $('yslider');
  box.style.paddingTop = `${PLOT_MARGIN.t}px`;
  box.style.height = `${h}px`;
  $('ydual').style.width = `${h}px`;
  $('ydual').style.transform = `translateY(${h}px) rotate(-90deg)`;
  $('xdual').style.marginLeft = `${PLOT_MARGIN.l}px`;
  $('xdual').style.marginRight = `${rvRight}px`;
}

function fromSliders(ax) {
  if (!view) return;
  let lo = +$(`${ax}lo`).value, hi = +$(`${ax}hi`).value;
  if (hi - lo < 5) {
    if (document.activeElement && document.activeElement.id === `${ax}lo`) lo = hi - 5; else hi = lo + 5;
  }
  view[ax] = ax === 'x' ? [xFromSlider(lo), xFromSlider(hi)] : [yFromSlider(lo), yFromSlider(hi)];
  if (ax === 'y') $('clip').checked = false;
  applyView();
}

// a zoom with the mouse moves the sliders too
function fromZoom(ev) {
  if (!view || !ev) return;
  if (ev['xaxis.range[0]'] !== undefined) view.x = [ev['xaxis.range[0]'], ev['xaxis.range[1]']];
  if (ev['yaxis.range[0]'] !== undefined) view.y = [ev['yaxis.range[0]'], ev['yaxis.range[1]']];
  if (ev['xaxis.autorange']) view.x = view.xdom.slice();
  if (ev['yaxis.autorange']) view.y = view.ydom.slice();
  syncSliders();
  if (ev['xaxis.range[0]'] !== undefined || ev['xaxis.autorange']) { dateAxis(); scheduleModel(); }
  if (ev['xaxis.range[0]'] !== undefined || ev['yaxis.range[0]'] !== undefined || ev['xaxis.autorange'] || ev['yaxis.autorange']) {
    rvOffscale();
    if (ev['yaxis.range[0]'] !== undefined || ev['yaxis.autorange']) foldY();
  }
  if (ev['yaxis.range[0]'] !== undefined || ev['yaxis.autorange']) foldY();
}

// the instruments left out: the field of the report is what counts; the
// boxes of the table and the legend of the plot write into it
let lastRV = null;
// the names of the field, separated by commas or spaces: a name with its
//   source in brackets is one ('HIRES (CLS)', 'HIRES (Teklu+ 2025)')
const namesOf = (text) => (String(text || '').match(/[^\s,(]+(?:\s*\([^)]*\))?/g) || []).map((s) => s.replace(/\s+/g, ' ').trim());
function excludedSet() {
  return new Set(namesOf($('exclude').value).map((s) => s.toUpperCase()));
}

function setExcluded(name, off) {
  const kept = namesOf($('exclude').value).filter((s) => s.toUpperCase() !== name.toUpperCase());
  if (off) kept.push(name);
  $('exclude').value = kept.join(', ');
  updateCommands();
  styleExcluded();
  checkStale();
}

function styleExcluded() {
  if (!lastRV) return;
  const ex = excludedSet();
  const off = lastRV.map((inst) => ex.has(inst.name.toUpperCase()));
  if (window.Plotly && $('rvplot').data) {
    Plotly.restyle('rvplot', { opacity: off.map((o) => (o ? 0.12 : 1)) }, off.map((_, i) => i));
    // their residuals too
    const resid = $('rvplot').data.map((d, k) => [String(d.uid || ''), k]).filter(([u]) => u.startsWith('resid-'));
    if (resid.length) Plotly.restyle('rvplot', { opacity: resid.map(([u]) => (off[+u.slice(6)] ? 0.12 : 1)) }, resid.map(([, k]) => k));
  }
  lastRV.forEach((inst, i) => {
    const row = document.querySelector(`tr[data-row="${CSS.escape(inst.name)}"]`);
    if (!row) return;
    row.classList.toggle('off', off[i]);
    row.querySelector('input').checked = !off[i];
  });
  // the time of the instruments kept, from end to end
  if (view) view.x = null;
  refitY();
}

// -----------------------------------------------------------------------------
// the quick look: a FIP of what is shown, made when it is plotted
// -----------------------------------------------------------------------------
let quick = null;   // the quick FIP: its id, state, result, and what it was of
let pview = null;   // the periods shown, and their domain

function shownOptions() {
  const asked = readOptions('detailed');
  return { files: filesNow(), target: $('target').value.trim(), root: $('root').value.trim(),
    dace: !!asked.dace, carmenes: !!asked.carmenes, vizier: !!asked.vizier, exclude: $('exclude').value,
    trend: !!asked.trend, curvature: !!asked.curvature, subtract: subtractList };
}

// several instruments: the quick FIP waits to be asked for
function quickPrompt() {
  if (quick && quick.status === 'running') api('/api/quickstop', { id: quick.id }).catch(() => {});
  clearTicks();
  quick = null;
  pview = null;
  foldShown = null; foldShownObj = null;
  $('fipcard').classList.add('on');
  if (window.Plotly) { Plotly.purge($('fipplot')); Plotly.purge($('foldplot')); }
  $('fipplot').classList.remove('on'); $('foldplot').classList.remove('on');
  $('foldbuttons').innerHTML = ''; $('foldnote').textContent = ''; $('fipeach').innerHTML = ''; $('fipstale').textContent = '';
  $('tsbuttons').innerHTML = ''; $('tsnote').textContent = ''; tsShown = null; tsOf = null;
  if (window.Plotly) Plotly.purge($('tsplot'));
  $('tsplot').classList.remove('on');
  $('remember').disabled = true; $('remstate').textContent = '';
  $('fipstatus').innerHTML = `<p class="hint">${esc(t('fip_prompt'))}</p>`
    + `<button type="button" class="go startfip">${esc(t('start_fip'))}</button>`;
}

async function startQuick() {
  clearTicks();
  $('fipcard').classList.add('on');
  $('fipstale').textContent = '';
  $('remember').disabled = true;
  $('remstate').textContent = '';
  $('fipstatus').innerHTML = `<p class="hint"><span class="spin"></span> ${esc(t('quick_noise'))}</p>`;
  try {
    const opts = shownOptions();
    quick = { ...(await api('/api/quickfip', { options: opts })), of: JSON.stringify(opts) };
    pollQuick(quick.id);
  } catch (err) {
    $('fipstatus').innerHTML = `<p class="hint bad">${esc(err.message)}</p>`;
  }
}

async function pollQuick(id) {
  if (!quick || quick.id !== id) return;   // a newer one took its place
  let state;
  try { state = await api(`/api/quickfip?id=${id}`); } catch (err) { return; }
  if (!quick || quick.id !== id) return;
  Object.assign(quick, state);
  if (state.status === 'failed') {
    $('fipstatus').innerHTML = `<p class="hint bad">${esc(state.error || 'failed')}</p>`;
    return;
  }
  // the joint FIP as soon as it is there, each instrument's as it comes
  const each = state.each || [];
  const drawn = `${state.result ? 1 : 0}/${each.length}`;
  if (state.result && quick.drawn !== drawn) {
    const first = !quick.drawn;
    drawFip(state.result, each);
    quick.drawn = drawn;
    if (first && zeroMode === 'fit' && state.result.trend_model) redrawSeries();
  }
  $('fipstatus').innerHTML = (state.result ? fipSummary(state.result, state.elapsed) : '') + quickRunning(state)
    + (state.status === 'stopped' ? `<p class="hint"><span class="bad">\u25a0</span> ${esc(t('fip_stopped'))} \u00b7 ${clock(state.elapsed)} `
      + `<button type="button" class="small startfip">${esc(t('start_fip'))}</button></p>` : '');
  $('remember').disabled = state.status !== 'done';
  if (state.status === 'running') setTimeout(() => pollQuick(id), 1000);
}

// what runs, with the bar of its sweeps
function quickRunning(state) {
  if (state.status !== 'running') return '';
  const step = { waiting: 'waiting', noise: 'quick_noise', fip1: 'quick_fip1', planets: 'quick_planets', fip2: 'quick_fip2', inst: 'quick_inst' }[state.step] || 'quick_noise';
  let bar = '';
  if (state.progress && ['fip1', 'fip2', 'inst'].includes(state.step)) {
    const p = state.progress;
    const frac = Math.min(1, p.done / Math.max(p.total, 1));
    const left = frac > 0.02 ? ` \u00b7 ~${clock(p.seconds * (1 - frac) / frac)} ${t('left_fip')}` : '';
    bar = `<span class="pbar"><span style="width:${Math.round(100 * frac)}%"></span></span><span class="ptext">${Math.round(100 * frac)} %${left}</span>`;
  }
  const detail = state.step === 'inst' && state.step_detail ? ` (${esc(state.step_detail)})` : '';
  return `<p class="hint"><span class="hourglass">\u23f3</span> ${esc(t(step))}${detail} \u00b7 ${clock(state.elapsed)}`
    + ` <button type="button" class="small stop stopfip">${esc(t('stop_fip'))}</button></p><div class="steps-bar">${bar}</div>`;
}

function fipSummary(r, elapsed) {
  const insts = Object.entries(r.instruments).map(([k, v]) => `${k} ${v}`).join(', ');
  const best = r.peaks[0];
  return `<p class="hint">${r.n} ${esc(t('nights'))} (${esc(insts)}); ${esc(t('no_gp'))}, ${r.settings.kmax} ${esc(t('signals_word'))}, `
    + `${r.settings.nsweep} ${esc(t('sweeps_word'))}, ${r.passes} ${esc(t('passes'))} \u00b7 ${clock(elapsed)}`
    + (best ? ` \u00b7 ${esc(t('strongest'))}: <b>${best.period.toFixed(4)} d</b>, FIP ${fipExp(best.family)}` : '') + '</p>'
    + `<p class="hint accel">${accelText(r)}</p>`
    + ((r.subtracted || []).length ? `<p class="resid"><b>${esc(t('residuals_of'))}</b> ${r.subtracted.map((s) => esc(s.label)).join(' ; ')}`
      + ` <button type="button" id="unsubtract" class="small">${esc(t('back_series'))}</button></p>` : '');
}

// the acceleration of the star the quick look measured (and its change)
function accelValue(v, unit) {
  const [val, lo, hi] = v;
  const err = pmErr(lo, hi, (x) => x.toPrecision(2));
  return `${val >= 0 ? '+' : '\u2212'}${Math.abs(val).toPrecision(3)}${err} ${unit}`;
}
function accelText(r) {
  const order = (r.settings || {}).trend;
  const acc = r.acceleration;
  if (order === 0) return esc(t('no_trend_fit'));
  if (!acc || !acc.accel) return '';
  let out = `${esc(t('accel_star'))}: <b>dv/dt = ${accelValue(acc.accel, 'm/s/yr')}</b> (${acc.accel_sigma.toFixed(1)}\u03c3)`;
  if (acc.jerk) out += ` \u00b7 d\u00b2v/dt\u00b2 = ${accelValue(acc.jerk, 'm/s/yr\u00b2')} (${acc.jerk_sigma.toFixed(1)}\u03c3)`;
  return out + ` \u00b7 <span class="hint">${esc(t('accel_note'))}</span>`;
}

// the colour of an instrument: its colour in the plot of the velocities
function instColour(name) {
  const i = (lastRV || []).findIndex((inst) => inst.name === name);
  return COLOURS[Math.max(0, i) % 8];
}

function drawFip(r, each) {
  each = each || [];
  const div = $('fipplot');
  div.classList.add('on');
  const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
  const lines = [], notes = [];
  for (const [name, per] of Object.entries(r.window)) {
    lines.push({ type: 'line', xref: 'x', yref: 'paper', x0: per, x1: per, y0: 0, y1: 1, line: { color: '#7a8597', width: 1, dash: 'dot' } });
    notes.push({ x: Math.log10(per), y: 1, xref: 'x', yref: 'paper', text: name, showarrow: false, yanchor: 'bottom', font: { size: 10, color: '#7a8597' } });
  }
  for (const pl of r.known) {
    lines.push({ type: 'line', xref: 'x', yref: 'paper', x0: pl.P, x1: pl.P, y0: 0, y1: 1, line: { color: '#e66767', width: 1.2, dash: 'dash' } });
    notes.push({ x: Math.log10(pl.P), y: 0.92, xref: 'x', yref: 'paper', text: pl.name, showarrow: false, xanchor: 'left', font: { size: 10, color: '#e66767' } });
  }
  lines.push({ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: 2, y1: 2, line: { color: '#a8b4ca', width: 1, dash: 'dot' } });
  // the numbered peaks; two labels too close to read, one above the other
  const placed = [];
  for (const pk of r.peak_list.filter((p) => p.named)) {
    const x = Math.log10(pk.period);
    let ay = -22;
    while (placed.some((q) => Math.abs(q.x - x) < 0.08 && q.ay === ay)) ay -= 16;
    placed.push({ x, ay });
    notes.push({ x, y: -Math.log10(Math.max(pk.family, 1e-15)), xref: 'x', yref: 'y',
      text: `<b>#${pk.id}</b>`, showarrow: true, arrowhead: 0, ax: 0, ay, font: { size: 12, color: '#ffffff' },
      arrowcolor: '#a8b4ca', hovertext: `#${pk.id}: ${pk.period.toFixed(4)} d, FIP ${fipExp(pk.family)}` });
  }
  const traces = [
    { x: r.period, y: r.alone, name: t('alone'), type: 'scatter', mode: 'lines', line: { color: '#8a93a3', width: 1 },
      legendrank: 1, hovertemplate: '%{x:.4f} d<br>-log10 FIP %{y:.2f}<extra></extra>' },
  ];
  // each instrument on its own: its FIP of the period or any of its aliases
  const valid = each.filter((x) => !x.skipped);
  $('fipv-each').disabled = !valid.length;
  document.querySelectorAll('[data-fipview]').forEach((b) => b.classList.toggle('on', b.dataset.fipview === (valid.length ? fipView : 'joint')));
  for (const one of valid) {
    traces.push({ x: one.period, y: one.family, name: `${one.name} ${t('inst_alone')}`, type: 'scatter', mode: 'lines',
      uid: `each-${safeId(one.name)}`, visible: fipView === 'each',
      line: { color: instColour(one.name), width: 1 }, opacity: 0.85, legendrank: 10 + traces.length,
      hovertemplate: `${esc(one.name)}<br>%{x:.4f} d<br>-log10 FIP %{y:.2f}<extra></extra>` });
  }
  // the joint FIP on top, in white (blue is the first instrument's colour)
  traces.push({ x: r.period, y: r.family, name: t('family'), type: 'scatter', mode: 'lines', line: { color: '#e8eef8', width: 1.6 },
    legendrank: 2, hovertemplate: '%{x:.4f} d<br>-log10 FIP %{y:.2f}<extra></extra>' });
  // the periods shown kept when the curves of the instruments come in
  const keepView = pview && pview.of === quick.id ? pview.p.slice() : null;
  const xaxis = { ...axis, type: 'log', title: t('period_axis') };
  if (keepView) xaxis.range = [Math.log10(keepView[0]), Math.log10(keepView[1])];
  Plotly.newPlot(div, traces, {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)', font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
    margin: { ...PLOT_MARGIN, t: 40 }, legend: { orientation: 'h', x: 0, y: 1.08, yanchor: 'bottom' }, shapes: lines, annotations: notes,
    xaxis, yaxis: { ...axis, title: t('fip_axis'), rangemode: 'tozero' },
  }, { responsive: true, displaylogo: false }).then(() => {
    div.removeAllListeners && div.removeAllListeners('plotly_click');
    // a click on the FIP: a fold at the dip of the period nearest
    div.on('plotly_click', (ev) => {
      const pt = ev.points && ev.points[0];
      if (pt && pt.x > 0) forceFold(+pt.x, true);
    });
    div.removeAllListeners && div.removeAllListeners('plotly_relayout');
    div.on('plotly_relayout', (ev) => {
      if (!pview) return;
      if (ev['xaxis.range[0]'] !== undefined) pview.p = [10 ** ev['xaxis.range[0]'], 10 ** ev['xaxis.range[1]']];
      if (ev['xaxis.autorange']) pview.p = pview.dom.slice();
      syncPeriods();
    });
    drawStages();
  });
  fipNotes = notes;
  const per = r.period;
  pview = { of: quick.id, dom: [per[0], per[per.length - 1]], p: keepView || [per[0], per[per.length - 1]] };
  syncPeriods();
  drawEach(each);
  if (keepView) return;
  // the folds at the numbered peaks
  renderFoldButtons();
  if (r.folds && r.folds.length) showFold(r.folds[0].id);
  // a transit in TESS at the peaks below a FIP of 1 %, and at the TOIs
  renderTransits(true);
}

// -----------------------------------------------------------------------------
// a transit in TESS: at each peak with a FIP (of the period or any alias)
//   below 1 %, about the conjunction of its fold, and at the ephemeris of
//   each transiting planet and TOI of the star; the light curve high-passed
//   and folded, its points and their medians in bins, the depths of a
//   1 Earth-radius and a 1 Jupiter-radius planet
// -----------------------------------------------------------------------------
const TRANSIT_FIP = 0.1;       // a transit looked for below this FIP
let tsShown = null;            // the candidate shown
let tsOf = null;               // the quick FIP its candidates are of
function transitCandidates() {
  const r = quick && quick.result;
  if (!r) return [];
  const times = (lastRV || []).flatMap((inst) => inst.time);
  const span = times.length ? Math.max(1, Math.max(...times) - Math.min(...times)) : 1000;
  const out = [];
  // a peak below a FIP of 10 %: a transit found would make it strong
  (r.peak_list || []).filter((pk) => pk.family < TRANSIT_FIP).forEach((pk) => {
    const f = (r.folds || []).find((x) => x.id === pk.id) || {};
    out.push({ key: `p${pk.id}`, label: `#${pk.id} \u00b7 ${pk.period.toFixed(4)} d`, period: f.period || pk.period, tc: f.tc, tc_err: f.tc_err,
      p_err: pk.period ** 2 / (4 * span), name: `#${pk.id}`, kind: 'rv' });
  });
  (r.transits || []).forEach((tr, i) => {
    if (tr.P && tr.tc !== null && tr.tc !== undefined) {
      out.push({ key: `t${i}`, label: `\u25d0 ${tr.name} \u00b7 ${(+tr.P).toFixed(4)} d`, period: +tr.P, tc: tr.tc, tc_err: tr.tc_err,
        p_err: tr.P_err, name: tr.name, kind: 'toi' });
    }
  });
  return out;
}
function renderTransits(auto) {
  const cands = transitCandidates();
  const box = $('tsbuttons');
  if (!cands.length) {
    box.innerHTML = ''; $('tsnote').textContent = quick && quick.result ? t('ts_none') : '';
    $('tsplot').classList.remove('on'); if (window.Plotly) Plotly.purge($('tsplot'));
    return;
  }
  box.innerHTML = cands.map((c) => `<button type="button" class="small${c.kind === 'toi' ? ' transit' : ''}" data-tsearch="${esc(c.key)}">${esc(c.label)}</button>`).join('');
  if (auto && tsOf !== quick.id) { tsOf = quick.id; showTransit(cands[0].key, true); }
  document.querySelectorAll('[data-tsearch]').forEach((b) => b.classList.toggle('on', b.dataset.tsearch === tsShown));
}
async function showTransit(key, fetch) {
  const c = transitCandidates().find((x) => x.key === key);
  if (!c) return;
  tsShown = key;
  document.querySelectorAll('[data-tsearch]').forEach((b) => b.classList.toggle('on', b.dataset.tsearch === key));
  const note = $('tsnote');
  note.innerHTML = `<span class="spin"></span> ${esc(t(fetch ? 'ts_fetching' : 'ts_searching'))}`;
  try {
    const res = await api('/api/transit', { options: { target: $('target').value.trim(), root: $('root').value.trim(), period: c.period,
      tc: c.tc, tc_err: c.tc_err, p_err: c.p_err, name: c.name, kind: c.kind, fetch: !!fetch, mstar: +$('mstar').value || null } });
    if (tsShown !== key) return;
    if (res.missing) {
      // TESS never looked at the star (nothing to fetch), has it in its
      //   full frames only, or its light curve is not on this machine yet
      note.innerHTML = res.reason === 'unobserved' ? `<span>${esc(t('ts_unobserved'))}</span>`
        : res.reason === 'frames' ? `<span>${esc(t('ts_frames'))} ${esc((res.sectors || []).join(', '))}</span>`
          : `${esc(t('ts_missing'))} (${esc(res.lc || '')}) <button type="button" class="small" id="tsfetch">${esc(t('ts_fetch'))}</button>`;
      $('tsplot').classList.remove('on');
      return;
    }
    drawTransit(res, c);
  } catch (err) {
    note.innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
}
// the TESS panel as a fold or as the time series of each sector
let tsView = 'fold';
try { tsView = localStorage.getItem('koloa-tsview') || 'fold'; } catch (err) { /* no storage */ }
let tsLast = null;
function drawTransit(res, c) {
  tsLast = { res, c };
  document.querySelectorAll('[data-tsview]').forEach((x) => x.classList.toggle('on', x.dataset.tsview === tsView));
  const b = res.best;
  // the light curve folded at the period found (the scan's best, when the
  //   period was scanned), the box and the window about it
  const P = res.fold_period || res.period;
  const fmt = (v, d) => (v === null || v === undefined ? '?' : (+v).toFixed(d));
  const verdict = res.plausible ? `<span class="ok">\u2691 ${esc(t('ts_plausible'))}</span>` : `<span class="hint">${esc(t('ts_not'))}</span>`;
  $('tsnote').innerHTML = `${verdict}: ${esc(res.why)}.`
    // the depth, duration and radius of the box fitted to the medians
    + (res.fit ? ` ${esc(t('ts_fit'))} ${fmt(res.fit.depth, 3)} \u00b1 ${fmt(res.fit.depth_err, 3)} ppt, ${fmt(res.fit.duration, 1)} h, Rp \u2248 ${fmt(res.fit.radius, 2)} R\u2295`
      : b ? ` ${esc(t('ts_box'))} ${fmt(b.depth, 3)} \u00b1 ${fmt(b.depth_err, 3)} ppt, ${fmt(b.duration, 1)} h, Rp \u2248 ${fmt(b.radius, 2)} R\u2295` : '')
    + ` \u00b7 R\u2605 ${fmt(res.star.radius, 3)} R\u2609 (${esc(res.star.radius_source)}), ${t('ts_expected')} ${fmt(res.duration, 1)} h`
    + (teqText(P, res.star) ? ` \u00b7${teqText(P, res.star)}` : '')
    + ` \u00b7 TESS ${esc((res.sectors || []).join(', '))} (${esc(res.lc)})`
    + (res.window ? ` \u00b7 ${t('ts_window')} \u00b1${fmt(res.window, 1)} h` : ` \u00b7 ${t('ts_whole')}`)
    // the period scanned: how many, and where the box is best
    + (res.scan ? ` \u00b7 ${t('ts_scan')} ${res.scan.n} ${t('ts_periods')} ${fmt(res.scan.low, 5)} - ${fmt(res.scan.high, 5)} d,`
      + ` ${t('ts_folded')} ${fmt(P, 5)} d${b && b.period_offset !== undefined ? ` (${b.period_offset >= 0 ? '+' : ''}${fmt(b.period_offset, 1)}\u03c3)` : ''}` : '');
  const div = $('tsplot');
  div.classList.add('on');
  div.classList.toggle('series', tsView === 'series');
  if (!window.Plotly) return;
  if (tsView === 'series') { drawTransitSeries(res); return; }
  div.style.height = '';
  // the phase at the bottom (0 at the transit expected, or found), the
  //   hours from it on top
  const ph = (h) => h / (P * 24);
  const traces = [{ x: res.hours.map(ph), y: res.flux, type: 'scattergl', mode: 'markers', name: t('ts_points'),
    marker: { size: 3, color: 'rgba(170,185,210,0.35)' }, hoverinfo: 'skip' }];
  const bins = res.bins || [];
  traces.push({ x: bins.map((v) => ph(v[0])), y: bins.map((v) => v[1]), type: 'scatter', mode: 'markers+lines', name: t('ts_bins'),
    error_y: { type: 'data', array: bins.map((v) => v[2]), visible: true, thickness: 1, width: 0, color: '#62c2ff' },
    marker: { size: 6, color: '#62c2ff' }, line: { color: '#62c2ff', width: 1 }, customdata: bins.map((v) => v[0]),
    hovertemplate: `${t('ts_phase')} %{x:.4f} \u00b7 %{customdata:.2f} h<br>%{y:.3f} ppt<extra></extra>` });
  // nothing to see on the axis of the hours: Plotly draws an axis a trace uses
  traces.push({ x: [-0.5, 0.5], y: [null, null], xaxis: 'x2', type: 'scatter', mode: 'markers', showlegend: false, hoverinfo: 'skip',
    marker: { opacity: 0 } });
  // the box fitted to the medians (else the box of the search), about the
  //   centre shown
  if (res.fit) {
    const f = res.fit, h = f.duration / 2;
    traces.push({ x: [f.centre - 3 * h, f.centre - h, f.centre - h, f.centre + h, f.centre + h, f.centre + 3 * h].map(ph),
      y: [f.level, f.level, f.level - f.depth, f.level - f.depth, f.level, f.level], type: 'scatter',
      mode: 'lines', name: t('ts_fitname'), line: { color: '#f5a524', width: 1.6 }, hoverinfo: 'skip' });
  } else if (b) {
    const off = (((b.centre - res.shown_centre) / P + 0.5) % 1 + 1) % 1 * P * 24 - 0.5 * P * 24;
    const h = b.duration / 2;
    traces.push({ x: [off - 3 * h, off - h, off - h, off + h, off + h, off + 3 * h].map(ph), y: [0, 0, -b.depth, -b.depth, 0, 0], type: 'scatter',
      mode: 'lines', name: t('ts_boxname'), line: { color: '#f5a524', width: 1.6 }, hoverinfo: 'skip' });
  }
  // the depths of a 1 Earth-radius and a 1 Jupiter-radius planet
  const span = Math.max(4 * res.duration, res.window ? 1.2 * res.window : 0, 6);
  const xr = [-Math.min(span, P * 12), Math.min(span, P * 12)];
  const vis = bins.filter((v) => v[0] >= xr[0] && v[0] <= xr[1]);
  const err = vis.length ? vis.map((v) => v[2]).sort((p, q) => p - q)[Math.floor(vis.length / 2)] : 0.1;
  let lo = Math.min(...vis.map((v) => v[1]), -(res.fit ? res.fit.depth : b ? b.depth : 0)) - 4 * err, hi = Math.max(...vis.map((v) => v[1]), 0) + 4 * err;
  lo = Math.min(lo, -1.25 * res.depth_earth);
  const lines = [[res.depth_earth, '1 R\u2295'], [res.depth_jupiter, '1 R\u2643']];
  const shapes = [], notes = [];
  lines.forEach(([d, name]) => {
    if (-d >= lo) {
      shapes.push({ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: -d, y1: -d, line: { color: '#e66767', width: 1, dash: 'dash' } });
      notes.push({ xref: 'paper', x: 1, y: -d, text: `${name}: ${fmt(d, d < 1 ? 3 : 1)} ppt`, showarrow: false, xanchor: 'right', yanchor: 'bottom', font: { size: 10, color: '#e66767' } });
    } else {
      notes.push({ xref: 'paper', yref: 'paper', x: 1, y: 0, text: `${name}: ${fmt(d, 1)} ppt \u2193`, showarrow: false, xanchor: 'right', yanchor: 'bottom', font: { size: 10, color: '#e66767' } });
    }
  });
  // where the transit is expected (the ephemeris carried to TESS), shaded
  if (res.window && res.t0 !== null && res.t0 !== undefined) {
    const at = (((res.t0 - res.shown_centre) / P + 0.5) % 1 + 1) % 1 * P * 24 - 0.5 * P * 24;
    shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: ph(at - res.window), x1: ph(at + res.window), y0: 0, y1: 1,
      fillcolor: 'rgba(98,194,255,0.06)', line: { width: 0 } });
  }
  const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
  Plotly.newPlot(div, traces, {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)', font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
    margin: { ...PLOT_MARGIN, t: 92 }, legend: { orientation: 'h', x: 0, y: 1, yref: 'container', yanchor: 'top' }, shapes, annotations: notes,
    xaxis: { ...axis, title: t('ts_phase'), range: xr.map(ph) },
    xaxis2: { ...axis, title: { text: t(res.plausible ? 'ts_axis_found' : 'ts_axis'), standoff: 4 }, overlaying: 'x', matches: 'x', side: 'top',
      showgrid: false, zeroline: false, tickmode: 'array', ...hourTicks(xr.map(ph), P) },
    yaxis: { ...axis, title: t('ts_flux'), range: [lo, hi] },
  }, { responsive: true, displaylogo: false }).then(() => {
    // the hours on top follow a zoom
    div.removeAllListeners && div.removeAllListeners('plotly_relayout');
    div.on('plotly_relayout', (ev) => {
      // a zoom (not the ticks set here)
      if (ev['xaxis2.tickvals'] !== undefined) return;
      if (!['xaxis.range', 'xaxis.range[0]', 'xaxis.autorange', 'xaxis2.range[0]', 'xaxis2.autorange'].some((key) => ev[key] !== undefined)) return;
      const r = div.layout.xaxis.range.map(Number);
      const tk = hourTicks(r, P);
      Plotly.relayout(div, { 'xaxis2.tickvals': tk.tickvals, 'xaxis2.ticktext': tk.ticktext });
    });
  });
}

// the high-passed light curve of each sector on its own, each putative
//   transit shaded with the fitted box on it and its own depth: whether
//   the transits agree, or one event makes the fold
function drawTransitSeries(res) {
  const div = $('tsplot');
  const sectors = [...new Set(res.sector || [])].sort((p, q) => p - q);
  if (!sectors.length) { Plotly.purge(div); return; }
  const cols = sectors.length > 1 ? 2 : 1, rows = Math.ceil(sectors.length / cols);
  div.style.height = `${170 * rows + 80}px`;
  const fit = res.fit, b = res.best;
  const depth = fit ? fit.depth : b ? b.depth : 0;
  const width = (fit ? fit.duration : b ? b.duration : res.duration) / 24;
  // one y range for all: the fold's
  const fl = (res.flux || []).slice().sort((p, q) => p - q);
  const q = (f) => fl[Math.min(fl.length - 1, Math.max(0, Math.floor(f * fl.length)))];
  const lo = Math.min(q(0.005), -1.6 * depth), hi = q(0.995);
  const traces = [], shapes = [], notes = [], layout = {};
  const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
  sectors.forEach((sec, k) => {
    const xa = k ? `x${k + 1}` : 'x', ya = k ? `y${k + 1}` : 'y';
    const idx = res.sector.map((v, i) => (v === sec ? i : -1)).filter((i) => i >= 0);
    const tt = idx.map((i) => res.time[i]), ff = idx.map((i) => res.flux[i]);
    traces.push({ x: tt, y: ff, xaxis: xa, yaxis: ya, type: 'scattergl', mode: 'markers', showlegend: k === 0, name: t('ts_points'),
      legendgroup: 'pts', marker: { size: 2, color: 'rgba(170,185,210,0.35)' }, hoverinfo: 'skip' });
    // medians in bins of an hour
    const bins = new Map();
    tt.forEach((x, i) => { const key = Math.floor(x * 24); if (!bins.has(key)) bins.set(key, []); bins.get(key).push([x, ff[i]]); });
    const mx = [], my = [];
    let prev = null;
    [...bins.keys()].sort((p, r) => p - r).forEach((key) => {
      const v = bins.get(key);
      if (v.length < 5) return;
      // the line broken across a gap (more than two hours without a bin)
      if (prev !== null && key - prev > 2) { mx.push(null); my.push(null); }
      prev = key;
      const ys = v.map((p) => p[1]).sort((p, r) => p - r);
      mx.push(v.reduce((acc, p) => acc + p[0], 0) / v.length); my.push(ys[Math.floor(ys.length / 2)]);
    });
    traces.push({ x: mx, y: my, xaxis: xa, yaxis: ya, type: 'scatter', mode: 'lines', showlegend: k === 0, name: t('ts_bins_hour'),
      legendgroup: 'bins', line: { color: '#62c2ff', width: 1 }, hovertemplate: '%{x:.3f}<br>%{y:.3f} ppt<extra></extra>' });
    const t0 = Math.min(...tt), t1 = Math.max(...tt);
    // each putative transit of the sector: shaded, the fitted box on it,
    //   its own depth
    (res.events || []).filter((e) => e.time >= t0 - width && e.time <= t1 + width).forEach((e, j) => {
      shapes.push({ type: 'rect', xref: xa, yref: `${ya} domain`, x0: e.time - width / 2, x1: e.time + width / 2, y0: 0, y1: 1,
        fillcolor: 'rgba(245,165,36,0.14)', line: { width: 0 } });
      const h = width / 2;
      traces.push({ x: [e.time - 3 * h, e.time - h, e.time - h, e.time + h, e.time + h, e.time + 3 * h],
        y: [e.level, e.level, e.level - depth, e.level - depth, e.level, e.level], xaxis: xa, yaxis: ya, type: 'scatter', mode: 'lines',
        showlegend: k === 0 && j === 0, legendgroup: 'box', name: t('ts_fitname'), line: { color: '#f5a524', width: 1.4 }, hoverinfo: 'skip' });
      notes.push({ xref: xa, yref: `${ya} domain`, x: e.time, y: 0.02, text: `${e.depth.toFixed(2)}`, showarrow: false, yanchor: 'bottom',
        font: { size: 9, color: e.depth > 0.5 * depth ? '#f5a524' : '#7a8597' } });
    });
    notes.push({ xref: `${xa} domain`, yref: `${ya} domain`, x: 0.01, y: 0.98, text: `${t('ts_sector')} ${sec}`, showarrow: false,
      xanchor: 'left', yanchor: 'top', font: { size: 10, color: '#a8b4ca' } });
    layout[k ? `xaxis${k + 1}` : 'xaxis'] = { ...axis, tickformat: '.0f', exponentformat: 'none', nticks: 6 };
    layout[k ? `yaxis${k + 1}` : 'yaxis'] = { ...axis, range: [lo, hi], title: k % cols === 0 ? 'ppt' : '' };
  });
  Plotly.newPlot(div, traces, {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)', font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
    margin: { l: 56, r: 12, t: 40, b: 30 }, legend: { orientation: 'h', x: 0, y: 1, yref: 'container', yanchor: 'top' },
    grid: { rows, columns: cols, pattern: 'independent', xgap: 0.06, ygap: 0.28 }, shapes, annotations: notes, ...layout,
  }, { responsive: true, displaylogo: false });
}

// the ticks of the hours on top of a fold in phase: a round step for the
//   hours shown, at their phase
function hourTicks(range, P) {
  const h0 = range[0] * P * 24, h1 = range[1] * P * 24;
  const step = [0.25, 0.5, 1, 2, 3, 6, 12, 24, 48, 96, 168, 336, 720].find((k) => (h1 - h0) / k <= 10) || 1440;
  const tickvals = [], ticktext = [];
  for (let h = Math.ceil(h0 / step) * step; h <= h1 + 1e-9; h += step) {
    tickvals.push(h / (P * 24));
    ticktext.push(step < 1 ? h.toFixed(2) : String(Math.round(h)));
  }
  return { tickvals, ticktext };
}

// the FIP: the period alone and with its aliases, or each instrument too
let fipView = 'joint';
try { fipView = localStorage.getItem('koloa-fipview') || 'joint'; } catch (err) { /* no storage */ }
function showFipView() {
  const div = $('fipplot');
  document.querySelectorAll('[data-fipview]').forEach((b) => b.classList.toggle('on', b.dataset.fipview === fipView));
  if (!window.Plotly || !div.data) return;
  const idx = div.data.map((d, k) => (String(d.uid || '').startsWith('each-') ? k : -1)).filter((k) => k >= 0);
  if (idx.length) Plotly.restyle(div, { visible: fipView === 'each' }, idx);
}

// a button per fold: the numbered peaks, then the periods asked
// a fold by its name: its transiting planet or TOI, its known planet, or
//   its peak
const foldLabel = (f) => (f.transit ? `\u25d0 ${esc(f.transit.name)}` : f.known ? `\u2605 ${esc(f.known.name)}` : `#${f.id}`);

function renderFoldButtons() {
  const folds = (quick && quick.result && quick.result.folds) || [];
  const known = (quick && quick.result && quick.result.known) || [];
  const transits = (quick && quick.result && quick.result.transits) || [];
  const label = foldLabel;
  // a tick takes the signal out: the colour of its subtraction
  const order = (f) => ticked.findIndex((tk) => tk.id === f.id);
  $('foldbuttons').innerHTML = folds.map((f) => `<span class="foldpick${order(f) >= 0 ? ' ticked' : ''}"${order(f) >= 0 ? ` style="--tick:${stageColour(order(f))}"` : ''}>`
    + `<input type="checkbox" data-tick="${f.id}"${order(f) >= 0 ? ' checked' : ''} title="${esc(t('tick_tip'))}">`
    + `<button type="button" class="small${f.forced ? ' forced' : ''}${f.known ? ' known' : ''}${f.transit ? ' transit' : ''}" data-fold="${f.id}"`
    + `${f.forced ? ` title="${esc(t('asked'))}"` : ''}>${label(f)} \u00b7 ${f.period.toFixed(4)} d</button></span>`).join('')
    // the known planets, at their published period whatever the FIP says
    + known.filter((pl) => !folds.some((f) => f.known && f.known.name === pl.name))
      .map((pl) => `<button type="button" class="small known" data-known="${esc(pl.name)}" title="${esc(t('known_fold'))}">`
        + `\u2605 ${esc(pl.name)} \u00b7 ${(+pl.P).toPrecision(6)} d</button>`).join('')
    // the transit ephemerides: phase 0 at the transit
    + transits.filter((tr) => !folds.some((f) => f.transit && f.transit.name === tr.name))
      .map((tr) => `<button type="button" class="small transit" data-transit="${esc(tr.name)}" title="${esc(t('transit_fold'))}">`
        + `\u25d0 ${esc(tr.name)} \u00b7 ${esc(t('transit_word'))} (${esc(tr.reference)})</button>`).join('');
  document.querySelectorAll('[data-fold]').forEach((b) => b.classList.toggle('on', +b.dataset.fold === foldShown));
}

// a fold at a period of one's own: clicked on the FIP (then at the dip of
//   the period nearest) or typed (as it is)
async function forceFold(period, snap, known, transit) {
  if (!quick || !quick.result || !(period > 0 || known || transit)) return;
  $('foldnote').innerHTML = `<span class="spin"></span> ${esc(t('folding'))}`;
  try {
    const last = lastStage();
    const res = await api('/api/fold', { quick: quick.id, period, snap: !!snap, kind: foldModel, options: shownOptions(),
      known: known || '', transit: transit || '', minus: ticked.map(tickModel).filter(Boolean), minus_tag: tickKey(ticked),
      snapq: last && last.job.id ? last.job.id : '' });
    const folds = quick.result.folds = quick.result.folds || [];
    const old = folds.findIndex((x) => x.id === res.fold.id);
    if (old >= 0) folds[old] = res.fold; else folds.push(res.fold);
    $('foldp').value = res.fold.period.toFixed(4);
    renderFoldButtons();
    showFold(res.fold.id);
  } catch (err) {
    $('foldnote').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
}

// the Keplerian orbit of a fold, fitted when first asked
let keplerBusy = null;
async function keplerOf(base) {
  if (keplerBusy === base.id) return;
  keplerBusy = base.id;
  $('foldnote').innerHTML = `<span class="spin"></span> ${esc(t('fitting_kepler'))}`;
  try {
    const res = await api('/api/fold', { quick: quick.id, id: base.transit ? undefined : base.id, kind: 'kepler', options: shownOptions(),
      transit: base.transit ? base.transit.name : '', ...minusOf(base.id) });
    base.kepler = res.fold.kepler;
    if (foldShown === base.id || foldShown === null) showFold(base.id);
  } catch (err) {
    $('foldnote').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  } finally {
    keplerBusy = null;
  }
}


// each instrument's nights and best peaks, under the FIP
function drawEach(each) {
  if (!each.length) { $('fipeach').innerHTML = ''; return; }
  $('fipeach').innerHTML = `<table class="mini"><tr><th>${esc(t('each_title'))}</th><th>${esc(t('nights'))}</th><th>${esc(t('best_peaks'))}</th></tr>`
    + each.map((one) => `<tr><td><span class="swatch" style="background:${instColour(one.name)}"></span>${esc(one.name)}</td>`
      + `<td class="num">${one.n}</td><td class="peaks">${one.skipped ? `<span class="hint">${esc(t('few_nights'))}</span>`
        : (one.peaks || []).map((pk) => `${pk.period.toFixed(4)} (${fipExp(pk.family)})`).join(' \u00b7 ') || esc(t('none_found'))}</td></tr>`).join('')
    + '</table>';
}

let foldShown = null;    // the id of the fold shown
let foldColour = 'inst'; // its points by instrument, date or BERV
try { foldColour = localStorage.getItem('koloa-foldcolour') || 'inst'; } catch (err) { /* no storage */ }
const RED = '#ff3b3b';
// a time as a calendar date and hour (of the barycentric time: within
//   minutes of UT)
const dateText = (r) => `${rjdDate(r).toISOString().slice(0, 16).replace('T', ' ')} (BJD)`;

let foldModel = 'sine';  // the fold's sinusoid, or its Keplerian orbit
try { foldModel = localStorage.getItem('koloa-foldmodel') || 'sine'; } catch (err) { /* no storage */ }
let subtractList = [];   // the solutions taken out of the series (a result recalled from before the ticks)
// the signals ticked in the list of folds: taken out one after the other in
//   the order ticked, their sum drawn on the series, and the FIP of what is
//   left after each drawn over the FIP of the series, each in its colour
let ticked = [];             // [{ id, kind }], kind: 'sine' or 'kepler'
const stages = new Map();    // the FIPs of the residuals: key -> { key, ids, label, job, error }
const STAGE_COLOURS = ['#f5a524', '#62c2ff', '#7ee787', '#ff7b72', '#d2a8ff', '#ffd866'];
const foldById = (id) => ((quick && quick.result && quick.result.folds) || []).find((x) => x.id === +id);
// a fold by its name, as text (a legend, a note)
const plainLabel = (f) => (!f ? '?' : f.transit ? `◐ ${f.transit.name}` : f.known ? `★ ${f.known.name}` : `#${f.id}`);
// the fold of a ticked signal as ticked: its Keplerian orbit, or its sinusoid
const tickFold = (tk) => { const f = foldById(tk.id); return !f ? null : tk.kind === 'kepler' && f.kepler ? { ...f.kepler, id: f.id } : f; };
// its solution, as the server takes it out (its signal, not its offsets and trend)
function tickModel(tk) {
  const f = tickFold(tk);
  if (!f || !f.model) return null;
  const kind = f.kind === 'kepler' ? `${t('fmodel_kepler_short')}, e ${f.e.toFixed(2)}` : t('fmodel_sine');
  return { ...f.model, period: f.model.period || f.period, label: `${plainLabel(foldById(tk.id))} ${f.period.toFixed(4)} d (${kind}, K ${f.K.toFixed(2)} m/s)` };
}
const tickKey = (list) => list.map((tk) => `${tk.id}${tk.kind === 'kepler' ? 'k' : 's'}`).join('+');
// what a fold is fitted without: the ticked signals other than itself
const othersOf = (id) => ticked.filter((tk) => tk.id !== +id);
const minusOf = (id) => { const o = othersOf(id); return { minus: o.map(tickModel).filter(Boolean), minus_tag: tickKey(o) }; };
// the consecutive subtractions: the first ticked, the first two, ...
const stageLists = () => ticked.map((_, k) => ticked.slice(0, k + 1));
const stageColour = (k) => STAGE_COLOURS[k % STAGE_COLOURS.length];
const lastStage = () => { const done = stageLists().map((l) => stages.get(tickKey(l))).filter((st) => st && st.job && st.job.result); return done.length ? done[done.length - 1] : null; };

// a fold fitted again on the series without the other ticked signals (the
//   server keeps what it was fitted without: nothing is done when the same)
async function refitFold(id, kind) {
  const base = foldById(id);
  if (!base || !quick) return;
  const want = tickKey(othersOf(id));
  if ((base.minus_tag || '') === want && (kind !== 'kepler' || base.kepler)) return;
  try {
    const res = await api('/api/fold', { quick: quick.id, id: +id, kind: kind === 'kepler' ? 'kepler' : 'sine', options: shownOptions(), ...minusOf(id) });
    const folds = quick.result.folds;
    const at = folds.findIndex((x) => x.id === +id);
    res.fold.minus_tag = want;
    if (at >= 0) folds[at] = res.fold;
  } catch (err) {
    base.minus_tag = want;   // not asked again and again
    $('foldnote').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
}
// each ticked signal fitted without the others (one pass, in the order ticked)
async function refitTicked() {
  for (const tk of ticked.slice()) await refitFold(tk.id, tk.kind);
}

// a signal ticked or unticked: the folds fitted again without the others,
//   the FIP of what is left after each, the sum on the series
let tickBusy = false;
async function toggleTick(id) {
  id = +id;
  if (tickBusy || !foldById(id)) { renderFoldButtons(); return; }
  tickBusy = true;
  const at = ticked.findIndex((tk) => tk.id === id);
  if (at >= 0) ticked.splice(at, 1);
  else ticked.push({ id, kind: foldModel === 'kepler' ? 'kepler' : 'sine' });
  renderFoldButtons();
  $('stagestatus').innerHTML = `<p class="hint"><span class="spin"></span> ${esc(t('tick_fitting'))}</p>`;
  try {
    await refitTicked();
    updateStages();
  } finally {
    tickBusy = false;
  }
  if (foldShown !== null && foldById(foldShown)) showFold(foldShown);
  else if ($('foldoverlay').checked) redrawSeries();
}

// the FIP of what is left after each consecutive subtraction: asked once
//   (kept by what was ticked), followed until done
function updateStages() {
  if (!quick || !quick.result) return;
  for (const list of stageLists()) {
    const key = tickKey(list);
    if (stages.has(key)) continue;
    const models = list.map(tickModel);
    if (models.some((m) => !m)) continue;
    const st = { key, ids: list.map((tk) => tk.id), label: list.map((tk) => plainLabel(foldById(tk.id))).join(' + '), job: null, error: null };
    stages.set(key, st);
    api('/api/quickfip', { options: { ...shownOptions(), subtract: models, stage: true } })
      .then((job) => { st.job = job; pollStage(st); })
      .catch((err) => { st.error = err.message; drawStages(); });
  }
  drawStages();
}
async function pollStage(st) {
  if (stages.get(st.key) !== st || !st.job) return;
  let state;
  try { state = await api(`/api/quickfip?id=${st.job.id}`); } catch (err) { return; }
  if (stages.get(st.key) !== st) return;
  Object.assign(st.job, state);
  drawStages();
  if (state.status === 'running') setTimeout(() => pollStage(st), 2000);
}
// the ticks and the FIPs of the residuals forgotten (a new FIP of the
//   series, another star): those still running are stopped
function clearTicks() {
  for (const st of stages.values()) {
    if (st.job && st.job.status === 'running' && st.job.id) api('/api/quickstop', { id: st.job.id }).catch(() => {});
  }
  stages.clear();
  ticked = [];
  if ($('stagestatus')) $('stagestatus').innerHTML = '';
}

// the FIPs of the residuals over the FIP of the series, each in the colour
//   of its subtraction, and where each stands (a line under the summary)
let fipNotes = [];   // the annotations of the FIP of the series
function drawStages() {
  const div = $('fipplot');
  const lists = stageLists();
  const lines = [];
  const traces = [];
  let notes = [];
  lists.forEach((list, k) => {
    const st = stages.get(tickKey(list));
    const colour = stageColour(k);
    const name = `${t('without')} ${st ? st.label : ''}`;
    const chip = `<span class="stagechip" style="--tick:${colour}">${esc(name)}</span>`;
    if (!st) return;
    const job = st.job || {};
    if (st.error || job.status === 'failed') { lines.push(`${chip} <span class="bad">${esc(st.error || job.error || 'failed')}</span>`); return; }
    const res = job.result;
    if (job.status === 'running' || !st.job) {
      const p = job.progress;
      const frac = p ? Math.min(1, p.done / Math.max(p.total, 1)) : 0;
      lines.push(`${chip} <span class="spin"></span> ${esc(t('stage_running'))}${p ? ` ${Math.round(100 * frac)} %` : ''}${job.elapsed ? ` · ${clock(job.elapsed)}` : ''}`);
    } else if (res) {
      const best = (res.peaks || [])[0];
      lines.push(`${chip} ${best ? `${esc(t('strongest'))}: <b>${best.period.toFixed(4)} d</b>, FIP ${fipExp(best.family)}` : esc(t('none_found'))}`);
    } else if (job.status === 'stopped') lines.push(`${chip} ${esc(t('fip_stopped'))}`);
    if (res && res.period) {
      traces.push({ x: res.period, y: res.family, name, type: 'scatter', mode: 'lines', uid: `stage-${k}`, line: { color: colour, width: 1.5 },
        legendrank: 3 + k, hovertemplate: `${esc(name)}<br>%{x:.4f} d<br>-log10 FIP %{y:.2f}<extra></extra>` });
      // the peaks of what is left after the last subtraction, at their periods
      if (k === lists.length - 1) {
        notes = (res.peak_list || []).filter((p) => p.named).map((pk) => ({ x: Math.log10(pk.period), y: -Math.log10(Math.max(pk.family, 1e-15)),
          xref: 'x', yref: 'y', text: `${pk.period.toFixed(pk.period < 100 ? 2 : 0)} d`, showarrow: true, arrowhead: 0, ax: 0, ay: -16,
          font: { size: 10, color: colour }, arrowcolor: colour, hovertext: `${pk.period.toFixed(4)} d, FIP ${fipExp(pk.family)}` }));
      }
    }
  });
  if ($('stagestatus')) $('stagestatus').innerHTML = lines.map((l) => `<p class="hint stageline">${l}</p>`).join('');
  if (!window.Plotly || !div.data) return;
  const old = div.data.map((d, k) => (String(d.uid || '').startsWith('stage-') ? k : -1)).filter((k) => k >= 0);
  const done = () => { if (traces.length) Plotly.addTraces(div, traces); Plotly.relayout(div, { annotations: fipNotes.concat(notes) }); };
  if (old.length) Plotly.deleteTraces(div, old).then(done); else done();
}
let foldShownObj = null; // the fold shown, as shown (sinusoid or Keplerian)

function currentFold() {
  return foldShownObj;
}

// a time span in minutes, hours or days
function span(days) {
  const a = Math.abs(days);
  return a < 2 / 24 ? `${(days * 1440).toFixed(1)} min` : a < 2 ? `${(days * 24).toFixed(2)} h` : `${days.toFixed(3)} d`;
}
// a fold on a transit ephemeris: the ephemeris at the velocities, and the
//   conjunction the velocities put on their own
function transitNote(f) {
  const e = f.transit;
  if (!e) return '';
  const sig = e.shift_sigma !== null && e.shift_sigma !== undefined ? ` (${e.shift_sigma.toFixed(1)}\u03c3)` : '';
  return `<span class="trnote">\u25d0 ${esc(e.name)} (${esc(e.reference || '')}), ${esc(t('on_transit'))}: P = ${pm(e.P, e.P_err, 7)} d, `
    + `T\u2080 = ${e.t0.toFixed(5)} \u00b1 ${span(e.t0_err)} (${e.cycles} ${esc(t('cycles_from'))} ${e.tc.toFixed(5)}; `
    + `${esc(t('phase_word'))} \u00b1 ${e.phase_err.toFixed(4)})`
    + (f.kind === 'kepler' ? '' : `; K = ${pm(f.K, f.K_err, 2)} m/s ${esc(t('phase_fixed'))}; ${esc(t('conj_alone'))} `
      + `${e.shift >= 0 ? '+' : ''}${span(e.shift)} \u00b1 ${span(e.shift_err)}${sig}, K ${pm(e.K_free, e.K_free_err, 2)} m/s`) + '</span>';
}

const pm = (v, e, d) => `${v.toFixed(d)}${Number.isFinite(e) ? ` \u00b1 ${e.toFixed(d)}` : ''}`;

function showFold(id) {
  const base = ((quick && quick.result && quick.result.folds) || []).find((x) => x.id === +id);
  if (!base) return;
  foldShown = base.id;
  document.querySelectorAll('[data-fold]').forEach((b) => b.classList.toggle('on', +b.dataset.fold === base.id));
  document.querySelectorAll('[data-fmodel]').forEach((b) => b.classList.toggle('on', b.dataset.fmodel === foldModel));
  // fitted without the other ticked signals: asked of the server when
  //   what it was fitted without is not that
  if ((base.minus_tag || '') !== tickKey(othersOf(base.id))) {
    $('foldnote').innerHTML = `<span class="spin"></span> ${esc(t('folding'))}`;
    refitFold(base.id, foldModel).then(() => { if (foldShown === base.id) showFold(base.id); });
    return;
  }
  if (foldModel === 'kepler' && !base.kepler) { keplerOf(base); return; }
  const f = foldModel === 'kepler' ? { ...base.kepler, id: base.id, forced: base.forced, transit: base.transit } : base;
  foldShownObj = f;
  const pubNote = f.published ? ` \u00b7 <span class="pubnote">\u2605 ${esc(t('published'))} (${esc(f.published.reference)}): `
    + `P = ${f.published.period.toFixed(4)} d, K = ${f.published.K.toFixed(2)} m/s, e = ${f.published.e.toFixed(2)}</span>` : '';
  $('foldnote').innerHTML = f.kind === 'kepler'
    ? `<b>${foldLabel(base)}</b> ${esc(t('fmodel_kepler_short'))}: P = ${pm(f.period, f.P_err, 4)} d, K = ${pm(f.K, f.K_err, 2)} m/s, `
      + `e = ${pm(f.e, f.e_err, 2)}, \u03c9 = ${f.omega.toFixed(0)}\u00b0, rms ${f.rms.toFixed(2)} m/s \u00b7 ${esc(t('fold_note_kep'))}`
    : `<b>${foldLabel(base)}</b>: P = ${f.P_err ? pm(f.period, f.P_err, Math.min(8, Math.max(4, 1 - Math.floor(Math.log10(f.P_err))))) : f.period.toFixed(4)} d, K = ${f.K.toFixed(2)} \u00b1 ${f.K_err.toFixed(2)} m/s, `
      + `rms ${f.rms.toFixed(2)} m/s \u00b7 ${esc(t('fold_note'))}`;
  $('foldnote').innerHTML += pubNote + transitNote(f)
    // the minimum mass of a Keplerian orbit; the mass of a transiting one
    // m sin i: of the sinusoid (a circular orbit), the Keplerian or the
    //   transit's orbit, with the mass of the star of the card; the
    //   equilibrium temperature
    + (othersOf(base.id).length ? `<span class="massnote hint">${esc(t('fold_minus'))} ${esc(othersOf(base.id).map((tk) => plainLabel(foldById(tk.id))).join(', '))}</span>` : '')
    + bicText(f) + massText(f) + teqText(f.period, starNow);
  // the colour of the points: their instrument, their date, or their BERV
  //   (when the series has it), on one scale for all
  const hasBerv = f.instruments.some((inst) => inst.berv);
  const hasTime = f.instruments.some((inst) => inst.time && inst.time.length);
  $('fcol-berv').disabled = !hasBerv;
  $('fcol-date').disabled = !hasTime;
  const mode = (foldColour === 'berv' && !hasBerv) || (foldColour === 'date' && !hasTime) ? 'inst' : foldColour;
  document.querySelectorAll('[data-fcol]').forEach((b) => b.classList.toggle('on', b.dataset.fcol === mode));
  let scale = null;
  if (mode === 'date') {
    const ts = f.instruments.flatMap((inst) => inst.time || []);
    const tk = dateTicks(Math.min(...ts), Math.max(...ts), 640);
    scale = { cmin: Math.min(...ts), cmax: Math.max(...ts), colorscale: 'Viridis',
      colorbar: { title: { text: t('fcol_date') }, tickvals: tk.vals, ticktext: tk.text } };
  } else if (mode === 'berv') {
    const top = Math.max(1e-3, ...f.instruments.flatMap((inst) => (inst.berv || []).filter((v) => v !== null).map(Math.abs)));
    // Plotly's RdBu runs from blue to red: a positive BERV red, as in the PDF
    scale = { cmin: -top, cmax: top, colorscale: 'RdBu', reversescale: false, colorbar: { title: { text: 'BERV [km/s]' } } };
  }
  const order = (lastRV || []).map((inst) => inst.name);
  let barShown = false;
  const traces = f.instruments.map((inst) => {
    const i = Math.max(0, order.indexOf(inst.name));
    const vals = mode === 'date' ? inst.time : mode === 'berv' ? inst.berv : null;
    let marker = { color: COLOURS[i % 8], symbol: SYMBOLS[i % 8], size: 7, line: { color: '#08111f', width: 1 } };
    if (scale && vals) {
      marker = { ...marker, color: vals, colorscale: scale.colorscale, reversescale: !!scale.reversescale, cmin: scale.cmin, cmax: scale.cmax,
        showscale: !barShown, colorbar: { ...scale.colorbar, thickness: 12, len: 0.9, x: 1.01, outlinewidth: 0, tickfont: { size: 10 } } };
      barShown = true;
    } else if (scale) marker = { ...marker, color: '#5a6476' };
    const custom = inst.phase.map((_, k) => [inst.valid ? inst.valid[k] : null, inst.time ? dateText(inst.time[k]) : '',
      inst.berv ? inst.berv[k] : null]);
    const hover = `${inst.name}<br>phase %{x:.3f}<br>%{y:.2f} m/s`
      + (inst.time ? '<br>%{customdata[1]}' : '') + (inst.berv ? '<br>BERV %{customdata[2]:.2f} km/s' : '')
      + (inst.valid ? `<br>${esc(t('p_valid'))} %{customdata[0]:.2f}` : '');
    return { x: inst.phase, y: inst.rv, name: inst.name, type: 'scatter', mode: 'markers', customdata: custom,
      error_y: { type: 'data', array: inst.err, visible: true, thickness: 1, width: 0, color: scale ? 'rgba(200,220,255,0.35)' : COLOURS[i % 8] },
      marker, hovertemplate: `${hover}<extra></extra>` };
  });
  // the 1-sigma envelope of the fit, light grey, under the points
  if (f.curve.lo) {
    traces.unshift({ x: f.curve.phase, y: f.curve.lo, type: 'scatter', mode: 'lines', line: { width: 0 }, hoverinfo: 'skip', showlegend: false },
      { x: f.curve.phase, y: f.curve.hi, type: 'scatter', mode: 'lines', line: { width: 0 }, fill: 'tonexty',
        fillcolor: 'rgba(215,222,235,0.22)', name: t('envelope'), hoverinfo: 'skip' });
  }
  traces.push({ x: f.curve.phase, y: f.curve.rv, name: `K = ${f.K.toFixed(2)} m/s`, type: 'scatter', mode: 'lines',
    line: { color: '#e8eef8', width: 2 }, hoverinfo: 'skip' });
  // a known planet: its published orbit, on the fold's phases
  if (f.published) {
    traces.push({ x: f.published.phase, y: f.published.rv, type: 'scatter', mode: 'lines', hoverinfo: 'skip',
      name: `${f.published.name} (${f.published.reference || t('published')})`, line: { color: '#f5a524', width: 2, dash: 'dash' } });
  }
  // a white circle around a night with less than an even chance of being
  //   valid (an outlier, as the FIP saw it)
  const low = { x: [], y: [], p: [] };
  for (const inst of f.instruments) {
    (inst.valid || []).forEach((p, k) => { if (p < 0.5) { low.x.push(inst.phase[k]); low.y.push(inst.rv[k]); low.p.push(p); } });
  }
  if (low.x.length) {
    traces.push({ x: low.x, y: low.y, customdata: low.p, name: t('p_valid_low'), type: 'scatter', mode: 'markers',
      marker: { symbol: 'circle-open', size: 17, color: '#ffffff', line: { width: 1.6 } },
      hovertemplate: `${esc(t('p_valid'))} %{customdata:.2f}<extra></extra>` });
  }
  traces.push(offscaleTrace());
  const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
  $('foldplot').classList.add('on');
  Plotly.newPlot('foldplot', traces, {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)', font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
    margin: { ...PLOT_MARGIN, t: 30, r: scale ? 100 : PLOT_MARGIN.r }, legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom' },
    xaxis: { ...axis, title: t('phase_axis'), range: [0, 1] },
    // the velocities on the range of the series above
    yaxis: { ...axis, title: 'RV [m/s]', ...(view && view.y ? { range: view.y.slice() } : {}) },
  }, { responsive: true, displaylogo: false }).then(foldOffscale);
  // the solution on the series: the instruments about the fold's when it
  //   is shown there
  if ($('foldoverlay').checked) redrawSeries(); else showModel();
}

// the fold follows the y range of the series
function foldY() {
  if (view && view.y && window.Plotly && $('foldplot').data) {
    Plotly.relayout('foldplot', { 'yaxis.range': view.y.slice() });
    foldOffscale();
  }
}

// a red triangle at the edge of a plot, pointing to a point beyond its y
//   range (hover: its velocity)
function offscaleTrace() {
  return { x: [], y: [], type: 'scatter', mode: 'markers', uid: 'offscale', showlegend: false, cliponaxis: false,
    marker: { color: RED, size: 10, symbol: [], line: { width: 0 } },
    hovertemplate: `${esc(t('offscale'))}: %{customdata:.1f} m/s<extra></extra>` };
}
function setOffscale(id, xs, ys, yr, xr) {
  const div = $(id);
  if (!window.Plotly || !div.data) return;
  const k = div.data.findIndex((d) => d.uid === 'offscale');
  if (k < 0) return;
  const out = { x: [], y: [], symbol: [], v: [] };
  if (yr) {
    const pad = 0.02 * (yr[1] - yr[0]);
    xs.forEach((x, j) => {
      const y = ys[j];
      if (y === null || y === undefined || (xr && (x < xr[0] || x > xr[1]))) return;
      if (y > yr[1]) { out.x.push(x); out.y.push(yr[1] - pad); out.symbol.push('triangle-up'); out.v.push(y); }
      if (y < yr[0]) { out.x.push(x); out.y.push(yr[0] + pad); out.symbol.push('triangle-down'); out.v.push(y); }
    });
  }
  Plotly.restyle(div, { x: [out.x], y: [out.y], customdata: [out.v], 'marker.symbol': [out.symbol] }, [k]);
}
function rvOffscale() {
  if (!view) return;
  const insts = kept();
  setOffscale('rvplot', insts.flatMap((inst) => inst.time), insts.flatMap((inst) => inst.rv), view.y, view.x);
}
function foldOffscale() {
  const f = currentFold();
  if (!f || !view || !view.y) return;
  setOffscale('foldplot', f.instruments.flatMap((inst) => inst.phase), f.instruments.flatMap((inst) => inst.rv), view.y, null);
}

// the velocity of a Keplerian orbit (radvel's convention)
function keplerRV(tt, m, period) {
  const M = 2 * Math.PI * (tt - m.tp) / period;
  if (!m.e) return m.K * Math.cos(M + m.omega);
  let E = M + m.e * Math.sin(M);
  for (let k = 0; k < 50; k++) {
    const d = (E - m.e * Math.sin(E) - M) / (1 - m.e * Math.cos(E));
    E -= d;
    if (Math.abs(d) < 1e-12) break;
  }
  const nu = 2 * Math.atan2(Math.sqrt(1 + m.e) * Math.sin(E / 2), Math.sqrt(1 - m.e) * Math.cos(E / 2));
  return m.K * (Math.cos(nu + m.omega) + m.e * Math.cos(m.omega));
}
let modelTimer = null;
function scheduleModel() {
  clearTimeout(modelTimer);
  if ($('foldoverlay').checked || seriesAligned) modelTimer = setTimeout(showModel, 150);
}
// a fold's solution without the offsets of the instruments, the same for
//   all of them: its trend, and (sig) its sinusoid or Keplerian orbit
function curveAt(m, period, tt, sig = true) {
  const span = (tt - (m.ttrend !== undefined ? m.ttrend : m.tref)) / (m.tscale || 365.25);
  let v = 0;
  if (sig && m.kind === 'kepler') v += keplerRV(tt, m, period);
  else if (sig) {
    const arg = 2 * Math.PI * (tt - m.tref) / period;
    v += m.cos * Math.cos(arg) + m.sin * Math.sin(arg);
  }
  (m.trend || []).forEach((c, k) => { v += c * span ** (k + 1); });
  return v;
}
// the fit of the quick look (or a draw of it): its trend, and (sig) the
//   signals of its second pass
function fitAt(tm, tt, sig = true) {
  const sp = (tt - tm.tref) / tm.tscale;
  let v = (tm.coefs || []).reduce((acc, c, d) => acc + c * sp ** (d + 1), 0);
  if (sig) v += (tm.signals || []).reduce((acc, s) => acc + keplerRV(tt, s, s.period), 0);
  return v;
}
// the solution the series is drawn about, the same for every instrument:
//   the fold shown on the series (its trend and its signal), else the fit
//   of the quick look (its trend and its signals); its draws (the full
//   covariance of its fit) for the 1-sigma envelope
function seriesCurve() {
  const f = $('foldoverlay').checked ? currentFold() : null;
  // the signals ticked: their sum, on the trend of the fold shown when it
  //   is one of them (else of the first ticked)
  const picks = $('foldoverlay').checked ? ticked.map(tickFold).filter((x) => x && x.model) : [];
  if (picks.length) {
    const ref = picks.find((x) => f && x.id === f.id) || picks[0];
    const per = (x, m) => m.period || x.model.period || x.period;
    const sigOf = (x, m, tt) => curveAt(m, per(x, m), tt, true) - curveAt(m, per(x, m), tt, false);
    const nd = Math.min(...picks.map((x) => (x.draws || []).length));
    const draws = [];
    for (let j = 0; j < nd; j++) {
      draws.push((tt, sig = true) => curveAt(ref.draws[j], per(ref, ref.draws[j]), tt, false)
        + (sig ? picks.reduce((acc, x) => acc + sigOf(x, x.draws[j], tt), 0) : 0));
    }
    return { name: picks.map((x) => `${plainLabel(foldById(x.id))} \u00b7 ${x.period.toFixed(4)} d`).join(' + '),
      shortest: Math.min(...picks.map((x) => x.period)),
      at: (tt, sig = true) => curveAt(ref.model, per(ref, ref.model), tt, false) + (sig ? picks.reduce((acc, x) => acc + sigOf(x, x.model, tt), 0) : 0),
      draws, pub: null };
  }
  if (f && f.model) {
    const m = f.model, p = m.period || f.period;
    const pub = f.published && f.published.model;
    return { name: `#${f.id} · ${f.period.toFixed(4)} d`, shortest: f.period,
      at: (tt, sig = true) => curveAt(m, p, tt, sig),
      draws: (f.draws || []).map((d) => (tt, sig = true) => curveAt(d, d.period || p, tt, sig)),
      pub: pub ? { name: `★ ${f.published.name} (${t('published')})`, at: (tt) => curveAt(pub, pub.period, tt) } : null };
  }
  const tm = quick && quick.result && quick.result.trend_model;
  if (!tm || !((tm.coefs || []).length || (tm.signals || []).length)) return null;
  return { name: t('fitted_model'), shortest: (tm.signals || []).length ? Math.min(...tm.signals.map((s) => s.period)) : Infinity,
    at: (tt, sig = true) => fitAt(tm, tt, sig),
    draws: (tm.draws || []).map((d) => (tt, sig = true) => fitAt({ ...tm, coefs: d.coefs, signals: d.signals }, tt, sig)),
    pub: null };
}

// the solution on the series: one curve for every instrument (each about
//   it), its 1-sigma envelope; over the time shown, its signal where it is
//   resolved (a few tens of cycles at most: beyond, a band that would hide
//   the points), else its trend alone, dashed
function showModel() {
  const div = $('rvplot');
  if (!window.Plotly || !div.data) return;
  const old = div.data.map((d, k) => (String(d.uid || '').startsWith('model') ? k : -1)).filter((k) => k >= 0);
  if (old.length) Plotly.deleteTraces(div, old);
  const curve = seriesAligned ? seriesCurve() : null;
  if (!curve || !lastRV || !lastRV.length) return;
  const times = lastRV.flatMap((inst) => inst.time);
  let lo = Math.min(...times), hi = Math.max(...times);
  const pad = 0.01 * (hi - lo + 1);
  lo -= pad; hi += pad;
  const xr = view && view.x ? view.x : null;
  if (xr) { lo = Math.max(lo, xr[0]); hi = Math.min(hi, xr[1]); }
  if (!(hi > lo)) return;
  const sig = (hi - lo) / curve.shortest <= MODEL_CYCLES;
  if (!sig && Number.isFinite(curve.shortest)) {
    // too many cycles to draw each: the range of the solution in each
    //   step of time (what its ups and downs would fill), its trend dashed
    const nb = 720, ns = 32, bx = [], blo = [], bhi = [], bmid = [];
    const dt = (hi - lo) / nb;
    // each step looked at over one cycle of the shortest period at least:
    //   the range the solution fills there, not a part of a cycle
    const wide = Math.max(dt, curve.shortest);
    for (let k = 0; k < nb; k++) {
      const mid = lo + dt * (k + 0.5);
      let mn = Infinity, mx = -Infinity;
      for (let j = 0; j < ns; j++) {
        // times spread without a regular spacing (which a period could match)
        const v = curve.at(mid + wide * (((0.5 + j * 0.6180339887) % 1) - 0.5), true);
        if (v < mn) mn = v;
        if (v > mx) mx = v;
      }
      bx.push(mid); blo.push(mn); bhi.push(mx); bmid.push(curve.at(mid, false));
    }
    Plotly.addTraces(div, [
      { x: bx, y: blo, type: 'scatter', mode: 'lines', uid: 'model-lo', line: { width: 0 }, showlegend: false, hoverinfo: 'skip' },
      { x: bx, y: bhi, type: 'scatter', mode: 'lines', uid: 'model-hi', line: { width: 0 }, fill: 'tonexty', fillcolor: 'rgba(232,238,248,0.20)',
        name: `${curve.name} (${t('band_note')})`, hoverinfo: 'skip' },
      { x: bx, y: bmid, type: 'scatter', mode: 'lines', uid: 'model-best', showlegend: false, line: { color: '#e8eef8', width: 1.2, dash: 'dash' },
        opacity: 0.85, hoverinfo: 'skip' }]);
    return;
  }
  const n = sig ? Math.round(Math.min(2500, Math.max(300, 30 * (hi - lo) / curve.shortest))) : 300;
  const draws = curve.draws.slice(0, 120);
  const x = [], y = [], ylo = [], yhi = [], yp = [];
  for (let k = 0; k < n; k++) {
    const tt = lo + (hi - lo) * k / (n - 1);
    const best = curve.at(tt, sig);
    x.push(tt); y.push(best);
    // the spread of the draws about the best model, not their median (the
    //   draws of a sharp eccentric peak have their peaks at different times)
    if (draws.length > 5) {
      const vals = draws.map((d) => d(tt, sig)).sort((a, b) => a - b);
      const mid = percentile(vals, 50);
      ylo.push(best - (mid - percentile(vals, 15.87)));
      yhi.push(best + (percentile(vals, 84.13) - mid));
    }
    if (curve.pub && sig) yp.push(curve.pub.at(tt));
  }
  const traces = [];
  if (ylo.length) {
    traces.push({ x, y: ylo, type: 'scatter', mode: 'lines', uid: 'model-lo', line: { width: 0 }, showlegend: false, hoverinfo: 'skip' },
      { x, y: yhi, type: 'scatter', mode: 'lines', uid: 'model-hi', line: { width: 0 }, fill: 'tonexty',
        fillcolor: 'rgba(215,222,235,0.22)', name: t('envelope'), showlegend: false, hoverinfo: 'skip' });
  }
  traces.push({ x, y, type: 'scatter', mode: 'lines', uid: 'model-best', name: sig ? curve.name : `${curve.name} (${t('trend_only')})`,
    line: { color: '#e8eef8', width: 1.4, dash: sig ? 'solid' : 'dash' }, opacity: 0.85, hoverinfo: 'skip' });
  if (yp.length) {
    traces.push({ x, y: yp, type: 'scatter', mode: 'lines', uid: 'model-pub', name: curve.pub.name,
      line: { color: '#f5a524', width: 1.2, dash: 'dash' }, hoverinfo: 'skip' });
  }
  Plotly.addTraces(div, traces);
}

// the period slider, on a log scale
const pSlider = (p) => Math.round(STEPS * Math.min(1, Math.max(0, Math.log(p / pview.dom[0]) / Math.log(pview.dom[1] / pview.dom[0]))));
const pFromSlider = (k) => pview.dom[0] * (pview.dom[1] / pview.dom[0]) ** (k / STEPS);
function syncPeriods() {
  if (!pview) return;
  $('plo').value = pSlider(pview.p[0]); $('phi').value = pSlider(pview.p[1]);
  const lo = +$('plo').value / STEPS, hi = +$('phi').value / STEPS;
  $('psel').style.left = `${100 * lo}%`; $('psel').style.width = `${100 * (hi - lo)}%`;
  $('pdual').style.marginLeft = `${PLOT_MARGIN.l}px`; $('pdual').style.marginRight = `${PLOT_MARGIN.r}px`;
}
function periodsFromSliders() {
  if (!pview) return;
  let lo = +$('plo').value, hi = +$('phi').value;
  if (hi - lo < 5) { if (document.activeElement && document.activeElement.id === 'plo') lo = hi - 5; else hi = lo + 5; }
  pview.p = [pFromSlider(lo), pFromSlider(hi)];
  Plotly.relayout('fipplot', { 'xaxis.range': [Math.log10(pview.p[0]), Math.log10(pview.p[1])] });
  syncPeriods();
}

// what is shown changed since the FIP was made
function checkStale() {
  if (quick && quick.status === 'done' && quick.of !== JSON.stringify(shownOptions())) $('fipstale').textContent = t('stale');
  else $('fipstale').textContent = '';
}

async function quicklookPdf() {
  const opts = shownOptions();
  try {
    const resp = await fetch('/api/quicklook_pdf', { method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ options: opts, x: view ? view.x : null, y: view ? view.y : null, p: pview ? pview.p : null,
        quick: quick && quick.result ? quick.id : '', command: $('cmd-detailed').textContent,
        fold_colour: foldColour, overlay: $('foldoverlay').checked && foldShown !== null ? foldShown : null,
        fold_model: foldModel, series_colour: seriesColour, fip_view: fipView,
        mstar: +$('mstar').value || null, mstar_err: +$('mstar_err').value || 0, series_zero: zeroMode }) });
    if (!resp.ok) throw new Error((await resp.json()).error || resp.statusText);
    const link = document.createElement('a');
    link.href = URL.createObjectURL(await resp.blob());
    link.download = `koloa_quicklook_${folderName(opts.target || 'star') || 'star'}.pdf`;
    link.click();
    URL.revokeObjectURL(link.href);
  } catch (err) {
    alert(err.message);
  }
}

// the analysis kit: the script of the analysis with the settings of the
//   report's card, the velocities and the star, as a .tar.gz
async function analysisScript() {
  const opts = options('detailed');
  const btn = $('script-detailed');
  btn.disabled = true;
  try {
    const resp = await fetch('/api/analysis_script', { method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ options: opts, lang, mstar: +$('mstar').value || null, mstar_err: +$('mstar_err').value || 0 }) });
    if (!resp.ok) throw new Error((await resp.json()).error || resp.statusText);
    const link = document.createElement('a');
    link.href = URL.createObjectURL(await resp.blob());
    link.download = `koloa_${folderName(opts.target || '') || 'series'}_analysis.tar.gz`;
    link.click();
    URL.revokeObjectURL(link.href);
  } catch (err) {
    alert(err.message);
  } finally {
    btn.disabled = false;
  }
}

// a new target: every field back to the page's own default, no file, no
//   star, no plot, no run that ended, as fresh as a new session (a run that
//   still runs stays: closing a page does not stop it either)
function resetPage() {
  document.querySelectorAll('#tab-analysis input, #tab-analysis select').forEach((el) => {
    if (el.type === 'checkbox') el.checked = el.defaultChecked;
    else if (el.tagName === 'SELECT') {
      const def = [...el.options].find((opt) => opt.defaultSelected) || el.options[0];
      el.value = def ? def.value : '';
    } else el.value = el.defaultValue;
  });
  fileRows = [{ path: '', label: '' }];
  renderFiles();
  $('ident').innerHTML = '';
  $('targetsrc').textContent = '';
  starNow = {};
  daceSeen.clear();
  $('rvnote').textContent = '';
  $('rvtable').innerHTML = '';
  if (window.Plotly) Plotly.purge($('rvplot'));
  $('rvplot').classList.remove('on');
  $('plotcard').classList.remove('on');
  $('fipcard').classList.remove('on');
  if (window.Plotly) { Plotly.purge($('fipplot')); Plotly.purge($('foldplot')); }
  $('foldbuttons').innerHTML = ''; $('foldnote').textContent = '';
  $('tsbuttons').innerHTML = ''; $('tsnote').textContent = ''; tsShown = null; tsOf = null;
  if (window.Plotly) Plotly.purge($('tsplot'));
  $('tsplot').classList.remove('on');
  $('fipeach').innerHTML = ''; $('fipstatus').innerHTML = '';
  $('fipstale').textContent = ''; $('remstate').textContent = ''; $('remember').disabled = true;
  syncMirrors();
  vizierSources();
  view = null;
  quick = null;
  pview = null;
  lastRV = null;
  lastRes = null;
  rvRight = PLOT_MARGIN.r;
  onDisk = null;
  foldShown = null;
  foldShownObj = null;
  subtractList = [];
  clearTicks();
  delete $('mstar').dataset.typed;
  $('mstarsrc').textContent = '';
  showDiskPoints();
}

function newTarget() {
  resetPage();
  showTab('analysis');
  // the runs that ended, the quick looks, the archive in memory: forgotten
  //   here and by the server (a run still running is kept: Stop stops it)
  for (const [id, job] of [...jobs]) if (job.status !== 'running') { jobs.delete(id); closedLogs.delete(id); }
  renderJobs();
  api('/api/forget', {}).catch(() => {});
  history.replaceState(null, '', location.pathname);
  updateCommands();
  checkArchives();
  $('target').focus();
}

// -----------------------------------------------------------------------------
// the results remembered: a quick look kept, listed in its tab, recalled as
//   it was
// -----------------------------------------------------------------------------
function showTab(name) {
  document.querySelectorAll('.tabs .tab').forEach((b) => b.classList.toggle('on', b.dataset.tab === name));
  $('tab-analysis').hidden = name !== 'analysis';
  $('tab-remembered').hidden = name !== 'remembered';
  $('tab-batch').hidden = name !== 'batch';
  if (name === 'remembered') loadRemembered();
  else if (view) setTimeout(() => { syncSliders(); syncPeriods(); }, 50);
}

// the page as it is: what a recall puts back
function pageState() {
  return { target: $('target').value.trim(), files: filesNow(), root: $('root').value.trim(), outdir: $('outdir').value.trim(),
    detailed: readOptions('detailed'), clip: $('clip').checked, subtract: subtractList,
    // the signals ticked, and the FIP of what is left after each (its curve and its peaks)
    ticks: ticked.slice(), stages: stageLists().map((list) => stages.get(tickKey(list))).filter((st) => st && st.job && st.job.result)
      .map((st) => ({ key: st.key, ids: st.ids, label: st.label, period: st.job.result.period, family: st.job.result.family,
        peaks: st.job.result.peaks, peak_list: st.job.result.peak_list })),
    view: view ? { x: view.x.slice(), y: view.y.slice() } : null, periods: pview ? pview.p.slice() : null };
}

async function rememberResult() {
  if (!quick || quick.status !== 'done') return;
  if (quick.of !== JSON.stringify(shownOptions())) { $('remstate').textContent = t('stale'); return; }
  $('remember').disabled = true;
  $('remstate').innerHTML = `<span class="spin"></span> ${esc(t('remembering'))}`;
  try {
    const res = await api('/api/remember', { page: pageState(), quick: quick.id, note: $('remnote').value.trim() });
    $('remstate').innerHTML = `<span class="ok">\u2713</span> ${esc(t('remembered_ok'))} (${esc(res.created)})`;
    loadRemembered();
  } catch (err) {
    $('remstate').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
    $('remember').disabled = false;
  }
}

let rememberedList = [];
async function loadRemembered() {
  try { rememberedList = await api('/api/remembered'); } catch (err) { return; }
  $('nrem').textContent = rememberedList.length ? `(${rememberedList.length})` : '';
  renderRemembered();
}

function renderRemembered() {
  const box = $('remlist');
  if (!rememberedList.length) { box.innerHTML = `<p class="hint">${esc(t('no_remembered'))}</p>`; return; }
  const exp = fipExp;
  box.innerHTML = `<div class="remwrap"><table class="mini"><tr><th>${esc(t('col_target'))}</th><th>${esc(t('col_when'))}</th>`
    + `<th>${esc(t('col_nights'))}</th><th>${esc(t('col_peaks'))}</th><th>dv/dt [m/s/yr]</th><th>${esc(t('col_known'))}</th>`
    + `<th>${esc(t('col_each'))}</th><th>${esc(t('col_data'))}</th><th></th></tr>`
    + rememberedList.map((e) => {
      const s = e.summary || {};
      const nights = `<b>${s.n || 0}</b><br>` + Object.entries(s.instruments || {}).map(([k, v]) => `${esc(k)} ${v}`).join('<br>');
      const peaks = (s.peaks || []).map((pk) => `#${pk.id} ${pk.period.toFixed(4)} (${exp(pk.family)})`).join('<br>') || esc(t('none'));
      const known = (s.known || []).map((pl) => `${esc(pl.name)} ${(+pl.P).toPrecision(5)}`).join('<br>') || esc(t('none'));
      const acc = s.acceleration && s.acceleration.accel
        ? `${accelValue(s.acceleration.accel, '')}<br>(${s.acceleration.accel_sigma.toFixed(1)}\u03c3)` : esc(t('none'));
      const each = (s.each || []).map((one) => `${esc(one.name)}: ${one.skipped ? '&lt; 10 n'
        : one.best ? `${one.best.period.toFixed(3)} (${exp(one.best.family)})` : esc(t('none'))}`).join('<br>') || esc(t('none'));
      const data = (s.files || []).map(esc).concat(s.archives || []).join('<br>')
        + (s.exclude ? `<br><span class="hint">\u2212 ${esc(s.exclude)}</span>` : '');
      return `<tr><td><b>${esc(e.target || '?')}</b>${e.note ? `<span class="note">${esc(e.note)}</span>` : ''}</td>`
        + `<td class="mono">${esc(e.created || '')}</td><td class="mono">${nights}</td><td class="mono">${peaks}</td><td class="mono">${acc}</td>`
        + `<td class="mono">${known}</td><td class="mono">${each}</td><td class="mono">${data}</td>`
        + `<td class="acts"><button type="button" class="small go" data-recall="${esc(e.id)}">${esc(t('recall'))}</button>`
        + `<button type="button" class="small stop" data-unremember="${esc(e.id)}">${esc(t('forget'))}</button></td></tr>`;
    }).join('') + '</table></div>';
}

async function recallResult(id) {
  let res;
  try { res = await api('/api/recall', { id }); } catch (err) { alert(err.message); return; }
  await applyRecall(res);
}

// a result put back in the page as it was (a result remembered, or a file
//   of a batch): the page, the velocities, its quick FIP, nothing computed
async function applyRecall(res) {
  showTab('analysis');
  resetPage();
  // a result recalled keeps the instruments it had
  (res.rv && res.rv.instruments || []).forEach((inst) => daceSeen.add(inst.name));
  const page = res.page;
  $('target').value = page.target || '';
  $('root').value = page.root || 'archives';
  $('outdir').value = page.outdir || '';
  fileRows = (page.files && page.files.length ? page.files : [{ path: '', label: '' }]).map((r) => ({ path: r.path, label: r.label || '' }));
  renderFiles();
  for (const [key, val] of Object.entries(page.detailed || {})) {
    const el = document.querySelector(`[data-for="detailed"][data-opt="${key}"]`);
    if (!el) continue;
    if (el.type === 'checkbox') el.checked = !!val; else el.value = val;
  }
  syncMirrors();
  $('clip').checked = !!page.clip;
  subtractList = page.subtract || [];
  clearTicks();
  ticked = page.ticks || [];
  for (const st of page.stages || []) {
    stages.set(st.key, { key: st.key, ids: st.ids, label: st.label, error: null,
      job: { status: 'done', result: { period: st.period, family: st.family, peaks: st.peaks, peak_list: st.peak_list } } });
  }
  if (page.target) resolveStar();
  checkArchives();
  updateCommands();
  // the velocities, the FIP and the folds as they were: nothing recomputed
  $('plotcard').classList.add('on');
  await drawVelocities(res.rv, res.notes);
  if (view && page.view) { view.x = page.view.x; view.y = page.view.y; applyView(); }
  quick = { ...res.quick, of: JSON.stringify(shownOptions()) };
  $('fipcard').classList.add('on');
  drawFip(quick.result, quick.each || []);
  quick.drawn = `1/${(quick.each || []).length}`;
  $('fipstatus').innerHTML = fipSummary(quick.result, quick.elapsed);
  if ((zeroMode === 'fit' || $('foldoverlay').checked) && quick.result.trend_model) await redrawSeries();
  if (pview && page.periods) {
    pview.p = page.periods;
    Plotly.relayout('fipplot', { 'xaxis.range': [Math.log10(pview.p[0]), Math.log10(pview.p[1])] });
    syncPeriods();
  }
  $('remnote').value = res.entry.note || '';
  $('remstate').textContent = res.entry.id ? `${t('recalled_from')} ${res.entry.created}` : t('from_batch');
  $('remember').disabled = !!res.entry.id;
  $('plotcard').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// -----------------------------------------------------------------------------
// the batch FIP: the quick FIP of many files, a table, each opened in the page
// -----------------------------------------------------------------------------
let batchFiles = [];
let batch = null;            // the batch running or run: its state
let batchSort = { key: 'fip', dir: 1 };

function renderBatchFiles() {
  const box = $('batchfiles');
  if (!batchFiles.length) { box.innerHTML = `<p class="hint">${esc(t('batch_none'))}</p>`; return; }
  box.innerHTML = `<p class="hint">${batchFiles.length} ${esc(t('batch_files'))}</p><div class="batchlist">`
    + batchFiles.map((path, i) => `<span class="bfile" title="${esc(path)}">${esc(path.split('/').pop())}`
      + `<button type="button" class="small" data-bdel="${i}">\u00d7</button></span>`).join('') + '</div>';
}
function addBatchFiles(paths) {
  for (const path of paths || []) if (!batchFiles.includes(path)) batchFiles.push(path);
  renderBatchFiles();
}

async function runBatch() {
  if (!batchFiles.length) return;
  try {
    batch = await api('/api/batch', { paths: batchFiles, options: readOptions('detailed'), archives: $('batcharchives').checked,
      regather: $('batchregather').checked, root: $('root').value.trim() });
    pollBatch(batch.id);
  } catch (err) {
    $('batchstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
}
async function pollBatch(id) {
  if (!batch || batch.id !== id) return;
  try { batch = await api(`/api/batch?id=${id}`); } catch (err) { return; }
  renderBatch();
  if (batch.status === 'running') setTimeout(() => pollBatch(id), 1500);
}

// a FIP: green below 1 %
const fipCell = (v) => (v === null || v === undefined ? ''
  : `<span class="${v < 0.01 ? 'ok' : ''}">${fipExp(v)}</span>`);
function batchRows() {
  const rows = batch.items.map((item, i) => ({ ...item, i, ...(item.summary || {}), object: objectName(item),
    // a plausible transit first, then by its signal-to-noise ratio
    transit_snr: item.summary && item.summary.transit && item.summary.transit.snr !== undefined && item.summary.transit.snr !== null
      ? -(item.summary.transit.snr + (item.summary.transit.plausible ? 1000 : 0)) : null }));
  const { key, dir } = batchSort;
  const val = (r) => (r[key] === null || r[key] === undefined ? Infinity * dir : r[key]);
  if (key === 'name' || key === 'object') return rows.sort((a, b) => dir * String(a[key]).localeCompare(String(b[key])));
  return rows.sort((a, b) => dir * (val(a) - val(b)));
}
// the object of a file of a batch: its APERO name, else the name read
const objectName = (item) => (item.star ? item.star.apero || item.star.raw || '' : '');
function objectCell(item) {
  const s = item.star;
  if (!s) return '';
  const from = s.source === 'OBJECT column' ? `OBJECT ${s.raw}` : `${t('file_name')} ${s.raw}`;
  const tip = s.apero ? `${from} \u2192 ${s.apero} (APERO${s.status ? `, ${s.status}` : ''}), ${s.target} ${t('for_simbad')}`
    : `${from} (${t('not_apero')})`;
  // its known planets (their periods) and candidate TOIs, in short; their
  //   names under the cursor
  const per = (p) => (p >= 100 ? p.toFixed(0) : p >= 10 ? p.toFixed(1) : p.toFixed(2));
  const pl = s.planets || [], ct = s.tois || [];
  const which = (pl.length ? `<span title="${esc(pl.map((p) => `${p.name}: ${per(p.P)} d`).join('\n'))}">\u2605 ${pl.map((p) => per(p.P)).join(', ')} d</span>` : '')
    + (pl.length && ct.length ? ' \u00b7 ' : '')
    + (ct.length ? `<span title="${esc(ct.map((c) => `TOI-${c.toi} (${c.disposition}): ${per(c.P)} d`).join('\n'))}">\u25d0 ${ct.map((c) => per(c.P)).join(', ')} d</span>` : '');
  return `<span title="${esc(tip)}">${esc(objectName(item))}</span>${s.spt ? ` <span class="hint">${esc(s.spt)}</span>` : ''}`
    + (which ? `<div class="hint planets">${which}</div>` : '');
}
// the best period of a line on a known planet's (\u2605) or a TOI's
//   (\u25d0): within the width of a peak (1 / the span, in frequency)
function periodMark(r) {
  const s = r.star || {};
  if (!r.period) return '';
  const near = (P) => Math.abs(1 / r.period - 1 / P) <= 1 / Math.max(r.baseline || 1, 1);
  if ((s.planets || []).some((p) => near(p.P))) return ' \u2605';
  if ((s.tois || []).some((c) => near(c.P))) return ' \u25d0';
  return '';
}
// the transit in TESS of a line (its best peak, when its FIP is below 1 %):
//   flagged when plausible, its depth, radius and signal-to-noise ratio
function transitCell(tr) {
  if (!tr) return '';
  const tip = esc(tr.why || '');
  if (tr.status === 'plausible') {
    return `<span class="ok" title="${tip}">\u2691 ${tr.depth.toFixed(2)} ppt \u00b7 ${tr.radius.toFixed(1)} R\u2295 \u00b7 ${tr.snr.toFixed(1)}\u03c3</span>`;
  }
  if (tr.status === 'none') return `<span class="hint" title="${tip}">${tr.snr !== null && tr.snr !== undefined ? `${tr.snr.toFixed(1)}\u03c3` : '-'}</span>`;
  return `<span class="hint" title="${tip}">${esc(tr.status === 'no TESS' ? t('ts_no_tess') : tr.status === 'not observed' ? t('ts_not_observed') : tr.status)}</span>`;
}
// the instruments of a line (nights of each), their sources in the tip
function instCell(r) {
  if (!r.instruments) return '';
  const src = r.sources || {};
  return Object.entries(r.instruments).map(([name, n]) => `<span title="${esc(src[name] || '')}">${esc(name)}\u00a0${n}</span>`).join('');
}
function renderBatch() {
  if (!batch) { $('batchtable').innerHTML = ''; $('batchstatus').textContent = ''; return; }
  const done = batch.items.filter((item) => item.status !== 'waiting' && item.status !== 'running').length;
  const running = batch.status === 'running';
  $('batchstatus').innerHTML = `${running ? '<span class="hourglass">\u23f3</span> ' : ''}${done} / ${batch.items.length} \u00b7 ${clock(batch.elapsed)}`
    + (running ? ` <button type="button" class="small stop" id="batchstop">${esc(t('stop_fip'))}</button>` : '')
    + (done ? ` <button type="button" class="small" id="batchcsv">CSV</button>` : '');
  const cols = [['name', t('col_file')], ['object', t('col_object')]].concat(batch.archives ? [['n', t('col_insts')]] : [['n', t('col_nights')]]).concat([['baseline', t('col_span')], ['period', 'P [d]'],
    ['fip', t('col_fip')], ['fip_alone', t('col_fip_alone')], ['K', 'K [m/s]'], ['rms', 'rms [m/s]'], ['accel', 'dv/dt [m/s/yr]'], ['accel_sigma', '\u03c3'],
    ['transit_snr', t('col_transit')]]);
  const head = cols.map(([k, label]) => `<th data-bsort="${k}" class="sortable${batchSort.key === k ? ' sorted' : ''}">${esc(label)}`
    + `${batchSort.key === k ? (batchSort.dir > 0 ? ' \u25b2' : ' \u25bc') : ''}</th>`).join('');
  const num = (v, d) => (v === null || v === undefined ? '' : v.toFixed(d));
  $('batchtable').innerHTML = `<div class="batchwrap"><table class="mini batch"><tr>${head}<th></th></tr>`
    + batchRows().map((r) => {
      let state = '';
      if (r.status === 'running' && r.stage === 'gather') state = `<span class="hourglass">\u23f3</span> ${esc(t('batch_gather'))}`;
      else if (r.stage === 'tess') state = `<span class="hourglass">\u23f3</span> ${esc(t('batch_tess'))}`;
      else if (r.status === 'running') {
        const p = r.progress;
        state = `<span class="hourglass">\u23f3</span> ${esc(t({ noise: 'quick_noise', fip1: 'quick_fip1', planets: 'quick_planets', fip2: 'quick_fip2' }[r.step] || 'quick_noise'))}`
          + (p ? ` ${Math.round(100 * p.done / Math.max(p.total, 1))} %` : '');
      } else if (r.status === 'failed' || r.status === 'stopped') state = `<span class="bad">${esc(r.error || t(r.status))}</span>`;
      else if (r.status === 'waiting') state = `<span class="hint">${esc(t('batch_waiting'))}</span>`;
      const acc = r.accel !== null && r.accel !== undefined ? `${r.accel >= 0 ? '+' : '\u2212'}${Math.abs(r.accel).toPrecision(3)} \u00b1 ${r.accel_err.toPrecision(2)}` : '';
      if (r.note) state += ` <span class="hint" title="${esc(r.note)}">\u24d8</span>`;
      const nights = batch.archives ? `<td class="insts">${r.n ? `<b>${r.n}</b>` : ''}${instCell(r)}</td>` : `<td class="num">${r.n ?? ''}</td>`;
      return `<tr><td class="bname" title="${esc(r.path)}">${esc(r.name).replace(/_/g, '_<wbr>')}</td><td>${objectCell(r)}</td>${nights}<td class="num">${num(r.baseline, 0)}</td>`
        + `<td class="num">${num(r.period, 4)}${periodMark(r)}</td><td class="num">${fipCell(r.fip)}</td><td class="num">${fipCell(r.fip_alone)}</td>`
        + `<td class="num">${r.K !== null && r.K !== undefined ? `${r.K.toFixed(2)} \u00b1 ${r.K_err.toFixed(2)}` : ''}</td><td class="num">${num(r.rms, 2)}</td>`
        + `<td class="num">${acc}</td><td class="num">${num(r.accel_sigma, 1)}</td><td>${transitCell(r.transit)}</td>`
        + `<td>${state}${r.status === 'done' ? `<button type="button" class="small go" data-bopen="${r.i}">${esc(t('open'))}</button>` : ''}</td></tr>`;
    }).join('') + '</table></div>';
}
function batchCsv() {
  const head = ['file', 'path', 'object', 'apero_name', 'simbad_name', 'known_planets_d', 'candidate_tois_d', 'nights', 'instruments', 'baseline_d', 'period_d', 'fip', 'fip_alone', 'K_ms', 'K_err_ms', 'rms_ms', 'dvdt_msyr', 'dvdt_err_msyr', 'dvdt_sigma',
    'transit_plausible', 'transit_snr', 'transit_depth_ppt', 'transit_radius_earth', 'transit_note', 'status'];
  const tr = (r, key) => (r.transit ? r.transit[key] : null);
  const insts = (r) => (r.instruments ? Object.entries(r.instruments).map(([k, n]) => `${k} ${n}`).join('; ') : '');
  const plist = (r, key, nm) => (r.star && r.star[key] ? r.star[key].map((p) => `${nm(p)} ${p.P}`).join('; ') : '');
  const lines = [head.join(',')].concat(batchRows().map((r) => [r.name, r.path, r.star ? r.star.raw : '', r.star ? r.star.apero : '', r.star ? r.star.target : '',
    plist(r, 'planets', (p) => p.name), plist(r, 'tois', (p) => `TOI-${p.toi}`),
    r.n, insts(r), r.baseline, r.period, r.fip, r.fip_alone, r.K, r.K_err,
    r.rms, r.accel, r.accel_err, r.accel_sigma, tr(r, 'plausible'), tr(r, 'snr'), tr(r, 'depth'), tr(r, 'radius'), tr(r, 'why'),
    r.status].map((v) => (v === null || v === undefined ? '' : `"${String(v).replace(/"/g, '""')}"`)).join(',')));
  const link = document.createElement('a');
  link.href = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], { type: 'text/csv' }));
  link.download = 'koloa_batch_fip.csv';
  link.click();
  URL.revokeObjectURL(link.href);
}

// -----------------------------------------------------------------------------
// start
// -----------------------------------------------------------------------------
document.addEventListener('input', (e) => {
  if (e.target.dataset && e.target.dataset.row !== undefined) {
    const row = fileRows[+e.target.dataset.row];
    if (e.target.classList.contains('fpath')) row.path = e.target.value;
    if (e.target.classList.contains('flabel')) row.label = e.target.value;
  }
  updateCommands();
  if (e.target.id === 'target' || e.target.id === 'root') checkArchives();
  if (e.target.id === 'target') $('targetsrc').textContent = '';
  if (e.target.id === 'exclude') styleExcluded();
});
// a file's path typed: its star, when the page has none
document.addEventListener('change', (e) => {
  if (e.target.classList && e.target.classList.contains('fpath')) starFromFile(e.target.value);
});
document.addEventListener('change', (e) => {
  if (e.target.dataset && e.target.dataset.inst) setExcluded(e.target.dataset.inst, !e.target.checked);
  if (e.target.dataset && e.target.dataset.mirror) {
    const twin = reportBox(e.target.dataset.mirror);
    if (twin) { twin.checked = e.target.checked; twin.dispatchEvent(new Event('change', { bubbles: true })); }
    return;
  }
  if (e.target.dataset && e.target.dataset.for === 'detailed' && ARCHIVES.includes(e.target.dataset.opt)) syncMirrors();
  // an archive ticked or not: the plot shows what the report will use
  if (e.target.dataset && e.target.dataset.for === 'detailed' && ARCHIVES.includes(e.target.dataset.opt)
      && $('plotcard').classList.contains('on')) plotVelocities();
  // a period of the literature ticked: the SHO at it; none: back to the bands
  if (e.target.name === 'prot') {
    $('rotation').value = e.target.value;
    $('fipgp').value = e.target.value ? 'sho' : 'banded';
    updateCommands();
  }
});
document.addEventListener('change', updateCommands);
// the sources of VizieR follow its box: off with it
function vizierSources() {
  const on = $('gather-vizier').checked;
  $('vzsources').classList.toggle('off', !on);
  $('vzsources').querySelectorAll('input').forEach((el) => { el.disabled = !on; });
}
$('gather-vizier').addEventListener('change', vizierSources);
// what is shown, or how the trend is fitted, changed since the quick FIP
document.addEventListener('change', checkStale);
$('resolve').addEventListener('click', resolveStar);
$('target').addEventListener('keydown', (e) => { if (e.key === 'Enter') resolveStar(); });
$('plot').addEventListener('click', plotVelocities);
for (const ax of ['x', 'y']) {
  for (const end of ['lo', 'hi']) $(`${ax}${end}`).addEventListener('input', () => fromSliders(ax));
}
$('clip').addEventListener('change', refitY);
$('plo').addEventListener('input', periodsFromSliders);
$('phi').addEventListener('input', periodsFromSliders);
$('refip').addEventListener('click', startQuick);
// the series by instrument or by BERV: drawn again, its ranges kept
document.addEventListener('click', async (e) => {
  const b = e.target.closest('[data-scol]');
  if (!b || b.disabled) return;
  seriesColour = b.dataset.scol;
  try { localStorage.setItem('koloa-seriescolour', seriesColour); } catch (err) { /* no storage */ }
  document.querySelectorAll('[data-scol]').forEach((x) => x.classList.toggle('on', x.dataset.scol === seriesColour));
  redrawSeries();
});
// each instrument about its fitted offset, or about its median
document.addEventListener('click', (e) => {
  const b = e.target.closest('[data-zero]');
  if (!b) return;
  zeroMode = b.dataset.zero;
  try { localStorage.setItem('koloa-zero', zeroMode); } catch (err) { /* no storage */ }
  redrawSeries();
});
// the quick FIP asked for, or stopped (while it waits for its turn too)
document.addEventListener('click', (e) => {
  if (e.target.closest('.startfip')) startQuick();
  if (e.target.closest('.stopfip') && quick) api('/api/quickstop', { id: quick.id }).catch(() => {});
});
$('pdf').addEventListener('click', quicklookPdf);
$('script-detailed').addEventListener('click', analysisScript);
document.addEventListener('click', (e) => { const b = e.target.closest('[data-fold]'); if (b) showFold(b.dataset.fold); });
// the colour of the fold: by instrument, date or BERV
document.addEventListener('click', (e) => {
  const b = e.target.closest('[data-fcol]');
  if (!b || b.disabled) return;
  foldColour = b.dataset.fcol;
  try { localStorage.setItem('koloa-foldcolour', foldColour); } catch (err) { /* no storage */ }
  if (foldShown !== null) showFold(foldShown);
  else document.querySelectorAll('[data-fcol]').forEach((x) => x.classList.toggle('on', x.dataset.fcol === foldColour));
});
$('foldoverlay').addEventListener('change', redrawSeries);
// the mass of the star typed: kept (the resolver no longer fills it), and
//   the masses of the fold shown again
for (const id of ['mstar', 'mstar_err']) {
  $(id).addEventListener('input', () => {
    $('mstar').dataset.typed = '1';
    $('mstarsrc').textContent = t('typed');
    if (foldShown !== null) showFold(foldShown);
  });
}
// the fold's sinusoid or its Keplerian orbit
document.addEventListener('click', (e) => {
  const b = e.target.closest('[data-fmodel]');
  if (!b) return;
  foldModel = b.dataset.fmodel;
  try { localStorage.setItem('koloa-foldmodel', foldModel); } catch (err) { /* no storage */ }
  document.querySelectorAll('[data-fmodel]').forEach((x) => x.classList.toggle('on', x.dataset.fmodel === foldModel));
  if (foldShown !== null) showFold(foldShown);
});
$('foldgo').addEventListener('click', () => forceFold(+$('foldp').value, false));
// a known planet: a fold at its published period
document.addEventListener('click', (e) => { const b = e.target.closest('[data-known]'); if (b) forceFold(0, false, b.dataset.known); });
document.addEventListener('click', (e) => { const b = e.target.closest('[data-transit]'); if (b) forceFold(0, false, '', b.dataset.transit); });
document.addEventListener('click', (e) => {
  const b = e.target.closest('[data-fipview]');
  if (!b || b.disabled) return;
  fipView = b.dataset.fipview;
  try { localStorage.setItem('koloa-fipview', fipView); } catch (err) { /* no storage */ }
  showFipView();
});
$('foldp').addEventListener('keydown', (e) => { if (e.key === 'Enter') forceFold(+$('foldp').value, false); });
document.addEventListener('change', (e) => { const b = e.target.closest('[data-tick]'); if (b) toggleTick(b.dataset.tick); });
document.addEventListener('click', (e) => { if (e.target.id === 'unsubtract') { subtractList = []; startQuick(); } });
$('fullrange').addEventListener('click', () => { $('clip').checked = false; view = null; refitY(); });
window.addEventListener('resize', () => { if (view) setTimeout(() => { syncSliders(); dateAxis(); }, 100); });
$('run-gather').addEventListener('click', () => run('gather', options('gather')));
$('run-detailed').addEventListener('click', () => run('detailed', options('detailed')));
$('run-archive').addEventListener('click', () => run('archive', options('archive')));
$('newtarget').addEventListener('click', newTarget);
document.addEventListener('click', async (e) => {
  const picker = e.target.closest('[data-pick]');
  if (picker) {
    const field = picker.dataset.row !== undefined
      ? document.querySelector(`input.fpath[data-row="${picker.dataset.row}"]`) : $(picker.dataset.into);
    const label = picker.textContent;
    picker.disabled = true; picker.textContent = t('picking');
    try {
      const res = await api(`/api/pick?${new URLSearchParams({ kind: picker.dataset.pick, start: field.value || '' })}`);
      if (res.path) {
        field.value = res.path;
        if (picker.dataset.row !== undefined) fileRows[+picker.dataset.row].path = res.path;
        updateCommands();
        if (picker.dataset.row !== undefined) starFromFile(res.path);
        if (picker.dataset.into === 'root') checkArchives();
      }
    } catch (err) {
      alert(err.message);
    } finally {
      picker.disabled = false; picker.textContent = label;
    }
  }
  if (e.target.id === 'addfile') {
    fileRows.push({ path: '', label: '' });
    renderFiles();
  }
  const del = e.target.closest('[data-del]');
  if (del) {
    fileRows.splice(+del.dataset.del, 1);
    if (!fileRows.length) fileRows.push({ path: '', label: '' });
    renderFiles();
    updateCommands();
  }
  const toiChip = e.target.closest('[data-toi]');
  if (toiChip) {
    $('toi_on').checked = true;
    const have = $('tois').value.split(/[\s,]+/).filter(Boolean);
    if (!have.includes(toiChip.dataset.toi)) have.push(toiChip.dataset.toi);
    $('tois').value = have.join(' ');
    updateCommands();
  }
  if (e.target.id === 'refresh-star') resolveStar(true);
  const tab = e.target.closest('[data-tab]');
  if (tab) showTab(tab.dataset.tab);
  const rec = e.target.closest('[data-recall]');
  if (rec) recallResult(rec.dataset.recall);
  const unrem = e.target.closest('[data-unremember]');
  if (unrem && confirm(t('forget_confirm'))) {
    try { await api('/api/unremember', { id: unrem.dataset.unremember }); } catch (err) { alert(err.message); }
    loadRemembered();
  }
});
$('remember').addEventListener('click', rememberResult);
// the batch FIP
$('batchpick').addEventListener('click', async () => {
  try {
    const res = await api(`/api/pick?${new URLSearchParams({ kind: 'files', start: batchFiles[0] || '' })}`);
    if (res.paths) addBatchFiles(res.paths);
  } catch (err) { alert(err.message); }
});
$('batchfolderadd').addEventListener('click', async () => {
  try {
    const res = await api(`/api/listfiles?${new URLSearchParams({ folder: $('batchfolder').value.trim(), pattern: $('batchpattern').value.trim() || '*.rdb' })}`);
    addBatchFiles(res.paths);
    if (!res.paths.length) alert(t('batch_nomatch'));
  } catch (err) { alert(err.message); }
});
$('batchclear').addEventListener('click', () => { batchFiles = []; renderBatchFiles(); });
$('batchrun').addEventListener('click', runBatch);
// the copy of APERO's names: what it is, and a new one
async function aperoState() {
  try {
    const s = await api('/api/apero_names');
    $('aperostate').textContent = s.objects ? `${t('apero_names')}: ${s.objects} ${t('objects')}${s.tarball ? `, ${s.tarball}` : ''}${s.fetched ? ` (${s.fetched})` : ''}` : t('apero_none');
  } catch (err) { $('aperostate').textContent = ''; }
}
$('aperorefresh').addEventListener('click', async () => {
  const btn = $('aperorefresh');
  btn.disabled = true;
  $('aperostate').innerHTML = `<span class="spin"></span> ${esc(t('apero_fetching'))}`;
  try { await api('/api/apero_refresh', {}); } catch (err) { alert(err.message); }
  btn.disabled = false;
  aperoState();
});
document.addEventListener('click', async (e) => {
  const del = e.target.closest('[data-bdel]');
  if (del) { batchFiles.splice(+del.dataset.bdel, 1); renderBatchFiles(); }
  const head = e.target.closest('[data-bsort]');
  if (head) {
    const key = head.dataset.bsort;
    batchSort = { key, dir: batchSort.key === key ? -batchSort.dir : 1 };
    renderBatch();
  }
  if (e.target.id === 'batchstop' && batch) api('/api/batchstop', { id: batch.id }).catch(() => {});
  const ts = e.target.closest('[data-tsearch]');
  if (ts) showTransit(ts.dataset.tsearch, true);
  const tv = e.target.closest('[data-tsview]');
  if (tv) {
    tsView = tv.dataset.tsview;
    try { localStorage.setItem('koloa-tsview', tsView); } catch (err) { /* no storage */ }
    if (tsLast) drawTransit(tsLast.res, tsLast.c);
    else document.querySelectorAll('[data-tsview]').forEach((x) => x.classList.toggle('on', x.dataset.tsview === tsView));
  }
  if (e.target.id === 'tsfetch' && tsShown) showTransit(tsShown, true);
  if (e.target.id === 'batchcsv') batchCsv();
  const open = e.target.closest('[data-bopen]');
  if (open && batch) {
    try { await applyRecall(await api('/api/batch_open', { id: batch.id, index: +open.dataset.bopen })); } catch (err) { alert(err.message); }
  }
});

(async () => {
  $('detailed-options').innerHTML = renderOptions('detailed');
  renderBatchFiles();
  aperoState();
  langHooks.push(aperoState, showCwd, renderJobs, renderFiles, checkArchives, renderRemembered, dateAxis, renderBatchFiles, renderBatch,
    () => { if (quick && quick.each) drawEach(quick.each); });
  loadRemembered();
  jobDoneHooks.push(async (job) => {
    if (job.action === 'archive') refreshInfo();
    if (job.action !== 'gather') return;
    checkArchives();
    if (job.status === 'done') archivesByDefault(await diskState());
  });
  jobStartHooks.push((action) => { if (action === 'gather') checkArchives(); });
  applyLang();
  await refreshInfo();
  try { (await api('/api/jobs')).forEach(keep); } catch (err) { /* nothing yet */ }
  renderJobs();
  // the fields from the address (?target=GJ%20436&file=...&root=...&plot=1)
  const params = new URLSearchParams(location.search);
  for (const key of ['target', 'root', 'outdir', 'exclude']) {
    if (params.get(key)) $(key).value = params.get(key);
  }
  const given = params.getAll('file');
  if (given.length) fileRows = given.map((path) => ({ path, label: '' }));
  renderFiles();
  updateCommands();
  checkArchives();
  if (params.get('target')) resolveStar();
  if (params.get('plot')) plotVelocities();
  setInterval(poll, 1500);
})();

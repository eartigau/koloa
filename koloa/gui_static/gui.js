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
const ARCHIVES = ['dace', 'carmenes'];
function reportBox(key) { return document.querySelector(`input[data-for="detailed"][data-opt="${key}"]`); }
function syncMirrors() {
  document.querySelectorAll('[data-mirror]').forEach((box) => { const twin = reportBox(box.dataset.mirror); if (twin) box.checked = twin.checked; });
}
function showDiskPoints() {
  const pts = (onDisk && onDisk.exists && onDisk.points) || {};
  for (const key of ARCHIVES) {
    $(`pts-${key}`).textContent = pts[key] ? `(${pts[key]} ${t('points_word')})` : onDisk && onDisk.exists ? '(0)' : '';
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
    if ((disk.points || {})[key] && reportBox(key)) { reportBox(key).checked = true; ticked = true; }
  }
  syncMirrors();
  updateCommands();
  return ticked;
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
      + `${text} <span class="hint">${esc(t('use_sho'))}</span></label>`;
    const rot = (id.variability || []).map((v) => chip(v.period, `${esc(v.type)} ${v.period} d <span class="hint">(${esc(v.bibcode || '')})</span>`));
    if (id.carmenes && id.carmenes.p_rot) rot.push(chip(id.carmenes.p_rot, `ROT ${esc(id.carmenes.p_rot)} d <span class="hint">(CARMENES, ${esc(id.carmenes.p_rot_source || '')})</span>`));
    if (id.archive_rotation) rot.push(chip(id.archive_rotation, `ROT ${esc(id.archive_rotation)} d <span class="hint">(NASA Exoplanet Archive)</span>`));
    const num = (v, d) => (v == null || v === '' ? '?' : (+v).toFixed(d));
    const planets = (id.planets || []).map((pl) => `<div class="planet"><b>${esc(pl.name)}</b> `
      + `<span class="pnum">P ${num(pl.P, 4)} d \u00b7 K ${num(pl.K, 2)} m/s`
      + (pl.mass_earth != null ? ` \u00b7 m sin i ${num(pl.mass_earth, 1)} M\u2295` : '') + `</span><br>`
      + `<span class="hint">${pl.reference_url ? `<a href="${esc(pl.reference_url)}" target="_blank">${esc(pl.reference || '')}</a>` : esc(pl.reference || '')}`
      + ` \u00b7 ${pl.solutions} ${esc(t(pl.solutions === 1 ? 'solution' : 'solutions'))}`
      + (pl.discovery ? ` \u00b7 ${esc(pl.discovery)}${pl.year ? ` ${esc(pl.year)}` : ''}` : '') + `</span></div>`).join('')
      || `<span class="hint">${esc(id.planets_error || t('none_known'))}</span>`;
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
    box.innerHTML = '<div class="stats-grid">'
      + tile(t('main'), esc(id.main)) + tile(t('tic'), esc((id.tic || '-').replace('TIC ', '')))
      + tile(t('pos'), pos) + tile(t('gaia'), esc((id.gaia_dr3 || '-').replace('Gaia DR3 ', '')))
      + tile(t('names'), others, true) + tile(t('rotation'), rot.join('<br>') || esc(t('none')), true)
      + tile(t('known'), planets, true)
      + tile(t('toi_title'), toiTile, true)
      + tile(t('carmenes'), carm, true) + '</div>'
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
    dace: asked.dace ? '1' : '', carmenes: asked.carmenes ? '1' : '' });
  try {
    const res = await api(`/api/rv?${q}`);
    if (await drawVelocities(res, auto ? [t('arch_auto')] : [])) startQuick();
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
  const traces = res.instruments.map((inst, i) => ({
    x: inst.time, y: inst.rv, name: `${inst.name} (${inst.source}, ${inst.n})`, type: 'scatter', mode: 'markers',
    error_y: { type: 'data', array: inst.err, visible: true, thickness: 1, width: 0, color: COLOURS[i % 8] },
    marker: { color: COLOURS[i % 8], symbol: SYMBOLS[i % 8], size: 7, line: { color: '#08111f', width: 1 } },
    hovertemplate: `${inst.name}<br>rjd %{x:.4f}<br>%{y:.2f} m/s<extra></extra>`,
  }));
  $('plotcard').classList.add('on');
  div.classList.add('on');
  const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
  lastRV = res.instruments;
  if (window.Plotly) {
    await Plotly.newPlot(div, traces, {
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)',
      font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
      margin: PLOT_MARGIN, showlegend: traces.length > 1,
      legend: { orientation: 'h', x: 0, y: 1.0, yanchor: 'bottom' },
      xaxis: { ...axis, title: 'BJD - 2400000', tickformat: '.0f', exponentformat: 'none' }, yaxis: { ...axis, title: 'RV - median [m/s]' },
    }, { responsive: true, displaylogo: false });
    div.removeAllListeners && div.removeAllListeners('plotly_relayout');
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
      const name = lastRV[ev.curveNumber].name;
      setExcluded(name, !excludedSet().has(name.toUpperCase()));
      return false;
    });
  }
  styleExcluded();
  return true;
}

// -----------------------------------------------------------------------------
// the ranges of the plot: the sliders along its axes (the y one in asinh,
//   fine around the bulk of the points, coarse toward the outliers), twice
//   the 3 to 97 percentile range, and the instruments kept
// -----------------------------------------------------------------------------
const PLOT_MARGIN = { l: 62, r: 12, t: 34, b: 48 };
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

// the domains of the sliders: every time, and the velocities kept
function domains() {
  const all = lastRV || [];
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
  const h = Math.max(60, div.clientHeight - PLOT_MARGIN.t - PLOT_MARGIN.b);
  const box = $('yslider');
  box.style.paddingTop = `${PLOT_MARGIN.t}px`;
  box.style.height = `${h}px`;
  $('ydual').style.width = `${h}px`;
  $('ydual').style.transform = `translateY(${h}px) rotate(-90deg)`;
  $('xdual').style.marginLeft = `${PLOT_MARGIN.l}px`;
  $('xdual').style.marginRight = `${PLOT_MARGIN.r}px`;
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
  if (ev['yaxis.range[0]'] !== undefined || ev['yaxis.autorange']) foldY();
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
  checkStale();
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
    dace: !!asked.dace, carmenes: !!asked.carmenes, exclude: $('exclude').value,
    trend: !!asked.trend, curvature: !!asked.curvature };
}

async function startQuick() {
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
  if (state.result && quick.drawn !== drawn) { drawFip(state.result, each); quick.drawn = drawn; }
  $('fipstatus').innerHTML = (state.result ? fipSummary(state.result, state.elapsed) : '') + quickRunning(state);
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
  return `<p class="hint"><span class="hourglass">\u23f3</span> ${esc(t(step))}${detail} \u00b7 ${clock(state.elapsed)}</p><div class="steps-bar">${bar}</div>`;
}

function fipSummary(r, elapsed) {
  const insts = Object.entries(r.instruments).map(([k, v]) => `${k} ${v}`).join(', ');
  const best = r.peaks[0];
  return `<p class="hint">${r.n} ${esc(t('nights'))} (${esc(insts)}); ${esc(t('no_gp'))}, ${r.settings.kmax} ${esc(t('signals_word'))}, `
    + `${r.settings.nsweep} ${esc(t('sweeps_word'))}, ${r.passes} ${esc(t('passes'))} \u00b7 ${clock(elapsed)}`
    + (best ? ` \u00b7 ${esc(t('strongest'))}: <b>${best.period.toFixed(4)} d</b>, FIP ${best.family.toExponential(1)}` : '') + '</p>'
    + `<p class="hint accel">${accelText(r)}</p>`;
}

// the acceleration of the star the quick look measured (and its change)
function accelValue(v, unit) {
  const [val, lo, hi] = v;
  const err = Math.abs(lo - hi) < 0.05 * Math.max(lo, hi) ? `\u00b1 ${(0.5 * (lo + hi)).toPrecision(2)}`
    : `\u2212${lo.toPrecision(2)} +${hi.toPrecision(2)}`;
  return `${val >= 0 ? '+' : '\u2212'}${Math.abs(val).toPrecision(3)} ${err} ${unit}`;
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
  // the numbered peaks
  for (const pk of r.peak_list.filter((p) => p.named)) {
    notes.push({ x: Math.log10(pk.period), y: -Math.log10(Math.max(pk.family, 1e-15)), xref: 'x', yref: 'y',
      text: `<b>#${pk.id}</b>`, showarrow: true, arrowhead: 0, ax: 0, ay: -22, font: { size: 12, color: '#ffffff' },
      arrowcolor: '#a8b4ca', hovertext: `#${pk.id}: ${pk.period.toFixed(4)} d, FIP ${pk.family.toExponential(1)}` });
  }
  const traces = [
    { x: r.period, y: r.alone, name: t('alone'), type: 'scatter', mode: 'lines', line: { color: '#8a93a3', width: 1 },
      legendrank: 1, hovertemplate: '%{x:.4f} d<br>-log10 FIP %{y:.2f}<extra></extra>' },
  ];
  // each instrument on its own: its FIP of the period or any of its aliases
  for (const one of each.filter((x) => !x.skipped)) {
    traces.push({ x: one.period, y: one.family, name: `${one.name} ${t('inst_alone')}`, type: 'scatter', mode: 'lines',
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
    div.removeAllListeners && div.removeAllListeners('plotly_relayout');
    div.on('plotly_relayout', (ev) => {
      if (!pview) return;
      if (ev['xaxis.range[0]'] !== undefined) pview.p = [10 ** ev['xaxis.range[0]'], 10 ** ev['xaxis.range[1]']];
      if (ev['xaxis.autorange']) pview.p = pview.dom.slice();
      syncPeriods();
    });
  });
  const per = r.period;
  pview = { of: quick.id, dom: [per[0], per[per.length - 1]], p: keepView || [per[0], per[per.length - 1]] };
  syncPeriods();
  drawEach(each);
  if (keepView) return;
  // the folds at the numbered peaks
  $('foldbuttons').innerHTML = (r.folds || []).map((f) => `<button type="button" class="small" data-fold="${f.id}">#${f.id} \u00b7 ${f.period.toFixed(4)} d</button>`).join('');
  if (r.folds && r.folds.length) showFold(r.folds[0].id);
}

// each instrument's nights and best peaks, under the FIP
function drawEach(each) {
  if (!each.length) { $('fipeach').innerHTML = ''; return; }
  $('fipeach').innerHTML = `<table class="mini"><tr><th>${esc(t('each_title'))}</th><th>${esc(t('nights'))}</th><th>${esc(t('best_peaks'))}</th></tr>`
    + each.map((one) => `<tr><td><span class="swatch" style="background:${instColour(one.name)}"></span>${esc(one.name)}</td>`
      + `<td class="num">${one.n}</td><td class="peaks">${one.skipped ? `<span class="hint">${esc(t('few_nights'))}</span>`
        : (one.peaks || []).map((pk) => `${pk.period.toFixed(4)} (${pk.family.toExponential(1)})`).join(' \u00b7 ') || esc(t('none_found'))}</td></tr>`).join('')
    + '</table>';
}

function showFold(id) {
  const f = ((quick && quick.result && quick.result.folds) || []).find((x) => x.id === +id);
  if (!f) return;
  document.querySelectorAll('[data-fold]').forEach((b) => b.classList.toggle('on', +b.dataset.fold === f.id));
  $('foldnote').innerHTML = `<b>#${f.id}</b>: P = ${f.period.toFixed(4)} d, K = ${f.K.toFixed(2)} \u00b1 ${f.K_err.toFixed(2)} m/s, `
    + `rms ${f.rms.toFixed(2)} m/s \u00b7 ${esc(t('fold_note'))}`;
  const order = (lastRV || []).map((inst) => inst.name);
  const traces = f.instruments.map((inst) => {
    const i = Math.max(0, order.indexOf(inst.name));
    const valid = inst.valid ? `<br>${esc(t('p_valid'))} %{customdata:.2f}` : '';
    return { x: inst.phase, y: inst.rv, name: inst.name, type: 'scatter', mode: 'markers', customdata: inst.valid,
      error_y: { type: 'data', array: inst.err, visible: true, thickness: 1, width: 0, color: COLOURS[i % 8] },
      marker: { color: COLOURS[i % 8], symbol: SYMBOLS[i % 8], size: 7, line: { color: '#08111f', width: 1 } },
      hovertemplate: `${inst.name}<br>phase %{x:.3f}<br>%{y:.2f} m/s${valid}<extra></extra>` };
  });
  traces.push({ x: f.curve.phase, y: f.curve.rv, name: `K = ${f.K.toFixed(2)} m/s`, type: 'scatter', mode: 'lines',
    line: { color: '#e8eef8', width: 2 }, hoverinfo: 'skip' });
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
  const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
  $('foldplot').classList.add('on');
  Plotly.newPlot('foldplot', traces, {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)', font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
    margin: { ...PLOT_MARGIN, t: 30 }, legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom' },
    xaxis: { ...axis, title: t('phase_axis'), range: [0, 1] },
    // the velocities on the range of the series above
    yaxis: { ...axis, title: 'RV [m/s]', ...(view && view.y ? { range: view.y.slice() } : {}) },
  }, { responsive: true, displaylogo: false });
}

// the fold follows the y range of the series
function foldY() {
  if (view && view.y && window.Plotly && $('foldplot').data) Plotly.relayout('foldplot', { 'yaxis.range': view.y.slice() });
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
        quick: quick && quick.result ? quick.id : '', command: $('cmd-detailed').textContent }) });
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

// a new target: every field back to the page's own default, no file, no
//   star, no plot, no run that ended, as fresh as a new session (a run that
//   still runs stays: closing a page does not stop it either)
function resetPage() {
  document.querySelectorAll('main input, main select').forEach((el) => {
    if (el.type === 'checkbox') el.checked = el.defaultChecked;
    else if (el.tagName === 'SELECT') {
      const def = [...el.options].find((opt) => opt.defaultSelected) || el.options[0];
      el.value = def ? def.value : '';
    } else el.value = el.defaultValue;
  });
  fileRows = [{ path: '', label: '' }];
  renderFiles();
  $('ident').innerHTML = '';
  $('rvnote').textContent = '';
  $('rvtable').innerHTML = '';
  if (window.Plotly) Plotly.purge($('rvplot'));
  $('rvplot').classList.remove('on');
  $('plotcard').classList.remove('on');
  $('fipcard').classList.remove('on');
  if (window.Plotly) { Plotly.purge($('fipplot')); Plotly.purge($('foldplot')); }
  $('foldbuttons').innerHTML = ''; $('foldnote').textContent = '';
  $('fipeach').innerHTML = ''; $('fipstatus').innerHTML = '';
  $('fipstale').textContent = ''; $('remstate').textContent = ''; $('remember').disabled = true;
  syncMirrors();
  view = null;
  quick = null;
  pview = null;
  lastRV = null;
  onDisk = null;
  showDiskPoints();
}

function newTarget() {
  resetPage();
  showTab('analysis');
  // the runs that ended, the quick looks, the archive in memory: forgotten
  //   here and by the server (a run still running is kept: Stop stops it)
  for (const [id, job] of [...jobs]) if (job.status !== 'running') { jobs.delete(id); openLogs.delete(id); }
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
  if (name === 'remembered') loadRemembered();
  else if (view) setTimeout(() => { syncSliders(); syncPeriods(); }, 50);
}

// the page as it is: what a recall puts back
function pageState() {
  return { target: $('target').value.trim(), files: filesNow(), root: $('root').value.trim(), outdir: $('outdir').value.trim(),
    detailed: readOptions('detailed'), clip: $('clip').checked,
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
  const exp = (v) => (+v).toExponential(1);
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
  showTab('analysis');
  resetPage();
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
  if (pview && page.periods) {
    pview.p = page.periods;
    Plotly.relayout('fipplot', { 'xaxis.range': [Math.log10(pview.p[0]), Math.log10(pview.p[1])] });
    syncPeriods();
  }
  $('remnote').value = res.entry.note || '';
  $('remstate').textContent = `${t('recalled_from')} ${res.entry.created}`;
  $('plotcard').scrollIntoView({ behavior: 'smooth', block: 'start' });
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
  if (e.target.id === 'exclude') styleExcluded();
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
  if (e.target.dataset && e.target.dataset.for === 'detailed' && ['dace', 'carmenes'].includes(e.target.dataset.opt)
      && $('plotcard').classList.contains('on')) plotVelocities();
  // a period of the literature ticked: the SHO at it; none: back to the bands
  if (e.target.name === 'prot') {
    $('rotation').value = e.target.value;
    $('fipgp').value = e.target.value ? 'sho' : 'banded';
    updateCommands();
  }
});
document.addEventListener('change', updateCommands);
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
$('pdf').addEventListener('click', quicklookPdf);
document.addEventListener('click', (e) => { const b = e.target.closest('[data-fold]'); if (b) showFold(b.dataset.fold); });
$('fullrange').addEventListener('click', () => { $('clip').checked = false; view = null; refitY(); });
window.addEventListener('resize', () => { if (view) setTimeout(syncSliders, 100); });
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

(async () => {
  $('detailed-options').innerHTML = renderOptions('detailed');
  langHooks.push(showCwd, renderJobs, renderFiles, checkArchives, renderRemembered,
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

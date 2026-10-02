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
    const where = id.disk ? `${t('from_disk')} (${id.disk}, ${id.disk_date})` : t('asked_now');
    if (rot.length) rot.push(`<label class="prot"><input type="radio" name="prot" value=""${now ? '' : ' checked'}> ${esc(t('no_sho'))}</label>`);
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
  // the archives the report will use, and only those
  const asked = readOptions('detailed');
  const q = new URLSearchParams({ files: JSON.stringify(filesNow()), target: $('target').value.trim(), root: $('root').value.trim(),
    dace: asked.dace ? '1' : '', carmenes: asked.carmenes ? '1' : '' });
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
    startQuick();
  } catch (err) {
    note.innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
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
    dace: !!asked.dace, carmenes: !!asked.carmenes, exclude: $('exclude').value };
}

async function startQuick() {
  $('fipcard').classList.add('on');
  $('fipstale').textContent = '';
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
  if (state.status === 'running') {
    const step = { waiting: 'waiting', noise: 'quick_noise', fip1: 'quick_fip1', planets: 'quick_planets', fip2: 'quick_fip2' }[state.step] || 'quick_noise';
    let bar = '';
    if (state.progress && ['fip1', 'fip2'].includes(state.step)) {
      const p = state.progress;
      const frac = Math.min(1, p.done / Math.max(p.total, 1));
      const left = frac > 0.02 ? ` \u00b7 ~${clock(p.seconds * (1 - frac) / frac)} ${t('left_fip')}` : '';
      bar = `<span class="pbar"><span style="width:${Math.round(100 * frac)}%"></span></span><span class="ptext">${Math.round(100 * frac)} %${left}</span>`;
    }
    $('fipstatus').innerHTML = `<p class="hint"><span class="hourglass">\u23f3</span> ${esc(t(step))} \u00b7 ${clock(state.elapsed)}</p><div class="steps-bar">${bar}</div>`;
    setTimeout(() => pollQuick(id), 1000);
  } else if (state.status === 'done') {
    drawFip(state.result, state.elapsed);
  } else {
    $('fipstatus').innerHTML = `<p class="hint bad">${esc(state.error || 'failed')}</p>`;
  }
}

function drawFip(r, elapsed) {
  const insts = Object.entries(r.instruments).map(([k, v]) => `${k} ${v}`).join(', ');
  const best = r.peaks[0];
  $('fipstatus').innerHTML = `<p class="hint">${r.n} ${esc(t('nights'))} (${esc(insts)}); ${esc(t('no_gp'))}, ${r.settings.kmax} ${esc(t('signals_word'))}, `
    + `${r.settings.nsweep} ${esc(t('sweeps_word'))}, ${r.passes} ${esc(t('passes'))} \u00b7 ${clock(elapsed)}`
    + (best ? ` \u00b7 ${esc(t('strongest'))}: <b>${best.period.toFixed(4)} d</b>, FIP ${best.family.toExponential(1)}` : '') + '</p>';
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
  Plotly.newPlot(div, [
    { x: r.period, y: r.alone, name: t('alone'), type: 'scatter', mode: 'lines', line: { color: '#8a93a3', width: 1 },
      hovertemplate: '%{x:.4f} d<br>-log10 FIP %{y:.2f}<extra></extra>' },
    { x: r.period, y: r.family, name: t('family'), type: 'scatter', mode: 'lines', line: { color: '#3987e5', width: 1.4 },
      hovertemplate: '%{x:.4f} d<br>-log10 FIP %{y:.2f}<extra></extra>' },
  ], {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)', font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
    margin: { ...PLOT_MARGIN, t: 40 }, legend: { orientation: 'h', x: 0, y: 1.08, yanchor: 'bottom' }, shapes: lines, annotations: notes,
    xaxis: { ...axis, type: 'log', title: t('period_axis') }, yaxis: { ...axis, title: t('fip_axis'), rangemode: 'tozero' },
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
  pview = { dom: [per[0], per[per.length - 1]], p: [per[0], per[per.length - 1]] };
  syncPeriods();
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
        quick: quick && quick.status === 'done' ? quick.id : '' }) });
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
//   star, no plot (the runs stay: closing a page does not stop them)
function newTarget() {
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
  if (window.Plotly) Plotly.purge($('fipplot'));
  view = null;
  quick = null;
  pview = null;
  lastRV = null;
  onDisk = null;
  history.replaceState(null, '', location.pathname);
  updateCommands();
  checkArchives();
  $('target').focus();
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
});

(async () => {
  $('detailed-options').innerHTML = renderOptions('detailed');
  langHooks.push(showCwd, renderJobs, renderFiles, checkArchives);
  jobDoneHooks.push((job) => { if (job.action === 'gather') checkArchives(); if (job.action === 'archive') refreshInfo(); });
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

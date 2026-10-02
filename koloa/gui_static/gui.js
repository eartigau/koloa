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

async function resolveStar(refresh) {
  const name = $('target').value.trim();
  if (!name) return;
  const box = $('ident');
  box.innerHTML = `<p class="hint"><span class="spin"></span> ${esc(t('resolving'))}</p>`;
  try {
    const id = await api(`/api/resolve?${new URLSearchParams({ name, root: $('root').value.trim(), refresh: refresh === true ? '1' : '' })}`);
    const pos = id.ra != null ? `${id.ra.toFixed(5)}, ${id.dec.toFixed(5)}` : '-';
    const others = [id.gj, id.hd, id.hip].filter(Boolean).map(esc).join(', ') || '-';
    const chip = (per, text) => `<button type="button" class="chip" data-prot="${esc(per)}">${text}</button>`;
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
    if (rot.length) rot.push(`<span class="hint">${esc(t('pick'))}</span>`);
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
  const q = new URLSearchParams({ files: JSON.stringify(filesNow()), target: $('target').value.trim(), root: $('root').value.trim() });
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
    div.classList.add('on');
    const axis = { gridcolor: 'rgba(200,220,255,0.10)', zerolinecolor: 'rgba(200,220,255,0.25)', color: '#7a8597' };
    if (window.Plotly) {
      Plotly.newPlot(div, traces, {
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(4,8,16,0.35)',
        font: { family: 'Space Grotesk, sans-serif', color: '#e8eef8' },
        margin: { l: 60, r: 10, t: 10, b: 50 }, legend: { orientation: 'h', y: -0.2 }, showlegend: traces.length > 1,
        xaxis: { ...axis, title: 'BJD - 2400000', tickformat: '.0f', exponentformat: 'none' }, yaxis: { ...axis, title: 'RV - median [m/s]' },
      }, { responsive: true, displaylogo: false });
    }
    lastRV = res.instruments;
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
});
document.addEventListener('change', updateCommands);
$('resolve').addEventListener('click', resolveStar);
$('target').addEventListener('keydown', (e) => { if (e.key === 'Enter') resolveStar(); });
$('plot').addEventListener('click', plotVelocities);
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
  const prot = e.target.closest('[data-prot]');
  if (prot) {
    $('rotation').value = prot.dataset.prot;
    $('fipgp').value = 'sho';
    updateCommands();
  }
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

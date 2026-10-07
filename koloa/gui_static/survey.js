// koloa's GUI: the survey tab. A sample asked of SIMBAD, what the archives
//   have of each star, the files put with them, the batch of those ticked
//   (here, or packed for another machine, which its README tells how to run).

// -----------------------------------------------------------------------------
// the sample
// -----------------------------------------------------------------------------
let survey = null;            // the sample: its id, its stars, its check, its packing
const svTicked = new Set();   // the stars ticked, by name
let svSort = { key: 'distance', dir: 1 };
let svPacked = null;          // the batch packed last: its folder, its tar

const svHas = (star) => {
  const a = star.archives || {};
  return !!(a.dace || a.carmenes || (a.surveys || []).length || (star.files || []).length);
};
// what the archives have of a star: its velocities counted (DACE, CARMENES
//   DR1), its spectrographs
const svKnown = (star) => (star.summary || []).reduce((sum, one) => sum + (one.n || 0), 0);
const svShown = () => {
  const words = $('sv-filter').value.trim().toUpperCase();
  // some data: not the stars cross-matched that have nothing (no velocity
  //   in an archive, no file); a star not cross-matched yet stays
  const some = $('sv-somedata').checked;
  const rows = (survey ? survey.stars : []).filter((star) => (!words
    || `${star.name} ${star.main} ${star.sptype} ${(star.summary || []).map((one) => one.name).join(' ')}`.toUpperCase().includes(words))
    && (!some || svHas(star) || !star.archives));
  const { key, dir } = svSort;
  const val = (star) => {
    if (key === 'known') return star.archives ? svKnown(star) : -1;
    if (key === 'summary') return star.archives ? (star.summary || []).length : -1;
    if (key === 'files') return (star.files || []).length;
    if (key === 'tick') return svTicked.has(star.name) ? 1 : 0;
    if (key === 'rotation') return (star.rotation || []).length ? star.rotation[0].period : null;
    return star[key];
  };
  return rows.sort((a, b) => {
    const [x, y] = [val(a), val(b)];
    if (typeof x === 'string' || typeof y === 'string') return dir * String(x ?? '').localeCompare(String(y ?? ''));
    return dir * ((x ?? Infinity) - (y ?? Infinity));
  });
};

// the rotation period of a star: the latest SIMBAD lists (its table of
//   variability), else that of CARMENES DR1; the others under the cursor
function svRotation(star) {
  const all = star.rotation || [];
  if (!all.length) return '';
  const fmt = (p) => (p < 10 ? p.toFixed(2) : p < 100 ? p.toFixed(1) : p.toFixed(0));
  const tip = all.map((one) => `${+one.period.toPrecision(5)} d \u00b7 ${one.source}`).join('\n');
  return `<span title="${esc(tip)}">${fmt(all[0].period)}${all.length > 1 ? ` <span class="hint">+${all.length - 1}</span>` : ''}</span>`;
}

// the summary of the cross-match: how many stars have velocities, and
//   spectrograph by spectrograph (the stars, the velocities counted, the
//   archives that have them)
function svOverview() {
  const o = survey && survey.overview;
  if (!o || !o.checked) { $('sv-overview').innerHTML = ''; return; }
  const rows = o.spectrographs.map((row) => `<tr><td><b>${esc(row.name)}</b></td><td class="num">${row.stars}</td>`
    + `<td class="num">${row.velocities || '<span class="hint">?</span>'}</td>`
    + `<td class="src">${esc(Object.entries(row.where).map(([where, num]) => `${where} ${num}`).join(' · '))}</td></tr>`).join('');
  $('sv-overview').innerHTML = `<p><b>${o.data}</b> ${esc(t('sv_ov_data'))} ${o.checked} ${esc(t('sv_ov_checked'))}`
    + ` (${o.several} ${esc(t('sv_ov_several'))}, ${o.single} ${esc(t('sv_ov_single'))}); <b>${o.none}</b> ${esc(t('sv_ov_none'))}.`
    + ` ${o.velocities} ${esc(t('sv_ov_velocities'))}`
    + (o.files ? ` · ${esc(t('sv_ov_files'))} ${o.nfiles}, ${esc(t('sv_ov_filestars'))} ${o.files}` : '')
    + (o.checked < o.n ? ` · ${o.n - o.checked} ${esc(t('sv_ov_todo'))}` : '') + '</p>'
    + (rows ? `<table class="mini"><tr><th>${esc(t('sv_ov_spec'))}</th><th>${esc(t('sv_stars'))}</th><th>${esc(t('sv_ov_counted'))}</th>`
      + `<th>${esc(t('sv_ov_where'))}</th></tr>${rows}</table>` : '');
}

function svRender() {
  $('sv-card').hidden = !survey;
  $('sv-batchcard').hidden = !survey;
  if (!survey) return;
  svOverview();
  const cols = [['tick', ''], ['name', t('sv_col_name')], ['main', 'SIMBAD'], ['spnum', t('sv_col_type')], ['distance', 'd [pc]'],
    ['V', 'V'], ['rotation', t('sv_col_rot')], ['summary', t('sv_col_arch')], ['known', t('sv_col_known')], ['files', t('sv_col_files')]];
  const head = cols.map(([key, label]) => `<th data-svsort="${key}" class="sortable svcol-${key}${svSort.key === key ? ' sorted' : ''}">${esc(label)}`
    + `${svSort.key === key ? (svSort.dir > 0 ? ' ▲' : ' ▼') : ''}</th>`).join('');
  const rows = svShown();
  const num = (v, d) => (v === null || v === undefined ? '' : (+v).toFixed(d));
  $('sv-table').innerHTML = `<div class="batchwrap svwrap"><table class="mini batch"><tr>${head}</tr>`
    + rows.map((star) => {
      const a = star.archives;
      // spectrograph by spectrograph: its velocities when an archive
      //   counted them, the archives that have it
      const arch = !a ? `<span class="hint">${esc(t('sv_not_checked'))}</span>`
        : (star.summary || []).length ? star.summary.map((one) => `<span class="svspec"><b>${esc(one.name)}</b>${one.n ? ` ${one.n}` : ''}`
          + ` <span class="hint">${esc(one.where.join(', '))}</span></span>`).join('')
          : `<span class="hint" title="${esc(a.dace_error || '')}">${esc(t(a.dace_error ? 'sv_dace_error' : 'sv_nothing'))}</span>`;
      // its files: how many, and their names (the first three; all of them under the cursor)
      const nf = (star.files || []).length;
      const files = nf ? `<span title="${esc(star.files.join('\n'))}"><b>${nf}</b> ${esc(star.files.slice(0, 3).join(', '))}`
        + `${nf > 3 ? ` <span class="hint">+${nf - 3}</span>` : ''}</span>` : '';
      return `<tr class="${a && !svHas(star) ? 'nodata' : ''}"><td><input type="checkbox" data-svtick="${esc(star.name)}"${svTicked.has(star.name) ? ' checked' : ''}></td>`
        + `<td><b>${esc(star.name)}</b>${star.near ? ` <span class="hint" title="${esc(t('sv_near'))} ${esc(star.near)}">⧉</span>` : ''}</td>`
        + `<td class="src">${esc(star.main)}</td><td>${esc(star.sptype)}</td><td class="num">${num(star.distance, 2)}</td>`
        + `<td class="num">${num(star.V ?? star.G, 1)}${star.V === null && star.G !== null ? '<span class="hint">G</span>' : ''}</td>`
        + `<td class="num">${svRotation(star)}</td>`
        + `<td class="svarch">${arch}</td><td class="num">${a ? (svKnown(star) || '') : ''}</td><td class="svfilecell">${files}</td></tr>`;
    }).join('') + '</table></div>';
  const withData = survey.stars.filter(svHas).length;
  $('sv-count').textContent = `${survey.stars.length} ${t('sv_stars')}${survey.check ? ` · ${withData} ${t('sv_with_data')}` : ''} · ${svTicked.size} ${t('sv_ticked')}`
    + (rows.length !== survey.stars.length ? ` · ${rows.length} ${t('sv_shown')}` : '');
  svCheckStatus();
  svPackStatus();
}

function svCheckStatus() {
  const c = survey && survey.check;
  if (!c) { $('sv-checkstatus').textContent = ''; return; }
  $('sv-checkstatus').innerHTML = c.status === 'running' ? `<span class="hourglass">⏳</span> ${c.done} / ${c.total}`
    : c.status === 'failed' ? `<span class="bad">${esc(c.error)}</span>` : `${c.total} ${esc(t('sv_checked'))}`;
  $('sv-check').disabled = c.status === 'running';
}

function svPackStatus() {
  const p = survey && survey.pack;
  $('sv-pack').disabled = !!p && p.status === 'running';
  if (!p) { $('sv-packstatus').textContent = ''; return; }
  if (p.status === 'running') {
    $('sv-packstatus').innerHTML = `<p class="hint"><span class="hourglass">⏳</span> ${p.done} / ${p.total} ${esc(p.star || '')}</p>`;
  } else if (p.status === 'failed') {
    $('sv-packstatus').innerHTML = `<p class="hint"><span class="bad">${esc(p.error)}</span></p>`;
  } else {
    const r = p.result;
    svPacked = r;
    // where it is on this machine, in full: the tar is what is copied to the server
    $('sv-packstatus').innerHTML = `<p class="svpacked"><b>${esc(t('sv_packed'))}</b> ${r.n} ${esc(t('sv_stars'))}, ${r.files} ${esc(t('sv_col_files'))}`
      + (r.missing.length ? ` · <span class="bad">${esc(t('sv_missing'))} ${esc(r.missing.join(', '))}</span>` : '') + '</p>'
      + (r.tar ? `<p class="svpacked">${esc(t('sv_packed_tar'))} <span class="mono sel">${esc(r.tar)}</span> (${(r.size / 1e6).toFixed(1)} MB)</p>` : '')
      + `<p class="svpacked">${esc(t('sv_packed_folder'))} <span class="mono sel">${esc(r.folder)}</span></p>`
      + `<p class="svpacked">${esc(t('sv_packed_root'))} <span class="mono sel">${esc(r.root)}</span></p>`
      + `<p class="hint">${esc(t('sv_packed_readme'))}</p>`;
    if (!$('sv-results').value) $('sv-results').value = r.folder_shown || r.folder;
  }
}

async function svAsk() {
  $('sv-status').innerHTML = `<span class="spin"></span> SIMBAD...`;
  $('sv-ask').disabled = true;
  try {
    survey = await api('/api/survey/sample', { sp_from: $('sv-from').value, sp_to: $('sv-to').value, dmax: $('sv-dmax').value,
      dec_min: $('sv-decmin').value, dec_max: $('sv-decmax').value, vmax: $('sv-vmax').value, dwarfs: $('sv-dwarfs').checked });
    svTicked.clear();
    $('sv-status').textContent = `${survey.stars.length} ${t('sv_stars')}`;
    svRender();
  } catch (err) {
    $('sv-status').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
  }
  $('sv-ask').disabled = false;
}

// the check of the archives, and the packing: followed until they end
let svTimer = null;
async function svFollow() {
  clearTimeout(svTimer);
  if (!survey) return;
  try {
    const now = await api(`/api/survey?id=${survey.id}&stars=0`);
    const wasChecking = survey.check && survey.check.status === 'running';
    survey.check = now.check; survey.pack = now.pack;
    const busy = (now.check && now.check.status === 'running') || (now.pack && now.pack.status === 'running');
    if (wasChecking || !busy) {
      // the archives of the stars checked so far
      const full = await api(`/api/survey?id=${survey.id}`);
      survey.stars = full.stars; survey.overview = full.overview;
      svRender();
    } else { svCheckStatus(); svPackStatus(); }
    if (busy) svTimer = setTimeout(svFollow, 1500);
  } catch (err) { svTimer = setTimeout(svFollow, 4000); }
}

async function svCheck() {
  if (!survey) return;
  try {
    const now = await api('/api/survey/check', { id: survey.id, root: $('root').value.trim(), dace: $('sv-dace').checked,
      carmenes: $('sv-carmenes').checked, vizier: $('sv-vizier').checked, names: svTicked.size ? [...svTicked] : null });
    survey.check = now.check;
    svCheckStatus();
    svFollow();
  } catch (err) { $('sv-checkstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

async function svMatch() {
  if (!survey) return;
  $('sv-matchstatus').innerHTML = '<span class="spin"></span>';
  try {
    const now = await api('/api/survey/files', { id: survey.id, folders: $('sv-folders').value, pattern: $('sv-pattern').value });
    survey.stars = now.stars; survey.unmatched = now.unmatched; survey.overview = now.overview;
    const other = now.unmatched.length ? ` · <span title="${esc(now.unmatched.map((u) => `${u.name}: ${u.target || u.raw || '?'}`).join('\n'))}">${now.unmatched.length} ${esc(t('sv_unmatched'))}</span>` : '';
    $('sv-matchstatus').innerHTML = `${now.matched} ${esc(t('sv_matched'))}${other}`;
    svRender();
  } catch (err) { $('sv-matchstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

function svTick(which) {
  if (!survey) return;
  if (which === 'none') svTicked.clear();
  else {
    for (const star of svShown()) {
      if (which === 'all' || (which === 'data' && svHas(star)) || (which === 'files' && (star.files || []).length)) svTicked.add(star.name);
    }
  }
  svRender();
}

async function svRun() {
  if (!survey || !svTicked.size) { $('sv-runstatus').textContent = t('sv_none_ticked'); return; }
  try {
    batch = await api('/api/survey/run', { id: survey.id, names: [...svTicked], root: $('root').value.trim(), rules: $('sv-rules').checked,
      report: $('sv-report').checked, report_fip: $('sv-reportfip').value, options: readOptions('detailed') });
    $('batcharchives').checked = true;
    showTab('batch');
    pollBatch(batch.id);
  } catch (err) { $('sv-runstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

async function svPack() {
  if (!survey || !svTicked.size) { $('sv-packstatus').textContent = t('sv_none_ticked'); return; }
  try {
    const now = await api('/api/survey/pack', { id: survey.id, names: [...svTicked], name: $('sv-name').value, out: $('sv-out').value,
      root: $('root').value.trim(), server_root: $('sv-root').value, jobs: $('sv-jobs').value, rules: $('sv-rules').checked,
      gather: $('sv-gather').checked, tess: $('sv-tess').checked, report: $('sv-report').checked, report_fip: $('sv-reportfip').value,
      options: readOptions('detailed') });
    survey.pack = now.pack;
    svPackStatus();
    svFollow();
  } catch (err) { $('sv-packstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

// the summary of a batch folder as one PDF: made (from what the batch kept) when it is not there, then opened
async function svSummary() {
  const root = $('sv-results').value.trim();
  $('sv-openstatus').innerHTML = `<span class="hourglass">⏳</span> ${esc(t('sv_summary_wait'))}`;
  try {
    const made = await api('/api/survey/summary', { root });
    $('sv-openstatus').innerHTML = `<a target="_blank" href="/api/batchsummary?root=${encodeURIComponent(made.root)}">${esc(made.path)}</a>`;
    window.open(`/api/batchsummary?root=${encodeURIComponent(made.root)}`, '_blank');
  } catch (err) { $('sv-openstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

async function svOpenResults() {
  try {
    batch = await api('/api/survey/results', { root: $('sv-results').value.trim() });
    showTab('batch');
    renderBatch();
    if (batch.todo) $('batchstatus').innerHTML += ` · ${batch.todo} ${esc(t('sv_todo'))}`;
  } catch (err) { $('sv-openstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

document.addEventListener('click', (e) => {
  const id = e.target.id;
  if (id === 'sv-ask') svAsk();
  else if (id === 'sv-check') svCheck();
  else if (id === 'sv-match') svMatch();
  else if (id === 'sv-tickdata') svTick('data');
  else if (id === 'sv-tickfiles') svTick('files');
  else if (id === 'sv-tickall') svTick('all');
  else if (id === 'sv-ticknone') svTick('none');
  else if (id === 'sv-run') svRun();
  else if (id === 'sv-pack') svPack();
  else if (id === 'sv-open') svOpenResults();
  else if (id === 'sv-summary') svSummary();
  const sort = e.target.closest('[data-svsort]');
  if (sort) {
    const key = sort.dataset.svsort;
    svSort = { key, dir: svSort.key === key ? -svSort.dir : (['known', 'summary', 'files', 'tick'].includes(key) ? -1 : 1) };
    svRender();
  }
});
document.addEventListener('change', (e) => {
  if (e.target.dataset && e.target.dataset.svtick !== undefined) {
    if (e.target.checked) svTicked.add(e.target.dataset.svtick); else svTicked.delete(e.target.dataset.svtick);
    svRender();
  }
  if (e.target.id === 'sv-somedata') svRender();
});
document.addEventListener('input', (e) => {
  if (e.target.id === 'sv-filter') svRender();
});
langHooks.push(() => { if (survey) svRender(); });

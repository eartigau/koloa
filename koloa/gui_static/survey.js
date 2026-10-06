// koloa's GUI: the survey tab. A sample asked of SIMBAD, what the archives
//   have of each star, the files put with them, the batch of those ticked
//   (here, or packed for another machine), and a terminal to that machine.

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
const svShown = () => {
  const words = $('sv-filter').value.trim().toUpperCase();
  const rows = (survey ? survey.stars : []).filter((star) => !words
    || `${star.name} ${star.main} ${star.sptype}`.toUpperCase().includes(words));
  const { key, dir } = svSort;
  const val = (star) => {
    if (key === 'dace') return ((star.archives || {}).dace || {}).n ?? -1;
    if (key === 'carmenes') return (star.archives || {}).carmenes ?? -1;
    if (key === 'surveys') return ((star.archives || {}).surveys || []).length;
    if (key === 'files') return (star.files || []).length;
    if (key === 'tick') return svTicked.has(star.name) ? 1 : 0;
    return star[key];
  };
  return rows.sort((a, b) => {
    const [x, y] = [val(a), val(b)];
    if (typeof x === 'string' || typeof y === 'string') return dir * String(x ?? '').localeCompare(String(y ?? ''));
    return dir * ((x ?? Infinity) - (y ?? Infinity));
  });
};

function svRender() {
  $('sv-card').hidden = !survey;
  $('sv-batchcard').hidden = !survey;
  if (!survey) return;
  const cols = [['tick', ''], ['name', t('sv_col_name')], ['main', 'SIMBAD'], ['spnum', t('sv_col_type')], ['distance', 'd [pc]'],
    ['V', 'V'], ['dace', 'DACE'], ['carmenes', 'CARMENES DR1'], ['surveys', t('sv_col_surveys')], ['files', t('sv_col_files')]];
  const head = cols.map(([key, label]) => `<th data-svsort="${key}" class="sortable${svSort.key === key ? ' sorted' : ''}">${esc(label)}`
    + `${svSort.key === key ? (svSort.dir > 0 ? ' ▲' : ' ▼') : ''}</th>`).join('');
  const rows = svShown();
  const num = (v, d) => (v === null || v === undefined ? '' : (+v).toFixed(d));
  $('sv-table').innerHTML = `<div class="batchwrap svwrap"><table class="mini batch"><tr>${head}</tr>`
    + rows.map((star) => {
      const a = star.archives;
      const dace = !a ? '' : a.dace ? `<b>${a.dace.n}</b> <span class="hint">${esc(Object.entries(a.dace.instruments).map(([k, v]) => `${k} ${v}`).join(', '))}</span>`
        : `<span class="hint" title="${esc(a.dace_error || '')}">${a.dace_error ? '?' : '-'}</span>`;
      const files = (star.files || []).length ? `<span title="${esc(star.files.join('\n'))}"><b>${star.files.length}</b> <span class="hint">${esc(star.files[0])}${star.files.length > 1 ? ', ...' : ''}</span></span>` : '';
      return `<tr class="${a && !svHas(star) ? 'nodata' : ''}"><td><input type="checkbox" data-svtick="${esc(star.name)}"${svTicked.has(star.name) ? ' checked' : ''}></td>`
        + `<td><b>${esc(star.name)}</b>${star.near ? ` <span class="hint" title="${esc(t('sv_near'))} ${esc(star.near)}">⧉</span>` : ''}</td>`
        + `<td class="src">${esc(star.main)}</td><td>${esc(star.sptype)}</td><td class="num">${num(star.distance, 2)}</td>`
        + `<td class="num">${num(star.V ?? star.G, 1)}${star.V === null && star.G !== null ? '<span class="hint">G</span>' : ''}</td>`
        + `<td>${dace}</td><td class="num">${a ? (a.carmenes || '-') : ''}</td>`
        + `<td class="src">${a ? esc((a.surveys || []).join(', ') || '-') : ''}</td><td>${files}</td></tr>`;
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
    $('sv-packstatus').innerHTML = `<span class="hourglass">⏳</span> ${p.done} / ${p.total} ${esc(p.star || '')}`;
  } else if (p.status === 'failed') {
    $('sv-packstatus').innerHTML = `<span class="bad">${esc(p.error)}</span>`;
  } else {
    const r = p.result;
    svPacked = r;
    $('sv-packstatus').innerHTML = `${esc(t('sv_stars'))}: ${r.n}, ${esc(t('sv_col_files'))}: ${r.files}· <span class="mono">${esc(r.tar_shown || r.tar || r.folder_shown || r.folder)}</span>`
      + (r.size ? ` (${(r.size / 1e6).toFixed(1)} MB)` : '') + ` · ROOT <span class="mono">${esc(r.root)}</span>`
      + (r.missing.length ? ` · <span class="bad">${esc(t('sv_missing'))} ${esc(r.missing.join(', '))}</span>` : '');
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
      survey.stars = full.stars;
      svRender();
    } else { svCheckStatus(); svPackStatus(); }
    if (busy) svTimer = setTimeout(svFollow, 1500);
  } catch (err) { svTimer = setTimeout(svFollow, 4000); }
}

async function svCheck() {
  if (!survey) return;
  try {
    const now = await api('/api/survey/check', { id: survey.id, root: $('root').value.trim(), dace: $('sv-dace').checked,
      names: svTicked.size ? [...svTicked] : null });
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
    survey.stars = now.stars; survey.unmatched = now.unmatched;
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
      options: readOptions('detailed') });
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
      gather: $('sv-gather').checked, tess: $('sv-tess').checked, options: readOptions('detailed') });
    survey.pack = now.pack;
    svPackStatus();
    svFollow();
  } catch (err) { $('sv-packstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

async function svOpenResults() {
  try {
    batch = await api('/api/survey/results', { root: $('sv-results').value.trim() });
    showTab('batch');
    renderBatch();
    if (batch.todo) $('batchstatus').innerHTML += ` · ${batch.todo} ${esc(t('sv_todo'))}`;
  } catch (err) { $('sv-openstatus').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

// where the batch will be on the server: the folder of the route, and its name
function svServerRoot() {
  const route = termRoutes[$('term-route').value];
  if (route && route.folder) $('sv-root').value = `${route.folder.replace(/\/+$/, '')}/${$('sv-name').value.trim() || 'koloa_batch'}`;
}

// -----------------------------------------------------------------------------
// the terminal: shells of this machine in the page, to go to the server
// -----------------------------------------------------------------------------
let termKey = '';
try {
  termKey = new URLSearchParams(location.search).get('key') || sessionStorage.getItem('koloa-key') || '';
  if (termKey) sessionStorage.setItem('koloa-key', termKey);
} catch (err) { /* no storage */ }
const terms = [];          // the terminals: { id, term, fit, since, alive, box, label }
let termNow = null;        // the one shown
let termRoutes = {};       // the routes kept, by name

async function termApi(what, body) {
  const resp = await fetch(`/api/term/${what}`, { method: 'POST',
    headers: { 'Content-Type': 'application/json', 'X-Koloa-Key': termKey }, body: JSON.stringify(body || {}) });
  const data = await resp.json();
  if (!resp.ok || data.error) throw new Error(data.error || resp.statusText);
  return data;
}

function termRenderTabs() {
  $('term-tabs').innerHTML = terms.map((one, i) => `<button type="button" class="small${one === termNow ? ' on' : ''}" data-term="${i}">`
    + `${esc(one.label)}${one.alive ? '' : ' †'}</button>`).join('')
    + (terms.length ? ` <button type="button" class="small" id="term-close">×</button>` : '');
  terms.forEach((one) => { one.box.hidden = one !== termNow; });
  if (termNow) { try { termNow.fit.fit(); termNow.term.focus(); } catch (err) { /* not shown yet */ } }
}

async function termRead(one) {
  while (one.alive) {
    try {
      const res = await termApi('read', { id: one.id, since: one.since, wait: 20 });
      if (res.data) one.term.write(Uint8Array.from(atob(res.data), (c) => c.charCodeAt(0)));
      one.since = res.next;
      one.alive = res.alive;
    } catch (err) {
      // the server asleep, or started again: once more, then given up
      one.misses = (one.misses || 0) + 1;
      if (one.misses > 5) one.alive = false;
      await new Promise((done) => setTimeout(done, 2000));
    }
  }
  termRenderTabs();
}

async function termOpen(label) {
  if (!window.Terminal) { $('term-status').textContent = t('term_no_lib'); return null; }
  const box = document.createElement('div');
  box.className = 'termbox';
  $('term-box').appendChild(box);
  const term = new window.Terminal({ fontFamily: '"IBM Plex Mono", monospace', fontSize: 13, cursorBlink: true, scrollback: 5000,
    theme: { background: '#05080f', foreground: '#e8eef8', cursor: '#7fd1ff' } });
  const fit = new window.FitAddon.FitAddon();
  term.loadAddon(fit);
  term.open(box);
  fit.fit();
  let made;
  try {
    made = await termApi('open', { cols: term.cols, rows: term.rows });
  } catch (err) {
    box.remove();
    $('term-status').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
    return null;
  }
  const one = { id: made.id, term, fit, since: 0, alive: true, box, label: label || `${t('term_tab')} ${terms.length + 1}` };
  term.onData((data) => { if (one.alive) termApi('write', { id: one.id, data }).catch(() => {}); });
  term.onResize(({ cols, rows }) => { if (one.alive) termApi('resize', { id: one.id, cols, rows }).catch(() => {}); });
  terms.push(one);
  termNow = one;
  termRenderTabs();
  termRead(one);
  $('term-status').textContent = '';
  return one;
}

// a command put in a terminal, not run: it is read, changed if need be, and Enter runs it
async function termType(one, text) {
  if (!one) return;
  termNow = one;
  termRenderTabs();
  await termApi('write', { id: one.id, data: text });
}
// the terminal of this machine (for rsync, which starts here), made when there is none
async function termLocal() {
  return terms.find((one) => one.local && one.alive) || (async () => {
    const one = await termOpen(t('term_here'));
    if (one) one.local = true;
    // its prompt first: what is typed before it is shown twice
    await new Promise((done) => setTimeout(done, 1200));
    return one;
  })();
}
const quoted = (text) => `'${String(text).replace(/'/g, "'\\''")}'`;

async function termLoadRoutes() {
  try {
    const state = await termApi('state');
    if (!state.available) { $('term-status').textContent = t('term_none'); $('term-card').classList.add('off'); return; }
    if (!state.allowed) { $('term-status').textContent = t('term_key'); $('term-card').classList.add('off'); return; }
    termRoutes = (await termApi('routes')).routes;
  } catch (err) { $('term-status').innerHTML = `<span class="bad">${esc(err.message)}</span>`; return; }
  const kept = $('term-route').value;
  $('term-route').innerHTML = Object.keys(termRoutes).map((name) => `<option value="${esc(name)}">${esc(name)}`
    + `${termRoutes[name].host ? ` (${esc(termRoutes[name].host)})` : ''}</option>`).join('') || `<option value="">${esc(t('term_no_route'))}</option>`;
  if (kept && termRoutes[kept]) $('term-route').value = kept;
  svServerRoot();
}

async function termGo() {
  const name = $('term-route').value;
  if (!termRoutes[name]) return;
  const one = (termNow && termNow.alive && !termNow.local) ? termNow : await termOpen(name);
  if (!one) return;
  one.label = name; one.route = name;
  termRenderTabs();
  await termApi('go', { id: one.id, name });
}

async function termRemember() {
  // what was typed in this terminal (what it showed: never a password), to be corrected
  let lines = [];
  if (termNow) { try { lines = (await termApi('typed', { id: termNow.id })).lines; } catch (err) { /* none */ } }
  const known = termRoutes[$('term-route').value] || {};
  $('route-name').value = (termNow && termNow.route) || $('term-route').value || '';
  $('route-host').value = known.host || (lines.map((l) => (l.match(/^ssh\s+(?:-\S+\s+)*([^\s-]\S*)/) || [])[1]).filter(Boolean)[0] || '');
  $('route-folder').value = known.folder || '';
  $('route-lines').value = lines.join('\n');
  $('route-dialog').showModal();
}

async function termSaveRoute() {
  try {
    termRoutes = (await termApi('route_save', { name: $('route-name').value, lines: $('route-lines').value,
      host: $('route-host').value, folder: $('route-folder').value })).routes;
    $('route-dialog').close();
    const name = $('route-name').value.trim();
    await termLoadRoutes();
    $('term-route').value = name;
    svServerRoot();
  } catch (err) { $('route-error').textContent = err.message; }
}

// the batch packed last, and the route chosen: what the buttons type
function termBatch() {
  const route = termRoutes[$('term-route').value];
  if (!route || !route.host || !route.folder) { $('term-status').textContent = t('term_need_route'); return null; }
  const name = (svPacked && svPacked.folder.split('/').pop()) || $('sv-name').value.trim();
  if (!name) { $('term-status').textContent = t('term_need_batch'); return null; }
  $('term-status').textContent = '';
  return { route, name, folder: route.folder.replace(/\/+$/, ''), tar: svPacked && (svPacked.tar_shown || svPacked.tar),
    local: svPacked && (svPacked.folder_shown || svPacked.folder) };
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
  else if (id === 'term-new') termOpen();
  else if (id === 'term-go') termGo();
  else if (id === 'term-remember') termRemember();
  else if (id === 'route-save') termSaveRoute();
  else if (id === 'route-cancel') $('route-dialog').close();
  else if (id === 'term-forget') {
    const name = $('term-route').value;
    if (name && confirm(`${t('term_forget')} ${name}?`)) termApi('route_delete', { name }).then(termLoadRoutes);
  } else if (id === 'term-close' && termNow) {
    const one = termNow;
    termApi('close', { id: one.id }).catch(() => {});
    one.alive = false; one.term.dispose(); one.box.remove();
    terms.splice(terms.indexOf(one), 1);
    termNow = terms[terms.length - 1] || null;
    termRenderTabs();
  } else if (id === 'term-copy') {
    const b = termBatch();
    if (b && !b.tar) $('term-status').textContent = t('term_need_batch');
    else if (b) termLocal().then((one) => termType(one, `rsync -av --progress ${quoted(b.tar)} ${b.route.host}:${quoted(`${b.folder}/`)}`));
  } else if (id === 'term-launch') {
    const b = termBatch();
    if (b && termNow && !termNow.local) {
      termType(termNow, `cd ${quoted(b.folder)} && tar xzf ${quoted(`${b.name}.tar.gz`)} && cd ${quoted(b.name)} && nohup python -u run_batch.py > run.log 2>&1 &`);
    } else if (b) $('term-status').textContent = t('term_need_server');
  } else if (id === 'term-follow') {
    const b = termBatch();
    if (b && termNow && !termNow.local) termType(termNow, `tail -n 30 -f ${quoted(`${b.folder}/${b.name}/run.log`)}`);
    else if (b) $('term-status').textContent = t('term_need_server');
  } else if (id === 'term-fetch') {
    const b = termBatch();
    if (b) {
      const local = b.local || `${$('sv-out').value.trim() || '.'}/${b.name}`;
      $('sv-results').value = local;
      termLocal().then((one) => termType(one, `rsync -av ${b.route.host}:${quoted(`${b.folder}/${b.name}/results/`)} ${quoted(`${local}/results/`)}`));
    }
  }
  const tab = e.target.closest('[data-term]');
  if (tab) { termNow = terms[+tab.dataset.term]; termRenderTabs(); }
  const sort = e.target.closest('[data-svsort]');
  if (sort) {
    const key = sort.dataset.svsort;
    svSort = { key, dir: svSort.key === key ? -svSort.dir : (['dace', 'carmenes', 'surveys', 'files', 'tick'].includes(key) ? -1 : 1) };
    svRender();
  }
});
document.addEventListener('change', (e) => {
  if (e.target.dataset && e.target.dataset.svtick !== undefined) {
    if (e.target.checked) svTicked.add(e.target.dataset.svtick); else svTicked.delete(e.target.dataset.svtick);
    svRender();
  }
  if (e.target.id === 'term-route') svServerRoot();
});
document.addEventListener('input', (e) => {
  if (e.target.id === 'sv-filter') svRender();
  if (e.target.id === 'sv-name') svServerRoot();
});
window.addEventListener('resize', () => { if (termNow) { try { termNow.fit.fit(); } catch (err) { /* hidden */ } } });
langHooks.push(() => { if (survey) svRender(); });
// the terminal's routes, once the survey tab is shown
let termLoaded = false;
function surveyShown() {
  if (!termLoaded) { termLoaded = true; termLoadRoutes(); }
  if (termNow) setTimeout(() => { try { termNow.fit.fit(); } catch (err) { /* hidden */ } }, 50);
}

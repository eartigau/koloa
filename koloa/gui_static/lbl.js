// The LBL tab: the public spectra of a star fetched from the archives, and
// LBL (lbl-rapide) run on them, a run of its own (koloa --lbl); its rdb
// kept in the folder of the star, to be given to the analysis as a file
// (koloa.lblrun).

let lblNow = null;        // what the server says: the choices, the instruments
let lblHave = {};         // what the archives have, per instrument (once asked)
const lblTicked = new Set(['HARPSN']);

function lblOptions(fetchOnly) {
  const opts = { target: $('lbl-target').value.trim(), root: $('root').value.trim(),
    lbl_dir: $('lbl-data').value.trim(), lbl_python: $('lbl-python').value.trim(),
    nproc: $('lbl-nproc').value, teff: $('lbl-teff').value, tellu: $('lbl-tellu').checked,
    fetch_only: !!fetchOnly };
  for (const key of lblTicked) opts[`lbl_${key}`] = true;
  return opts;
}

async function lblLoad() {
  const q = new URLSearchParams({ target: $('lbl-target').value.trim(), root: $('root').value.trim(),
    data: $('lbl-data').value.trim() });
  try { lblNow = await api(`/api/lbl?${q}`); } catch (err) { $('lbl-pystate').textContent = err.message; return; }
  if (!$('lbl-data').value.trim()) $('lbl-data').value = lblNow.data;
  if (!$('lbl-python').value.trim()) $('lbl-python').value = lblNow.python;
  if (!$('lbl-nproc').value) $('lbl-nproc').value = lblNow.nproc;
  const ck = lblNow.check;
  $('lbl-pystate').innerHTML = ck.ok ? `<span class="ok">✓</span> LBL ${esc(ck.version)} (${esc(ck.where)})`
    : `<span class="bad">${esc(t('lbl_nopython'))}: ${esc(ck.error)}</span>`;
  lblDraw();
  lblCommand();
}

function lblDraw() {
  if (!lblNow) return;
  const rows = lblNow.instruments.map((one) => {
    const have = lblHave[one.key];
    const inArchive = !have ? '<span class="hint">?</span>' : have.error ? `<span class="bad">${esc(have.error)}</span>`
      : `${have.n} <span class="hint">(~${have.mb} MB)</span>`;
    const result = one.rdb ? `${one.npoints} ${esc(t('points'))} <span class="hint">(${esc(one.made)})</span> `
      + `<button type="button" class="small" data-lbluse="${esc(one.key)}">${esc(t('lbl_use'))}</button>` : '<span class="hint">-</span>';
    return `<tr><td><label><input type="checkbox" data-lbltick="${esc(one.key)}"${lblTicked.has(one.key) ? ' checked' : ''}> ${esc(one.label)}</label></td>`
      + `<td class="hint">${esc((TEXT[lang] || {})[`lbl_archive_${one.key}`] || one.archive)}</td><td class="num">${inArchive}</td><td class="num">${one.there || '-'}</td><td>${result}</td></tr>`;
  }).join('');
  $('lbl-insts').innerHTML = `<table class="mini"><tr><th>${esc(t('lbl_inst'))}</th><th>${esc(t('lbl_from'))}</th>`
    + `<th>${esc(t('lbl_inarchive'))}</th><th>${esc(t('lbl_here'))}</th><th>${esc(t('lbl_result'))}</th></tr>${rows}</table>`;
}

async function lblCommand() {
  try {
    const res = await api('/api/command', { action: 'lbl', options: lblOptions(false) });
    $('cmd-lbl').textContent = res.line;
  } catch (err) { $('cmd-lbl').textContent = err.message; }
}

async function lblCheck() {
  const target = $('lbl-target').value.trim();
  if (!target) { $('lbl-checkstate').textContent = t('lbl_nostar'); return; }
  $('lbl-checkstate').innerHTML = `<span class="spin"></span> ${esc(t('lbl_asking'))}`;
  try {
    lblHave = await api('/api/lbl/check', { target, root: $('root').value.trim() });
    $('lbl-checkstate').textContent = '';
  } catch (err) { $('lbl-checkstate').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
  lblDraw();
}

async function lblSave() {
  try {
    lblNow = await api('/api/lbl/settings', { python: $('lbl-python').value.trim(), data: $('lbl-data').value.trim(),
      nproc: $('lbl-nproc').value, target: $('lbl-target').value.trim(), root: $('root').value.trim() });
  } catch (err) { $('lbl-pystate').textContent = err.message; return; }
  lblLoad();
}

async function lblRun(fetchOnly) {
  $('lbl-runstate').textContent = '';
  try {
    const state = await api('/api/run', { action: 'lbl', options: lblOptions(fetchOnly) });
    keep(state);
    lblJob();
  } catch (err) { $('lbl-runstate').innerHTML = `<span class="bad">${esc(err.message)}</span>`; }
}

// the run of LBL, the last one: its steps, the spectra as they come, its log
function lblJob() {
  const mine = [...jobs.values()].filter((job) => job.action === 'lbl');
  const job = mine[mine.length - 1];
  $('lbl-jobcard').hidden = !job;
  if (!job) return;
  const running = job.status === 'running';
  const steps = job.steps.map((st, i) => {
    const last = i === job.steps.length - 1;
    const icon = last && running ? '<span class="hourglass">⏳</span>' : last && job.status !== 'done'
      ? '<span class="bad">✗</span>' : '<span class="ok">✓</span>';
    return `<li><span class="icon">${icon}</span><span>${esc(st.name)}</span><span class="time">${clock(st.elapsed)}</span>`
      + `${st.result ? `<span class="result">${esc(st.result)}</span>` : ''}</li>`;
  }).join('');
  const bars = Object.entries(job.rows || {}).filter(([, row]) => running && row.total).map(([, row]) => {
    const pct = Math.round(100 * row.done / row.total);
    return `<p><span class="pbar"><span style="width:${pct}%"></span></span><span class="ptext">${esc(t('lbl_spectra'))}: ${row.done} / ${row.total}</span></p>`;
  }).join('');
  $('lbl-job').innerHTML = `<div class="job"><div class="job-head"><span class="what">${esc(t('lbl_job'))}</span>`
    + `<span class="status ${job.status}">${esc(t(job.status))}</span><span class="when">${clock(job.elapsed)}</span><span class="spacer"></span>`
    + `${running ? `<button class="stop" type="button" data-stop="${job.id}">${esc(t('stop'))}</button>` : ''}</div>`
    + `<pre class="codeblock">${esc(job.line)}</pre><ol class="steps">${steps}</ol>${bars}`
    + `<details open><summary>${esc(t('log'))} (${job.lines.length})</summary><pre class="codeblock log">${esc(job.lines.slice(-200).join('\n'))}</pre></details></div>`;
  const pre = $('lbl-job').querySelector('pre.log');
  if (pre) pre.scrollTop = pre.scrollHeight;
}

// the rdb of LBL given to the analysis: a file of the star, its dataset named
//   so that koloa knows its spectrograph (HARPS-N (LBL) and the CCFs of the
//   TNG are two releases of the same spectra)
function lblUse(key) {
  const one = lblNow && lblNow.instruments.find((inst) => inst.key === key);
  if (!one || !one.rdb) return;
  fileRows = fileRows.filter((row) => row.path.trim());
  if (!fileRows.some((row) => row.path === one.rdb)) fileRows.push({ path: one.rdb, label: one.dataset });
  renderFiles();
  if (!$('target').value.trim()) $('target').value = $('lbl-target').value.trim();
  showTab('analysis');
  updateCommands();
  checkArchives();
}

// the tab shown: the star of the analysis, unless one is typed here
function lblShow() {
  if (!$('lbl-target').value.trim() && $('target').value.trim()) $('lbl-target').value = $('target').value.trim();
  lblLoad();
  lblJob();
}
window.lblShow = lblShow;

$('lbl-check').addEventListener('click', lblCheck);
$('lbl-run').addEventListener('click', () => lblRun(false));
$('lbl-fetch').addEventListener('click', () => lblRun(true));
$('lbl-target').addEventListener('change', () => { lblHave = {}; lblLoad(); });
$('lbl-target').addEventListener('keydown', (e) => { if (e.key === 'Enter') { lblHave = {}; lblLoad(); lblCheck(); } });
for (const id of ['lbl-data', 'lbl-python', 'lbl-nproc']) $(id).addEventListener('change', lblSave);
for (const id of ['lbl-teff', 'lbl-tellu']) $(id).addEventListener('change', lblCommand);
document.addEventListener('change', (e) => {
  const box = e.target.closest('[data-lbltick]');
  if (!box) return;
  if (box.checked) lblTicked.add(box.dataset.lbltick); else lblTicked.delete(box.dataset.lbltick);
  lblCommand();
});
document.addEventListener('click', (e) => {
  const use = e.target.closest('[data-lbluse]');
  if (use) lblUse(use.dataset.lbluse);
});
// the run followed while the tab is open; the results read again when it ends
setInterval(() => { if (!$('tab-lbl').hidden) lblJob(); }, 1500);
jobDoneHooks.push((job) => { if (job.action === 'lbl') lblLoad(); });
langHooks.push(() => { lblDraw(); lblJob(); });

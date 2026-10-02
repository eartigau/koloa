// koloa's GUI, the batch page: a list of targets, the options of all, and
// the bash script that makes all their reports (koloa.gui serves it)
'use strict';

const KEPT = 'koloa-batch';
let targets = [{ target: '', files: [] }];
let made = null;

// -----------------------------------------------------------------------------
// what is kept in this browser
// -----------------------------------------------------------------------------
function runSettings() {
  return {
    mode: document.querySelector('input[name="mode"]:checked').value,
    host: $('host').value.trim(), workdir: $('workdir').value.trim(),
    koloa: $('koloa').value.trim(), batch: $('batchdir').value.trim(),
    jobs: $('jobs').value, bashrc: $('bashrc').checked,
  };
}

function save() {
  try {
    localStorage.setItem(KEPT, JSON.stringify({ targets, run: runSettings(), options: readOptions('batch') }));
  } catch (err) { /* no storage: nothing kept */ }
}

function restore() {
  let kept = null;
  try { kept = JSON.parse(localStorage.getItem(KEPT) || 'null'); } catch (err) { kept = null; }
  if (!kept) return;
  if (Array.isArray(kept.targets) && kept.targets.length) targets = kept.targets;
  const run = kept.run || {};
  const radio = document.querySelector(`input[name="mode"][value="${run.mode || 'local'}"]`);
  if (radio) radio.checked = true;
  for (const [id, key] of [['host', 'host'], ['workdir', 'workdir'], ['koloa', 'koloa'], ['batchdir', 'batch'], ['jobs', 'jobs']]) {
    if (run[key] !== undefined && run[key] !== null) $(id).value = run[key];
  }
  if (run.bashrc !== undefined) $('bashrc').checked = !!run.bashrc;
  for (const [key, val] of Object.entries(kept.options || {})) {
    const el = document.querySelector(`[data-for="batch"][data-opt="${key}"]`);
    if (!el) continue;
    if (el.type === 'checkbox') el.checked = !!val; else el.value = val;
  }
}

// -----------------------------------------------------------------------------
// the list of targets
// -----------------------------------------------------------------------------
function renderTargets() {
  $('targetlist').innerHTML = `<div class="trow thead"><span>#</span><span>${esc(t('name_col'))} ${info('target')}</span>`
    + `<span class="flabel-head">${esc(t('files_col'))} ${info('files')}<span class="spacer"></span>${esc(t('inst_head'))} ${info('file_label')}<span class="pad"></span></span><span></span></div>`
    + targets.map((tg, i) => {
      const files = (tg.files.length ? tg.files : []).map((f, j) => `<div class="tfile">
          <span class="picker"><input type="text" class="tpath" data-t="${i}" data-f="${j}" value="${esc(f.path)}" autocomplete="off">
          <button type="button" data-browse="${i}" data-f="${j}">${esc(t('browse'))}</button></span>
          <input type="text" class="tlabel" data-t="${i}" data-f="${j}" value="${esc(f.label || '')}" autocomplete="off">
          <button type="button" class="small" data-delfile="${i}" data-f="${j}">&times;</button></div>`).join('');
      return `<div class="trow"><span class="tnum">${i + 1}</span>
        <span><input type="text" class="tname" data-t="${i}" value="${esc(tg.target)}" autocomplete="off"></span>
        <span class="tfiles">${files}<button type="button" class="small" data-addfile="${i}">${esc(t('add_file_short'))}</button></span>
        <button type="button" class="small" data-deltarget="${i}">&times;</button></div>`;
    }).join('');
  const n = targets.filter((tg) => tg.target.trim() || tg.files.some((f) => f.path.trim())).length;
  $('ntargets').textContent = `(${n})`;
}

// one line of the pasted list: NAME | file1 file2=LABEL
function parseLine(text) {
  const [name, rest] = text.includes('|') ? [text.slice(0, text.indexOf('|')), text.slice(text.indexOf('|') + 1)] : [text, ''];
  const files = rest.trim().split(/\s+/).filter(Boolean).map((word) => {
    const cut = word.lastIndexOf('=');
    return cut > 0 && !word.slice(cut + 1).includes('/') ? { path: word.slice(0, cut), label: word.slice(cut + 1) } : { path: word, label: '' };
  });
  return { target: name.trim(), files };
}

// -----------------------------------------------------------------------------
// where it runs
// -----------------------------------------------------------------------------
function showMode() {
  const ssh = runSettings().mode === 'ssh';
  document.querySelectorAll('.sshonly').forEach((el) => { el.style.display = ssh ? '' : 'none'; });
  document.querySelectorAll('.localonly').forEach((el) => { el.style.display = ssh ? 'none' : ''; });
}

// a file (or a folder): a dialog of this machine, or the browser of the page
//   on the server
async function pickPath(kind, start) {
  const run = runSettings();
  if (run.mode === 'ssh') {
    if (!run.host) { alert(t('need_host')); return null; }
    return browse({ host: run.host, start: start || run.workdir, kind });
  }
  try {
    const res = await api(`/api/pick?${new URLSearchParams({ kind, start: start || run.workdir || '' })}`);
    return res.path || null;
  } catch (err) {
    // no dialog on this machine (a server seen through a tunnel): the page's
    return browse({ start: start || run.workdir, kind });
  }
}

// -----------------------------------------------------------------------------
// the script
// -----------------------------------------------------------------------------
async function makeScript() {
  try {
    made = await api('/api/batch_script', { targets, options: readOptions('batch'), run: runSettings() });
    $('script').textContent = made.script;
    $('make-note').textContent = `${made.targets.length} ${t('n_targets')}: ${made.batch}/${made.targets.join(', ')}`;
    const run = runSettings();
    const where = run.workdir || '~';
    $('follow').textContent = run.mode === 'ssh'
      ? `ssh ${run.host} 'cd ${where} && ${made.start}'\nssh ${run.host} tail -f ${where}/${made.name.replace(/\.sh$/, '.out')}`
      : '';
    $('script-note').textContent = '';
    return made;
  } catch (err) {
    $('make-note').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
    return null;
  }
}

function download() {
  if (!made) return;
  const link = document.createElement('a');
  link.href = URL.createObjectURL(new Blob([made.script], { type: 'text/x-sh' }));
  link.download = made.name;
  link.click();
  URL.revokeObjectURL(link.href);
}

// -----------------------------------------------------------------------------
// the page
// -----------------------------------------------------------------------------
document.addEventListener('input', (e) => {
  const el = e.target;
  if (el.classList.contains('tname')) targets[+el.dataset.t].target = el.value;
  if (el.classList.contains('tpath')) targets[+el.dataset.t].files[+el.dataset.f].path = el.value;
  if (el.classList.contains('tlabel')) targets[+el.dataset.t].files[+el.dataset.f].label = el.value;
  if (el.classList.contains('tname') || el.classList.contains('tpath')) {
    const n = targets.filter((tg) => tg.target.trim() || tg.files.some((f) => f.path.trim())).length;
    $('ntargets').textContent = `(${n})`;
  }
  save();
});
document.addEventListener('change', (e) => {
  if (e.target.name === 'mode') showMode();
  save();
});
document.addEventListener('click', async (e) => {
  const el = e.target;
  if (el.id === 'addtarget') { targets.push({ target: '', files: [] }); renderTargets(); save(); }
  if (el.id === 'cleartargets' && confirm(t('confirm_clear'))) { targets = [{ target: '', files: [] }]; renderTargets(); save(); }
  if (el.dataset.addfile !== undefined) { targets[+el.dataset.addfile].files.push({ path: '', label: '' }); renderTargets(); save(); }
  if (el.dataset.delfile !== undefined) { targets[+el.dataset.delfile].files.splice(+el.dataset.f, 1); renderTargets(); save(); }
  if (el.dataset.deltarget !== undefined) {
    targets.splice(+el.dataset.deltarget, 1);
    if (!targets.length) targets.push({ target: '', files: [] });
    renderTargets(); save();
  }
  if (el.dataset.browse !== undefined) {
    const file = targets[+el.dataset.browse].files[+el.dataset.f];
    const path = await pickPath('file', file.path ? file.path.replace(/\/[^/]*$/, '') : '');
    if (path) { file.path = path; renderTargets(); save(); }
  }
  if (el.id === 'pick-workdir') {
    const path = await pickPath('folder', $('workdir').value);
    if (path) { $('workdir').value = path; save(); }
  }
  if (el.id === 'addpasted') {
    const rows = $('paste').value.split('\n').map((line) => line.trim()).filter(Boolean).map(parseLine);
    targets = targets.filter((tg) => tg.target.trim() || tg.files.some((f) => f.path.trim())).concat(rows);
    if (!targets.length) targets.push({ target: '', files: [] });
    $('paste').value = '';
    renderTargets(); save();
  }
  if (el.id === 'ssh-test') {
    $('ssh-result').innerHTML = `<span class="spin"></span> ${esc(t('testing'))}`;
    try {
      const res = await api(`/api/ssh_test?${new URLSearchParams({ host: $('host').value.trim(), koloa: $('koloa').value.trim() })}`);
      $('ssh-result').innerHTML = `<span class="${res.ok ? 'ok' : 'bad'}">${esc(res.text)}</span>`;
    } catch (err) {
      $('ssh-result').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
    }
  }
  if (el.id === 'make') await makeScript();
  if (el.id === 'download') { if (made || await makeScript()) download(); }
  if (el.id === 'runhere') {
    if (!(await makeScript())) return;
    try {
      keep(await api('/api/batch_run', { workdir: runSettings().workdir, name: made.name, script: made.script, batch: made.batch }));
      renderJobs();
      $('jobs').scrollIntoView({ behavior: 'smooth', block: 'start' });
    } catch (err) { alert(err.message); }
  }
  if (el.id === 'send' || el.id === 'start') {
    if (!(await makeScript())) return;
    const run = runSettings();
    if (el.id === 'start' && !confirm(`${t('confirm_start')} ${run.host}?`)) return;
    try {
      const sent = await api('/api/batch_send', { host: run.host, workdir: run.workdir, name: made.name, script: made.script });
      let note = `${t('sent')} ${run.host}:${sent.path}`;
      if (el.id === 'start') {
        await api('/api/batch_start', { host: run.host, workdir: run.workdir, name: made.name });
        note += `; ${t('started')}`;
      }
      $('script-note').innerHTML = `<span class="ok">${esc(note)}</span>`;
    } catch (err) {
      $('script-note').innerHTML = `<span class="bad">${esc(err.message)}</span>`;
    }
  }
});

(async () => {
  $('batch-options').innerHTML = renderOptions('batch');
  restore();
  langHooks.push(renderTargets, renderJobs);
  applyLang();
  showMode();
  renderTargets();
  try { (await api('/api/jobs')).filter((job) => job.action === 'batch').forEach(keep); } catch (err) { /* nothing yet */ }
  renderJobs();
  setInterval(poll, 2000);
})();

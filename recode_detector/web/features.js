/* Saved runs, exports, quality checks and the raw-spectrum viewer. Loaded after app.js and dbcompare.js. */

/* ---------- export links ---------- */
const exLink = (name, fmt, label) => `<a href="/api/download/${name}?fmt=${fmt}" download>${label || fmt.toUpperCase()}</a>`;
function fillExportBars() {
  $$('[data-export]').forEach(el => {
    const n = el.dataset.export;
    el.innerHTML = exLink(n, 'csv') + exLink(n, 'xlsx', 'Excel');
  });
}
const Q_CARDS = [
  ['results', 'Differential results', 'Every tested protein: fold change, p-value, FDR and call.'],
  ['processed', 'Processed matrix', 'Normalised, imputed log2 intensities.'],
  ['qc_checks', 'Quality checks', 'Each check with its result and rule.'],
  ['sample_qc', 'Per-sample QC', 'Proteins, missing values, correlation and flags per sample.'],
  ['enrichment', 'Enrichment', 'Gene-set results (after running enrichment).'],
];
const D_CARDS = [
  ['db_sites', 'Sites', 'One row per differing position with its verdict, peptides and detectability.'],
  ['db_peptides', 'Peptides', 'Every peptide that spans a site, with start, end, charge and intensity per run.'],
  ['db_coverage', 'Coverage', 'Every reference protein: length, peptides, percent covered.'],
  ['db_qc_checks', 'Quality checks', 'Each check with its result and rule.'],
  ['db_qc_runs', 'Per-run QC', 'Precursors, peptides and signal for each run.'],
  ['db_spectra', 'Spectrum checks', 'One row per checked peptide: MS1 peak, fragments matched, chance probability, verdict.'],
  ['db_spectra_ions', 'Matched fragment ions', 'Every matched b and y ion with its mass error.'],
];
function renderExports() {
  const q = $('#x-grid'), d = $('#dx-grid');
  const card = (n, t, desc, fmt, only) => `<div class="card dl-card"><div class="kicker">${only || fmt.toUpperCase()}</div><h3>${t}</h3><p>${desc}</p>
    <span class="xp">${only ? exLink(n, fmt, 'Download') : exLink(n, fmt) + exLink(n, 'tsv') + exLink(n, 'xlsx', 'Excel')}</span></div>`;
  if (q) {
    const f = $('#fmt-q').value;
    q.innerHTML = Q_CARDS.map(c => `<div class="card dl-card"><div class="kicker">TABLE</div><h3>${c[1]}</h3><p>${c[2]}</p><span class="xp">${exLink(c[0], f)}</span></div>`).join('') +
      `<div class="card dl-card"><div class="kicker">JSON</div><h3>Analysis settings</h3><p>Groups, parameters and each processing step, for your methods section.</p><span class="xp">${exLink('settings', 'json', 'Download')}</span></div>` +
      `<div class="card dl-card"><div class="kicker">ZIP</div><h3>Everything</h3><p>All tables above plus the settings file.</p><span class="xp">${exLink('zip_quant', f, 'Download ZIP')}</span></div>`;
  }
  if (d) {
    const f = $('#fmt-d').value;
    d.innerHTML = `<div class="card dl-card"><div class="kicker">PDF</div><h3>Final report</h3><p>The whole analysis explained in plain English, with figures, results, quality checks, spectrum evidence, caveats and a glossary.</p><span class="xp">${exLink('db_report', 'pdf', 'Download PDF')}</span></div>` +
      D_CARDS.map(c => `<div class="card dl-card"><div class="kicker">TABLE</div><h3>${c[1]}</h3><p>${c[2]}</p><span class="xp">${exLink(c[0], f)}</span></div>`).join('') +
      `<div class="card dl-card"><div class="kicker">ZIP</div><h3>Everything</h3><p>All tables plus the PDF report.</p><span class="xp">${exLink('zip_db', f, 'Download ZIP')}</span></div>`;
  }
}
$('#fmt-q').onchange = $('#fmt-d').onchange = renderExports;

/* ---------- quality-check rows ---------- */
const CHIP = { pass: 'Pass', warn: 'Check', fail: 'Problem', info: 'Info' };
function renderChecks(el, checks) {
  el.innerHTML = (checks || []).map(c => `<div class="chk"><span class="chip s-${c.status}">${CHIP[c.status]}</span><span>${esc(c.name)}<br><b>${esc(c.value)}</b></span><span class="rule">${esc(c.rule)}</span><span class="why">${esc(c.why)}</span></div>`).join('') || '<p class="small">No checks yet.</p>';
}
function summarise(checks) {
  const n = k => checks.filter(c => c.status === k).length;
  return `${n('pass')} pass · ${n('warn')} to look at · ${n('fail')} problem${n('fail') === 1 ? '' : 's'}`;
}
window.renderQuantQc = function () {
  if (!S.checks) return;
  renderChecks($('#q-checks'), S.checks); $('#q-sum').textContent = summarise(S.checks);
  $('#t-sq tbody').innerHTML = (S.sampleQc || []).map(r => `<tr><td>${esc(r.sample)}</td><td>${esc(r.group)}</td><td class="num">${r.proteins_identified}</td><td class="num">${r.missing_pct}</td><td class="num">${r.median_log2}</td><td class="num">${r.correlation_to_group ?? '–'}</td><td>${r.flag ? `<span class="chip s-warn">${esc(r.flag)}</span>` : ''}</td></tr>`).join('');
};
window.drawDbQc = function () {
  const r = DBS.result; if (!r || !r.qc || !r.qc.checks) return;
  const sm = r.qc.summary;
  $('#dq-read').innerHTML = ro('overall', sm.overall === 'pass' ? 'Pass' : sm.overall === 'warn' ? 'Look closer' : 'Problem') + ro('pass', sm.pass) + ro('to look at', sm.warn) + ro('problems', sm.fail);
  renderChecks($('#dq-checks'), r.qc.checks);
  $('#t-dqr tbody').innerHTML = r.qc.runs.map(x => `<tr><td>${esc(x.run)}</td><td class="num">${x.precursors.toLocaleString()}</td><td class="num">${x.peptides.toLocaleString()}</td><td class="num">${x.total_intensity.toExponential(2)}</td></tr>`).join('');
};

/* ---------- saved runs ---------- */
const fmtDate = t => new Date(t * 1000).toLocaleString(undefined, { day: 'numeric', month: 'short', year: 'numeric', hour: 'numeric', minute: '2-digit' });
const fmtMB = b => b > 1e6 ? (b / 1e6).toFixed(1) + ' MB' : Math.max(1, Math.round(b / 1e3)) + ' kB';
function runSummary(r) {
  const s = r.summary || {};
  if (r.kind === 'db') return `${s.sites ?? '?'} sites: ${s.alternative ?? 0} alternative, ${s.reference ?? 0} reference, ${s.both ?? 0} both, ${s.none ?? 0} not covered`;
  return `${s.samples ?? '?'} samples, ${(s.proteins_kept ?? 0).toLocaleString()} proteins kept${s.groups ? ', groups ' + s.groups.join(', ') : ''}`;
}
async function loadRuns() {
  try {
    const d = await api('/api/runs');
    $('#runs-root').textContent = 'Stored in ' + d.root;
    $('#runs-list').innerHTML = d.runs.length ? d.runs.map(r => `<div class="run" data-id="${r.id}">
      <div><input class="name" value="${esc(r.name)}" aria-label="Run name">
        <div class="meta"><span class="badge b-${r.kind === 'db' ? 'alternative' : 'reference'}">${r.kind === 'db' ? 'database check' : 'quantification'}</span>
          ${r.has_spectra ? '<span class="badge b-both">spectra</span>' : ''} ${esc(runSummary(r))}</div>
        <div class="meta">Saved ${fmtDate(r.updated)} · ${Object.values(r.files).map(esc).join(', ')} · ${fmtMB(r.disk_bytes)}</div>
        <textarea placeholder="Notes about this run…" aria-label="Notes">${esc(r.note || '')}</textarea></div>
      <div class="acts"><button class="btn" data-act="open">Open</button><button class="ghost" data-act="del">Delete</button></div></div>`).join('')
      : '<p class="small" style="padding:10px 0">No saved runs yet. Run an analysis and it appears here automatically.</p>';
  } catch (e) { toast(e.message); }
}
window.loadRuns = loadRuns;
$('#runs-list').addEventListener('click', async e => {
  const b = e.target.closest('[data-act]'); if (!b) return;
  const id = b.closest('.run').dataset.id;
  if (b.dataset.act === 'del') {
    if (!confirm('Delete this saved run? Its saved input files are removed. This cannot be undone.')) return;
    try { await api('/api/runs/delete?id=' + id, {}); loadRuns(); } catch (x) { toast(x.message); }
  } else busy(b, async () => { applyOpenedRun(await api('/api/runs/open?id=' + id, {})); });
});
$('#runs-list').addEventListener('change', async e => {
  const id = e.target.closest('.run')?.dataset.id; if (!id) return;
  try {
    if (e.target.classList.contains('name')) await api('/api/runs/rename', { id, name: e.target.value });
    else if (e.target.tagName === 'TEXTAREA') await api('/api/runs/note', { id, note: e.target.value });
    toast('Saved', true);
  } catch (x) { toast(x.message); }
});

function applyOpenedRun(r) {
  if (r.mode === 'quant') {
    S.st = r.state; S.comp = null; S.qc = r.analysis.qc; S.checks = r.analysis.checks; S.sampleQc = r.analysis.sample_qc; S.levels = r.analysis.groups;
    S.st.analysed = true; S.st.loaded = true;
    const p = (r.params && r.params.processing) || {};
    if (p.normalization) $('#pr-norm').value = p.normalization; if (p.imputation) $('#pr-imp').value = p.imputation;
    if (p.min_valid) $('#pr-min').value = p.min_valid; if (p.valid_rule) $('#pr-rule').value = p.valid_rule;
    $('#pr-log').textContent = (r.analysis.log || []).join('\n');
    const opts = r.analysis.groups.map(g => `<option>${esc(g)}</option>`).join('');
    $('#d-a').innerHTML = opts; $('#d-b').innerHTML = opts; $('#d-a').selectedIndex = Math.max(0, r.analysis.groups.length - 1);
    afterLoad(true);
    if (r.comparison) {
      S.comp = r.comparison; const m = r.comparison.meta;
      $('#d-a').value = m.a; $('#d-b').value = m.b; $('#d-method').value = m.method; $('#d-use').value = m.use; $('#d-fdr').value = m.fdr; $('#d-lfc').value = m.lfc;
      $('#d-out').classList.remove('hidden'); S.sel = null;
    }
    setMode('quant'); updateNav(); show(r.comparison ? 'diff' : 'qc');
    toast('Opened ' + S.st.filename, true);
  } else {
    DBS.st = r.state; DBS.result = r.state.result; $('#db-marker').value = (r.params && r.params.marker) || $('#db-marker').value;
    if (window.afterDbChange) afterDbChange(true);
    renderDb(); setMode('db'); updateNav(); show('dbsites');
    toast('Opened the saved run', true);
  }
}

/* ---------- spectra ---------- */
const SP = { poll: null, sel: null, sites: null };
const stored = k => { try { return localStorage.getItem(k) || ''; } catch (e) { return ''; } };
const verdictWord = v => ({ supported: 'supported', weak: 'weak', 'not seen': 'not seen', 'no MS1 signal': 'no MS1 signal', 'no MS/MS': 'no MS/MS', error: 'error' }[v] || v);
const vcls = v => v === 'supported' ? 'supported' : v === 'weak' ? 'weak' : 'no';
const specData = () => (DBS.result && DBS.result.spectra) || null;

$('#sp-pick').onclick = () => busy($('#sp-pick'), async () => { const r = await api('/api/pick', {}); $('#sp-path').value = r.path; try { localStorage.setItem('rd-spectra-path', r.path); } catch (e) {} });
$('#sp-path').addEventListener('change', e => { try { localStorage.setItem('rd-spectra-path', e.target.value); } catch (x) {} });
$('#sp-run').onclick = () => startSpectra(null);
$('#sp-cancel').onclick = () => api('/api/db/spectra/cancel', {}).catch(() => {});
$('#sp-one').onclick = async () => {
  const path = $('#sp-path').value.trim(); if (!path) return toast('Choose the raw data file first.');
  try { await api('/api/db/spectra/peptide', { path, sequence: $('#sp-seq').value, charge: +$('#sp-z').value, ppm: +$('#sp-ppm').value, half_window: +$('#sp-win').value }); watchJob(); } catch (e) { toast(e.message); }
};
async function startSpectra(sites) {
  const path = $('#sp-path').value.trim(); if (!path) { show('dbspectra'); $('#sp-path').focus(); return toast('Choose the raw data file first.'); }
  try { await api('/api/db/spectra/start', { path, sites, ppm: +$('#sp-ppm').value, half_window: +$('#sp-win').value }); show('dbspectra'); watchJob(); } catch (e) { toast(e.message); }
}
document.addEventListener('click', e => { const b = e.target.closest('[data-check-site]'); if (b) startSpectra([+b.dataset.checkSite]); });

function watchJob() {
  clearInterval(SP.poll);
  const tick = async () => {
    let st; try { st = await api('/api/db/spectra/status'); } catch (e) { return; }
    const j = st.job; if (!j) return;
    $('#sp-prog').classList.toggle('hidden', !(j.status === 'running' || j.status === 'queued'));
    $('#sp-cancel').hidden = !(j.status === 'running' || j.status === 'queued');
    $('#sp-prog .progbar i').style.width = (j.total ? 100 * j.done / j.total : 5) + '%';
    $('#sp-prog span').textContent = j.message || '';
    if (['done', 'error', 'cancelled'].includes(j.status)) {
      clearInterval(SP.poll);
      if (j.status === 'error') toast(j.error); else if (j.status === 'done') toast('Spectrum check finished', true);
      if (st.spectra && DBS.result) { DBS.result.spectra = st.spectra; renderSpectra(); }
    }
  };
  SP.poll = setInterval(tick, 1200); tick();
}

window.renderSpectra = function () {
  const r = DBS.result; if (!r) return;
  if (!$('#sp-path').value) $('#sp-path').value = (r.spectra && r.spectra.path) || stored('rd-spectra-path');
  const st = DBS.st || {};
  $('#sp-note').textContent = st.raw_reader === false ? 'Thermo .raw files need ThermoRawFileParser (installed by MassSpec Bench). Without it, convert the file to mzML first.' : '';
  const sp = specData();
  $('#sp-out').classList.toggle('hidden', !sp || !(sp.results.length || (sp.custom || []).length));
  if (!sp) return;
  const res = sp.results, n = k => res.filter(x => x.verdict === k).length;
  $('#sp-read').innerHTML = ro('peptides checked', res.length) + ro('supported', n('supported'), 'up') + ro('weak', n('weak')) + ro('not seen', res.filter(x => !['supported', 'weak'].includes(x.verdict)).length, 'down') +
    (sp.source ? ro('instrument', sp.source.instrument || sp.source.name || 'file') : '');
  const sites = r.sites;
  $('#t-sp tbody').innerHTML = res.map((x, i) => {
    const s = x.site != null ? sites[x.site] : null, bm = (x.ms2 || {}).best || {};
    return `<tr data-i="${i}" class="${SP.sel === i ? 'sel' : ''}" style="cursor:pointer"><td>${s ? esc(shortId(s.ref_protein)) : ''}</td><td class="num">${s ? s.position : ''}</td>
      <td><span class="badge b-${x.version === 'alternative' ? 'alternative' : 'reference'}">${esc(x.version)}</span></td><td class="mono">${esc(x.modified)}</td><td class="num">${x.charge}+</td>
      <td class="num">${x.ms1 && x.ms1.apex_rt != null ? x.ms1.apex_rt.toFixed(2) : '–'}</td><td class="num">${bm.n_matched ?? '–'} (${bm.n_disc ?? '–'})</td>
      <td class="num">${bm.p_chance != null ? (bm.p_chance <= 1e-12 ? '<1e-12' : bm.p_chance.toExponential(0)) : '–'}</td><td><span class="chip s-${vcls(x.verdict)}">${esc(verdictWord(x.verdict))}</span></td></tr>`;
  }).join('');
  if (SP.sel != null && res[SP.sel]) renderSpecDetail(); else if (res.length && SP.sel == null) { SP.sel = 0; renderSpecDetail(); }
  const cu = sp.custom || [];
  $('#sp-custom-card').hidden = !cu.length;
  if (cu.length) { $('#sp-custom-detail').innerHTML = ''; cu.slice(0, 2).forEach(x => specColumn($('#sp-custom-detail'), x, 'c' + cu.indexOf(x))); }
};
$('#t-sp tbody').addEventListener('click', e => { const tr = e.target.closest('tr'); if (!tr) return; SP.sel = +tr.dataset.i; $$('#t-sp tbody tr').forEach(t => t.classList.toggle('sel', t === tr)); renderSpecDetail(); });
function renderSpecDetail() {
  const sp = specData(); if (!sp || SP.sel == null) return;
  const x = sp.results[SP.sel], box = $('#sp-detail'); box.innerHTML = '';
  const partner = sp.results.find((y, i) => i !== SP.sel && y.pair && y.pair === x.pair);
  specColumn(box, x, 'a'); if (partner) specColumn(box, partner, 'b');
}
function specColumn(box, x, tag) {
  const col = document.createElement('div'); col.className = 'card spec-col';
  const bm = (x.ms2 || {}).best, id = 'sp' + tag + Math.random().toString(36).slice(2, 6);
  col.innerHTML = `<h4>${esc(x.modified)} <span class="chip s-${vcls(x.verdict)}">${esc(verdictWord(x.verdict))}</span> <span class="small">${esc(x.version)}${x.charge ? ', ' + x.charge + '+' : ''}${x.precursor_mz ? ', m/z ' + x.precursor_mz : ''}</span></h4>
    <div class="small">${(x.reasons || []).map(esc).join('<br>')}</div>
    <div class="kicker">Intact peptide (MS1)</div><canvas id="${id}m" height="150"></canvas>
    <div class="kicker">Fragment spectrum${bm ? ` (scan ${bm.scan}, ${bm.rt.toFixed(2)} min)` : ''}</div><canvas id="${id}s" height="230"></canvas>
    <div class="kicker">Fragment ions over time</div><canvas id="${id}x" height="170"></canvas>`;
  box.appendChild(col);
  if (x.ms1 && x.ms1.rt.length) drawMs1(col.querySelector('#' + id + 'm'), x); else blank(col.querySelector('#' + id + 'm'), 'no MS1 signal at this mass');
  if (bm) { drawSpec(col.querySelector('#' + id + 's'), x); drawFragXic(col.querySelector('#' + id + 'x'), x); }
  else { blank(col.querySelector('#' + id + 's'), 'no fragmentation scan'); blank(col.querySelector('#' + id + 'x'), ''); }
}
function blank(cv, msg) { const { c, w } = setup(cv); c.fillStyle = MUTE; c.textAlign = 'center'; c.fillText(msg, w / 2, 40); }
function drawMs1(cv, x) {
  const m = x.ms1, fr = frame(cv, { xr: ext(m.rt, 0), yr: [0, Math.max(...m.intensity, 1) * 1.1], xl: 'retention time (min)', box: [60, 8, 12, 34] });
  fr.c.beginPath(); m.rt.forEach((t, i) => { const px = fr.X(t), py = fr.Y(m.intensity[i]); i ? fr.c.lineTo(px, py) : fr.c.moveTo(px, py); });
  fr.c.strokeStyle = TH.a; fr.c.lineWidth = 1.6; fr.c.stroke();
  fr.c.strokeStyle = TH.b; fr.c.setLineDash([4, 4]); fr.c.beginPath(); fr.c.moveTo(fr.X(m.apex_rt), fr.T); fr.c.lineTo(fr.X(m.apex_rt), fr.h - fr.B); fr.c.stroke(); fr.c.setLineDash([]);
  fr.c.fillStyle = INK; fr.c.textAlign = 'left'; fr.c.fillText(`apex ${m.apex_rt.toFixed(2)} min`, fr.X(m.apex_rt) + 6, fr.T + 12);
}
function drawSpec(cv, x) {
  const b = x.ms2.best, mx = Math.max(...b.int, 1), lo = Math.min(...b.mz) - 8, hi = Math.max(...b.mz) + 8;
  const fr = frame(cv, { xr: [lo, hi], yr: [0, 1.12], xl: 'm/z', box: [44, 8, 12, 34] }), hit = new Map(b.matches.map(m => [m.index, m]));
  b.mz.forEach((m, i) => { if (!hit.has(i)) { fr.c.strokeStyle = TH.ns; fr.c.lineWidth = 1; fr.c.beginPath(); fr.c.moveTo(fr.X(m), fr.Y(0)); fr.c.lineTo(fr.X(m), fr.Y(b.int[i] / mx)); fr.c.stroke(); } });
  fr.c.font = '10px ui-monospace,Menlo,monospace'; fr.c.textAlign = 'center';
  b.matches.forEach(m => {
    const px = fr.X(m.mz), py = fr.Y(m.int / mx); fr.c.strokeStyle = m.disc ? TH.b : TH.a; fr.c.lineWidth = 2.4; fr.c.beginPath(); fr.c.moveTo(px, fr.Y(0)); fr.c.lineTo(px, py); fr.c.stroke();
    fr.c.fillStyle = m.disc ? TH.b : TH.a; fr.c.fillText(m.label, px, py - 4);
  });
  fr.c.textAlign = 'right'; fr.c.fillStyle = TH.b; fr.c.fillText('● contains the changed residue', fr.w - 16, fr.T + 14); fr.c.fillStyle = TH.a; fr.c.fillText('● other matched ions', fr.w - 16, fr.T + 28);
}
function drawFragXic(cv, x) {
  const tr = x.ms2.traces; if (!tr.length) return blank(cv, 'no matched fragments');
  const all = tr.flatMap(t => t.intensity), rt = tr[0].rt, fr = frame(cv, { xr: ext(rt, 0), yr: [0, Math.max(...all, 1) * 1.1], xl: 'retention time (min)', box: [60, 8, 12, 34] });
  tr.forEach((t, k) => { fr.c.beginPath(); t.rt.forEach((v, i) => { const px = fr.X(v), py = fr.Y(t.intensity[i]); i ? fr.c.lineTo(px, py) : fr.c.moveTo(px, py); }); fr.c.strokeStyle = PAL[k % PAL.length]; fr.c.lineWidth = t.disc ? 2.2 : 1.2; fr.c.stroke(); });
  fr.c.font = '10px ui-monospace,Menlo,monospace'; fr.c.textAlign = 'left'; tr.forEach((t, k) => { fr.c.fillStyle = PAL[k % PAL.length]; fr.c.fillText(t.label, fr.L + 6 + k * 46, fr.T + 12); });
  if (x.ms2.coelution != null) { fr.c.fillStyle = MUTE; fr.c.textAlign = 'right'; fr.c.fillText('co-elution r = ' + x.ms2.coelution.toFixed(2), fr.w - 14, fr.T + 12); }
}
const _redraw = window.onLookChange;
window.onLookChange = () => { _redraw && _redraw(); if (S.page === 'dbspectra') renderSpectra(); };

/* ---------- start-up ---------- */
fillExportBars(); renderExports();
(async () => {
  try {
    const st = await api('/api/db/state');
    if (st.job && ['running', 'queued'].includes(st.job.status)) watchJob();
  } catch (e) {}
})();

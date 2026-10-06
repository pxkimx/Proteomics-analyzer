/* Recode Detector front end. No dependencies; plots are drawn on canvases. */
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const S = { mode: 'quant', last: { quant: 'data', db: 'dbfiles' }, st: null, qc: null, comp: null, sel: null, boxWhich: 'raw', pc: [0, 1], gmt: '', sort: { k: 'p_value', dir: 1 }, hover: null };
const DBS = { st: null, result: null };
let PAL, INK, MUTE, GRID;
function syncColors() { PAL = TH.pal; INK = TH.ink; MUTE = TH.mute; GRID = TH.grid; }
syncColors();
const lerp3 = (stops, t) => { t = Math.max(0, Math.min(1, t)) * (stops.length - 1); const i = Math.min(stops.length - 2, Math.floor(t)), f = t - i; return `rgb(${stops[i].map((v, j) => Math.round(v + (stops[i + 1][j] - v) * f)).join(',')})`; };
const heatCol = t => lerp3(TH.heat, t);                                                  // low -> high
const divCol = z => z >= 0 ? lerp3([TH.div[0], TH.div[1]], z) : lerp3([TH.div[0], TH.div[2]], -z);
const esc = s => String(s ?? '').replace(/[&<>"']/g, m => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[m]));
const fmt = (v, d = 2) => v == null ? '–' : (Math.abs(v) < 0.001 && v !== 0 ? v.toExponential(1) : (+v).toFixed(d));

/* ---------- plumbing ---------- */
function toast(msg, ok) {
  const t = $('#toast'); t.textContent = msg; t.className = 'show' + (ok ? ' ok' : '');
  clearTimeout(toast.t); toast.t = setTimeout(() => t.className = '', ok ? 2500 : 7000);
}
async function api(path, body, raw) {
  const opt = body === undefined ? {} : { method: 'POST', body: raw ? body : JSON.stringify(body) };
  const r = await fetch(path, opt); const j = await r.json();
  if (!r.ok) throw new Error(j.error || r.statusText);
  return j;
}
async function busy(btn, fn) {
  const label = btn && btn.textContent; if (btn) { btn.disabled = true; btn.textContent = 'Working…'; }
  try { return await fn(); } catch (e) { toast(e.message); } finally { if (btn) { btn.disabled = false; btn.textContent = label; } }
}
const levels = () => [...new Set(Object.values(S.st.groups).filter(Boolean))].sort();
const gcolor = g => PAL[levels().indexOf(g) % PAL.length] || PAL[7];

/* ---------- navigation ---------- */
const order = ['data', 'groups', 'process', 'qc', 'diff', 'enrich', 'export'];
function reach() {
  const st = S.st, r = DBS.result;
  return { data: true, groups: st.loaded, process: st.loaded, qc: st.analysed, diff: st.analysed, enrich: !!S.comp, export: st.analysed,
    runs: true, dbfiles: true, dbsites: !!r, dbspectra: !!r, dbcover: !!r, dbqc: !!r, dbscatter: !!(r && r.summary.has_b), dbexport: !!r };
}
function updateNav() {
  const r = reach();
  $$('#nav button').forEach(b => { b.hidden = b.dataset.mode !== S.mode && b.dataset.mode !== 'both'; b.disabled = !r[b.dataset.page]; });
  $$('#mode button').forEach(b => b.classList.toggle('on', b.dataset.mode === S.mode));
}
function setMode(m) {
  S.mode = m; updateNav();
  show(S.last[m]);
}
$('#mode').addEventListener('click', e => { const b = e.target.closest('button'); if (b && b.dataset.mode !== S.mode) setMode(b.dataset.mode); });
function show(page) {
  if (!reach()[page]) return;
  $$('.page').forEach(p => p.classList.toggle('on', p.id === 'p-' + page));
  $$('#nav button').forEach(b => b.classList.toggle('on', b.dataset.page === page));
  S.page = page; if (page !== 'runs') S.last[page.startsWith('db') ? 'db' : 'quant'] = page; scrollTo(0, 0); requestAnimationFrame(() => redraw(page));
}
$('#nav').addEventListener('click', e => { const b = e.target.closest('button'); if (b) show(b.dataset.page); });
document.addEventListener('click', e => { const g = e.target.closest('[data-go]'); if (g) show(g.dataset.go); });

/* ---------- 01 data ---------- */
function wireDrop(el, onFile) {
  ['dragenter', 'dragover'].forEach(ev => el.addEventListener(ev, e => { e.preventDefault(); el.classList.add('over'); }));
  ['dragleave', 'drop'].forEach(ev => el.addEventListener(ev, e => { e.preventDefault(); el.classList.remove('over'); }));
  el.addEventListener('drop', e => { const f = e.dataTransfer.files[0]; if (f) onFile(f); });
}
async function loadFile(f, kind) {
  S.file = f; S.comp = null; S.qc = null;
  toast('Reading ' + f.name + '…', true);
  const url = '/api/load?name=' + encodeURIComponent(f.name) + (kind ? '&kind=' + encodeURIComponent(kind) : '');
  try { S.st = await api(url, f, true); afterLoad(); } catch (e) { toast(e.message); }
}
$('#browse').onclick = () => $('#file').click();
$('#file').onchange = e => e.target.files[0] && loadFile(e.target.files[0]);
wireDrop($('#drop'), loadFile);
$('#example').onclick = ev => busy(ev.target, async () => { S.comp = null; S.qc = null; S.st = await api('/api/example', {}); afterLoad(); });
$('#kind').onchange = e => busy(null, async () => { S.st = await api('/api/kind?kind=' + encodeURIComponent(e.target.value), {}); afterLoad(true); });

function ro(label, val, cls = '') { return `<div class="ro ${cls}"><b>${esc(val)}</b><span>${esc(label)}</span></div>`; }
function afterLoad(keepPage) {
  const st = S.st;
  $('#loaded').classList.remove('hidden');
  $('#l-name').textContent = st.filename + (st.is_example ? '  (simulated)' : '');
  $('#l-read').innerHTML = ro('proteins', st.n_proteins.toLocaleString()) + ro('samples', st.samples.length) +
    ro('missing values', st.missing_pct.toFixed(1) + '%') + ro('format', st.format);
  const removed = Object.entries(st.removed || {}).map(([k, v]) => `${v.toLocaleString()} ${k} removed`);
  $('#l-notes').innerHTML = [...st.notes, ...removed].map(n => `<li>${esc(n)}</li>`).join('');
  const k = $('#kindrow'); k.classList.toggle('hidden', !(st.kinds && st.kinds.length > 1));
  if (st.kinds) $('#kind').innerHTML = st.kinds.map(x => `<option ${st.notes[0].includes('“' + x + '”') ? 'selected' : ''}>${esc(x)}</option>`).join('');
  buildGroups(); updateNav(); if (!keepPage) toast('Loaded ' + st.n_proteins.toLocaleString() + ' proteins × ' + st.samples.length + ' samples', true);
}

/* ---------- 02 groups ---------- */
function buildGroups() {
  const st = S.st;
  $('#gtable tbody').innerHTML = st.samples.map((s, i) =>
    `<tr><td>${esc(s)}</td><td><input list="glist" data-i="${i}" value="${esc(st.groups[s] || '')}" placeholder="(excluded)"></td></tr>`).join('');
  refreshGlist();
}
function refreshGlist() { $('#glist').innerHTML = levels().map(g => `<option value="${esc(g)}">`).join(''); }
$('#gtable').addEventListener('input', e => { const i = e.target.dataset.i; if (i != null) { S.st.groups[S.st.samples[i]] = e.target.value.trim(); refreshGlist(); } });
$('#guess').onclick = () => busy($('#guess'), async () => {
  const st = S.st; st.groups = guessGroups(st.samples); buildGroups();
});
function guessGroups(samples) {
  const out = {}; samples.forEach(s => out[s] = s.replace(/[\s_\-.]*(rep(licate)?|r|n|bio|tech)?[\s_\-.]*\d+$/i, '').replace(/^[ _\-.]+|[ _\-.]+$/g, '') || s);
  return new Set(Object.values(out)).size === samples.length ? Object.fromEntries(samples.map(s => [s, ''])) : out;
}
$('#gsave').onclick = () => busy($('#gsave'), async () => {
  const lv = levels(); if (lv.length < 2) throw new Error('Give at least two different group names.');
  const small = lv.filter(g => Object.values(S.st.groups).filter(x => x === g).length < 2);
  if (small.length) throw new Error('Each group needs at least two samples. Too few in: ' + small.join(', '));
  S.st = await api('/api/groups', { groups: S.st.groups }); S.comp = null; S.qc = null; updateNav(); show('process');
});

/* ---------- 03 process ---------- */
const syncImp = () => $('#pr-ds').classList.toggle('hidden', $('#pr-imp').value !== 'downshift');
$('#pr-imp').onchange = syncImp; syncImp();
$('#run').onclick = () => busy($('#run'), async () => {
  const params = { normalization: $('#pr-norm').value, imputation: $('#pr-imp').value, min_valid: +$('#pr-min').value,
    valid_rule: $('#pr-rule').value, downshift: +$('#pr-shift').value, width: +$('#pr-width').value, seed: +$('#pr-seed').value };
  const r = await api('/api/analyze', { params });
  S.qc = r.qc; S.checks = r.checks; S.sampleQc = r.sample_qc; S.comp = null; S.st.analysed = true; S.levels = r.groups;
  if (r.run_id) toast('Saved to Runs', true);
  $('#pr-log').textContent = r.log.join('\n');
  const opts = r.groups.map(g => `<option>${esc(g)}</option>`).join('');
  $('#d-a').innerHTML = opts; $('#d-b').innerHTML = opts; if (r.groups.length > 1) $('#d-a').selectedIndex = r.groups.length - 1;
  $('#d-out').classList.add('hidden'); updateNav(); show('qc');
});

/* ---------- plot toolkit ---------- */
function setup(canvas) {
  if (!canvas.dataset.h) canvas.dataset.h = canvas.getAttribute('height');
  const d = devicePixelRatio || 1, w = canvas.clientWidth, h = +canvas.dataset.h;
  canvas.width = w * d; canvas.height = h * d; canvas.style.height = h + 'px';
  const c = canvas.getContext('2d'); c.setTransform(d, 0, 0, d, 0, 0); c.clearRect(0, 0, w, h);
  c.font = '11px ui-monospace,Menlo,monospace'; return { c, w, h };
}
function ticks(lo, hi, n = 6) {
  const span = hi - lo || 1, raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) || raw, out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}
function frame(canvas, { xr, yr, box = [54, 12, 14, 38], xl, yl, noX, noY } = {}) {
  const { c, w, h } = setup(canvas), [L, T, R, B] = box;
  const X = v => L + (v - xr[0]) / (xr[1] - xr[0]) * (w - L - R), Y = v => h - B - (v - yr[0]) / (yr[1] - yr[0]) * (h - T - B);
  c.fillStyle = MUTE; c.strokeStyle = GRID; c.lineWidth = 1; c.textAlign = 'right'; c.textBaseline = 'middle';
  if (!noY) for (const v of ticks(...yr, 5)) { c.beginPath(); c.moveTo(L, Y(v)); c.lineTo(w - R, Y(v)); c.stroke(); c.fillText(+v.toFixed(2), L - 6, Y(v)); }
  c.textAlign = 'center'; c.textBaseline = 'top';
  if (!noX) for (const v of ticks(...xr, 7)) { c.beginPath(); c.moveTo(X(v), T); c.lineTo(X(v), h - B); c.stroke(); c.fillText(+v.toFixed(2), X(v), h - B + 5); }
  c.strokeStyle = TH.border; c.strokeRect(L, T, w - L - R, h - T - B);
  c.fillStyle = INK;
  if (xl) { c.textAlign = 'center'; c.fillText(xl, (L + w - R) / 2, h - 15); }
  if (yl) { c.save(); c.translate(13, (T + h - B) / 2); c.rotate(-Math.PI / 2); c.textAlign = 'center'; c.textBaseline = 'middle'; c.fillText(yl, 0, 0); c.restore(); }
  return { c, w, h, X, Y, L, T, R, B };
}
const ext = (arr, pad = 0.05) => { let lo = Infinity, hi = -Infinity; for (const v of arr) if (v != null && isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } const p = (hi - lo || 1) * pad; return [lo - p, hi + p]; };
function wrapText(c, txt, max) { while (c.measureText(txt).width > max && txt.length > 3) txt = txt.slice(0, -2); return txt; }
function xLabels(c, labels, X0, step, y, rot) {
  c.fillStyle = MUTE; c.textBaseline = 'top';
  labels.forEach((l, i) => { c.save(); c.translate(X0 + i * step, y); c.rotate(rot); c.textAlign = 'right'; c.fillText(wrapText(c, l, 90), 0, 0); c.restore(); });
}

/* ---------- 04 qc plots ---------- */
function drawIdent() {
  const q = S.qc, n = q.per_sample.length, cvs = $('#c-ident'), { c, w, h } = setup(cvs);
  const L = 54, B = 78, T = 12, R = 14, mx = Math.max(...q.per_sample.map(s => s.identified)) * 1.08;
  const fr = frame(cvs, { xr: [0, n], yr: [0, mx], box: [L, T, R, B], noX: true });
  const bw = (w - L - R) / n;
  q.per_sample.forEach((s, i) => {
    const x = L + i * bw + bw * .18, ww = bw * .64, col = gcolor(s.group);
    fr.c.globalAlpha = .35; fr.c.fillStyle = col; fr.c.fillRect(x, fr.Y(s.identified), ww, fr.Y(0) - fr.Y(s.identified));
    fr.c.globalAlpha = 1; fr.c.fillRect(x, fr.Y(s.after_filter), ww, fr.Y(0) - fr.Y(s.after_filter));
  });
  xLabels(fr.c, q.per_sample.map(s => s.sample), L + bw * .8, bw, h - B + 8, -Math.PI / 4);
  fr.c.fillStyle = MUTE; fr.c.textAlign = 'left'; fr.c.fillText('pale = identified, solid = kept after filter', L + 6, T + 4);
}
function drawBox() {
  const q = S.qc, key = S.boxWhich === 'raw' ? 'box_raw' : 'box_norm', n = q.per_sample.length, cvs = $('#c-box');
  const all = q.per_sample.flatMap(s => [s.box_raw.min, s.box_raw.max, s.box_norm.min, s.box_norm.max]);
  const L = 54, B = 78, T = 12, R = 14, fr = frame(cvs, { xr: [0, n], yr: ext(all, .04), box: [L, T, R, B], noX: true, yl: 'log2 intensity' });
  const bw = (fr.w - L - R) / n;
  q.per_sample.forEach((s, i) => {
    const b = s[key], cx = L + (i + .5) * bw, ww = bw * .52, col = gcolor(s.group), c = fr.c;
    c.strokeStyle = col; c.fillStyle = col + '55'; c.lineWidth = 1.5;
    c.beginPath(); c.moveTo(cx, fr.Y(b.min)); c.lineTo(cx, fr.Y(b.q1)); c.moveTo(cx, fr.Y(b.q3)); c.lineTo(cx, fr.Y(b.max)); c.stroke();
    c.fillRect(cx - ww / 2, fr.Y(b.q3), ww, fr.Y(b.q1) - fr.Y(b.q3)); c.strokeRect(cx - ww / 2, fr.Y(b.q3), ww, fr.Y(b.q1) - fr.Y(b.q3));
    c.lineWidth = 2.5; c.beginPath(); c.moveTo(cx - ww / 2, fr.Y(b.med)); c.lineTo(cx + ww / 2, fr.Y(b.med)); c.stroke();
  });
  xLabels(fr.c, q.per_sample.map(s => s.sample), L + bw * .8, bw, fr.h - B + 8, -Math.PI / 4);
}
function drawPca() {
  const q = S.qc, [a, b] = S.pc, sc = q.pca.scores, cvs = $('#c-pca');
  const xs = sc.map(r => r[a]), ys = sc.map(r => r[b]), fr = frame(cvs, { xr: ext(xs, .22), yr: ext(ys, .22), box: [54, 12, 14, 38], xl: `PC${a + 1} (${q.pca.explained[a].toFixed(1)}%)`, yl: `PC${b + 1} (${q.pca.explained[b].toFixed(1)}%)` });
  q.samples.forEach((s, i) => {
    const g = q.per_sample[i].group, x = fr.X(xs[i]), y = fr.Y(ys[i]);
    fr.c.fillStyle = gcolor(g); fr.c.beginPath(); fr.c.arc(x, y, 6.5, 0, 6.3); fr.c.fill();
    fr.c.strokeStyle = TH.bg1; fr.c.lineWidth = 1.5; fr.c.stroke();
    fr.c.fillStyle = INK; fr.c.textAlign = 'left'; fr.c.textBaseline = 'middle'; fr.c.fillText(s, x + 10, y);
  });
  let lx = fr.L + 10; levels().forEach(g => { fr.c.fillStyle = gcolor(g); fr.c.fillRect(lx, fr.T + 8, 10, 10); fr.c.fillStyle = INK; fr.c.textAlign = 'left'; fr.c.fillText(g, lx + 15, fr.T + 14); lx += 30 + g.length * 7; });
}
function drawCorr() {
  const q = S.qc, cvs = $('#c-corr'), { c, w, h } = setup(cvs), n = q.samples.length, ord = q.corr.order, M = q.corr.matrix;
  const L = 96, T = 10, size = Math.min(w - L - 70, h - T - 96), cell = size / n;
  let lo = 1; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (i !== j) lo = Math.min(lo, M[i][j]);
  lo = Math.min(lo, .99);
  ord.forEach((si, r) => ord.forEach((sj, k) => {
    const v = M[si][sj], t = (v - lo) / (1 - lo);
    c.fillStyle = heatCol(t); c.fillRect(L + k * cell, T + r * cell, cell - 1, cell - 1);
    if (cell >= 34) { c.fillStyle = t > .55 ? TH.bg1 : INK; c.textAlign = 'center'; c.textBaseline = 'middle'; c.fillText(v.toFixed(2), L + (k + .5) * cell, T + (r + .5) * cell); }
  }));
  c.fillStyle = MUTE; c.textBaseline = 'middle'; c.textAlign = 'right';
  ord.forEach((si, r) => c.fillText(wrapText(c, q.samples[si], L - 8), L - 6, T + (r + .5) * cell));
  ord.forEach((si, k) => { c.save(); c.translate(L + (k + .5) * cell, T + size + 8); c.rotate(Math.PI / 3); c.textAlign = 'left'; c.fillText(wrapText(c, q.samples[si], 80), 0, 0); c.restore(); });
  const gx = L + size + 18; for (let i = 0; i < 60; i++) { const t = 1 - i / 59; c.fillStyle = heatCol(t); c.fillRect(gx, T + i * size / 60, 12, size / 60 + 1); }
  c.fillStyle = MUTE; c.textAlign = 'left'; c.fillText('1.00', gx + 17, T + 5); c.fillText(lo.toFixed(2), gx + 17, T + size - 5);
}
function drawCv() {
  const q = S.qc, gs = Object.keys(q.cv), cvs = $('#c-cv');
  const mx = Math.max(1, ...gs.flatMap(g => q.cv[g].hist)), fr = frame(cvs, { xr: [0, 100], yr: [0, mx * 1.1], xl: 'CV (%), proteins with ≥3 values', yl: 'proteins' });
  gs.forEach(g => {
    const h = q.cv[g].hist, c = fr.c; c.strokeStyle = gcolor(g); c.lineWidth = 2; c.beginPath();
    h.forEach((v, i) => { const x = fr.X(i * 5 + 2.5), y = fr.Y(v); i ? c.lineTo(x, y) : c.moveTo(x, y); }); c.stroke();
    c.fillStyle = gcolor(g) + '22'; c.lineTo(fr.X(97.5), fr.Y(0)); c.lineTo(fr.X(2.5), fr.Y(0)); c.fill();
    c.strokeStyle = gcolor(g); c.setLineDash([4, 4]); c.beginPath(); c.moveTo(fr.X(q.cv[g].median), fr.T); c.lineTo(fr.X(q.cv[g].median), fr.h - fr.B); c.stroke(); c.setLineDash([]);
  });
  let ly = fr.T + 14; gs.forEach(g => { fr.c.fillStyle = gcolor(g); fr.c.textAlign = 'right'; fr.c.fillText(`${g}: median ${q.cv[g].median.toFixed(1)}%`, fr.w - fr.R - 8, ly); ly += 15; });
}
function drawRank() {
  const q = S.qc, r = q.rank, fr = frame($('#c-rank'), { xr: [0, r.length * q.rank_step], yr: ext(r, .05), xl: 'protein rank', yl: 'mean log2' });
  fr.c.strokeStyle = TH.a; fr.c.lineWidth = 2; fr.c.beginPath(); r.forEach((v, i) => { const x = fr.X(i * q.rank_step), y = fr.Y(v); i ? fr.c.lineTo(x, y) : fr.c.moveTo(x, y); }); fr.c.stroke();
}
function drawQc() {
  const q = S.qc; if (!q) return;
  if (window.renderQuantQc) renderQuantQc();
  $('#qc-read').innerHTML = ro('proteins kept', q.n_proteins.toLocaleString()) + ro('samples', q.samples.length) + ro('missing after filter', q.missing_total_pct.toFixed(1) + '%') +
    ro('median CV', Object.values(q.cv).length ? (Object.values(q.cv).reduce((a, b) => a + b.median, 0) / Object.values(q.cv).length).toFixed(1) + '%' : '–');
  const np = q.pca.explained.length, pairs = []; for (let i = 0; i < np; i++) for (let j = i + 1; j < np; j++) pairs.push([i, j]);
  const sel = $('#pcsel'); if (sel.options.length !== pairs.length) sel.innerHTML = pairs.map(([i, j], k) => `<option value="${k}">PC${i + 1} vs PC${j + 1}</option>`).join('');
  S.pcPairs = pairs; S.pc = pairs[+sel.value || 0] || [0, 1];
  drawIdent(); drawBox(); drawPca(); drawCorr(); drawCv(); drawRank();
}
$('#pcsel').onchange = () => drawQc();
$('#boxwhich').addEventListener('click', e => { const b = e.target.closest('button'); if (!b) return; S.boxWhich = b.dataset.w; $$('#boxwhich button').forEach(x => x.classList.toggle('on', x === b)); drawBox(); });

/* ---------- 05 differential ---------- */
const colorOf = call => call === 'up' ? TH.up : call === 'down' ? TH.down : TH.ns;
function runComparison(btn) {
  return busy(btn, async () => {
    const a = $('#d-a').value, b = $('#d-b').value; if (a === b) throw new Error('Pick two different groups.');
    S.comp = await api('/api/compare', { a, b, method: $('#d-method').value, use: $('#d-use').value, fdr: +$('#d-fdr').value, lfc: +$('#d-lfc').value });
    S.sel = null; $('#d-out').classList.remove('hidden'); updateNav(); drawDiff();
  });
}
$('#d-run').onclick = e => runComparison(e.target);
['d-fdr', 'd-lfc'].forEach(id => $('#' + id).addEventListener('change', () => S.comp && runComparison(null)));
function drawDiff() {
  const m = S.comp; if (!m) return;
  $('#d-read').innerHTML = ro(`up in ${m.meta.a}`, m.counts.up, 'up') + ro(`up in ${m.meta.b}`, m.counts.down, 'down') + ro('tested', m.n_tested.toLocaleString()) +
    ro('prior df', m.meta.prior_df == null ? 'n/a' : (isFinite(m.meta.prior_df) ? m.meta.prior_df.toFixed(1) : '∞'));
  drawVolcano(); drawHeat(); drawSel(); fillTable();
}
function volcanoPts() {
  return S.comp.rows.filter(r => r.p_value != null && r.log2fc != null).map(r => ({ r, x: r.log2fc, y: -Math.log10(r.p_value) }));
}
function drawVolcano() {
  const pts = S.vp = volcanoPts(), m = S.comp.meta, cvs = $('#c-volc');
  const xm = Math.max(1.5, ...pts.map(p => Math.abs(p.x))) * 1.05, fr = S.vf = frame(cvs, { xr: [-xm, xm], yr: [0, Math.max(2, ...pts.map(p => p.y)) * 1.06], xl: `log2 fold change   (← higher in ${m.b}   |   higher in ${m.a} →)`, yl: '−log10 p' });
  const c = fr.c;
  c.strokeStyle = TH.cAlpha; c.setLineDash([5, 5]); [-m.lfc, m.lfc].forEach(v => { c.beginPath(); c.moveTo(fr.X(v), fr.T); c.lineTo(fr.X(v), fr.h - fr.B); c.stroke(); });
  const sigP = pts.filter(p => p.r.call !== 'ns').map(p => p.y); if (sigP.length) { const yv = Math.min(...sigP); c.beginPath(); c.moveTo(fr.L, fr.Y(yv)); c.lineTo(fr.w - fr.R, fr.Y(yv)); c.stroke(); } c.setLineDash([]);
  for (const pass of ['ns', 'sig']) for (const p of pts) {
    if ((p.r.call === 'ns') !== (pass === 'ns')) continue;
    c.fillStyle = colorOf(p.r.call); c.beginPath(); c.arc(fr.X(p.x), fr.Y(p.y), p.r.call === 'ns' ? 2.2 : 3.2, 0, 6.3); c.fill();
  }
  // label the strongest hits without overlap
  const placed = [], labelled = pts.filter(p => p.r.call !== 'ns').sort((a, b) => b.y - a.y).slice(0, 14);
  c.font = '11px ui-monospace,Menlo,monospace'; c.textBaseline = 'middle';
  for (const p of labelled) {
    const t = p.r.gene || p.r.protein, x = fr.X(p.x), y = fr.Y(p.y), wd = c.measureText(t).width, left = p.x < 0, bx = left ? x - wd - 8 : x + 7;
    if (placed.some(b => bx < b[0] + b[2] && bx + wd > b[0] && Math.abs(y - b[1]) < 12)) continue;
    placed.push([bx, y, wd]); c.fillStyle = INK; c.textAlign = 'left'; c.fillText(t, bx, y);
  }
  if (S.sel) { const p = pts.find(q => q.r.protein === S.sel); if (p) { c.strokeStyle = '#fff'; c.lineWidth = 2; c.beginPath(); c.arc(fr.X(p.x), fr.Y(p.y), 8, 0, 6.3); c.stroke(); } }
  const f = $('#d-find').value.trim().toUpperCase();
  if (f) for (const p of pts) if ((p.r.gene || '').toUpperCase() === f || p.r.protein.toUpperCase() === f) { c.strokeStyle = TH.c; c.lineWidth = 2.5; c.beginPath(); c.arc(fr.X(p.x), fr.Y(p.y), 9, 0, 6.3); c.stroke(); }
}
function nearest(e) {
  const cv = $('#c-volc'), b = cv.getBoundingClientRect(), mx = e.clientX - b.left, my = e.clientY - b.top; let best = null, bd = 144;
  for (const p of S.vp || []) { const dx = S.vf.X(p.x) - mx, dy = S.vf.Y(p.y) - my, d = dx * dx + dy * dy; if (d < bd) { bd = d; best = p; } }
  return { best, mx, my };
}
$('#c-volc').addEventListener('mousemove', e => {
  const { best, mx, my } = nearest(e), tip = $('#tip');
  if (!best) return tip.classList.add('hidden');
  const r = best.r; tip.innerHTML = `<b>${esc(r.gene || r.protein)}</b> ${esc(r.protein)}<br>log2 FC ${fmt(r.log2fc)} · FDR ${fmt(r.adj_p, 3)}`;
  tip.classList.remove('hidden'); tip.style.left = Math.min(mx + 14, $('#c-volc').clientWidth - 220) + 'px'; tip.style.top = (my + 14) + 'px';
});
$('#c-volc').addEventListener('mouseleave', () => $('#tip').classList.add('hidden'));
$('#c-volc').addEventListener('click', e => { const { best } = nearest(e); if (best) selectProtein(best.r.protein); });
$('#d-find').addEventListener('input', () => drawVolcano());

async function selectProtein(id) {
  S.sel = id; drawVolcano();
  try { S.selData = await api('/api/protein?id=' + encodeURIComponent(id)); drawSel(); } catch (e) { toast(e.message); }
}
function drawSel() {
  const d = S.selData, cvs = $('#c-sel');
  if (!S.sel || !d) { setup(cvs); $('#sel-title').textContent = 'Click a protein in the volcano'; $('#sel-info').textContent = ''; return; }
  const row = S.comp.rows.find(r => r.protein === S.sel);
  $('#sel-title').textContent = (row && row.gene) || S.sel; $('#sel-info').innerHTML = row ? `${esc(row.protein)}<br>log2 FC ${fmt(row.log2fc)} · p ${fmt(row.p_value, 3)} · FDR ${fmt(row.adj_p, 3)}<br>Open dots are imputed values.` : '';
  const gs = [...new Set(d.groups)], n = gs.length, fr = frame(cvs, { xr: [0, n], yr: ext(d.values, .12), box: [54, 12, 14, 36], noX: true, yl: 'log2 intensity' }), bw = (fr.w - fr.L - fr.R) / n;
  gs.forEach((g, gi) => {
    const idx = d.groups.map((x, i) => x === g ? i : -1).filter(i => i >= 0), cx = fr.L + (gi + .5) * bw;
    idx.forEach((i, k) => { const x = cx + (k - (idx.length - 1) / 2) * 14, y = fr.Y(d.values[i]); fr.c.beginPath(); fr.c.arc(x, y, 5, 0, 6.3); if (d.imputed[i]) { fr.c.strokeStyle = gcolor(g); fr.c.lineWidth = 2; fr.c.stroke(); } else { fr.c.fillStyle = gcolor(g); fr.c.fill(); } });
    const mean = idx.reduce((s, i) => s + d.values[i], 0) / idx.length; fr.c.strokeStyle = INK; fr.c.lineWidth = 2; fr.c.beginPath(); fr.c.moveTo(cx - bw * .3, fr.Y(mean)); fr.c.lineTo(cx + bw * .3, fr.Y(mean)); fr.c.stroke();
    fr.c.fillStyle = INK; fr.c.textAlign = 'center'; fr.c.textBaseline = 'top'; fr.c.fillText(g, cx, fr.h - fr.B + 8);
  });
}
function drawHeat() {
  const h = S.comp.heatmap, cvs = $('#c-heat');
  if (!h.z.length) { cvs.dataset.h = 60; const { c, w } = setup(cvs); c.fillStyle = MUTE; c.textAlign = 'left'; c.fillText('No significant proteins at these cut-offs.', 12, 30); return; }
  const rows = h.row_order.length, cols = h.col_order.length, L = 8, R = 140, T = 20, B = 70, rh = Math.max(3, Math.min(14, 520 / rows));
  cvs.dataset.h = Math.round(T + rows * rh + B);
  const { c, w } = setup(cvs), cw = (w - L - R) / cols, showL = rh >= 8;
  h.col_order.forEach((ci, k) => { c.fillStyle = gcolor(h.groups[ci]); c.fillRect(L + k * cw, 4, cw - 1, 10); });
  h.row_order.forEach((ri, r) => h.col_order.forEach((ci, k) => {
    const z = Math.max(-2, Math.min(2, h.z[ri][ci])) / 2, t = Math.abs(z);
    c.fillStyle = divCol(z);
    c.fillRect(L + k * cw, T + r * rh, cw - 1, rh - (rh >= 8 ? 1 : 0));
  }));
  if (showL) { c.fillStyle = INK; c.textAlign = 'left'; c.textBaseline = 'middle'; h.row_order.forEach((ri, r) => c.fillText(wrapText(c, h.labels[ri], R - 10), w - R + 6, T + (r + .5) * rh)); }
  c.fillStyle = MUTE; h.col_order.forEach((ci, k) => { c.save(); c.translate(L + (k + .5) * cw, T + rows * rh + 8); c.rotate(Math.PI / 4); c.textAlign = 'left'; c.textBaseline = 'middle'; c.fillText(wrapText(c, h.samples[ci], 80), 0, 0); c.restore(); });
}
function fillTable() {
  const m = S.comp, f = $('#t-filter').value, q = $('#t-search').value.trim().toUpperCase(), { k, dir } = S.sort;
  let rows = m.rows.filter(r => r.p_value != null && (f === 'all' || (f === 'sig' ? r.call !== 'ns' : r.call === f)) && (!q || (r.gene || '').toUpperCase().includes(q) || r.protein.toUpperCase().includes(q)));
  rows.sort((a, b) => ((a[k] ?? Infinity) > (b[k] ?? Infinity) ? 1 : -1) * dir);
  $('#t-res tbody').innerHTML = rows.slice(0, 400).map(r => `<tr data-id="${esc(r.protein)}"><td>${esc(r.gene)}</td><td>${esc(r.protein)}</td><td class="num ${r.call === 'ns' ? '' : r.call}">${fmt(r.log2fc)}</td><td class="num">${fmt(r.mean_a)}</td><td class="num">${fmt(r.mean_b)}</td><td class="num">${fmt(r.p_value, 3)}</td><td class="num">${fmt(r.adj_p, 3)}</td><td class="num">${r.n_a_observed}</td><td class="num">${r.n_b_observed}</td></tr>`).join('');
  $('#t-foot').textContent = `${Math.min(rows.length, 400)} of ${rows.length} shown. FC is ${m.meta.a} minus ${m.meta.b}. "Obs" counts values actually measured (the rest are imputed). The CSV export has everything.`;
}
$('#t-filter').onchange = $('#t-search').oninput = fillTable;
$('#t-res thead').addEventListener('click', e => { const th = e.target.closest('th'); if (!th) return; const k = th.dataset.k; S.sort = { k, dir: S.sort.k === k ? -S.sort.dir : 1 }; fillTable(); });
$('#t-res tbody').addEventListener('click', e => { const tr = e.target.closest('tr'); if (tr) { selectProtein(tr.dataset.id); scrollTo({ top: $('#c-volc').getBoundingClientRect().top + scrollY - 90, behavior: 'smooth' }); } });
document.addEventListener('click', e => { const b = e.target.closest('[data-save]'); if (!b) return; const a = document.createElement('a'); a.download = b.dataset.save.slice(2) + '.png'; a.href = $('#' + b.dataset.save).toDataURL('image/png'); a.click(); });

/* ---------- 06 enrichment ---------- */
async function loadGmt(f) { S.gmt = await f.text(); $('#gmtname').textContent = f.name + ' (' + S.gmt.split('\n').filter(Boolean).length + ' sets)'; }
$('#gmtbrowse').onclick = () => $('#gmtfile').click(); $('#gmtfile').onchange = e => e.target.files[0] && loadGmt(e.target.files[0]); wireDrop($('#gmtdrop'), loadGmt);
$('#en-run').onclick = e => busy(e.target, async () => {
  if (!S.gmt) throw new Error('Load a .gmt gene-set file first.');
  const r = await api('/api/enrich', { gmt: S.gmt, direction: $('#en-dir').value });
  $('#en-out').classList.remove('hidden'); $('#en-info').textContent = `${r.n_hits} proteins tested against ${r.n_sets} sets (sets need ≥3 members among the quantified proteins).`;
  $('#t-en tbody').innerHTML = r.rows.map(x => `<tr><td>${esc(x.set)}</td><td class="num">${x.set_size}</td><td class="num">${x.hits}</td><td class="num">${fmt(x.fold_enrichment, 1)}</td><td class="num">${fmt(x.p_value, 4)}</td><td class="num">${fmt(x.adj_p, 4)}</td><td>${esc(x.genes)}</td></tr>`).join('') || '<tr><td colspan="7">No set had enough members.</td></tr>';
});

/* ---------- redraw / lifecycle ---------- */
function redraw(page) { if (page === 'runs' && window.loadRuns) loadRuns(); if (page === 'export' || page === 'dbexport') window.renderExports && renderExports(); if (page === 'dbspectra') window.renderSpectra && renderSpectra(); if (page === 'dbqc') window.drawDbQc && drawDbQc(); if (page === 'qc') drawQc(); if (page === 'diff') drawDiff(); if (page === 'dbscatter' && window.drawDbScatter) drawDbScatter(); if (page === 'dbcover' && window.drawDbCover) drawDbCover(); }
let rt; addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => redraw(S.page), 120); });
$('#anim').onchange = e => window.setBackgroundAnimation(e.target.checked);
$('#quit').onclick = async () => { if (!confirm('Quit Recode Detector?')) return; try { await api('/api/quit', {}); } catch (e) {} document.body.innerHTML = '<p style="padding:40px;color:#8fb5b4;font:16px sans-serif">Recode Detector has quit. You can close this tab.</p>'; };

// Holding this stream open is how the server knows a window is open; closing the tab drops it.
function holdWindow() { const es = new EventSource('/api/window'); es.onerror = () => { es.close(); setTimeout(holdWindow, 1500); }; }
if (/[?&]shot=1/.test(location.search)) document.documentElement.classList.add('shot'); else holdWindow();       // ?shot=1: screenshots without holding the program open
(async () => {
  S.st = await api('/api/state'); $('#ver').textContent = 'v' + S.st.version;
  if (S.st.loaded) {
    afterLoad(true);
    if (S.st.analysed) {
      S.qc = S.st.qc; S.checks = S.st.checks; S.sampleQc = S.st.sample_qc; S.levels = S.st.levels; $('#pr-log').textContent = S.st.log.join('\n');
      const opts = S.levels.map(g => `<option>${esc(g)}</option>`).join('');
      $('#d-a').innerHTML = opts; $('#d-b').innerHTML = opts; $('#d-a').selectedIndex = S.levels.length - 1;
    }
  }
  try { DBS.st = await api('/api/db/state'); if (DBS.st.has_result) DBS.result = DBS.st.result; if (window.afterDbChange) afterDbChange(true); } catch (e) {}
  updateNav();
  const hp = new URLSearchParams(location.hash.slice(1)).get('page');
  if (hp && reach()[hp]) { S.mode = hp.startsWith('db') ? 'db' : 'quant'; updateNav(); show(hp); const sc = new URLSearchParams(location.hash.slice(1)).get('scroll'); if (sc) setTimeout(() => { const el = document.getElementById(sc); if (el) { el.scrollIntoView({ block: 'start' }); scrollBy(0, -100); } }, 900); if (new URLSearchParams(location.hash.slice(1)).get('open')) setTimeout(() => { const r = document.querySelector('#t-db tbody tr'); if (r) r.click(); }, 300); }
})();

window.onLookChange = () => { syncColors(); redraw(S.page); if (window.restartBackground) restartBackground(); };

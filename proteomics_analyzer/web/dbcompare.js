/* Database check mode: which version of a protein did the search see? */
const DB_SLOTS = [
  { slot: 'report', title: 'Peptide report from your search', need: true, hint: 'DIA-NN precursor matrix, report.pr_matrix.tsv', accept: '.tsv,.txt,.csv' },
  { slot: 'fasta_ref', title: 'Protein database (FASTA)', need: true, hint: 'The standard proteins, or one file holding standard and variant proteins', accept: '.fasta,.fa,.faa,.txt' },
  { slot: 'fasta_alt', title: 'Alternative-coding FASTA', need: false, hint: 'Optional. Only the variants, or the whole alternative proteome', accept: '.fasta,.fa,.faa,.txt' },
  { slot: 'report_b', title: 'Report from a standard-database search', need: false, hint: 'Optional. The same data searched without the alternative proteins', accept: '.tsv,.txt,.csv' },
];

function buildSlots() {
  const files = (DBS.st && DBS.st.files) || {};
  $('#db-slots').innerHTML = DB_SLOTS.map(s => {
    const f = files[s.slot];
    return `<div class="card slot" data-slot="${s.slot}"><div class="kicker">${s.need ? 'Required' : 'Optional'}</div><h3>${esc(s.title)}</h3>
      <div class="small">${esc(s.hint)}</div>
      ${f ? `<div class="got">✓ ${esc(f.name)} · ${esc(f.detail)}<button class="ghost sm" data-drop="${s.slot}">Remove</button></div>`
        : `<div class="drop slim" tabindex="0" data-zone="${s.slot}"><b>Drop file</b><span>or</span><button class="btn" data-pick="${s.slot}">Browse…</button></div>`}
      <input type="file" hidden accept="${s.accept}" data-file="${s.slot}"></div>`;
  }).join('');
  $$('#db-slots [data-zone]').forEach(z => wireDrop(z, f => dbUpload(z.dataset.zone, f)));
}
$('#db-slots').addEventListener('click', e => {
  const pick = e.target.closest('[data-pick]'), drop = e.target.closest('[data-drop]');
  if (pick) $(`#db-slots [data-file="${pick.dataset.pick}"]`).click();
  if (drop) busy(null, async () => { DBS.st = await api('/api/db/drop?slot=' + drop.dataset.drop, {}); DBS.result = null; afterDbChange(); });
});
$('#db-slots').addEventListener('change', e => { const i = e.target.closest('[data-file]'); if (i && i.files[0]) dbUpload(i.dataset.file, i.files[0]); });

async function dbUpload(slot, f) {
  toast('Reading ' + f.name + '…', true);
  try { DBS.st = await api(`/api/db/load?slot=${slot}&name=${encodeURIComponent(f.name)}`, f, true); DBS.result = null; afterDbChange(); }
  catch (e) { toast(e.message); }
}
function afterDbChange(quiet) {
  buildSlots(); updateNav();
  if (DBS.result) renderDb();
}

$('#db-example').onclick = e => busy(e.target, async () => {
  DBS.st = await api('/api/db/example', {}); DBS.result = null; afterDbChange();
  await runDb(null);
});
$('#db-run').onclick = e => runDb(e.target);
function runDb(btn) {
  return busy(btn, async () => {
    DBS.result = await api('/api/db/run', { marker: $('#db-marker').value });
    DBS.st.has_result = true; updateNav(); renderDb(); show('dbsites');
  });
}

/* ---------- sites ---------- */
const VERDICT = { alternative: 'alternative seen', reference: 'reference seen', both: 'both seen', none: 'no coverage' };
const fmtI = v => v >= 1e6 ? (v / 1e6).toFixed(1) + 'M' : v >= 1e3 ? (v / 1e3).toFixed(0) + 'k' : v.toFixed(0);
const shortId = id => { const t = id.split('|').filter(x => x && !/^(CANONICAL|CGG2W|ALT|VARIANT)$/i.test(x)); return t.length ? t[t.length - 1] : id; };

function renderDb() {
  const r = DBS.result; if (!r) return;
  const s = r.summary, c = s.status;
  $('#db-read').innerHTML = ro('differing positions', s.n_sites) + ro('alternative seen', c.alternative, 'up') + ro('reference seen', c.reference, 'down') +
    ro('both seen', c.both) + ro('no coverage', c.none) + ro('peptides in report', s.n_peptides.toLocaleString()) +
    (s.n_unpaired ? ro('unpaired variants', s.n_unpaired) : '');
  fillSites(); drawDbScatter();
}
function peptideLine(p, version, samples) {
  const mx = Math.max(...p.intensity, 1);
  const bars = p.intensity.map((v, i) => `<div class="bar2" title="${esc(samples[i])}: ${fmtI(v)}" style="width:${Math.max(2, v / mx * 100)}%"></div>`).join('');
  const hl = (p.sequence);
  return `<div class="pep ${version === 'alt' ? 'alt' : ''} ${p.also_elsewhere ? 'amb' : ''}"><span>${esc(hl)}</span>
    <em>${p.start}–${p.end} · z ${esc(p.charges.join('/') || '?')}</em><em>${p.missed_cleavages} missed</em>
    <span>${bars}<em>${fmtI(p.total_intensity)}${p.also_elsewhere ? ' · also elsewhere (not counted)' : ''}${p.seen_in_b ? ' · in standard search' : ''}</em></span></div>`;
}
function fillSites() {
  const r = DBS.result, f = $('#db-filter').value, q = $('#db-search').value.trim().toUpperCase();
  const rows = r.sites.filter(s => (f === 'all' || s.status === f) && (!q || s.ref_protein.toUpperCase().includes(q) || s.description.toUpperCase().includes(q)));
  const tb = $('#t-db tbody'); tb.innerHTML = '';
  rows.forEach((s, i) => {
    const ref = s.ref_peptides.filter(p => !p.also_elsewhere).map(p => p.sequence), alt = s.alt_peptides.filter(p => !p.also_elsewhere).map(p => p.sequence);
    const tr = document.createElement('tr'); tr.dataset.i = i; tr.style.cursor = 'pointer';
    tr.innerHTML = `<td title="${esc(s.description)}">${esc(shortId(s.ref_protein))}</td><td class="num">${s.position}</td><td class="num">${esc(s.ref_aa)} → ${esc(s.alt_aa)}</td>
      <td><span class="badge b-${s.status}">${VERDICT[s.status]}</span></td><td class="mono">${esc(ref.slice(0, 2).join(', '))}${ref.length > 2 ? ' +' + (ref.length - 2) : ''}</td>
      <td class="mono">${esc(alt.slice(0, 2).join(', '))}${alt.length > 2 ? ' +' + (alt.length - 2) : ''}</td>
      <td class="mono">${esc(s.context_ref.replace(/\[(.)\]/, '·$1·'))}<br>${esc(s.context_alt.replace(/\[(.)\]/, '·$1·'))}</td>`;
    tb.appendChild(tr);
    tr._site = s;
  });
  $('#db-foot').textContent = `${rows.length} of ${r.sites.length} sites. Click a row for the peptides, intensities and context.`;
}
$('#t-db tbody').addEventListener('click', e => {
  const tr = e.target.closest('tr'); if (!tr || tr.classList.contains('detail') || !tr._site) return;
  const nxt = tr.nextElementSibling;
  if (nxt && nxt.classList.contains('detail')) return nxt.remove();
  const s = tr._site, samples = DBS.result.samples, d = document.createElement('tr'); d.className = 'detail';
  const blk = (title, list, v) => `<div class="kicker" style="margin:8px 0 4px">${title}</div>` + (list.length ? list.map(p => peptideLine(p, v, samples)).join('') : '<div class="small">none identified</div>');
  d.innerHTML = `<td colspan="7"><div class="small">${esc(s.description)}</div>
    ${blk(`Reference peptides spanning ${s.ref_aa}${s.position}`, s.ref_peptides, 'ref')}${blk(`Alternative peptides spanning ${s.alt_aa}${s.position}`, s.alt_peptides, 'alt')}
    ${s.ref_peptides_b.length ? `<div class="kicker" style="margin:8px 0 4px">Reference peptides in the standard-database search</div>` + s.ref_peptides_b.map(p => `<div class="pep"><span>${esc(p.sequence)}</span><em>${p.start}–${p.end}</em><em></em><span><em>${fmtI(p.total_intensity)}${p.seen_in_a ? ' · also in main search' : ' · not in main search'}</em></span></div>`).join('') : ''}
    </td>`;
  tr.after(d);
});
$('#db-filter').onchange = $('#db-search').oninput = fillSites;

/* ---------- standard-search comparison ---------- */
function drawDbScatter() {
  const r = DBS.result; if (!r || !r.summary.has_b) return;
  const sc = r.scatter, cvs = $('#c-dbsc'); if (!cvs.clientWidth) return;
  $('#db-sread').innerHTML = ro('shared peptides', sc.n_shared.toLocaleString()) + ro('only in main search', sc.n_only_a.toLocaleString()) +
    ro('only in standard search', sc.n_only_b.toLocaleString()) + ro('log2 correlation', sc.r == null ? '–' : sc.r.toFixed(3));
  if (sc.a.length) {
    const lo = Math.min(...sc.a, ...sc.b) - .5, hi = Math.max(...sc.a, ...sc.b) + .5;
    const fr = frame(cvs, { xr: [lo, hi], yr: [lo, hi], xl: 'main search (with alternative proteins)', yl: 'standard-database search' });
    fr.c.strokeStyle = 'rgba(251,191,36,.5)'; fr.c.setLineDash([5, 5]); fr.c.beginPath(); fr.c.moveTo(fr.X(lo), fr.Y(lo)); fr.c.lineTo(fr.X(hi), fr.Y(hi)); fr.c.stroke(); fr.c.setLineDash([]);
    fr.c.fillStyle = 'rgba(45,212,191,.35)'; sc.a.forEach((x, i) => { fr.c.beginPath(); fr.c.arc(fr.X(x), fr.Y(sc.b[i]), 2.2, 0, 6.3); fr.c.fill(); });
  }
  const hits = r.sites.filter(s => s.ref_peptides_b.length);
  $('#db-bsites').innerHTML = hits.length ? hits.map(s => `<div style="margin:6px 0"><b>${esc(shortId(s.ref_protein))}</b> position ${s.position} <span class="badge b-${s.status}">${VERDICT[s.status]}</span><br><span class="mono">${esc(s.ref_peptides_b.map(p => p.sequence + (p.seen_in_a ? ' (also main)' : ' (main: absent)')).join(', '))}</span></div>`).join('')
    : 'No standard-search peptide spans any of the differing positions.';
}

buildSlots();

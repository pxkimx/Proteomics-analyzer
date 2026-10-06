"""Database comparison: which version of a protein did a search actually see?

Given a reference protein database and an alternative one (for example the same genome translated
with a codon read differently), find every position where a paired protein differs, then look at the
identified peptides that span each position. A peptide that spans the position can only come from one
of the two versions, so it is the evidence for that version. Peptides elsewhere in the protein cannot tell
the versions apart and are ignored.
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

MIN_PEP, MAX_PEP = 5, 70
SEQ_COLS = ("stripped.sequence", "peptide sequence", "peptide", "sequence", "stripped_sequence")
META_COLS = {"protein.group", "protein.ids", "protein.names", "genes", "proteotypic", "stripped.sequence",
             "modified.sequence", "precursor.charge", "precursor.id", "first.protein.description",
             "peptide sequence", "peptide", "sequence", "stripped_sequence", "protein", "charge"}


# ---------------------------------------------------------------- FASTA

def parse_fasta(text: str) -> dict[str, tuple[str, str]]:
    """id -> (full header line without '>', sequence). The id is the first whitespace-free word."""
    out: dict[str, tuple[str, str]] = {}
    head, seq = None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if head is not None:
                out[head.split()[0] if head.split() else head] = (head, "".join(seq).replace("*", "").upper())
            head, seq = line[1:].strip(), []
        elif head is not None and line:
            seq.append(line)
    if head is not None:
        out[head.split()[0] if head.split() else head] = (head, "".join(seq).replace("*", "").upper())
    return out


def _tokens(pid: str, marker: str) -> set[str]:
    t = {x for x in re.split(r"[|;,\s]+", pid) if x}
    if marker and marker in pid:
        t.add(pid.replace(marker, "").strip("_|-:. "))
        t.update(x.strip("_|-:. ") for x in pid.split(marker) if x.strip("_|-:. "))
    return {x for x in t if x}


@dataclass
class Pair:
    ref_id: str
    alt_id: str
    ref_header: str
    ref_seq: str
    alt_seq: str
    positions: list[int]               # 0-based


@dataclass
class Pairing:
    pairs: list[Pair]
    unpaired: list[str] = field(default_factory=list)       # variant ids with no equal-length partner
    identical: int = 0                                      # alt entries equal to their reference
    ref_pool: dict[str, str] = field(default_factory=dict)
    alt_pool: dict[str, str] = field(default_factory=dict)
    ref_headers: dict[str, str] = field(default_factory=dict)


def pair_databases(ref: dict, alt: dict | None = None, marker: str = "") -> Pairing:
    """Pair every alternative protein with its reference counterpart.

    `ref` and `alt` are parse_fasta() outputs. If `alt` is None, `ref` is a combined file and the entries
    whose id contains `marker` are the alternatives.
    """
    if alt is None:
        if not marker:
            raise ValueError("One FASTA file was given, so say which entries are the alternatives "
                             "(text that appears in their IDs, for example CGG2W).")
        alts = {k: v for k, v in ref.items() if marker in k}
        refs = {k: v for k, v in ref.items() if marker not in k}
        if not alts:
            raise ValueError(f"No entry in the FASTA has “{marker}” in its ID.")
    else:
        refs, alts = ref, alt
    by_token: dict[str, list[str]] = {}
    for rid in refs:
        for tok in _tokens(rid, ""):
            by_token.setdefault(tok, []).append(rid)
    by_ends: dict[str, list[str]] = {}
    for rid, (_, s) in refs.items():
        if len(s) >= 12:
            by_ends.setdefault("P" + s[:12], []).append(rid)
            by_ends.setdefault("S" + s[-12:], []).append(rid)

    pairs, unpaired, identical = [], [], 0
    for aid, (_, aseq) in alts.items():
        cands: list[str] = []
        sets = [by_token[t] for t in _tokens(aid, marker) if t in by_token]
        if sets:
            cands = min(sets, key=len)
        if not cands and len(aseq) >= 12:
            cands = by_ends.get("P" + aseq[:12], []) + by_ends.get("S" + aseq[-12:], [])
        same_len = [c for c in cands if len(refs[c][1]) == len(aseq)]
        if not same_len:
            unpaired.append(aid)
            continue
        best, best_d = None, None
        for c in same_len:
            rs = refs[c][1]
            d = [i for i, (x, y) in enumerate(zip(rs, aseq)) if x != y]
            if best_d is None or len(d) < len(best_d):
                best, best_d = c, d
        if not best_d:
            identical += 1
            continue
        pairs.append(Pair(best, aid, refs[best][0], refs[best][1], aseq, best_d))
    return Pairing(pairs, unpaired, identical, {k: v[1] for k, v in refs.items()},
                   {k: v[1] for k, v in alts.items()}, {k: v[0] for k, v in refs.items()})


# ---------------------------------------------------------------- peptide reports

@dataclass
class PeptideTable:
    samples: list[str]
    peptides: dict[str, dict]          # sequence -> {"charges", "intensity", "groups", "precursors", "rt", "q"}
    n_precursors: int
    prec_counts: list[int] = field(default_factory=list)       # precursors with signal, per run
    pep_counts: list[int] = field(default_factory=list)        # unique peptides with signal, per run
    total_intensity: list[float] = field(default_factory=list)
    proteotypic: float | None = None
    has_rt: bool = False
    has_q: bool = False


def _new_pep(n):
    return {"charges": set(), "intensity": np.zeros(n), "groups": set(), "precursors": {}, "rt": None, "q": None}


def _clean_run_name(c: str) -> str:
    return re.sub(r"\.(raw|mzml|d|wiff|dia)$", "", str(c).split("/")[-1].split("\\")[-1], flags=re.I)


def read_peptides(raw: bytes, filename: str = "") -> PeptideTable:
    """DIA-NN precursor matrix (wide, one column per run) or the long-format report (one row per precursor and run)."""
    text = raw.decode("utf-8-sig", errors="replace")
    delim = "\t" if text[:5000].count("\t") >= text[:5000].count(",") else ","
    df = pd.read_csv(io.StringIO(text), sep=delim, dtype=str, keep_default_na=False, low_memory=False)
    lower = {c.lower().strip(): c for c in df.columns}
    seq_col = next((lower[c] for c in SEQ_COLS if c in lower), None)
    if seq_col is None:
        raise ValueError("No peptide-sequence column found (expected DIA-NN's “Stripped.Sequence”). "
                         "Load the precursor matrix, report.pr_matrix.tsv. Columns seen: "
                         + ", ".join(list(df.columns)[:8]))
    mod_col, charge_col, group_col = lower.get("modified.sequence"), lower.get("precursor.charge") or lower.get("charge"), lower.get("protein.group")
    proto_col = lower.get("proteotypic")
    seqs = df[seq_col].str.upper().str.replace(r"[^A-Z]", "", regex=True).tolist()
    long_fmt = "run" in lower and any(k in lower for k in ("precursor.quantity", "precursor.normalised"))

    if long_fmt:
        q_col = lower.get("precursor.normalised") or lower.get("precursor.quantity")
        runs = sorted(df[lower["run"]].unique())
        ridx = {r: i for i, r in enumerate(runs)}
        peps: dict[str, dict] = {}
        vals = pd.to_numeric(df[q_col], errors="coerce").fillna(0.0).to_numpy(float)
        rts = pd.to_numeric(df[lower["rt"]], errors="coerce").to_numpy(float) if "rt" in lower else None
        qv = pd.to_numeric(df[lower["q.value"]], errors="coerce").to_numpy(float) if "q.value" in lower else None
        prec_by_run = np.zeros(len(runs))
        pep_by_run = [set() for _ in runs]
        total = np.zeros(len(runs))
        rt_lists: dict[str, list] = {}
        for i, sq in enumerate(seqs):
            if not sq:
                continue
            r = ridx[df[lower["run"]].iat[i]]
            p = peps.setdefault(sq, _new_pep(len(runs)))
            p["intensity"][r] += vals[i]
            ch = str(df[charge_col].iat[i]) if charge_col else ""
            if ch:
                p["charges"].add(ch)
            ms = df[mod_col].iat[i] if mod_col else sq
            pr = p["precursors"].setdefault((ms, ch), np.zeros(len(runs)))
            pr[r] += vals[i]
            if group_col:
                p["groups"].add(df[group_col].iat[i])
            if rts is not None and rts[i] == rts[i]:
                rt_lists.setdefault(sq, []).append(rts[i])
            if qv is not None and qv[i] == qv[i]:
                p["q"] = qv[i] if p["q"] is None else min(p["q"], qv[i])
            if vals[i] > 0:
                prec_by_run[r] += 1
                pep_by_run[r].add(sq)
                total[r] += vals[i]
        for sq, lst in rt_lists.items():
            peps[sq]["rt"] = float(np.median(lst))
        proteo = float(pd.to_numeric(df[proto_col], errors="coerce").mean()) if proto_col else None
        return PeptideTable([_clean_run_name(r) for r in runs], peps, len(df), prec_by_run.astype(int).tolist(),
                            [len(x) for x in pep_by_run], total.tolist(), proteo, rts is not None, qv is not None)

    sample_cols = [c for c in df.columns if c.lower().strip() not in META_COLS
                   and pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.2]
    if not sample_cols:
        raise ValueError("No intensity columns found next to the peptide sequences.")
    inten = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(float)
    peps = {}
    for i, sq in enumerate(seqs):
        if not sq:
            continue
        p = peps.setdefault(sq, _new_pep(len(sample_cols)))
        p["intensity"] += inten[i]
        ch = str(df[charge_col].iat[i]) if charge_col else ""
        if ch:
            p["charges"].add(ch)
        ms = df[mod_col].iat[i] if mod_col else sq
        pr = p["precursors"].setdefault((ms, ch), np.zeros(len(sample_cols)))
        pr += inten[i]
        if group_col:
            p["groups"].add(df[group_col].iat[i])
    keep = np.array([bool(x) for x in seqs])
    prec_counts = (inten[keep] > 0).sum(axis=0).astype(int).tolist()
    pep_counts = [sum(1 for p in peps.values() if p["intensity"][j] > 0) for j in range(len(sample_cols))]
    proteo = float(pd.to_numeric(df[proto_col], errors="coerce").mean()) if proto_col else None
    return PeptideTable([_clean_run_name(c) for c in sample_cols], peps, len(df), prec_counts, pep_counts,
                        inten[keep].sum(axis=0).tolist(), proteo)


# ---------------------------------------------------------------- site evidence

def tryptic_cuts(seq: str) -> list[int]:
    return [0] + [i + 1 for i in range(len(seq) - 1) if seq[i] in "KR" and seq[i + 1] != "P"] + [len(seq)]


def expected_peptides(seq: str, q: int, lo: int = 7, hi: int = 30) -> dict:
    """The tryptic peptides (up to one missed cleavage) that would span position q, and whether a search could see them."""
    cuts = tryptic_cuts(seq)
    k = max(i for i in range(len(cuts) - 1) if cuts[i] <= q)
    opts = []
    for a_ in (k - 1, k):
        for b_ in (k + 1, k + 2):
            if a_ >= 0 and b_ < len(cuts) and cuts[a_] <= q < cuts[b_]:
                opts.append((cuts[a_], cuts[b_]))
    opts = sorted(set(opts), key=lambda t: t[1] - t[0])
    cleaved = [(cuts[k], cuts[k + 1])]
    rows = [{"sequence": seq[a_:b_], "start": a_ + 1, "end": b_, "length": b_ - a_,
             "detectable": lo <= b_ - a_ <= hi, "missed_cleavages": _missed(seq[a_:b_])} for a_, b_ in opts]
    full = rows[0] if rows else None
    best = next((r for r in rows if r["detectable"]), None)
    first = {"sequence": seq[cleaved[0][0]:cleaved[0][1]], "length": cleaved[0][1] - cleaved[0][0]}
    if best:
        note = f"{best['sequence']} ({best['length']} residues) can be detected"
    elif full:
        note = (f"the peptide here would be {first['sequence']} ({first['length']} residues): "
                f"{'too short' if first['length'] < lo else 'too long'} to be detected"
                + (f"; with a missed cleavage {full['sequence']} ({full['length']})" if full['length'] != first['length'] else ""))
    else:
        note = "no peptide of a detectable length spans this position"
    return {"candidates": rows, "detectable": best is not None, "first": first, "note": note}


def _missed(pep: str) -> int:
    return sum(1 for i in range(len(pep) - 1) if pep[i] in "KR" and pep[i + 1] != "P")


def _spanning(seq: str, q: int, peps: dict) -> list[tuple[str, int]]:
    """Identified peptides that occupy seq[start:start+len] with start <= q < start+len."""
    out = []
    for start in range(max(0, q - MAX_PEP + 1), q + 1):
        for ln in range(max(MIN_PEP, q - start + 1), min(MAX_PEP, len(seq) - start) + 1):
            sub = seq[start:start + ln]
            if sub in peps:
                out.append((sub, start))
    return out


def _count(blob: str, pep: str) -> int:
    n, i = 0, blob.find(pep)
    while i != -1:
        n += 1
        i = blob.find(pep, i + 1)
    return n


LENGTH_BINS = [(0, 150, "under 150"), (150, 300, "150-300"), (300, 600, "300-600"), (600, 10 ** 9, "over 600")]


def coverage(pairing: Pairing, table: PeptideTable) -> dict:
    """How much of the reference proteome the identified peptides cover."""
    import bisect
    ids = list(pairing.ref_pool)
    seqs = [pairing.ref_pool[i] for i in ids]
    offsets, blob, pos = [], [], 0
    for sq in seqs:
        offsets.append(pos)
        blob.append(sq)
        pos += len(sq) + 1
    blob = "\n".join(blob)
    cov = [bytearray(len(sq)) for sq in seqs]
    peps_of = [set() for _ in ids]
    unmapped = shared = 0
    for pep in table.peptides:
        hits, i = set(), blob.find(pep)
        while i != -1:
            j = bisect.bisect_right(offsets, i) - 1
            start = i - offsets[j]
            if start + len(pep) <= len(seqs[j]):
                hits.add(j)
                cov[j][start:start + len(pep)] = b"\x01" * len(pep)
                peps_of[j].add(pep)
            i = blob.find(pep, i + 1)
        if not hits:
            unmapped += 1
        elif len(hits) > 1:
            shared += 1
    variant_ids = {p.ref_id for p in pairing.pairs}
    rows = []
    for j, pid in enumerate(ids):
        n = len(seqs[j])
        rows.append({"protein": pid, "description": pairing.ref_headers.get(pid, "")[len(pid):].strip()[:140],
                     "length": n, "n_peptides": len(peps_of[j]),
                     "coverage_pct": round(100 * sum(cov[j]) / n, 1) if n else 0.0,
                     "has_variant": pid in variant_ids})
    det = [r for r in rows if r["n_peptides"] > 0]
    pct = [r["coverage_pct"] for r in det]
    total_res = sum(r["length"] for r in rows)
    covered = sum(sum(c) for c in cov)
    by_len = []
    for lo, hi, label in LENGTH_BINS:
        grp = [r for r in rows if lo <= r["length"] < hi]
        by_len.append({"label": label, "n": len(grp), "detected": sum(r["n_peptides"] > 0 for r in grp)})
    lengths = {}
    for pep in table.peptides:
        lengths[len(pep)] = lengths.get(len(pep), 0) + 1
    charges: dict[str, int] = {}
    for v in table.peptides.values():
        for c in v["charges"]:
            charges[c] = charges.get(c, 0) + 1
    missed: dict[int, int] = {}
    for pep in table.peptides:
        m = _missed(pep)
        missed[m] = missed.get(m, 0) + 1
    ints = sorted((float(v["intensity"].sum()) for v in table.peptides.values() if v["intensity"].sum() > 0))
    var_rows = [r for r in rows if r["has_variant"]]
    return {
        "n_proteins": len(rows), "n_residues": total_res, "n_detected": len(det),
        "n_ge2": sum(r["n_peptides"] >= 2 for r in rows), "n_ge3": sum(r["n_peptides"] >= 3 for r in rows),
        "n_ge5": sum(r["n_peptides"] >= 5 for r in rows),
        "median_cov": float(np.median(pct)) if pct else 0.0, "mean_cov": float(np.mean(pct)) if pct else 0.0,
        "n_cov50": sum(p >= 50 for p in pct), "overall_cov": 100 * covered / total_res if total_res else 0.0,
        "hist": np.histogram(pct, bins=10, range=(0, 100))[0].tolist() if pct else [0] * 10,
        "by_length": by_len, "pep_lengths": {str(k): v for k, v in sorted(lengths.items())},
        "charges": dict(sorted(charges.items())), "missed": {str(k): v for k, v in sorted(missed.items())},
        "n_unmapped": unmapped, "n_shared": shared, "n_peptides": len(table.peptides),
        "intensity_orders": float(np.log10(ints[-1] / ints[0])) if len(ints) > 1 else 0.0,
        "variant_proteins": len(var_rows), "variant_proteins_detected": sum(r["n_peptides"] > 0 for r in var_rows),
        "proteins": rows,
    }


@dataclass
class Result:
    samples: list[str]
    sites: list[dict]
    summary: dict
    scatter: dict
    coverage: dict = field(default_factory=dict)


def analyse(pairing: Pairing, table: PeptideTable, table_b: PeptideTable | None = None) -> Result:
    ref_blob = "\n".join(pairing.ref_pool.values())
    alt_blob = "\n".join(pairing.alt_pool.values())
    pep_a, pep_b = table.peptides, (table_b.peptides if table_b else {})
    sites: list[dict] = []
    for pr in pairing.pairs:
        for q in pr.positions:
            row = {"ref_protein": pr.ref_id, "alt_protein": pr.alt_id, "position": q + 1,
                   "ref_aa": pr.ref_seq[q], "alt_aa": pr.alt_seq[q],
                   "description": pr.ref_header[len(pr.ref_id):].strip()[:160],
                   "context_ref": pr.ref_seq[max(0, q - 7):q] + "[" + pr.ref_seq[q] + "]" + pr.ref_seq[q + 1:q + 8],
                   "context_alt": pr.alt_seq[max(0, q - 7):q] + "[" + pr.alt_seq[q] + "]" + pr.alt_seq[q + 1:q + 8],
                   "ref_peptides": [], "alt_peptides": [], "ref_peptides_b": [],
                   "expected_ref": expected_peptides(pr.ref_seq, q), "expected_alt": expected_peptides(pr.alt_seq, q)}
            own_ref = pr.ref_seq
            for version, seq in (("ref", pr.ref_seq), ("alt", pr.alt_seq)):
                for pep, start in _spanning(seq, q, pep_a):
                    if version == "ref":          # the same peptide in another reference protein
                        other = _count(ref_blob, pep) - _count(own_ref, pep)
                    else:                         # a variant peptide that also exists elsewhere
                        other = _count(ref_blob, pep) + _count(alt_blob, pep) - _count(pr.alt_seq, pep)
                    rec = {"sequence": pep, "start": start + 1, "end": start + len(pep),
                           "charges": sorted(pep_a[pep]["charges"]), "missed_cleavages": _missed(pep),
                           "intensity": pep_a[pep]["intensity"].tolist(), "rt": pep_a[pep]["rt"], "q": pep_a[pep]["q"],
                           "precursors": [[m, z] for (m, z) in pep_a[pep]["precursors"]],
                           "total_intensity": float(pep_a[pep]["intensity"].sum()),
                           "also_elsewhere": bool(other > 0),
                           "seen_in_b": pep in pep_b}
                    row[f"{version}_peptides"].append(rec)
            if pep_b:
                for pep, start in _spanning(pr.ref_seq, q, pep_b):
                    row["ref_peptides_b"].append({"sequence": pep, "start": start + 1, "end": start + len(pep),
                                                  "total_intensity": float(pep_b[pep]["intensity"].sum()),
                                                  "seen_in_a": pep in pep_a})
            for key in ("ref_peptides", "alt_peptides"):
                row[key].sort(key=lambda r: -r["total_intensity"])
            n_ref = sum(not p["also_elsewhere"] for p in row["ref_peptides"])
            n_alt = sum(not p["also_elsewhere"] for p in row["alt_peptides"])
            row["n_ref"], row["n_alt"] = n_ref, n_alt
            row["status"] = ("both" if n_ref and n_alt else "alternative" if n_alt else
                             "reference" if n_ref else "none")
            row["ambiguous"] = (len(row["ref_peptides"]) - n_ref) + (len(row["alt_peptides"]) - n_alt)
            sites.append(row)
    counts = {k: sum(1 for s in sites if s["status"] == k) for k in ("alternative", "reference", "both", "none")}
    summary = {"n_sites": len(sites), "n_proteins": len(pairing.pairs), "status": counts,
               "unpaired": pairing.unpaired[:50], "n_unpaired": len(pairing.unpaired),
               "identical_alt": pairing.identical, "n_peptides": len(pep_a), "n_precursors": table.n_precursors,
               "has_b": table_b is not None, "n_peptides_b": len(pep_b)}
    scatter = {}
    if table_b is not None:
        shared = [p for p in pep_a if p in pep_b]
        a = np.array([pep_a[p]["intensity"].sum() for p in shared])
        b = np.array([pep_b[p]["intensity"].sum() for p in shared])
        ok = (a > 0) & (b > 0)
        r = float(np.corrcoef(np.log2(a[ok]), np.log2(b[ok]))[0, 1]) if ok.sum() > 2 else None
        idx = np.linspace(0, ok.sum() - 1, min(ok.sum(), 4000)).astype(int) if ok.sum() else []
        scatter = {"a": np.log2(a[ok])[idx].round(3).tolist() if len(idx) else [],
                   "b": np.log2(b[ok])[idx].round(3).tolist() if len(idx) else [],
                   "n_shared": len(shared), "n_only_a": len(pep_a) - len(shared),
                   "n_only_b": len(pep_b) - len(shared), "r": r}
    return Result(table.samples, sites, summary, scatter, coverage(pairing, table))


def sites_frame(res: Result) -> pd.DataFrame:
    rows = []
    for s in res.sites:
        rows.append({
            "reference_protein": s["ref_protein"], "alternative_protein": s["alt_protein"], "position": s["position"],
            "reference_aa": s["ref_aa"], "alternative_aa": s["alt_aa"], "status": s["status"],
            "reference_peptides": ";".join(p["sequence"] for p in s["ref_peptides"] if not p["also_elsewhere"]),
            "alternative_peptides": ";".join(p["sequence"] for p in s["alt_peptides"] if not p["also_elsewhere"]),
            "ambiguous_peptides": ";".join(p["sequence"] for k in ("ref_peptides", "alt_peptides")
                                           for p in s[k] if p["also_elsewhere"]),
            "reference_intensity": sum(p["total_intensity"] for p in s["ref_peptides"] if not p["also_elsewhere"]),
            "alternative_intensity": sum(p["total_intensity"] for p in s["alt_peptides"] if not p["also_elsewhere"]),
            "reference_peptides_in_standard_search": ";".join(p["sequence"] for p in s["ref_peptides_b"]),
            "reference_version_detectable": s["expected_ref"]["detectable"], "alternative_version_detectable": s["expected_alt"]["detectable"],
            "reference_detectability_note": s["expected_ref"]["note"], "alternative_detectability_note": s["expected_alt"]["note"],
            "context_reference": s["context_ref"], "context_alternative": s["context_alt"],
            "description": s["description"]})
    return pd.DataFrame(rows)


def peptides_frame(res: Result) -> pd.DataFrame:
    rows = []
    for s in res.sites:
        for version in ("ref", "alt"):
            for p in s[f"{version}_peptides"]:
                row = {"site_protein": s["ref_protein"], "site_position": s["position"],
                       "version": "reference" if version == "ref" else "alternative",
                       "residue_at_site": s["ref_aa"] if version == "ref" else s["alt_aa"],
                       "peptide": p["sequence"], "start": p["start"], "end": p["end"],
                       "charges": ";".join(p["charges"]), "missed_cleavages": p["missed_cleavages"],
                       "also_elsewhere_in_database": p["also_elsewhere"], "seen_in_standard_search": p["seen_in_b"]}
                for name, v in zip(res.samples, p["intensity"]):
                    row[f"intensity_{name}"] = v
                rows.append(row)
    return pd.DataFrame(rows)


def coverage_frame(res: Result) -> pd.DataFrame:
    return pd.DataFrame(res.coverage["proteins"]).rename(columns={"has_variant": "has_alternative_version"})

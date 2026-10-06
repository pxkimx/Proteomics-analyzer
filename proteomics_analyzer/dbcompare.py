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
                   {k: v[1] for k, v in alts.items()})


# ---------------------------------------------------------------- peptide reports

@dataclass
class PeptideTable:
    samples: list[str]
    peptides: dict[str, dict]          # sequence -> {"charges": set, "intensity": np.ndarray, "groups": set}
    n_precursors: int


def read_peptides(raw: bytes, filename: str = "") -> PeptideTable:
    text = raw.decode("utf-8-sig", errors="replace")
    delim = "\t" if text[:5000].count("\t") >= text[:5000].count(",") else ","
    df = pd.read_csv(io.StringIO(text), sep=delim, dtype=str, keep_default_na=False, low_memory=False)
    lower = {c.lower().strip(): c for c in df.columns}
    seq_col = next((lower[c] for c in SEQ_COLS if c in lower), None)
    if seq_col is None:
        raise ValueError("No peptide-sequence column found (expected DIA-NN's “Stripped.Sequence”). "
                         "Load the precursor matrix, report.pr_matrix.tsv. Columns seen: "
                         + ", ".join(list(df.columns)[:8]))
    sample_cols = [c for c in df.columns if c.lower().strip() not in META_COLS
                   and pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.2]
    if not sample_cols:
        raise ValueError("No intensity columns found next to the peptide sequences.")
    inten = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(float)
    charge_col = lower.get("precursor.charge") or lower.get("charge")
    group_col = lower.get("protein.group")
    peps: dict[str, dict] = {}
    seqs = df[seq_col].str.upper().str.replace(r"[^A-Z]", "", regex=True).tolist()
    for i, s in enumerate(seqs):
        if not s:
            continue
        p = peps.setdefault(s, {"charges": set(), "intensity": np.zeros(len(sample_cols)), "groups": set()})
        p["intensity"] += inten[i]
        if charge_col:
            p["charges"].add(str(df[charge_col].iat[i]))
        if group_col:
            p["groups"].add(df[group_col].iat[i])
    names = [re.sub(r"\.(raw|mzml|d|wiff|dia)$", "", c.split("/")[-1].split("\\")[-1], flags=re.I) for c in sample_cols]
    return PeptideTable(names, peps, len(df))


# ---------------------------------------------------------------- site evidence

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


@dataclass
class Result:
    samples: list[str]
    sites: list[dict]
    summary: dict
    scatter: dict


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
                   "ref_peptides": [], "alt_peptides": [], "ref_peptides_b": []}
            own_ref = pr.ref_seq
            for version, seq in (("ref", pr.ref_seq), ("alt", pr.alt_seq)):
                for pep, start in _spanning(seq, q, pep_a):
                    if version == "ref":          # the same peptide in another reference protein
                        other = _count(ref_blob, pep) - _count(own_ref, pep)
                    else:                         # a variant peptide that also exists elsewhere
                        other = _count(ref_blob, pep) + _count(alt_blob, pep) - _count(pr.alt_seq, pep)
                    rec = {"sequence": pep, "start": start + 1, "end": start + len(pep),
                           "charges": sorted(pep_a[pep]["charges"]), "missed_cleavages": _missed(pep),
                           "intensity": pep_a[pep]["intensity"].tolist(),
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
    return Result(table.samples, sites, summary, scatter)


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

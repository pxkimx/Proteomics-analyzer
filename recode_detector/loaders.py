"""Read protein-level quantification tables (MaxQuant, DIA-NN, generic) into one shape."""
from __future__ import annotations

import io
import os
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

DIANN_META = {
    "protein.group", "protein.ids", "protein.names", "genes", "first.protein.description",
    "n.sequences", "n.proteotypic.sequences", "proteotypic", "precursor.charge", "precursor.id",
    "stripped.sequence", "modified.sequence", "protein.q.value", "pg.q.value",
}
DIANN_PRECURSOR_MARKERS = {"precursor.id", "modified.sequence", "stripped.sequence", "precursor.charge"}
GENE_NAMES = ("gene names", "gene name", "genes", "gene", "gene symbol", "symbol", "pg.genes")
ID_NAMES = ("majority protein ids", "protein ids", "protein.group", "protein group", "protein",
            "accession", "uniprot", "pg.proteingroups", "protein id", "id")


@dataclass
class Dataset:
    """Protein x sample matrix of raw (linear) intensities; missing values are NaN."""
    source_format: str
    matrix: pd.DataFrame                       # index = protein id, columns = sample names
    genes: pd.Series                           # protein id -> gene label
    notes: list[str] = field(default_factory=list)
    removed: dict[str, int] = field(default_factory=dict)

    @property
    def samples(self) -> list[str]:
        return list(self.matrix.columns)


def _read_text(raw: bytes) -> pd.DataFrame:
    text = raw.decode("utf-8-sig", errors="replace")
    head = text[:20000]
    delim = "\t" if head.count("\t") >= head.count(",") else ","
    return pd.read_csv(io.StringIO(text), sep=delim, low_memory=False, dtype=str,
                       keep_default_na=False)


def read_any(raw: bytes, filename: str) -> pd.DataFrame:
    ext = os.path.splitext(filename.lower())[1]
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(io.BytesIO(raw), dtype=str).fillna("")
    return _read_text(raw)


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.replace(",", "", regex=False).str.strip(),
                         errors="coerce")


def clean_sample_name(name: str) -> str:
    n = re.sub(r"^(LFQ intensity|Intensity|iBAQ|Reporter intensity corrected \d+|Reporter intensity)\s+", "",
               str(name).strip(), flags=re.I)
    n = n.replace("\\", "/").split("/")[-1]
    n = re.sub(r"\.(raw|mzml|d|wiff|dia|parquet|tsv|txt)$", "", n, flags=re.I)
    return n.strip() or str(name)


def _first_token(value: str) -> str:
    return str(value).split(";")[0].strip()


def _find(df: pd.DataFrame, wanted: tuple[str, ...]) -> str | None:
    lower = {c.lower().strip(): c for c in df.columns}
    for w in wanted:
        if w in lower:
            return lower[w]
    return None


def detect_format(df: pd.DataFrame) -> str:
    cols = [c.lower() for c in df.columns]
    if "protein ids" in cols and any(c.startswith(("lfq intensity ", "intensity ")) for c in cols):
        return "maxquant"
    if "protein.group" in cols:
        return "diann"
    return "generic"


def maxquant_intensity_kinds(df: pd.DataFrame) -> list[str]:
    kinds = []
    for kind in ("LFQ intensity", "Intensity", "iBAQ"):
        if any(c.startswith(kind + " ") and not (kind == "Intensity" and c.startswith("Intensity L"))
               for c in df.columns):
            kinds.append(kind)
    return kinds


def load(raw: bytes, filename: str, intensity_kind: str | None = None) -> Dataset:
    df = read_any(raw, filename)
    if df.shape[0] == 0 or df.shape[1] < 3:
        raise ValueError("The table has no data rows or fewer than three columns.")
    fmt = detect_format(df)
    notes: list[str] = []
    removed: dict[str, int] = {}

    if fmt == "maxquant":
        kinds = maxquant_intensity_kinds(df)
        kind = intensity_kind if intensity_kind in kinds else ("LFQ intensity" if "LFQ intensity" in kinds
                                                                 else kinds[0])
        sample_cols = [c for c in df.columns if c.startswith(kind + " ") and c != kind]
        if kind == "Intensity":
            sample_cols = [c for c in sample_cols if not re.match(r"Intensity [HLM]\b", c)]
        notes.append(f"MaxQuant proteinGroups table; using “{kind}” columns "
                     f"({len(sample_cols)} samples).")
        keep = pd.Series(True, index=df.index)
        for flag, label in (("Reverse", "reverse hits"), ("Potential contaminant", "contaminants"),
                            ("Only identified by site", "only identified by site")):
            if flag in df.columns:
                hit = df[flag].astype(str).str.strip() == "+"
                removed[label] = int((hit & keep).sum())
                keep &= ~hit
        df = df[keep]
        id_col = "Majority protein IDs" if "Majority protein IDs" in df.columns else "Protein IDs"
        gene_col = "Gene names" if "Gene names" in df.columns else None
        ids = df[id_col].map(_first_token)
        genes = df[gene_col].map(_first_token) if gene_col else ids
    elif fmt == "diann":
        if DIANN_PRECURSOR_MARKERS & {c.lower() for c in df.columns}:
            raise ValueError(
                "This looks like DIA-NN's precursor matrix (pr_matrix: one row per peptide ion). This program "
                "works on protein-level tables here: load the protein group matrix (report.pg_matrix.tsv) instead. "
                "To look at peptides, for example which version of a protein was seen against two databases, "
                "use the “Database check” mode (switch at the top of the page).")
        sample_cols = [c for c in df.columns if c.lower() not in DIANN_META
                       and _num(df[c]).notna().mean() > 0.2]
        notes.append(f"DIA-NN protein group matrix ({len(sample_cols)} runs).")
        ids = df[_find(df, ("protein.group",))].map(_first_token)
        g = _find(df, ("genes",))
        genes = df[g].map(_first_token) if g else ids
    else:
        id_col = _find(df, ID_NAMES) or df.columns[0]
        gene_col = _find(df, GENE_NAMES)
        used = {id_col, gene_col}
        sample_cols = [c for c in df.columns if c not in used and _num(df[c]).notna().mean() > 0.5]
        if len(sample_cols) < 2:
            raise ValueError("Could not find at least two numeric intensity columns.")
        notes.append(f"Generic table; ID column “{id_col}”, "
                     f"{'gene column “' + gene_col + '”' if gene_col else 'no gene column'}, "
                     f"{len(sample_cols)} numeric sample columns.")
        ids = df[id_col].map(_first_token)
        genes = df[gene_col].map(_first_token) if gene_col else ids

    if len(sample_cols) < 2:
        raise ValueError(
            f"Only {len(sample_cols)} sample column found"
            f"{' (' + clean_sample_name(sample_cols[0]) + ')' if sample_cols else ''}. Comparing groups needs "
            "replicate runs: process several runs together and load the combined protein table.")

    mat = pd.DataFrame({c: _num(df[c]) for c in sample_cols})
    mat = mat.where(mat > 0)                       # 0 means "not quantified"
    names = [clean_sample_name(c) for c in sample_cols]
    if len(set(names)) != len(names):
        names = [str(c) for c in sample_cols]
    mat.columns = names
    ids = ids.reset_index(drop=True)
    genes = genes.reset_index(drop=True)
    mat = mat.reset_index(drop=True)

    blank = (ids == "")
    if blank.any():
        removed["rows without an ID"] = int(blank.sum())
    ok = ~blank & mat.notna().any(axis=1)
    removed["rows with no quantified value"] = int((~blank & ~mat.notna().any(axis=1)).sum())
    mat, ids, genes = mat[ok], ids[ok], genes[ok]

    # make IDs unique (isoform groups can repeat a leading accession)
    ids = ids.where(~ids.duplicated(keep=False), ids + "_" + (ids.groupby(ids).cumcount() + 1).astype(str))
    genes = genes.where(genes != "", ids)
    mat.index = ids.values
    gene_series = pd.Series(genes.values, index=ids.values)
    return Dataset(fmt, mat, gene_series, notes, {k: v for k, v in removed.items() if v})


def guess_groups(samples: list[str]) -> dict[str, str]:
    """Group = sample name without its trailing replicate number (Ctrl_1, ctrl-2, KO rep3 ...)."""
    out = {}
    for s in samples:
        g = re.sub(r"[\s_\-\.]*(rep(licate)?|r|n|bio|tech)?[\s_\-\.]*\d+$", "", s, flags=re.I).strip(" _-.")
        out[s] = g or s
    counts = pd.Series(out).value_counts()
    if len(counts) == len(samples):                # nothing repeats: the guess is useless
        return {s: "" for s in samples}
    return out

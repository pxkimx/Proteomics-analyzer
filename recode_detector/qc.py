"""Quality-control checks. Thresholds are rules of thumb, shown next to each result so they can be judged."""
from __future__ import annotations

import numpy as np
import pandas as pd

PASS, WARN, FAIL, INFO = "pass", "warn", "fail", "info"


def _check(cid, name, value, status, rule, why):
    return {"id": cid, "name": name, "value": value, "status": status, "rule": rule, "why": why}


def _grade(v, ok, warn, higher_is_better=True):
    if higher_is_better:
        return PASS if v >= ok else WARN if v >= warn else FAIL
    return PASS if v <= ok else WARN if v <= warn else FAIL


# ---------------------------------------------------------------- database-check mode

def db_checks(table, cov: dict) -> tuple[list[dict], pd.DataFrame]:
    peps = table.peptides
    n = max(len(peps), 1)
    checks = []
    mc = cov["missed"]
    mc0 = mc.get("0", 0) / n * 100
    checks.append(_check("missed", "Peptides with no missed cleavage", f"{mc0:.0f}%", _grade(mc0, 80, 65),
                         "pass at 80% or more", "A high share of missed cuts suggests incomplete digestion."))
    ch = cov["charges"]
    tot = max(sum(ch.values()), 1)
    c23 = (ch.get("2", 0) + ch.get("3", 0)) / tot * 100
    checks.append(_check("charge", "Precursors with charge 2+ or 3+", f"{c23:.0f}%", _grade(c23, 85, 70),
                         "pass at 85% or more", "Tryptic peptides are mostly 2+ and 3+; a very different mix points to a problem."))
    lens = {int(k): v for k, v in cov["pep_lengths"].items()}
    med = _median_from_counts(lens)
    checks.append(_check("length", "Median peptide length", f"{med} amino acids", PASS if 8 <= med <= 16 else WARN,
                         "pass between 8 and 16", "Typical tryptic peptides are 8 to 16 residues long."))
    det = cov["n_detected"]
    single = det - cov["n_ge2"]
    sh = single / max(det, 1) * 100
    checks.append(_check("single", "Proteins identified from one peptide only", f"{sh:.0f}%", _grade(sh, 25, 40, False),
                         "pass at 25% or less", "One-peptide identifications are the least certain."))
    det_pct = 100 * det / max(cov["n_proteins"], 1)
    checks.append(_check("detected", "Share of the database detected", f"{det_pct:.0f}%", INFO,
                         "informative only", "Depends on the organism, sample and depth; compare runs of the same kind."))
    if table.proteotypic is not None:
        checks.append(_check("proteotypic", "Precursors unique to one protein", f"{table.proteotypic * 100:.0f}%",
                             _grade(table.proteotypic * 100, 90, 75), "pass at 90% or more",
                             "Peptides shared between proteins cannot tell them apart."))
    orders = cov["intensity_orders"]
    checks.append(_check("range", "Dynamic range of peptide signal", f"{orders:.1f} orders of magnitude",
                         _grade(orders, 3, 2), "pass at 3 or more", "A narrow range means weak or saturated data."))
    ns = len(table.samples)
    runs = pd.DataFrame({"run": table.samples, "precursors": table.prec_counts,
                         "peptides": table.pep_counts, "total_intensity": table.total_intensity})
    if ns >= 3 and runs["precursors"].mean() > 0:
        cv = runs["precursors"].std(ddof=1) / runs["precursors"].mean() * 100
        checks.append(_check("runs", "Spread of identifications across runs (CV)", f"{cv:.0f}%", _grade(cv, 15, 30, False),
                             "pass at 15% or less", "Runs of the same sample should identify similar numbers of precursors."))
    elif ns == 1:
        checks.append(_check("runs", "Replicate runs", "1 run", INFO, "informative only",
                             "With one run, repeatability cannot be judged."))
    return checks, runs


def _median_from_counts(counts: dict[int, int]) -> int:
    total, run = sum(counts.values()), 0
    for k in sorted(counts):
        run += counts[k]
        if run >= total / 2:
            return k
    return 0


# ---------------------------------------------------------------- quantification mode

def quant_checks(raw_matrix: pd.DataFrame, proc, qc: dict) -> tuple[list[dict], pd.DataFrame]:
    samples = list(proc.log_raw.columns)
    lin = raw_matrix[samples]
    ident = lin.notna().sum()
    spread = (ident.max() - ident.min()) / max(ident.median(), 1) * 100
    miss = proc.log_raw.isna().mean() * 100
    cvs = [v["median"] for v in qc["cv"].values()]
    # correlation of each sample to the median profile of its own group (leave-one-out)
    corr_to_group = {}
    groups = proc.groups
    for s in samples:
        mates = [t for t in samples if groups[t] == groups[s] and t != s]
        if mates:
            ref = proc.log_imp[mates].median(axis=1)
            corr_to_group[s] = float(np.corrcoef(proc.log_imp[s], ref)[0, 1])
        else:
            corr_to_group[s] = float("nan")
    cs = pd.Series(corr_to_group)
    flagged = []
    ok = cs.dropna()
    if len(ok) >= 4:
        med, mad = ok.median(), max((ok - ok.median()).abs().median() * 1.4826, 1e-3)
        flagged = [s for s, v in ok.items() if (v - med) / mad < -3 or v < 0.85]
    elif len(ok):
        flagged = [s for s, v in ok.items() if v < 0.85]
    rows = pd.DataFrame({
        "sample": samples, "group": [groups[s] for s in samples], "proteins_identified": [int(ident[s]) for s in samples],
        "missing_pct": [round(float(miss[s]), 1) for s in samples], "median_log2": [round(float(proc.log_raw[s].median()), 2) for s in samples],
        "correlation_to_group": [round(corr_to_group[s], 3) if corr_to_group[s] == corr_to_group[s] else None for s in samples],
        "flag": ["possible outlier" if s in flagged else "" for s in samples]})
    checks = [
        _check("ident", "Proteins identified per sample: spread", f"{spread:.0f}% of the median", _grade(spread, 20, 35, False),
               "pass at 20% or less", "Samples of one experiment should have similar depth."),
        _check("missing", "Worst missing-value share in a sample", f"{miss.max():.0f}%", _grade(miss.max(), 30, 50, False),
               "pass at 30% or less", "Many missing values weaken statistics and depend on imputation."),
        _check("cv", "Median CV within groups", f"{np.mean(cvs):.0f}%" if cvs else "n/a",
               _grade(np.mean(cvs), 20, 30, False) if cvs else INFO, "pass at 20% or less", "Replicates of one condition should agree."),
        _check("corr", "Lowest replicate correlation", f"{np.nanmin(cs):.3f}" if cs.notna().any() else "n/a",
               _grade(np.nanmin(cs), .9, .8) if cs.notna().any() else INFO, "pass at 0.90 or more",
               "A sample far from its replicates may be a failed run."),
        _check("outliers", "Possible outlier samples", ", ".join(flagged) if flagged else "none", WARN if flagged else PASS,
               "none flagged", "Flagged when correlation to the group is unusually low (robust z below -3) or under 0.85."),
    ]
    sh = qc["pca"]["explained"][0] if qc["pca"]["explained"] else 0
    checks.append(_check("pc1", "Variance on the first principal component", f"{sh:.0f}%", INFO, "informative only",
                         "A strong first component usually reflects the main difference between groups (or a batch effect)."))
    return checks, rows


def summarise(checks: list[dict]) -> dict:
    c = {k: sum(1 for x in checks if x["status"] == k) for k in (PASS, WARN, FAIL, INFO)}
    c["overall"] = FAIL if c[FAIL] else WARN if c[WARN] else PASS
    return c

"""Filtering, normalisation, imputation and quality-control summaries."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy

NORMALIZATIONS = ("none", "median", "quantile")
IMPUTATIONS = ("none", "downshift", "halfmin")


@dataclass
class Params:
    normalization: str = "median"
    imputation: str = "downshift"
    min_valid: int = 2                 # valid values needed per group
    valid_rule: str = "any"            # "any" group meets min_valid, or "all" groups do
    downshift: float = 1.8             # Perseus convention, in SD of each sample
    width: float = 0.3
    seed: int = 1

    @classmethod
    def from_dict(cls, d: dict) -> "Params":
        p = cls()
        for k, v in (d or {}).items():
            if hasattr(p, k):
                setattr(p, k, type(getattr(p, k))(v))
        if p.normalization not in NORMALIZATIONS:
            raise ValueError(f"Unknown normalization: {p.normalization}")
        if p.imputation not in IMPUTATIONS:
            raise ValueError(f"Unknown imputation: {p.imputation}")
        if p.valid_rule not in ("any", "all"):
            raise ValueError("valid_rule must be 'any' or 'all'")
        return p


@dataclass
class Processed:
    log_raw: pd.DataFrame              # log2 intensities, NaN = missing, kept proteins only, grouped samples only
    log_norm: pd.DataFrame             # after normalisation, NaN kept
    log_imp: pd.DataFrame              # after imputation (== log_norm when imputation is "none")
    imputed: pd.DataFrame              # bool mask of imputed cells
    groups: dict[str, str]             # sample -> group
    filtered_out: int
    params: Params
    log: list[str] = field(default_factory=list)


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalisation that tolerates missing values (ranks are taken over valid values only)."""
    grid = np.linspace(0, 1, 201)
    curves = []
    for c in df.columns:
        v = np.sort(df[c].dropna().values)
        curves.append(np.interp(grid, np.linspace(0, 1, len(v)), v) if len(v) else np.full(len(grid), np.nan))
    ref = np.nanmean(np.vstack(curves), axis=0)
    out = df.copy()
    for c in df.columns:
        v = df[c]
        valid = v.notna()
        if valid.sum() < 2:
            continue
        ranks = v[valid].rank(method="average").values
        pct = (ranks - 1) / max(valid.sum() - 1, 1)
        out.loc[valid, c] = np.interp(pct, grid, ref)
    return out


def process(matrix: pd.DataFrame, groups: dict[str, str], params: Params) -> Processed:
    keep = [s for s in matrix.columns if groups.get(s, "").strip()]
    if len(keep) < 2:
        raise ValueError("Assign at least two samples to groups first.")
    gmap = {s: groups[s].strip() for s in keep}
    levels = sorted(set(gmap.values()))
    if len(levels) < 2:
        raise ValueError("At least two different groups are needed.")
    log: list[str] = []

    raw = np.log2(matrix[keep])
    n0 = len(raw)

    counts = pd.DataFrame({g: raw[[s for s in keep if gmap[s] == g]].notna().sum(axis=1) for g in levels})
    ok = (counts >= params.min_valid).any(axis=1) if params.valid_rule == "any" \
        else (counts >= params.min_valid).all(axis=1)
    raw = raw[ok]
    log.append(f"Valid-value filter ({params.min_valid} per group, {params.valid_rule} group"
               f"{'s' if params.valid_rule == 'all' else ''}): kept {len(raw)} of {n0} proteins.")
    if len(raw) < 10:
        raise ValueError(f"Only {len(raw)} proteins pass the valid-value filter; relax it.")

    if params.normalization == "median":
        shift = raw.median(axis=0) - raw.median(axis=0).mean()
        norm = raw - shift
        log.append("Median normalisation: each sample's median set to the mean of the sample medians.")
    elif params.normalization == "quantile":
        norm = quantile_normalize(raw)
        log.append("Quantile normalisation (missing values excluded from the ranks).")
    else:
        norm = raw.copy()
        log.append("No normalisation.")

    mask = norm.isna()
    imp = norm.copy()
    if params.imputation == "downshift":
        rng = np.random.default_rng(params.seed)
        for c in imp.columns:
            v = norm[c].dropna()
            m = mask[c]
            if m.any():
                imp.loc[m, c] = rng.normal(v.mean() - params.downshift * v.std(),
                                           max(params.width * v.std(), 1e-6), int(m.sum()))
        log.append(f"Imputed {int(mask.values.sum())} missing values from a normal distribution "
                   f"shifted {params.downshift} SD down, width {params.width} SD (seed {params.seed}).")
    elif params.imputation == "halfmin":
        for c in imp.columns:
            m = mask[c]
            if m.any():
                imp.loc[m, c] = norm[c].min() - 1.0          # log2: half of the minimum
        log.append(f"Imputed {int(mask.values.sum())} missing values with half the sample minimum.")
    else:
        log.append("No imputation: statistics use the observed values only.")
        mask = pd.DataFrame(False, index=norm.index, columns=norm.columns)
    return Processed(raw, norm, imp, mask, gmap, n0 - len(raw), params, log)


# ---------------------------------------------------------------- QC summaries

def _box(values: np.ndarray) -> dict:
    v = values[~np.isnan(values)]
    if len(v) == 0:
        return {}
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    iqr = q3 - q1
    return {"min": float(max(v.min(), q1 - 1.5 * iqr)), "q1": float(q1), "med": float(med),
            "q3": float(q3), "max": float(min(v.max(), q3 + 1.5 * iqr))}


def pca(imp: pd.DataFrame, n: int = 4) -> dict:
    X = imp.values.T                       # samples x proteins
    X = X - X.mean(axis=0)
    n = min(n, X.shape[0] - 1)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    var = S ** 2
    total = var.sum() or 1.0
    scores = (U * S)[:, :n]
    return {"scores": scores.tolist(), "explained": (var[:n] / total * 100).tolist()}


def qc_summary(raw_matrix: pd.DataFrame, proc: Processed) -> dict:
    samples = list(proc.log_raw.columns)
    lin = raw_matrix[samples]
    per_sample = []
    for s in samples:
        per_sample.append({
            "sample": s, "group": proc.groups[s],
            "identified": int(lin[s].notna().sum()),
            "after_filter": int(proc.log_raw[s].notna().sum()),
            "missing_pct": float(proc.log_raw[s].isna().mean() * 100),
            "box_raw": _box(proc.log_raw[s].values), "box_norm": _box(proc.log_norm[s].values),
        })
    # CV on the linear scale, proteins with >=3 valid values in the group
    cv = {}
    for g in sorted(set(proc.groups.values())):
        cols = [s for s in samples if proc.groups[s] == g]
        sub = 2 ** proc.log_norm[cols]
        ok = sub.notna().sum(axis=1) >= 3
        if ok.any():
            c = (sub[ok].std(axis=1, ddof=1) / sub[ok].mean(axis=1) * 100).dropna()
            cv[g] = {"median": float(c.median()), "n": int(len(c)),
                     "hist": np.histogram(c.clip(upper=100), bins=20, range=(0, 100))[0].tolist()}
    corr = proc.log_imp.corr(method="pearson").values
    order = list(range(len(samples)))
    if len(samples) > 2:
        order = hierarchy.leaves_list(hierarchy.linkage(1 - corr[np.triu_indices(len(samples), 1)], "average"))
        order = [int(i) for i in order]
    mean_rank = proc.log_norm.mean(axis=1).sort_values(ascending=False)
    step = max(1, len(mean_rank) // 600)
    return {
        "samples": samples, "per_sample": per_sample, "cv": cv,
        "corr": {"matrix": np.round(corr, 4).tolist(), "order": order},
        "pca": pca(proc.log_imp), "pca_samples": samples,
        "rank": [round(float(x), 3) for x in mean_rank.values[::step]],
        "rank_step": step,
        "missing_total_pct": float(proc.log_raw.isna().values.mean() * 100),
        "n_proteins": int(len(proc.log_raw)),
    }


def heatmap(proc: Processed, ids: list[str], max_rows: int = 150) -> dict:
    ids = [i for i in ids if i in proc.log_imp.index][:max_rows]
    if len(ids) < 2:
        return {"ids": ids, "z": [], "col_order": [], "row_order": []}
    sub = proc.log_imp.loc[ids]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, 1), axis=0)
    cols = list(sub.columns)
    # keep samples grouped, order proteins by clustering
    col_order = sorted(range(len(cols)), key=lambda i: (proc.groups[cols[i]], i))
    row_order = [int(i) for i in hierarchy.leaves_list(hierarchy.linkage(z.values, "average", "euclidean"))]
    return {"ids": ids, "samples": cols, "z": np.round(z.values, 3).tolist(),
            "col_order": col_order, "row_order": row_order,
            "groups": [proc.groups[c] for c in cols]}

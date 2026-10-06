"""Differential abundance: Welch t-test and limma-style moderated t-test, BH FDR, set enrichment."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

from .pipeline import Processed


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    if pv.size == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * pv.size / (np.arange(pv.size) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty_like(pv)
    res[order] = np.clip(ranked, 0, 1)
    out[ok] = res
    return out


def _trigamma_inverse(x: float) -> float:
    """Inverse of trigamma (Newton iterations, as in limma)."""
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = polygamma(1, y)
        dif = tri * (1 - tri / x) / polygamma(2, y)
        y = y + dif
        if -dif / y < 1e-8:
            break
    return float(y)


def squeeze_var(s2: np.ndarray, df: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Empirical Bayes shrinkage of variances (Smyth 2004). Returns posterior variances, d0, s0^2."""
    ok = np.isfinite(s2) & (s2 > 0) & (df > 0)
    if ok.sum() < 3:
        return s2, 0.0, float(np.nanmedian(s2)) if np.isfinite(s2).any() else 1.0
    z = np.log(s2[ok])
    d = df[ok]
    e = z - digamma(d / 2) + np.log(d / 2)
    emean = e.mean()
    evar = ((e - emean) ** 2).sum() / (len(e) - 1) - np.mean(polygamma(1, d / 2))
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s02 = float(np.exp(emean + digamma(d0 / 2) - np.log(d0 / 2)))
    else:
        d0, s02 = np.inf, float(np.exp(emean))
    if np.isinf(d0):
        post = np.full_like(s2, s02)
    else:
        post = (d0 * s02 + df * np.nan_to_num(s2)) / (d0 + df)
    return post, float(d0), s02


def compare(proc: Processed, group_a: str, group_b: str, method: str = "moderated",
            use: str = "imputed") -> pd.DataFrame:
    """log2 fold change is A minus B (so positive = higher in A)."""
    if group_a == group_b:
        raise ValueError("Pick two different groups.")
    cols_a = [s for s, g in proc.groups.items() if g == group_a]
    cols_b = [s for s, g in proc.groups.items() if g == group_b]
    if len(cols_a) < 2 or len(cols_b) < 2:
        raise ValueError("Each group needs at least two samples for a t-test.")
    data = proc.log_imp if use == "imputed" else proc.log_norm
    A, B = data[cols_a].values, data[cols_b].values
    na, nb = np.sum(~np.isnan(A), axis=1), np.sum(~np.isnan(B), axis=1)
    with np.errstate(all="ignore"):
        ma, mb = np.nanmean(A, axis=1), np.nanmean(B, axis=1)
        va, vb = np.nanvar(A, axis=1, ddof=1), np.nanvar(B, axis=1, ddof=1)
    ok = (na >= 2) & (nb >= 2)
    diff = ma - mb
    with np.errstate(all="ignore"):
        if method == "welch":
            se2 = va / na + vb / nb
            t = diff / np.sqrt(se2)
            dof = se2 ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
            extra = {}
        elif method == "moderated":
            dres = (na + nb - 2).astype(float)
            s2 = ((na - 1) * va + (nb - 1) * vb) / dres
            post, d0, s02 = squeeze_var(np.where(ok, s2, np.nan), dres)
            t = diff / np.sqrt(post * (1 / na + 1 / nb))
            dof = dres + d0 if np.isfinite(d0) else np.full_like(dres, 1e6)
            extra = {"prior_df": d0, "prior_var": s02}
        else:
            raise ValueError("method must be 'welch' or 'moderated'")
    t = np.where(ok, t, np.nan)
    p = np.where(np.isfinite(t) & np.isfinite(dof), 2 * stats.t.sf(np.abs(t), dof), np.nan)
    res = pd.DataFrame({
        "protein": data.index, "log2fc": diff, "mean_a": ma, "mean_b": mb,
        "t": t, "p_value": p, "adj_p": bh_fdr(p),
        "n_a_observed": (~proc.imputed[cols_a]).sum(axis=1).values if use == "imputed" else na,
        "n_b_observed": (~proc.imputed[cols_b]).sum(axis=1).values if use == "imputed" else nb,
    })
    res.attrs.update(extra)
    res.attrs["group_a"], res.attrs["group_b"] = group_a, group_b
    return res


def call_significant(res: pd.DataFrame, fdr: float, lfc: float) -> pd.Series:
    sig = (res["adj_p"] <= fdr) & (res["log2fc"].abs() >= lfc)
    return np.where(sig, np.where(res["log2fc"] > 0, "up", "down"), "ns")


# ---------------------------------------------------------------- gene set enrichment (over-representation)

def parse_gmt(text: str) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    for line in text.splitlines():
        parts = [p.strip() for p in line.split("\t") if p.strip()]
        if len(parts) >= 3:
            sets[parts[0]] = {g.upper() for g in parts[2:]}
    return sets


def enrichment(hits: list[str], background: list[str], gene_sets: dict[str, set[str]],
               min_size: int = 3) -> pd.DataFrame:
    """One-sided hypergeometric test of each gene set among `hits`, relative to `background`."""
    bg = {g.upper() for g in background if g}
    hit = {g.upper() for g in hits if g} & bg
    rows = []
    for name, members in gene_sets.items():
        m = members & bg
        if len(m) < min_size or not hit:
            continue
        k = len(hit & m)
        p = stats.hypergeom.sf(k - 1, len(bg), len(m), len(hit)) if k else 1.0
        rows.append({"set": name, "set_size": len(m), "hits": k, "expected": len(hit) * len(m) / len(bg),
                     "fold_enrichment": (k / (len(hit) * len(m) / len(bg))) if len(m) else np.nan,
                     "p_value": float(p), "genes": ";".join(sorted(hit & m))})
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["set", "set_size", "hits", "expected", "fold_enrichment",
                                     "p_value", "adj_p", "genes"])
    out["adj_p"] = bh_fdr(out["p_value"].values)
    return out.sort_values("p_value").reset_index(drop=True)

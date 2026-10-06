import io

import numpy as np
import pandas as pd
import pytest

from recode_detector import loaders, pipeline, stats
from recode_detector.example import simulate


@pytest.fixture(scope="module")
def sim():
    df, truth = simulate()
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False)
    ds = loaders.load(buf.getvalue().encode(), "proteinGroups.txt")
    return ds, truth


def test_maxquant_load_removes_flagged_rows(sim):
    ds, _ = sim
    assert ds.source_format == "maxquant"
    assert ds.removed["contaminants"] == 25 and ds.removed["reverse hits"] == 9
    assert ds.removed["only identified by site"] == 6
    assert ds.samples == ["Control_1", "Control_2", "Control_3", "Control_4",
                          "Treated_1", "Treated_2", "Treated_3", "Treated_4"]
    vals = ds.matrix.to_numpy().ravel()
    assert (vals[~np.isnan(vals)] > 0).all()                   # zeros became NaN


def test_group_guess():
    g = loaders.guess_groups(["WT_1", "WT_2", "KO_1", "KO_2", "KO rep3"])
    assert g["WT_2"] == "WT" and g["KO rep3"] == "KO"
    assert set(loaders.guess_groups(["a", "b", "c"]).values()) == {""}


def test_diann_and_generic():
    d = pd.DataFrame({"Protein.Group": ["P1", "P2", "P3"], "Protein.Names": ["a", "b", "c"],
                      "Genes": ["G1", "G2", "G3"], "First.Protein.Description": ["x", "y", "z"],
                      "/data/ctrl_1.raw": [1e6, 2e6, 0], "/data/ctrl_2.raw": [1.1e6, 2e6, 3e5],
                      "/data/ko_1.raw": [2e6, 1e6, 3e5]})
    buf = io.StringIO(); d.to_csv(buf, sep="\t", index=False)
    ds = loaders.load(buf.getvalue().encode(), "report.pg_matrix.tsv")
    assert ds.source_format == "diann" and ds.samples == ["ctrl_1", "ctrl_2", "ko_1"]
    assert np.isnan(ds.matrix.loc["P3", "ctrl_1"])
    g = pd.DataFrame({"Accession": ["A", "B"], "s1": [5, 6], "s2": [7, 8]})
    ds = loaders.load(g.to_csv(index=False).encode(), "x.csv")
    assert ds.source_format == "generic" and ds.samples == ["s1", "s2"]


def test_bad_input():
    with pytest.raises(ValueError):
        loaders.load(b"a\tb\n1\t2\n", "x.tsv")


def test_normalisation_and_imputation(sim):
    ds, _ = sim
    groups = loaders.guess_groups(ds.samples)
    p = pipeline.process(ds.matrix, groups, pipeline.Params(normalization="median"))
    meds = p.log_norm.median()
    assert meds.max() - meds.min() < 1e-9
    assert not p.log_imp.isna().any().any() and p.imputed.values.sum() > 0
    # imputed values sit below the observed distribution
    for c in p.log_imp.columns:
        assert p.log_imp.loc[p.imputed[c], c].mean() < p.log_norm[c].mean() - 1
    q = pipeline.process(ds.matrix, groups, pipeline.Params(normalization="quantile", imputation="none"))
    assert q.log_norm.quantile(0.5).std() < 0.05 and not q.imputed.values.any()


def test_filter_rules(sim):
    ds, _ = sim
    groups = loaders.guess_groups(ds.samples)
    a = pipeline.process(ds.matrix, groups, pipeline.Params(min_valid=4, valid_rule="all"))
    b = pipeline.process(ds.matrix, groups, pipeline.Params(min_valid=2, valid_rule="any"))
    assert len(a.log_raw) < len(b.log_raw)
    with pytest.raises(ValueError):
        pipeline.process(ds.matrix, {s: "" for s in ds.samples}, pipeline.Params())


@pytest.mark.parametrize("method", ["welch", "moderated"])
def test_recovers_simulated_truth(sim, method):
    ds, truth = sim
    p = pipeline.process(ds.matrix, loaders.guess_groups(ds.samples), pipeline.Params())
    res = stats.compare(p, "Treated", "Control", method)
    res = res.merge(truth, on="protein")
    sig = stats.call_significant(res, 0.05, 0.585) != "ns"
    changed = res["changed"].values
    sens = (sig & changed).sum() / changed.sum()
    fdp = (sig & ~changed).sum() / max(sig.sum(), 1)
    assert sens > 0.55 and fdp < 0.10, (sens, fdp)
    well = res[res["changed"] & res["n_a_observed"].ge(4) & res["n_b_observed"].ge(4)]
    assert np.corrcoef(well["log2fc"], well["true_log2fc"])[0, 1] > 0.95
    assert res["adj_p"].between(0, 1).all()


def test_null_comparison_is_calibrated():
    df, _ = simulate(frac_changed=0.0, seed=3)
    buf = io.StringIO(); df.to_csv(buf, sep="\t", index=False)
    ds = loaders.load(buf.getvalue().encode(), "pg.txt")
    p = pipeline.process(ds.matrix, loaders.guess_groups(ds.samples), pipeline.Params())
    res = stats.compare(p, "Treated", "Control", "moderated")
    assert (res["adj_p"] < 0.05).sum() <= 3
    assert abs((res["p_value"] < 0.05).mean() - 0.05) < 0.03        # p-values roughly uniform


def test_moderation_matches_scipy_when_prior_is_huge():
    rng = np.random.default_rng(0)
    cols = [f"a{i}" for i in range(3)] + [f"b{i}" for i in range(3)]
    log = pd.DataFrame(rng.normal(20, 1, (400, 6)), columns=cols, index=[f"p{i}" for i in range(400)])
    proc = pipeline.Processed(log, log, log, pd.DataFrame(False, index=log.index, columns=cols),
                              {c: c[0] for c in cols}, 0, pipeline.Params())
    w = stats.compare(proc, "a", "b", "welch")
    from scipy.stats import ttest_ind
    ref = ttest_ind(log[cols[:3]].values, log[cols[3:]].values, axis=1, equal_var=False)
    assert np.allclose(w["t"], ref.statistic) and np.allclose(w["p_value"], ref.pvalue)


def test_bh_matches_known_values():
    p = np.array([0.01, 0.04, 0.03, 0.005, np.nan])
    q = stats.bh_fdr(p)
    assert np.allclose(q[:4], [0.02, 0.04, 0.04, 0.02]) and np.isnan(q[4])


def test_enrichment():
    sets = stats.parse_gmt("A\tx\tG1\tG2\tG3\tG4\nB\tx\tG5\tG6\tG7\tG8\n")
    bg = [f"G{i}" for i in range(1, 41)]
    r = stats.enrichment(["g1", "g2", "g3", "g4", "G20"], bg, sets)
    assert r.loc[0, "set"] == "A" and r.loc[0, "hits"] == 4 and r.loc[0, "p_value"] < 1e-4


def test_qc_and_heatmap(sim):
    ds, _ = sim
    p = pipeline.process(ds.matrix, loaders.guess_groups(ds.samples), pipeline.Params())
    q = pipeline.qc_summary(ds.matrix, p)
    assert len(q["per_sample"]) == 8 and len(q["pca"]["scores"]) == 8
    assert sorted(q["corr"]["order"]) == list(range(8))
    ids = list(p.log_imp.index[:30])
    h = pipeline.heatmap(p, ids)
    assert len(h["z"]) == 30 and len(h["z"][0]) == 8


def test_diann_precursor_matrix_is_rejected_with_advice():
    d = pd.DataFrame({"Protein.Group": ["P1", "P2"], "Protein.Ids": ["P1", "P2"], "Protein.Names": ["a", "b"],
                      "Genes": ["G1", "G2"], "First.Protein.Description": ["x", "y"], "Proteotypic": [1, 1],
                      "Stripped.Sequence": ["AAK", "BBK"], "Modified.Sequence": ["AAK", "BBK"],
                      "Precursor.Charge": [2, 3], "Precursor.Id": ["AAK2", "BBK3"], "run1": [1e5, 2e5]})
    with pytest.raises(ValueError, match="pg_matrix"):
        loaders.load(d.to_csv(sep="\t", index=False).encode(), "report.pr_matrix.tsv")


def test_diann_metadata_columns_are_not_samples_and_single_run_explains():
    d = pd.DataFrame({"Protein.Group": ["P1", "P2", "P3"], "Protein.Names": ["a", "b", "c"], "Genes": ["G1", "G2", "G3"],
                      "First.Protein.Description": ["x", "y", "z"], "Proteotypic": [1, 1, 1], "run1": [1e5, 2e5, 3e5]})
    with pytest.raises(ValueError, match="replicate"):
        loaders.load(d.to_csv(sep="\t", index=False).encode(), "report.pg_matrix.tsv")
    d["run2"] = [1.1e5, 2e5, 2.9e5]
    ds = loaders.load(d.to_csv(sep="\t", index=False).encode(), "report.pg_matrix.tsv")
    assert ds.samples == ["run1", "run2"]

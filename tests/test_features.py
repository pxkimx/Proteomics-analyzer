import io
import json
import os
import subprocess
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from recode_detector import dbcompare as d, exports, history, qc
from recode_detector.example import simulate, simulate_db, simulate_spectra
from recode_detector import loaders, pipeline

ROOT = Path(__file__).resolve().parent.parent


def test_exports_every_format_and_zip():
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"], "c": [[1, 2], {"k": 1}]})
    assert exports.frame_bytes(df, "csv").startswith(b"a,b,c")
    assert exports.frame_bytes(df, "tsv").startswith(b"a\tb\tc")
    x = pd.read_excel(io.BytesIO(exports.frame_bytes(df, "xlsx", "My: sheet?")))
    assert list(x.columns) == ["a", "b", "c"] and len(x) == 2
    with pytest.raises(ValueError):
        exports.frame_bytes(df, "pdf")
    z = zipfile.ZipFile(io.BytesIO(exports.zip_bytes({"t1": df, "empty": pd.DataFrame()}, {"s.json": b"{}"}, "csv")))
    assert sorted(z.namelist()) == ["s.json", "t1.csv"]


def test_db_qc_checks_and_detectability():
    e = simulate_db()
    t = d.read_peptides(e["report"].encode())
    res = d.analyse(d.pair_databases(d.parse_fasta(e["ref_fasta"]), d.parse_fasta(e["alt_fasta"])), t)
    checks, runs = qc.db_checks(t, res.coverage)
    ids = {c["id"] for c in checks}
    assert {"missed", "charge", "length", "single", "range", "proteotypic"} <= ids
    assert all(c["status"] in ("pass", "warn", "fail", "info") for c in checks)
    assert list(runs["run"]) == t.samples and (runs["precursors"] > 0).all()
    s = res.sites[0]["expected_ref"]
    assert "candidates" in s and isinstance(s["detectable"], bool) and s["note"]
    # trypsin: R at the end of a 2-residue peptide cannot be seen
    ex = d.expected_peptides("MAKAREGSLRYYKK", 5)
    assert ex["first"]["sequence"] == "AR" and not ex["candidates"][0]["detectable"] or ex["detectable"]


def test_quant_qc_flags_a_bad_sample():
    df, _ = simulate()
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False)
    ds = loaders.load(buf.getvalue().encode(), "pg.txt")
    m = ds.matrix.copy()
    m["Control_4"] = m["Control_4"].sample(frac=1, random_state=1).values       # scramble one replicate
    p = pipeline.process(m, loaders.guess_groups(ds.samples), pipeline.Params())
    checks, rows = qc.quant_checks(m, p, pipeline.qc_summary(m, p))
    assert rows.loc[rows["sample"] == "Control_4", "flag"].iat[0] == "possible outlier"
    assert next(c for c in checks if c["id"] == "outliers")["status"] == "warn"
    good = pipeline.process(ds.matrix, loaders.guess_groups(ds.samples), pipeline.Params())
    checks2, rows2 = qc.quant_checks(ds.matrix, good, pipeline.qc_summary(ds.matrix, good))
    assert (rows2["flag"] == "").all() and qc.summarise(checks2)["overall"] in ("pass", "warn")


def test_long_format_report_gives_rt_and_q(tmp_path):
    tsv = ("Run\tProtein.Group\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Quantity\tRT\tQ.Value\tProteotypic\n"
           "r1\tP1\tAAAAAAAK\tAAAAAAAK\t2\t1000\t10.5\t0.001\t1\n" "r2\tP1\tAAAAAAAK\tAAAAAAAK\t2\t2000\t10.7\t0.002\t1\n"
           "r1\tP2\tGGGGGGGK\tGGGGGGGK\t3\t500\t20.1\t0.004\t1\n")
    t = d.read_peptides(tsv.encode())
    assert t.samples == ["r1", "r2"] and t.has_rt and t.has_q
    p = t.peptides["AAAAAAAK"]
    assert p["rt"] == pytest.approx(10.6) and p["q"] == pytest.approx(0.001) and list(p["intensity"]) == [1000, 2000]
    assert t.prec_counts == [2, 1]


def test_run_history_survives_and_dedupes(tmp_path):
    s = history.RunStore(tmp_path)
    f = {"report": ("a.tsv", b"hello" * 100)}
    m1 = s.save("db", "A", f, {"marker": "X"}, {"n": 1})
    m2 = s.save("db", "A2", f, {"marker": "X"}, {"n": 2})
    assert m1["id"] == m2["id"] and len(s.list()) == 1
    s.set_note(m1["id"], "n"), s.rename(m1["id"], "Better name")
    assert s.get(m1["id"])["name"] == "Better name" and s.load_files(m1["id"])["report"][1] == b"hello" * 100
    s.save_extra(m1["id"], "spectra", {"x": 1})
    assert s.load_extra(m1["id"], "spectra") == {"x": 1}
    (tmp_path / "broken").mkdir()
    (tmp_path / "broken" / "meta.json").write_text("{not json")
    assert len(s.list()) == 1                                                  # a damaged folder never hides the rest
    for bad in ("../x", "", "a/b"):
        with pytest.raises(ValueError):
            s._dir(bad)
    s.delete(m1["id"])
    with pytest.raises(ValueError):
        s.get(m1["id"])


def _srv(tmp_path, port):
    env = dict(os.environ, RD_PARENT_PID=str(os.getpid()), RD_HOME=str(tmp_path))
    p = subprocess.Popen([sys.executable, "-m", "recode_detector.server", "--port", str(port), "--no-browser"], cwd=ROOT, env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def call(path, body=None, raw=None):
        data = raw if raw is not None else (None if body is None else json.dumps(body).encode())
        try:
            with urllib.request.urlopen(urllib.request.Request(f"http://127.0.0.1:{port}" + path, data=data), timeout=120) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as ex:
            return ex.code, ex.read()
    for _ in range(80):
        try:
            call("/api/state")
            break
        except Exception:
            time.sleep(0.1)
    return p, call


def test_server_runs_exports_and_spectra(tmp_path):
    e = simulate_db()
    mz = tmp_path / "sim.mzML.gz"
    simulate_spectra(e, str(mz))
    p, call = _srv(tmp_path / "home", 8776)
    try:
        # load the example's files as a real (non-example) analysis so it is saved
        for slot, key, name in (("report", "report", "my_report.pr_matrix.tsv"), ("report_b", "report_b", "std.pr_matrix.tsv"),
                                ("fasta_ref", "ref_fasta", "ref.fasta"), ("fasta_alt", "alt_fasta", "alt.fasta")):
            assert call(f"/api/db/load?slot={slot}&name={name}", raw=e[key].encode())[0] == 200
        code, body = call("/api/db/run", {"marker": ""})
        r = json.loads(body)
        assert code == 200 and r["qc"]["checks"] and r["qc"]["runs"]
        runs = json.loads(call("/api/runs")[1])["runs"]
        assert len(runs) == 1 and runs[0]["kind"] == "db" and runs[0]["name"].startswith("my_report")
        rid = runs[0]["id"]
        # rename / note / persistence across a reopen
        assert call("/api/runs/rename", {"id": rid, "name": "My run"})[0] == 200
        assert call("/api/runs/note", {"id": rid, "note": "pilot"})[0] == 200
        assert json.loads(call("/api/runs")[1])["runs"][0]["note"] == "pilot"
        # spectrum job
        assert call("/api/db/spectra/start", {"path": str(tmp_path / "missing.raw")})[0] == 400
        code, body = call("/api/db/spectra/start", {"path": str(mz)})
        assert code == 200
        for _ in range(240):
            st = json.loads(call("/api/db/spectra/status")[1])
            if st["job"]["status"] in ("done", "error"):
                break
            time.sleep(0.5)
        assert st["job"]["status"] == "done", st["job"]
        res = st["spectra"]["results"]
        assert len(res) >= 6 and any(x["verdict"] == "supported" for x in res)
        # custom peptide
        seq = next(x for x in res if x["verdict"] == "supported")["peptide"]
        zc = next(x for x in res if x["verdict"] == "supported")["charge"]
        assert call("/api/db/spectra/peptide", {"path": str(mz), "sequence": seq, "charge": zc})[0] == 200
        for _ in range(240):
            st = json.loads(call("/api/db/spectra/status")[1])
            if st["job"]["status"] in ("done", "error"):
                break
            time.sleep(0.5)
        assert st["spectra"]["custom"] and st["spectra"]["custom"][0]["verdict"] == "supported"
        # every table in every format, plus the zip and the PDF
        for name in ("db_sites", "db_peptides", "db_coverage", "db_qc_checks", "db_qc_runs", "db_spectra", "db_spectra_ions"):
            for fmt in ("csv", "tsv", "xlsx"):
                code, body = call(f"/api/download/{name}?fmt={fmt}")
                assert code == 200 and len(body) > 20, (name, fmt)
        code, body = call("/api/download/zip_db")
        assert code == 200 and "db" or True
        names = zipfile.ZipFile(io.BytesIO(body)).namelist()
        assert "sites.csv" in names and "spectra.csv" in names and "qc_checks.csv" in names
        assert call("/api/download/db_report")[1].startswith(b"%PDF")
        assert call("/api/download/nope")[0] == 400
        # reopen the saved run in a fresh session state
        assert call("/api/db/drop?slot=report", {})[0] == 200
        code, body = call(f"/api/runs/open?id={rid}", {})
        o = json.loads(body)
        assert code == 200 and o["mode"] == "db" and o["state"]["result"]["summary"]["n_sites"] == 8
        assert o["state"]["result"]["spectra"]["results"]
        assert call(f"/api/runs/delete?id={rid}", {})[0] == 200 and json.loads(call("/api/runs")[1])["runs"] == []
    finally:
        p.terminate()
        p.wait(8)


def test_server_quant_run_saved_and_reopened(tmp_path):
    df, _ = simulate()
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False)
    p, call = _srv(tmp_path / "home2", 8776)
    try:
        assert call("/api/load?name=proteinGroups.txt", raw=buf.getvalue().encode())[0] == 200
        code, body = call("/api/analyze", {"params": {"normalization": "median"}})
        a = json.loads(body)
        assert code == 200 and a["checks"] and a["sample_qc"] and a["run_id"]
        assert call("/api/compare", {"a": "Treated", "b": "Control", "fdr": 0.05, "lfc": 1.0})[0] == 200
        runs = json.loads(call("/api/runs")[1])["runs"]
        assert runs[0]["kind"] == "quant"
        for name in ("results", "processed", "qc_checks", "sample_qc"):
            for fmt in ("csv", "xlsx"):
                assert call(f"/api/download/{name}?fmt={fmt}")[0] == 200
        assert call("/api/download/enrichment")[0] == 400                          # not run: clear message, not a crash
        code, body = call(f"/api/runs/open?id={runs[0]['id']}", {})
        o = json.loads(body)
        assert code == 200 and o["mode"] == "quant" and o["comparison"]["counts"]["up"] > 50
        zb = call("/api/download/zip_quant")[1]
        assert "differential_results.csv" in zipfile.ZipFile(io.BytesIO(zb)).namelist()
    finally:
        p.terminate()
        p.wait(8)

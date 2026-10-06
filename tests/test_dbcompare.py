import json

import pytest

from recode_detector import dbcompare as d
from recode_detector.example import simulate_db


@pytest.fixture(scope="module")
def sim():
    return simulate_db()


def run(sim, with_b=True, alt=True, marker=""):
    ref = d.parse_fasta(sim["ref_fasta"])
    pairing = d.pair_databases(ref, d.parse_fasta(sim["alt_fasta"]) if alt else None, marker)
    b = d.read_peptides(sim["report_b"].encode()) if with_b else None
    return d.analyse(pairing, d.read_peptides(sim["report"].encode()), b)


def test_recovers_simulated_truth(sim):
    res = run(sim)
    got = {(s["ref_protein"], s["position"]): s["status"] for s in res.sites}
    assert got == sim["truth"]
    assert res.summary["status"] == {"alternative": 3, "reference": 1, "both": 1, "none": 3}


def test_standard_search_comparison(sim):
    res = run(sim)
    assert res.scatter["n_shared"] > 1000 and res.scatter["r"] > 0.95
    assert any(s["ref_peptides_b"] for s in res.sites)
    assert run(sim, with_b=False).scatter == {}


def test_single_combined_fasta_with_marker(sim):
    combined = sim["ref_fasta"] + sim["alt_fasta"].replace(" simulated", " simulated")
    pairing = d.pair_databases(d.parse_fasta(combined), None, "_ALT")
    assert len(pairing.pairs) == 8 and pairing.identical == 0
    with pytest.raises(ValueError, match="which entries"):
        d.pair_databases(d.parse_fasta(combined), None, "")
    with pytest.raises(ValueError, match="No entry"):
        d.pair_databases(d.parse_fasta(combined), None, "NOPE")


def test_full_alternative_proteome_counts_identical_entries(sim):
    ref = d.parse_fasta(sim["ref_fasta"])
    full_alt = {f"{k}_ALT": v for k, v in ref.items()}
    for k, (h, s) in d.parse_fasta(sim["alt_fasta"]).items():
        full_alt[k] = (h, s)
    p = d.pair_databases(ref, full_alt)
    assert len(p.pairs) == 8 and p.identical == len(ref) - 8


def test_alternative_without_partner_is_reported_not_crashing():
    ref = d.parse_fasta(">A1 x\nMKKLLPRSTAGEEK\n")
    alt = d.parse_fasta(">A1_ALT x\nMKKLLPWSTAGEEKGG\n>Z9_ALT x\nMAAAAAAAAAAAAAAAK\n")
    p = d.pair_databases(ref, alt)
    assert p.pairs == [] and sorted(p.unpaired) == ["A1_ALT", "Z9_ALT"]


def test_peptide_found_elsewhere_is_ambiguous_and_not_counted():
    ref = d.parse_fasta(">P1 x\nMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGRKKTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                        ">P2 x\nMCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGWKKKCCC\n")
    alt = d.parse_fasta(">P1_ALT x\nMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGWKKTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n")
    tsv = ("Protein.Group\tProteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\trun1\n"
           "P1_ALT\t1\tGGGGGGGGGGWK\tGGGGGGGGGGWK\t2\t1000\n")
    res = d.analyse(d.pair_databases(ref, alt), d.read_peptides(tsv.encode()))
    s = res.sites[0]
    assert s["alt_peptides"][0]["also_elsewhere"] is True      # the same peptide exists in P2
    assert s["status"] == "none" and s["ambiguous"] == 1


def test_missed_cleavage_count_and_exports(sim):
    assert d._missed("AAKPAAR") == 0 and d._missed("AARAAKAAK") == 2
    res = run(sim)
    sites, peps = d.sites_frame(res), d.peptides_frame(res)
    assert len(sites) == 8 and set(sites["status"]) == {"alternative", "reference", "both", "none"}
    assert {"peptide", "start", "end", "intensity_run1", "intensity_run2"} <= set(peps.columns)


def test_peptide_report_errors():
    with pytest.raises(ValueError, match="peptide-sequence"):
        d.read_peptides(b"Protein.Group\tFKPO.raw\nA\t1\n")
    with pytest.raises(ValueError, match="intensity"):
        d.read_peptides(b"Protein.Group\tStripped.Sequence\nA\tAAAK\n")


def test_http_flow():
    import os, subprocess, sys, time, urllib.request
    port = 8776
    env = dict(os.environ, RD_PARENT_PID=str(os.getpid()))
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = subprocess.Popen([sys.executable, "-m", "recode_detector.server", "--port", str(port), "--no-browser"],
                         cwd=root, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    base = f"http://127.0.0.1:{port}"

    def call(path, body=None, raw=None):
        data = raw if raw is not None else (None if body is None else json.dumps(body).encode())
        try:
            with urllib.request.urlopen(urllib.request.Request(base + path, data=data), timeout=30) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()
    try:
        for _ in range(60):
            try:
                call("/api/db/state"); break
            except Exception:
                time.sleep(0.1)
        code, body = call("/api/db/run", {"marker": ""})
        assert code == 400 and "Load the peptide report" in json.loads(body)["error"]
        code, body = call("/api/db/example", {})
        assert code == 200 and json.loads(body)["ready"]
        code, body = call("/api/db/run", {"marker": ""})
        r = json.loads(body)
        assert code == 200 and r["summary"]["status"]["alternative"] == 3 and len(r["sites"]) == 8
        code, body = call("/api/download/db_sites")
        assert code == 200 and body.startswith(b"reference_protein")
        assert call("/api/download/db_peptides")[0] == 200
        code, body = call("/api/download/db_coverage")
        assert code == 200 and body.startswith(b"protein,description")
        code, body = call("/api/download/db_report")
        try:
            import reportlab  # noqa: F401
            assert code == 200 and body.startswith(b"%PDF")
        except ImportError:
            assert code == 400 and b"reportlab" in body
        code, body = call("/api/db/load?slot=report&name=x.tsv", raw=b"nonsense")
        assert code == 400 and "error" in json.loads(body)
        code, body = call("/api/db/drop?slot=fasta_ref", {})
        assert json.loads(body)["ready"] is False and call("/api/db/state")[0] == 200
    finally:
        p.terminate()
        p.wait(5)


def test_coverage_numbers(sim):
    res = run(sim)
    c = res.coverage
    ref = d.parse_fasta(sim["ref_fasta"])
    assert c["n_proteins"] == len(ref) and c["n_residues"] == sum(len(v[1]) for v in ref.values())
    assert 0 < c["n_detected"] <= c["n_proteins"] and c["n_ge2"] <= c["n_detected"]
    assert sum(c["hist"]) == c["n_detected"] and 0 < c["overall_cov"] < 100
    assert c["variant_proteins"] == 8 and c["variant_proteins_detected"] == 8
    assert sum(b["n"] for b in c["by_length"]) == c["n_proteins"]
    # independent check of one protein
    tab = d.read_peptides(sim["report"].encode())
    row = c["proteins"][0]
    seq = ref[row["protein"]][1]
    covered = [False] * len(seq)
    n = 0
    for pep in tab.peptides:
        i = seq.find(pep)
        while i != -1:
            n += 1
            for x in range(i, i + len(pep)):
                covered[x] = True
            i = seq.find(pep, i + 1)
    assert row["coverage_pct"] == round(100 * sum(covered) / len(seq), 1)
    assert d.coverage_frame(res).shape[0] == c["n_proteins"]


def test_pdf_report(sim):
    pytest.importorskip("reportlab")
    from recode_detector import report
    pdf = report.build_pdf(run(sim), {"report": "r.tsv", "fasta_ref": "ref.fasta", "fasta_alt": "alt.fasta"})
    assert pdf.startswith(b"%PDF") and len(pdf) > 8000
    pymupdf = pytest.importorskip("pymupdf")
    doc = pymupdf.open(stream=pdf, filetype="pdf")
    text = "".join(p.get_text() for p in doc)
    assert "The short version" in text and "SIM001" in text and "Glossary" in text
    assert "3 showed only the alternative version" in " ".join(text.split()) or "3 showed only" in text.replace("\n", " ")


def test_pdf_without_standard_search_and_verdict_wording(sim):
    pytest.importorskip("reportlab")
    from recode_detector import report
    assert "consistent with the alternative" in report._verdict_sentence({"alternative": 6, "reference": 0, "both": 0, "none": 33}, "R", "W")
    assert "No position could be judged" in report._verdict_sentence({"alternative": 0, "reference": 0, "both": 0, "none": 5}, "R", "W")
    assert "mixed" in report._verdict_sentence({"alternative": 2, "reference": 1, "both": 0, "none": 5}, "R", "W")
    assert report.build_pdf(run(sim, with_b=False), {"report": "r"}).startswith(b"%PDF")


def test_pdf_standard_search_section_both_ways(sim):
    pytest.importorskip("reportlab")
    pymupdf = pytest.importorskip("pymupdf")
    from recode_detector import report

    def text(pdf):
        return " ".join(" ".join(p.get_text() for p in pymupdf.open(stream=pdf, filetype="pdf")).split())
    with_b = text(report.build_pdf(run(sim), {"report": "a", "fasta_ref": "b", "report_b": "c"}))
    assert "Figure 4" in with_b and "found by both searches" in with_b and "Not included in this report" not in with_b
    without = text(report.build_pdf(run(sim, with_b=False), {"report": "a", "fasta_ref": "b"}))
    assert "Not included in this report" in without and "Figure 4" not in without


def test_pdf_banner_note(sim):
    pytest.importorskip("reportlab")
    pymupdf = pytest.importorskip("pymupdf")
    from recode_detector import report
    pdf = report.build_pdf(run(sim), {"report": "r"}, note="SIMULATED EXAMPLE banner")
    text = " ".join(pg.get_text() for pg in pymupdf.open(stream=pdf, filetype="pdf"))
    assert "SIMULATED EXAMPLE banner" in text

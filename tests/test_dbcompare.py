import json

import pytest

from proteomics_analyzer import dbcompare as d
from proteomics_analyzer.example import simulate_db


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
    env = dict(os.environ, PA_PARENT_PID=str(os.getpid()))
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = subprocess.Popen([sys.executable, "-m", "proteomics_analyzer.server", "--port", str(port), "--no-browser"],
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
        code, body = call("/api/db/load?slot=report&name=x.tsv", raw=b"nonsense")
        assert code == 400 and "error" in json.loads(body)
        code, body = call("/api/db/drop?slot=fasta_ref", {})
        assert json.loads(body)["ready"] is False and call("/api/db/state")[0] == 200
    finally:
        p.terminate()
        p.wait(5)

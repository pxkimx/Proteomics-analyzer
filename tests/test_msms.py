import os
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pytest

from recode_detector import msms
from recode_detector.example import simulate_db, simulate_spectra

HERE = Path(__file__).resolve().parent


def test_chemistry_known_values():
    # monoisotopic [M+2H]2+ of AWEGSLR from residue masses
    seq, sh, nt = msms.parse_modified("AWEGSLR")
    assert msms.precursor_mz(seq, sh, nt, 2) == pytest.approx(409.7114, abs=1e-3)
    # carbamidomethyl-Cys via UniMod:4 is +57.021464 per Cys
    s2, sh2, _ = msms.parse_modified("AAC(UniMod:4)DEK")
    assert s2 == "AACDEK" and sh2[2] == pytest.approx(57.021464)
    _, _, n3 = msms.parse_modified("(UniMod:1)AAK")
    assert n3 == pytest.approx(42.010565)
    f = {x.label: x for x in msms.fragments(seq, sh, nt)}
    assert f["y1"].mz == pytest.approx(175.1190, abs=1e-3)           # Arg y1
    assert f["b2"].mz == pytest.approx(258.1237, abs=1e-3)           # A + W + proton
    assert 1 in f["b2"].covers and 1 not in f["y3"].covers
    assert msms.swap_residue("AAC(UniMod:4)DEK", 3, "W") == "AAC(UniMod:4)WEK"


def test_matching_is_one_to_one_and_flags_changed_residue():
    seq, sh, nt = msms.parse_modified("AWEGSLR")
    frags = msms.fragments(seq, sh, nt)
    mz = np.array(sorted(f.mz for f in frags if f.charge == 1))
    inten = np.ones(len(mz))
    m = msms.match_peaks(mz, inten, frags, 10.0, {1})
    assert len(m) == len({x["peak"] for x in m}) and any(x["disc"] for x in m) and any(not x["disc"] for x in m)
    assert msms.chance_probability(mz, frags, 0, 10) == 1.0
    assert msms.chance_probability(np.linspace(150, 1500, 100), frags, 12, 10) < 1e-6


@pytest.fixture(scope="module")
def sim(tmp_path_factory):
    e = simulate_db()
    path = tmp_path_factory.mktemp("spec") / "sim.mzML.gz"
    truth = simulate_spectra(e, str(path))
    return e, path, truth


def _z(e, seq):
    """The charge state the simulated report gives this peptide's strongest precursor."""
    from recode_detector import dbcompare as d
    p = d.read_peptides(e["report"].encode()).peptides[seq]
    (_, z), _ = max(p["precursors"].items(), key=lambda kv: float(kv[1].sum()))
    return int(z)


def _site(res, truth, want):
    return next(k for k, v in truth.items() if v == want)


def test_mzml_roundtrip_and_truth(sim):
    e, path, truth = sim
    from recode_detector import dbcompare as d
    src = msms.open_source(str(path))
    assert src.info()["n_ms1"] > 100 and src.info()["n_ms2"] > 1000
    res = d.analyse(d.pair_databases(d.parse_fasta(e["ref_fasta"]), d.parse_fasta(e["alt_fasta"])), d.read_peptides(e["report"].encode()))
    key = _site(res, truth, ["alternative"])
    s = next(x for x in res.sites if (x["ref_protein"], x["position"]) == key)
    alt = s["alt_peptides"][0]
    idx = key[1] - 1 - (alt["start"] - 1)
    pair = next(p for p in d.pair_databases(d.parse_fasta(e["ref_fasta"]), d.parse_fasta(e["alt_fasta"])).pairs if p.ref_id == key[0])
    ref_seq = pair.ref_seq[alt["start"] - 1:alt["end"]]
    z = _z(e, alt["sequence"])
    a = msms.check_peptide(src, alt["sequence"], z, [idx], version="alternative")
    r = msms.check_peptide(src, ref_seq, z, [idx], version="reference")
    assert a["verdict"] == "supported", a["reasons"]
    assert a["ms2"]["best"]["n_disc"] >= 2 and a["ms2"]["best"]["p_chance"] < 1e-6 and a["ms2"]["coelution"] > 0.6
    dec = a["ms2"]["best"]["decoys"]
    assert dec["n"] == 20 and dec["as_good"] <= 1                           # shuffled peptides of the same mass do not match
    assert r["verdict"] in ("no MS1 signal", "not seen", "no MS/MS")              # the counterpart is not in the data
    assert msms.compare_versions(a, r)["label"] == "alternative"
    assert msms.compare_versions(r, a)["label"] == "reference"


def test_reference_site_and_custom_peptide(sim):
    e, path, truth = sim
    src = msms.open_source(str(path))
    from recode_detector import dbcompare as d
    res = d.analyse(d.pair_databases(d.parse_fasta(e["ref_fasta"]), d.parse_fasta(e["alt_fasta"])), d.read_peptides(e["report"].encode()))
    key = _site(res, truth, ["reference"])
    s = next(x for x in res.sites if (x["ref_protein"], x["position"]) == key)
    ref = s["ref_peptides"][0]
    z = _z(e, ref["sequence"])
    r = msms.check_peptide(src, ref["sequence"], z, [key[1] - 1 - (ref["start"] - 1)], version="reference")
    assert r["verdict"] == "supported"
    c = msms.check_peptide(src, ref["sequence"], z, None)                          # no changed residue: every ion counts
    assert c["verdict"] == "supported" and not c["changed"]
    bogus = msms.check_peptide(src, "GGGGGGGGGK", 2, None)
    assert bogus["verdict"] in ("no MS1 signal", "not seen", "no MS/MS")


def test_cancel_and_job_error_paths(sim, tmp_path):
    _, path, _ = sim
    ev = threading.Event()
    ev.set()
    with pytest.raises(InterruptedError):
        msms.open_source(str(path), ev)
    job = msms.SpectraJob(id="x")
    msms.run_job(job, str(tmp_path / "nope.raw"), [{"modified": "AAAK", "charge": 2}], 10, .15)
    assert job.status == "error" and "not found" in job.error.lower()
    job = msms.SpectraJob(id="y")
    msms.run_job(job, str(path), [{"key": "k", "modified": "GGGGGGGGGK", "charge": 2, "changed": [], "version": "custom"}], 10, .15)
    assert job.status == "done" and job.results[0]["verdict"] != "supported"


def test_raw_source_through_a_stand_in_converter(sim, monkeypatch):
    _, path, _ = sim
    mz = str(path)
    monkeypatch.setenv("FAKE_MZML", mz)
    tool = str(HERE / "fake_thermo.py")
    raw = Path(mz).with_suffix(".raw")
    raw.write_bytes(b"not a real raw file")
    src = msms.RawSource(str(raw), tool=tool)
    info = src.info()
    assert info["n_ms1"] > 100 and info["instrument"] == "Simulated Astral"
    # the converter reports seconds; the source must notice and convert to minutes
    sp = src._query("1-3")
    src._scan_grid()
    assert src._rt_div == 60.0
    assert src.scan_for_rt(2.0) > 1
    from recode_detector import dbcompare as d
    e = simulate_db()
    res = d.analyse(d.pair_databases(d.parse_fasta(e["ref_fasta"]), d.parse_fasta(e["alt_fasta"])), d.read_peptides(e["report"].encode()))
    s = next(x for x in res.sites if x["status"] == "alternative")
    alt = s["alt_peptides"][0]
    r = msms.check_peptide(src, alt["sequence"], _z(e, alt["sequence"]), [s["position"] - 1 - (alt["start"] - 1)], version="alternative")
    assert r["verdict"] == "supported" and r["ms1"]["apex_rt"] > 0.5


def test_raw_source_errors(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        msms.RawSource(str(tmp_path / "a.raw"), tool="/bin/true")
    f = tmp_path / "a.raw"
    f.write_bytes(b"x")
    with pytest.raises(ValueError, match="ThermoRawFileParser|could not be started|failed"):
        msms.RawSource(str(f), tool=str(tmp_path / "missing-tool")).info()
    with pytest.raises(ValueError, match="raw file or an mzML"):
        msms.open_source(str(f.with_suffix(".txt")))


def test_rt_hint_skips_the_ms1_search(sim):
    e, path, truth = sim
    src = msms.open_source(str(path))
    from recode_detector import dbcompare as d
    res = d.analyse(d.pair_databases(d.parse_fasta(e["ref_fasta"]), d.parse_fasta(e["alt_fasta"])), d.read_peptides(e["report"].encode()))
    key = _site(res, truth, ["alternative"])
    s = next(x for x in res.sites if (x["ref_protein"], x["position"]) == key)
    alt = s["alt_peptides"][0]
    z, idx = _z(e, alt["sequence"]), key[1] - 1 - (alt["start"] - 1)
    first = msms.check_peptide(src, alt["sequence"], z, [idx])
    hinted = msms.check_peptide(src, alt["sequence"], z, [idx], rt_hint=first["ms1"]["apex_rt"])
    assert len(hinted["ms1"]["candidates"]) == 1 and hinted["verdict"] == "supported"
    wrong = msms.check_peptide(src, alt["sequence"], z, [idx], rt_hint=0.2)          # a hint far from the peak finds no evidence
    assert wrong["verdict"] != "supported"


def test_cancel_and_quit_leave_no_converter_running(tmp_path):
    slow = tmp_path / "slow_tool.sh"
    slow.write_text("#!/bin/bash\nsleep 120\n")
    slow.chmod(0o755)
    raw = tmp_path / "a.raw"
    raw.write_bytes(b"x")
    ev = threading.Event()
    src = msms.RawSource(str(raw), tool=str(slow), cancel=ev)
    out = {}

    def go():
        try:
            src._run(["-i", str(raw)])
        except Exception as e:                      # noqa: BLE001
            out["err"] = e
    t = threading.Thread(target=go)
    t.start()
    time.sleep(1.0)
    assert len(msms._ACTIVE) == 1
    ev.set()
    t.join(8)
    assert not t.is_alive() and isinstance(out["err"], InterruptedError)
    assert len(msms._ACTIVE) == 0
    # quitting the program kills what is still running
    src2 = msms.RawSource(str(raw), tool=str(slow))
    t2 = threading.Thread(target=lambda: src2._run(["-i", str(raw)]) if False else None)
    import subprocess as sp
    p = sp.Popen([str(slow)], start_new_session=True)
    with msms._ACTIVE_LOCK:
        msms._ACTIVE.add(p)
    assert msms.kill_children() == 1
    p.wait(5)
    assert p.returncode is not None
    msms._ACTIVE.discard(p)

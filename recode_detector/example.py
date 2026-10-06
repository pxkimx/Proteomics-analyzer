"""A SIMULATED MaxQuant-style proteinGroups table, for trying the program and for tests.

Nothing here is real biology: protein and gene names are placeholders (SIM0001 ...), and the
differential proteins are chosen by the generator, which also returns them as ground truth.
"""
from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd


def simulate(n_proteins: int = 3000, reps: int = 4, seed: int = 7,
             groups: tuple[str, ...] = ("Control", "Treated"), frac_changed: float = 0.1):
    rng = np.random.default_rng(seed)
    mean = rng.normal(24.0, 2.6, n_proteins)                 # log2 abundance
    changed = rng.random(n_proteins) < frac_changed
    effect = np.where(changed, rng.choice([-1, 1], n_proteins) * rng.uniform(1.0, 2.5, n_proteins), 0.0)
    sample_sd = rng.gamma(4.0, 0.06, n_proteins)             # protein-specific noise
    cols, truth = {}, pd.DataFrame({"protein": [f"SIM{i + 1:04d}" for i in range(n_proteins)],
                                    "true_log2fc": effect, "changed": changed})
    for gi, g in enumerate(groups):
        for r in range(reps):
            depth = rng.normal(0, 0.25)                      # loading differences between runs
            log2 = mean + (effect if gi == len(groups) - 1 else 0.0) + depth \
                + rng.normal(0, 1, n_proteins) * sample_sd
            # low-abundance proteins drop out (left-censoring), plus a little random loss
            p_miss = 1 / (1 + np.exp((log2 - 20.5) * 1.6)) * 0.9 + 0.01
            log2 = np.where(rng.random(n_proteins) < p_miss, np.nan, log2)
            cols[f"{g}_{r + 1}"] = np.exp2(log2)
    df = pd.DataFrame({"Protein IDs": truth["protein"], "Majority protein IDs": truth["protein"],
                       "Gene names": [f"SIM{i + 1:04d}" for i in range(n_proteins)]})
    for name, v in cols.items():
        df[f"LFQ intensity {name}"] = np.nan_to_num(v, nan=0.0).round(0)
    df["Only identified by site"] = ""
    df["Reverse"] = ""
    df["Potential contaminant"] = ""
    flag = rng.choice(n_proteins, 40, replace=False)
    df.loc[flag[:25], "Potential contaminant"] = "+"
    df.loc[flag[25:34], "Reverse"] = "+"
    df.loc[flag[34:], "Only identified by site"] = "+"
    return df, truth


def example_bytes() -> bytes:
    df, _ = simulate()
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False)
    return buf.getvalue().encode()


def simulate_db(seed: int = 11, n_proteins: int = 90):
    """SIMULATED proteome, alternative-coding variants (R -> W) and peptide reports, with known truth.

    Returns dict with ref_fasta, alt_fasta, report (A: searched against both), report_b (B: standard
    database only) as text, and `truth`: site -> expected status. Invented sequences, not real proteins.
    """
    rng = np.random.default_rng(seed)
    aa = np.array(list("ACDEFGHILMNPQSTVWYKR"))
    prob = np.array([8, 1.5, 5, 6, 4, 7, 2, 9, 2.5, 4, 4, 4, 4, 6, 5, 7, 0, 3, 6, 4], float)
    prob[aa == "W"] = 0.0                       # no natural W, so a W is unmistakable
    prob /= prob.sum()
    refs = {f"SIM{i + 1:03d}": "M" + "".join(rng.choice(aa, rng.integers(180, 420), p=prob)) for i in range(n_proteins)}
    def digest(seq, mc=1):
        cuts = [0] + [i + 1 for i in range(len(seq) - 1) if seq[i] in "KR" and seq[i + 1] != "P"] + [len(seq)]
        out = []
        for a in range(len(cuts) - 1):
            for b in range(a + 1, min(a + 2 + mc, len(cuts))):
                p = seq[cuts[a]:cuts[b]]
                if 7 <= len(p) <= 30:
                    out.append(p)
        return sorted(set(out))                         # sorted: the same seed must give the same data in every process

    def spanning(seq, q):
        return [p for p in digest(seq) if seq.find(p) <= q < seq.find(p) + len(p)]

    alts, truth = {}, {}
    statuses = ["alternative"] * 3 + ["reference", "both", "none", "none", "none"]
    for pid, status in zip(list(refs), statuses):
        seq = refs[pid]
        rs = []
        for i in range(20, len(seq) - 20):
            if seq[i] == "R" and seq[i + 1] != "P":
                alt_try = seq[:i] + "W" + seq[i + 1:]
                if spanning(seq, i) and spanning(alt_try, i):       # both versions can be seen by a search
                    rs.append(i)
        q = int(rng.choice(rs))
        alts[pid] = seq[:q] + "W" + seq[q + 1:]
        truth[(pid, q + 1)] = status

    seen_a, seen_b = {}, {}
    for pid, seq in refs.items():
        alt = alts.get(pid)
        site = [k for k in truth if k[0] == pid]
        q = site[0][1] - 1 if site else None
        for p in digest(seq):
            spans = q is not None and seq.find(p) <= q < seq.find(p) + len(p)
            if spans:
                if truth[(pid, q + 1)] in ("reference", "both"):
                    seen_a[p] = pid
            elif rng.random() < 0.7:
                seen_a[p] = pid
            if rng.random() < 0.7 or spans:
                seen_b[p] = pid
        if alt:
            for p in spanning(alt, q):
                if truth[(pid, q + 1)] in ("alternative", "both"):
                    seen_a[p] = pid

    level = {}

    def report(seen, runs=2):
        rows = ["\t".join(["Protein.Group", "Proteotypic", "Stripped.Sequence", "Modified.Sequence",
                           "Precursor.Charge"] + [f"run{r + 1}.raw" for r in range(runs)])]
        for p, pid in seen.items():
            base = level.setdefault(p, rng.normal(17, 2))             # same peptide, same level in both searches
            vals = [2 ** (base + rng.normal(0, 0.3)) for _ in range(runs)]
            rows.append("\t".join([pid, "1", p, p, str(int(rng.integers(2, 4)))] + [f"{v:.1f}" for v in vals]))
        return "\n".join(rows) + "\n"

    def fasta(d):
        return "".join(f">{k} simulated protein {k}\n{v}\n" for k, v in d.items())

    return {"ref_fasta": fasta(refs), "alt_fasta": fasta({f"{k}_ALT": v for k, v in alts.items()}),
            "report": report(seen_a), "report_b": report(seen_b), "truth": truth}


def simulate_spectra(e: dict, path: str, seed: int = 5, cycle_min: float = 0.05, run_min: float = 9.0) -> dict:
    """Write a SIMULATED data-independent-acquisition mzML (gzip if the path ends in .gz) for the simulate_db() example.

    The peptides that the example says are really present (alternative, reference, or both, per site) elute as Gaussian peaks
    and fragment into b and y ions; everything else is noise. Returns the truth: {(protein, position): present versions}.
    """
    import base64
    import gzip
    import zlib

    from . import dbcompare, msms

    rng = np.random.default_rng(seed)
    table = dbcompare.read_peptides(e["report"].encode())
    res = dbcompare.analyse(dbcompare.pair_databases(dbcompare.parse_fasta(e["ref_fasta"]), dbcompare.parse_fasta(e["alt_fasta"])),
                            table, dbcompare.read_peptides(e["report_b"].encode()))
    present = []                                               # (modified sequence, charge, rt, truth key, version)
    truth = {}
    for i, s in enumerate(res.sites):
        key = (s["ref_protein"], s["position"])
        vers = {"alternative": ["alternative"], "reference": ["reference"], "both": ["alternative", "reference"]}.get(
            e["truth"][key], [])
        truth[key] = vers
        pr = next(p for p in dbcompare.pair_databases(dbcompare.parse_fasta(e["ref_fasta"]),
                                                      dbcompare.parse_fasta(e["alt_fasta"])).pairs if p.ref_id == key[0])
        for v in vers:
            lst = s["alt_peptides"] if v == "alternative" else s["ref_peptides"]
            if lst:
                seq0 = lst[0]["sequence"]
                (_, z0), _ = max(table.peptides[seq0]["precursors"].items(), key=lambda kv: float(kv[1].sum()))
                present.append((seq0, int(z0), float(rng.uniform(1.5, run_min - 1.5)), key, v))
    # a few unrelated background peptides so the windows are not empty
    n_cycles = int(run_min / cycle_min)
    edges = np.arange(300, 1800, 100)
    spectra = []
    scan = 0
    for c in range(n_cycles):
        t = c * cycle_min
        scan += 1
        mz1, in1 = list(rng.uniform(300, 1700, 40)), list(rng.uniform(1e3, 5e3, 40))
        for seq, z, rt, _, _ in present:
            g = np.exp(-((t - rt) ** 2) / (2 * 0.04 ** 2))
            sq, sh, nt = msms.parse_modified(seq)
            m = msms.precursor_mz(sq, msms.with_fixed_cam(sq, sh, False), nt, z)
            if g > 0.01:
                for k, f in enumerate((1.0, 0.55, 0.2)):
                    mz1.append(m + k * 1.00335 / z)
                    in1.append(5e6 * g * f)
        o = np.argsort(mz1)
        spectra.append({"scan": scan, "level": 1, "rt": t, "mz": np.array(mz1)[o], "int": np.array(in1)[o]})
        for w in edges:
            scan += 1
            mz2, in2 = list(rng.uniform(150, 1700, 120)), list(rng.uniform(5e2, 4e3, 120))
            for seq, z, rt, _, _ in present:
                sq, sh, nt = msms.parse_modified(seq)
                sh = msms.with_fixed_cam(sq, sh, False)
                m = msms.precursor_mz(sq, sh, nt, z)
                if w - 1 <= m <= w + 99 and abs(t - rt) < 0.2:
                    g = np.exp(-((t - rt) ** 2) / (2 * 0.04 ** 2))
                    for f in msms.fragments(sq, sh, nt):
                        if f.charge == 1 and rng.random() < 0.8:
                            mz2.append(f.mz * (1 + rng.normal(0, 2e-6)))
                            in2.append(2e5 * g * rng.uniform(0.2, 1.0))
            o = np.argsort(mz2)
            spectra.append({"scan": scan, "level": 2, "rt": t + 0.001, "mz": np.array(mz2)[o], "int": np.array(in2)[o],
                            "target": w + 49.0, "lo": 50.0, "hi": 50.0})

    def b64(a):
        return base64.b64encode(zlib.compress(np.asarray(a, dtype="<f8").tobytes())).decode()
    parts = ['<?xml version="1.0" encoding="utf-8"?>\n<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="%d">' % len(spectra)]
    for i, s_ in enumerate(spectra):
        cv = [f'<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{s_["level"]}"/>',
              f'<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{s_["rt"]:.5f}" unitName="minute"/>',
              f'<cvParam cvRef="MS" accession="MS:1003057" name="scan number" value="{s_["scan"]}"/>']
        if s_["level"] == 2:
            cv += [f'<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{s_["target"]}"/>',
                   f'<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{s_["lo"]}"/>',
                   f'<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{s_["hi"]}"/>']
        parts.append(f'<spectrum index="{i}" id="scan={s_["scan"]}" defaultArrayLength="{len(s_["mz"])}">' + "".join(cv)
                     + '<binaryDataArrayList count="2">'
                     + f'<binaryDataArray><cvParam accession="MS:1000523" name="64-bit float"/><cvParam accession="MS:1000574" name="zlib compression"/><cvParam accession="MS:1000514" name="m/z array"/><binary>{b64(s_["mz"])}</binary></binaryDataArray>'
                     + f'<binaryDataArray><cvParam accession="MS:1000523" name="64-bit float"/><cvParam accession="MS:1000574" name="zlib compression"/><cvParam accession="MS:1000515" name="intensity array"/><binary>{b64(s_["int"])}</binary></binaryDataArray>'
                     + "</binaryDataArrayList></spectrum>")
    parts.append("</spectrumList></run></mzML>")
    data = "".join(parts).encode()
    if str(path).endswith(".gz"):
        data = gzip.compress(data, 5)
    Path(path).write_bytes(data)
    return truth

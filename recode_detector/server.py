"""Local web server: serves the page and a small JSON API. Standard library only."""
from __future__ import annotations

import io
import json
import math
import os
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np
import pandas as pd

from . import __version__, dbcompare, exports, loaders, msms, pipeline, qc as qcmod, report, stats
from .example import example_bytes, simulate_db
from .history import RunStore
from .lifecycle import Lifecycle

WEB = Path(__file__).resolve().parent / "web"
DEFAULT_PORT = 8775
MAX_UPLOAD = 400 * 1024 * 1024
STORE = RunStore()
TYPES = {".html": "text/html; charset=utf-8", ".js": "text/javascript", ".css": "text/css",
         ".svg": "image/svg+xml", ".png": "image/png", ".ico": "image/x-icon"}


def clean(o):
    """Make numpy / pandas values JSON-safe (NaN and inf become null)."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, pd.Series):
        return clean(o.tolist())
    return o


class Session:
    """The one dataset the open window is working on."""

    def __init__(self):
        self.lock = threading.RLock()
        self.reset()

    def reset(self):
        self.filename = ""
        self.raw_bytes: bytes | None = None
        self.dataset: loaders.Dataset | None = None
        self.groups: dict[str, str] = {}
        self.proc: pipeline.Processed | None = None
        self.qc: dict | None = None
        self.result: pd.DataFrame | None = None
        self.result_meta: dict = {}
        self.is_example = False
        self.kinds: list[str] = []
        self.kind: str | None = None
        self.run_id: str | None = None
        self.qc_checks: list[dict] = []
        self.sample_qc: pd.DataFrame = pd.DataFrame()
        self.enrich_df: pd.DataFrame = pd.DataFrame()

    # ---- operations (each returns a JSON-ready dict)
    def state(self) -> dict:
        d = self.dataset
        out = {"version": __version__, "loaded": d is not None, "filename": self.filename,
               "is_example": self.is_example, "analysed": self.proc is not None}
        if d is not None:
            out.update({"format": d.source_format, "samples": d.samples, "groups": self.groups,
                        "notes": d.notes, "removed": d.removed, "n_proteins": int(len(d.matrix)),
                        "kinds": self._kinds(), "missing_pct": float(d.matrix.isna().values.mean() * 100)})
        if self.proc is not None:
            out.update({"qc": self.qc, "levels": sorted(set(self.proc.groups.values())), "log": self.proc.log,
                        "checks": self.qc_checks, "check_summary": qcmod.summarise(self.qc_checks),
                        "sample_qc": self.sample_qc.to_dict("records")})
        out["run_id"] = self.run_id
        return out

    def _kinds(self) -> list[str]:
        return self.kinds

    def reload_kind(self, kind: str) -> dict:
        with self.lock:
            self._need_data()
            return self.load(self.raw_bytes, self.filename, kind, self.is_example)

    def protein(self, pid: str) -> dict:
        with self.lock:
            self._need_analysis()
            if pid not in self.proc.log_imp.index:
                raise ValueError("Protein not found.")
            cols = list(self.proc.log_imp.columns)
            return {"id": pid, "samples": cols, "groups": [self.proc.groups[c] for c in cols],
                    "values": self.proc.log_imp.loc[pid].tolist(),
                    "imputed": [bool(x) for x in self.proc.imputed.loc[pid].tolist()]}

    def load(self, raw: bytes, name: str, kind: str | None = None, example: bool = False) -> dict:
        ds = loaders.load(raw, name, kind)
        with self.lock:
            self.reset()
            self.kind = kind
            self.filename, self.raw_bytes, self.dataset, self.is_example = name, raw, ds, example
            self.groups = loaders.guess_groups(ds.samples)
            self.kinds = []
            if ds.source_format == "maxquant":
                try:
                    self.kinds = loaders.maxquant_intensity_kinds(loaders.read_any(raw, name))
                except Exception:
                    pass
            return self.state()

    def set_groups(self, groups: dict) -> dict:
        with self.lock:
            self._need_data()
            self.groups = {s: str(groups.get(s, "")).strip() for s in self.dataset.samples}
            self.proc = self.qc = self.result = None
            return self.state()

    def analyse(self, params: dict) -> dict:
        with self.lock:
            self._need_data()
            p = pipeline.Params.from_dict(params)
            self.proc = pipeline.process(self.dataset.matrix, self.groups, p)
            self.qc = pipeline.qc_summary(self.dataset.matrix, self.proc)
            self.qc_checks, self.sample_qc = qcmod.quant_checks(self.dataset.matrix, self.proc, self.qc)
            self.result = None
            self.enrich_df = pd.DataFrame()
            levels = sorted(set(self.proc.groups.values()))
            self._autosave(params)
            return {"qc": self.qc, "log": self.proc.log, "groups": levels, "filtered_out": self.proc.filtered_out,
                    "checks": self.qc_checks, "check_summary": qcmod.summarise(self.qc_checks),
                    "sample_qc": self.sample_qc.to_dict("records"), "run_id": self.run_id}

    def _autosave(self, processing: dict):
        if self.is_example or self.raw_bytes is None:
            return
        name = Path(self.filename).stem + " (quantification)"
        meta = STORE.save("quant", name, {"data": (self.filename, self.raw_bytes)},
                          {"kind": self.kind, "groups": self.groups, "processing": processing, "comparison": None},
                          {"samples": len(self.proc.groups), "proteins_kept": int(len(self.proc.log_raw)),
                           "groups": sorted(set(self.proc.groups.values()))}, fp_keys=["kind", "groups", "processing"])
        self.run_id = meta["id"]

    def compare(self, a: str, b: str, method: str, use: str, fdr: float, lfc: float) -> dict:
        with self.lock:
            self._need_analysis()
            res = stats.compare(self.proc, a, b, method, use)
            res["gene"] = self.dataset.genes.reindex(res["protein"]).values
            res["call"] = stats.call_significant(res, fdr, lfc)
            self.result = res
            self.result_meta = {"a": a, "b": b, "method": method, "use": use, "fdr": fdr, "lfc": lfc,
                                "prior_df": res.attrs.get("prior_df")}
            if self.run_id:
                try:
                    STORE.update(self.run_id, params={"comparison": {"a": a, "b": b, "method": method, "use": use, "fdr": fdr, "lfc": lfc}})
                except ValueError:
                    pass
            order = res.sort_values("p_value", na_position="last")
            counts = res["call"].value_counts().to_dict()
            sig = order[order["call"] != "ns"]
            hm = pipeline.heatmap(self.proc, list(sig["protein"].head(100)))
            hm["labels"] = [str(self.dataset.genes.get(i, i)) for i in hm["ids"]]
            cols = ["protein", "gene", "log2fc", "mean_a", "mean_b", "t", "p_value", "adj_p",
                    "n_a_observed", "n_b_observed", "call"]
            return {"meta": self.result_meta, "counts": {k: int(counts.get(k, 0)) for k in ("up", "down", "ns")},
                    "rows": order[cols].to_dict("records"), "heatmap": hm, "n_tested": int(res["p_value"].notna().sum())}

    def enrich(self, gmt: str, direction: str) -> dict:
        with self.lock:
            if self.result is None:
                raise ValueError("Run a comparison first.")
            sets = stats.parse_gmt(gmt)
            if not sets:
                raise ValueError("No gene sets found. Expected GMT: name, description, then genes, tab-separated.")
            r = self.result
            hits = r[r["call"].isin(["up", "down"] if direction == "both" else [direction])]["gene"]
            out = stats.enrichment(list(hits), list(r["gene"]), sets)
            self.enrich_df = out
            return {"n_hits": int(len(hits)), "n_sets": len(sets), "rows": out.head(200).to_dict("records")}

    def open_run(self, meta: dict, files: dict) -> dict:
        name, raw = files["data"]
        prm = meta["params"]
        st = self.load(raw, name, prm.get("kind"))
        self.set_groups(prm.get("groups", {}))
        an = self.analyse(prm.get("processing", {}))
        cmp_ = None
        c = prm.get("comparison")
        if c:
            try:
                cmp_ = self.compare(c["a"], c["b"], c["method"], c["use"], c["fdr"], c["lfc"])
            except ValueError:
                cmp_ = None
        self.run_id = meta["id"]
        return {"mode": "quant", "state": self.state(), "analysis": an, "comparison": cmp_, "params": prm}

    def tables(self) -> dict[str, pd.DataFrame]:
        out: dict[str, pd.DataFrame] = {}
        if self.result is not None:
            out["differential_results"] = self.result.sort_values("p_value")
        if self.proc is not None:
            m = self.proc.log_imp.copy()
            m.insert(0, "gene", self.dataset.genes.reindex(m.index).values)
            m.index.name = "protein"
            out["processed_matrix"] = m.reset_index()
            out["qc_checks"] = pd.DataFrame(self.qc_checks)
            out["sample_qc"] = self.sample_qc
        if len(self.enrich_df):
            out["enrichment"] = self.enrich_df
        return out

    def settings_json(self) -> bytes:
        with self.lock:
            self._need_analysis()
            return json.dumps({"version": __version__, "file": self.filename, "groups": self.groups,
                               "parameters": self.proc.params.__dict__, "steps": self.proc.log,
                               "comparison": self.result_meta}, indent=2, default=str).encode()

    def _need_data(self):
        if self.dataset is None:
            raise ValueError("Load a protein table first.")

    def _need_analysis(self):
        self._need_data()
        if self.proc is None:
            raise ValueError("Run the processing step first.")


class DbSession:
    """Inputs and result of the database-check mode."""

    SLOTS = ("report", "report_b", "fasta_ref", "fasta_alt")

    def __init__(self):
        self.lock = threading.RLock()
        self.job: msms.SpectraJob | None = None
        self.reset()

    def reset(self):
        self.files: dict[str, dict] = {}            # slot -> {"name", "raw", "table"/"fasta", "detail"}
        self.result: dbcompare.Result | None = None
        self.pairing = None
        self.marker = ""
        self.run_id: str | None = None
        self.is_example = False
        self.qc: dict = {}
        self.spec: dict | None = None

    # ---- state
    def state(self) -> dict:
        info = {k: {"name": v["name"], "detail": v["detail"]} for k, v in self.files.items()}
        out = {"files": info, "ready": "report" in self.files and "fasta_ref" in self.files,
               "is_example": self.is_example, "has_result": self.result is not None, "run_id": self.run_id,
               "raw_reader": bool(msms.find_thermo_parser()), "job": self.job.public() if self.job else None}
        if self.result is not None:
            out["result"] = self._payload()
        return out

    def _payload(self) -> dict:
        r = self.result
        return {"samples": r.samples, "sites": r.sites, "summary": r.summary, "scatter": r.scatter,
                "coverage": r.coverage, "qc": self.qc, "spectra": self.spec}

    # ---- inputs
    def load(self, slot: str, raw: bytes, name: str) -> dict:
        if slot not in self.SLOTS:
            raise ValueError("Unknown slot.")
        with self.lock:
            if slot in ("report", "report_b"):
                t = dbcompare.read_peptides(raw, name)
                kind = "long-format report" if t.has_rt else "precursor matrix"
                self.files[slot] = {"name": name, "raw": raw, "table": t,
                                    "detail": f"{len(t.peptides):,} peptides, {len(t.samples)} run(s), {kind}"}
            else:
                fa = dbcompare.parse_fasta(raw.decode("utf-8-sig", errors="replace"))
                if not fa:
                    raise ValueError("No FASTA entries found (lines starting with “>”).")
                self.files[slot] = {"name": name, "raw": raw, "fasta": fa, "detail": f"{len(fa):,} proteins"}
            self.result, self.run_id, self.spec, self.qc = None, None, None, {}
            return self.state()

    def drop(self, slot: str) -> dict:
        with self.lock:
            self.files.pop(slot, None)
            self.result, self.run_id, self.spec, self.qc = None, None, None, {}
            return self.state()

    def example(self) -> dict:
        e = simulate_db()
        with self.lock:
            self.reset()
            self.is_example = True
        for slot, key, name in (("report", "report", "simulated_report.pr_matrix.tsv"),
                                ("report_b", "report_b", "simulated_standard_only.pr_matrix.tsv"),
                                ("fasta_ref", "ref_fasta", "simulated_standard.fasta"),
                                ("fasta_alt", "alt_fasta", "simulated_variants.fasta")):
            self.load(slot, e[key].encode(), name)
        return self.state()

    # ---- analysis
    def run(self, marker: str) -> dict:
        with self.lock:
            if "report" not in self.files or "fasta_ref" not in self.files:
                raise ValueError("Load the peptide report and the protein FASTA first.")
            ref = self.files["fasta_ref"]["fasta"]
            alt = self.files["fasta_alt"]["fasta"] if "fasta_alt" in self.files else None
            pairing = dbcompare.pair_databases(ref, alt, marker.strip())
            if not pairing.pairs:
                raise ValueError("No protein differs between the two databases"
                                 + (f" ({len(pairing.unpaired)} alternative proteins could not be paired: "
                                    "they have no reference protein of the same length)." if pairing.unpaired else "."))
            b = self.files["report_b"]["table"] if "report_b" in self.files else None
            table = self.files["report"]["table"]
            self.result = dbcompare.analyse(pairing, table, b)
            self.pairing, self.marker = pairing, marker.strip()
            checks, runs = qcmod.db_checks(table, self.result.coverage)
            self.qc = {"checks": checks, "runs": runs.to_dict("records"), "summary": qcmod.summarise(checks)}
            self.spec = None
            if not self.is_example:
                st = self.result.summary["status"]
                meta = STORE.save("db", Path(self.files["report"]["name"]).stem + " (database check)",
                                  {slot: (v["name"], v["raw"]) for slot, v in self.files.items()}, {"marker": self.marker},
                                  {"sites": self.result.summary["n_sites"], "alternative": st["alternative"],
                                   "reference": st["reference"], "both": st["both"], "none": st["none"],
                                   "peptides": self.result.summary["n_peptides"]}, fp_keys=["marker"])
                self.run_id = meta["id"]
                self.spec = STORE.load_extra(self.run_id, "spectra")
            return self._payload()

    def open_run(self, meta: dict, files: dict) -> dict:
        with self.lock:
            self.reset()
            for slot, (name, raw) in files.items():
                self.load(slot, raw, name)
            self.run(meta["params"].get("marker", ""))
            self.run_id = meta["id"]
            self.spec = STORE.load_extra(self.run_id, "spectra")
        return {"mode": "db", "state": self.state(), "params": meta["params"]}

    # ---- spectrum evidence
    def _best_precursor(self, seq: str):
        table = self.files["report"]["table"]
        p = table.peptides.get(seq)
        if not p or not p["precursors"]:
            return seq, 2
        (modseq, z), _ = max(p["precursors"].items(), key=lambda kv: float(np.sum(kv[1])))
        return modseq, int(z) if str(z).isdigit() else 2

    def default_targets(self, site_indexes: list[int] | None = None) -> list[dict]:
        out, seen = [], set()
        sites = self.result.sites
        for i, s in enumerate(sites):
            if site_indexes is not None and i not in site_indexes:
                continue
            if site_indexes is None and s["status"] == "none":
                continue
            pr = next(p for p in self.pairing.pairs if p.ref_id == s["ref_protein"] and p.alt_id == s["alt_protein"])
            q = s["position"] - 1
            for version, key in (("alternative", "alt_peptides"), ("reference", "ref_peptides")):
                cand = sorted((p for p in s[key] if not p["also_elsewhere"]), key=lambda p: -p["total_intensity"])
                if not cand:
                    continue
                pep = cand[0]
                modseq, z = self._best_precursor(pep["sequence"])
                idx = q - (pep["start"] - 1)
                other_seq = (pr.ref_seq if version == "alternative" else pr.alt_seq)[pep["start"] - 1:pep["end"]]
                pair = f"{i}|{pep['start']}-{pep['end']}"
                for ver, ms, role in ((version, modseq, "observed"),
                                      ("reference" if version == "alternative" else "alternative",
                                       msms.swap_residue(modseq, idx, other_seq[idx]), "counterpart")):
                    k = f"{pair}|{ver}"
                    if k in seen:
                        continue
                    seen.add(k)
                    rt_hint = None
                    if role == "observed" and self.files["report"]["table"].has_rt:
                        rt_hint = self.files["report"]["table"].peptides[pep["sequence"]].get("rt")   # DIA-NN's own retention time
                    out.append({"key": k, "site": i, "version": ver, "role": role, "modified": ms, "charge": z,
                                "changed": [idx], "pair": pair, "rt": rt_hint})
        return out[:60]

    def spectra_start(self, path: str, ppm: float, half_window: float, sites: list[int] | None = None) -> dict:
        with self.lock:
            if self.result is None:
                raise ValueError("Run the comparison first.")
            if self.job is not None and self.job.status in ("queued", "running"):
                raise ValueError("A spectrum check is already running.")
            path = str(path).strip().strip('"').strip("'")
            if not path or not Path(path).is_file():
                raise ValueError(f"File not found: {path or '(no path given)'}")
            if not path.lower().endswith(msms.RAW_EXT + msms.MZML_EXT):
                raise ValueError("Use a Thermo .raw file or an mzML / mzML.gz file.")
            targets = self.default_targets(sites)
            if not targets:
                raise ValueError("No site has an identified peptide to check yet.")
            job = msms.SpectraJob(id=str(int(__import__("time").time() * 1000)))
            self.job = job
        threading.Thread(target=self._run_job, args=(job, path, targets, ppm, half_window), daemon=True).start()
        return {"job": job.public(), "n_targets": len(targets)}

    def _run_job(self, job, path, targets, ppm, hw):
        msms.run_job(job, path, targets, ppm, hw)
        if job.status == "done":
            with self.lock:
                self.spec = {"path": path, "params": {"ppm": ppm, "half_window": hw}, "source": job.source_info,
                             "results": job.results, "custom": (self.spec or {}).get("custom", [])}
                if self.run_id:
                    try:
                        STORE.save_extra(self.run_id, "spectra", self.spec)
                        STORE.update(self.run_id, params={"spectra_path": path})
                    except ValueError:
                        pass

    def spectra_peptide(self, path: str, sequence: str, charge: int, ppm: float, half_window: float) -> dict:
        with self.lock:
            if self.job is not None and self.job.status in ("queued", "running"):
                raise ValueError("A spectrum check is already running.")
            path = str(path).strip().strip('"').strip("'")
            if not path or not Path(path).is_file():
                raise ValueError(f"File not found: {path or '(no path given)'}")
            seq = sequence.strip().replace(" ", "")
            if not seq:
                raise ValueError("Enter a peptide sequence.")
            msms.parse_modified(seq)                      # fails early on unknown letters
            if not 1 <= int(charge) <= 6:
                raise ValueError("Charge must be between 1 and 6.")
            job = msms.SpectraJob(id=str(int(__import__("time").time() * 1000)))
            self.job = job
        target = {"key": "custom|" + seq, "site": None, "version": "custom", "modified": seq, "charge": int(charge), "changed": []}

        def go():
            msms.run_job(job, path, [target], ppm, half_window)
            if job.status == "done" and job.results:
                with self.lock:
                    self.spec = self.spec or {"path": path, "params": {"ppm": ppm, "half_window": half_window},
                                              "source": job.source_info, "results": [], "custom": []}
                    self.spec["custom"] = ([job.results[0]] + self.spec.get("custom", []))[:12]
                    if self.run_id:
                        try:
                            STORE.save_extra(self.run_id, "spectra", self.spec)
                        except ValueError:
                            pass
        threading.Thread(target=go, daemon=True).start()
        return {"job": job.public()}

    def spectra_status(self) -> dict:
        with self.lock:
            return {"job": self.job.public() if self.job else None, "spectra": self.spec}

    def spectra_cancel(self) -> dict:
        with self.lock:
            if self.job is not None:
                self.job.cancel.set()
            msms.kill_children()
            return {"ok": True}

    # ---- tables and downloads
    def tables(self) -> dict[str, pd.DataFrame]:
        if self.result is None:
            return {}
        out = {"sites": dbcompare.sites_frame(self.result), "peptides": dbcompare.peptides_frame(self.result),
               "coverage": dbcompare.coverage_frame(self.result)}
        if self.qc:
            out["qc_checks"] = pd.DataFrame(self.qc["checks"])
            out["qc_runs"] = pd.DataFrame(self.qc["runs"])
        if self.spec and self.spec.get("results"):
            out["spectra"], out["spectra_ions"] = msms.spectra_frames(self.spec["results"], self.result.sites)
        return out

    def pdf(self) -> bytes:
        with self.lock:
            if self.result is None:
                raise ValueError("Run the comparison first.")
            names = {k: v["name"] for k, v in self.files.items()}
            kw = {"note": "SIMULATED EXAMPLE: invented proteins and peptides, not real results."} if self.is_example else {}
            return report.build_pdf(self.result, names, qc=self.qc, spectra=self.spec, **kw)


SESSION = Session()
DB = DbSession()
LIFE: Lifecycle | None = None

DB_TABLES = {"db_sites": "sites", "db_peptides": "peptides", "db_coverage": "coverage", "db_qc_checks": "qc_checks",
             "db_qc_runs": "qc_runs", "db_spectra": "spectra", "db_spectra_ions": "spectra_ions"}
Q_TABLES = {"results": "differential_results", "processed": "processed_matrix", "qc_checks": "qc_checks",
            "sample_qc": "sample_qc", "enrichment": "enrichment"}


def runs_payload() -> dict:
    runs = []
    for m in STORE.list():
        runs.append({k: m.get(k) for k in ("id", "kind", "name", "created", "updated", "note", "summary", "disk_bytes")}
                    | {"files": {k: v.get("name") for k, v in m.get("files", {}).items()},
                       "has_spectra": bool(m.get("params", {}).get("spectra_path"))})
    return {"runs": runs, "root": str(STORE.root)}


def open_run(run_id: str) -> dict:
    meta = STORE.get(run_id)
    files = STORE.load_files(run_id)
    return SESSION.open_run(meta, files) if meta["kind"] == "quant" else DB.open_run(meta, files)


def pick_file() -> str:
    if sys.platform != "darwin":
        raise ValueError("The file picker is only available on macOS. Type or paste the path instead.")
    import subprocess
    r = subprocess.run(["osascript", "-e", 'POSIX path of (choose file with prompt "Choose a Thermo .raw or mzML file")'],
                       capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        raise ValueError("No file was chosen.")
    return r.stdout.strip()


def build_download(which: str, fmt: str) -> tuple[str, bytes, str]:
    fmt = (fmt or "csv").lower()
    if which == "db_report":
        return "recode_detector_report.pdf", DB.pdf(), "application/pdf"
    if which == "settings":
        return "analysis_settings.json", SESSION.settings_json(), "application/json"
    if which == "runs":
        df = pd.DataFrame([{"name": r["name"], "kind": r["kind"], "created": pd.to_datetime(r["created"], unit="s"),
                            "updated": pd.to_datetime(r["updated"], unit="s"), "note": r["note"],
                            "files": ", ".join(r["files"].values())} for r in runs_payload()["runs"]])
        return f"saved_runs.{fmt}", exports.frame_bytes(df, fmt, "Runs"), exports.FORMATS[fmt]
    if which in ("zip_db", "zip_quant"):
        if which == "zip_db":
            tables, extras = DB.tables(), {}
            if not tables:
                raise ValueError("Run the comparison first.")
            try:
                extras["recode_detector_report.pdf"] = DB.pdf()
            except ValueError:
                pass
            return "recode_detector_database_check.zip", exports.zip_bytes(tables, extras, fmt), "application/zip"
        tables = SESSION.tables()
        if not tables:
            raise ValueError("Run the processing step first.")
        extras = {}
        try:
            extras["analysis_settings.json"] = SESSION.settings_json()
        except ValueError:
            pass
        return "recode_detector_quantification.zip", exports.zip_bytes(tables, extras, fmt), "application/zip"
    if which in DB_TABLES:
        df = DB.tables().get(DB_TABLES[which])
        label = DB_TABLES[which]
    elif which in Q_TABLES:
        df = SESSION.tables().get(Q_TABLES[which])
        label = Q_TABLES[which]
    else:
        raise ValueError("Unknown download.")
    if df is None:
        raise ValueError("That table is not available yet. Run the analysis first.")
    return f"{label}.{fmt}", exports.frame_bytes(df, fmt, label), exports.FORMATS[fmt]


class Handler(BaseHTTPRequestHandler):
    server_version = "RecodeDetector"
    protocol_version = "HTTP/1.1"

    def log_message(self, *a):
        pass

    # ---- helpers
    def _send(self, code, body: bytes, ctype="application/json", extra=None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(clean(obj), allow_nan=False).encode())

    def _body(self) -> bytes:
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_UPLOAD:
            raise ValueError("File is larger than 400 MB.")
        return self.rfile.read(n) if n else b""

    def _guard(self, fn):
        try:
            self._json(fn())
        except ValueError as e:
            self._json({"error": str(e)}, 400)
        except Exception as e:                           # keep the server alive, show the reason
            import traceback
            traceback.print_exc()
            self._json({"error": f"{type(e).__name__}: {e}"}, 500)

    # ---- routes
    def do_GET(self):
        u = urlparse(self.path)
        if u.path == "/api/state":
            return self._json(SESSION.state())
        if u.path == "/api/window":
            return self._window()
        if u.path == "/api/db/state":
            return self._json(DB.state())
        if u.path == "/api/protein":
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            return self._guard(lambda: SESSION.protein(q.get("id", "")))
        if u.path == "/api/runs":
            return self._guard(runs_payload)
        if u.path == "/api/db/spectra/status":
            return self._guard(DB.spectra_status)
        if u.path.startswith("/api/download/"):
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            try:
                name, data, ctype = build_download(u.path.rsplit("/", 1)[1], q.get("fmt", "csv"))
            except ValueError as e:
                return self._json({"error": str(e)}, 400)
            except Exception as e:
                return self._json({"error": f"{type(e).__name__}: {e}"}, 500)
            return self._send(200, data, ctype, {"Content-Disposition": f'attachment; filename="{name}"'})
        rel = "index.html" if u.path in ("/", "") else u.path.lstrip("/")
        f = (WEB / rel).resolve()
        if WEB.resolve() not in f.parents and f != WEB.resolve() or not f.is_file():
            return self._send(404, b"not found", "text/plain")
        self._send(200, f.read_bytes(), TYPES.get(f.suffix, "application/octet-stream"))

    def do_POST(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            body = self._body()
        except ValueError as e:
            return self._json({"error": str(e)}, 413)
        js = lambda: json.loads(body or b"{}")             # noqa: E731
        if u.path == "/api/load":
            return self._guard(lambda: SESSION.load(body, q.get("name", "table.txt"), q.get("kind")))
        if u.path == "/api/example":
            return self._guard(lambda: SESSION.load(example_bytes(), "simulated_proteinGroups.txt",
                                                    None, example=True))
        if u.path == "/api/db/load":
            return self._guard(lambda: DB.load(q.get("slot", ""), body, q.get("name", "file")))
        if u.path == "/api/db/drop":
            return self._guard(lambda: DB.drop(q.get("slot", "")))
        if u.path == "/api/db/example":
            return self._guard(DB.example)
        if u.path == "/api/db/run":
            return self._guard(lambda: DB.run(str(js().get("marker", ""))))
        if u.path == "/api/kind":
            return self._guard(lambda: SESSION.reload_kind(q.get("kind", "")))
        if u.path == "/api/groups":
            return self._guard(lambda: SESSION.set_groups(js().get("groups", {})))
        if u.path == "/api/analyze":
            return self._guard(lambda: SESSION.analyse(js().get("params", {})))
        if u.path == "/api/compare":
            d = js()
            return self._guard(lambda: SESSION.compare(
                d["a"], d["b"], d.get("method", "moderated"), d.get("use", "imputed"),
                float(d.get("fdr", 0.05)), float(d.get("lfc", 1.0))))
        if u.path == "/api/enrich":
            d = js()
            return self._guard(lambda: SESSION.enrich(d.get("gmt", ""), d.get("direction", "both")))
        if u.path == "/api/runs/open":
            return self._guard(lambda: open_run(q.get("id", "")))
        if u.path == "/api/runs/rename":
            d = js()
            return self._guard(lambda: (STORE.rename(d.get("id", ""), d.get("name", "")), runs_payload())[1])
        if u.path == "/api/runs/note":
            d = js()
            return self._guard(lambda: (STORE.set_note(d.get("id", ""), d.get("note", "")), runs_payload())[1])
        if u.path == "/api/runs/delete":
            return self._guard(lambda: (STORE.delete(q.get("id", "")), runs_payload())[1])
        if u.path == "/api/pick":
            return self._guard(lambda: {"path": pick_file()})
        if u.path == "/api/db/spectra/start":
            d = js()
            return self._guard(lambda: DB.spectra_start(d.get("path", ""), float(d.get("ppm", 10)),
                                                        float(d.get("half_window", 0.15)), d.get("sites")))
        if u.path == "/api/db/spectra/peptide":
            d = js()
            return self._guard(lambda: DB.spectra_peptide(d.get("path", ""), d.get("sequence", ""), int(d.get("charge", 2)),
                                                          float(d.get("ppm", 10)), float(d.get("half_window", 0.15))))
        if u.path == "/api/db/spectra/cancel":
            return self._guard(DB.spectra_cancel)
        if u.path == "/api/quit":
            self._json({"ok": True})
            if LIFE:
                LIFE.quit_now()
            return
        self._send(404, b"not found", "text/plain")

    def _window(self):
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Connection", "close")
        self.end_headers()
        if LIFE:
            LIFE.window_opened()
        try:
            while True:
                self.wfile.write(b": hold\n\n")
                self.wfile.flush()
                threading.Event().wait(2.0)
        except (BrokenPipeError, ConnectionResetError, OSError):
            pass
        finally:
            if LIFE:
                LIFE.window_closed()
        self.close_connection = True


def serve(port: int = DEFAULT_PORT, open_browser: bool = True) -> None:
    global LIFE
    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    httpd.daemon_threads = True

    def quit_():
        sys.stdout.flush()
        try:
            msms.kill_children()                 # never leave a converter running after the program quits
        except Exception:
            pass
        os._exit(0)

    LIFE = Lifecycle(on_quit=quit_)
    LIFE.watch()
    url = f"http://localhost:{port}"
    print(f"Recode Detector {__version__} at {url}", flush=True)
    if open_browser:
        import webbrowser
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(prog="recode-detector", description="Recode Detector local app")
    ap.add_argument("--port", type=int, default=int(os.environ.get("RD_PORT", DEFAULT_PORT)))
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()
    serve(a.port, not a.no_browser)


if __name__ == "__main__":
    main()

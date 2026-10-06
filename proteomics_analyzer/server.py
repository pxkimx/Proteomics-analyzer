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

from . import __version__, loaders, pipeline, stats
from .example import example_bytes
from .lifecycle import Lifecycle

WEB = Path(__file__).resolve().parent / "web"
DEFAULT_PORT = 8775
MAX_UPLOAD = 400 * 1024 * 1024
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
            out.update({"qc": self.qc, "levels": sorted(set(self.proc.groups.values())), "log": self.proc.log})
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
            self.result = None
            levels = sorted(set(self.proc.groups.values()))
            return {"qc": self.qc, "log": self.proc.log, "groups": levels,
                    "filtered_out": self.proc.filtered_out}

    def compare(self, a: str, b: str, method: str, use: str, fdr: float, lfc: float) -> dict:
        with self.lock:
            self._need_analysis()
            res = stats.compare(self.proc, a, b, method, use)
            res["gene"] = self.dataset.genes.reindex(res["protein"]).values
            res["call"] = stats.call_significant(res, fdr, lfc)
            self.result = res
            self.result_meta = {"a": a, "b": b, "method": method, "use": use, "fdr": fdr, "lfc": lfc,
                                "prior_df": res.attrs.get("prior_df")}
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
            return {"n_hits": int(len(hits)), "n_sets": len(sets), "rows": out.head(200).to_dict("records")}

    def csv(self, which: str) -> tuple[str, bytes]:
        with self.lock:
            buf = io.StringIO()
            if which == "results":
                if self.result is None:
                    raise ValueError("Run a comparison first.")
                m = self.result_meta
                self.result.sort_values("p_value").to_csv(buf, index=False)
                return f"differential_{m['a']}_vs_{m['b']}.csv", buf.getvalue().encode()
            if which == "processed":
                self._need_analysis()
                out = self.proc.log_imp.copy()
                out.insert(0, "gene", self.dataset.genes.reindex(out.index).values)
                out.index.name = "protein"
                out.to_csv(buf)
                return "processed_log2_intensities.csv", buf.getvalue().encode()
            if which == "settings":
                self._need_analysis()
                return "analysis_settings.json", json.dumps(
                    {"version": __version__, "file": self.filename, "groups": self.groups,
                     "parameters": self.proc.params.__dict__, "steps": self.proc.log,
                     "comparison": self.result_meta}, indent=2, default=str).encode()
        raise ValueError("unknown download")

    def _need_data(self):
        if self.dataset is None:
            raise ValueError("Load a protein table first.")

    def _need_analysis(self):
        self._need_data()
        if self.proc is None:
            raise ValueError("Run the processing step first.")


SESSION = Session()
LIFE: Lifecycle | None = None


class Handler(BaseHTTPRequestHandler):
    server_version = "ProteomicsAnalyzer"
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
        if u.path == "/api/protein":
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            return self._guard(lambda: SESSION.protein(q.get("id", "")))
        if u.path.startswith("/api/download/"):
            try:
                name, data = SESSION.csv(u.path.rsplit("/", 1)[1])
            except ValueError as e:
                return self._json({"error": str(e)}, 400)
            ctype = "application/json" if name.endswith(".json") else "text/csv"
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
        os._exit(0)

    LIFE = Lifecycle(on_quit=quit_)
    LIFE.watch()
    url = f"http://localhost:{port}"
    print(f"Proteomics Analyzer {__version__} at {url}", flush=True)
    if open_browser:
        import webbrowser
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(prog="proteomics-analyzer", description="Proteomics Analyzer local app")
    ap.add_argument("--port", type=int, default=int(os.environ.get("PA_PORT", DEFAULT_PORT)))
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()
    serve(a.port, not a.no_browser)


if __name__ == "__main__":
    main()

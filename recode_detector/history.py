"""Saved runs: every analysis is stored on disk so it can be reopened later.

A run = the input files (gzipped), the settings, a short summary and a note. Results are not stored: they are
recomputed from the inputs when a run is opened, which is deterministic and keeps the store small.
Layout: <root>/<id>/meta.json and <root>/<id>/files/<slot>.gz  (plus optional <root>/<id>/extra.json).
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import secrets
import shutil
import sys
import time
from pathlib import Path


def default_root() -> Path:
    env = os.environ.get("RD_HOME")
    if env:
        return Path(env) / "runs"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "RecodeDetector" / "runs"
    return Path.home() / ".recode_detector" / "runs"


class RunStore:
    def __init__(self, root: Path | None = None):
        self.root = Path(root) if root else default_root()

    # ---- helpers
    def _dir(self, run_id: str) -> Path:
        if not run_id or "/" in run_id or ".." in run_id:
            raise ValueError("Bad run id.")
        return self.root / run_id

    def _read(self, run_id: str) -> dict:
        p = self._dir(run_id) / "meta.json"
        if not p.exists():
            raise ValueError("That saved run no longer exists.")
        return json.loads(p.read_text())

    def _write(self, run_id: str, meta: dict) -> None:
        d = self._dir(run_id)
        d.mkdir(parents=True, exist_ok=True)
        tmp = d / "meta.json.tmp"
        tmp.write_text(json.dumps(meta, indent=1, default=str))
        tmp.replace(d / "meta.json")                       # atomic: a crash never leaves a half-written file

    @staticmethod
    def fingerprint(kind: str, fhash: dict, params: dict) -> str:
        h = hashlib.sha256(kind.encode())
        for slot in sorted(fhash):
            h.update(slot.encode() + b"\0" + fhash[slot].encode())
        h.update(json.dumps(params, sort_keys=True, default=str).encode())
        return h.hexdigest()[:20]

    # ---- operations
    def save(self, kind: str, name: str, files: dict[str, tuple[str, bytes]], params: dict, summary: dict,
             fp_keys: list[str] | None = None) -> dict:
        """Create a run, or refresh the existing one with identical inputs and settings.
        Only `fp_keys` of the params define a run's identity; later additions (a comparison, a spectra path) do not."""
        fhash = {slot: hashlib.sha256(data).hexdigest()[:20] for slot, (_, data) in files.items()}
        keys = fp_keys if fp_keys is not None else sorted(params)
        fp = self.fingerprint(kind, fhash, {k: params.get(k) for k in keys})
        for m in self.list():
            if m.get("fingerprint") == fp:
                m["updated"] = time.time()
                m["summary"] = summary
                self._write(m["id"], m)
                return m
        run_id = time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(2)
        d = self._dir(run_id)
        (d / "files").mkdir(parents=True, exist_ok=True)
        meta_files = {}
        for slot, (fname, data) in files.items():
            with gzip.open(d / "files" / f"{slot}.gz", "wb", compresslevel=6) as fh:
                fh.write(data)
            meta_files[slot] = {"name": fname, "size": len(data)}
        meta = {"id": run_id, "kind": kind, "name": name, "created": time.time(), "updated": time.time(), "note": "",
                "params": params, "summary": summary, "files": meta_files, "fhash": fhash, "fp_keys": keys, "fingerprint": fp}
        self._write(run_id, meta)
        return meta

    def list(self) -> list[dict]:
        out = []
        if not self.root.exists():
            return out
        for d in self.root.iterdir():
            try:
                if (d / "meta.json").exists():
                    m = json.loads((d / "meta.json").read_text())
                    m["disk_bytes"] = sum(f.stat().st_size for f in d.rglob("*") if f.is_file())
                    out.append(m)
            except (OSError, ValueError):
                continue                                    # a damaged folder never hides the others
        return sorted(out, key=lambda m: m.get("updated", 0), reverse=True)

    def get(self, run_id: str) -> dict:
        return self._read(run_id)

    def load_files(self, run_id: str) -> dict[str, tuple[str, bytes]]:
        meta, out = self._read(run_id), {}
        for slot, info in meta["files"].items():
            with gzip.open(self._dir(run_id) / "files" / f"{slot}.gz", "rb") as fh:
                out[slot] = (info["name"], fh.read())
        return out

    def update(self, run_id: str, **fields) -> dict:
        meta = self._read(run_id)
        for k, v in fields.items():
            meta[k] = {**meta.get("params", {}), **v} if k == "params" else v
        meta["updated"] = time.time()
        if "params" in fields:
            meta["fingerprint"] = self.fingerprint(meta["kind"], meta.get("fhash", {}),
                                                   {k: meta["params"].get(k) for k in meta.get("fp_keys", sorted(meta["params"]))})
        self._write(run_id, meta)
        return meta

    def rename(self, run_id: str, name: str) -> dict:
        name = name.strip()
        if not name:
            raise ValueError("A run needs a name.")
        meta = self._read(run_id)
        meta["name"] = name[:120]
        self._write(run_id, meta)
        return meta

    def set_note(self, run_id: str, note: str) -> dict:
        meta = self._read(run_id)
        meta["note"] = note[:5000]
        self._write(run_id, meta)
        return meta

    def save_extra(self, run_id: str, key: str, data: dict) -> None:
        d = self._dir(run_id)
        p = d / "extra.json"
        extra = json.loads(p.read_text()) if p.exists() else {}
        extra[key] = data
        p.write_text(json.dumps(extra, default=str))

    def load_extra(self, run_id: str, key: str):
        p = self._dir(run_id) / "extra.json"
        return json.loads(p.read_text()).get(key) if p.exists() else None

    def delete(self, run_id: str) -> None:
        d = self._dir(run_id)
        if not d.exists():
            raise ValueError("That saved run no longer exists.")
        shutil.rmtree(d)

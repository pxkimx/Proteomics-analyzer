"""Spectrum evidence: look at the raw data behind a peptide call.

For a peptide (and its counterpart with the other residue at a changed position) this finds the precursor in the MS1
scans, the fragmentation (MS/MS) scans that isolate it around its elution peak, matches theoretical b and y ions,
and measures how well the fragments rise and fall together. Two sources are supported: Thermo .raw files (read in place
with ThermoRawFileParser's query and xic modes, so a multi-gigabyte file is never converted) and mzML (small files).
"""
from __future__ import annotations

import base64
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import threading
import zlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

PROTON = 1.007276
H2O = 18.010565
MONO = {"A": 71.037114, "R": 156.101111, "N": 114.042927, "D": 115.026943, "C": 103.009185, "E": 129.042593,
        "Q": 128.058578, "G": 57.021464, "H": 137.058912, "I": 113.084064, "L": 113.084064, "K": 128.094963,
        "M": 131.040485, "F": 147.068414, "P": 97.052764, "S": 87.032028, "T": 101.047679, "W": 186.079313,
        "Y": 163.063329, "V": 99.068414}
UNIMOD = {4: 57.021464, 35: 15.994915, 1: 42.010565, 21: 79.966331, 7: 0.984016, 5: 43.005814, 28: -17.026549,
          27: -18.010565, 34: 14.01565, 36: 28.0313, 37: 42.04695, 121: 114.042927, 188: 6.020129}
RAW_EXT, MZML_EXT = (".raw",), (".mzml", ".mzml.gz")


# ---------------------------------------------------------------- chemistry

def parse_modified(modseq: str) -> tuple[str, list[float], float]:
    """'AAC(UniMod:4)DEK' -> ('AACDEK', per-residue mass shifts, N-terminal shift)."""
    seq, shifts, nterm = [], [], 0.0
    i = 0
    while i < len(modseq):
        ch = modseq[i]
        if ch == "(":
            j = modseq.find(")", i)
            tag = modseq[i + 1:j] if j != -1 else ""
            m = re.search(r"(\d+)", tag)
            delta = UNIMOD.get(int(m.group(1)), 0.0) if m and "unimod" in tag.lower() else 0.0
            if not seq:
                nterm += delta
            else:
                shifts[-1] += delta
            i = (j if j != -1 else i) + 1
            continue
        if ch.isalpha() and ch.isupper():
            seq.append(ch)
            shifts.append(0.0)
        i += 1
    return "".join(seq), shifts, nterm


def with_fixed_cam(seq: str, shifts: list[float], explicit: bool) -> list[float]:
    """Carbamidomethyl-Cys is the usual fixed modification; add it when no modification string was given."""
    if explicit:
        return shifts
    return [s + (57.021464 if aa == "C" else 0.0) for aa, s in zip(seq, shifts)]


def residue_masses(seq: str, shifts: list[float]) -> list[float]:
    try:
        return [MONO[a] + shifts[i] for i, a in enumerate(seq)]
    except KeyError as e:
        raise ValueError(f"Unknown amino acid {e} in {seq}.")


def precursor_mz(seq: str, shifts: list[float], nterm: float, z: int) -> float:
    return (sum(residue_masses(seq, shifts)) + nterm + H2O + z * PROTON) / z


@dataclass
class Fragment:
    label: str
    kind: str          # "b" or "y"
    index: int         # number of residues in the fragment
    charge: int
    mz: float
    covers: tuple      # residue positions (0-based) inside the fragment


def fragments(seq: str, shifts: list[float], nterm: float = 0.0, max_charge: int = 2) -> list[Fragment]:
    m = residue_masses(seq, shifts)
    n = len(seq)
    out = []
    cum = nterm
    for i in range(1, n):
        cum += m[i - 1]
        for z in range(1, max_charge + 1):
            out.append(Fragment(f"b{i}" + "+" * z if z > 1 else f"b{i}", "b", i, z, (cum + z * PROTON) / z, tuple(range(i))))
    cum = H2O
    for j in range(1, n):
        cum += m[n - j]
        for z in range(1, max_charge + 1):
            out.append(Fragment(f"y{j}" + "+" * z if z > 1 else f"y{j}", "y", j, z, (cum + z * PROTON) / z, tuple(range(n - j, n))))
    return out


def swap_residue(modseq: str, idx: int, new: str) -> str:
    """Replace the idx-th residue (0-based) of a modified sequence, keeping every modification tag."""
    out, n, depth = [], 0, 0
    for ch in modseq:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            out.append(ch)
            continue
        if depth == 0 and ch.isalpha() and ch.isupper():
            out.append(new if n == idx else ch)
            n += 1
        else:
            out.append(ch)
    return "".join(out)


# ---------------------------------------------------------------- spectra and sources

@dataclass
class Spectrum:
    scan: int
    rt: float                       # minutes
    level: int
    target: float | None
    lo: float
    hi: float
    mz: np.ndarray
    inten: np.ndarray


def _spectrum_from_attrs(attrs: dict, mz, inten) -> Spectrum:
    g = lambda k, d=None: attrs.get(k, d)             # noqa: E731
    f = lambda v, d=0.0: float(v) if v not in (None, "") else d   # noqa: E731
    tgt = g("isolation window target m/z", g("selected ion m/z"))
    return Spectrum(int(float(g("scan number", 0))), f(g("scan start time")), int(float(g("ms level", 1))),
                    f(tgt) if tgt is not None else None, f(g("isolation window lower offset"), 1.0),
                    f(g("isolation window upper offset"), 1.0), np.asarray(mz, float), np.asarray(inten, float))


_ACTIVE: set = set()
_ACTIVE_LOCK = threading.Lock()


def _kill(proc) -> None:
    if proc.poll() is None:
        try:
            os.killpg(os.getpgid(proc.pid), 9)
        except (OSError, ProcessLookupError):
            try:
                proc.kill()
            except OSError:
                pass


def kill_children() -> int:
    """Stop every converter process this program started. Called when the program quits."""
    with _ACTIVE_LOCK:
        procs = list(_ACTIVE)
    for p in procs:
        _kill(p)
    return len(procs)


class SpectrumSource:
    kind = "?"
    name = ""

    def info(self) -> dict:
        raise NotImplementedError

    def xic(self, targets: list[tuple[float, float]]) -> list[tuple[np.ndarray, np.ndarray]]:
        raise NotImplementedError

    def ms2_around(self, rt0: float, rt1: float, mz: float) -> list[Spectrum]:
        raise NotImplementedError


def find_thermo_parser() -> str | None:
    cands = [os.environ.get("RD_THERMO_PARSER", ""),
             str(Path.home() / "Library/Application Support/MassSpecBench/tools/ThermoRawFileParser/ThermoRawFileParser"),
             str(Path.home() / "Library/Application Support/RecodeDetector/tools/ThermoRawFileParser/ThermoRawFileParser"),
             shutil.which("ThermoRawFileParser") or ""]
    for c in cands:
        if c and Path(c).is_file() and os.access(c, os.X_OK):
            return c
    return None


class RawSource(SpectrumSource):
    """Thermo .raw read in place. Needs ThermoRawFileParser (found automatically, or set RD_THERMO_PARSER)."""
    kind = "raw"

    def __init__(self, path: str, tool: str | None = None, cancel: threading.Event | None = None):
        self.path = str(path)
        self.name = Path(path).name
        self.tool = tool or find_thermo_parser()
        self.cancel = cancel
        if not Path(self.path).is_file():
            raise ValueError(f"File not found: {self.path}")
        if not self.tool:
            raise ValueError("ThermoRawFileParser was not found. Open MassSpec Bench once (it installs it), or set "
                             "RD_THERMO_PARSER to the program's path, or convert the file to mzML first.")
        self._info = None
        self._grid = None
        self._rt_div = 1.0           # the converter's scan query reports retention time in seconds; detected from the run length
        self._blocks: dict[tuple[int, int], list[Spectrum]] = {}

    def _run(self, args: list[str], timeout: int = 900) -> bytes:
        """Run the converter, watching for cancel and registering the process so that quitting the program kills it."""
        import time
        if self.cancel is not None and self.cancel.is_set():
            raise InterruptedError("cancelled")
        try:
            proc = subprocess.Popen([self.tool, *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
        except OSError as e:
            raise ValueError(f"ThermoRawFileParser could not be started: {e}")
        with _ACTIVE_LOCK:
            _ACTIVE.add(proc)
        out_chunks, err = [], b""
        try:
            t0 = time.time()
            reader = threading.Thread(target=lambda: out_chunks.append(proc.stdout.read()), daemon=True)
            reader.start()
            while reader.is_alive():
                reader.join(0.4)
                if self.cancel is not None and self.cancel.is_set():
                    _kill(proc)
                    raise InterruptedError("cancelled")
                if time.time() - t0 > timeout:
                    _kill(proc)
                    raise ValueError("Reading the RAW file took too long.")
            err = proc.stderr.read()
            proc.wait()
        finally:
            _kill(proc)
            with _ACTIVE_LOCK:
                _ACTIVE.discard(proc)
        out = b"".join(out_chunks)
        if proc.returncode != 0 and not out:
            msg = err.decode(errors="replace").strip().splitlines()[-1:] or ["no message"]
            raise ValueError(f"ThermoRawFileParser failed: {msg[0][:200]}")
        return out

    def info(self) -> dict:
        if self._info is None:
            with tempfile.TemporaryDirectory() as t:
                self._run(["-i", self.path, "-f", "4", "-m", "0", "-o", t, "-l", "0"], 600)
                files = list(Path(t).glob("*metadata*.json"))
                if not files:
                    raise ValueError("Could not read the RAW file's header.")
                d = json.loads(files[0].read_text())
            flat = {x.get("name"): x.get("value") for sec in d.values() if isinstance(sec, list) for x in sec if isinstance(x, dict)}
            self._info = {"n_scans": int(float(flat.get("Number of scans", 0))), "n_ms1": int(float(flat.get("Number of MS1 spectra", 0))),
                          "n_ms2": int(float(flat.get("Number of MS2 spectra", 0))), "rt_max": float(flat.get("MS max RT", 0)),
                          "instrument": flat.get("Thermo Scientific instrument model", ""), "name": self.name}
            if self._info["n_scans"] < 10:
                raise ValueError("The RAW file has no scans.")
        return self._info

    def _query(self, spec: str) -> list[Spectrum]:
        out = self._run(["query", "-i", self.path, "-n", spec, "-s"], 900)
        try:
            data = json.loads(out)
        except ValueError:
            raise ValueError("ThermoRawFileParser returned unreadable data for the scan query.")
        res = []
        for d in data:
            attrs = {a.get("name"): a.get("value") for a in d.get("attributes", [])}
            sp = _spectrum_from_attrs(attrs, d.get("mzs", []), d.get("intensities", []))
            sp.rt /= self._rt_div
            res.append(sp)
        return res

    def _scan_grid(self):
        if self._grid is None:
            n = self.info()["n_scans"]
            pts = np.unique(np.linspace(1, n, min(400, n)).astype(int))
            sp = self._query(",".join(map(str, pts)))
            last, rt_max = max(s.rt for s in sp), self.info()["rt_max"]
            if rt_max > 0 and last > rt_max * 30:           # seconds, not minutes
                self._rt_div = 60.0
                for s_ in sp:
                    s_.rt /= 60.0
            self._grid = (np.array([s.rt for s in sp]), np.array([s.scan for s in sp], float))
        return self._grid

    def scan_for_rt(self, rt: float) -> int:
        rts, scans = self._scan_grid()
        return int(round(float(np.interp(rt, rts, scans))))

    def xic(self, targets):
        with tempfile.TemporaryDirectory() as t:
            jp = Path(t) / "xic.json"
            jp.write_text(json.dumps([{"mz": round(mz, 6), "tolerance": ppm, "tolerance_unit": "ppm"} for mz, ppm in targets]))
            out = self._run(["xic", "-i", self.path, "-j", str(jp), "-s"], 900)
        content = json.loads(out)["Content"]
        return [(np.array(c["RetentionTimes"], float), np.array(c["Intensities"], float)) for c in content]

    def _block(self, rt0: float, rt1: float) -> list[Spectrum]:
        n = self.info()["n_scans"]
        s0 = max(1, self.scan_for_rt(max(rt0, 0)) - 150)
        s1 = min(n, self.scan_for_rt(rt1) + 150)
        key = (s0 // 200, s1 // 200)
        if key not in self._blocks:
            self._blocks[key] = self._query(f"{s0}-{s1}")
            if len(self._blocks) > 6:
                self._blocks.pop(next(iter(self._blocks)))
        return self._blocks[key]

    def ms2_around(self, rt0, rt1, mz):
        return [s for s in self._block(rt0, rt1) if s.level == 2 and s.target is not None
                and rt0 <= s.rt <= rt1 and s.target - s.lo <= mz <= s.target + s.hi]


class MzmlSource(SpectrumSource):
    """mzML or mzML.gz, read fully into memory (for small files)."""
    kind = "mzml"
    MAX_BYTES = 1_500_000_000

    def __init__(self, path: str, cancel: threading.Event | None = None):
        self.path = str(path)
        self.name = Path(path).name
        if not Path(self.path).is_file():
            raise ValueError(f"File not found: {self.path}")
        self.ms1: list[Spectrum] = []
        self.ms2: list[Spectrum] = []
        self.cancel = cancel
        self._load()

    @staticmethod
    def _decode(el, ns_strip) -> tuple[str, np.ndarray]:
        kinds = {ns_strip(c.tag): c for c in el}
        accs = {c.get("accession") for c in el if ns_strip(c.tag) == "cvParam"}
        raw = (el.find(".//binary") if el.find(".//binary") is not None else None)
        for c in el.iter():
            if ns_strip(c.tag) == "binary":
                raw = c
        data = base64.b64decode(raw.text or "") if raw is not None else b""
        if "MS:1000574" in accs:
            data = zlib.decompress(data)
        dt = np.float64 if "MS:1000523" in accs else np.float32
        arr = np.frombuffer(data, dtype=dt).astype(float)
        name = "mz" if "MS:1000514" in accs else "int" if "MS:1000515" in accs else "other"
        return name, arr

    def _load(self):
        import xml.etree.ElementTree as ET
        if Path(self.path).stat().st_size > self.MAX_BYTES:
            raise ValueError("This mzML file is very large. Load the .raw file instead, or convert only the MS levels you need.")
        opener = gzip.open if self.path.lower().endswith(".gz") else open
        strip = lambda t: t.split("}")[-1]       # noqa: E731
        with opener(self.path, "rb") as fh:
            for _, el in ET.iterparse(fh, events=("end",)):
                if strip(el.tag) != "spectrum":
                    continue
                if self.cancel is not None and self.cancel.is_set():
                    raise InterruptedError("cancelled")
                attrs, mz, inten = {}, np.array([]), np.array([])
                for sub in el.iter():
                    t = strip(sub.tag)
                    if t == "cvParam":
                        v = sub.get("value")
                        nm = sub.get("name")
                        if nm == "scan start time" and (sub.get("unitName") or "").lower().startswith("sec"):
                            v = str(float(v) / 60.0)
                        if nm and nm not in attrs:
                            attrs[nm] = v
                    elif t == "binaryDataArray":
                        name, arr = self._decode(sub, strip)
                        if name == "mz":
                            mz = arr
                        elif name == "int":
                            inten = arr
                idx = el.get("index") or "0"
                attrs.setdefault("scan number", str(int(idx) + 1))
                sp = _spectrum_from_attrs(attrs, mz, inten)
                (self.ms1 if sp.level == 1 else self.ms2).append(sp)
                el.clear()
        if not self.ms1 and not self.ms2:
            raise ValueError("No spectra found in this mzML file.")
        self.ms2.sort(key=lambda s: s.rt)
        self.ms1.sort(key=lambda s: s.rt)

    def info(self):
        allsp = self.ms1 + self.ms2
        return {"n_scans": len(allsp), "n_ms1": len(self.ms1), "n_ms2": len(self.ms2),
                "rt_max": max((s.rt for s in allsp), default=0.0), "instrument": "", "name": self.name}

    def xic(self, targets):
        rts = np.array([s.rt for s in self.ms1])
        out = []
        for mz, ppm in targets:
            tol = mz * ppm * 1e-6
            vals = np.zeros(len(self.ms1))
            for i, s in enumerate(self.ms1):
                lo, hi = np.searchsorted(s.mz, [mz - tol, mz + tol])
                if hi > lo:
                    vals[i] = s.inten[lo:hi].sum()
            out.append((rts, vals))
        return out

    def ms2_around(self, rt0, rt1, mz):
        return [s for s in self.ms2 if rt0 <= s.rt <= rt1 and s.target is not None and s.target - s.lo <= mz <= s.target + s.hi]


def open_source(path: str, cancel: threading.Event | None = None) -> SpectrumSource:
    p = str(path).strip().strip('"').strip("'")
    low = p.lower()
    if low.endswith(RAW_EXT):
        return RawSource(p, cancel=cancel)
    if low.endswith(MZML_EXT):
        return MzmlSource(p, cancel=cancel)
    raise ValueError("Use a Thermo .raw file or an mzML / mzML.gz file.")


# ---------------------------------------------------------------- matching and scoring

def match_peaks(sp_mz: np.ndarray, sp_int: np.ndarray, frags: list[Fragment], ppm: float, changed: set[int]):
    """Each fragment takes its nearest peak inside the tolerance; a peak is used once (closest fragment wins)."""
    if len(sp_mz) == 0:
        return []
    cand = []
    for f in frags:
        tol = f.mz * ppm * 1e-6
        i = int(np.searchsorted(sp_mz, f.mz))
        for j in (i - 1, i):
            if 0 <= j < len(sp_mz) and abs(sp_mz[j] - f.mz) <= tol:
                cand.append((abs(sp_mz[j] - f.mz), f, j))
    cand.sort(key=lambda t: t[0])
    used_p, used_f, out = set(), set(), []
    for _, f, j in cand:
        if j in used_p or f.label in used_f:
            continue
        used_p.add(j)
        used_f.add(f.label)
        out.append({"mz": float(sp_mz[j]), "int": float(sp_int[j]), "label": f.label, "kind": f.kind, "theory": f.mz,
                    "ppm": float((sp_mz[j] - f.mz) / f.mz * 1e6), "disc": bool(changed & set(f.covers)), "peak": int(j)})
    return out


def chance_probability(sp_mz: np.ndarray, frags: list[Fragment], k: int, ppm: float) -> float:
    """Probability that at least k fragments would match by chance, given how crowded this spectrum is around each
    fragment (Poisson approximation: each fragment hits a peak with probability 1 - exp(-density * window))."""
    import math
    if k <= 0 or len(sp_mz) == 0:
        return 1.0
    lam = 0.0
    for f in frags:
        a, b = np.searchsorted(sp_mz, [f.mz - 25.0, f.mz + 25.0])
        density = max(b - a, 1) / 50.0                       # peaks per Da around the fragment
        lam += 1.0 - math.exp(-density * 2.0 * f.mz * ppm * 1e-6)
    term, tail = math.exp(-lam), 0.0
    cdf = 0.0
    for i in range(k):
        cdf += term
        term *= lam / (i + 1)
    tail = max(1.0 - cdf, 1e-15)                 # below this the arithmetic is not meaningful
    return float(tail)


def fmt_p(p: float) -> str:
    return "<1e-12" if p <= 1e-12 else f"{p:.0e}"


def _smooth(y: np.ndarray, k: int = 5) -> np.ndarray:
    if len(y) < k:
        return y
    return np.convolve(y, np.ones(k) / k, mode="same")


def _down(rt: np.ndarray, y: np.ndarray, lo: float, hi: float, n: int = 300):
    m = (rt >= lo) & (rt <= hi)
    rt, y = rt[m], y[m]
    if len(rt) > n:
        idx = np.linspace(0, len(rt) - 1, n).astype(int)
        rt, y = rt[idx], y[idx]
    return rt.round(4).tolist(), y.round(1).tolist()


N_DECOYS = 20


def _decoys(seq: str, shifts: list[float], n: int = N_DECOYS) -> list[tuple[str, list[float]]]:
    """Shuffled versions with the same composition (so the same precursor mass) and the same C-terminal residue."""
    rng = np.random.default_rng(sum(ord(c) * (i + 1) for i, c in enumerate(seq)))
    out, seen = [], {seq}
    body = list(range(len(seq) - 1))
    for _ in range(n * 6):
        if len(out) >= n or len(body) < 2:
            break
        idx = list(rng.permutation(body))
        d = "".join(seq[i] for i in idx) + seq[-1]
        if d in seen:
            continue
        seen.add(d)
        out.append((d, [shifts[i] for i in idx] + [shifts[-1]]))
    return out


def _local_peaks(rt: np.ndarray, sm: np.ndarray, n: int) -> list[int]:
    """Indices of the n strongest separate peaks in a smoothed chromatogram (each at least 0.4 min from a stronger one)."""
    order = np.argsort(-sm)
    picked: list[int] = []
    floor = sm.max() * 0.15
    for i in order:
        if sm[i] < floor or len(picked) >= n:
            break
        if all(abs(rt[i] - rt[j]) >= 0.4 for j in picked):
            picked.append(int(i))
    return picked


def _evaluate_apex(source, seq, shifts, nterm, frags, changed, mz, apex_rt, ppm, half_window, cancel):
    """MS/MS evidence at one elution peak. Returns (ms2 dict or None, sort key, n_scans)."""
    scans = source.ms2_around(apex_rt - half_window, apex_rt + half_window, mz)
    scans.sort(key=lambda s: s.rt)
    if not scans:
        return None, (-1, -1.0, -1.0), 0
    disc_frags = [f for f in frags if changed & set(f.covers)]
    best, best_key, per_scan = None, (-1.0, -1.0, -1.0), []
    for s in scans:
        if cancel is not None and cancel.is_set():
            raise InterruptedError("cancelled")
        m = match_peaks(s.mz, s.inten, frags, ppm, changed)
        n_d = sum(1 for x in m if x["disc"])
        tic = float(s.inten.sum()) or 1.0
        expl = sum(x["int"] for x in m) / tic
        p_all = chance_probability(s.mz, frags, len(m), ppm)
        p_disc = chance_probability(s.mz, disc_frags, n_d, ppm) if disc_frags else 1.0
        per_scan.append((s, m))
        key = (-np.log10(p_disc) if disc_frags else 0.0, -np.log10(p_all), expl)
        if key > best_key:
            best, best_key = (s, m, expl, p_all, p_disc), key
    s_best, m_best, expl, p_all, p_disc = best
    # negative control: shuffled peptides of the same mass, each given the same chance to find its best scan
    decoy_p = []
    for dseq, dsh in _decoys(seq, shifts):
        dfr = fragments(dseq, dsh, nterm)
        pm = 1.0
        for s_, _ in per_scan:
            mm = match_peaks(s_.mz, s_.inten, dfr, ppm, set())
            pm = min(pm, chance_probability(s_.mz, dfr, len(mm), ppm))
        decoy_p.append(pm)
    as_good = int(sum(1 for q in decoy_p if q <= p_all))
    top = sorted(m_best, key=lambda x: -x["int"])[:8]
    tr_rt = np.array([s.rt for s, _ in per_scan])
    traces, mats = [], []
    for x in top:
        tol = x["theory"] * ppm * 1e-6
        y = np.zeros(len(per_scan))
        for n_, (s, _) in enumerate(per_scan):
            a, b = np.searchsorted(s.mz, [x["theory"] - tol, x["theory"] + tol])
            if b > a:
                y[n_] = s.inten[a:b].max()
        mats.append(y)
        traces.append({"label": x["label"], "disc": x["disc"], "rt": tr_rt.round(4).tolist(), "intensity": y.round(1).tolist()})
    co = None
    usable = [y for y in mats if (y > 0).sum() >= 3]
    if len(usable) >= 2:
        cs = [np.corrcoef(a, b)[0, 1] for i, a in enumerate(usable) for b in usable[i + 1:] if a.std() > 0 and b.std() > 0]
        co = float(np.mean(cs)) if cs else None
    n_m, n_d = len(m_best), sum(1 for x in m_best if x["disc"])
    tight = sum(1 for x in m_best if abs(x["ppm"]) <= ppm / 2)
    order = np.argsort(-s_best.inten)[:160]
    keep = sorted(set(order.tolist()) | {x["peak"] for x in m_best})
    ms2 = {"n_scans": len(scans), "apex_rt": round(apex_rt, 3),
           "best": {"scan": s_best.scan, "rt": round(s_best.rt, 4), "target": s_best.target,
                    "mz": s_best.mz[keep].round(4).tolist(), "int": s_best.inten[keep].round(1).tolist(),
                    "matches": [{k_: v for k_, v in x.items() if k_ != "peak"} | {"index": keep.index(x["peak"])} for x in m_best],
                    "n_matched": n_m, "n_disc": n_d, "n_tight": tight, "explained": round(expl, 4),
                    "n_fragments_possible": len(frags), "p_chance": p_all, "p_chance_disc": p_disc,
                    "decoys": {"n": len(decoy_p), "as_good": as_good}},
           "traces": traces, "coelution": co}
    return ms2, best_key, len(scans)


def check_peptide(source: SpectrumSource, modified: str, charge: int, changed: list[int] | None = None, ppm: float = 10.0,
                  half_window: float = 0.15, rt_hint: float | None = None, explicit_mods: bool = True,
                  version: str = "", cancel: threading.Event | None = None, ms1_ppm: float = 10.0, n_apex: int = 2) -> dict:
    """Evidence for one peptide in the raw data. `changed` = residue positions (0-based) that differ between versions."""
    seq, shifts, nterm = parse_modified(modified)
    shifts = with_fixed_cam(seq, shifts, explicit_mods and ("(" in modified))
    changed = set(changed or [])
    mz = precursor_mz(seq, shifts, nterm, charge)
    res = {"peptide": seq, "modified": modified, "charge": charge, "precursor_mz": round(mz, 5), "version": version,
           "changed": sorted(changed), "verdict": "not seen", "reasons": [], "ms1": None, "ms2": None}
    (rt, inten), = source.xic([(mz, ms1_ppm)])
    sm = _smooth(inten)
    nz = inten[inten > 0]
    if len(nz) < 3 or sm.max() <= 0:
        res["reasons"].append("no MS1 signal at the expected mass")
        res["verdict"] = "no MS1 signal"
        res["ms1"] = {"rt": [], "intensity": [], "apex_rt": None, "apex_intensity": 0.0, "snr": 0.0, "candidates": []}
        return res
    idxs = [int(np.argmin(np.abs(rt - rt_hint)))] if rt_hint is not None else _local_peaks(rt, sm, n_apex)
    frags = fragments(seq, shifts, nterm)
    best_ms2, best_key, best_i, tried, trials = None, (-2.0, -2.0, -2.0), idxs[0], [], 0
    for k in idxs:
        ms2, key, n_sc = _evaluate_apex(source, seq, shifts, nterm, frags, changed, mz, float(rt[k]), ppm, half_window, cancel)
        tried.append(round(float(rt[k]), 3))
        trials += n_sc
        if ms2 is not None and key > best_key:
            best_ms2, best_key, best_i = ms2, key, k
    k = best_i
    apex_rt, apex_int = float(rt[k]), float(inten[max(0, k - 2):k + 3].max())
    r_, i_ = _down(rt, inten, max(rt.min(), apex_rt - 1.0), apex_rt + 1.0)
    res["ms1"] = {"rt": r_, "intensity": i_, "apex_rt": round(apex_rt, 3), "apex_intensity": apex_int,
                  "snr": round(apex_int / float(np.median(nz)), 1), "candidates": tried}
    if best_ms2 is None:
        res["reasons"].append("no fragmentation scan isolated this mass at the peak")
        res["verdict"] = "no MS/MS"
        return res
    res["ms2"] = best_ms2
    bm = best_ms2["best"]
    # we kept the best of `trials` scans, so each chance probability is multiplied by that number (Bonferroni)
    bm["p_chance"], bm["p_chance_disc"] = min(1.0, bm["p_chance"] * trials), min(1.0, bm["p_chance_disc"] * trials)
    bm["n_trials"] = trials
    n_m, n_d, p_all, p_disc, co = bm["n_matched"], bm["n_disc"], bm["p_chance"], bm["p_chance_disc"], best_ms2["coelution"]
    # the verdict rests on the fragments that contain the changed residue: fragments from the unchanged part are
    # shared by both versions and cannot tell them apart
    dec = bm.get("decoys", {"n": 0, "as_good": 0})
    if not changed:                                   # a peptide checked on its own: every ion counts
        res["verdict"] = "supported" if (p_all <= 1e-4 and n_m >= 5 and dec["as_good"] <= 1) else "weak" if (p_all <= 1e-2 and n_m >= 3) else "not seen"
    elif n_d >= 2 and p_disc <= 1e-3 and p_all <= 1e-3 and dec["as_good"] <= 1:
        res["verdict"] = "supported"
    elif n_d >= 1 and p_disc <= 0.05 and p_all <= 1e-2:
        res["verdict"] = "weak"
    res["reasons"] = [f"{n_m} fragment ions matched ({bm['n_tight']} within half the tolerance); {n_d} contain the changed residue",
                      f"chance of this many matches by accident (corrected for choosing the best of {trials} scans): {fmt_p(p_all)}; "
                      f"for the changed-residue ions alone: {fmt_p(p_disc)}",
                      f"{dec['as_good']} of {dec['n']} shuffled decoy peptides of the same mass matched this well on the same scans"
                      + (" (a good sign: decoys should not match)" if dec["as_good"] == 0 else ""),
                      f"fragment traces co-elute (r = {co:.2f})" if co is not None else "too few fragment traces to judge co-elution"]
    return res


def compare_versions(alt: dict, ref: dict) -> dict:
    """Which version do the spectra favour? A plain label plus the numbers it rests on."""
    def strength(r):
        return {"supported": 2, "weak": 1}.get(r["verdict"], 0)
    a, b = strength(alt), strength(ref)
    if a > b:
        label, text = "alternative", "The spectra favour the alternative version."
    elif b > a:
        label, text = "reference", "The spectra favour the reference version."
    elif a == 2:
        label, text = "both", "Both versions look supported; the peptides may co-elute or share fragments."
    else:
        label, text = "neither", "Neither version is clearly supported by the spectra."
    return {"label": label, "text": text}


# ---------------------------------------------------------------- job for the app

@dataclass
class SpectraJob:
    id: str
    status: str = "queued"           # queued, running, done, error, cancelled
    message: str = ""
    done: int = 0
    total: int = 0
    results: list = field(default_factory=list)
    error: str = ""
    cancel: threading.Event = field(default_factory=threading.Event)
    source_info: dict = field(default_factory=dict)

    def public(self) -> dict:
        return {"id": self.id, "status": self.status, "message": self.message, "done": self.done, "total": self.total,
                "error": self.error, "results": self.results, "source": self.source_info}


def run_job(job: SpectraJob, path: str, targets: list[dict], ppm: float, half_window: float) -> None:
    """targets: [{"key", "site", "version", "modified", "charge", "changed", "pair"}]"""
    try:
        job.status, job.message = "running", "Opening the file…"
        src = open_source(path, job.cancel)
        job.source_info = src.info()
        job.total = len(targets)
        for t in targets:
            if job.cancel.is_set():
                raise InterruptedError("cancelled")
            job.message = f"Checking {t.get('modified')} ({t.get('version')}), {job.done + 1} of {job.total}…"
            try:
                r = check_peptide(src, t["modified"], int(t["charge"]), t.get("changed"), ppm, half_window,
                                  rt_hint=t.get("rt"), version=t.get("version", ""), cancel=job.cancel)
            except InterruptedError:
                raise
            except Exception as e:                       # one bad peptide must not sink the rest
                r = {"peptide": t.get("modified"), "charge": t.get("charge"), "version": t.get("version", ""),
                     "verdict": "error", "reasons": [str(e)], "ms1": None, "ms2": None, "precursor_mz": None}
            r.update({"key": t.get("key"), "site": t.get("site"), "pair": t.get("pair")})
            job.results.append(r)
            job.done += 1
        job.status, job.message = "done", f"Checked {job.done} peptides."
    except InterruptedError:
        job.status, job.message = "cancelled", "Stopped."
    except Exception as e:
        job.status, job.error = "error", str(e)


def spectra_frames(results: list[dict], sites: list[dict]):
    """Two tables for export: one row per checked peptide, and one row per matched fragment ion."""
    import pandas as pd
    rows, ions = [], []
    for r in results:
        site = sites[r["site"]] if r.get("site") is not None and r["site"] < len(sites) else None
        m1, m2 = r.get("ms1") or {}, r.get("ms2") or {}
        bm = (m2.get("best") or {})
        row = {"protein": site["ref_protein"] if site else "", "position": site["position"] if site else None,
               "version": r.get("version"), "peptide": r.get("modified"), "charge": r.get("charge"),
               "precursor_mz": r.get("precursor_mz"), "ms1_apex_rt_min": m1.get("apex_rt"),
               "ms1_apex_intensity": m1.get("apex_intensity"), "ms1_signal_to_median": m1.get("snr"),
               "ms2_scans_at_peak": m2.get("n_scans"), "best_scan": bm.get("scan"), "fragments_matched": bm.get("n_matched"),
               "changed_residue_fragments_matched": bm.get("n_disc"), "chance_probability": bm.get("p_chance"),
               "chance_probability_changed_residue_ions": bm.get("p_chance_disc"), "decoys_matching_as_well": (bm.get("decoys") or {}).get("as_good"), "decoys_tested": (bm.get("decoys") or {}).get("n"), "fragment_coelution_r": m2.get("coelution"),
               "verdict": r.get("verdict"), "notes": " | ".join(r.get("reasons", []))}
        rows.append(row)
        for x in bm.get("matches", []):
            ions.append({"protein": row["protein"], "position": row["position"], "version": row["version"], "peptide": row["peptide"],
                         "charge": row["charge"], "scan": bm.get("scan"), "ion": x["label"], "observed_mz": x["mz"],
                         "theoretical_mz": x["theory"], "error_ppm": round(x["ppm"], 2), "intensity": x["int"],
                         "contains_changed_residue": x["disc"]})
    return pd.DataFrame(rows), pd.DataFrame(ions)

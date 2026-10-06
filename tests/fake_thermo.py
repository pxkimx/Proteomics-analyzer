#!/usr/bin/env python3
"""A stand-in for ThermoRawFileParser used in tests. It serves a simulated mzML (path in FAKE_MZML) through the same
command line as the real tool, including the real tool's habit of reporting scan start times in seconds from `query`."""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np  # noqa: E402

from recode_detector import msms  # noqa: E402

src = msms.MzmlSource(os.environ["FAKE_MZML"])
allsp = sorted(src.ms1 + src.ms2, key=lambda s: s.scan)
args = sys.argv[1:]


def opt(flag):
    return args[args.index(flag) + 1] if flag in args else None


if args and args[0] == "xic":
    req = json.loads(Path(opt("-j")).read_text())
    out = [{"Meta": {}, "RetentionTimes": r.tolist(), "Intensities": i.tolist()}
           for r, i in src.xic([(q["mz"], q["tolerance"]) for q in req])]
    print(json.dumps({"OutputMeta": {"base64": False, "timeunit": "minutes"}, "Content": out}))
elif args and args[0] == "query":
    wanted = []
    for part in opt("-n").split(","):
        if "-" in part:
            a, b = part.split("-")
            wanted += range(int(a), int(b) + 1)
        else:
            wanted.append(int(part))
    by = {s.scan: s for s in allsp}
    res = []
    for n in wanted:
        s = by.get(n)
        if s is None:
            continue
        at = [{"name": "scan number", "value": str(s.scan)}, {"name": "scan start time", "value": str(s.rt * 60.0)},   # seconds!
              {"name": "ms level", "value": str(s.level)}]
        if s.level == 2:
            at += [{"name": "isolation window target m/z", "value": str(s.target)},
                   {"name": "isolation window lower offset", "value": str(s.lo)}, {"name": "isolation window upper offset", "value": str(s.hi)}]
        res.append({"mzs": s.mz.tolist(), "intensities": s.inten.tolist(), "attributes": at})
    print(json.dumps(res))
else:                                                       # header: -i raw -f 4 -m 0 -o dir
    out = Path(opt("-o"))
    rt_max = max(s.rt for s in allsp)
    meta = {"MsData": [{"name": "Number of MS1 spectra", "value": str(len(src.ms1))},
                       {"name": "Number of MS2 spectra", "value": str(len(src.ms2))},
                       {"name": "MS max RT", "value": str(rt_max)}],
            "ScanSettings": [{"name": "Number of scans", "value": str(len(allsp))}],
            "InstrumentProperties": [{"name": "Thermo Scientific instrument model", "value": "Simulated Astral"}]}
    (out / "x-metadata.json").write_text(json.dumps(meta))

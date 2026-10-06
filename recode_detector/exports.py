"""Turn any table into a download: CSV, TSV or Excel, or everything at once as a ZIP."""
from __future__ import annotations

import io
import json
import re
import zipfile

import pandas as pd

FORMATS = {"csv": "text/csv", "tsv": "text/tab-separated-values",
           "xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}


def _plain(df: pd.DataFrame) -> pd.DataFrame:
    """Lists and dicts become text so every format can hold them."""
    out = df.copy()
    for c in out.columns:
        if out[c].map(lambda v: isinstance(v, (list, dict, tuple, set))).any():
            out[c] = out[c].map(lambda v: json.dumps(v, default=str) if isinstance(v, (list, dict, tuple, set)) else v)
    return out


def frame_bytes(df: pd.DataFrame, fmt: str = "csv", sheet: str = "Table") -> bytes:
    fmt = (fmt or "csv").lower()
    if fmt not in FORMATS:
        raise ValueError(f"Unknown format “{fmt}”. Use csv, tsv or xlsx.")
    df = _plain(df)
    if fmt == "xlsx":
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as xw:
            df.to_excel(xw, index=False, sheet_name=re.sub(r"[\[\]:*?/\\]", "_", sheet)[:31] or "Table")
            ws = xw.sheets[next(iter(xw.sheets))]
            for i, col in enumerate(df.columns, 1):                      # readable column widths
                width = min(60, max(len(str(col)), *(len(str(v)) for v in df[col].head(200))) + 2) if len(df) else len(str(col)) + 2
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = width
            ws.freeze_panes = "A2"
        return buf.getvalue()
    return df.to_csv(index=False, sep="\t" if fmt == "tsv" else ",").encode()


def zip_bytes(tables: dict[str, pd.DataFrame], extras: dict[str, bytes] | None = None, fmt: str = "csv") -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, df in tables.items():
            if df is not None and len(df.columns):
                z.writestr(f"{name}.{fmt}", frame_bytes(df, fmt, name))
        for name, data in (extras or {}).items():
            z.writestr(name, data)
    return buf.getvalue()

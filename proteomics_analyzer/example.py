"""A SIMULATED MaxQuant-style proteinGroups table, for trying the program and for tests.

Nothing here is real biology: protein and gene names are placeholders (SIM0001 ...), and the
differential proteins are chosen by the generator, which also returns them as ground truth.
"""
from __future__ import annotations

import io

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

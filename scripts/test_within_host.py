#!/usr/bin/env python3
"""Permutation test: within-host vs between-host H. pylori structural divergence.

Uses the cohort pair table (syn2b_structural_pairs_raw.tsv). Within-host pairs
are those with the same participant (same_group == 1); all others are
between-host. For each structural metric we report medians, the Mann-Whitney U
statistic, a label-permutation p-value (10,000 permutations of the
within/between assignment across pairs), and Cliff's delta as effect size.

Output: results/metric_validation/h_pylori_host_permutation.tsv
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
PAIRS = ROOT / "data" / "syntracker_validation" / "syn2b_structural_raw" / "syn2b_structural_pairs_raw.tsv"
OUT = ROOT / "results" / "metric_validation" / "h_pylori_host_permutation.tsv"

METRICS = ["syn2b_breakpoints", "syn2b_scj_distance", "syn2b_observable_fraction"]
N_PERM = 10_000
SEED = 42


def cliffs_delta(x, y):
    """Cliff's delta: P(x > y) - P(x < y) for x = within-host, y = between-host."""
    x = np.asarray(x)
    y = np.asarray(y)
    gt = (x[:, None] > y[None, :]).mean()
    lt = (x[:, None] < y[None, :]).mean()
    return gt - lt


def main():
    df = pd.read_csv(PAIRS, sep="\t")
    df = df[(df["cohort"] == "Helicobacter_pylori") & (df["status"] == "ok") & (df["is_self"] != 1)].copy()

    rows = []
    rng = np.random.default_rng(SEED)
    n_within = int(df["same_group"].sum())

    for metric in METRICS:
        within = df.loc[df["same_group"] == 1, metric].values
        between = df.loc[df["same_group"] == 0, metric].values

        u, p_mw = stats.mannwhitneyu(within, between, alternative="two-sided")

        observed = np.median(within) - np.median(between)
        values = df[metric].values
        perm_diffs = np.empty(N_PERM)
        for i in range(N_PERM):
            idx = rng.permutation(len(values))[:n_within]
            mask = np.zeros(len(values), dtype=bool)
            mask[idx] = True
            perm_diffs[i] = np.median(values[mask]) - np.median(values[~mask])
        p_perm = (np.sum(np.abs(perm_diffs) >= abs(observed)) + 1) / (N_PERM + 1)

        rows.append({
            "metric": metric,
            "n_within": len(within),
            "n_between": len(between),
            "median_within": np.median(within),
            "median_between": np.median(between),
            "median_diff": observed,
            "cliffs_delta": cliffs_delta(within, between),
            "mannwhitney_u": u,
            "p_mannwhitney": p_mw,
            "p_permutation": p_perm,
            "n_permutations": N_PERM,
        })

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()

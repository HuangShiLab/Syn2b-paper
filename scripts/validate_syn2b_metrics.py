#!/usr/bin/env python3
"""
validate_syn2b_metrics.py

Comprehensive validation of Syn2b structural metrics across GTDB-R207,
high-ANI subsets, SynTracker isolate cohorts, and closed genomes.
Writes summary tables to results/metric_validation/.
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats


def ci_r(r, n, alpha=0.05):
    """Fisher z-transform confidence interval for Pearson r."""
    if n <= 3:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    lo, hi = z - z_crit * se, z + z_crit * se
    return float(np.tanh(lo)), float(np.tanh(hi))


def pearson_summary(x, y, label=""):
    """Return dict with Pearson r, CI, n, MAE, RMSE, slope, intercept."""
    mask = x.notna() & y.notna()
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 3:
        return {"label": label, "n": n}
    r, p = stats.pearsonr(x, y)
    lo, hi = ci_r(r, n)
    mae = (y - x).abs().mean()
    rmse = np.sqrt(((y - x) ** 2).mean())
    slope, intercept, _, _, _ = stats.linregress(x, y)
    return {
        "label": label,
        "n": n,
        "r": r,
        "r_lo": lo,
        "r_hi": hi,
        "p": p,
        "mae": mae,
        "rmse": rmse,
        "slope": slope,
        "intercept": intercept,
        "mean_x": x.mean(),
        "mean_y": y.mean(),
    }


def analyze_gtdb_four():
    """Main GTDB-R207 held-out set, four-enzyme panel vs dnadiff."""
    df = pd.read_csv("results/gtdb50k/inverted_fraction_truth_four.tsv", sep="\t")
    records = []
    records.append(pearson_summary(
        df["dnadiff_inverted_fraction"], df["syn2b_raw_inverted_fraction"],
        "raw_inverted_fraction_all"
    ))
    records.append(pearson_summary(
        df["dnadiff_inverted_fraction"], df["syn2b_inverted_fraction"],
        "min_inverted_fraction_all"
    ))
    for band, sub in df.groupby("band"):
        if len(sub) < 30:
            continue
        records.append(pearson_summary(
            sub["dnadiff_inverted_fraction"], sub["syn2b_raw_inverted_fraction"],
            f"raw_{band}"
        ))
    return pd.DataFrame(records)


def analyze_high_ani():
    """High-ANI subset: demonstrate raw vs min convention."""
    syn = pd.read_csv("results/gtdb50k/syn2b_inverted_fraction_high_ani_all.tsv", sep="\t")
    dna = pd.read_csv("results/gtdb50k/dnadiff_inverted_fraction_high_ani_all.tsv", sep="\t")
    merged = syn.merge(dna, on="pairid")
    records = []
    records.append(pearson_summary(
        merged["dnadiff_inverted_fraction"], merged["syn2b_raw_inverted_fraction"],
        "high_ani_all_raw"
    ))
    records.append(pearson_summary(
        merged["dnadiff_inverted_fraction"], merged["syn2b_inverted_fraction"],
        "high_ani_all_min"
    ))
    # Subset by dnadiff inverted fraction to show saturation effect
    saturated = merged[merged["dnadiff_inverted_fraction"] > 0.5]
    if len(saturated) > 30:
        records.append(pearson_summary(
            saturated["dnadiff_inverted_fraction"], saturated["syn2b_raw_inverted_fraction"],
            "high_ani_all_raw_saturated"
        ))
        records.append(pearson_summary(
            saturated["dnadiff_inverted_fraction"], saturated["syn2b_inverted_fraction"],
            "high_ani_all_min_saturated"
        ))
    return pd.DataFrame(records)


def analyze_syntracker_cohorts():
    """SynTracker isolate cohorts: self-control and structural signal."""
    pairs = pd.read_csv(
        "data/syntracker_validation/syn2b_structural_raw/syn2b_structural_pairs_raw.tsv",
        sep="\t",
    )
    records = []
    for cohort, sub in pairs.groupby("cohort"):
        # SynTracker cohort files use same_group/same_patient to mark
        # within-host comparisons; is_self is not populated.
        group_col = None
        for col in ("same_group", "same_patient", "is_self"):
            if col in sub.columns and sub[col].astype(bool).sum() > 0:
                group_col = col
                break
        if group_col is not None:
            self_mask = sub[group_col].astype(bool)
        else:
            self_mask = sub["status"] == "self" if "status" in sub.columns else pd.Series(False, index=sub.index)
        self_rows = sub[self_mask]
        cross_rows = sub[~self_mask]
        records.append({
            "label": f"{cohort}_self_floor_bp",
            "n": len(self_rows),
            "median_bp": self_rows["syn2b_breakpoints"].median(),
            "median_scj": self_rows["syn2b_scj_distance"].median(),
            "median_obs_frac": self_rows["syn2b_observable_fraction"].median(),
        })
        records.append({
            "label": f"{cohort}_cross_median_bp",
            "n": len(cross_rows),
            "median_bp": cross_rows["syn2b_breakpoints"].median(),
            "median_scj": cross_rows["syn2b_scj_distance"].median(),
            "median_obs_frac": cross_rows["syn2b_observable_fraction"].median(),
        })
    return pd.DataFrame(records)


def analyze_closed_genomes():
    """Closed-genome inversion validation."""
    f = "results/closed_inversions/position_agreement_allpairs.tsv"
    if not os.path.exists(f):
        return pd.DataFrame()
    df = pd.read_csv(f, sep="\t")
    records = []
    for col in df.select_dtypes(include=[np.number]).columns:
        records.append({
            "label": f"closed_{col}",
            "n": df[col].notna().sum(),
            "median": df[col].median(),
            "mean": df[col].mean(),
        })
    return pd.DataFrame(records)


def analyze_fracminhash():
    """Compare enzyme panel to FracMinHash at multiple densities."""
    records = []
    for label, f in [
        ("bcgI", "results/gtdb50k/inverted_fraction_truth_bcgI.tsv"),
        ("four", "results/gtdb50k/inverted_fraction_truth_four.tsv"),
        ("fmh250", "results/gtdb50k/inverted_fraction_truth_fmh250.tsv"),
        ("fmh750", "results/gtdb50k/inverted_fraction_truth_fmh750.tsv"),
        ("fmh2000", "results/gtdb50k/inverted_fraction_truth_fmh2000.tsv"),
        ("fmh6000", "results/gtdb50k/inverted_fraction_truth_fmh6000.tsv"),
    ]:
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, sep="\t")
        rec = pearson_summary(
            df["dnadiff_inverted_fraction"], df["syn2b_raw_inverted_fraction"],
            label
        )
        rec["median_shared_tags"] = df["syn2b_shared_tags"].median()
        records.append(rec)
    return pd.DataFrame(records)


def main():
    out_dir = "results/metric_validation"
    os.makedirs(out_dir, exist_ok=True)

    gtdb = analyze_gtdb_four()
    gtdb.to_csv(f"{out_dir}/gtdb_four_validation.tsv", sep="\t", index=False)
    print("GTDB four:", gtdb[["label", "n", "r", "mae"]].to_string(index=False))

    high = analyze_high_ani()
    high.to_csv(f"{out_dir}/high_ani_validation.tsv", sep="\t", index=False)
    print("\nHigh ANI:", high[["label", "n", "r", "mae"]].to_string(index=False))

    cohorts = analyze_syntracker_cohorts()
    cohorts.to_csv(f"{out_dir}/syntracker_cohort_summary.tsv", sep="\t", index=False)
    print("\nSynTracker cohorts:", cohorts.to_string(index=False))

    closed = analyze_closed_genomes()
    if not closed.empty:
        closed.to_csv(f"{out_dir}/closed_genome_summary.tsv", sep="\t", index=False)
        print("\nClosed genomes:", closed.to_string(index=False))

    fmh = analyze_fracminhash()
    fmh.to_csv(f"{out_dir}/fracminhash_comparison.tsv", sep="\t", index=False)
    print("\nFracMinHash:", fmh[["label", "n", "r", "mae", "median_shared_tags"]].to_string(index=False))

    print(f"\nWrote summary tables to {out_dir}/")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate Supplementary Figure 2: closed-genome inversion validation."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "closed_inversions" / "position_agreement_allpairs.tsv"
OUT_PNG = ROOT / "figures" / "supplementary" / "fig2_closed_genome_validation.png"
OUT_PDF = ROOT / "figures" / "supplementary" / "fig2_closed_genome_validation.pdf"


def main():
    parser = argparse.ArgumentParser(description="Generate Supplementary Figure 2.")
    parser.add_argument("--input", default=str(DATA))
    parser.add_argument("--out-png", default=str(OUT_PNG))
    parser.add_argument("--out-pdf", default=str(OUT_PDF))
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    df = df.dropna(subset=["n_syn2b", "n_dnadiff", "n_matched"])

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)

    # (a) Syn2b vs dnadiff junction counts
    ax = axes[0]
    max_n = max(df["n_dnadiff"].max(), df["n_syn2b"].max())
    ax.plot([0, max_n], [0, max_n], color="gray", ls="--", lw=1)
    ax.scatter(df["n_dnadiff"], df["n_syn2b"], alpha=0.5, s=20, color="steelblue")
    ax.set_xlabel("dnadiff junctions")
    ax.set_ylabel("Syn2b junctions")
    ax.set_title("a  Junction count agreement")
    ax.set_xlim(0, max_n + 1)
    ax.set_ylim(0, max_n + 1)

    # (b) recovery fraction
    ax = axes[1]
    recovery = df["n_matched"] / df["n_dnadiff"].replace(0, np.nan)
    recovery = recovery.dropna()
    ax.hist(recovery, bins=np.linspace(0, 1, 21), color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(recovery.median(), color="darkred", ls="--", lw=1.5,
               label=f"median = {recovery.median():.2f}")
    ax.set_xlabel("Matched / dnadiff junctions")
    ax.set_ylabel("Pairs")
    ax.set_title("b  Recovery of dnadiff junctions")
    ax.legend(loc="upper left")

    # (c) median distance of matched junctions
    ax = axes[2]
    median_dists = df["median_dist"].dropna()
    median_dists = median_dists[median_dists > 0]
    bins = np.logspace(np.log10(median_dists.min()), np.log10(median_dists.max()), 40)
    ax.hist(median_dists, bins=bins, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(median_dists.median(), color="darkred", ls="--", lw=1.5,
               label=f"median = {median_dists.median():.0f} bp")
    ax.set_xscale("log")
    ax.set_xlabel("Median distance (bp)")
    ax.set_ylabel("Pairs")
    ax.set_title("c  Junction position error")
    ax.legend(loc="upper right")

    fig.suptitle("Supplementary Figure 2 | Closed-genome inversion validation (n = "
                 f"{len(df)} pairs)", fontsize=12)

    out_png = Path(args.out_png)
    out_pdf = Path(args.out_pdf)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate Supplementary Figure 1: tag spacing distribution for the production panel."""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "ecoli_k12_MG1655.fasta"
OUT_PNG = ROOT / "figures" / "supplementary" / "fig1_tag_spacing.png"
OUT_PDF = ROOT / "figures" / "supplementary" / "fig1_tag_spacing.pdf"
PANEL = "BcgI,AlfI,AloI,FalI"


def digest(syn2b, fasta, tgt):
    subprocess.run([str(syn2b), "digest", "--enzymes", PANEL,
                    "--input", str(fasta), "--output", str(tgt)],
                   check=True, capture_output=True)


def parse_distances(tgt_path):
    distances = []
    with open(tgt_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">") or line.startswith("#") or not line.strip():
                continue
            # distances are encoded as " -<int>- " between consecutive landmarks
            distances.extend(int(x) for x in re.findall(r" -(\d+)- ", line))
    return np.array(distances, dtype=float)


def main():
    parser = argparse.ArgumentParser(description="Generate Supplementary Figure 1.")
    parser.add_argument("--syn2b", default=str(ROOT.parent / "Syn2b" / "target" / "release" / "syn2b"))
    parser.add_argument("--input", default=str(DATA))
    parser.add_argument("--out-png", default=str(OUT_PNG))
    parser.add_argument("--out-pdf", default=str(OUT_PDF))
    args = parser.parse_args()

    syn2b = Path(args.syn2b)
    if not syn2b.exists():
        sys.exit(f"syn2b binary not found: {syn2b}")

    with tempfile.NamedTemporaryFile(suffix=".tgt", delete=False) as tmp:
        tgt_path = tmp.name
    try:
        digest(syn2b, args.input, tgt_path)
        dists = parse_distances(tgt_path)
    finally:
        os.unlink(tgt_path)

    median = np.median(dists)
    mean = np.mean(dists)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)

    ax = axes[0]
    bins = np.logspace(np.log10(max(1, dists.min())), np.log10(dists.max()), 60)
    ax.hist(dists, bins=bins, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(median, color="darkred", ls="--", lw=1.5, label=f"median = {median:.0f} bp")
    ax.axvline(mean, color="darkorange", ls="--", lw=1.5, label=f"mean = {mean:.0f} bp")
    ax.set_xscale("log")
    ax.set_xlabel("Inter-tag distance (bp)")
    ax.set_ylabel("Count")
    ax.set_title("a  Tag spacing histogram")
    ax.legend(loc="upper right")

    ax = axes[1]
    sorted_dists = np.sort(dists)
    cdf = np.arange(1, len(sorted_dists) + 1) / len(sorted_dists)
    ax.plot(sorted_dists, cdf, color="steelblue", lw=1.5)
    ax.axvline(median, color="darkred", ls="--", lw=1.5)
    ax.set_xscale("log")
    ax.set_xlabel("Inter-tag distance (bp)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("b  Cumulative spacing distribution")
    ax.set_ylim(0, 1)

    fig.suptitle("Supplementary Figure 1 | Tag spacing distribution for BcgI+AlfI+AloI+FalI in *E. coli* K-12",
                 fontsize=12)

    out_png = Path(args.out_png)
    out_pdf = Path(args.out_pdf)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")
    print(f"n = {len(dists)}, median = {median:.1f} bp, mean = {mean:.1f} bp")


if __name__ == "__main__":
    main()

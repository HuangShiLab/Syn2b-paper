#!/usr/bin/env python3
"""Generate Supplementary Figure 4: inversion detection-size resolution.

Plots the per-size detection rate (both junctions recovered) for the production
four-enzyme panel and BcgI alone, from results/detection_size_metrics.csv
(produced by scripts/measure_detection_size.py).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "results" / "detection_size_metrics.csv"
OUT_PNG = ROOT / "figures" / "supplementary" / "fig4_detection_size.png"
OUT_PDF = ROOT / "figures" / "supplementary" / "fig4_detection_size.pdf"

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "axes.linewidth": 0.7,
})


def main():
    df = pd.read_csv(CSV)
    fig, ax = plt.subplots(figsize=(3.6, 3.0), layout="constrained")

    styles = {
        "panel": ("BcgI+AlfI+AloI+FalI", "#1f77b4", "o"),
        "BcgI": ("BcgI alone", "#d62728", "s"),
    }
    for arm, (label, color, marker) in styles.items():
        sub = df[df["enzyme_arm"] == arm].sort_values("event_size_bp")
        ax.plot(sub["event_size_bp"] / 1000, sub["detection_rate"],
                marker=marker, color=color, lw=1.8, markersize=6, label=label)
        for _, row in sub.iterrows():
            ax.annotate(f"{int(row['detected'])}/{int(row['replicates'])}",
                        (row["event_size_bp"] / 1000, row["detection_rate"]),
                        textcoords="offset points", xytext=(0, 6),
                        ha="center", fontsize=6.5, color=color)

    ax.axhline(0.5, color="gray", ls="--", lw=0.8)
    ax.text(2.05, 0.52, "50% detection", fontsize=7, color="gray")
    ax.set_xlabel("Inverted segment size (kb)")
    ax.set_ylabel("Detection rate (both junctions)")
    ax.set_ylim(-0.05, 1.12)
    ax.set_xticks([2, 4, 8, 16, 32])
    ax.legend(loc="lower right", frameon=False)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()

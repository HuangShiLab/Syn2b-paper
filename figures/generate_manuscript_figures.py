#!/usr/bin/env python3
"""
generate_manuscript_figures.py

Generate publication-quality figures for the Syn2b manuscript.
"""

import os
import random
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "figures")
os.makedirs(OUTDIR, exist_ok=True)


def save(fig, name):
    path = os.path.join(OUTDIR, name)
    fig.savefig(path, dpi=300)
    print(f"Saved {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 1: Algorithm schematic
# ---------------------------------------------------------------------------

def fig1_algorithm():
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.3))

    # (a) Type IIB enzyme cut
    ax = axes[0]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("a  Type IIB/IIG tag generation")

    ax.plot([0.1, 0.9], [0.65, 0.65], color="steelblue", lw=2.5)
    ax.plot([0.1, 0.9], [0.55, 0.55], color="steelblue", lw=2.5)
    ax.add_patch(mpatches.Rectangle((0.35, 0.48), 0.3, 0.25, fc="gold", ec="black", lw=0.6))
    ax.text(0.5, 0.60, "motif", ha="center", va="center", fontsize=8)
    ax.annotate("", xy=(0.35, 0.45), xytext=(0.35, 0.35),
                arrowprops=dict(arrowstyle="->", color="red", lw=1))
    ax.annotate("", xy=(0.65, 0.45), xytext=(0.65, 0.35),
                arrowprops=dict(arrowstyle="->", color="red", lw=1))
    ax.text(0.35, 0.26, "27–32 bp tag", ha="center", fontsize=7)
    ax.text(0.65, 0.26, "27–32 bp tag", ha="center", fontsize=7)

    # (b) Inversion flips tag order/orientation
    ax = axes[1]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("b  Inversion reverses tag order")

    y = 0.75
    for i, c in enumerate(["A", "B", "C", "D", "E", "F"]):
        x = 0.12 + i * 0.13
        ax.add_patch(mpatches.Rectangle((x - 0.045, y - 0.045), 0.09, 0.09, fc="lightgreen", ec="black", lw=0.6))
        ax.text(x, y, c, ha="center", va="center", fontsize=8)
    ax.text(0.5, y + 0.13, "Reference", ha="center", fontsize=8)

    y = 0.32
    tags = [("A", False), ("D", True), ("C", True), ("B", True), ("E", False), ("F", False)]
    for i, (c, inv) in enumerate(tags):
        x = 0.12 + i * 0.13
        fc = "salmon" if inv else "lightgreen"
        ax.add_patch(mpatches.Rectangle((x - 0.045, y - 0.045), 0.09, 0.09, fc=fc, ec="black", lw=0.6))
        ax.text(x, y, c, ha="center", va="center", fontsize=8)
    ax.text(0.5, y - 0.13, "Query (inverted segment shaded)", ha="center", fontsize=8)

    # (c) Fragmentation invariance
    ax = axes[2]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("c  Fragmentation invariance")

    ax.add_patch(mpatches.Rectangle((0.1, 0.70), 0.8, 0.10, fc="lightgray", ec="black", lw=0.6))
    ax.add_patch(mpatches.Rectangle((0.1, 0.70), 0.4, 0.10, fc="salmon", ec="black", lw=0.6))
    ax.text(0.5, 0.88, "True: inverted fraction = 0.4", ha="center", fontsize=8)

    ax.add_patch(mpatches.Rectangle((0.1, 0.40), 0.20, 0.10, fc="lightgray", ec="black", lw=0.6))
    ax.add_patch(mpatches.Rectangle((0.30, 0.40), 0.15, 0.10, fc="salmon", ec="black", lw=0.6))
    ax.add_patch(mpatches.Rectangle((0.45, 0.40), 0.20, 0.10, fc="lightgray", ec="black", lw=0.6))
    ax.add_patch(mpatches.Rectangle((0.65, 0.40), 0.25, 0.10, fc="salmon", ec="black", lw=0.6))
    ax.text(0.5, 0.58, "Draft assembly: ratio still 0.4", ha="center", fontsize=8)
    ax.text(0.5, 0.26, "Transition count: inflated by K−1", ha="center", fontsize=8, color="darkred")

    plt.tight_layout()
    save(fig, "fig1_algorithm_schematic.png")


# ---------------------------------------------------------------------------
# Figure 2: SV sensitivity
# ---------------------------------------------------------------------------

def fig2_sv_sensitivity():
    df = pd.read_csv(os.path.join(ROOT, "data", "enzyme_comparison.csv"))

    enzymes = ["BcgI", "AlfI", "BplI", "CjePI"]
    sv_types = ["control", "inversion_500kb", "translocation_500kb", "insertion_10kb", "deletion_10kb"]
    labels = ["Control", "Inv\n500 kb", "Tra\n500 kb", "Ins\n10 kb", "Del\n10 kb"]

    fig = plt.figure(figsize=(10.0, 5.5))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.55)

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    x = np.arange(len(sv_types))

    # (a) Breakpoint count
    ax = fig.add_subplot(gs[0, 0])
    width = 0.2
    for i, enzyme in enumerate(enzymes):
        sub = df[df["enzyme"] == enzyme].set_index("sv_type").loc[sv_types]
        ax.bar(x + i * width, sub["breakpoints"], width, label=enzyme, color=colors[i])
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Breakpoint count")
    ax.set_title("a  Breakpoints")
    ax.legend(loc="upper left", frameon=False)
    ax.set_ylim(0, 3500)

    # (b) Kendall tau
    ax = fig.add_subplot(gs[0, 1])
    for i, enzyme in enumerate(enzymes):
        sub = df[df["enzyme"] == enzyme].set_index("sv_type").loc[sv_types]
        ax.plot(x, sub["kendall_tau"], marker="o", label=enzyme, color=colors[i], lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Kendall tau")
    ax.set_title("b  Global order")
    ax.set_ylim(0.75, 1.02)
    ax.axhline(1.0, color="gray", ls="--", lw=0.6)

    # (c) Mash proxy distance
    ax = fig.add_subplot(gs[0, 2])
    mash_vals = [df[(df["enzyme"] == e) & (df["sv_type"] == "control")]["mash"].values[0] for e in enzymes]
    inv_mash = [df[(df["enzyme"] == e) & (df["sv_type"] == "inversion_500kb")]["mash"].values[0] for e in enzymes]
    xpos = np.arange(len(enzymes))
    ax.bar(xpos - 0.2, mash_vals, 0.35, label="Control", color="lightgray", ec="black", lw=0.6)
    ax.bar(xpos + 0.2, inv_mash, 0.35, label="500-kb inversion", color="steelblue", ec="black", lw=0.6)
    ax.set_xticks(xpos)
    ax.set_xticklabels(enzymes)
    ax.set_ylabel("Mash proxy distance")
    ax.set_title("c  Mash is blind to SV")
    ax.legend(loc="upper right", frameon=False)

    # (d) 10-kb indel breakpoint counts
    ax = fig.add_subplot(gs[1, 0])
    ins_bp = [df[(df["enzyme"] == e) & (df["sv_type"] == "insertion_10kb")]["breakpoints"].values[0] for e in enzymes]
    del_bp = [df[(df["enzyme"] == e) & (df["sv_type"] == "deletion_10kb")]["breakpoints"].values[0] for e in enzymes]
    xpos = np.arange(len(enzymes))
    width = 0.35
    ax.bar(xpos - width / 2, ins_bp, width, label="10-kb insertion", color="#9467bd", ec="black", lw=0.6)
    ax.bar(xpos + width / 2, del_bp, width, label="10-kb deletion", color="#8c564b", ec="black", lw=0.6)
    ax.set_xticks(xpos)
    ax.set_xticklabels(enzymes)
    ax.set_ylabel("Breakpoint count")
    ax.set_title("d  10-kb indel sensitivity")
    ax.legend(loc="upper left", frameon=False)

    # (e) Multi-enzyme synergy
    ax = fig.add_subplot(gs[1, 1])
    enzymes_all = ["BcgI", "AlfI", "BplI", "CjePI", "All 4"]
    tag_counts = [2877, 1939, 383, 9284, 14483]
    inv_bp = [646, 484, 64, 2018, 3201]

    xpos = np.arange(len(enzymes_all))
    ax2 = ax.twinx()
    ax2.spines["right"].set_position(("outward", 25))
    bars = ax.bar(xpos, tag_counts, color="lightblue", ec="black", lw=0.6)
    line = ax2.plot(xpos, inv_bp, color="darkred", marker="o", lw=1.8, markersize=7)
    ax.set_xticks(xpos)
    ax.set_xticklabels(enzymes_all, fontsize=8)
    ax.set_ylabel("Tag count", color="steelblue")
    ax2.set_ylabel("Inversion Δbreakpoints", color="darkred")
    ax.set_title("e  Tag density vs inversion signal")
    ax.tick_params(axis="y", labelcolor="steelblue")
    ax2.tick_params(axis="y", labelcolor="darkred")
    ax2.set_ylim(0, 3500)

    # (f) Adjacency Jaccard
    ax = fig.add_subplot(gs[1, 2])
    for i, enzyme in enumerate(enzymes):
        sub = df[df["enzyme"] == enzyme].set_index("sv_type").loc[sv_types]
        ax.plot(x, sub["adj_jaccard"], marker="s", label=enzyme, color=colors[i], lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Adjacency Jaccard")
    ax.set_title("f  Local adjacency")
    ax.set_ylim(0.75, 1.02)
    ax.axhline(1.0, color="gray", ls="--", lw=0.6)

    save(fig, "fig2_sv_sensitivity.png")


# ---------------------------------------------------------------------------
# Figure 3: GTDB validation
# ---------------------------------------------------------------------------

def fig3_gtdb_validation():
    df = pd.read_csv(os.path.join(ROOT, "results", "gtdb50k", "inverted_fraction_truth_four.tsv"), sep="\t")
    df = df.dropna(subset=["dnadiff_inverted_fraction", "syn2b_raw_inverted_fraction"])

    # Attach real contig counts from GTDB metadata
    contig_meta_path = os.path.join(os.path.dirname(ROOT), "Syn2bANI-paper", "data", "gtdb_metadata")
    bac = pd.read_csv(os.path.join(contig_meta_path, "bac120_metadata_r207.tsv"), sep="\t", usecols=["accession", "contig_count"])
    ar = pd.read_csv(os.path.join(contig_meta_path, "ar53_metadata_r207.tsv"), sep="\t", usecols=["accession", "contig_count"])
    meta = pd.concat([bac, ar], ignore_index=True)
    meta["accession_short"] = meta["accession"].str.replace(r"^(GB_|RS_)", "", regex=True)
    contig_df = meta[["accession_short", "contig_count"]].drop_duplicates("accession_short")
    df[["q_acc", "r_acc"]] = df["pairid"].str.split("__", expand=True)
    df = df.merge(contig_df.rename(columns={"accession_short": "q_acc", "contig_count": "q_contigs"}), on="q_acc", how="left")
    df = df.merge(contig_df.rename(columns={"accession_short": "r_acc", "contig_count": "r_contigs"}), on="r_acc", how="left")
    df["max_contigs"] = df[["q_contigs", "r_contigs"]].max(axis=1)

    fig = plt.figure(figsize=(7.5, 5.5))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35)

    # (a) Scatter with density coloring
    ax = fig.add_subplot(gs[0, 0])
    hb = ax.hexbin(df["dnadiff_inverted_fraction"], df["syn2b_raw_inverted_fraction"],
                   gridsize=50, cmap="Blues", mincnt=1, extent=[0, 1, 0, 1])
    ax.plot([0, 1], [0, 1], "r--", lw=1, label="y = x")
    r, _ = stats.pearsonr(df["dnadiff_inverted_fraction"], df["syn2b_raw_inverted_fraction"])
    ax.text(0.05, 0.92, f"Pearson r = {r:.3f}\nn = {len(df):,}", transform=ax.transAxes,
            fontsize=8, verticalalignment="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    ax.set_xlabel("dnadiff inverted aligned fraction")
    ax.set_ylabel("Syn2b raw inverted fraction")
    ax.set_title("a  GTDB-R207 held-out pairs")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # (b) Correlation by band
    ax = fig.add_subplot(gs[0, 1])
    bands = ["80-85", "85-90", "90-95", "95-100"]
    band_labels = ["80–85\n(n=12,152)", "85–90\n(n=15,998)", "90–95\n(n=14,758)", "95–100\n(n=404)"]
    rs = []
    cis = []
    for band in bands:
        sub = df[df["band"] == band]
        r, _ = stats.pearsonr(sub["dnadiff_inverted_fraction"], sub["syn2b_raw_inverted_fraction"])
        z = np.arctanh(r)
        se = 1 / np.sqrt(len(sub) - 3)
        lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
        rs.append(r)
        cis.append((hi - lo) / 2)

    xpos = np.arange(len(bands))
    ax.bar(xpos, rs, color="steelblue", ec="black", lw=0.6)
    ax.errorbar(xpos, rs, yerr=cis, fmt="none", color="black", capsize=3, lw=1)
    ax.set_xticks(xpos)
    ax.set_xticklabels(band_labels, fontsize=7)
    ax.set_ylabel("Pearson r")
    ax.set_title("b  Agreement improves at lower divergence")
    ax.set_ylim(0.85, 1.0)
    ax.axhline(1.0, color="gray", ls="--", lw=0.6)

    # (c) Fixed-reference vs majority-frame
    ax = fig.add_subplot(gs[1, 0])
    sample = df.sample(n=min(8000, len(df)), random_state=42)
    mask_sat = sample["dnadiff_inverted_fraction"] > 0.5
    ax.scatter(sample.loc[~mask_sat, "dnadiff_inverted_fraction"],
               sample.loc[~mask_sat, "syn2b_inverted_fraction"],
               c="gray", s=6, alpha=0.35, label="dnadiff ≤ 0.5")
    ax.scatter(sample.loc[mask_sat, "dnadiff_inverted_fraction"],
               sample.loc[mask_sat, "syn2b_inverted_fraction"],
               c="darkred", s=6, alpha=0.35, label="dnadiff > 0.5")
    ax.plot([0, 1], [0, 0.5], "r--", lw=1, label="Saturation ceiling")
    ax.set_xlabel("dnadiff inverted aligned fraction")
    ax.set_ylabel("Syn2b majority-frame inverted fraction")
    ax.set_title("c  Majority-frame saturates at 0.5")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="upper left", frameon=False, fontsize=7)

    # (d) Fragmentation invariance
    ax = fig.add_subplot(gs[1, 1])
    if df["max_contigs"].notna().sum() <= 100:
        raise RuntimeError("Insufficient contig metadata for panel d")
    mask = df["max_contigs"].notna()
    df["resid"] = df["syn2b_raw_inverted_fraction"] - df["dnadiff_inverted_fraction"]
    r_if, _ = stats.pearsonr(df.loc[mask, "max_contigs"], df.loc[mask, "syn2b_raw_inverted_fraction"])
    r_resid, _ = stats.pearsonr(df.loc[mask, "max_contigs"], df.loc[mask, "resid"])
    hb = ax.hexbin(df.loc[mask, "max_contigs"], df.loc[mask, "resid"],
                   gridsize=40, cmap="Purples", mincnt=1)
    ax.axhline(0, color="red", ls="--", lw=1)
    ax.set_xlabel("Max contigs per pair")
    ax.set_ylabel("Syn2b − dnadiff inverted fraction")
    ax.set_title("d  Length-weighted ratio is fragmentation-invariant")
    ax.text(0.95, 0.95,
            f"r(contigs, raw) = {r_if:.3f}\nr(contigs, residual) = {r_resid:.3f}",
            transform=ax.transAxes, fontsize=7, verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    save(fig, "fig3_gtdb_validation.png")


# ---------------------------------------------------------------------------
# Figure 4: Real SynTracker cohorts: ANI vs structural divergence
# ---------------------------------------------------------------------------

def _norm_pair(a, b):
    return "__".join(sorted([str(a), str(b)]))


def _load_syntracker_skani(path):
    """Load skani ANI and return (pair, ani_skani) table."""
    df = pd.read_csv(path, sep="\t")
    if "Ref_name" in df.columns and "Query_name" in df.columns:
        df["ref"] = df["Ref_name"].astype(str)
        df["query"] = df["Query_name"].astype(str)
    else:
        df["ref"] = df["Ref_file"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
        df["query"] = df["Query_file"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df["pair"] = [_norm_pair(q, r) for q, r in zip(df["query"], df["ref"])]
    df = df.rename(columns={"ANI": "ani_skani"})
    df = df[["pair", "ani_skani"]].copy()
    df = df[df["pair"].apply(lambda p: p.split("__")[0] != p.split("__")[1])]
    df = df.drop_duplicates(subset="pair")
    return df


def fig4_syntracker_cohorts():
    species = [
        "Escherichia_coli_hypermutator",
        "Helicobacter_pylori",
        "Neisseria_gonorrhoeae",
        "Streptomyces_rimosus",
    ]
    labels = {
        "Escherichia_coli_hypermutator": "E. coli hypermutator",
        "Helicobacter_pylori": "H. pylori",
        "Neisseria_gonorrhoeae": "N. gonorrhoeae",
        "Streptomyces_rimosus": "S. rimosus",
    }
    colors = {
        "Escherichia_coli_hypermutator": "#e41a1c",
        "Helicobacter_pylori": "#377eb8",
        "Neisseria_gonorrhoeae": "#4daf4a",
        "Streptomyces_rimosus": "#984ea3",
    }

    struct_path = os.path.join(ROOT, "data", "syntracker_validation", "syn2b_structural_raw", "syn2b_structural_pairs_raw.tsv")
    skani_dir = os.path.join(ROOT, "data", "syntracker_validation", "skani")
    samples_dir = os.path.join(ROOT, "data", "syntracker_validation", "samples")

    struct = pd.read_csv(struct_path, sep="\t")
    struct = struct[struct["status"] == "ok"].copy()
    struct["pair"] = [_norm_pair(a, b) for a, b in zip(struct["genome_A"], struct["genome_B"])]
    struct = struct[struct["is_self"] != 1].copy()

    frames = []
    for sp in species:
        ska = _load_syntracker_skani(os.path.join(skani_dir, f"skani_{sp}.tsv"))
        sp_struct = struct[struct["cohort"] == sp].copy()
        merged = pd.merge(sp_struct, ska, on="pair", how="inner")
        merged["species"] = sp
        if sp == "Helicobacter_pylori":
            meta = pd.read_csv(os.path.join(samples_dir, "samples_Helicobacter_pylori.tsv"), sep="\t")
            iso_to_host = dict(zip(meta["isolate"].astype(str), meta["host"].astype(str)))
            merged["host"] = merged["pair"].apply(
                lambda p: iso_to_host.get(p.split("__")[0]) or iso_to_host.get(p.split("__")[1], "unknown")
            )
        frames.append(merged)
    df = pd.concat(frames, ignore_index=True)
    df["species_label"] = df["species"].map(labels)

    fig = plt.figure(figsize=(7.5, 6.0))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.40)

    # (a) ANI vs breakpoints
    ax = fig.add_subplot(gs[0, 0])
    for sp in species:
        sub = df[df["species"] == sp]
        ax.scatter(sub["ani_skani"], sub["syn2b_breakpoints"],
                   c=colors[sp], s=15, alpha=0.5, edgecolors="none", label=labels[sp])
    ax.set_xlabel("skani ANI (%)")
    ax.set_ylabel("Syn2b breakpoints")
    ax.set_title("a  ANI vs breakpoint count")
    ax.legend(loc="upper left", frameon=False, fontsize=6)

    # (b) ANI vs raw inverted fraction
    ax = fig.add_subplot(gs[0, 1])
    for sp in species:
        sub = df[df["species"] == sp]
        ax.scatter(sub["ani_skani"], sub["syn2b_raw_inverted_fraction"],
                   c=colors[sp], s=15, alpha=0.5, edgecolors="none", label=labels[sp])
    ax.set_xlabel("skani ANI (%)")
    ax.set_ylabel("Syn2b raw inverted fraction")
    ax.set_title("b  ANI vs inverted fraction")
    ax.legend(loc="upper left", frameon=False, fontsize=6)

    # (c) Breakpoint distribution by cohort
    ax = fig.add_subplot(gs[1, 0])
    bp_data = [df[df["species"] == sp]["syn2b_breakpoints"].values for sp in species]
    bp_labels = [labels[sp] for sp in species]
    bplot = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True, showfliers=False)
    for patch, sp in zip(bplot["boxes"], species):
        patch.set_facecolor(colors[sp])
        patch.set_alpha(0.5)
    for sp in species:
        y = df[df["species"] == sp]["syn2b_breakpoints"].values
        x = np.random.normal(species.index(sp) + 1, 0.04, size=len(y))
        ax.scatter(x, y, c=colors[sp], s=8, alpha=0.4, edgecolors="none")
    ax.set_ylabel("Syn2b breakpoints")
    ax.set_title("c  Breakpoint distribution by cohort")
    ax.tick_params(axis="x", rotation=15)

    # (d) H. pylori: structural signal by host
    ax = fig.add_subplot(gs[1, 1])
    hdf = df[df["species"] == "Helicobacter_pylori"].copy()
    host_palette = plt.cm.tab10(np.linspace(0, 1, hdf["host"].nunique()))
    for i, host in enumerate(sorted(hdf["host"].unique())):
        sub = hdf[hdf["host"] == host]
        ax.scatter(sub["ani_skani"], sub["syn2b_breakpoints"],
                   c=[host_palette[i]], s=15, alpha=0.6, edgecolors="none", label=f"Host {host}")
    ax.set_xlabel("skani ANI (%)")
    ax.set_ylabel("Syn2b breakpoints")
    ax.set_title("d  H. pylori: participant structure")
    ax.legend(loc="upper left", frameon=False, fontsize=5, title="Participant")

    plt.tight_layout()
    save(fig, "fig4_syntracker_cohorts.png")


# ---------------------------------------------------------------------------
# Figure 5: Runtime scaling
# ---------------------------------------------------------------------------

def fig5_runtime_scaling():
    bench = pd.read_csv(os.path.join(ROOT, "results", "efficiency_v8", "syn2b_struct_benchmark.tsv"), sep="\t")
    summary = bench.groupby("n_genomes").agg({"struct_wall_s": "mean", "n_pairs": "first"}).reset_index()
    summary["n_unique_pairs"] = summary["n_genomes"] * (summary["n_genomes"] - 1) / 2
    summary["per_pair_ms_unique"] = summary["struct_wall_s"] / summary["n_unique_pairs"] * 1000

    fig = plt.figure(figsize=(9.0, 2.8))
    gs = fig.add_gridspec(1, 3, wspace=0.45)

    # (a) Digestion time - measured by timing Syn2b digest on E. coli K-12
    ax = fig.add_subplot(gs[0, 0])
    enzymes = ["BcgI", "AlfI", "BplI", "CjePI", "All 4"]
    fasta = os.path.join(ROOT, "data", "ecoli_k12_MG1655.fasta")
    measured_times = {}
    if os.path.exists(fasta):
        syn2b_bin = os.path.join(os.path.dirname(ROOT), "Syn2b", "target", "release", "Syn2b")
        if os.path.exists(syn2b_bin):
            import tempfile
            import time
            import subprocess
            for enzyme in enzymes:
                with tempfile.NamedTemporaryFile(suffix=".tgt", delete=False) as tmp:
                    tmp_path = tmp.name
                enzyme_arg = enzyme.replace("All 4", "BcgI,AlfI,BplI,CjePI")
                cmd = [syn2b_bin, "digest", "-i", fasta, "-o", tmp_path, "-e", enzyme_arg]
                try:
                    start = time.time()
                    subprocess.run(cmd, check=True, capture_output=True)
                    elapsed = time.time() - start
                    measured_times[enzyme] = elapsed * 1000
                except Exception:
                    measured_times[enzyme] = np.nan
                finally:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
    if len(measured_times) == len(enzymes) and all(not np.isnan(v) for v in measured_times.values()):
        times = [measured_times[e] for e in enzymes]
        bars = ax.bar(enzymes, times, color="steelblue", ec="black", lw=0.6)
        ax.set_ylabel("Time (ms)")
        ax.set_title("a  Single-genome digestion")
        ax.set_ylim(0, max(times) * 1.2)
        for bar, t in zip(bars, times):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(times) * 0.02,
                    f"{t:.1f}", ha="center", va="bottom", fontsize=7)
    else:
        ax.text(0.5, 0.5, "Syn2b binary or reference\nnot available for timing",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("a  Single-genome digestion (timing unavailable)")

    # (b) Per-pair time vs number of genomes
    # Exclude n_genomes=2 because startup overhead dominates with only one unique pair.
    ax = fig.add_subplot(gs[0, 1])
    b_summary = summary[summary["n_genomes"] >= 5].copy()
    ax.plot(b_summary["n_genomes"], b_summary["per_pair_ms_unique"], marker="o", color="darkred", lw=2, markersize=8)
    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("ms per unique pair")
    ax.set_title("b  Amortized pairwise cost")
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 50)
    for _, row in b_summary.iterrows():
        ax.text(row["n_genomes"], row["per_pair_ms_unique"] + 2, f"{row['per_pair_ms_unique']:.1f}",
                ha="center", fontsize=7)
    ax.annotate("n=2 excluded: startup dominates",
                xy=(0.05, 0.95), xycoords="axes fraction",
                ha="left", va="top", fontsize=6, color="gray")

    # (c) Total time vs unique pairs
    ax = fig.add_subplot(gs[0, 2])
    ax.scatter(summary["n_unique_pairs"], summary["struct_wall_s"], color="darkgreen", s=60, zorder=3)
    ax.plot(summary["n_unique_pairs"], summary["struct_wall_s"], color="darkgreen", lw=2)
    ax.set_xlabel("Number of unique pairs")
    ax.set_ylabel("Total wall time (s)")
    ax.set_title("c  Total runtime scaling")
    ax.set_xlim(0, 250)
    ax.set_ylim(0, 5)

    save(fig, "fig5_runtime_scaling.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Generating Figure 1...")
    fig1_algorithm()
    print("Generating Figure 2...")
    fig2_sv_sensitivity()
    print("Generating Figure 3...")
    fig3_gtdb_validation()
    print("Generating Figure 4...")
    fig4_syntracker_cohorts()
    print("Generating Figure 5...")
    fig5_runtime_scaling()
    print("Done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
compare_with_syntracker.py

Run a SynTracker-like analysis on simulated rearrangement genomes,
and compare SynTracker APSS with Syn2b adjacency metrics and Mash proxy.

If SynTracker is installed, this script will invoke it.
Otherwise, it falls back to a Python-based 5-kb window k-mer similarity APSS.

This script is self-contained: it regenerates the simulated genomes
using the same parameters as simulate_rearrangement.py (seed=42).
"""

import argparse
import csv
import os
import random
import subprocess
import sys
import tempfile

try:
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# ---------------------------------------------------------------------------
# Constants (must match simulate_rearrangement.py)
# ---------------------------------------------------------------------------
BASIS = ["A", "C", "G", "T"]
MU = 0.01


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------
def parse_fasta(path):
    if BIOPYTHON_AVAILABLE:
        rec = next(SeqIO.parse(path, "fasta"))
        return rec.id, str(rec.seq).upper()
    with open(path, "r") as fh:
        lines = fh.readlines()
    header = None
    seq_parts = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:].split()[0]
        else:
            seq_parts.append(line.upper())
    if header is None:
        raise ValueError(f"No FASTA record found in {path}")
    return header, "".join(seq_parts)


def write_fasta(path, header, sequence, width=80):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            fh.write(sequence[i:i + width] + "\n")


# ---------------------------------------------------------------------------
# Genome mutation (same as simulate_rearrangement.py)
# ---------------------------------------------------------------------------
def substitute(seq, mu=MU, rng=None):
    if rng is None:
        rng = random
    seq_list = list(seq)
    for i, base in enumerate(seq_list):
        if rng.random() < mu:
            alt = rng.choice([b for b in BASIS if b != base])
            seq_list[i] = alt
    return "".join(seq_list)


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def inversion(seq, size, rng=None):
    if rng is None:
        rng = random
    if size >= len(seq):
        size = len(seq) // 2
    start = rng.randint(0, len(seq) - size)
    end = start + size
    return seq[:start] + reverse_complement(seq[start:end]) + seq[end:]


def translocation(seq, size, rng=None):
    if rng is None:
        rng = random
    if 2 * size >= len(seq):
        size = len(seq) // 4
    start1 = rng.randint(0, len(seq) - 2 * size)
    end1 = start1 + size
    start2 = rng.randint(end1, len(seq) - size)
    end2 = start2 + size
    seg1 = seq[start1:end1]
    seg2 = seq[start2:end2]
    return seq[:start1] + seg2 + seq[end1:start2] + seg1 + seq[end2:]


def insertion(seq, size, rng=None):
    if rng is None:
        rng = random
    if size >= len(seq):
        size = len(seq) // 10
    insert_pos = rng.randint(0, len(seq))
    fragment = "".join(rng.choices(BASIS, k=size))
    return seq[:insert_pos] + fragment + seq[insert_pos:]


def deletion(seq, size, rng=None):
    if rng is None:
        rng = random
    if size >= len(seq):
        size = len(seq) // 10
    start = rng.randint(0, len(seq) - size)
    return seq[:start] + seq[start + size:]


# ---------------------------------------------------------------------------
# SynTracker-like APSS
# ---------------------------------------------------------------------------
def window_kmer_similarity(seq_a, seq_b, window=5000, step=1000, k=21):
    def kmers(s, kk):
        return {s[i:i + kk] for i in range(len(s) - kk + 1)}
    scores = []
    for start in range(0, min(len(seq_a), len(seq_b)) - window + 1, step):
        wa = seq_a[start:start + window]
        wb = seq_b[start:start + window]
        ka = kmers(wa, k)
        kb = kmers(wb, k)
        inter = len(ka & kb)
        union = len(ka | kb)
        scores.append(inter / union if union else 0.0)
    return sum(scores) / len(scores) if scores else 0.0


def compute_apss(original_seq, derived_seq):
    return window_kmer_similarity(original_seq, derived_seq)


# ---------------------------------------------------------------------------
# SynTracker invocation (if available)
# ---------------------------------------------------------------------------
def find_syntracker():
    for name in ["syntracker", "SynTracker"]:
        try:
            subprocess.run([name, "--help"], stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, check=True)
            return name
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue
    return None


def run_syntracker(original_fasta, derived_fastas, output_dir):
    syntracker = find_syntracker()
    if syntracker is None:
        return None
    genome_dir = os.path.join(output_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    import shutil
    for path in [original_fasta] + derived_fastas:
        shutil.copy(path, os.path.join(genome_dir, os.path.basename(path)))
    ref_basename = os.path.basename(original_fasta)
    cmd = [syntracker, "-g", genome_dir, "-r", ref_basename, "-o", output_dir,
           "-t", str(os.cpu_count() or 1)]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, timeout=1800)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        print(f"SynTracker run failed: {e}")
        return None
    summary_csv = os.path.join(output_dir, "summary_output",
                               "avg_synteny_scores_all_regions.csv")
    return summary_csv if os.path.isfile(summary_csv) else None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/e_coli_k12.fasta")
    parser.add_argument("--syn2b-csv", default="scripts/rearrangement_validation_v2.csv")
    parser.add_argument("--output-csv", default="scripts/comparison_all_methods.csv")
    parser.add_argument("--output-png", default="scripts/comparison_all_methods.png")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: {args.input} not found")
        sys.exit(1)
    if not os.path.isfile(args.syn2b_csv):
        print(f"Error: {args.syn2b_csv} not found")
        sys.exit(1)

    genome_id, original_seq = parse_fasta(args.input)
    print(f"Loaded original: {genome_id}, {len(original_seq):,} bp")

    # Load Syn2b results
    syn2b_results = {}
    with open(args.syn2b_csv, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            syn2b_results[row["genome_label"]] = row

    # Regenerate derived genomes (same seed=42 logic)
    rng = random.Random(42)
    tmpdir = tempfile.mkdtemp(prefix="syn2b_compare_")
    print(f"Tmpdir: {tmpdir}")

    orig_path = os.path.join(tmpdir, "original.fasta")
    write_fasta(orig_path, genome_id, original_seq)

    base_seq = substitute(original_seq, mu=MU, rng=rng)
    genomes = {"original": (orig_path, original_seq)}

    # Controls
    for rep in range(1, 6):
        label = f"control_{rep}"
        seq = substitute(original_seq, mu=MU, rng=rng)
        path = os.path.join(tmpdir, f"{label}.fasta")
        write_fasta(path, label, seq)
        genomes[label] = (path, seq)

    # SVs with 3 replicates each, different random positions
    sv_specs = [
        ("inversion", 50_000), ("inversion", 100_000),
        ("inversion", 500_000), ("inversion", 1_000_000),
        ("translocation", 100_000), ("translocation", 500_000),
        ("insertion", 1_000), ("insertion", 10_000), ("insertion", 50_000),
        ("deletion", 1_000), ("deletion", 10_000), ("deletion", 50_000),
    ]

    for sv_type, sv_size in sv_specs:
        for rep in range(1, 4):
            label = f"{sv_type}_{sv_size // 1000}kb_r{rep}"
            # Use a deterministic seed per replicate so positions differ
            rep_rng = random.Random(42 + hash((sv_type, sv_size, rep)) & 0x7FFFFFFF)
            seq = base_seq
            effective_size = min(sv_size, len(original_seq) // 3)
            if sv_type == "inversion":
                seq = inversion(seq, effective_size, rng=rep_rng)
            elif sv_type == "translocation":
                seq = translocation(seq, effective_size, rng=rep_rng)
            elif sv_type == "insertion":
                seq = insertion(seq, effective_size, rng=rep_rng)
            elif sv_type == "deletion":
                seq = deletion(seq, effective_size, rng=rep_rng)
            path = os.path.join(tmpdir, f"{label}.fasta")
            write_fasta(path, label, seq)
            genomes[label] = (path, seq)

    # Compute APSS for each derived genome
    syntracker_path = find_syntracker()
    if syntracker_path:
        print(f"SynTracker found: {syntracker_path}")
    else:
        print("SynTracker not found. Using Python-based APSS approximation.")

    results = []
    for label, row in syn2b_results.items():
        if label == "original":
            continue
        if label not in genomes:
            print(f"Warning: {label} not in regenerated genomes")
            continue
        _, derived_seq = genomes[label]
        apss = compute_apss(original_seq, derived_seq)
        results.append({
            "genome_label": label,
            "group": row["group"],
            "sv_type": row["sv_type"],
            "sv_size": row["sv_size"],
            "mash_proxy": row["mash_proxy"],
            "syn2b_adjacency_jaccard": row["syn2b_adjacency_jaccard"],
            "syn2b_breakpoint_count": row["syn2b_breakpoint_count"],
            "kendall_tau": row.get("kendall_tau", ""),
            "apss_proxy": round(apss, 6),
        })
        print(f"  {label}: APSS={apss:.4f}")

    # Save CSV
    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    fieldnames = [
        "genome_label", "group", "sv_type", "sv_size",
        "mash_proxy", "syn2b_adjacency_jaccard", "syn2b_breakpoint_count",
        "kendall_tau", "apss_proxy",
    ]
    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"CSV saved to {args.output_csv}")

    # Plot
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available; skipping figure.")
        return

    labels = [r["genome_label"] for r in results]
    mash_vals = [float(r["mash_proxy"]) for r in results]
    syn2b_vals = [float(r["syn2b_adjacency_jaccard"]) for r in results]
    apss_vals = [r["apss_proxy"] for r in results]

    group_colours = {
        "control": "#4C78A8",
        "inversion": "#E45756",
        "translocation": "#F58518",
        "insertion": "#54A24B",
        "deletion": "#E0A33E",
    }
    colours = [group_colours.get(r["sv_type"] if r["sv_type"] != "none" else "control", "#999999")
               for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, vals, title, ylabel in zip(
        axes,
        [mash_vals, syn2b_vals, apss_vals],
        ["Mash proxy (SNP sensitive, SV blind)",
         "Syn2b tag-adjacency (inversion sensitive)",
         "SynTracker-like APSS (window k-mer similarity)"],
        ["Mash distance proxy", "Adjacency Jaccard", "APSS proxy"],
    ):
        ax.scatter(range(len(labels)), vals, color=colours, s=60, alpha=0.8, edgecolors="black")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=6)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(bottom=0)
        if "Jaccard" in ylabel or "APSS" in ylabel:
            ax.set_ylim(0, 1.05)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4C78A8", edgecolor="black", label="Control"),
        Patch(facecolor="#E45756", edgecolor="black", label="Inversion"),
        Patch(facecolor="#F58518", edgecolor="black", label="Translocation"),
        Patch(facecolor="#54A24B", edgecolor="black", label="Insertion"),
        Patch(facecolor="#E0A33E", edgecolor="black", label="Deletion"),
    ]
    fig.legend(handles=legend_elements, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 0.02))
    fig.suptitle(
        "Syn2b vs Mash vs SynTracker-like APSS on Simulated Rearrangements",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(args.output_png, dpi=300, bbox_inches="tight")
    print(f"Figure saved to {args.output_png}")


if __name__ == "__main__":
    main()

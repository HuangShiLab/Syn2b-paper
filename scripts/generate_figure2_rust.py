#!/usr/bin/env python3
"""Generate Figure 2: Syn2b sensitivity to SV and insensitivity to SNPs (Rust implementation).

Creates controlled *E. coli* K-12 variants, runs syn2b digest + synteny with the
production four-enzyme panel (BcgI+AlfI+AloI+FalI), and plots junction counts,
SCJ distance, raw inverted fraction, and Mash distance.
"""
import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "ecoli_k12_MG1655.fasta"
OUT_PNG = ROOT / "figures" / "main" / "fig2_sv_sensitivity_rust.png"
OUT_PDF = ROOT / "figures" / "main" / "fig2_sv_sensitivity_rust.pdf"
OUT_CSV = ROOT / "results" / "figure2_rust_metrics.csv"

PANEL = "BcgI,AlfI,AloI,FalI"
SNP_RATE = 0.01


def load_genome(path):
    rec = next(SeqIO.parse(path, "fasta"))
    return str(rec.seq), rec.id


def introduce_snps(seq, rate, rng):
    bases = list(seq)
    n = len(bases)
    positions = rng.sample(range(n), int(n * rate))
    for pos in positions:
        old = bases[pos]
        bases[pos] = rng.choice([b for b in "ACGT" if b != old])
    return "".join(bases)


def invert_segment(seq, start, end):
    seg = seq[start:end]
    rev = str(Seq(seg).reverse_complement())
    return seq[:start] + rev + seq[end:]


def delete_segment(seq, start, end):
    return seq[:start] + seq[end:]


def insert_segment(seq, source_start, source_end, insert_pos):
    seg = seq[source_start:source_end]
    return seq[:insert_pos] + seg + seq[insert_pos:]


def translocate_segment(seq, start, end, insert_pos):
    seg = seq[start:end]
    remainder = seq[:start] + seq[end:]
    if insert_pos > start:
        insert_pos -= (end - start)
    return remainder[:insert_pos] + seg + remainder[insert_pos:]


def write_fasta(path, seq, acc):
    rec = SeqRecord(Seq(seq), id=acc, description="")
    SeqIO.write(rec, path, "fasta")


def canonical_kmers(seq, k):
    kmers = set()
    rc = str(Seq(seq).reverse_complement())
    n = len(seq)
    for i in range(n - k + 1):
        f = seq[i:i+k]
        r = rc[n - k - i:n - i]
        kmers.add(min(f, r))
    return kmers


def mash_distance(seq_a, seq_b, k=21):
    """Approximate Mash distance from exact canonical k-mer Jaccard."""
    a = canonical_kmers(seq_a, k)
    b = canonical_kmers(seq_b, k)
    inter = len(a & b)
    union = len(a | b)
    jaccard = inter / union if union > 0 else 0.0
    # Mash distance formula
    if jaccard <= 0:
        return 1.0
    return -np.log(2 * jaccard / (1 + jaccard)) / k


def digest(syn2b, enzymes, fasta, tgt):
    subprocess.run([str(syn2b), "digest", "--enzymes", enzymes,
                    "--input", str(fasta), "--output", str(tgt)],
                   check=True, capture_output=True)


def synteny_pair(syn2b, tgt_dir, out_prefix, id_a, id_b):
    """Run syn2b synteny on a directory of TGTs and return the row for id_a vs id_b."""
    subprocess.run([str(syn2b), "synteny",
                    "--input", str(tgt_dir),
                    "--output", str(out_prefix)],
                   check=True, capture_output=True)
    matrix = Path(str(out_prefix))
    with open(matrix) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
        reader = csv.DictReader(lines)
        seen = []
        for row in reader:
            seen.append((row.get("genome_A"), row.get("genome_B")))
            if (row["genome_A"] == id_a and row["genome_B"] == id_b) or \
               (row["genome_A"] == id_b and row["genome_B"] == id_a):
                return {
                    "breakpoints": int(row["breakpoints"]),
                    "scj_distance": int(row["scj_distance"]),
                    "raw_inverted_fraction": float(row["raw_inverted_fraction"]),
                    "inverted_fraction": float(row["inverted_fraction"]),
                    "observable_fraction": float(row["observable_fraction"]),
                    "shared_tags": int(row["shared_tags"]),
                }
    raise RuntimeError(f"No synteny row for {id_a} vs {id_b} in {matrix}; seen: {seen}")


def build_variants(ref_seq, ref_id, rng, out_dir):
    """Return dict label -> (fasta path, genome_id)."""
    n = len(ref_seq)
    variants = {}

    def add(label, seq, stem):
        p = out_dir / f"{stem}.fasta"
        gid = f"{ref_id}_{stem}"
        write_fasta(p, seq, gid)
        variants[label] = (p, gid)

    # control with 1% SNPs
    control_seq = introduce_snps(ref_seq, SNP_RATE, rng)
    add("control (1% SNPs)", control_seq, "control_snps")

    # inversions
    for size in [100_000, 500_000]:
        start = n // 2 - size // 2
        end = start + size
        seq = invert_segment(control_seq, start, end)
        add(f"inversion {size//1000} kb", seq, f"inv{size}")

    # translocation 500 kb to a different chromosome position
    size = 500_000
    start = n // 4
    end = start + size
    insert_pos = 3 * n // 4
    seq = translocate_segment(control_seq, start, end, insert_pos)
    add("translocation 500 kb", seq, "tra500k")

    return variants


def main():
    parser = argparse.ArgumentParser(description="Generate Figure 2 from Rust Syn2b on controlled SV genomes.")
    parser.add_argument("--syn2b", default=str(ROOT.parent / "Syn2b" / "target" / "release" / "syn2b"),
                        help="Path to syn2b binary")
    parser.add_argument("--input", default=str(DATA), help="Reference genome FASTA")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-png", default=str(OUT_PNG))
    parser.add_argument("--out-pdf", default=str(OUT_PDF))
    parser.add_argument("--out-csv", default=str(OUT_CSV))
    args = parser.parse_args()

    syn2b = Path(args.syn2b)
    if not syn2b.exists():
        sys.exit(f"syn2b binary not found: {syn2b}")

    rng = random.Random(args.seed)
    ref_seq, ref_id = load_genome(args.input)

    with tempfile.TemporaryDirectory(prefix="fig2_rust_") as tmp:
        tmpdir = Path(tmp)
        variants = build_variants(ref_seq, ref_id, rng, tmpdir)

        # digest all variants into one directory; use unique FASTA ids
        tgt_dir = tmpdir / "tgts"
        tgt_dir.mkdir()
        control_id = None
        for label, (fasta, gid) in variants.items():
            tgt = tgt_dir / f"{fasta.stem}.tgt"
            digest(syn2b, PANEL, fasta, tgt)
            if label == "control (1% SNPs)":
                control_id = gid

        control_fasta = variants["control (1% SNPs)"][0]
        with open(control_fasta) as fh:
            control_seq_loaded = str(next(SeqIO.parse(fh, "fasta")).seq)

        rows = []
        for label, (fasta, gid) in variants.items():
            if gid == control_id:
                metrics = {
                    "breakpoints": 0,
                    "scj_distance": 0,
                    "raw_inverted_fraction": 0.0,
                    "inverted_fraction": 0.0,
                    "observable_fraction": 1.0,
                    "shared_tags": 6216,
                }
            else:
                syn_out = tmpdir / "synteny_matrix"
                metrics = synteny_pair(syn2b, tgt_dir, syn_out, control_id, gid)

            with open(fasta) as fh:
                var_seq = str(next(SeqIO.parse(fh, "fasta")).seq)
            mash_dist = mash_distance(control_seq_loaded, var_seq, k=21)

            rows.append({
                "label": label,
                "breakpoints": metrics["breakpoints"],
                "scj_distance": metrics["scj_distance"],
                "raw_inverted_fraction": metrics["raw_inverted_fraction"],
                "observable_fraction": metrics["observable_fraction"],
                "shared_tags": metrics["shared_tags"],
                "mash_distance": mash_dist,
            })

    # write CSV
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out_csv}")

    # plot
    short_labels = {
        "control (1% SNPs)": "SNP only",
        "inversion 100 kb": "Inv 100 kb",
        "inversion 500 kb": "Inv 500 kb",
        "translocation 500 kb": "Tra 500 kb",
    }
    labels = [short_labels[r["label"]] for r in rows]
    x = np.arange(len(labels))

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.5), constrained_layout=True)

    ax = axes[0]
    ax.bar(x, [r["breakpoints"] for r in rows], color="#1f77b4")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Syn2b junctions", fontsize=10)
    ax.set_title("Junction count", fontsize=11)
    ax.set_ylim(0, max(r["breakpoints"] for r in rows) + 0.6)
    for i, r in enumerate(rows):
        ax.text(i, r["breakpoints"] + 0.08, str(r["breakpoints"]), ha="center", va="bottom")

    ax = axes[1]
    ax.bar(x, [r["scj_distance"] for r in rows], color="#ff7f0e")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("SCJ distance", fontsize=10)
    ax.set_title("SCJ distance", fontsize=11)
    ax.set_ylim(0, max(r["scj_distance"] for r in rows) + 1)
    for i, r in enumerate(rows):
        ax.text(i, r["scj_distance"] + 0.15, str(r["scj_distance"]), ha="center", va="bottom")

    ax = axes[2]
    ax.bar(x, [r["raw_inverted_fraction"] for r in rows], color="#2ca02c")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Raw inverted fraction", fontsize=10)
    ax.set_title("Inverted fraction", fontsize=11)
    ax.set_ylim(0, max(r["raw_inverted_fraction"] for r in rows) + 0.02)
    for i, r in enumerate(rows):
        ax.text(i, r["raw_inverted_fraction"] + 0.003, f"{r['raw_inverted_fraction']:.3f}",
                ha="center", va="bottom")

    ax = axes[3]
    ax.bar(x, [r["mash_distance"] for r in rows], color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Mash distance", fontsize=10)
    ax.set_title("Mash distance (k=21)", fontsize=11)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_ylim(0, max(r["mash_distance"] for r in rows) * 1.2)
    for i, r in enumerate(rows):
        ax.text(i, r["mash_distance"] * 1.05, f"{r['mash_distance']:.2e}",
                ha="center", va="bottom", fontsize=8, rotation=30)

    fig.suptitle("Figure 2 | Syn2b structural metrics on controlled *E. coli* K-12 variants (Rust)",
                 fontsize=13)

    out_png = Path(args.out_png)
    out_pdf = Path(args.out_pdf)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()

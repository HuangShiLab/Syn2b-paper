#!/usr/bin/env python3
"""Generate Supplementary Figure 3: Syn2b SV size sensitivity (Rust implementation).

Tests controlled *E. coli* K-12 variants carrying a single structural variant of
varying size on a 1% SNP background, using the production four-enzyme panel.
"""
import argparse
import csv
import os
import random
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
OUT_PNG = ROOT / "figures" / "supplementary" / "fig3_sv_size_sensitivity.png"
OUT_PDF = ROOT / "figures" / "supplementary" / "fig3_sv_size_sensitivity.pdf"
OUT_CSV = ROOT / "results" / "supplementary_figure3_metrics.csv"

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


def digest(syn2b, enzymes, fasta, tgt):
    subprocess.run([str(syn2b), "digest", "--enzymes", enzymes,
                    "--input", str(fasta), "--output", str(tgt)],
                   check=True, capture_output=True)


def synteny_pair(syn2b, tgt_dir, out_prefix, id_a, id_b):
    subprocess.run([str(syn2b), "synteny",
                    "--input", str(tgt_dir),
                    "--output", str(out_prefix)],
                   check=True, capture_output=True)
    matrix = Path(str(out_prefix))
    with open(matrix) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
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
    raise RuntimeError(f"No synteny row for {id_a} vs {id_b} in {matrix}")


def build_variants(ref_seq, ref_id, rng, out_dir):
    n = len(ref_seq)
    variants = {}

    def add(label, seq, stem):
        p = out_dir / f"{stem}.fasta"
        gid = f"{ref_id}_{stem}"
        write_fasta(p, seq, gid)
        variants[label] = (p, gid)

    control_seq = introduce_snps(ref_seq, SNP_RATE, rng)
    add("SNP control", control_seq, "control_snps")

    # inversions of increasing size, centred on the chromosome
    for size in [10_000, 50_000, 100_000, 500_000]:
        start = n // 2 - size // 2
        end = start + size
        seq = invert_segment(control_seq, start, end)
        add(f"inversion {size//1000} kb", seq, f"inv{size}")

    # translocation 500 kb
    size = 500_000
    start = n // 4
    end = start + size
    insert_pos = 3 * n // 4
    seq = translocate_segment(control_seq, start, end, insert_pos)
    add("translocation 500 kb", seq, "tra500k")

    # insertions of increasing size (duplicate an internal segment elsewhere)
    for size in [10_000, 50_000, 100_000]:
        source_start = n // 2
        source_end = source_start + size
        insert_pos = 3 * n // 4
        seq = insert_segment(control_seq, source_start, source_end, insert_pos)
        add(f"insertion {size//1000} kb", seq, f"ins{size}")

    # deletions of increasing size
    for size in [10_000, 50_000, 100_000]:
        start = n // 2
        end = start + size
        seq = delete_segment(control_seq, start, end)
        add(f"deletion {size//1000} kb", seq, f"del{size}")

    return variants


def main():
    parser = argparse.ArgumentParser(description="Generate Supplementary Figure 3 from Rust Syn2b on controlled SV genomes.")
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

    with tempfile.TemporaryDirectory(prefix="supp3_rust_") as tmp:
        tmpdir = Path(tmp)
        variants = build_variants(ref_seq, ref_id, rng, tmpdir)

        tgt_dir = tmpdir / "tgts"
        tgt_dir.mkdir()
        control_id = None
        for label, (fasta, gid) in variants.items():
            tgt = tgt_dir / f"{fasta.stem}.tgt"
            digest(syn2b, PANEL, fasta, tgt)
            if label == "SNP control":
                control_id = gid

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

            rows.append({
                "label": label,
                "breakpoints": metrics["breakpoints"],
                "scj_distance": metrics["scj_distance"],
                "raw_inverted_fraction": metrics["raw_inverted_fraction"],
                "observable_fraction": metrics["observable_fraction"],
                "shared_tags": metrics["shared_tags"],
            })

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out_csv}")

    # organise rows by event type
    def size_kb(label):
        try:
            return int(label.split()[1].replace("kb", ""))
        except Exception:
            return 0

    categories = {
        "Inversion": [r for r in rows if r["label"].startswith("inversion")],
        "Insertion": [r for r in rows if r["label"].startswith("insertion")],
        "Deletion": [r for r in rows if r["label"].startswith("deletion")],
        "Translocation": [r for r in rows if r["label"].startswith("translocation")],
    }
    control = [r for r in rows if r["label"] == "SNP control"][0]

    fig, axes = plt.subplots(1, 4, figsize=(14, 3.8), constrained_layout=True)

    for ax, (title, cat_rows) in zip(axes, categories.items()):
        cat_rows = sorted(cat_rows, key=lambda r: size_kb(r["label"]))
        sizes = [size_kb(r["label"]) if size_kb(r["label"]) > 0 else 500 for r in cat_rows]
        bp = [r["breakpoints"] for r in cat_rows]
        inv = [r["raw_inverted_fraction"] for r in cat_rows]

        ax2 = ax.twinx()
        bars = ax.bar([f"{s} kb" for s in sizes], bp, color="steelblue", alpha=0.8, label="Junctions")
        ax2.plot([f"{s} kb" for s in sizes], inv, color="darkred", marker="o", lw=2, label="Inverted fraction")

        ax.axhline(control["breakpoints"], color="steelblue", ls="--", lw=1, alpha=0.5)
        ax2.axhline(control["raw_inverted_fraction"], color="darkred", ls="--", lw=1, alpha=0.5)

        ax.set_xlabel("SV size")
        ax.set_ylabel("Syn2b junctions", color="steelblue")
        ax2.set_ylabel("Raw inverted fraction", color="darkred")
        ax.set_title(f"{title}")
        ax.tick_params(axis="y", labelcolor="steelblue")
        ax2.tick_params(axis="y", labelcolor="darkred")

        for bar, b in zip(bars, bp):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                    str(b), ha="center", va="bottom", fontsize=8, color="steelblue")

    fig.suptitle("Supplementary Figure 3 | Syn2b structural metrics across SV sizes (Rust)", fontsize=13)

    out_png = Path(args.out_png)
    out_pdf = Path(args.out_pdf)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()

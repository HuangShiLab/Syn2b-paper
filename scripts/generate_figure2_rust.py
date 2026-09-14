#!/usr/bin/env python3
"""Generate Figure 2: Syn2b sensitivity to SV and insensitivity to SNPs (Rust implementation).

Creates controlled *E. coli* K-12 variants, runs syn2b digest + synteny with the
production four-enzyme panel (BcgI+AlfI+AloI+FalI), and plots junction counts,
SCJ distance, raw inverted fraction, and Mash distance.

All metrics are measured by running the Rust binary; nothing is hardcoded. The
SNP-only control is compared against the unmutated reference, and every SV
variant (which already carries the 1% SNP background) is compared against the
control, so the Mash panel contrasts the SNP signal (~1e-2) with the SV
increment (~1e-7).

Also writes:
  - results/snp_sweep_metrics.csv: junction behaviour at 0.5/1/2/5% substitutions
  - results/table2_sv_metrics.csv: BcgI-only vs four-enzyme panel on the same
    controlled variants (source of manuscript Table 2)
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
OUT_PNG = ROOT / "figures" / "main" / "fig2_sv_sensitivity_rust.png"
OUT_PDF = ROOT / "figures" / "main" / "fig2_sv_sensitivity_rust.pdf"
OUT_CSV = ROOT / "results" / "figure2_rust_metrics.csv"
OUT_SWEEP = ROOT / "results" / "snp_sweep_metrics.csv"
OUT_TABLE2 = ROOT / "results" / "table2_sv_metrics.csv"

PANEL = "BcgI,AlfI,AloI,FalI"
BCGI = "BcgI"
BASE_SNP_RATE = 0.01

FIGURE_ROWS = [
    ("snp_1pct", "SNP 1%\n(vs ref)"),
    ("inv_100kb", "Inv\n100 kb"),
    ("inv_500kb", "Inv\n500 kb"),
    ("tra_500kb", "Tra\n500 kb"),
]
SWEEP_RATES = [0.005, 0.02, 0.05]


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


def translocate_segment(seq, start, end, insert_pos):
    seg = seq[start:end]
    remainder = seq[:start] + seq[end:]
    if insert_pos > start:
        insert_pos -= (end - start)
    return remainder[:insert_pos] + seg + remainder[insert_pos:]


def write_fasta(path, seq, acc):
    rec = SeqRecord(Seq(seq), id=acc, description="")
    SeqIO.write(rec, path, "fasta")


def read_fasta_seq(path):
    return str(next(SeqIO.parse(str(path), "fasta")).seq)


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
    """Approximate Mash distance from exact canonical k-mer Jaccard (stride 1)."""
    a = canonical_kmers(seq_a, k)
    b = canonical_kmers(seq_b, k)
    inter = len(a & b)
    union = len(a | b)
    jaccard = inter / union if union > 0 else 0.0
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
    """Return (variants, comparisons).

    variants: label -> (fasta path, genome_id). All variants carry the 1% SNP
    background except the sweep controls, which are built directly on the
    reference at 0.5/2/5% to isolate the substitution effect.
    """
    n = len(ref_seq)
    variants = {}

    def add(label, seq, stem):
        p = out_dir / f"{stem}.fasta"
        gid = f"{ref_id}_{stem}"
        write_fasta(p, seq, gid)
        variants[label] = (p, gid)

    # SNP background control (1%)
    control_seq = introduce_snps(ref_seq, BASE_SNP_RATE, rng)
    add("control (1% SNPs)", control_seq, "control_snps")

    # inversions on the SNP background
    for size in [100_000, 500_000]:
        start = n // 2 - size // 2
        end = start + size
        seq = invert_segment(control_seq, start, end)
        add(f"inversion {size//1000} kb", seq, f"inv{size}")

    # translocation 500 kb on the SNP background
    size = 500_000
    start = n // 4
    end = start + size
    insert_pos = 3 * n // 4
    seq = translocate_segment(control_seq, start, end, insert_pos)
    add("translocation 500 kb", seq, "tra500k")

    return variants


def run_arm(syn2b, enzymes, fasta_map, gid_of, comparisons, tmpdir, arm_tag):
    """Digest all fastas with `enzymes` and measure every comparison.

    fasta_map: label -> fasta path; comparisons: list of (row_key, label_a, label_b).
    Returns {row_key: metrics dict} and {label: shared_tags vs nothing...}.
    """
    tgt_dir = tmpdir / f"tgts_{arm_tag}"
    tgt_dir.mkdir()
    for label, fasta in fasta_map.items():
        digest(syn2b, enzymes, fasta, tgt_dir / f"{Path(fasta).stem}.tgt")

    rows = {}
    for row_key, label_a, label_b in comparisons:
        syn_out = tmpdir / f"synteny_matrix_{arm_tag}"
        rows[row_key] = synteny_pair(syn2b, tgt_dir, syn_out,
                                     gid_of[label_a], gid_of[label_b])
    return rows


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

        # reference itself as a digest target (comparison baseline)
        ref_fasta = tmpdir / "reference.fasta"
        write_fasta(ref_fasta, ref_seq, f"{ref_id}_reference")

        variants = build_variants(ref_seq, ref_id, rng, tmpdir)

        # sweep controls at 0.5/2/5% built directly on the reference
        sweep_labels = []
        for rate in SWEEP_RATES:
            seq = introduce_snps(ref_seq, rate, rng)
            stem = f"sweep_{rate}"
            p = tmpdir / f"{stem}.fasta"
            write_fasta(p, seq, f"{ref_id}_{stem}")
            variants[f"sweep {rate}"] = (p, f"{ref_id}_{stem}")
            sweep_labels.append((rate, f"sweep {rate}"))

        fasta_map = {"reference": ref_fasta}
        gid_of = {"reference": f"{ref_id}_reference"}
        for label, (fasta, gid) in variants.items():
            fasta_map[label] = fasta
            gid_of[label] = gid

        comparisons = (
            [("snp_1pct", "reference", "control (1% SNPs)")]
            + [(key, "control (1% SNPs)", lbl) for key, lbl in
               [("inv_100kb", "inversion 100 kb"),
                ("inv_500kb", "inversion 500 kb"),
                ("tra_500kb", "translocation 500 kb")]]
        )

        panel_rows = run_arm(syn2b, PANEL, fasta_map, gid_of, comparisons, tmpdir, "panel")
        bcgI_rows = run_arm(syn2b, BCGI, fasta_map, gid_of, comparisons, tmpdir, "bcgI")

        # sweep comparisons (reference vs each sweep control), panel only
        sweep_rows = []
        for rate, label in sweep_labels:
            m = synteny_pair(syn2b, tmpdir / "tgts_panel",
                             tmpdir / "synteny_matrix_panel",
                             gid_of["reference"], gid_of[label])
            m["snp_rate"] = rate
            seq_b = read_fasta_seq(fasta_map[label])
            m["mash_distance"] = mash_distance(ref_seq, seq_b, k=21)
            sweep_rows.append(m)

        # Mash distances for the figure rows
        seqs = {label: read_fasta_seq(p) for label, p in fasta_map.items()}
        mash = {
            "snp_1pct": mash_distance(seqs["reference"], seqs["control (1% SNPs)"], k=21),
            "inv_100kb": mash_distance(seqs["control (1% SNPs)"], seqs["inversion 100 kb"], k=21),
            "inv_500kb": mash_distance(seqs["control (1% SNPs)"], seqs["inversion 500 kb"], k=21),
            "tra_500kb": mash_distance(seqs["control (1% SNPs)"], seqs["translocation 500 kb"], k=21),
        }

    # write figure CSV (panel arm, measured)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["label", "breakpoints", "scj_distance",
                                                "raw_inverted_fraction", "observable_fraction",
                                                "shared_tags", "mash_distance"])
        writer.writeheader()
        for key, label in FIGURE_ROWS:
            m = panel_rows[key]
            writer.writerow({"label": label.replace("\n", " "), "breakpoints": m["breakpoints"],
                             "scj_distance": m["scj_distance"],
                             "raw_inverted_fraction": m["raw_inverted_fraction"],
                             "observable_fraction": m["observable_fraction"],
                             "shared_tags": m["shared_tags"],
                             "mash_distance": mash[key]})
    print(f"Wrote {out_csv}")

    # write SNP sweep CSV
    out_sweep = ROOT / "results" / "snp_sweep_metrics.csv"
    with open(out_sweep, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["snp_rate", "breakpoints", "scj_distance",
                                                "raw_inverted_fraction", "observable_fraction",
                                                "shared_tags", "mash_distance"])
        writer.writeheader()
        for m in sweep_rows:
            writer.writerow({"snp_rate": m["snp_rate"], "breakpoints": m["breakpoints"],
                             "scj_distance": m["scj_distance"],
                             "raw_inverted_fraction": m["raw_inverted_fraction"],
                             "observable_fraction": m["observable_fraction"],
                             "shared_tags": m["shared_tags"],
                             "mash_distance": m["mash_distance"]})
    print(f"Wrote {out_sweep}")

    # write Table 2 source CSV (both enzyme arms)
    out_t2 = ROOT / "results" / "table2_sv_metrics.csv"
    with open(out_t2, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["enzyme_arm", "condition", "junctions",
                                                "scj_distance", "raw_inverted_fraction",
                                                "shared_tags"])
        writer.writeheader()
        for key, _label in FIGURE_ROWS:
            for arm, rows in [("BcgI", bcgI_rows), ("BcgI+AlfI+AloI+FalI", panel_rows)]:
                m = rows[key]
                writer.writerow({"enzyme_arm": arm, "condition": key,
                                 "junctions": m["breakpoints"],
                                 "scj_distance": m["scj_distance"],
                                 "raw_inverted_fraction": m["raw_inverted_fraction"],
                                 "shared_tags": m["shared_tags"]})
    print(f"Wrote {out_t2}")

    # plot
    labels = [lbl for _, lbl in FIGURE_ROWS]
    x = np.arange(len(labels))

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.5), constrained_layout=True)

    ax = axes[0]
    vals = [panel_rows[k]["breakpoints"] for k, _ in FIGURE_ROWS]
    ax.bar(x, vals, color="#1f77b4")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0, fontsize=8)
    ax.set_ylabel("Syn2b junctions", fontsize=10)
    ax.set_title("Junction count", fontsize=11)
    ax.set_ylim(0, max(vals) + 0.6)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.08, str(v), ha="center", va="bottom")

    ax = axes[1]
    vals = [panel_rows[k]["scj_distance"] for k, _ in FIGURE_ROWS]
    ax.bar(x, vals, color="#ff7f0e")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("SCJ distance", fontsize=10)
    ax.set_title("SCJ distance", fontsize=11)
    ax.set_ylim(0, max(vals) + 1)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.15, str(v), ha="center", va="bottom")

    ax = axes[2]
    vals = [panel_rows[k]["raw_inverted_fraction"] for k, _ in FIGURE_ROWS]
    ax.bar(x, vals, color="#2ca02c")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Raw inverted fraction", fontsize=10)
    ax.set_title("Inverted fraction", fontsize=11)
    ax.set_ylim(0, max(vals) + 0.02)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.003, f"{v:.3f}", ha="center", va="bottom")

    ax = axes[3]
    vals = [mash[k] for k, _ in FIGURE_ROWS]
    ax.bar(x, vals, color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Mash distance (k=21)", fontsize=10)
    ax.set_title("Mash distance", fontsize=11)
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 0.05)
    for i, v in enumerate(vals):
        ax.text(i, v * 1.6, f"{v:.1e}", ha="center", va="bottom", fontsize=8)
    ax.axhline(0, color="black", lw=0.5)
    ax.annotate("SV conditions add nothing\nabove the SNP background\n(~1e-7 vs ~1e-2)",
                xy=(0.03, 0.03), xycoords="axes fraction", fontsize=7, color="black",
                va="bottom")

    out_png = Path(args.out_png)
    out_pdf = Path(args.out_pdf)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()

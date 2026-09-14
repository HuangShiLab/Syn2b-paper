#!/usr/bin/env python3
"""Measure Syn2b SV detection-size resolution for the production panel vs BcgI.

For each inversion size, random positions (10 replicates) are inverted on an
*E. coli* K-12 genome carrying a 1% SNP background. An event is "detected" when
the Rust implementation reports both junctions (breakpoints >= 2). The reported
resolution is the smallest size with >=50% detection, for the production
four-enzyme panel (BcgI+AlfI+AloI+FalI) and for BcgI alone.

Output: results/detection_size_metrics.csv
"""
import argparse
import csv
import random
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "ecoli_k12_MG1655.fasta"
OUT_CSV = ROOT / "results" / "detection_size_metrics.csv"

PANEL = "BcgI,AlfI,AloI,FalI"
BCGI = "BcgI"
SIZES = [2_000, 4_000, 8_000, 16_000, 32_000]
REPLICATES = 10
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
    return seq[:start] + str(Seq(seg).reverse_complement()) + seq[end:]


def write_fasta(path, seq, acc):
    SeqIO.write(SeqRecord(Seq(seq), id=acc, description=""), str(path), "fasta")


def digest(syn2b, enzymes, fasta, tgt):
    subprocess.run([str(syn2b), "digest", "--enzymes", enzymes,
                    "--input", str(fasta), "--output", str(tgt)],
                   check=True, capture_output=True)


def get_pair_row(syn2b, tgt_dir, out_prefix, id_a, id_b):
    subprocess.run([str(syn2b), "synteny", "--input", str(tgt_dir),
                    "--output", str(out_prefix)], check=True, capture_output=True)
    with open(out_prefix) as fh:
        reader = csv.DictReader([ln for ln in fh if not ln.startswith("#")])
        for row in reader:
            if {row["genome_A"], row["genome_B"]} == {id_a, id_b}:
                return row
    raise RuntimeError(f"pair {id_a} vs {id_b} not found")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--syn2b", default=str(ROOT.parent / "Syn2b" / "target" / "release" / "syn2b"))
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    syn2b = Path(args.syn2b)
    if not syn2b.exists():
        sys.exit(f"syn2b binary not found: {syn2b}")

    rng = random.Random(args.seed)
    ref_seq, ref_id = load_genome(DATA)
    n = len(ref_seq)
    control_seq = introduce_snps(ref_seq, SNP_RATE, rng)

    rows = []
    for arm, enzymes in [("panel", PANEL), ("BcgI", BCGI)]:
        for size in SIZES:
            detected = 0
            inv_frac = []
            with tempfile.TemporaryDirectory(prefix="det_size_") as tmp:
                tmpdir = Path(tmp)
                tgt_dir = tmpdir / "tgts"
                tgt_dir.mkdir()
                write_fasta(tmpdir / "control.fasta", control_seq, f"{ref_id}_control")
                digest(syn2b, enzymes, tmpdir / "control.fasta", tgt_dir / "control.tgt")
                gids = []
                for rep in range(REPLICATES):
                    start = rng.randint(n // 10, n - n // 10 - size)
                    var_seq = invert_segment(control_seq, start, start + size)
                    gid = f"{ref_id}_inv{size}_r{rep}"
                    p = tmpdir / f"inv{size}_r{rep}.fasta"
                    write_fasta(p, var_seq, gid)
                    digest(syn2b, enzymes, p, tgt_dir / f"inv{size}_r{rep}.tgt")
                    gids.append(gid)
                for gid in gids:
                    row = get_pair_row(syn2b, tgt_dir, tmpdir / "matrix",
                                       f"{ref_id}_control", gid)
                    bp = int(row["breakpoints"])
                    detected += bp >= 2
                    inv_frac.append(float(row["raw_inverted_fraction"]))
            rate = detected / REPLICATES
            rows.append({"enzyme_arm": arm, "event_size_bp": size,
                         "detected": detected, "replicates": REPLICATES,
                         "detection_rate": rate,
                         "mean_raw_inverted_fraction": np.mean(inv_frac)})
            print(f"{arm} size={size}: {detected}/{REPLICATES} "
                  f"(mean inv frac {np.mean(inv_frac):.5f})")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()

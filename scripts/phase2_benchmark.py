#!/usr/bin/env python3
"""
phase2_benchmark.py

Benchmark Syn2b vs BLAST-based approach (SynTracker-like) on E. coli genomes.
Compares:
  - Syn2b: multi-enzyme digest (O(N) time)
  - BLAST: database creation + search (SynTracker Step 1-2)
"""

import sys
import os
import time
import random
import tempfile
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from simulate_rearrangement import (
    parse_fasta, write_fasta, substitute, inversion, translocation,
    digest_multi
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def generate_test_genomes(reference_seq, n_genomes=20, rng=None):
    """Generate N diverse E. coli genomes for benchmarking."""
    if rng is None:
        rng = random.Random(42)
    
    genomes = [("ref", reference_seq)]
    for i in range(n_genomes):
        seq = substitute(reference_seq, mu=0.001, rng=rng)
        # Add occasional SV
        if rng.random() < 0.3:
            seq = inversion(seq, 100000, rng=rng)
        if rng.random() < 0.3:
            seq = translocation(seq, 50000, rng=rng)
        genomes.append((f"genome_{i+1}", seq))
    return genomes


def benchmark_syn2b(genomes):
    """Benchmark Syn2b multi-enzyme digestion."""
    times = []
    for label, seq in genomes:
        t0 = time.time()
        tags = digest_multi(seq, include_cjepi=True)
        t1 = time.time()
        times.append(t1 - t0)
    return times


def benchmark_blast(genomes, blastn_path, tmpdir):
    """
    Benchmark BLAST-based approach (SynTracker Step 1-2).
    1. Build BLAST db from all genomes
    2. Search with 1-kb central regions (4kb apart) from reference
    """
    ref_label, ref_seq = genomes[0]
    
    # Write all genomes to multi-FASTA
    multi_fasta = os.path.join(tmpdir, "all_genomes.fasta")
    with open(multi_fasta, "w") as fh:
        for label, seq in genomes:
            fh.write(f">{label}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
    
    # Step 1: Build BLAST db
    t0 = time.time()
    db_path = os.path.join(tmpdir, "blast_db")
    cmd = [blastn_path, "-dbtype", "nucl", "-in", multi_fasta, "-out", db_path]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        print(f"BLAST db build failed: {e}")
        return None, None
    t_db = time.time() - t0
    
    # Step 2: Create 1-kb central regions (4kb apart) from reference
    central_regions = []
    region_size = 1000
    spacing = 4000
    for start in range(0, len(ref_seq) - region_size, spacing):
        central_regions.append(ref_seq[start:start + region_size])
    
    # Write central regions to query FASTA
    query_fasta = os.path.join(tmpdir, "queries.fasta")
    with open(query_fasta, "w") as fh:
        for i, region in enumerate(central_regions):
            fh.write(f">region_{i}\n")
            fh.write(region + "\n")
    
    # Step 3: BLAST search
    t0 = time.time()
    output_file = os.path.join(tmpdir, "blast_out.txt")
    cmd = [
        blastn_path, "-query", query_fasta, "-db", db_path,
        "-outfmt", "6 qseqid sseqid pident length qlen slen",
        "-perc_identity", "97", "-qcov_hsp_perc", "70",
        "-out", output_file, "-num_threads", "1"
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=300)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        print(f"BLAST search failed: {e}")
        return t_db, None
    t_search = time.time() - t0
    
    return t_db, t_search


def plot_benchmark(syn2b_times, blast_times, n_genomes_list, output_png):
    """Generate benchmark comparison figure."""
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: per-genome digest time
    ax = axes[0]
    ax.bar(["Syn2b\n(digest)", "BLAST\n(db + search)"],
           [sum(syn2b_times) / len(syn2b_times),
            sum(blast_times) / len(blast_times) if blast_times else 0],
           color=["#4C78A8", "#E45756"],
           edgecolor="black", linewidth=1)
    ax.set_ylabel("Time per genome (seconds)", fontsize=12)
    ax.set_title("Single-genome processing time", fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")
    
    # Right: scaling with N genomes
    ax = axes[1]
    ax.plot(n_genomes_list, [sum(syn2b_times[:n]) for n in n_genomes_list],
            'o-', color="#4C78A8", linewidth=2, markersize=8, label="Syn2b")
    ax.plot(n_genomes_list, [sum(blast_times[:n]) for n in n_genomes_list],
            's-', color="#E45756", linewidth=2, markersize=8, label="BLAST-based")
    ax.set_xlabel("Number of genomes", fontsize=12)
    ax.set_ylabel("Cumulative time (seconds)", fontsize=12)
    ax.set_title("Scaling with dataset size", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    fig.suptitle("Syn2b vs BLAST-based (SynTracker-like) Speed Benchmark",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Benchmark figure saved to {output_png}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/e_coli_k12.fasta")
    parser.add_argument("--n-genomes", type=int, default=20)
    parser.add_argument("--blastn", default="/Users/shihuang/ncbi-blast/bin/blastn")
    parser.add_argument("--csv", default="/Users/shihuang/Documents/kimi/workspace/phase2_benchmark.csv")
    parser.add_argument("--png", default="/Users/shihuang/Documents/kimi/workspace/phase2_benchmark.png")
    args = parser.parse_args()
    
    # Check BLAST exists
    if not os.path.isfile(args.blastn):
        # Try makeblastdb instead
        blast_dir = os.path.dirname(args.blastn)
        makeblastdb = os.path.join(blast_dir, "makeblastdb")
        if os.path.isfile(makeblastdb):
            print(f"Using makeblastdb at {makeblastdb}")
        else:
            print(f"ERROR: BLAST not found at {args.blastn}")
            print("Please install BLAST or provide correct path with --blastn")
            return
    
    _, reference_seq = parse_fasta(args.input)
    print(f"Reference: {len(reference_seq):,} bp")
    print(f"Generating {args.n_genomes} test genomes...")
    
    rng = random.Random(42)
    genomes = generate_test_genomes(reference_seq, n_genomes=args.n_genomes, rng=rng)
    print(f"Generated {len(genomes)} genomes")
    
    # Benchmark Syn2b
    print("\n=== Benchmarking Syn2b ===")
    syn2b_times = benchmark_syn2b(genomes)
    print(f"Syn2b digest: mean={sum(syn2b_times)/len(syn2b_times):.4f}s, total={sum(syn2b_times):.2f}s")
    
    # Benchmark BLAST
    print("\n=== Benchmarking BLAST (SynTracker-like) ===")
    tmpdir = tempfile.mkdtemp(prefix="syn2b_bench_")
    
    # Use makeblastdb for db creation, blastn for search
    blast_dir = os.path.dirname(args.blastn)
    makeblastdb = os.path.join(blast_dir, "makeblastdb")
    blastn = args.blastn
    
    # Write all genomes
    multi_fasta = os.path.join(tmpdir, "all_genomes.fasta")
    with open(multi_fasta, "w") as fh:
        for label, seq in genomes:
            fh.write(f">{label}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
    
    # Build BLAST db
    t0 = time.time()
    db_path = os.path.join(tmpdir, "blast_db")
    cmd = [makeblastdb, "-dbtype", "nucl", "-in", multi_fasta, "-out", db_path]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        t_db = time.time() - t0
        print(f"BLAST db build: {t_db:.2f}s")
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        print(f"BLAST db build failed: {e}")
        t_db = None
    
    # Create central regions
    ref_label, ref_seq = genomes[0]
    central_regions = []
    region_size = 1000
    spacing = 4000
    for start in range(0, len(ref_seq) - region_size, spacing):
        central_regions.append(ref_seq[start:start + region_size])
    print(f"Created {len(central_regions)} central regions")
    
    query_fasta = os.path.join(tmpdir, "queries.fasta")
    with open(query_fasta, "w") as fh:
        for i, region in enumerate(central_regions):
            fh.write(f">region_{i}\n")
            fh.write(region + "\n")
    
    # BLAST search
    if t_db is not None:
        t0 = time.time()
        output_file = os.path.join(tmpdir, "blast_out.txt")
        cmd = [
            blastn, "-query", query_fasta, "-db", db_path,
            "-outfmt", "6 qseqid sseqid pident length qlen slen",
            "-perc_identity", "97", "-qcov_hsp_perc", "70",
            "-out", output_file, "-num_threads", "1"
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=300)
            t_search = time.time() - t0
            print(f"BLAST search: {t_search:.2f}s")
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            print(f"BLAST search failed: {e}")
            t_search = None
    else:
        t_search = None
    
    # Calculate per-genome BLAST time (approximate)
    # BLAST time scales with number of genomes, but we ran all at once
    # Approximate: total BLAST time / N genomes for comparison
    if t_db is not None and t_search is not None:
        total_blast = t_db + t_search
        blast_per_genome = total_blast / len(genomes)
        print(f"\nTotal BLAST time: {total_blast:.2f}s")
        print(f"Per-genome equivalent: {blast_per_genome:.4f}s")
        
        # For scaling plot, assume linear scaling (rough approximation)
        blast_times = [blast_per_genome] * len(genomes)
    else:
        blast_times = []
    
    # Save CSV
    with open(args.csv, "w", newline="") as fh:
        import csv
        writer = csv.writer(fh)
        writer.writerow(["tool", "genome_index", "time_seconds"])
        for i, t in enumerate(syn2b_times):
            writer.writerow(["syn2b", i, f"{t:.6f}"])
        if blast_times:
            for i, t in enumerate(blast_times):
                writer.writerow(["blast", i, f"{t:.6f}"])
    print(f"\nCSV saved to {args.csv}")
    
    # Plot
    n_genomes_list = list(range(1, len(genomes) + 1, max(1, len(genomes) // 10)))
    if n_genomes_list[-1] != len(genomes):
        n_genomes_list.append(len(genomes))
    
    if blast_times:
        plot_benchmark(syn2b_times, blast_times, n_genomes_list, args.png)
    else:
        print("BLAST benchmarking failed, skipping figure")
    
    # Cleanup
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)
    
    print("\n=== Summary ===")
    print(f"Syn2b:  {sum(syn2b_times):.2f}s total for {len(genomes)} genomes")
    if blast_times:
        print(f"BLAST:  {total_blast:.2f}s total for {len(genomes)} genomes")
        speedup = total_blast / sum(syn2b_times)
        print(f"Speedup: {speedup:.1f}x")


if __name__ == "__main__":
    main()

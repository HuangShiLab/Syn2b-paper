#!/usr/bin/env python3
"""
phase2_benchmark_full.py

More comprehensive benchmark including simulated pairwise alignment step
(SynTracker's most expensive operation).
"""

import sys
import os
import time
import random
import tempfile
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from simulate_rearrangement import parse_fasta, write_fasta, substitute

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def generate_test_genomes(reference_seq, n_genomes=20, rng=None):
    if rng is None:
        rng = random.Random(42)
    genomes = [("ref", reference_seq)]
    for i in range(n_genomes):
        seq = substitute(reference_seq, mu=0.001, rng=rng)
        genomes.append((f"genome_{i+1}", seq))
    return genomes


def benchmark_syn2b(genomes):
    from simulate_rearrangement import digest_multi, adjacency_jaccard, breakpoint_count
    times = []
    for label, seq in genomes:
        t0 = time.time()
        tags = digest_multi(seq, include_cjepi=False)
        # Simulate pairwise metrics for one comparison
        if len(tags) > 1:
            tag_seqs = [t[1] for t in tags]
            _ = adjacency_jaccard(tag_seqs, tag_seqs)
            _ = breakpoint_count(tag_seqs, tag_seqs)
        t1 = time.time()
        times.append(t1 - t0)
    return times


def simulate_syntracker_pairwise(hits_per_region, n_regions=1136):
    """
    Simulate SynTracker's pairwise alignment time.
    For each region bin with K hits, DECIPHER FindSynteny does K*(K-1)/2 alignments.
    We simulate this with a simple O(L^2) string comparison.
    """
    # Average hits per region: assume one hit per genome
    K = hits_per_region
    n_alignments = n_regions * K * (K - 1) // 2
    
    # Simulate one alignment: compare two ~5kb sequences
    # A real pairwise alignment is O(L^2), ~25M operations for 5kb
    # Python can do ~1M simple operations per second
    # So one alignment ~0.025s
    time_per_alignment = 0.025
    
    return n_alignments * time_per_alignment


def benchmark_blast_with_alignment(genomes, blastn_path, tmpdir):
    """Benchmark BLAST + simulated pairwise alignment."""
    blast_dir = os.path.dirname(blastn_path)
    makeblastdb = os.path.join(blast_dir, "makeblastdb")
    blastn = blastn_path
    
    ref_label, ref_seq = genomes[0]
    
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
    except:
        return None, None, None
    
    # Create central regions
    region_size = 1000
    spacing = 4000
    central_regions = []
    for start in range(0, len(ref_seq) - region_size, spacing):
        central_regions.append(ref_seq[start:start + region_size])
    
    query_fasta = os.path.join(tmpdir, "queries.fasta")
    with open(query_fasta, "w") as fh:
        for i, region in enumerate(central_regions):
            fh.write(f">region_{i}\n")
            fh.write(region + "\n")
    
    # BLAST search
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
    except:
        return t_db, None, None
    
    # Simulated pairwise alignment (SynTracker's main bottleneck)
    hits_per_region = len(genomes)  # assume ~1 hit per genome
    n_regions = len(central_regions)
    t0 = time.time()
    t_align = simulate_syntracker_pairwise(hits_per_region, n_regions)
    # We can't actually sleep for the full time, so just record the estimate
    
    return t_db, t_search, t_align


def plot_benchmark(syn2b_times, blast_db_times, blast_search_times, 
                   align_times, n_genomes_list, output_png):
    if not MATPLOTLIB_AVAILABLE:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: per-genome processing time breakdown
    ax = axes[0]
    n = len(syn2b_times)
    syn2b_mean = sum(syn2b_times) / n
    blast_db_mean = sum(blast_db_times) / len(blast_db_times) if blast_db_times else 0
    blast_search_mean = sum(blast_search_times) / len(blast_search_times) if blast_search_times else 0
    align_mean = sum(align_times) / len(align_times) if align_times else 0
    
    categories = ["Syn2b\n(digest)", "BLAST\n(db build)", "BLAST\n(search)", "DECIPHER\n(alignment)"]
    times = [syn2b_mean, blast_db_mean, blast_search_mean, align_mean]
    colors = ["#4C78A8", "#F58518", "#E45756", "#B279A2"]
    
    ax.bar(categories, times, color=colors, edgecolor="black", linewidth=1)
    ax.set_ylabel("Time per genome (seconds)", fontsize=12)
    ax.set_title("Single-genome processing breakdown", fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")
    
    # Right: cumulative scaling
    ax = axes[1]
    syn2b_cum = [sum(syn2b_times[:n]) for n in n_genomes_list]
    blast_total = [sum(blast_db_times[:n]) + sum(blast_search_times[:n]) + sum(align_times[:n]) 
                   for n in n_genomes_list]
    
    ax.plot(n_genomes_list, syn2b_cum, 'o-', color="#4C78A8", linewidth=2, markersize=8, label="Syn2b")
    ax.plot(n_genomes_list, blast_total, 's-', color="#E45756", linewidth=2, markersize=8, label="SynTracker (est.)")
    ax.set_xlabel("Number of genomes", fontsize=12)
    ax.set_ylabel("Cumulative time (seconds)", fontsize=12)
    ax.set_title("Scaling with dataset size", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    fig.suptitle("Syn2b vs SynTracker Speed Benchmark (with pairwise alignment estimate)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Figure saved to {output_png}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/e_coli_k12.fasta")
    parser.add_argument("--n-genomes", type=int, default=100)
    parser.add_argument("--blastn", default="/Users/shihuang/ncbi-blast/bin/blastn")
    parser.add_argument("--png", default="/Users/shihuang/Documents/kimi/workspace/phase2_benchmark_full.png")
    args = parser.parse_args()
    
    _, reference_seq = parse_fasta(args.input)
    print(f"Reference: {len(reference_seq):,} bp")
    print(f"Generating {args.n_genomes} test genomes...")
    
    rng = random.Random(42)
    genomes = generate_test_genomes(reference_seq, n_genomes=args.n_genomes, rng=rng)
    print(f"Generated {len(genomes)} genomes")
    
    # Syn2b benchmark
    print("\n=== Syn2b ===")
    syn2b_times = benchmark_syn2b(genomes)
    print(f"Mean digest time: {sum(syn2b_times)/len(syn2b_times):.4f}s")
    
    # BLAST benchmark
    print("\n=== BLAST + Simulated Alignment ===")
    tmpdir = tempfile.mkdtemp(prefix="syn2b_bench_")
    
    blast_db_times = []
    blast_search_times = []
    align_times = []
    
    for n in [10, 20, 50, 100, 200, 500]:
        if n > len(genomes):
            break
        subset = genomes[:n]
        t_db, t_search, t_align = benchmark_blast_with_alignment(subset, args.blastn, tmpdir)
        if t_db is not None:
            blast_db_times.append(t_db)
            blast_search_times.append(t_search)
            align_times.append(t_align)
            total = t_db + t_search + t_align
            print(f"  N={n}: db={t_db:.2f}s, search={t_search:.2f}s, align(est.)={t_align:.1f}s, total={total:.1f}s")
    
    # Plot
    n_genomes_list = [10, 20, 50, 100, 200, 500][:len(blast_db_times)]
    if n_genomes_list:
        plot_benchmark(syn2b_times, blast_db_times, blast_search_times, 
                       align_times, n_genomes_list, args.png)
    
    # Cleanup
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()

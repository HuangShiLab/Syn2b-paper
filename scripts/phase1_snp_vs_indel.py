#!/usr/bin/env python3
"""
phase1_snp_vs_indel.py

Reproduce SynTracker Figure 2a-b:
(a) SNP-only evolution: adjacency Jaccard insensitive to SNPs
(b) Indel-only evolution: adjacency Jaccard sensitive to indels

Simulates bacterial population evolution over N generations,
sampling 20 genomes at each time point for pairwise comparison.
"""

import sys
import os
import random
import math
import csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from simulate_rearrangement import (
    parse_fasta, write_fasta, substitute, insertion, deletion,
    digest_multi, adjacency_jaccard, breakpoint_count, kendall_tau_rank
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def simulate_population_evolution(reference_seq, mode="snp", n_generations=50,
                                   mutation_rate=1e-6, sample_size=20,
                                   region_length=20000, rng=None):
    """
    Simulate population evolution.
    
    mode: "snp" or "indel"
    mutation_rate: per-nucleotide per-generation (1e-6 for SNP, 1e-7 for indel)
    sample_size: number of genomes to sample at each time point
    region_length: length of genomic region to analyze (default 20kb like SynTracker)
    
    Returns: list of (generation, list_of_sampled_sequences)
    """
    if rng is None:
        rng = random.Random(42)
    
    # Start with a single ancestor
    ancestor = reference_seq[:region_length]
    
    # For simplicity: evolve a single lineage for N generations,
    # then branch into sample_size descendants at each time point
    # (this approximates the Bacmeta simulation in SynTracker)
    
    results = []
    current_pop = [ancestor] * sample_size
    
    for gen in range(n_generations + 1):
        # Sample current population
        sampled = []
        for i in range(sample_size):
            # Each individual diverges slightly from ancestor
            seq = current_pop[i]
            if mode == "snp":
                # Introduce SNPs only
                seq = substitute(seq, mu=mutation_rate, rng=rng)
            elif mode == "indel":
                # Introduce indels only (no SNPs)
                # Random indel events
                n_events = rng.poissonvariate(len(seq) * mutation_rate) if hasattr(rng, 'poissonvariate') else max(1, int(len(seq) * mutation_rate * rng.random() * 2))
                if hasattr(rng, 'poissonvariate'):
                    n_events = rng.poissonvariate(len(seq) * mutation_rate)
                else:
                    # Approximate Poisson
                    lam = len(seq) * mutation_rate
                    n_events = 0
                    p = 1.0
                    while p > math.exp(-lam):
                        p *= rng.random()
                        n_events += 1
                    n_events -= 1
                
                for _ in range(n_events):
                    if rng.random() < 0.5:
                        # Insertion
                        ins_size = rng.randint(1, 50)
                        seq = insertion(seq, ins_size, rng=rng)
                    else:
                        # Deletion
                        del_size = rng.randint(1, 50)
                        seq = deletion(seq, del_size, rng=rng)
            sampled.append(seq)
        
        results.append((gen, sampled))
        
        # Update population for next generation
        if mode == "snp":
            current_pop = [substitute(seq, mu=mutation_rate, rng=rng) for seq in current_pop]
        elif mode == "indel":
            new_pop = []
            for seq in current_pop:
                lam = len(seq) * mutation_rate
                n_events = 0
                p = 1.0
                while p > math.exp(-lam):
                    p *= rng.random()
                    n_events += 1
                n_events = max(0, n_events - 1)
                for _ in range(n_events):
                    if rng.random() < 0.5:
                        seq = insertion(seq, rng.randint(1, 50), rng=rng)
                    else:
                        seq = deletion(seq, rng.randint(1, 50), rng=rng)
                new_pop.append(seq)
            current_pop = new_pop
    
    return results


def compute_pairwise_metrics(sampled_sequences):
    """Compute all pairwise metrics for a set of sampled genomes."""
    n = len(sampled_sequences)
    
    # Pre-compute tags
    all_tags = []
    for seq in sampled_sequences:
        tags = digest_multi(seq)
        all_tags.append([t[1] for t in tags])
    
    adj_jaccards = []
    break_points = []
    kendall_taus = []
    
    for i in range(n):
        for j in range(i + 1, n):
            aj = adjacency_jaccard(all_tags[i], all_tags[j])
            bp = breakpoint_count(all_tags[i], all_tags[j])
            kt = kendall_tau_rank(all_tags[i], all_tags[j])
            
            adj_jaccards.append(aj)
            break_points.append(bp)
            if kt is not None:
                kendall_taus.append(kt)
    
    return adj_jaccards, break_points, kendall_taus


def plot_figure2a_b(snp_results, indel_results, output_png):
    """Generate Figure 2a-b style plot."""
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Top row: SNP-only
    ax = axes[0, 0]
    gens = [r[0] for r in snp_results]
    adj_medians = [sum(r[2]) / len(r[2]) for r in snp_results]  # mean adj_jaccard
    ax.plot(gens, adj_medians, 'o-', color='#4C78A8', linewidth=2, markersize=6)
    ax.set_xlabel('Generation', fontsize=12)
    ax.set_ylabel('Mean Adjacency Jaccard', fontsize=12)
    ax.set_title('(a) SNP-only evolution', fontsize=13, fontweight='bold')
    ax.set_ylim(0.2, 1.0)
    ax.grid(True, alpha=0.3)
    
    ax = axes[0, 1]
    bp_medians = [sum(r[3]) / len(r[3]) for r in snp_results]
    ax.plot(gens, bp_medians, 'o-', color='#4C78A8', linewidth=2, markersize=6)
    ax.set_xlabel('Generation', fontsize=12)
    ax.set_ylabel('Mean Breakpoint Count', fontsize=12)
    ax.set_title('(a) SNP-only evolution', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Bottom row: Indel-only
    ax = axes[1, 0]
    gens = [r[0] for r in indel_results]
    adj_medians = [sum(r[2]) / len(r[2]) for r in indel_results]
    ax.plot(gens, adj_medians, 'o-', color='#E45756', linewidth=2, markersize=6)
    ax.set_xlabel('Generation', fontsize=12)
    ax.set_ylabel('Mean Adjacency Jaccard', fontsize=12)
    ax.set_title('(b) Indel-only evolution', fontsize=13, fontweight='bold')
    ax.set_ylim(0.2, 1.0)
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 1]
    bp_medians = [sum(r[3]) / len(r[3]) for r in indel_results]
    ax.plot(gens, bp_medians, 'o-', color='#E45756', linewidth=2, markersize=6)
    ax.set_xlabel('Generation', fontsize=12)
    ax.set_ylabel('Mean Breakpoint Count', fontsize=12)
    ax.set_title('(b) Indel-only evolution', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    fig.suptitle('Syn2b Sensitivity Profile (SynTracker Figure 2a-b reproduction)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_png}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/e_coli_k12.fasta")
    parser.add_argument("--generations", type=int, default=30)
    parser.add_argument("--snp-rate", type=float, default=1e-6)
    parser.add_argument("--indel-rate", type=float, default=1e-7)
    parser.add_argument("--sample-size", type=int, default=20)
    parser.add_argument("--region-length", type=int, default=20000)
    parser.add_argument("--csv", default="/Users/shihuang/Documents/kimi/workspace/phase1_results.csv")
    parser.add_argument("--png", default="/Users/shihuang/Documents/kimi/workspace/phase1_figure2ab.png")
    args = parser.parse_args()
    
    _, reference_seq = parse_fasta(args.input)
    print(f"Reference: {len(reference_seq):,} bp")
    
    rng = random.Random(42)
    
    # Phase 1a: SNP-only evolution
    print("\n=== SNP-only evolution ===")
    snp_results = []
    snp_evolution = simulate_population_evolution(
        reference_seq, mode="snp", n_generations=args.generations,
        mutation_rate=args.snp_rate, sample_size=args.sample_size,
        region_length=args.region_length, rng=rng
    )
    
    for gen, sampled in snp_evolution:
        if gen % 5 == 0 or gen == args.generations:
            adj, bp, kt = compute_pairwise_metrics(sampled)
            mean_adj = sum(adj) / len(adj)
            mean_bp = sum(bp) / len(bp)
            mean_kt = sum(kt) / len(kt) if kt else None
            snp_results.append((gen, sampled, adj, bp, kt, mean_adj, mean_bp, mean_kt))
            print(f"  Gen {gen:3d}: n_pairs={len(adj)}, mean_adj={mean_adj:.4f}, mean_bp={mean_bp:.1f}")
    
    # Phase 1b: Indel-only evolution
    print("\n=== Indel-only evolution ===")
    indel_results = []
    indel_evolution = simulate_population_evolution(
        reference_seq, mode="indel", n_generations=args.generations,
        mutation_rate=args.indel_rate, sample_size=args.sample_size,
        region_length=args.region_length, rng=rng
    )
    
    for gen, sampled in indel_evolution:
        if gen % 5 == 0 or gen == args.generations:
            adj, bp, kt = compute_pairwise_metrics(sampled)
            mean_adj = sum(adj) / len(adj)
            mean_bp = sum(bp) / len(bp)
            mean_kt = sum(kt) / len(kt) if kt else None
            indel_results.append((gen, sampled, adj, bp, kt, mean_adj, mean_bp, mean_kt))
            print(f"  Gen {gen:3d}: n_pairs={len(adj)}, mean_adj={mean_adj:.4f}, mean_bp={mean_bp:.1f}")
    
    # Save CSV
    with open(args.csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["mode", "generation", "n_pairs", "mean_adj_jaccard", "mean_breakpoints", "mean_kendall_tau"])
        for gen, _, _, _, _, mean_adj, mean_bp, mean_kt in snp_results:
            writer.writerow(["snp", gen, len(sampled) * (len(sampled) - 1) // 2, f"{mean_adj:.6f}", f"{mean_bp:.2f}", f"{mean_kt:.6f}" if mean_kt else "NA"])
        for gen, _, _, _, _, mean_adj, mean_bp, mean_kt in indel_results:
            writer.writerow(["indel", gen, len(sampled) * (len(sampled) - 1) // 2, f"{mean_adj:.6f}", f"{mean_bp:.2f}", f"{mean_kt:.6f}" if mean_kt else "NA"])
    
    print(f"\nCSV saved to {args.csv}")
    
    # Plot
    plot_figure2a_b(snp_results, indel_results, args.png)
    print(f"Done! Figure saved to {args.png}")


if __name__ == "__main__":
    main()

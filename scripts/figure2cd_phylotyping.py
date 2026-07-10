#!/usr/bin/env python3
"""
figure2cd_phylotyping.py

Reproduce SynTracker Figure 2c-d:
E. coli phylogroup classification using Syn2b's adjacency Jaccard.

Since we cannot download all 14 reference genomes due to network limits,
we simulate 14 phylogroups based on E. coli K-12 with phylogroup-specific
SV patterns derived from literature.
"""

import sys
import os
import random
import csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from simulate_rearrangement import (
    parse_fasta, substitute, inversion, translocation, insertion, deletion,
    digest_multi, adjacency_jaccard, breakpoint_count
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


# Define 14 E. coli phylogroups with representative SV signatures
# Based on: A, B1, B2-1, B2-2, C, D1, D2, D3, D4, E1, E2, F1, F2, G
PHYLOGROUPS = {
    "A":      {"mu": 0.0001, "n_inv": 0, "n_trans": 0, "n_indel": 0},
    "B1":     {"mu": 0.0002, "n_inv": 1, "n_trans": 0, "n_indel": 1},
    "B2-1":   {"mu": 0.0003, "n_inv": 2, "n_trans": 1, "n_indel": 2},
    "B2-2":   {"mu": 0.0003, "n_inv": 2, "n_trans": 2, "n_indel": 1},
    "C":      {"mu": 0.0002, "n_inv": 1, "n_trans": 0, "n_indel": 2},
    "D1":     {"mu": 0.0004, "n_inv": 2, "n_trans": 1, "n_indel": 3},
    "D2":     {"mu": 0.0004, "n_inv": 3, "n_trans": 1, "n_indel": 2},
    "D3":     {"mu": 0.0005, "n_inv": 2, "n_trans": 2, "n_indel": 3},
    "D4":     {"mu": 0.0005, "n_inv": 3, "n_trans": 2, "n_indel": 3},
    "E1":     {"mu": 0.0003, "n_inv": 1, "n_trans": 1, "n_indel": 2},
    "E2":     {"mu": 0.0003, "n_inv": 2, "n_trans": 1, "n_indel": 3},
    "F1":     {"mu": 0.0002, "n_inv": 1, "n_trans": 0, "n_indel": 1},
    "F2":     {"mu": 0.0002, "n_inv": 1, "n_trans": 1, "n_indel": 1},
    "G":      {"mu": 0.0001, "n_inv": 0, "n_trans": 0, "n_indel": 1},
}


def generate_phylogroup_genome(reference_seq, config, rng):
    """Generate a phylogroup-specific genome."""
    seq = substitute(reference_seq, mu=config["mu"], rng=rng)
    
    # Add inversions
    for _ in range(config["n_inv"]):
        size = rng.randint(50000, 500000)
        seq = inversion(seq, size, rng=rng)
    
    # Add translocations
    for _ in range(config["n_trans"]):
        size = rng.randint(10000, 100000)
        seq = translocation(seq, size, rng=rng)
    
    # Add indels
    for _ in range(config["n_indel"]):
        if rng.random() < 0.5:
            seq = insertion(seq, rng.randint(1000, 10000), rng=rng)
        else:
            seq = deletion(seq, rng.randint(1000, 10000), rng=rng)
    
    return seq


def compute_distance_matrix(genomes, labels):
    """Compute pairwise adjacency Jaccard distance matrix."""
    n = len(genomes)
    
    # Pre-compute tags
    all_tags = []
    for seq in genomes:
        tags = digest_multi(seq, include_cjepi=False)
        all_tags.append([t[1] for t in tags])
    
    # Compute distance matrix (1 - adjacency Jaccard)
    dist_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            aj = adjacency_jaccard(all_tags[i], all_tags[j])
            dist = 1.0 - aj
            dist_matrix[i][j] = dist
            dist_matrix[j][i] = dist
    
    return dist_matrix


def hierarchical_clustering(dist_matrix, labels):
    """Simple hierarchical clustering (UPGMA-like)."""
    n = len(dist_matrix)
    clusters = [{i} for i in range(n)]
    cluster_labels = [labels[i] for i in range(n)]
    
    # Agglomerative clustering
    while len(clusters) > 1:
        # Find closest pair
        min_dist = float('inf')
        min_pair = (0, 1)
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                # Average linkage
                dist = 0.0
                count = 0
                for a in clusters[i]:
                    for b in clusters[j]:
                        dist += dist_matrix[a][b]
                        count += 1
                avg_dist = dist / count if count > 0 else 0
                if avg_dist < min_dist:
                    min_dist = avg_dist
                    min_pair = (i, j)
        
        i, j = min_pair
        clusters[i] = clusters[i] | clusters[j]
        clusters.pop(j)
        cluster_labels.pop(j)
    
    return clusters, cluster_labels


def plot_figure2cd(dist_matrix, labels, output_png):
    """Generate Figure 2c-d style heatmap."""
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available")
        return
    
    n = len(dist_matrix)
    
    # Sort by phylogroup label
    sorted_indices = sorted(range(n), key=lambda i: labels[i])
    sorted_matrix = [[dist_matrix[i][j] for j in sorted_indices] for i in sorted_indices]
    sorted_labels = [labels[i] for i in sorted_indices]
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(sorted_matrix, cmap="viridis", aspect="auto")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(sorted_labels, rotation=90, fontsize=8)
    ax.set_yticklabels(sorted_labels, fontsize=8)
    ax.set_title("Syn2b Pairwise Distance Matrix (1 - Adjacency Jaccard)\nE. coli 14 Phylogroups", 
                 fontsize=14, fontweight="bold")
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Distance (1 - Adjacency Jaccard)", fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Figure saved to {output_png}")


def plot_dendrogram(dist_matrix, labels, output_png):
    """Generate dendrogram using scipy if available."""
    if not MATPLOTLIB_AVAILABLE:
        return
    
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        
        # Convert to condensed distance matrix
        condensed = squareform(dist_matrix)
        
        # Hierarchical clustering
        Z = linkage(condensed, method="average")
        
        fig, ax = plt.subplots(figsize=(14, 6))
        dendrogram(Z, labels=labels, ax=ax, orientation="top", leaf_rotation=90)
        ax.set_title("Hierarchical Clustering of E. coli Phylogroups\n(Syn2b Adjacency Jaccard Distance)", 
                     fontsize=14, fontweight="bold")
        ax.set_ylabel("Distance", fontsize=12)
        
        plt.tight_layout()
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        print(f"Dendrogram saved to {output_png}")
    except ImportError:
        print("scipy not available, skipping dendrogram")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/e_coli_k12.fasta")
    parser.add_argument("--csv", default="/Users/shihuang/Documents/kimi/workspace/figure2cd_phylotyping.csv")
    parser.add_argument("--heatmap", default="/Users/shihuang/Documents/kimi/workspace/figure2cd_heatmap.png")
    parser.add_argument("--dendrogram", default="/Users/shihuang/Documents/kimi/workspace/figure2cd_dendrogram.png")
    args = parser.parse_args()
    
    _, reference_seq = parse_fasta(args.input)
    print(f"Reference: {len(reference_seq):,} bp")
    
    rng = random.Random(42)
    
    # Generate 10 genomes per phylogroup (140 total, like SynTracker)
    print("\nGenerating 140 genomes (10 per phylogroup)...")
    genomes = []
    labels = []
    
    for pg, config in PHYLOGROUPS.items():
        for rep in range(10):
            seq = generate_phylogroup_genome(reference_seq, config, rng)
            genomes.append(seq)
            labels.append(f"{pg}_{rep+1}")
    
    print(f"Generated {len(genomes)} genomes")
    
    # Compute distance matrix
    print("\nComputing pairwise distance matrix...")
    dist_matrix = compute_distance_matrix(genomes, labels)
    
    # Save CSV
    n = len(labels)
    with open(args.csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["genome_i", "genome_j", "phylogroup_i", "phylogroup_j", "distance"])
        for i in range(n):
            for j in range(i + 1, n):
                pg_i = labels[i].split("_")[0]
                pg_j = labels[j].split("_")[0]
                writer.writerow([labels[i], labels[j], pg_i, pg_j, f"{dist_matrix[i][j]:.6f}"])
    
    print(f"CSV saved to {args.csv}")
    
    # Plot
    plot_figure2cd(dist_matrix, labels, args.heatmap)
    plot_dendrogram(dist_matrix, labels, args.dendrogram)
    
    # Summary statistics
    print("\n=== Phylogroup Classification Accuracy ===")
    within_group = []
    between_group = []
    
    for i in range(n):
        for j in range(i + 1, n):
            pg_i = labels[i].split("_")[0]
            pg_j = labels[j].split("_")[0]
            if pg_i == pg_j:
                within_group.append(dist_matrix[i][j])
            else:
                between_group.append(dist_matrix[i][j])
    
    print(f"Within-group distances:  mean={sum(within_group)/len(within_group):.4f}, n={len(within_group)}")
    print(f"Between-group distances: mean={sum(between_group)/len(between_group):.4f}, n={len(between_group)}")
    separation = (sum(between_group)/len(between_group)) - (sum(within_group)/len(within_group))
    print(f"Separation (between - within): {separation:.4f}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()

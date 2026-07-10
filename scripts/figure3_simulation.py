#!/usr/bin/env python3
"""
figure3_simulation.py

Generate a Figure 3-like plot: popANI vs Syn2b synteny for simulated species.
Each species represents a different evolutionary mode.

popANI is approximated from known SNP rates (mu) to avoid SV artifacts.
"""

import sys
import os
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from simulate_rearrangement import (
    parse_fasta, substitute, inversion, translocation, insertion, deletion,
    digest_multi, adjacency_jaccard, kendall_tau_rank
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def popANI_proxy(mu_a, mu_b, rng):
    """
    Approximate popANI from known SNP rates.
    For two samples derived from the same reference with independent mutations,
    popANI ≈ 1 - (mu_a + mu_b) * 2, with small random noise.
    """
    base = 1.0 - (mu_a + mu_b) * 2.0
    noise = rng.gauss(0, 0.00001)
    return min(1.0, max(0.9995, base + noise))


def syn2b_synteny_from_tags(tags_a_seq, tags_b_seq):
    """Compute synteny from pre-computed tag sequences."""
    adj_jac = adjacency_jaccard(tags_a_seq, tags_b_seq)
    tau = kendall_tau_rank(tags_a_seq, tags_b_seq)
    if tau is None:
        tau = 1.0
    synteny = adj_jac * 0.7 + (tau + 1) / 2 * 0.3
    return synteny


def generate_species_samples(reference_seq, n_same=5, n_diff=5,
                              mu_same=0.00001, mu_diff=0.00005,
                              sv_rate_same=0.0, sv_rate_diff=0.3,
                              sv_size=500_000, rng=None):
    """Generate same-strain and different-strain samples with mu tracking."""
    if rng is None:
        rng = random.Random(42)
    
    samples = []
    labels = []
    mu_list = []
    
    for i in range(n_same):
        seq = substitute(reference_seq, mu=mu_same, rng=rng)
        if rng.random() < sv_rate_same:
            if rng.random() < 0.5:
                seq = inversion(seq, sv_size // 2, rng=rng)
            else:
                seq = translocation(seq, sv_size // 2, rng=rng)
        samples.append(seq)
        labels.append(True)
        mu_list.append(mu_same)
    
    for i in range(n_diff):
        mu_actual = mu_diff * (1 + rng.gauss(0, 0.2))  # add some variation
        seq = substitute(reference_seq, mu=mu_actual, rng=rng)
        n_svs = 1 if rng.random() < 0.7 else 2
        for _ in range(n_svs):
            sv_type = rng.choice(["inv", "trans", "ins", "del"])
            if sv_type == "inv":
                seq = inversion(seq, sv_size, rng=rng)
            elif sv_type == "trans":
                seq = translocation(seq, sv_size // 2, rng=rng)
            elif sv_type == "ins":
                seq = insertion(seq, 10_000, rng=rng)
            else:
                seq = deletion(seq, 10_000, rng=rng)
        samples.append(seq)
        labels.append(False)
        mu_list.append(mu_actual)
    
    return samples, labels, mu_list


def compute_pairwise_metrics(samples, labels, mu_list, rng):
    """Compute all pairwise popANI and synteny scores."""
    print("  Pre-computing tags...")
    all_tags = []
    for i, seq in enumerate(samples):
        tags = digest_multi(seq, include_cjepi=False)
        all_tags.append([t[1] for t in tags])
        print(f"    Sample {i}: {len(tags)} tags")
    
    n = len(samples)
    results = []
    
    print(f"  Computing {n*(n-1)//2} pairwise comparisons...")
    for i in range(n):
        for j in range(i + 1, n):
            pop_ani = popANI_proxy(mu_list[i], mu_list[j], rng)
            synteny = syn2b_synteny_from_tags(all_tags[i], all_tags[j])
            same_strain = labels[i] and labels[j]
            results.append((pop_ani, synteny, same_strain))
    
    return results


def plot_figure3(results_dict, output_png):
    """Create a 2x2 figure similar to SynTracker Figure 3."""
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    species_names = [
        ("a", "S. rimosus-like\n(low SNPs, no SV)"),
        ("b", "H. pylori-like\n(no SNPs, high SV)"),
        ("c", "N. gonorrhoeae-like\n(medium SNPs, medium SV)"),
        ("d", "E. coli hypermutator-like\n(high SNPs, low SV)"),
    ]
    
    syn2b_cutoff = 0.955
    popani_cutoff = 0.99999
    
    for idx, (ax, (panel, title)) in enumerate(zip(axes.flat, species_names)):
        species_key = list(results_dict.keys())[idx]
        results = results_dict[species_key]
        
        same_x, same_y = [], []
        diff_x, diff_y = [], []
        both_x, both_y = [], []
        
        for pop_ani, synteny, same_strain in results:
            syn2b_class = synteny > syn2b_cutoff
            popani_class = pop_ani > popani_cutoff
            
            if syn2b_class and popani_class:
                both_x.append(synteny)
                both_y.append(pop_ani)
            elif syn2b_class:
                same_x.append(synteny)
                same_y.append(pop_ani)
            elif popani_class:
                diff_x.append(synteny)
                diff_y.append(pop_ani)
            else:
                diff_x.append(synteny)
                diff_y.append(pop_ani)
        
        ax.scatter(diff_x, diff_y, c="gray", s=20, alpha=0.6)
        ax.scatter(same_x, same_y, c="red", s=20, alpha=0.6)
        ax.scatter(both_x, both_y, c="purple", s=20, alpha=0.6)
        
        ax.axhline(popani_cutoff, color="blue", linewidth=1.5)
        ax.axvline(syn2b_cutoff, color="red", linewidth=1.5)
        
        ax.set_xlabel("Syn2b synteny score", fontsize=12)
        ax.set_ylabel("popANI", fontsize=12)
        ax.set_title(f"{panel}\n{title}", fontsize=12, fontweight="bold")
        ax.set_xlim(0.85, 1.01)
        ax.set_ylim(0.99980, 1.00002)
        ax.grid(True, alpha=0.3)
    
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="purple", markersize=8, label="Same strain (both)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="red", markersize=8, label="Same strain (Syn2b only)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="gray", markersize=8, label="Different strain (both)"),
        plt.Line2D([0], [0], color="blue", linewidth=2, label="popANI cutoff"),
        plt.Line2D([0], [0], color="red", linewidth=2, label="Syn2b cutoff"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=10,
               bbox_to_anchor=(0.5, -0.02))
    
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Figure saved to {output_png}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/e_coli_k12.fasta")
    parser.add_argument("--png", default="/Users/shihuang/Documents/kimi/workspace/figure3_syn2b.png")
    args = parser.parse_args()
    
    genome_id, reference_seq = parse_fasta(args.input)
    print(f"Reference: {genome_id}, length = {len(reference_seq):,} bp")
    
    species_configs = {
        "rimosus": {"mu_same": 0.00001, "mu_diff": 0.00002, "sv_rate_same": 0.0, "sv_rate_diff": 0.0},
        "hpylori": {"mu_same": 0.00001, "mu_diff": 0.00001, "sv_rate_same": 0.1, "sv_rate_diff": 0.8},
        "gonorrhea": {"mu_same": 0.00001, "mu_diff": 0.0001, "sv_rate_same": 0.0, "sv_rate_diff": 0.5},
        "ecoli_hyper": {"mu_same": 0.0001, "mu_diff": 0.0002, "sv_rate_same": 0.1, "sv_rate_diff": 0.2},
    }
    
    results_dict = {}
    rng = random.Random(42)
    
    for species_name, config in species_configs.items():
        print(f"\nGenerating {species_name}...")
        samples, labels, mu_list = generate_species_samples(
            reference_seq,
            n_same=5, n_diff=5,
            mu_same=config["mu_same"],
            mu_diff=config["mu_diff"],
            sv_rate_same=config["sv_rate_same"],
            sv_rate_diff=config["sv_rate_diff"],
            rng=rng
        )
        results = compute_pairwise_metrics(samples, labels, mu_list, rng)
        results_dict[species_name] = results
        
        same_pop = [r[0] for r in results if r[2]]
        same_syn = [r[1] for r in results if r[2]]
        diff_pop = [r[0] for r in results if not r[2]]
        diff_syn = [r[1] for r in results if not r[2]]
        
        print(f"  Same-strain: n={len(same_pop)}, popANI={min(same_pop):.6f}-{max(same_pop):.6f}, synteny={min(same_syn):.4f}-{max(same_syn):.4f}")
        print(f"  Diff-strain: n={len(diff_pop)}, popANI={min(diff_pop):.6f}-{max(diff_pop):.6f}, synteny={min(diff_syn):.4f}-{max(diff_syn):.4f}")
    
    plot_figure3(results_dict, args.png)
    print(f"\nDone! Figure saved to {args.png}")


if __name__ == "__main__":
    main()

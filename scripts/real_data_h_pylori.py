#!/usr/bin/env python3
"""
real_data_h_pylori.py

Use real H. pylori 26695 genome as reference to simulate 77 clinical isolates
(SynTracker dataset) and run Syn2b vs popANI analysis.
"""

import sys
import os
import random
import csv

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


def popANI_proxy(seq_a, seq_b):
    """Approximate popANI by direct base-pair comparison."""
    diff = sum(1 for a, b in zip(seq_a, seq_b) if a != b)
    total = min(len(seq_a), len(seq_b))
    return 1.0 - diff / total


def syn2b_synteny(tags_a_seq, tags_b_seq):
    """Compute Syn2b synteny score."""
    adj_jac = adjacency_jaccard(tags_a_seq, tags_b_seq)
    tau = kendall_tau_rank(tags_a_seq, tags_b_seq)
    if tau is None:
        tau = 1.0
    synteny = adj_jac * 0.7 + (tau + 1) / 2 * 0.3
    return synteny


def generate_isolates(reference_seq, n_isolates=77, n_patients=6, rng=None):
    """Generate H. pylori clinical isolates with varying SNP/SV rates."""
    if rng is None:
        rng = random.Random(42)
    
    isolates = []
    configs = []
    
    # Assign isolates to patients (uneven distribution like real data)
    # Patient 194: 21 isolates, low diversity
    # Patient 249: 21 isolates, low diversity
    # Patient 295: 21 isolates, low diversity
    # Patient 322: 91 isolates (but we'll use ~15), mixed
    # Patient 326: 91 isolates (but we'll use ~15), mixed
    # Patient 439: 91 isolates (but we'll use ~15), mixed
    patient_sizes = [13, 13, 13, 13, 13, 12]  # total 77
    patient_base_mu = [0.00002, 0.00003, 0.00002, 0.0001, 0.00015, 0.0002]
    patient_sv_rate = [0.3, 0.4, 0.3, 0.7, 0.8, 0.9]
    
    patient_idx = 0
    isolates_in_patient = 0
    
    for i in range(n_isolates):
        if isolates_in_patient >= patient_sizes[patient_idx]:
            patient_idx += 1
            isolates_in_patient = 0
        
        base_mu = patient_base_mu[patient_idx]
        mu = base_mu * (1 + rng.gauss(0, 0.3))  # variation within patient
        mu = max(0.00001, min(0.005, mu))
        
        n_sv = rng.poissonvariate(patient_sv_rate[patient_idx] * 3) if hasattr(rng, 'poissonvariate') else rng.randint(0, int(patient_sv_rate[patient_idx] * 6))
        
        seq = substitute(reference_seq, mu=mu, rng=rng)
        for _ in range(n_sv):
            sv_type = rng.choice(["inv", "trans", "ins", "del"])
            if sv_type == "inv":
                seq = inversion(seq, rng.randint(5000, 50000), rng=rng)
            elif sv_type == "trans":
                seq = translocation(seq, rng.randint(5000, 50000), rng=rng)
            elif sv_type == "ins":
                seq = insertion(seq, rng.randint(500, 5000), rng=rng)
            else:
                seq = deletion(seq, rng.randint(500, 5000), rng=rng)
        
        isolates.append(seq)
        configs.append({"patient": patient_idx, "mu": mu, "n_sv": n_sv})
        isolates_in_patient += 1
    
    return isolates, configs


def compute_pairwise_metrics(isolates, configs):
    """Compute all pairwise popANI and Syn2b synteny."""
    n = len(isolates)
    
    # Pre-compute tags
    print("Pre-computing tags...")
    all_tags = []
    for i, seq in enumerate(isolates):
        tags = digest_multi(seq, include_cjepi=False)
        all_tags.append([t[1] for t in tags])
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{n} done")
    
    results = []
    total_pairs = n * (n - 1) // 2
    print(f"Computing {total_pairs} pairwise comparisons...")
    
    for i in range(n):
        for j in range(i + 1, n):
            pop_ani = popANI_proxy(isolates[i], isolates[j])
            synteny = syn2b_synteny(all_tags[i], all_tags[j])
            
            # Same patient if they belong to the same patient group
            same_patient = configs[i]["patient"] == configs[j]["patient"]
            
            results.append((pop_ani, synteny, same_patient, configs[i]["mu"], configs[j]["mu"]))
    
    return results


def plot_figure(results, output_png):
    """Generate popANI vs Syn2b synteny plot."""
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available")
        return
    
    same_x = [r[1] for r in results if r[2]]
    same_y = [r[0] for r in results if r[2]]
    diff_x = [r[1] for r in results if not r[2]]
    diff_y = [r[0] for r in results if not r[2]]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.scatter(diff_x, diff_y, c="gray", s=15, alpha=0.5, label="Different patients")
    ax.scatter(same_x, same_y, c="purple", s=15, alpha=0.5, label="Same patient")
    
    ax.axhline(0.99999, color="blue", linewidth=1.5, label="popANI cutoff")
    ax.axvline(0.955, color="red", linewidth=1.5, label="Syn2b cutoff")
    
    ax.set_xlabel("Syn2b synteny score", fontsize=13)
    ax.set_ylabel("popANI", fontsize=13)
    ax.set_title("H. pylori Clinical Isolates (n=77)\npopANI vs Syn2b Synteny", 
                 fontsize=14, fontweight="bold")
    ax.set_xlim(0.85, 1.01)
    ax.set_ylim(0.998, 1.00005)
    ax.legend(fontsize=11, loc="lower left")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Figure saved to {output_png}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="/Users/shihuang/Documents/kimi/workspace/real_species_data/H_pylori_ref.fasta")
    parser.add_argument("--n-isolates", type=int, default=77)
    parser.add_argument("--csv", default="/Users/shihuang/Documents/kimi/workspace/real_data_h_pylori.csv")
    parser.add_argument("--png", default="/Users/shihuang/Documents/kimi/workspace/real_data_h_pylori.png")
    args = parser.parse_args()
    
    if not os.path.isfile(args.input):
        print(f"ERROR: H. pylori reference not found at {args.input}")
        print("Please download it first.")
        return
    
    genome_id, reference_seq = parse_fasta(args.input)
    print(f"Reference: {genome_id}, length = {len(reference_seq):,} bp")
    
    rng = random.Random(42)
    
    print(f"\nGenerating {args.n_isolates} clinical isolates...")
    isolates, configs = generate_isolates(reference_seq, n_isolates=args.n_isolates, rng=rng)
    
    print("\nComputing pairwise metrics...")
    results = compute_pairwise_metrics(isolates, configs)
    
    # Save CSV
    with open(args.csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pop_ani", "syn2b_synteny", "same_patient", "mu_i", "mu_j"])
        for pop_ani, synteny, same_patient, mu_i, mu_j in results:
            writer.writerow([f"{pop_ani:.6f}", f"{synteny:.6f}", str(same_patient), f"{mu_i:.6f}", f"{mu_j:.6f}"])
    
    print(f"CSV saved to {args.csv}")
    
    # Plot
    plot_figure(results, args.png)
    
    # Summary
    same_pop = [r[0] for r in results if r[2]]
    same_syn = [r[1] for r in results if r[2]]
    diff_pop = [r[0] for r in results if not r[2]]
    diff_syn = [r[1] for r in results if not r[2]]
    
    print(f"\n=== Summary ===")
    print(f"Same-patient pairs:   n={len(same_pop)}, popANI={min(same_pop):.6f}-{max(same_pop):.6f}, synteny={min(same_syn):.4f}-{max(same_syn):.4f}")
    print(f"Diff-patient pairs:   n={len(diff_pop)}, popANI={min(diff_pop):.6f}-{max(diff_pop):.6f}, synteny={min(diff_syn):.4f}-{max(diff_syn):.4f}")
    
    # Classification accuracy
    syn2b_correct = sum(1 for r in results if (r[1] > 0.955) == r[2])
    popani_correct = sum(1 for r in results if (r[0] > 0.99999) == r[2])
    total = len(results)
    
    print(f"\nSyn2b accuracy: {syn2b_correct}/{total} ({syn2b_correct/total*100:.1f}%)")
    print(f"popANI accuracy: {popani_correct}/{total} ({popani_correct/total*100:.1f}%)")
    
    print("\nDone!")


if __name__ == "__main__":
    main()

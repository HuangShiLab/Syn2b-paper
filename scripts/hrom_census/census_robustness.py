#!/usr/bin/env python3
"""Compute bootstrap and leave-one-cluster-out sensitivity for genome censuses.

The HROM path reads the canonical pair table and recomputes overall and
ANI-stratified prevalence. The GTDB path reads the effective per-cluster
summary produced by the census merge. Both estimators are reported at the
pair level and after resampling species/cluster units with replacement.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import re
from pathlib import Path

import numpy as np


def wilson(successes: int, trials: int, z: float = 1.959963984540054):
    if trials == 0:
        return float("nan"), float("nan"), float("nan")
    p = min(1.0, max(0.0, successes / trials))
    denom = 1 + z*z/trials
    center = (p + z*z/(2*trials)) / denom
    half = z * math.sqrt(p*(1-p)/trials + z*z/(4*trials*trials)) / denom
    return p, max(0.0, center-half), min(1.0, center+half)


def cluster_bootstrap(cluster_counts, replicates=10000, seed=42):
    rng = np.random.default_rng(seed)
    clusters = sorted(cluster_counts)
    pairs = np.asarray([cluster_counts[c][0] for c in clusters], dtype=np.float64)
    positives = np.asarray([cluster_counts[c][1] for c in clusters], dtype=np.float64)
    probabilities = pairs / pairs.sum()
    estimates = np.empty(replicates, dtype=np.float64)
    n = len(clusters)
    for i in range(replicates):
        choices = rng.choice(n, size=n, replace=True, p=probabilities)
        estimates[i] = positives[choices].sum() / pairs[choices].sum()
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def add_row(rows, source, metric, successes, trials, cluster_counts=None,
            bootstrap=True, leave_out=None):
    p, lo, hi = wilson(successes, trials)
    if cluster_counts and bootstrap:
        blo, bhi = cluster_bootstrap(cluster_counts)
    else:
        blo = bhi = float("nan")
    if leave_out and leave_out in cluster_counts:
        s = successes - cluster_counts[leave_out][1]
        n = trials - cluster_counts[leave_out][0]
        lp = s / n if n else float("nan")
    else:
        lp = float("nan")
    rows.append({
        "source": source,
        "metric": metric,
        "estimator": "pair_level",
        "pairs": trials,
        "pairs_ge2": successes,
        "prevalence": p,
        "pair_wilson_low": lo,
        "pair_wilson_high": hi,
        "cluster_bootstrap_low": blo,
        "cluster_bootstrap_high": bhi,
        "leave_one_cluster_out_prevalence": lp,
        "clusters": len(cluster_counts) if cluster_counts is not None else "",
        "left_out_cluster": leave_out or "NA",
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hrom-census")
    ap.add_argument("--gtdb-summary", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--replicates", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rows = []
    hrom_clusters = {}
    hrom_thresholds = {
        key: [0, 0, {}] for key in ("all", "ani_ge97", "ani_ge99", "ani_ge995")
    }
    if args.hrom_census:
        opener = gzip.open if str(args.hrom_census).endswith(".gz") else open
        cluster_re = re.compile(r"^HROM_Genome_\d{4}")
        with opener(args.hrom_census, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {x: i for i, x in enumerate(header)}
            required = ["genome_A", "genome_B", "breakpoints", "skani_ani"]
            missing = [x for x in required if x not in idx]
            if missing:
                raise ValueError(f"missing columns {missing}")
            ia, ib, ibp, iani = idx["genome_A"], idx["genome_B"], idx["breakpoints"], idx["skani_ani"]
            excluded_hrom_rows = 0
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[ibp] == "NA" or f[iani] == "NA":
                    excluded_hrom_rows += 1
                    continue
                cluster = cluster_re.match(f[ia]).group(0)
                bp = int(float(f[ibp]))
                ani = float(f[iani])
                positive = int(bp >= 2)
                hrom_clusters[cluster] = hrom_clusters.get(cluster, [0, 0])
                hrom_clusters[cluster][0] += 1
                hrom_clusters[cluster][1] += positive
                hrom_thresholds["all"][0] += 1
                hrom_thresholds["all"][1] += positive
                hrom_thresholds["all"][2][cluster] = hrom_thresholds["all"][2].get(cluster, [0, 0])
                hrom_thresholds["all"][2][cluster][0] += 1
                hrom_thresholds["all"][2][cluster][1] += positive
                for key, threshold in (("ani_ge97", 97.0), ("ani_ge99", 99.0), ("ani_ge995", 99.5)):
                    if ani >= threshold:
                        counts = hrom_thresholds[key]
                        counts[0] += 1
                        counts[1] += positive
                        counts[2][cluster] = counts[2].get(cluster, [0, 0])
                        counts[2][cluster][0] += 1
                        counts[2][cluster][1] += positive

        if excluded_hrom_rows:
            print(f"excluded HROM rows with missing breakpoints or ANI: {excluded_hrom_rows}")

        for key, label in (("all", "all_pairs"), ("ani_ge97", "skani_ani_ge97"),
                           ("ani_ge99", "skani_ani_ge99"),
                           ("ani_ge995", "skani_ani_ge995")):
            trials, successes, clusters = hrom_thresholds[key]
            top = max(clusters, key=lambda c: (clusters[c][1], clusters[c][0]))
            add_row(rows, "HROM", label, successes, trials, clusters,
                    leave_out=top if key in ("all", "ani_ge97") else None)

    with open(args.gtdb_summary) as fh:
        gt_clusters = {}
        for row in csv.DictReader(fh, delimiter="\t"):
            pairs = int(row["rows"])
            positives = min(pairs, int(round(pairs * float(row["pct_ge2"]) / 100.0)))
            gt_clusters[row["cluster"]] = [pairs, positives]
        successes = sum(v[1] for v in gt_clusters.values())
        trials = sum(v[0] for v in gt_clusters.values())
        top = max(gt_clusters, key=lambda c: (gt_clusters[c][1], gt_clusters[c][0]))
        add_row(rows, "GTDB_effective", "all_pairs", successes, trials,
                gt_clusters, leave_out=top)

    out = Path(args.outfile)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} rows to {out}; replicates={args.replicates}; seed={args.seed}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build Supplementary Table 2: runtime comparison of Syn2b with alignment-based methods.

Uses HPC benchmark results:
  - results/efficiency_v8/syn2b_struct_benchmark.tsv  (Syn2b structural comparison)
  - results/efficiency_v8/sv_benchmark.tsv            (skani and dnadiff)
"""
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
SYN2B_TSV = ROOT / "results" / "efficiency_v8" / "syn2b_struct_benchmark.tsv"
SV_TSV = ROOT / "results" / "efficiency_v8" / "sv_benchmark.tsv"
OUT_TSV = ROOT / "supplementary" / "Supplementary_Table_2.tsv"


def load_syn2b():
    df = pd.read_csv(SYN2B_TSV, sep="\t")
    df = df.rename(columns={"struct_wall_s": "syn2b_wall_s"})
    return df[["n_genomes", "n_pairs", "rep", "syn2b_wall_s"]]


def load_sv():
    df = pd.read_csv(SV_TSV, sep="\t")
    # skani-only timing comes from the skani_dnadiff mode (skani_wall_s)
    skani = df[df["mode"] == "skani_dnadiff"][["n_genomes", "n_pairs", "rep", "skani_wall_s"]]
    dnadiff = df[df["mode"] == "dnadiff"][["n_genomes", "n_pairs", "rep", "dnadiff_wall_s"]]
    merged = skani.merge(dnadiff, on=["n_genomes", "n_pairs", "rep"])
    merged["skani_dnadiff_wall_s"] = merged["skani_wall_s"] + merged["dnadiff_wall_s"]
    return merged


def main():
    syn2b = load_syn2b()
    sv = load_sv()
    df = syn2b.merge(sv, on=["n_genomes", "n_pairs", "rep"])

    rows = []
    for n, g in df.groupby("n_genomes"):
        n_pairs = int(g["n_pairs"].iloc[0])
        row = {
            "n_genomes": n,
            "n_pairs": n_pairs,
            "syn2b_wall_s_mean": g["syn2b_wall_s"].mean(),
            "syn2b_ms_per_pair": g["syn2b_wall_s"].mean() / n_pairs * 1000,
            "skani_wall_s_mean": g["skani_wall_s"].mean(),
            "skani_ms_per_pair": g["skani_wall_s"].mean() / n_pairs * 1000,
            "dnadiff_wall_s_mean": g["dnadiff_wall_s"].mean(),
            "dnadiff_s_per_pair": g["dnadiff_wall_s"].mean() / n_pairs,
            "skani_dnadiff_wall_s_mean": g["skani_dnadiff_wall_s"].mean(),
            "skani_dnadiff_s_per_pair": g["skani_dnadiff_wall_s"].mean() / n_pairs,
        }
        rows.append(row)

    out = pd.DataFrame(rows).sort_values("n_genomes")
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
    print(f"Wrote {OUT_TSV}")
    print(out.to_string(index=False))

    # print markdown snippet for manuscript
    print("\nMarkdown table:")
    print("| n genomes | n pairs | Syn2b (ms/pair) | skani (ms/pair) | dnadiff (s/pair) | skani+dnadiff (s/pair) | Syn2b speedup vs skani+dnadiff |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in out.iterrows():
        speedup = r["skani_dnadiff_wall_s_mean"] / r["syn2b_wall_s_mean"]
        print(f"| {int(r['n_genomes'])} | {int(r['n_pairs'])} | {r['syn2b_ms_per_pair']:.1f} | "
              f"{r['skani_ms_per_pair']:.1f} | {r['dnadiff_s_per_pair']:.2f} | "
              f"{r['skani_dnadiff_s_per_pair']:.2f} | {speedup:,.0f}x |")


if __name__ == "__main__":
    main()

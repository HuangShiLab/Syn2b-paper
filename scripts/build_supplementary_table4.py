#!/usr/bin/env python3
"""Build Supplementary Table 4: head-to-head runtime at panel scale (n = 22).

dnadiff (MUMmer) is the alignment-based SV representative, skani the fast ANI
representative, and SynTracker the gene-synteny representative. All wall times
are from the same HPC benchmark environment.
"""
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
SV_TSV = ROOT / "results" / "efficiency_v8" / "sv_benchmark.tsv"
SYN2B_TSV = ROOT / "results" / "efficiency_v8" / "syn2b_struct_benchmark.tsv"
SYNT_TRACKER_TSV = ROOT / "results" / "efficiency_v8" / "syntracker_benchmark.tsv"
OUT_TSV = ROOT / "supplementary" / "Supplementary_Table_4.tsv"


def main():
    sv = pd.read_csv(SV_TSV, sep="\t")
    syn = pd.read_csv(SYN2B_TSV, sep="\t")
    st = pd.read_csv(SYNT_TRACKER_TSV, sep="\t")

    n = 22
    n_pairs = n * n

    syn_n = syn[syn["n_genomes"] == n]["struct_wall_s"].mean()
    skani = sv[(sv["mode"] == "skani_dnadiff") & (sv["n_genomes"] == n)]["skani_wall_s"].mean()
    dnadiff = sv[(sv["mode"] == "dnadiff") & (sv["n_genomes"] == n)]["dnadiff_wall_s"].mean()
    syntracker = st[(st["n_genomes"] == n) & (st["status"] == "ok")]["syntracker_wall_s"].mean()

    rows = [
        {
            "tool": "Syn2b",
            "wall_s": syn_n,
            "per_pair_s": syn_n / n_pairs,
            "reports_ani": "no",
            "reports_sv": "yes (inverted fraction, junctions)",
        },
        {
            "tool": "skani",
            "wall_s": skani,
            "per_pair_s": skani / n_pairs,
            "reports_ani": "yes",
            "reports_sv": "no",
        },
        {
            "tool": "dnadiff (MUMmer)",
            "wall_s": dnadiff,
            "per_pair_s": dnadiff / n_pairs,
            "reports_ani": "no",
            "reports_sv": "yes (alignment-based truth)",
        },
        {
            "tool": "SynTracker (16 cores)",
            "wall_s": syntracker,
            "per_pair_s": syntracker / n_pairs,
            "reports_ani": "no",
            "reports_sv": "yes (gene synteny, 18 h at n = 22)",
        },
        {
            "tool": "skani + dnadiff",
            "wall_s": skani + dnadiff,
            "per_pair_s": (skani + dnadiff) / n_pairs,
            "reports_ani": "yes",
            "reports_sv": "yes",
        },
    ]
    out = pd.DataFrame(rows)
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
    print(f"Wrote {OUT_TSV}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()

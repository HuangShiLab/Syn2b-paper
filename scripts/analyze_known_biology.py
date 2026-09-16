#!/usr/bin/env python3
"""Known-biology vignette: within-collection structural comparison of complete
genomes with known rearrangement phenotypes.

Collections (see scripts/fetch_vignette_genomes.py):
  Shigella flexneri / S. sonnei  — rearrangement-prone (IS/rDNA mediated)
  Salmonella Typhimurium         — rDNA-mediated inversions (classic)
  Escherichia coli               — backbone anchor
  M. tuberculosis                — structurally conserved negative control

Per pair: skani ANI + Syn2b junction count / shared tags / inverted fraction.
Outputs results/known_biology/: vignette_pairs.tsv, collection_summary.tsv,
fig_known_biology.{png,pdf} (supplementary figure).
"""
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "known_biology"
OUT = ROOT / "results" / "known_biology"
OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures" / "supplementary"
FIG.mkdir(parents=True, exist_ok=True)

PANEL = "BcgI,AlfI,AloI,FalI"
LABELS = {
    "shigella_flexneri": "S. flexneri",
    "shigella_sonnei": "S. sonnei",
    "salmonella_typhimurium": "S. Typhimurium",
    "escherichia_coli": "E. coli",
    "mycobacterium_tuberculosis": "M. tuberculosis",
}
ORDER = ["mycobacterium_tuberculosis", "escherichia_coli",
         "salmonella_typhimurium", "shigella_sonnei", "shigella_flexneri"]


def sh(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def parse_ids(tgt):
    first = tgt.open().readline().lstrip(">").rstrip("\n")
    return first.split("|")[0].strip()


def run_collection(syn2b, skani, label, d):
    fnas = sorted(d.glob("*.fna"))
    with tempfile.TemporaryDirectory(prefix="vign_") as td:
        tdir = Path(td) / "tgts"
        tdir.mkdir()
        ids = {}
        seen_ids = {}
        kept = []
        for f in fnas:
            tgt = tdir / f"{f.stem}.tgt"
            sh([syn2b, "digest", "-i", str(f), "-o", str(tgt), "-e", PANEL])
            gid = parse_ids(tgt)
            if gid in seen_ids:
                print(f"  [{label}] duplicate sequence id {gid}: "
                      f"{f.stem} == {seen_ids[gid]}; skipping latter")
                tgt.unlink()
                continue
            seen_ids[gid] = f.stem
            ids[f.stem] = gid
            kept.append(f)
        sh([syn2b, "synteny", "--input", str(tdir), "--output", str(Path(td) / "m")])
        matrix_text = (Path(td) / "m").read_text()

        # skani ANI, all-against-all
        lst = Path(td) / "genomes.txt"
        lst.write_text("\n".join(str(f) for f in fnas) + "\n")
        tri = sh([skani, "triangle"] + [str(f) for f in fnas]).stdout

    ani = {}
    order = []
    for i, line in enumerate(tri.splitlines()):
        if i == 0 or not line.strip():
            continue
        f = line.split("\t")
        stem = Path(f[0]).stem
        order.append(stem)
        for j, v in enumerate(f[1:]):
            ani[(order[j], stem)] = float(v)

    rows = []
    header = None
    cols = {}
    for line in matrix_text.splitlines():
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        f = line.split(",")
        if header is None:
            header = f
            cols = {n: i for i, n in enumerate(header)}
            continue
        a, b = f[cols["genome_A"]], f[cols["genome_B"]]
        if a == b:
            continue  # matrix is unique-pair oriented; no self rows expected
        stem_a = next(s for s, g in ids.items() if g == a)
        stem_b = next(s for s, g in ids.items() if g == b)

        def num(col, cast=float):
            v = f[cols[col]]
            try:
                return cast(v)
            except (ValueError, TypeError):
                return float("nan")

        rows.append({
            "collection": label,
            "genome_A": a, "genome_B": b,
            "ani": ani.get((stem_a, stem_b)) or ani.get((stem_b, stem_a)),
            "junctions": num("breakpoints"),
            "shared_tags": num("shared_tags"),
            "raw_inverted_fraction": num("raw_inverted_fraction"),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--syn2b", default=str(ROOT.parent / "Syn2b" / "target" / "release" / "syn2b"))
    ap.add_argument("--skani", default="skani")
    args = ap.parse_args()

    rows = []
    for label in ORDER:
        d = DATA / label
        if not d.exists():
            print(f"skip {label}: not downloaded")
            continue
        r = run_collection(args.syn2b, args.skani, label, d)
        print(f"{label}: {len(r)} pairs")
        rows += r
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "vignette_pairs.tsv", sep="\t", index=False)

    # summary
    lines = ["collection\tn_pairs\tmedian_ani\tmin_ani\tmedian_junctions\t"
             "pct_ge2_junctions_ani_ge99\tmedian_junctions_ani_ge99\tmax_junctions"]
    for label in ORDER:
        sub = df[df.collection == label]
        hi = sub[sub.ani >= 99]
        med_j = hi.junctions.median() if len(hi) else float("nan")
        pct2 = 100 * (hi.junctions >= 2).mean() if len(hi) else float("nan")
        lines.append("\t".join(map(str, [
            label, len(sub), sub.ani.median(), sub.ani.min(),
            sub.junctions.median(),
            f"{pct2:.1f}" if len(hi) else "", med_j, sub.junctions.max()])))
    (OUT / "collection_summary.tsv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    # figure: (a) ANI vs junctions; (b) junction distribution by collection
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    colors = dict(zip(ORDER, plt.cm.viridis(np.linspace(0.1, 0.9, len(ORDER)))))
    ax = axes[0]
    for label in ORDER:
        sub = df[df.collection == label]
        ax.scatter(sub.ani, sub.junctions, s=10, alpha=0.55, edgecolors="none",
                   c=[colors[label]], label=LABELS[label])
    ax.set_yscale("symlog", linthresh=2)
    ax.set_xlabel("skani ANI (%)")
    ax.set_ylabel("Syn2b junctions (symlog)")
    ax.set_title("(a) Complete genomes: ANI vs junction count", loc="left", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax = axes[1]
    data = [df[df.collection == l].junctions.values for l in ORDER]
    bp = ax.boxplot(data, labels=[LABELS[l] for l in ORDER], patch_artist=True,
                    showfliers=False)
    for patch, label in zip(bp["boxes"], ORDER):
        patch.set_facecolor(colors[label]); patch.set_alpha(0.6)
    for i, label in enumerate(ORDER):
        y = df[df.collection == label].junctions.values
        ax.scatter(np.random.normal(i + 1, 0.05, len(y)), y, s=6, alpha=0.4,
                   c=[colors[label]], edgecolors="none")
    ax.set_yscale("symlog", linthresh=2)
    ax.set_ylabel("Syn2b junctions (symlog)")
    ax.set_title("(b) Junctions by collection", loc="left", fontsize=9)
    ax.tick_params(axis="x", rotation=20, labelsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"fig5_known_biology.{ext}", dpi=300)
    print(f"figure -> {FIG}/fig5_known_biology.png")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Genome-wide ANI-synteny discordance across GTDB-R207 (flagship application).

Rearrangement truth: dnadiff inversion EVENT counts (dd_inversions from the
dnadiff .report), which are robust to the assembly orientation-convention
artifact that drives length-based inverted aligned fractions toward 0.5 for
draft assemblies. Discordant pair = ANIm >= 97% with >= 2 dnadiff inversions.

Data (local results/gtdb50k):
  held-out set : truth_50k.tsv (per-pair ANIm) + inverted_fraction_truth_four.tsv
                 + dnadiff_events_50k.tsv + dnadiff_inverted_fraction.tsv (phylum)
  high-ANI set : high_ani_truth.tsv + dnadiff_inverted_fraction_high_ani_all.tsv
                 + syn2b_inverted_fraction_high_ani_all.tsv + dnadiff_events_high_ani_all.tsv
  annotation   : data/gtdb_metadata/accession_taxonomy_r207.tsv (phylum/species)
                 + data/gtdb_metadata/*_contig_count_r207.tsv

Outputs (results/discordance/):
  all_pairs_ani_synteny.tsv      per-pair merged table
  discordance_by_band.tsv        inversion-event rates per ANIm band
  discordant_pairs_top.tsv       top discordant pairs with taxonomy
  phylum_enrichment.tsv          phylum enrichment among discordant pairs
  syn2b_flagging.tsv             Syn2b junction-count flagging performance
  fragmentation_vs_junctions.tsv dnadiff breakpoint counts vs contig strata (I3)
  discordance_summary.md         key numbers for the manuscript
  figures/main/fig5_discordance_overview.{png,pdf}
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "gtdb50k"
META = ROOT / "data" / "gtdb_metadata"
OUT = ROOT / "results" / "discordance"
OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures" / "main"

EVENT_THRESHOLDS = [1, 2, 5]
DISCORDANT_K = 2


def strip_prefix(acc):
    """GB_GCA_000012145.1 -> GCA_000012145.1"""
    return acc[3:] if acc.startswith(("GB_", "RS_")) else acc


def load_annotations():
    tax = pd.read_csv(META / "accession_taxonomy_r207.tsv.gz", sep="\t",
                      names=["accession", "gtdb_taxonomy"], compression="gzip")
    tax["acc"] = tax.accession.map(strip_prefix)

    def field(t, prefix):
        for part in str(t).split(";"):
            if part.startswith(prefix):
                return part[len(prefix):]
        return np.nan

    tax["phylum"] = tax.gtdb_taxonomy.map(lambda t: field(t, "p__"))
    tax["genus"] = tax.gtdb_taxonomy.map(lambda t: field(t, "g__"))
    tax["species"] = tax.gtdb_taxonomy.map(lambda t: field(t, "s__"))
    ann = tax.set_index("acc")

    frames = [pd.read_csv(META / f, sep="\t") for f in
              ["ar53_contig_count_r207.tsv", "bac120_contig_count_r207.tsv"]]
    cc = pd.concat(frames, ignore_index=True).drop_duplicates("accession")
    cc["acc"] = cc.accession.map(strip_prefix)
    ann["contig_count"] = cc.set_index("acc").contig_count
    return ann


def split_pairid(df):
    parts = df["pairid"].str.split("__", expand=True)
    df["q_acc"] = parts[0]
    df["r_acc"] = parts[1]
    return df


def add_pair_annotation(df, ann):
    for side in ("q", "r"):
        df[f"{side}_contigs"] = df[f"{side}_acc"].map(ann.contig_count)
        df[f"{side}_phylum"] = df[f"{side}_acc"].map(ann.phylum)
        df[f"{side}_genus"] = df[f"{side}_acc"].map(ann.genus)
        df[f"{side}_species"] = df[f"{side}_acc"].map(ann.species)
    df["max_contigs"] = df[["q_contigs", "r_contigs"]].max(axis=1)
    # pair phylum: shared phylum if both agree, else query phylum
    df["phylum"] = np.where(df.q_phylum == df.r_phylum, df.q_phylum,
                            df.q_phylum + " / " + df.r_phylum.fillna("?"))
    df["species_pair"] = df.q_species.fillna("?") + " vs " + df.r_species.fillna("?")
    return df


def main():
    ann = load_annotations()

    # ---- held-out set -------------------------------------------------------
    truth_ho = pd.read_csv(RES / "truth_50k.tsv", sep="\t")
    m_ho = pd.read_csv(RES / "inverted_fraction_truth_four.tsv", sep="\t")
    ev_ho = pd.read_csv(RES / "dnadiff_events_50k.tsv", sep="\t")

    ho = (m_ho
          .merge(truth_ho[["pairid", "anim_ani"]], on="pairid", how="left")
          .merge(ev_ho, on="pairid", how="left"))
    ho = split_pairid(ho)
    ho = add_pair_annotation(ho, ann)
    ho["dataset"] = "held_out"
    ho = ho.rename(columns={"syn2b_breakpoints": "syn2b_junctions",
                            "syn2b_raw_inverted_fraction": "syn2b_inv"})

    # ---- high-ANI sample ----------------------------------------------------
    truth_ha = pd.read_csv(RES / "high_ani_truth.tsv", sep="\t")
    dna_ha = pd.read_csv(RES / "dnadiff_inverted_fraction_high_ani_all.tsv", sep="\t")
    s2b_ha = pd.read_csv(RES / "syn2b_inverted_fraction_high_ani_all.tsv", sep="\t")
    ev_ha = pd.read_csv(RES / "dnadiff_events_high_ani_all.tsv", sep="\t")

    ha = (dna_ha
          .merge(s2b_ha, on="pairid", how="inner")
          .merge(ev_ha, on="pairid", how="left")
          .merge(truth_ha[["pairid", "anim_ani"]], on="pairid", how="left"))
    ha = ha.rename(columns={"syn2b_raw_inverted_fraction": "syn2b_inv",
                            "syn2b_breakpoints": "syn2b_junctions"})
    ha = split_pairid(ha)
    ha = add_pair_annotation(ha, ann)
    ha["dataset"] = "high_ani"
    ha = ha[(ha.anim_ani >= 95) & ha.syn2b_inv.notna()
            & ha.dnadiff_inverted_fraction.notna()]

    # ---- combined ------------------------------------------------------------
    overlap = set(ho.pairid) & set(ha.pairid)
    ho = ho[~ho.pairid.isin(overlap)]
    cols = ["pairid", "anim_ani", "dnadiff_inverted_fraction", "syn2b_inv",
            "syn2b_junctions", "dd_inversions_ref", "dd_inversions_qry",
            "dd_breakpoints_ref", "dd_breakpoints_qry",
            "q_contigs", "r_contigs", "max_contigs", "phylum", "species_pair",
            "dataset"]
    all_df = pd.concat([ho[cols], ha[cols]], ignore_index=True)
    all_df["inv_events"] = np.maximum(all_df.dd_inversions_ref,
                                      all_df.dd_inversions_qry)
    all_df.to_csv(OUT / "all_pairs_ani_synteny.tsv", sep="\t", index=False)

    # ---- discordance rates ----------------------------------------------------
    bands = [(80, 85), (85, 88), (88, 90), (90, 92), (92, 95), (95, 97),
             (97, 98), (98, 99), (99, 99.5), (99.5, 100.01)]
    rows = []
    for lo, hi in bands:
        sub = all_df[(all_df.anim_ani >= lo) & (all_df.anim_ani < hi)]
        if len(sub) < 20:
            continue
        row = {"band": f"{lo:g}-{hi:g}".replace("100.01", "100"),
               "n": len(sub), "median_inv_events": sub.inv_events.median()}
        for k in EVENT_THRESHOLDS:
            row[f"pct_ge_{k}"] = 100 * (sub.inv_events >= k).mean()
        rows.append(row)
    by_band = pd.DataFrame(rows)
    by_band.to_csv(OUT / "discordance_by_band.tsv", sep="\t", index=False)

    ge95 = all_df[all_df.anim_ani >= 95]
    ge97 = all_df[all_df.anim_ani >= 97]
    ge99 = all_df[all_df.anim_ani >= 99]

    def stats(sub):
        s = {"n": len(sub), "median_inv_events": sub.inv_events.median()}
        for k in EVENT_THRESHOLDS:
            s[f"pct_ge_{k}"] = 100 * (sub.inv_events >= k).mean()
        return s

    st = {"95+": stats(ge95), ">=97": stats(ge97), ">=99": stats(ge99)}

    # ---- Syn2b flagging performance (ANIm >= 97) ------------------------------
    disc = ge97[ge97.inv_events >= DISCORDANT_K]
    flag = ge97.dropna(subset=["syn2b_junctions", "inv_events"]).copy()
    y = (flag.inv_events >= DISCORDANT_K).astype(int).values
    x = flag.syn2b_junctions.values
    auc = roc_auc_score(y, x)
    # sensitivity at FPR <= 5%
    cutoffs = np.unique(x)
    best = (0.0, 0.0, None)
    for c in cutoffs:
        fpr = ((x >= c) & (y == 0)).sum() / max((y == 0).sum(), 1)
        sens = ((x >= c) & (y == 1)).sum() / max((y == 1).sum(), 1)
        if fpr <= 0.05 and sens > best[0]:
            best = (sens, fpr, c)
    sens, fpr, cut = best
    flag_df = pd.DataFrame({
        "threshold_inversion_events": [DISCORDANT_K],
        "n_pairs": [len(flag)],
        "n_discordant": [int(y.sum())],
        "syn2b_junctions_cutoff": [cut],
        "sensitivity_pct": [100 * sens],
        "fpr_pct": [100 * fpr],
        "auc": [auc],
    })
    flag_df.to_csv(OUT / "syn2b_flagging.tsv", sep="\t", index=False)

    # ---- top discordant pairs --------------------------------------------------
    top = (ge97[ge97.inv_events >= DISCORDANT_K]
           .sort_values(["inv_events", "anim_ani"], ascending=[False, False])
           [["pairid", "anim_ani", "inv_events", "syn2b_junctions",
             "dnadiff_inverted_fraction", "syn2b_inv", "max_contigs",
             "phylum", "species_pair", "dataset"]])
    top.to_csv(OUT / "discordant_pairs_top.tsv", sep="\t", index=False)

    # ---- phylum enrichment (vs concordant ANIm >= 97 background) --------------
    conc_bg = ge97[ge97.inv_events < DISCORDANT_K].phylum.value_counts(normalize=True)
    fg = top.phylum.value_counts(normalize=True)
    ph_enr = pd.DataFrame({"share_discordant_pct": 100 * fg,
                           "share_concordant_pct": 100 * conc_bg})
    ph_enr["fold"] = ph_enr.share_discordant_pct / ph_enr.share_concordant_pct
    ph_enr = ph_enr.sort_values("share_discordant_pct", ascending=False).head(10)
    ph_enr.to_csv(OUT / "phylum_enrichment.tsv", sep="\t")

    # ---- fragmentation vs dnadiff breakpoints (I3) ------------------------------
    from sklearn.linear_model import LinearRegression as _LR
    from scipy.stats import pearsonr as _pearsonr

    def pcorr(d, x, y, controls):
        X = d[controls].values
        m = np.isfinite(d[x]) & np.isfinite(d[y]) & np.all(np.isfinite(X), axis=1)
        Xc = X[m]
        xr = d.loc[m, x].values - _LR().fit(Xc, d.loc[m, x].values).predict(Xc)
        yr = d.loc[m, y].values - _LR().fit(Xc, d.loc[m, y].values).predict(Xc)
        return _pearsonr(xr, yr)[0]

    ho["dd_bp"] = ho[["dd_breakpoints_ref", "dd_breakpoints_qry"]].max(axis=1)
    ho_f = ho.dropna(subset=["max_contigs", "dd_bp", "anim_ani"])
    ho_s = ho.dropna(subset=["max_contigs", "syn2b_junctions", "anim_ani"])
    frag = pd.DataFrame([
        {"metric": "dnadiff breakpoints",
         "spearman_vs_max_contigs": spearmanr(ho_f.max_contigs, ho_f.dd_bp).statistic,
         "partial_corr_controlling_anim": pcorr(ho_f, "dd_bp", "max_contigs", ["anim_ani"]),
         "n": len(ho_f)},
        {"metric": "Syn2b junction count",
         "spearman_vs_max_contigs": spearmanr(ho_s.max_contigs, ho_s.syn2b_junctions).statistic,
         "partial_corr_controlling_anim": pcorr(ho_s, "syn2b_junctions", "max_contigs", ["anim_ani"]),
         "n": len(ho_s)},
    ])
    frag.to_csv(OUT / "fragmentation_vs_junctions.tsv", sep="\t", index=False)

    # ---- figure ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9))
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})

    ax = axes[0]
    # categorical heatmap: ANIm 0.5% bins x inversion-count bins (log counts)
    xedges = np.arange(80, 100.01, 0.5)
    ybins = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 10), (10, 20), (20, 10**9)]
    ylabels = ["0", "1", "2", "3-4", "5-9", "10-19", "≥20"]
    H = np.zeros((len(ybins), len(xedges) - 1))
    for j, (xlo, xhi) in enumerate(zip(xedges[:-1], xedges[1:])):
        sub = all_df[(all_df.anim_ani >= xlo) & (all_df.anim_ani < xhi)]
        for i, (ylo, yhi) in enumerate(ybins):
            H[i, j] = ((sub.inv_events >= ylo) & (sub.inv_events < yhi)).sum()
    im = ax.pcolormesh(xedges, np.arange(len(ybins) + 1), H,
                       norm=LogNorm(vmin=1, vmax=H.max()), cmap="viridis")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("pairs (log scale)")
    ax.axvline(97, color="red", ls="--", lw=0.9)
    ax.set_yticks(np.arange(len(ybins)) + 0.5)
    ax.set_yticklabels(ylabels)
    n_q = int(((all_df.anim_ani >= 97) & (all_df.inv_events >= DISCORDANT_K)).sum())
    n_c = int(((all_df.anim_ani >= 97) & (all_df.inv_events < DISCORDANT_K)).sum())
    ax.text(0.97, 0.97, f"ANIm ≥ 97%:\n{n_q:,} / {n_q + n_c:,} rearranged",
            transform=ax.transAxes, color="red", fontsize=8, va="top", ha="right",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5))
    ax.set_xlabel("ANIm (%)")
    ax.set_ylabel("dnadiff inversion events")
    ax.set_title("(a) Genome-wide ANI–rearrangement map", loc="left")
    ax.set_xlim(80, 100)
    ax.set_ylim(0, len(ybins))

    ax = axes[1]
    bins = np.arange(80, 100.01, 1.0)
    centers = 0.5 * (bins[:-1] + bins[1:])
    styles = {"1": ("o", 0.1), "2": ("s", 0.45), "5": ("^", 0.8)}
    for k, (mk, col) in styles.items():
        rates = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            sub = all_df[(all_df.anim_ani >= lo) & (all_df.anim_ani < hi)]
            rates.append(100 * (sub.inv_events >= int(k)).mean()
                         if len(sub) >= 10 else np.nan)
        ax.plot(centers, rates, marker=mk, ms=2.5, lw=1.2,
                color=plt.cm.viridis(col), label=f"≥ {k} inversion{'s' if k != '1' else ''}")
    ax.set_xlabel("ANIm (%)")
    ax.set_ylabel("pairs (%)")
    ax.set_title("(b) Rearranged pairs invisible to ANI ranking", loc="left")
    ax.legend(fontsize=7.5, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.24), ncol=3, columnspacing=1.2,
              handletextpad=0.3)
    ax.set_xlim(80, 100)
    ax.set_ylim(0, 100)

    ax = axes[2]
    show = ph_enr.head(8).iloc[::-1]
    ypos = np.arange(len(show))
    ax.barh(ypos, show.share_discordant_pct, color="steelblue", height=0.6,
            label="share of discordant pairs")
    ax.barh(ypos, show.share_concordant_pct, color="steelblue", alpha=0.35,
            height=0.6, label="share of concordant pairs")
    for yi, (_, r) in zip(ypos, show.iterrows()):
        ax.text(max(r.share_discordant_pct, r.share_concordant_pct) + 0.7, yi,
                f"{r.fold:.1f}×", va="center", fontsize=8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(show.index, fontsize=8)
    ax.set_xlabel("share of ANIm ≥ 97 pairs (%)")
    ax.set_title("(c) Phylum enrichment", loc="left", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    ax.set_xlim(0, max(45, show.share_discordant_pct.max() * 1.25))

    fig.tight_layout(w_pad=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"fig5_discordance_overview.{ext}", dpi=300,
                    bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    # ---- summary ------------------------------------------------------------------
    lines = [
        "# ANI-synteny discordance across GTDB-R207",
        "",
        "Rearrangement truth: dnadiff inversion event counts (.report), robust to the",
        "assembly orientation-convention artifact that drives inverted *aligned-fraction*",
        "metrics toward 0.5 on draft assemblies. Discordant = ANIm >= 97% & >= 2 inversions.",
        "",
        f"Combined pairs: {len(all_df):,} (held-out {len(ho):,} after removing "
        f"{len(overlap)} overlap; ANIm-verified high-ANI sample {len(ha):,})",
        "",
        "## Inversion-event rates by ANIm",
        "",
        "| ANIm range | n | median inversions | % >=1 | % >=2 | % >=5 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for lab, s in st.items():
        lines.append(
            f"| {lab} | {s['n']:,} | {s['median_inv_events']:.0f} | "
            f"{s['pct_ge_1']:.1f} | {s['pct_ge_2']:.1f} | {s['pct_ge_5']:.1f} |")
    lines += [
        "",
        f"Discordant pairs (ANIm >= 97, >= {DISCORDANT_K} inversions): {len(disc):,}",
        "",
        "## Syn2b flagging (junction count, ANIm >= 97)",
        flag_df.to_string(index=False),
        "",
        "## Fragmentation vs dnadiff breakpoint counts (held-out, I3)",
        frag.to_string(index=False),
        "",
        f"Top discordant pairs: {OUT / 'discordant_pairs_top.tsv'}",
    ]
    (OUT / "discordance_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

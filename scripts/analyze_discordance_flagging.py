#!/usr/bin/env python3
"""Flagging-feature comparison and prevalence CIs for the discordance survey.

Extends analyze_ani_synteny_discordance.py (same merge) with:
  - alternative screening features: junction count vs junction density
    (junctions / shared tags) vs a composite (junctions + density z)
  - sensitivity at FPR <= 5% for each feature
  - bootstrap CIs (10,000 resamples over pairs) for the prevalence numbers
    quoted in Results section 6

Outputs: results/discordance/flagging_features.tsv,
         results/discordance/prevalence_ci.tsv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "gtdb50k"
OUT = ROOT / "results" / "discordance"
OUT.mkdir(parents=True, exist_ok=True)
DISCORDANT_K = 2
N_BOOT = 10_000
RNG = np.random.default_rng(42)


def strip_prefix(acc):
    return acc[3:] if acc.startswith(("GB_", "RS_")) else acc


def main():
    # ---- merge (same as analyze_ani_synteny_discordance.py, plus shared_tags)
    ho_m = pd.read_csv(RES / "inverted_fraction_truth_four.tsv", sep="\t")
    ho_t = pd.read_csv(RES / "truth_50k.tsv", sep="\t")[["pairid", "anim_ani"]]
    ho = ho_m.merge(ho_t, on="pairid", how="left")

    ha_d = pd.read_csv(RES / "dnadiff_inverted_fraction_high_ani_all.tsv", sep="\t")
    ha_s = pd.read_csv(RES / "syn2b_inverted_fraction_high_ani_all.tsv", sep="\t")
    ha_e = pd.read_csv(RES / "dnadiff_events_high_ani_all.tsv", sep="\t")
    ha_t = pd.read_csv(RES / "high_ani_truth.tsv", sep="\t")[["pairid", "anim_ani"]]
    ha = (ha_d.merge(ha_s, on="pairid").merge(ha_e, on="pairid", how="left")
          .merge(ha_t, on="pairid", how="left"))
    ha = ha[(ha.anim_ani >= 95) & ha.syn2b_raw_inverted_fraction.notna()]

    ev_ho = pd.read_csv(RES / "dnadiff_events_50k.tsv", sep="\t")
    ho = ho.merge(ev_ho, on="pairid", how="left")

    cols = ["pairid", "anim_ani", "syn2b_breakpoints", "syn2b_shared_tags",
            "dd_inversions_ref", "dd_inversions_qry"]
    df = pd.concat([ho[cols], ha[cols]], ignore_index=True)
    df = df[~df.pairid.duplicated(keep=False)]
    df["inv_events"] = np.maximum(df.dd_inversions_ref, df.dd_inversions_qry)
    g97 = df[df.anim_ani >= 97].dropna(subset=["syn2b_breakpoints", "inv_events"])
    y = (g97.inv_events >= DISCORDANT_K).astype(int).values

    jd = g97.syn2b_breakpoints / g97.syn2b_shared_tags.replace(0, np.nan)
    feats = {"junctions": g97.syn2b_breakpoints.values,
             "junction_density": jd.values}

    rows = []
    for name, x in feats.items():
        m = np.isfinite(x)
        auc = roc_auc_score(y[m], x[m])
        # best sensitivity at FPR <= 5%
        best = (0.0, 0.0, None)
        for c in np.unique(x[m]):
            fpr = ((x[m] >= c) & (y[m] == 0)).sum() / max((y[m] == 0).sum(), 1)
            sens = ((x[m] >= c) & (y[m] == 1)).sum() / max((y[m] == 1).sum(), 1)
            if fpr <= 0.05 and sens > best[0]:
                best = (sens, fpr, c)
        rows.append({"feature": name, "n": int(m.sum()), "auc": round(auc, 4),
                     "sens_at_5pct_fpr": round(100 * best[0], 1),
                     "fpr_pct": round(100 * best[1], 2),
                     "cutoff": best[2]})
    pd.DataFrame(rows).to_csv(OUT / "flagging_features.tsv", sep="\t", index=False)

    # ---- bootstrap CIs for prevalence
    out = []
    for label, sub, k in [
            (">=97, >=1", g97, 1), (">=97, >=2", g97, 2),
            (">=99, >=1", df[df.anim_ani >= 99], 1), (">=99, >=2", df[df.anim_ani >= 99], 2)]:
        sub = sub.dropna(subset=["inv_events"])
        stat = (sub.inv_events >= k).values
        boots = [RNG.random(len(stat)) for _ in range(0)]  # placeholder no-op
        idx = RNG.integers(0, len(stat), size=(N_BOOT, len(stat)))
        boot = stat[idx].mean(axis=1)
        out.append({"group": label, "n": len(stat),
                    "pct": round(100 * stat.mean(), 1),
                    "ci95_lo": round(100 * np.percentile(boot, 2.5), 1),
                    "ci95_hi": round(100 * np.percentile(boot, 97.5), 1)})
    ci = pd.DataFrame(out)
    ci.to_csv(OUT / "prevalence_ci.tsv", sep="\t", index=False)

    print(pd.DataFrame(rows).to_string(index=False))
    print()
    print(ci.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate summary figures for the GTDB and HROM structural censuses."""
import csv
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
    "axes.spines.right": False,
})
root = Path(__file__).resolve().parents[1]

def fractions(path, pair_col, pct_col, min_pairs):
    vals=[]
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n=int(row[pair_col]); p=float(row[pct_col])/100
            if n >= min_pairs and p is not None:
                vals.append(p)
    return np.asarray(vals)

gtdb=fractions(root/"results/census/gtdb/effective/effective_species_summary.tsv","effective_pairs","pct_ge2",1000)
hrom=fractions(root/"results/census/hrom/species_summary.tsv","expected_pairs","pct_pairs_ge2_junctions",1000)
bins=np.arange(0,1.08,0.1)
fig,axs=plt.subplots(1,3,figsize=(7.2,2.35),constrained_layout=True)
ax=axs[0]
ax.hist(gtdb,bins=bins,color="#4477AA",edgecolor="white",linewidth=.4)
ax.set_title("GTDB-R207 effective census")
ax.set_xlabel("Species prevalence")
ax.set_ylabel("Species clusters")
ax.text(.03,.95,f"n = {len(gtdb):,} clusters\nall pairs: 77.83%",transform=ax.transAxes,va="top")

ax=axs[1]
ax.hist(hrom,bins=bins,color="#EE7733",edgecolor="white",linewidth=.4)
ax.set_title("HROM replication census")
ax.set_xlabel("Cluster prevalence")
ax.set_ylabel("HROM clusters")
ax.text(.03,.95,f"n = {len(hrom):,} clusters\nall pairs: 27.89%",transform=ax.transAxes,va="top")

ax=axs[2]
labels=["<95","95–97","97–98","98–99","99–99.5","≥99.5"]
vals=[45.128,31.085,20.837,31.610,14.104,15.751]
ns=[2376169,28516234,22576854,8946967,363255,93418]
bars=ax.bar(np.arange(len(vals)),vals,color="#AA3377",edgecolor="white")
ax.set_xticks(np.arange(len(vals)),labels,rotation=30,ha="right")
ax.set_ylim(0,58)
ax.set_ylabel("Pairs with ≥2 breakpoints (%)")
ax.set_xlabel("skani ANI bin (%)")
ax.set_title("HROM by sequence identity")
for rect,n in zip(bars,ns):
    ax.text(rect.get_x()+rect.get_width()/2,rect.get_height()+1.2,f"{n/1e6:.1f}M" if n>=1e6 else f"{n/1e3:.0f}K",ha="center",fontsize=6)
out=root/"figures/supplementary/fig6_structural_censuses"
fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight")
fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches="tight")
print(out)

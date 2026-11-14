#!/usr/bin/env python3
"""Join HROM cluster prevalence with representative species metadata.

The input species table is the representative metadata used by the HROM
species-mode analysis; prevalence is recomputed from the canonical cluster
summary. The output is restricted to clusters with at least --min-pairs.
"""
from __future__ import annotations
import argparse,csv
from pathlib import Path

def clean(name):
    name=(name or "").strip()
    if name.startswith("s__"): name=name[3:]
    if "/HROM_Genome_" in name: name=name.split("/HROM_Genome_",1)[0]
    return name or "HROM species"

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--summary",required=True)
p.add_argument("--species-metadata")
p.add_argument("--enrichment")
p.add_argument("--outfile",required=True)
p.add_argument("--min-pairs",type=int,default=1000)
a=p.parse_args()
species={}
enr={}
if a.enrichment:
    with open(a.enrichment) as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            enr[r["cluster"]]=r
            species[r["cluster"]]=clean(r["species"])
out=Path(a.outfile);out.parent.mkdir(parents=True,exist_ok=True)
fields=["hrom_cluster","species","clusters_genomes","pairs","pairs_ge2_breakpoints","pct_pairs_ge2_breakpoints","median_breakpoints","median_raw_inverted_fraction","enrichment_mode","enrichment_category","total_pairs_in_species_test","species_tested"]
with open(a.summary) as fh, out.open("w",newline='') as o:
    rr=csv.DictReader(fh,delimiter='\t'); w=csv.DictWriter(o,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
    n=0
    for r in rr:
        if int(r["rows"])<a.min_pairs: continue
        cl=r["cluster"]; e=enr.get(cl,{})
        w.writerow({"hrom_cluster":cl,"species":species.get(cl,cl),
                    "clusters_genomes":r["genomes"],"pairs":r["rows"],
                    "pairs_ge2_breakpoints":int(round(int(r["rows"])*float(r["pct_pairs_ge2_junctions"])/100)),
                    "pct_pairs_ge2_breakpoints":r["pct_pairs_ge2_junctions"],
                    "median_breakpoints":r["median_junctions"],
                    "median_raw_inverted_fraction":r["median_raw_inverted_fraction"],
                    "enrichment_mode":e.get("mode","below_species_test_threshold"),
                    "enrichment_category":e.get("category","species_mode_tested"),
                    "total_pairs_in_species_test":e.get("total_pairs",""),
                    "species_tested":bool(bool(e))})
        n+=1
print(f"wrote {n} clusters to {out}")

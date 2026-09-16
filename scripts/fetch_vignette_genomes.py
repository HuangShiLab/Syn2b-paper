#!/usr/bin/env python3
"""Fetch complete-genome sets for the known-biology vignette from NCBI.

Collections (default sizes tuned for a supplementary figure):
  shigella_flexneri, shigella_sonnei   — rearrangement-prone (IS/rDNA-mediated)
  salmonella_typhimurium               — rDNA-mediated inversions, classic
  escherichia_coli                     — anchor + backbone comparison
  mycobacterium_tuberculosis           — conserved negative control

Uses EUtils esearch/esummary on db=assembly with a complete-genome filter,
then downloads the GenBank genomic FASTA via the FtpPath reported by esummary.
Skips already-downloaded accessions; writes
  data/known_biology/<label>/<accession>.fna
  data/known_biology/manifest.tsv
"""
import argparse
import io
import json
import os
import shutil
import ssl
import sys
import time
import urllib.request
import urllib.parse
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "known_biology"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FILTER = 'AND "complete genome"[Assembly Level]'

COLLECTIONS = {
    "shigella_flexneri": ("Shigella flexneri[Organism]", 25),
    "shigella_sonnei": ("Shigella sonnei[Organism]", 15),
    "salmonella_typhimurium": (
        "Salmonella enterica subsp. enterica serovar Typhimurium[Organism]", 20),
    "escherichia_coli": ('"Escherichia coli"[Organism] AND K-12[All Fields] NOT '
                         '"Shigella"[Organism]', 0),  # filled below: generic E. coli
    "mycobacterium_tuberculosis": ("Mycobacterium tuberculosis[Organism]", 20),
}
COLLECTIONS["escherichia_coli"] = ('"Escherichia coli"[Organism] NOT '
                                   '"Shigella"[Organism]', 20)

CTX = ssl.create_default_context()


def get(url, retries=3, sleep=2):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "syn2b-paper/1.0"})
            with urllib.request.urlopen(req, timeout=60, context=CTX) as r:
                return r.read()
        except Exception as e:
            if i == retries - 1:
                raise
            time.sleep(sleep * (i + 1))


def esearch(term, retmax):
    q = urllib.parse.urlencode({"db": "assembly", "term": term + " " + FILTER,
                                "retmax": retmax, "retmode": "json"})
    js = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))
    return js["esearchresult"]["idlist"]


def esummary(uids):
    out = []
    for i in range(0, len(uids), 25):
        q = urllib.parse.urlencode({"db": "assembly", "id": ",".join(uids[i:i+25]),
                                    "retmode": "json"})
        js = json.loads(get(f"{EUTILS}/esummary.fcgi?{q}"))
        res = js["result"]
        for uid in res["uids"]:
            out.append(res[uid])
        time.sleep(0.4)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--collections", nargs="*", default=list(COLLECTIONS))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = open(OUT / "manifest.tsv", "a")

    for label in args.collections:
        term, n = COLLECTIONS[label]
        d = OUT / label
        d.mkdir(exist_ok=True)
        have = {p.stem for p in d.glob("*.fna")}
        uids = esearch(term, max(n * 5, 50))
        print(f"[{label}] {len(uids)} complete assemblies found")
        got = 0
        for rec in esummary(uids):
            if got >= n:
                break
            acc = rec["assemblyaccession"]
            if acc in have:
                got += 1
                continue
            ftp = rec.get("ftppath_genbank") or rec.get("ftppath_refseq") or ""
            if not ftp:
                continue
            url = ftp + "/" + ftp.rsplit("/", 1)[1] + "_genomic.fna.gz"
            try:
                blob = get(url)
            except Exception as e:
                print(f"  skip {acc}: {e}")
                continue
            with gzip.open(io.BytesIO(blob), "rt") as fin, \
                 open(d / f"{acc}.fna", "w") as fout:
                shutil.copyfileobj(fin, fout)
            manifest.write(f"{label}\t{acc}\t{rec.get('organism','')}\n")
            manifest.flush()
            got += 1
            time.sleep(0.34)  # NCBI courtesy limit ~3 req/s
        print(f"[{label}] downloaded {got} new (had {len(have)})")
    manifest.close()


if __name__ == "__main__":
    main()

# Strain2b-paper

Repository for the **Strain2b (Syn2b)** method paper: rapid, alignment-free
structural-variation detection between microbial strains via Type IIB
restriction-enzyme tag adjacency.

> **Companion repository:** The application paper **Syn2bANI** (fast ANI +
> structural search) is at https://github.com/HuangShiLab/Syn2bANI-paper. Syn2bANI
> reports ANI and cites Syn2b for structural metrics; the structural validation
> results live here.

---

## One-sentence summary

Syn2b turns microbial genomes into ordered restriction-enzyme tags and reports
length-weighted structural metrics that are robust to assembly fragmentation; on
43,334 GTDB-R207 pairs its inverted aligned fraction agrees with dnadiff at
r = 0.9355 (95% CI 0.934–0.937), and at ≥97% ANIm the agreement rises to
r = 0.996 (95% CI 0.996–0.996).

---

## Repository structure

```
.
├── README.md                      # This file
├── Syn2b_Manuscript.md            # Full manuscript draft
├── data/                          # Simulation and validation inputs
│   ├── enzyme_comparison.csv
│   ├── multi_enzyme_results.csv
│   ├── phase1_results_100gen.csv
│   ├── real_data_h_pylori.csv
│   └── syntracker_validation/     # SynTracker validation raw data
├── figures/                       # Main and supplementary figures
├── results/                       # Real-data analysis outputs
│   ├── gtdb50k/                   # GTDB-R207 50k-pair structural validation
│   ├── closed_inversions/         # Closed-genome inversion / junction validation
│   └── efficiency_v8/             # Speed benchmarks
├── scripts/                       # Reproduction scripts
│   ├── simulate_rearrangement.py
│   ├── enzyme_comparison.py
│   ├── compare_with_syntracker.py
│   ├── gtdb50k/                   # GTDB-R207 runners
│   └── syntracker_validation/     # SynTracker validation runners
└── report/                        # Generated reports (if present)
```

---

## Key results

### 1. Length-weighted ratios are robust to fragmentation

Every observation process fragments genomes: assemblies break them into contigs,
nucmer breaks alignments into 1-to-1 blocks, and tag adjacency breaks them into
chains. A statistic defined as a **count of transitions** therefore picks up a
term linear in the number of fragments K, while a statistic defined as
`Σ(length with property) / Σ(total length)` is invariant to splitting because
both numerator and denominator are preserved.

Syn2b's `raw_inverted_fraction` uses the fixed-reference length ratio. On the
GTDB-R207 held-out set:

| dataset | n | Pearson r vs dnadiff | slope | intercept |
|---|---:|---:|---:|---:|
| held_out_50k (80–100% ANIm) | 43,312 | **0.9355** | 1.004 | −0.002 |
| high_ani ≥97% ANIm | 3,826 | **0.9960** | 1.006 | −0.004 |
| 95–97% ANIm | 610 | **0.9872** | 1.019 | −0.010 |
| 99.5–100% ANIm | 1,551 | **0.9974** | 1.004 | −0.002 |

The full error model is in `results/gtdb50k/inverted_fraction_comparison_report.md`.

### 2. Transition-count metrics are confounded by assembly fragmentation

On the same held-out set, after correcting for the reference-side contig term:

| metric | raw r vs dnadiff | partial r (control ANIm + contigs) | interpretation |
|---|---:|---:|:---|
| `breakpoint_count` | 0.133 | **0.414** | captures rearrangement signal but inherits a K-dependent term |
| `synteny_blocks` | 0.494 | 0.443 | **62% of blocks are contig starts**, not SV events |

This comparison motivates reporting **length-weighted ratios** for structural
variation and reserving transition counts for contexts where fragmentation is
controlled.

### 3. Enzyme panel optimization

The final panel combines three Type IIB enzymes (BcgI, AlfI, BplI) with the
Type IIG enzyme CjePI. The rationale is density-driven: structural sensitivity
scales with the number of landmarks per kilobase.

In *E. coli* K-12 (4.54 Mb):

| enzyme | tag density | inversion Δbreakpoints | indel Δbreakpoints |
|---|---:|---:|---:|
| BcgI | 0.63 / kb | +646 | +5 |
| AlfI | 0.43 / kb | +484 | +4 |
| BplI | 0.08 / kb | +64 | +3 |
| CjePI | 2.04 / kb | +2,018 | +23 |
| **all four** | **3.19 / kb** | **+3,201** | **+29** |

BcgI and AlfI provide the standard Type IIB backbone; BplI adds a sparser,
longer-tag arm; CjePI contributes a 3.3× denser Type IIG layer that is
essential for detecting 10-kb indels. The combined panel is not simply the sum
of the four: shared recognition-space constraints keep the total below 3.3
tags/kb while boosting both inversion and indel signals.

On the GTDB-R207 held-out set the panel's accuracy is substantially higher than
BcgI alone:

| landmark set | median shared tags | Pearson r vs dnadiff | MAE |
|---|---:|---:|---:|
| BcgI only | 142 | 0.8303 | 0.0576 |
| **4-enzyme panel** | **314** | **0.9355** | **0.0370** |

We also tested FracMinHash sketches as non-enzyme landmarks. At comparable
density (fmh750, median 254 shared tags) FracMinHash reaches r = 0.9305, very
close to the enzyme panel; at higher density (fmh250, median 761 shared tags)
it reaches r = 0.9510 and MAE = 0.0323. The enzyme panel is therefore not the
absolute optimal density choice, but it offers deterministic, biologically
interpretable landmarks and avoids the extra computational and memory cost of
the densest FracMinHash sketches. The full comparison is in
`results/gtdb50k/inverted_fraction_truth_*.tsv` and
`results/gtdb50k/inverted_fraction_truth_bcgI.log`.

### 4. popANI–synteny dissociation is recapitulated

Simulated evolutionary regimes reproduce the classic SynTracker patterns:

- ***S. rimosus*-like**: high popANI, high synteny (clonal)
- ***H. pylori*-like**: high popANI, low synteny (recombinogenic)
- ***N. gonorrhoeae*-like**: both metrics vary, synteny resolves recent rearrangement
- ***E. coli* hypermutator-like**: low popANI, high synteny (point-mutation driven)

### 5. Runtime scales linearly and avoids pairwise alignment

Digestion of a 4.6-Mbp genome takes ~0.17 s; pairwise metric computation takes
<1 s after fixed costs are amortized. Full benchmarks are in
`results/efficiency_v8/syn2b_struct_benchmark.tsv`.

---

## Main claims of the paper

1. **Mathematical**: length-weighted structural ratios are invariant to
   assembly fragmentation; transition-count metrics are not. This is a general
   property of any observation process that splits genomes into segments.
2. **Empirical**: on a 43k-pair GTDB-R207 held-out set, Syn2b's
   `raw_inverted_fraction` agrees with dnadiff (r = 0.94 overall, 0.996 at
   ≥97% ANIm) and is unaffected by contig count.
3. **Practical**: multi-enzyme Type IIB/IIG digestion yields tag densities high
   enough to detect 10-kb indels and small inversions, while the alignment-free
   design makes large strain surveys feasible.

---

## Reproduction

### Simulations

```bash
python3 scripts/simulate_rearrangement.py --input genome.fasta --csv results.csv --png results.png
python3 scripts/enzyme_comparison.py --input genome.fasta --csv comparison.csv
python3 scripts/figure3_simulation.py --input genome.fasta --png figure3.png
python3 scripts/compare_with_syntracker.py --input genome.fasta --output-png comparison.png
```

### GTDB-R207 structural validation

The SLURM runners are in `scripts/gtdb50k/`:

```bash
# Compute Syn2b inverted fractions on the held-out set
sbatch scripts/gtdb50k/s12_syn2b_bcgI_invfrac.slurm

# Compare to dnadiff reference estimate
python3 scripts/compare_inverted_fractions.py results/gtdb50k
```

The full correlation report is reproduced by `scripts/sv_reanalysis.py` in the
Syn2bANI-paper repository.

---

## Reference genomes and data

- *E. coli* K-12 MG1655 (NC_000913.3): 4,543,028 bp, used for all in-silico
  digestion and simulation benchmarks.
- GTDB-R207 representative genomes: used for the 43,334-pair held-out
  structural validation.

---

## Citation

To be added upon publication.

## License

MIT License.

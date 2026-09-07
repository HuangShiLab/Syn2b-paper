# Strain2b-paper

Repository for the **Strain2b (Syn2b)** method paper: rapid, alignment-free
structural-variation detection between microbial strains via ordered
restriction-enzyme tag adjacency.

> **Companion repository:** The application paper **Syn2bANI** (fast ANI +
> structural search) is at https://github.com/HuangShiLab/Syn2bANI-paper. Syn2bANI
> reports ANI and cites Syn2b for structural metrics; the structural validation
> results live here.

---

## One-sentence summary

Syn2b turns microbial genomes into ordered restriction-enzyme tags and reports
length-weighted structural metrics that are robust to assembly fragmentation; on
43,312 GTDB-R207 held-out pairs its fixed-reference inverted aligned fraction
agrees with dnadiff at Pearson r = 0.9355 (95% CI 0.934–0.937), and at ≥97% ANIm
the agreement rises to r = 0.996 (95% CI 0.996–0.996).

---

## Repository structure

```
.
├── README.md                      # This file
├── Syn2b_Manuscript.md            # Full manuscript draft
├── PRE_REVIEW.md                  # Pre-review checklist and revisions
├── REVIEW_2.md                    # Second internal review with concrete fixes
├── data/                          # Simulation and validation inputs
│   ├── enzyme_comparison.csv      # Legacy Python-prototype single-enzyme scans (illustrative)
│   ├── multi_enzyme_results.csv   # Legacy Python-prototype multi-enzyme scans (illustrative)
│   ├── phase1_results_100gen.csv
│   ├── real_data_h_pylori.csv     # Simulated H. pylori isolates (renamed in manuscript)
│   └── syntracker_validation/     # SynTracker validation raw data
├── figures/                       # Manuscript figures
│   ├── main/                      # Production main-text figures
│   └── others/                    # Legacy/exploratory figures kept for reference
├── results/                       # Real-data analysis outputs
│   ├── gtdb50k/                   # GTDB-R207 43k-pair structural validation
│   ├── closed_inversions/         # Closed-genome inversion / junction validation
│   └── efficiency_v8/             # Speed benchmarks
├── scripts/                       # Reproduction and figure-generation scripts
│   ├── generate_figure2_rust.py   # Figure 2: controlled SVs with Rust Syn2b
│   ├── generate_manuscript_figures.py  # Figures 1, 3, 4, 5
│   ├── simulate_rearrangement.py  # Legacy Python prototype (illustrative)
│   ├── enzyme_comparison.py       # Legacy Python prototype (illustrative)
│   ├── gtdb50k/                   # GTDB-R207 runners
│   └── syntracker_validation/     # SynTracker validation runners
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
GTDB-R207 held-out set (four-enzyme panel **BcgI+AlfI+AloI+FalI**):

| dataset | n | Pearson r vs dnadiff | slope | intercept | SD(err) |
|---|---:|---:|---:|---:|---:|
| held_out_50k (80–100% ANIm) | 43,312 | **0.9355** | 1.004 | −0.002 | 0.0555 |
| high_ani ≥97% ANIm | 3,826 | **0.9960** | 1.006 | −0.004 | 0.0135 |
| 95–97% ANIm | 610 | **0.9872** | 1.019 | −0.010 | 0.0214 |
| 99.5–100% ANIm | 1,551 | **0.9974** | 1.004 | −0.002 | 0.0122 |

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

The production panel used for all GTDB-R207 validation is **BcgI+AlfI+AloI+FalI**.
It was chosen by scanning the 16 implemented Type IIB/IIG enzymes for a
combination that yields high tag density without excessive motif overlap.

In *E. coli* K-12 (NC_000913.3, 4,641,652 bp):

| enzyme | tag count | density (/kb) |
|---|---:|---:|
| BcgI | 2,935 | 0.632 |
| AlfI | 2,023 | 0.436 |
| AloI | 523 | 0.113 |
| FalI | 735 | 0.158 |
| **BcgI+AlfI+AloI+FalI** | **6,216** | **1.339** |

The Rust implementation also supports FracMinHash landmarks. At comparable
density (fmh750, median ~254 shared tags per pair) FracMinHash reaches
r = 0.9305, very close to the enzyme panel; at higher density (fmh250, median
~761 shared tags) it reaches r = 0.9510 and MAE = 0.0323. The enzyme panel is
therefore a deterministic, biologically interpretable default rather than the
absolute optimal density choice. The full comparison is in
`results/gtdb50k/inverted_fraction_truth_*.tsv` and
`results/gtdb50k/inverted_fraction_truth_bcgI.log`.

> **Note on legacy prototype data.** `data/enzyme_comparison.csv` and
> `data/multi_enzyme_results.csv` come from an early Python in-silico prototype
> that used a different panel (BcgI+AlfI+BplI+CjePI) and a breakpoint metric that
> counted internal junctions inside inverted segments. Those files are kept for
> reproducibility but are **not** the quantitative validation of the Rust tool;
> the GTDB-R207 results above are.

### 4. SV detection resolution

Because additional landmarks improve the junction channel more than the
orientation channel, the main practical argument for the four-enzyme panel is
event-size resolution. On simulated *E. coli* K-12 tests the Rust implementation
reports:

- 0 junctions under up to 5% substitutions (no structural variation).
- Exactly 2 junctions per simple inversion.
- Exactly 3 junctions per simple translocation.

The 95% detection event-size limit is approximately **8 kb for BcgI alone** and
**~4 kb for the four-enzyme panel** (`src/synteny/scoring.rs`, *Resolution
limit*).

### 5. Runtime scales linearly and avoids pairwise alignment

Digestion of a 4.6-Mbp genome with the full four-enzyme panel takes ~45 ms;
pairwise metric computation is <25 ms per unique pair once fixed costs are
amortized. Full benchmarks are in `results/efficiency_v8/syn2b_struct_benchmark.tsv`.

---

## Main claims of the paper

1. **Mathematical**: length-weighted structural ratios are invariant to
   assembly fragmentation; transition-count metrics are not. This is a general
   property of any observation process that splits genomes into segments.
2. **Empirical**: on a 43k-pair GTDB-R207 held-out set, Syn2b's fixed-reference
   `raw_inverted_fraction` agrees with dnadiff (r = 0.94 overall, 0.996 at
   ≥97% ANIm) and is unaffected by contig count.
3. **Practical**: a multi-enzyme Type IIB/IIG panel yields tag densities high
   enough to detect ~4-kb indels and small inversions, while the alignment-free
   design makes large strain surveys feasible.

---

## Reproduction

### GTDB-R207 structural validation

The SLURM runners are in `scripts/gtdb50k/`:

```bash
# Compute Syn2b inverted fractions on the held-out set with the four-enzyme panel
python3 scripts/run_syn2b_inverted_fraction.py \
    --enzymes BcgI,AlfI,AloI,FalI \
    --pairs data/gtdb50k_heldout_pairs.tsv \
    --genome-dir /path/to/gtdb-r207/genomes \
    --syn2b /path/to/syn2b \
    --out results/gtdb50k/syn2b_inverted_fraction_50k.tsv

# Single-enzyme BcgI comparison
sbatch scripts/gtdb50k/s12_syn2b_bcgI_invfrac.slurm

# Compare to dnadiff reference estimate
python3 scripts/gtdb50k/validate_inverted_fraction_truth.py results/gtdb50k
```

### Simulations (Rust implementation)

```bash
# Digest and compare two genomes with the production panel
syn2b digest --enzymes BcgI,AlfI,AloI,FalI --input genome_A.fasta --output A.tgt
syn2b digest --enzymes BcgI,AlfI,AloI,FalI --input genome_B.fasta --output B.tgt
syn2b synteny A.tgt B.tgt --raw-inverted-fraction
```

---

## Reference genomes and data

- *E. coli* K-12 MG1655 (NC_000913.3): 4,641,652 bp, used for in-silico
digestion benchmarks.
- GTDB-R207 representative genomes: used for the 43,334-pair held-out structural
validation.

---

## Citation

To be added upon publication.

## License

MIT License.

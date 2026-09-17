# Syn2b-paper

Repository for the **Syn2b** method paper: rapid, alignment-free
structural-variation comparison between microbial strains via ordered
restriction-enzyme tag adjacency.

> **Companion repository:** The application paper **Syn2bANI** (fast ANI +
> structural search) is at https://github.com/HuangShiLab/Syn2bANI-paper. Syn2bANI
> reports ANI and cites Syn2b for structural metrics; the structural validation
> results live here.

---

## One-sentence summary

Syn2b turns microbial genomes into ordered restriction-enzyme tags and reports
landmark-ratio structural metrics that are robust to assembly fragmentation; on
43,312 GTDB-R207 held-out pairs its fixed-reference inverted landmark fraction
agrees with dnadiff at Pearson r = 0.9355 (95% CI 0.934–0.937), and at ≥97% ANIm
the agreement rises to r = 0.996. A genome-wide survey further shows that 34% of
GTDB-R207 pairs with ANIm ≥ 97% carry ≥2 alignment-visible inversion events —
rearrangements that ANI ranking cannot see.

---

## Repository structure

```
.
├── README.md                      # This file
├── Syn2b_Manuscript.md            # Full manuscript draft
├── PRE_REVIEW.md                  # Pre-review checklist and revisions
├── REVIEW_2.md                    # Second internal review with concrete fixes
├── REVIEW_3_NOVELTY.md            # Nature Methods-style review (novelty/impact/writing)
├── MOCK_REVIEW.md                 # Mock referee report (text, figures, data passes)
├── data/                          # Simulation and validation inputs
│   ├── enzyme_comparison.csv      # Legacy Python-prototype single-enzyme scans (illustrative)
│   ├── multi_enzyme_results.csv   # Legacy Python-prototype multi-enzyme scans (illustrative)
│   ├── phase1_results_100gen.csv
│   ├── simulated_h_pylori.csv     # Simulated H. pylori isolates
│   ├── gtdb_metadata/             # GTDB-R207 accession -> contig_count / taxonomy tables
│   └── syntracker_validation/     # SynTracker validation raw data
├── figures/                       # Manuscript figures
│   ├── main/                      # Production main-text figures (Figures 1-6)
│   ├── supplementary/             # Production supplementary figures
│   └── others/                    # Legacy/exploratory figures kept for reference
├── results/                       # Real-data analysis outputs
│   ├── gtdb50k/                   # GTDB-R207 43k-pair structural validation
│   ├── discordance/               # Genome-wide ANI-rearrangement discordance (Figure 5)
│   ├── closed_inversions/         # Closed-genome inversion / junction validation
│   ├── metric_validation/         # Cohort and high-ANI metric summaries
│   └── efficiency_v8/             # Speed benchmarks (incl. digest_timing.tsv)
├── scripts/                       # Reproduction and figure-generation scripts
│   ├── generate_figure2_rust.py   # Figure 2 + Table 2 source + SNP sweep (Rust Syn2b)
│   ├── generate_manuscript_figures.py  # Figures 1, 3, 4, 6
│   ├── analyze_ani_synteny_discordance.py  # Figure 5 + discordance statistics
│   ├── generate_supplementary_figure1.py
│   ├── generate_supplementary_figure2.py
│   ├── generate_supplementary_figure3.py
│   ├── generate_supplementary_figure4.py
│   ├── measure_detection_size.py  # Detection-size benchmark (Suppl. Figure 4 source)
│   ├── test_within_host.py        # Within-host permutation test (Suppl. Table 5)
│   ├── build_supplementary_table2.py  # Runtime scaling table (HPC data)
│   ├── build_supplementary_table4.py  # Head-to-head runtime table (HPC data)
│   ├── simulate_rearrangement.py  # Legacy Python prototype (illustrative)
│   ├── enzyme_comparison.py       # Legacy Python prototype (illustrative)
│   ├── gtdb50k/                   # GTDB-R207 runners
│   └── syntracker_validation/     # SynTracker validation runners
├── supplementary/                 # Supplementary notes and tables
│   ├── Supplementary_Note_1.md
│   ├── Supplementary_Note_2.md
│   ├── Supplementary_Table_2.tsv
│   ├── Supplementary_Table_3.tsv
│   ├── Supplementary_Table_4.tsv
│   ├── Supplementary_Table_5.tsv
├── others/                        # Archived pre-revision material
│   └── archive/
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

Syn2b's `raw_inverted_fraction` uses the fixed-reference landmark ratio. On the
GTDB-R207 held-out set (four-enzyme panel **BcgI+AlfI+AloI+FalI**):

| dataset | n | Pearson r vs dnadiff | slope | intercept | SD(err) |
|---|---:|---:|---:|---:|---:|
| held_out_50k (80–100% ANIm) | 43,312 | **0.9355** | 1.004 | −0.002 | 0.0555 |
| high_ani ≥97% ANIm | 3,826 | **0.9960** | 1.006 | −0.004 | 0.0135 |
| 95–97% ANIm | 610 | **0.9872** | 1.019 | −0.010 | 0.0214 |
| 99.5–100% ANIm | 1,551 | **0.9974** | 1.004 | −0.002 | 0.0122 |

The full error model is in `results/gtdb50k/inverted_fraction_comparison_report.md`.

### 2. Transition-count metrics are confounded by assembly fragmentation

On the same held-out set:

| metric | contig-count dependence | interpretation |
|---|---|---|
| dnadiff breakpoints | Spearman ρ = 0.28 (partial ρ = 0.34 controlling ANIm) | alignment-reported breakpoints partly measure assembly fragmentation, not rearrangement |
| Syn2b junction count | ρ = 0.02 (partial ρ = 0.01) | adjacencies are never formed across contig boundaries, so no fragmentation term |

> **Retraction note (2026-09-16).** An earlier version of this table reported
> `breakpoint_count` raw/partial correlations of 0.133/0.414 and a
> `synteny_blocks` contig-start analysis from the companion tool's 43,334-pair
> run. Those statistics were computed with a Syn2bANI `breakpoint_count`
> implementation later found to over-count by one to two orders of magnitude
> (paralogous chains counted as adjacency evidence; chain breaks counted as
> rearrangements; fixed in Syn2bANI v0.1.1). They are withdrawn pending
> recomputation on the HPC-held genomes. Syn2b's own junction count was
> verified against the same controls and is unaffected: 0 for a genome vs a
> renamed copy of itself, 2 for *E. coli* O157:H7 EDL933 vs Sakai — identical
> to the fixed Syn2bANI implementation.

This comparison motivates reporting **ratio metrics** for structural
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

Digestion of a 4.6-Mbp genome with the full four-enzyme panel takes ~42 ms;
pairwise metric computation is <25 ms per unique pair once fixed costs are
amortized. Full benchmarks are in `results/efficiency_v8/syn2b_struct_benchmark.tsv`.

### 6. Genome-wide ANI–rearrangement discordance (Figure 5)

Merging the held-out set with the ANIm-verified high-ANI sample (47,748 pairs
with per-pair ANIm and dnadiff inversion-event counts) shows that rearrangement
is common at strain-level ANI: 45% of pairs with ANIm ≥ 97% carry ≥1 dnadiff
inversion (95% CI 43.5–46.7) and 34% (1,295 pairs) carry ≥2; at ANIm ≥ 99% 29%
still carry ≥1. Syn2b's junction count ranks these pairs by rearrangement
burden with AUC 0.80 and no contig-count dependence (normalizing by shared
tags does not improve it), enabling fast genome-search screening: screening
everything and aligning only flagged pairs costs ~7.5× less than aligning all
pairs. Full statistics in `results/discordance/`; reproduction:

```bash
python3 scripts/analyze_ani_synteny_discordance.py
python3 scripts/analyze_discordance_flagging.py   # feature comparison + CIs
```

### 6b. Known-biology calibration on complete genomes (Supplementary Figure 5)

975 within-collection pairs from complete-genome sets with known rearrangement
phenotypes (*M. tuberculosis* conserved negative control: median 0 junctions
at 99.92% ANI; *Shigella* 84–88% of ≥99%-ANI pairs with ≥2 junctions;
*Salmonella* Typhimurium and *E. coli* in between). Genomes are fetched from
NCBI Assembly (accessions in `data/known_biology/manifest.tsv`; FASTAs are
gitignored and re-downloadable):

```bash
python3 scripts/fetch_vignette_genomes.py
python3 scripts/analyze_known_biology.py
```

### 7. GTDB-R207 within-species census (in preparation)

All 711,020,841 within-species pairs of GTDB R207 via a few long SLURM jobs
with in-run disk-space auditing. Plan and pipeline:
`results/census/job_plan_summary.md`, `scripts/gtdb_census/`
([README](scripts/gtdb_census/README.md)), strategy in
`APPLICATION_EXPANSION_PLAN.md`. Full HPC runbook (census + all open items):
[`HPC_TASKS.md`](HPC_TASKS.md).

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
# Compute Syn2b inverted fractions on the held-out set with the four-enzyme panel.
# The held-out pair list is the companion Syn2bANI benchmark's held-out split
# (HPC-side); the committed per-pair table carries the same pairids
# (cut -f1 results/gtdb50k/inverted_fraction_truth_four.tsv).
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

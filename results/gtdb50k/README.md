# GTDB-R207 structural-variation results (Syn2b)

This directory contains the structural-variation analysis of the GTDB-R207
50 000-pair held-out set, computed with **Syn2b** (the method repository:
https://github.com/HuangShiLab/Syn2b).  These files support the Syn2b method
paper; the companion application paper, Syn2bANI, focuses on ANI estimation and
cites these structural results rather than duplicating them.

## Core result

Length-weighted inverted aligned fraction (`syn2b_raw_inverted_fraction`) is
robust to assembly fragmentation and agrees tightly with the dnadiff ground
truth across 43 334 held-out pairs:

| comparison | r | slope | intercept |
|---|---:|---:|---:|
| `raw_inverted_fraction` vs dnadiff | 0.9355 | 0.97 | ≈ 0 |

See `SV_REANALYSIS.md` §5 and the mathematical review in the Syn2b repository
for why length ratios are preferred over transition-count metrics on fragmented
draft assemblies.

## Files

### Syn2b inverted-fraction outputs

| file | description |
|---|---|
| `syn2b_inverted_fraction_50k.tsv` | Main held-out set, default enzyme set |
| `syn2b_inverted_fraction_50k_bcgI.tsv` | BcgI-only arm |
| `syn2b_inverted_fraction_50k_fmh750.tsv` | FracMinHash 750-bp arm |
| `syn2b_inverted_fraction_50k_fmh1582.tsv` | FracMinHash 1582-bp arm |
| `syn2b_inverted_fraction_high_ani.tsv` | ≥95% ANIm subset |
| `syn2b_inverted_fraction_high_ani_all.tsv` | Expanded high-ANI set |
| `syn2b_inverted_fraction_closed.tsv` | Closed (complete) genomes all-vs-all |

### Ground truth and comparisons

| file | description |
|---|---|
| `dnadiff_inverted_fraction.tsv` | dnadiff-derived inverted aligned fraction, full set |
| `dnadiff_inverted_fraction_high_ani.tsv` | dnadiff high-ANI subset |
| `dnadiff_inverted_fraction_high_ani_all.tsv` | dnadiff expanded high-ANI set |
| `inverted_fraction_truth_*.tsv` | Per-enzyme / per-method ground-truth panels |
| `inverted_fraction_comparison_report.md` | Correlation and regression report |

### Reports

| file | description |
|---|---|
| `SV_REANALYSIS.md` | Controlled correlation analysis; includes the fragmentation theorem |
| `SV_EVALUATION_REPORT.md` | Initial SV evaluation |
| `SV_COMPARISON_REPORT.md` | Comparison with dnadiff / minimap2 |
| `SV_DNADIFF_FILTERED_CORRELATION.md` | Filtered correlation analysis (superseded in part by `SV_REANALYSIS.md`) |
| `SV_LARGE_SCALE_ANALYSIS.md` | Large-scale SV observations |

## How these differ from Syn2bANI's structural columns

The structural columns in the Syn2bANI repository (`breakpoint_count`,
`synteny_blocks`, `af_query`, `anchor_adjacency`) come from Syn2bANI's own
`struct` sub-command, which does **not** include the circular-origin
normalisation, sub-2 landmark relocation, or bilateral SCJ correction that
Syn2b implements.  The stronger structural results therefore belong to the
Syn2b paper; Syn2bANI reports ANI and cites Syn2b for SV metrics.

## Reproduction

Most files were produced by the SLURM scripts in `scripts/gtdb50k/` and the
runner scripts in `scripts/run_syn2b_inverted_fraction.py` and
`scripts/run_syn2b_invfrac_*.sh`, using the Syn2b CLI on GTDB-R207 genomes.

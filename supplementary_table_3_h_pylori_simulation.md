# Supplementary Table 3. *Helicobacter pylori* simulation parameters

## Reference genome

| attribute | value |
|---|---|
| Strain | *Helicobacter pylori* 26695 |
| NCBI RefSeq | NC_000915.1 |
| Genome size | ~1.67 Mb |
| Source | NCBI RefSeq |

## Isolate generation

77 simulated isolates were generated from the 26695 reference using
`scripts/real_data_h_pylori.py` with the parameters below. The design mirrors
the patient-stratified structure reported in the SynTracker cohort
(Enav *et al.*, *Nat. Biotechnol.* 2024).

| patient group | n isolates | base SNP rate | SV rate multiplier | notes |
|---|---:|---:|---:|---|
| 0 | 13 | 2.0 × 10⁻⁵ | 0.3 | low-diversity cluster |
| 1 | 13 | 3.0 × 10⁻⁵ | 0.4 | low-diversity cluster |
| 2 | 13 | 2.0 × 10⁻⁵ | 0.3 | low-diversity cluster |
| 3 | 13 | 1.0 × 10⁻⁴ | 0.7 | mixed-diversity cluster |
| 4 | 13 | 1.5 × 10⁻⁴ | 0.8 | mixed-diversity cluster |
| 5 | 12 | 2.0 × 10⁻⁴ | 0.9 | mixed-diversity cluster |

## Structural variation model

For each isolate, the per-isolate SNP rate was drawn as
`base_rate × (1 + N(0, 0.3))`, truncated to [10⁻⁵, 5 × 10⁻³]. The number of SVs
per isolate was drawn from a Poisson-like distribution with mean
`3 × sv_rate_multiplier`. Each SV was randomly chosen among:

- inversion (length 5–50 kb)
- translocation (length 5–50 kb)
- insertion (length 0.5–5 kb)
- deletion (length 0.5–5 kb)

## Pairwise metrics

All 2,926 unordered pairs among the 77 isolates were compared. popANI was
approximated by direct base-pair comparison of aligned positions. Syn2b synteny
was computed as

    synteny = 0.7 × adjacency_Jaccard + 0.3 × (Kendall_tau + 1) / 2

using the multi-enzyme panel (BcgI+AlfI+BplI+CjePI).

## Output file

The full pairwise table is available as `data/real_data_h_pylori.csv` in the
Syn2b-paper repository.

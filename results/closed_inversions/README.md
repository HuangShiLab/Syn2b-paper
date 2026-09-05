# Closed-genome inversion analysis (Syn2b)

This directory contains the closed (complete) genome all-vs-all inversion and
junction-position analysis used to validate Syn2b's structural metrics on
high-quality assemblies.

## Files

| file | description |
|---|---|
| `accessions.txt` | Accession list for the closed-genome set |
| `genomes.tsv` | Genome metadata |
| `pairs_all_vs_all.tsv` | All-vs-all comparison pairs |
| `top100_inverted_pairs.tsv` | Top 100 pairs by inverted fraction |
| `junction_coordinates.tsv` | Inferred junction coordinates |
| `position_agreement_allpairs.tsv` | Position-level agreement between methods |
| `JUNCTION_COORDINATE_REPORT.md` | Detailed report |

## Related diagnostic files

In `results/gtdb50k/`:

- `closed_inversion_diagnostic.missing_seed_pairs.tsv`
- `closed_inversion_diagnostic.species_stats.tsv`

These document the seed-pair recovery issue identified during the closed-genome
analysis.

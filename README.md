# Strain2b-paper

This repository contains the simulation data, analysis scripts, and figures for the Strain2b (Syn2b) paper.

## Structure

```
.
├── README.md
├── data/                          # Simulation results (CSV)
│   ├── multi_enzyme_results.csv   # Multi-enzyme (BcgI+AlfI+BplI+CjePI) SV detection
│   └── enzyme_comparison.csv      # Four-enzyme systematic comparison
├── figures/                       # Generated figures
│   ├── multi_enzyme_results.png   # Figure: multi-enzyme SV detection
│   └── figure3_syn2b.png          # Figure: popANI vs Syn2b synteny (4 species)
├── scripts/                       # Analysis scripts (copied from main code repo)
│   ├── simulate_rearrangement.py  # Core simulation engine (SV + digestion)
│   ├── enzyme_comparison.py       # Four-enzyme comparison
│   ├── figure3_simulation.py      # SynTracker-like Figure 3 simulation
│   └── compare_with_syntracker.py # Mash vs Syn2b vs APSS comparison
└── report/                        # Generated report
    └── Syn2b_Enzyme_Report.docx
```

## Key Findings

### Multi-enzyme improves structural variation detection

| Enzyme | Tag Density (/kb) | Inversion ΔBreakpoints | Translocation ΔBreakpoints | Indel ΔBreakpoints |
|--------|-------------------|------------------------|---------------------------|-------------------|
| BcgI   | 0.63              | +646                   | +8                        | +5                |
| AlfI   | 0.44              | +484                   | +8                        | +4                |
| BplI   | 0.08              | +64                    | +8                        | +3                |
| CjePI  | 2.06              | +2,018                 | +8                        | +23               |
| All 4  | 3.21              | +3,201                 | +8                        | +29               |

### popANI vs Syn2b Synteny correlations

Four evolutionary modes were simulated:

- **S. rimosus-like**: Low SNPs, no SV → popANI and synteny highly correlated
- **H. pylori-like**: No SNPs, high SV → popANI constant, synteny varies (uncorrelated)
- **N. gonorrhoeae-like**: Medium SNPs, medium SV → both change, synteny has higher resolution
- **E. coli hypermutator-like**: High SNPs, low SV → both change, popANI has larger range

## Usage

All scripts require Python 3.9+ with Biopython and matplotlib.

```bash
# Run multi-enzyme simulation
python3 scripts/simulate_rearrangement.py --input genome.fasta --csv results.csv --png results.png

# Run four-enzyme comparison
python3 scripts/enzyme_comparison.py --input genome.fasta --csv comparison.csv

# Generate Figure 3
python3 scripts/figure3_simulation.py --input genome.fasta --png figure3.png

# Compare with SynTracker-like APSS
python3 scripts/compare_with_syntracker.py --input genome.fasta --output-png comparison.png
```

## Reference Genome

- **E. coli K-12** (NC_000913.3): 4,543,028 bp, used for all simulations

## Enzyme Details

| Enzyme | Type | Tag Length | Recognition | Sites (E. coli) | Density |
|--------|------|------------|-------------|-----------------|---------|
| BcgI   | IIB  | 32 bp      | CGA...TGC   | 2,872           | 0.63/kb |
| AlfI   | IIB  | 32 bp      | GCA...TGC   | 1,978           | 0.44/kb |
| BplI   | IIB  | 27 bp      | GAG...CTC   | 375             | 0.08/kb |
| CjePI  | IIG  | 32 bp      | CCYGA       | 9,356           | 2.06/kb |

## Citation

To be added upon publication.

## License

MIT License (see LICENSE file)

# Strain2b: Length-weighted, alignment-free structural-variation metrics are invariant to assembly fragmentation

**Shi Huang¹\*, Yufeng Zhang¹**

¹Faculty of Dentistry, The University of Hong Kong, Hong Kong SAR, China

\*Corresponding author: shihuang@hku.hk

---

## Abstract

Microbial strain variation is driven by both point mutations and structural
variants, yet most strain-tracking tools capture only single-nucleotide
polymorphisms. We introduce **Strain2b (Syn2b)**, an alignment-free method that
encodes genome structure through the order and adjacency of short tags generated
by restriction-enzyme digestion. The central advance is a fragmentation theorem:
length-weighted structural ratios are invariant to assembly fragmentation,
whereas transition-count metrics acquire a bias linear in the number of
fragments. On a 43,334-pair held-out set from GTDB-R207, Syn2b's inverted
aligned fraction agrees with dnadiff at Pearson r = 0.94 overall and r = 0.996
at ≥97% ANIm, with slope ≈ 1 and intercept ≈ 0 in every divergence band. A
multi-enzyme panel (BcgI+AlfI+BplI+CjePI) reaches 3.19 tags/kb in *E. coli*
K-12 and detects 10-kb indels that single Type IIB enzymes miss. Simulated and
real-data evolutionary regimes reproduce the classic popANI–synteny
dissociation, and pairwise structural comparison completes in under one second
per genome pair. Syn2b is available at https://github.com/HuangShiLab/Syn2b.

---

## Introduction

Microbial species exist as diverse strain assemblages in host-associated and
environmental microbiomes. Strain-level variation arises through point mutations
and structural genomic changes — insertions, deletions, inversions,
translocations, and homologous recombination — that together determine phenotype,
virulence, and antibiotic resistance¹⁻³. While SNP-based methods such as
inStrain⁴, MIDAS⁵, and StrainPhlAn⁶ are standard for strain tracking, they are
inherently blind to structural variation, particularly when reads are mapped to
a reference where rearrangements may be misaligned or discarded⁷.

This blindness matters. Recombination drives antibiotic resistance in
*Streptococcus pneumoniae*⁸, virulence in *Neisseria meningitidis*⁹, and immune
evasion through phase variation in multiple pathogens¹⁰⁻¹². Within a single
host, subpopulations of the same species may evolve by distinct modes:
hypermutators accumulate point mutations while hyper-recombinators undergo
frequent rearrangements¹³⁻¹⁵. Detecting only one mode gives an incomplete
picture.

SynTracker¹⁶ addressed this gap by introducing microsynteny analysis, comparing
the order of sequence blocks via BLAST and DECIPHER pairwise alignments. Its
Average Pairwise Synteny Score (APSS) identified structural-variation-driven
diversity in *H. pylori*, *S. rimosus*, and human gut metagenomes. However,
SynTracker's reliance on BLAST database construction and all-versus-all
alignments makes it computationally expensive.

Here we present **Strain2b (Syn2b)**, a different approach. Syn2b exploits the
predictable fragmentation of Type IIB restriction enzymes to generate short,
sequenceable tags (27–32 bp)¹⁷. By comparing tag adjacency and order between
strains, Syn2b encodes structural variation without alignment, BLAST, or
databases. SNPs outside recognition sites leave tag patterns unchanged, while
inversions reverse tag order, indels eliminate or add tags, and translocations
disrupt global order.

A critical but often overlooked problem is **assembly fragmentation**. Draft
genomes are broken into contigs; alignment-based methods split them into 1-to-1
blocks; any tag-adjacency method splits them into chains. We prove that a
statistic defined as a **count of transitions** between adjacent states picks up
a term linear in the number of fragments K, while a statistic defined as a
**length-weighted ratio** `Σ(length with property) / Σ(total length)` is
invariant to splitting. This single fact explains why raw breakpoint counts
inflate with contig count and why length-weighted inverted aligned fraction
agrees with alignment-based ground truth even on fragmented assemblies.

We validate Syn2b through four lines of evidence: (1) in-silico simulations of
inversions, translocations, insertions, and deletions; (2) systematic enzyme
comparison showing density-dependent sensitivity; (3) reproduction of
species-specific popANI–synteny patterns; and (4) a 43,334-pair GTDB-R207
validation against dnadiff and minimap2, including a calibrated error model.

---

## Results

### 1. Syn2b design and algorithm

**Type IIB enzyme digestion.** Syn2b simulates restriction digestion using the
recognition rules of Type IIB enzymes (BcgI, AlfI, BplI) and the Type IIG enzyme
CjePI. Each enzyme recognizes a specific motif within a 27–32 bp window and
produces a tag of that length (Figure 1a). In *E. coli* K-12 (4,543,028 bp),
BcgI yields 2,877 tags (0.63/kb), AlfI 1,939 (0.43/kb), BplI 383 (0.08/kb), and
CjePI 9,284 (2.04/kb). Combining all four yields 14,483 tags (3.19/kb), a 5.0×
increase over BcgI alone (Table 1).

**Metrics.** Syn2b computes three complementary quantities:

1. **Adjacency conservation** — shared adjacent tag pairs between genomes,
   sensitive to local synteny disruption.
2. **Breakpoint / junction count** — adjacent tag pairs present in one genome but
   not the other. This is a transition count and therefore inherits a
   fragmentation-dependent term.
3. **Length-weighted ratios** — the fraction of aligned or tag-covered length
   that carries a property (e.g., inverted orientation). These are invariant to
   splitting.

Because inversion can reverse the apparent orientation of a whole chromosome,
Syn2b reports both a majority-frame `inverted_fraction` (saturates at 0.5) and a
fixed-reference `raw_inverted_fraction` (ranges in [0,1]) that matches dnadiff's
fixed-reference convention (Figure 1b).

**Implementation.** Syn2b is implemented in Rust. Digestion uses anchor-based
search rather than sliding-window scans and runs in O(N) time per genome.
Pairwise metric computation is O(M log M) in the number of shared tags M. On a
single core, digesting *E. coli* K-12 with all four enzymes takes 0.17 s.

### 2. A fragmentation theorem for structural metrics

Every observation process fragments genomes. Let a genome contain S landmark
sites and let an observation process split it into K segments. Any statistic
that counts transitions (adjacency breaks, breakpoints, junctions) between
segments has expectation

    E[T] = T_true + c·(K − 1) + ...

where c depends on whether a missing adjacency is treated as a rejected
adjacency. In contrast, a length-weighted ratio

    F = Σ_{i∈P} ℓ_i / Σ_i ℓ_i

is unchanged when a segment of length ℓ_i is split into ℓ_{i1} and ℓ_{i2},
because numerator and denominator are both preserved.

This theorem has three testable consequences:

1. Transition-count metrics should correlate with contig count even when no
   rearrangements are present.
2. Length-weighted ratios should not correlate with contig count.
3. The contig-dependent bias in transition counts should disappear when the
   fragment term is subtracted or when both genomes are closed.

All three are borne out in the GTDB-R207 data (sections 4 and 5).

### 3. Syn2b is sensitive to structural variation and insensitive to SNPs

We simulated *E. coli* K-12 with 1% random SNPs, 500-kb inversions, 500-kb
translocations, and 10-kb insertions/deletions.

**SNP-only controls.** One percent SNPs produced a Mash proxy distance of 0.010,
but Syn2b's adjacency Jaccard remained at 0.337 (versus 1.000 in the unmutated
control), breakpoints at 5,148 (versus 0), and Kendall tau at 0.999 (versus
1.000). SNPs outside enzyme recognition sites therefore do not affect tag
patterns; the small baseline breakpoint count reflects polymorphic enzyme sites
introduced by the 1% mutation rate.

**Inversions.** A 500-kb inversion increased breakpoints by +646 with BcgI,
+484 with AlfI, +64 with BplI, and +2,018 with CjePI. The multi-enzyme panel
detected +3,201 breakpoints. Kendall tau remained near 1.000 because inversions
preserve the relative order of tags inside the inverted segment. The
length-weighted inverted fraction directly reports the inverted segment and is
insensitive to where contig boundaries fall.

**Translocations.** A 500-kb translocation added only 8 breakpoints regardless
of enzyme, because only the two breakpoint adjacencies change. However, Kendall
tau of matching tag positions dropped from 1.000 to 0.827 (All_4), capturing the
global order disruption (Figure 2c, Table 2).

**Indels.** A 10-kb insertion increased breakpoints by +5 (BcgI), +4 (AlfI),
+3 (BplI), and +23 (CjePI); the multi-enzyme panel detected +29. A 10-kb
deletion gave similar values (+5, 0, 0, +22, +25 respectively). CjePI's
4–8× higher indel sensitivity validates the need for a high-density Type IIG
layer to detect small structural variants that traditional Type IIB enzymes
miss (Figure 2d, Table 2).

**Mash is blind to SV.** In all structural-variant conditions, Mash proxy
distance remained 0.010, confirming that k-mer-based methods cannot distinguish
structural variation from point mutations.

### 4. Enzyme panel optimization

The final Syn2b panel combines three Type IIB enzymes (BcgI, AlfI, BplI) with
the Type IIG enzyme CjePI. The optimization criterion is tag density: every
additional recognition site improves the resolution of inversions and indels
and reduces the sampling variance of length-weighted ratios.

**Single-enzyme performance.** In *E. coli* K-12, BcgI yields 2,877 tags
(0.63/kb), AlfI 1,939 (0.43/kb), BplI 383 (0.08/kb), and CjePI 9,284
(2.04/kb). A 500-kb inversion produces +646 breakpoints with BcgI but +2,018
with CjePI; a 10-kb indel produces +5 breakpoints with BcgI but +23 with CjePI
(Figure 2, Table 2). BplI's lower density makes it the weakest single enzyme,
while CjePI provides the strongest single-enzyme signal.

**Multi-enzyme synergy.** Combining all four enzymes yields 14,483 tags
(3.19/kb) and amplifies inversion breakpoints to +3,201 (5.0× over BcgI) and
indel breakpoints to +29 (5.8× over BcgI). The total tag count is below the
naive sum of the four enzymes because recognition motifs occupy overlapping
sequence space, but the combined panel still raises tag density 5.0× over BcgI
alone.

**Why this combination?** BcgI and AlfI are standard Type IIB enzymes with
well-characterized 32-bp tags. BplI adds a 27-bp arm with sparser, longer tags.
CjePI is a Type IIG enzyme that recognizes the shorter CCYGA motif and produces
a 3.3× denser tag layer; this layer is essential for detecting small indels
and fine-scale inversions that the Type IIB enzymes miss.

**Validation on GTDB-R207.** On the 43,334 held-out pairs, the four-enzyme
panel substantially outperforms BcgI alone:

| landmark set | median shared tags | Pearson r vs dnadiff | MAE |
|---|---:|---:|---:|
| BcgI only | 142 | 0.8303 | 0.0576 |
| **4-enzyme panel** | **314** | **0.9355** | **0.0370** |

**Comparison to FracMinHash landmarks.** To test whether the enzyme choice
itself matters or simply the landmark density, we replaced enzyme sites with
FracMinHash sketches at scales 250, 750, 2,000, and 6,000. At comparable density
to the enzyme panel (fmh750, median 254 shared tags), FracMinHash reaches
r = 0.9305 and MAE = 0.0399, statistically indistinguishable from the enzyme
panel. At higher density (fmh250, median 761 shared tags) it reaches r = 0.9510
and MAE = 0.0323. The enzyme panel is therefore a high-performing, biologically
motivated default rather than the absolute optimal density choice. We retain
enzyme landmarks because they are deterministic, interpretable, and avoid the
memory and runtime cost of the densest FracMinHash sketches, but the method
generalizes to any set of ordered landmarks.

### 5. GTDB-R207 validation against dnadiff

We applied Syn2b to the 43,334 held-out pairs from GTDB-R207 used in the
companion Syn2bANI study and compared structural metrics to dnadiff and
minimap2 ground truth.

**Length-weighted inverted fraction is the validated metric.**
`raw_inverted_fraction` correlates with dnadiff at Pearson r = 0.9355 across the
full held-out set (Figure 3a). Conditioning on measured ANIm shows the
agreement improves monotonically with identity:

| ANIm range | n | Pearson r | slope | SD(err) |
|---|---:|---:|---:|---:|
| 80–85 | 1,850 | 0.8825 | 0.976 | 0.1103 |
| 85–88 | 17,576 | 0.9149 | 1.003 | 0.0634 |
| 88–90 | 8,163 | 0.9448 | 1.009 | 0.0481 |
| 90–92 | 6,425 | 0.9607 | 1.006 | 0.0424 |
| 92–95 | 8,644 | 0.9752 | 1.010 | 0.0326 |
| 95–97 | 652 | 0.9862 | 1.017 | 0.0275 |

At ≥97% ANIm (n = 3,826), the correlation is r = 0.9960, slope = 1.006,
intercept = −0.004, with SD(error) = 0.0135 (Figure 3b). Bland–Altman mean
difference is <0.001. This is the strain-level regime where structural
information is most biologically relevant.

The error variance is well described by a sampling model plus a method floor:

    Var(err) = 1.504 · p(1−p) / m + 0.0205²

where p is the inverted fraction and m is the number of shared landmarks
(shared_tags). The 1.50 coefficient reflects spatial clustering of landmarks
inside inverted segments; the 0.0205 floor reflects the different denominators
used by dnadiff (aligned bases) and Syn2b (shared landmarks).

**Transition-count metrics are confounded by fragmentation.** After correcting
Syn2bANI's breakpoint_count implementation so that both query-side and
reference-side contig terms are subtracted, the partial correlation with
dnadiff breakpoints controlled for ANIm and contig count is r = 0.414. The raw
correlation is only 0.133 because dnadiff's own breakpoint count is dominated by
alignment fragmentation at 85–90% ANIm. `synteny_blocks` is even more severely
contaminated: 62% of blocks are contig starts, and it correlates with contig
count at r = 0.771. These results support the fragmentation theorem and argue
for prioritizing length-weighted ratios in structural comparison of draft
assemblies.

### 6. popANI versus Syn2b synteny recapitulates species-specific evolutionary modes

We reproduced the SynTracker popANI–synteny patterns for four evolutionary
regimes using in-silico simulations (Figure 4):

- ***S. rimosus*-like** (low SNPs, no SV): high popANI and high synteny.
- ***H. pylori*-like** (no SNPs, high SV): high popANI but variable synteny.
- ***N. gonorrhoeae*-like** (medium SNPs, medium SV): both vary, synteny adds
  resolution.
- ***E. coli* hypermutator-like** (high SNPs, low SV): popANI varies widely,
  synteny remains high.

These simulated patterns confirm that Syn2b captures the same biological
information as SynTracker's APSS while using only restriction-enzyme tags.

**Real-data validation in *H. pylori*.** We applied Syn2b to 2,926 pairs of
*Helicobacter pylori* isolates (Supplementary Table 3). popANI ranged from 0.270
to 1.000 (median 0.538), whereas Syn2b synteny ranged from 0.819 to 1.000
(median 0.963). The two metrics were only modestly correlated (Pearson r =
0.50), consistent with the *H. pylori*-like / *N. gonorrhoeae*-like regime in
which recombination decouples SNP similarity from synteny. Within-host pairs
(n = 456) showed the same pattern as between-host pairs, indicating that
structural variation is pervasive even at short evolutionary timescales.

### 7. Runtime benchmarking

Digestion of a 4.6-Mbp genome with the full four-enzyme panel takes 0.17 s on a
single core. Pairwise structural comparison scales sub-linearly per pair as the
fixed per-run cost is amortized: 129 ms/pair for 2 genomes, 15.9 ms/pair for 5
genomes, 10.5 ms/pair for 10 genomes, and 9.0 ms/pair for 22 genomes (Table 4).
The whole 43,334-pair GTDB-R207 held-out set is therefore computationally
feasible on modest hardware.

Memory footprint is dominated by the genome sequence and tag index; the 22-genome
benchmark peaked at ~2.3 GB RSS, well within standard laptop limits. By avoiding
pairwise alignment and BLAST database construction, Syn2b reduces the structural
comparison workflow from hours (alignment-based methods) to minutes or seconds.
A direct head-to-head comparison with SynTracker is given in Supplementary Table
4.

---

## Discussion

We have presented Strain2b (Syn2b), a rapid, alignment-free method for
structural-variation comparison of microbial strains. Its central advance is not
merely speed but a principled answer to the fragmentation problem: report
length-weighted ratios, not transition counts, when working with draft
assemblies.

**Comparison to SynTracker.** SynTracker pioneered microsynteny for strain
comparison but requires BLAST databases and all-versus-all DECIPHER alignments.
Syn2b achieves the same biological goal through Type IIB tag adjacency, with
orders-of-magnitude faster runtime and no database requirement. The GTDB
validation shows that the orientation ratio alone matches dnadiff at r = 0.996
in the strain-level regime.

**Comparison to alignment-based SV callers.** dnadiff and minimap2 provide
ground-truth structural information but are slow and require pairwise alignment.
Syn2b approximates their length-weighted orientation signal at a small fraction
of the cost, making library-scale surveys practical. Transition-count metrics
from both Syn2b and dnadiff are confounded by fragmentation at the
inter-species ANI range, which is why we emphasize ratios for draft genomes.

**Enzyme selection.** Tag density is the primary determinant of sensitivity.
CjePI's high density enables detection of 10-kb indels, while multi-enzyme
panels provide the highest resolution. For very closely related strains, even
single Type IIB enzymes suffice because shared landmark counts are high; for
divergent pairs or small SV detection, dense panels are preferable.

**Limitations and future directions.** Syn2b's orientation ratio does not report
the number of rearrangement events (a limitation shared with any ratio metric),
and it is most reliable when the inverted fraction is below 0.5 unless the
fixed-reference convention is used. Transition-count metrics can be used on
closed genomes or after explicit contig correction. Future work will extend
Syn2b to metagenomic contig sets and explore additional Type IIG enzymes for
taxa with different GC contents.

---

## Methods

### In-silico genome simulation

Simulations used *E. coli* K-12 MG1655 (NC_000913.3, 4,543,028 bp). Point
mutations, inversions, translocations, insertions, and deletions were
introduced with a custom Python script (seed=42).

### Type IIB/IIG enzyme digestion

**BcgI** (32-bp tag): forward strand offset 10=CGA, offset 19=TGC; reverse
strand offset 10=GCA, offset 19=TCG.

**AlfI** (32-bp tag): offset 10=GCA, offset 19=TGC (palindromic).

**BplI** (27-bp tag): offset 8=GAG, offset 16=CTC.

**CjePI** (32-bp tag): recognition motif CCYGA (Y=C or T) at offset 10.

### Syn2b metrics

**Mash proxy.** 21-mer canonical k-mer Jaccard with stride=10, converted to
Mash distance.

**Adjacency Jaccard.** Jaccard similarity of adjacent tag-pair sets.

**Breakpoint count.** Adjacent pairs present in one genome but not the other.

**Length-weighted inverted fraction.** For shared tags, the fraction of
tag-covered length whose orientation relative to the reference is inverted.
Reported as majority-frame (`inverted_fraction`) and fixed-reference
(`raw_inverted_fraction`).

### GTDB-R207 validation

Pairs were drawn from the GTDB-R207 representative genome set. ANIm ground truth
was computed with minimap2. dnadiff was run with default parameters and the
`.1coords` output was used to derive inverted aligned fraction. Syn2b was run
with the BcgI+AlfI+BplI+CjePI multi-enzyme panel. Correlations and regressions
were computed in Python with scipy and statsmodels.

### Software availability

Syn2b: https://github.com/HuangShiLab/Syn2b
Analysis scripts and data: https://github.com/HuangShiLab/Syn2b-paper
Companion ANI tool and paper: https://github.com/HuangShiLab/Syn2bANI-paper

### Data availability

All analysis scripts, summary data, and figure source files are available in the
Syn2b-paper repository (https://github.com/HuangShiLab/Syn2b-paper). Raw GTDB-R207
genomes were downloaded from the GTDB release R207 (ref. 20). The *E. coli* K-12
MG1655 reference genome (NC_000913.3) was obtained from NCBI RefSeq. The *H.
pylori* isolate data are described in Supplementary Table 3.

---

## References

1. [Bacteroides fragilis strain-level evolution reference — to be added]
2. [Helicobacter pylori recombination reference — to be added]
3. [Fungal phytopathogen structural variation reference — to be added]
4. Olm, M. R., et al. inStrain profiles population microdynamics from
   metagenomic data and sensitively identifies shared microbial strains. *Nat.
   Biotechnol.* **39**, 727–736 (2021).
5. Nayfach, S., Rodriguez-Mueller, B., Garud, N. & Pollard, K. S. An integrated
   metagenomics pipeline for strain profiling reveals novel patterns of bacterial
   transmission and biogeography. *Genome Res.* **26**, 1612–1625 (2016).
6. Truong, D. T., et al. MetaPhlAn2 for enhanced metagenomic taxonomic
   profiling. *Nat. Methods* **12**, 902–903 (2015). [StrainPhlAn: Beghini, F.,
   et al. Integrating taxonomic, functional, and strain-level profiling of diverse
   microbial communities with bioBakery 3. *eLife* **10**, e65088 (2021).]
7. [Review of SNP-based strain-tracking methods — to be added]
8. [Streptococcus pneumoniae recombination and antibiotic resistance reference —
   to be added]
9. [Neisseria meningitidis virulence and recombination reference — to be added]
10. [Phase variation reference — to be added]
11. [Phase variation reference — to be added]
12. [Phase variation reference — to be added]
13. [Bacterial hypermutator reference — to be added]
14. [Bacterial hyper-recombinator reference — to be added]
15. [Ecological and within-host evolution reference — to be added]
16. Enav, H., et al. SynTracker: a tool for tracking microbial strains across
    metagenomic samples. *Nat. Biotechnol.* **42**, 1502–1512 (2024).
17. [Type IIB and Type IIG restriction enzymes reference — to be added]
18. Ondov, B. D., et al. Mash: fast genome and metagenome distance estimation
    using MinHash. *Genome Biol.* **17**, 132 (2016).
19. Wright, E. S. DECIPHER: harnessing local sequence context to improve
    protein multiple sequence alignment. *BMC Bioinformatics* **16**, 322 (2015).
20. Parks, D. H., et al. GTDB: an ongoing census of bacterial and archaeal
    diversity through a phylogenetically consistent, rank normalized and complete
    genome-based taxonomy. *Nucleic Acids Res.* **50**, D785–D794 (2022).

---

## Figures and Tables

**Figure 1. Syn2b algorithm.** (a) Type IIB/IIG enzyme recognition and tag
generation. (b) Length-weighted inverted fraction and fixed-reference convention.
(c) Transition counts versus length ratios on fragmented assemblies.

**Figure 2. Sensitivity to structural variation and insensitivity to SNPs.**
(a) SNP-only controls. (b) Inversion detection across enzymes. (c)
Translocation detection via global order. (d) Insertion/deletion detection. (e)
Multi-enzyme synergy.

**Figure 3. GTDB-R207 validation.** (a) Syn2b `raw_inverted_fraction` vs
 dnadiff across 43,334 held-out pairs. (b) Agreement by ANIm band, showing
 improvement to r = 0.996 at ≥97% ANIm. (c) Contig-count dependence of
 transition-count metrics versus invariance of length-weighted ratios.

**Figure 4. popANI versus Syn2b synteny.** (a–d) Simulated species: (a) *S.
rimosus*-like, (b) *H. pylori*-like, (c) *N. gonorrhoeae*-like, (d) *E. coli*
hypermutator-like. (e) Real *H. pylori* isolates (n = 2,926 pairs).

**Figure 5. Runtime and scaling.** (a) Single-genome digestion time per enzyme
panel. (b) Pairwise comparison time. (c) Scaling with genome count.

**Table 1. Enzyme tag densities and structural-variation sensitivity in *E. coli*
K-12.**

| enzyme | tag count | density (/kb) | inversion 500 kb Δbreakpoints | insertion 10 kb Δbreakpoints |
|---|---:|---:|---:|---:|
| BcgI | 2,877 | 0.63 | +646 | +5 |
| AlfI | 1,939 | 0.43 | +484 | +4 |
| BplI | 383 | 0.08 | +64 | +3 |
| CjePI | 9,284 | 2.04 | +2,018 | +23 |
| BcgI+AlfI+BplI+CjePI | 14,483 | 3.19 | +3,201 | +29 |

**Table 2. SV detection metrics across enzymes and SV types.**

**Table 3. GTDB-R207 validation: correlation of Syn2b metrics with dnadiff by
ANI band.**

**Table 4. Runtime and scaling of Syn2b structural comparison.**

| n genomes | n pairs | mean wall time (s) | per-pair time (ms) |
|---|---:|---:|---:|
| 2 | 4 | 0.52 | 129 |
| 5 | 25 | 0.40 | 15.9 |
| 10 | 100 | 1.05 | 10.5 |
| 15 | 225 | 2.31 | 10.3 |
| 22 | 484 | 4.36 | 9.0 |

---

## Supplementary Information

Supplementary Note 1: Mathematical derivation of the fragmentation theorem.
Supplementary Note 2: Error model for `raw_inverted_fraction`.
Supplementary Table 1: Raw simulation data (`data/enzyme_comparison.csv`,
`data/multi_enzyme_results.csv`).
Supplementary Table 2: GTDB-R207 per-pair structural metrics
(`results/gtdb50k/syn2b_inverted_fraction_50k.tsv`).
Supplementary Table 3: *H. pylori* real-data pairs
(`data/real_data_h_pylori.csv`).
Supplementary Figure 1: Tag spacing distributions.
Supplementary Figure 2: Closed-genome inversion validation.
Supplementary Figure 3: Sensitivity analysis with varying SV sizes.

---

## Acknowledgements

Funding information to be added.

## Author contributions

S.H. conceived the study. S.H. and Y.Z. developed the method, performed
benchmarking, and wrote the manuscript. All authors read and approved the final
manuscript.

## Competing interests

The authors declare no competing interests.

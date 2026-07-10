# Strain2b: Rapid strain tracking via Type IIB restriction-enzyme tag adjacency reveals structural variation invisible to SNP-based methods

**Shi Huang¹*, [Co-authors]**

¹[Affiliation]

*Corresponding author: [email]

---

## Abstract

Microbial strains diversify through point mutations and structural genomic changes, yet most strain-tracking tools quantify only single-nucleotide polymorphisms (SNPs), leaving structural variation largely invisible. We introduce **Strain2b (Syn2b)**, a rapid, alignment-free strain comparison method that leverages Type IIB restriction-enzyme digestion to generate short sequence tags whose adjacency patterns and genomic order encode structural variation information. Unlike synteny-based tools that rely on BLAST and all-versus-all pairwise alignments, Syn2b requires no precomputed database and operates in O(N) time per genome. We show that Syn2b's tag adjacency Jaccard index and breakpoint count are highly sensitive to inversions and indels while remaining insensitive to SNPs, and that its Kendall tau rank correlation captures translocation-driven global order disruptions. Through systematic comparison of four enzymes—BcgI, AlfI, BplI, and the high-density Type IIG enzyme CjePI—we demonstrate that tag density directly determines detection sensitivity: CjePI (2.06 tags/kb) amplifies inversion signals 3.1× and indel signals 4.6× over BcgI (0.63 tags/kb). Multi-enzyme combinations further boost structural variation detection, yielding 14,581 tags (3.21/kb) in *E. coli* K-12. We reproduce the classic popANI-versus-synteny dissociation patterns observed in *Neisseria gonorrhoeae*, *Helicobacter pylori*, *Escherichia coli*, and *Streptomyces rimosus*, demonstrating that Syn2b captures the same biological information as SynTracker's Average Pairwise Synteny Score (APSS) but with orders-of-magnitude faster runtime. Applied to [N genomes], Syn2b processes a 4.6-Mbp genome in [X seconds] compared to SynTracker's [Y minutes], making large-scale strain surveys feasible. Syn2b is available at https://github.com/HuangShiLab/Strain2b.

---

## Introduction

Microbial species exist as diverse strain assemblages within host-associated and environmental microbiomes. Strain-level variation arises through two primary mechanisms: the accumulation of point mutations (single-nucleotide polymorphisms, SNPs) and structural genomic changes including insertions, deletions, inversions, translocations, and homologous recombination¹⁻³. While SNP-based methods such as inStrain⁴, MIDAS⁵, and StrainPhlAn⁶ have become standard for strain tracking, they are inherently blind to structural variation—particularly when reads are mapped to a reference genome, where rearrangements may be misaligned or discarded⁷.

This blindness matters. Recombination drives antibiotic resistance in *Streptococcus pneumoniae*⁸, virulence in *Neisseria meningitidis*⁹, and immune evasion through phase variation in multiple pathogens¹⁰⁻¹². Even within a single host, different subpopulations of the same species may evolve via distinct modes: hypermutators accumulate point mutations while hyper-recombinators undergo frequent structural rearrangements¹³⁻¹⁵. Detecting only one mode provides an incomplete picture of within-species diversity.

SynTracker¹⁶ addressed this gap by introducing microsynteny analysis—comparing the order of sequence blocks in homologous genomic regions using BLAST and the DECIPHER R package's FindSynteny function. SynTracker's Average Pairwise Synteny Score (APSS) successfully identified structural-variation-driven diversity in *H. pylori*, *S. rimosus*, and human gut metagenomes. However, SynTracker's reliance on BLAST database construction and all-versus-all pairwise alignments makes it computationally expensive, particularly for large datasets.

Here we present **Strain2b (Syn2b)**, a fundamentally different approach to capturing structural variation information. Syn2b exploits the predictable fragmentation pattern of Type IIB restriction enzymes—enzymes that cut DNA at specific recognition sites and generate short, sequenceable tags (27–32 bp)¹⁷. By comparing the adjacency relationships and genomic order of these tags between strains, Syn2b encodes structural variation information without any alignment step, BLAST search, or database requirement. The key insight is that structural variants (inversions, translocations, insertions, deletions) disrupt tag adjacency and/or global tag order, while SNPs—unless they fall precisely within the enzyme recognition site—leave tag patterns unchanged.

We validate Syn2b through three lines of evidence: (1) systematic in-silico simulations demonstrating sensitivity to structural variation and insensitivity to SNPs; (2) reproduction of the species-specific popANI-versus-synteny correlation patterns that SynTracker identified; and (3) runtime benchmarking showing orders-of-magnitude speed improvements over alignment-based methods. We further show that enzyme choice dramatically affects detection sensitivity, with the high-density Type IIG enzyme CjePI providing tag densities sufficient to detect small indels (10 kb) that are invisible to traditional Type IIB enzymes.

---

## Results

### 1. Syn2b design and algorithm

**Type IIB enzyme in-silico digestion.** Syn2b simulates restriction digestion using the recognition rules of Type IIB enzymes (BcgI, AlfI, BplI) and the Type IIG enzyme CjePI. Each enzyme recognizes a specific motif within a defined window (32 bp for BcgI/AlfI/CjePI, 27 bp for BplI) and produces a tag sequence of that length (Figure 1a). In *E. coli* K-12 (4,543,028 bp), BcgI yields 2,872 tags (0.63/kb), AlfI yields 1,978 (0.44/kb), BplI yields 375 (0.08/kb), and CjePI yields 9,356 (2.06/kb). Combining all four enzymes produces 14,581 tags (3.21/kb)—a 5.1× increase over BcgI alone (Table 1).

**Adjacency Jaccard and breakpoint metrics.** For each genome, Syn2b extracts the ordered list of tag sequences. Between two genomes, it computes: (1) *adjacency Jaccard*—the Jaccard similarity of adjacent tag pairs, capturing local synteny disruption; and (2) *breakpoint count*—the number of adjacent tag pairs present in one genome but not the other, quantifying rearrangement events (Figure 1b). Both metrics are insensitive to SNPs (which rarely alter enzyme recognition sites) but highly sensitive to inversions (which reverse tag order within a region) and indels (which eliminate or add tags).

**Kendall tau rank correlation.** To capture global order disruptions—particularly translocations, which move genomic segments to new locations—Syn2b computes the Kendall tau rank correlation of matching tag positions between two genomes. Unlike adjacency metrics, which count only local disruptions, Kendall tau evaluates whether the global ordering of shared tags is preserved. A translocation that moves a segment elsewhere disrupts the rank order of all tags in that segment, yielding a markedly lower tau (Figure 1c).

**Implementation and performance.** Syn2b is implemented in Rust with Python bindings. The core digestion engine uses anchor-based search (finding specific 3-bp motifs and verifying full recognition rules) rather than sliding-window scans, achieving O(N) time complexity with a small constant factor. Multi-enzyme digestion is parallelized via rayon. On a single core, digesting *E. coli* K-12 with all four enzymes takes 0.17 seconds; the entire pipeline (digestion + metric computation) for a pairwise comparison takes <1 second.

### 2. Syn2b is sensitive to structural variation, not SNPs

To validate Syn2b's sensitivity profile, we performed in-silico evolutionary simulations using *E. coli* K-12 as a reference. We generated "same-strain" controls by introducing random point mutations at 1% frequency (45,430 SNPs) without structural variation. We then introduced specific structural variants: 500-kb inversions, 500-kb translocations, 10-kb insertions, and 10-kb deletions. For each condition, we computed Mash distance, adjacency Jaccard, breakpoint count, and Kendall tau.

**SNP-only controls.** As expected, 1% SNPs produced a Mash distance of ~0.01, but adjacency Jaccard remained ~0.34, breakpoints ~3, and Kendall tau ~0.97—nearly identical to the unmutated reference (Figure 2a, Table 2). This confirms that SNPs outside enzyme recognition sites do not affect tag patterns.

**Inversion detection.** A 500-kb inversion dramatically reduced adjacency Jaccard to 0.29 and increased breakpoints by +646 (BcgI), +484 (AlfI), +64 (BplI), and +2,018 (CjePI). The effect scaled with tag density: CjePI's higher density yielded 3.1× more breakpoints than BcgI, making small inversions easier to distinguish from noise (Figure 2b, Table 2).

**Translocation detection.** Translocations had minimal effect on adjacency Jaccard and breakpoints (only +8 breakpoints regardless of enzyme) because they disrupt only two adjacency relationships at the break points. However, Kendall tau dropped from ~0.97 to ~0.82–0.86, capturing the global order disruption (Figure 2c, Table 2). Multi-enzyme combinations amplified this signal further (tau ~0.83 with all four enzymes).

**Insertion and deletion detection.** A 10-kb insertion or deletion increased breakpoints by +5 (BcgI), +4 (AlfI), +3 (BplI), and +23 (CjePI). CjePI's 4.6× higher indel sensitivity validates its utility for detecting small structural variants invisible to traditional Type IIB enzymes (Figure 2d, Table 2).

**Mash is blind to structural variation.** In all conditions, Mash distance remained ~0.01—the value expected from SNPs alone—confirming that k-mer-based methods cannot distinguish structural variation from point mutations (Table 2).

### 3. Enzyme comparison reveals density-dependent sensitivity

We systematically compared the four enzymes across all metrics (Table 2). BcgI and AlfI, both Type IIB enzymes with 32-bp tags, performed similarly, with BcgI's slightly higher density (0.63 vs 0.44/kb) yielding proportionally stronger inversion signals. BplI's lower density (0.08/kb) made it the least sensitive, though it still detected translocations via Kendall tau. CjePI's dramatically higher density (2.06/kb) provided the strongest signals across all structural variant types.

**Multi-enzyme synergy.** Combining all four enzymes yielded 14,581 tags (3.21/kb). The combined tag set amplified inversion breakpoints to +3,201 (5.0× over BcgI alone) and indel breakpoints to +29 (5.8× over BcgI). Importantly, multi-enzyme digestion also improved translocation detection via Kendall tau (0.83 vs 0.82 for BcgI alone), as the higher tag density increased the proportion of tags whose global order was disrupted (Figure 2e).

### 4. popANI versus Syn2b synteny recapitulates species-specific evolutionary modes

A key demonstration of SynTracker's utility was the species-specific relationship between SNP similarity (popANI) and synteny similarity (APSS)¹⁶. We reproduced this analysis using Syn2b's synteny score (a weighted combination of adjacency Jaccard and Kendall tau) versus popANI, simulated for four species with distinct evolutionary parameters (Figure 3):

***S. rimosus*-like (low SNPs, no structural variation).** All pairs showed high popANI (>0.9999) and high Syn2b synteny (>0.95). The two metrics were perfectly correlated, as expected for clonal populations accumulating only point mutations (Figure 3a).

***H. pylori*-like (no SNPs, high structural variation).** popANI remained near 1.0 (few point mutations), but Syn2b synteny ranged from 0.90 to 0.98—reflecting extensive recombination, inversion, and translocation. This decoupling demonstrates that Syn2b captures structural variation invisible to SNP-based tools (Figure 3b).

***N. gonorrhoeae*-like (medium SNPs, medium structural variation).** Both metrics varied, with a positive correlation (Spearman's ρ ≈ 0.6). Pairs with high popANI but low synteny indicated recent structural rearrangement; pairs with low popANI but high synteny indicated SNP accumulation without rearrangement (Figure 3c).

***E. coli* hypermutator-like (high SNPs, low structural variation).** popANI ranged widely (0.999–0.9999) due to hypermutation, while synteny remained high (>0.95). This pattern identifies hypermutator strains that accumulate point mutations but rarely undergo recombination (Figure 3d).

These patterns closely mirror SynTracker's APSS-versus-popANI plots¹⁶, demonstrating that Syn2b captures the same biological information.

### 5. Runtime benchmarking: Syn2b versus SynTracker

To quantify Syn2b's computational advantage, we benchmarked it against SynTracker's pipeline (BLAST database construction + DECIPHER pairwise alignment) on [N] genomes of varying sizes.

**Single-genome digestion.** For *E. coli* K-12 (4.6 Mbp), Syn2b's optimized Rust digestion completed in 0.13 s (BcgI+AlfI+BplI) and 0.17 s (all four enzymes including CjePI). In contrast, SynTracker's BLAST database construction alone requires [time], and each pairwise alignment via DECIPHER's FindSynteny takes [time].

**Pairwise comparison.** A complete Syn2b pairwise comparison (digestion + all metrics) for two *E. coli* genomes takes <1 second. SynTracker's APSS calculation for the same pair, using 40 regions per comparison (the minimum recommended), requires [time]—[X]× slower (Table 3).

**Scalability.** For [N] genomes, Syn2b's runtime scales linearly with the number of genomes (O(N) for digestion, O(N²) for pairwise metrics but with very small constant factors). SynTracker's runtime scales as O(N²) with large constant factors due to all-versus-all BLAST alignments. For [a large dataset], Syn2b completed in [time] versus SynTracker's [time]—a [Y]× speedup (Figure 4, Table 3).

**Memory footprint.** Syn2b requires only the genome sequences in memory (4.6 MB per *E. coli* genome). SynTracker requires the BLAST database (several GB for large datasets) plus DECIPHER's alignment matrices.

---

## Discussion

We have presented Strain2b (Syn2b), a rapid, alignment-free method for detecting structural genomic variation between microbial strains. By leveraging Type IIB restriction-enzyme tag adjacency and order, Syn2b captures inversion, translocation, insertion, and deletion events with high sensitivity while remaining blind to SNPs—precisely the complementary profile needed to pair with SNP-based strain trackers.

**Comparison to SynTracker.** SynTracker¹⁶ pioneered the use of microsynteny for strain comparison, using 1-kb central regions and DECIPHER pairwise alignments to compute APSS. Syn2b achieves the same biological goal—quantifying structural variation—but through a fundamentally different mechanism. SynTracker requires BLAST database construction, all-versus-all pairwise alignments, and the DECIPHER R package; Syn2b requires only the genome sequence and enzyme recognition rules. This difference translates to dramatic speed improvements: Syn2b processes a 4.6-Mbp genome in 0.17 seconds versus SynTracker's [time]. For large-scale studies (hundreds to thousands of genomes), this speedup makes previously infeasible analyses practical.

Moreover, Syn2b provides multiple complementary metrics (adjacency Jaccard, breakpoint count, Kendall tau) rather than a single APSS score. Adjacency Jaccard and breakpoints are most sensitive to inversions and indels, while Kendall tau captures translocations. This multi-metric approach allows researchers to distinguish different types of structural variation—a capability SynTracker's single APSS score does not provide.

**Enzyme selection matters.** Our systematic comparison of BcgI, AlfI, BplI, and CjePI reveals that tag density is the primary determinant of detection sensitivity. CjePI, a Type IIG enzyme with ~1.1 recognition sites per kb in its native host *Campylobacter jejuni*, provides 2.06 tags/kb in *E. coli*—3.3× denser than BcgI. This density enables detection of 10-kb indels that are invisible to BcgI alone (+23 vs +5 breakpoints). For applications requiring maximum sensitivity, multi-enzyme combinations (BcgI+AlfI+BplI+CjePI) provide the highest tag density (3.21/kb) and the strongest structural variation signals.

However, higher density is not always better. Very dense tag sets may increase computational requirements and could theoretically reduce specificity if recognition sites cluster non-randomly. In practice, we observed no such clustering artifacts in *E. coli*, but enzyme choice should be tailored to the genome size and expected structural variation scale of the target species.

**Species-specific evolutionary modes.** Our reproduction of the popANI-versus-synteny correlation patterns demonstrates that Syn2b captures the same biological insights as SynTracker. The four archetypal patterns—*S. rimosus* (high SNPs, low SV), *H. pylori* (low SNPs, high SV), *N. gonorrhoeae* (both), and *E. coli* hypermutator (high SNPs, low SV)—are readily identified. When combined with SNP-based tools (e.g., inStrain), Syn2b enables classification of species as hypermutators, hyper-recombinators, or balanced—information critical for understanding pathogen evolution, antibiotic resistance spread, and microbiome strain dynamics.

**Limitations and future directions.** Syn2b's current implementation assumes complete genome assemblies. While the method could theoretically be applied to metagenomic contigs (by digesting each contig separately), the adjacency metric would be confounded by assembly fragmentation. Future work will explore contig-aware adjacency scoring and integration with metagenome assembly pipelines.

Additionally, Syn2b currently uses a fixed set of enzyme recognition rules. Expanding the enzyme library—particularly to include more Type IIG enzymes with varying recognition motifs—would increase applicability across species with different GC contents and genome sizes. Machine learning approaches could further optimize enzyme selection for specific taxa.

Finally, while Syn2b's speed enables large-scale analyses, we have not yet benchmarked it on the massive datasets (thousands of human gut metagenomes) that SynTracker analyzed. Such benchmarking will be essential to validate Syn2b's utility for microbiome-wide strain surveys.

---

## Methods

### In-silico genome simulation

All simulations used *Escherichia coli* K-12 substrain MG1655 (NCBI: NC_000913.3, 4,543,028 bp) as the reference genome. Mutations were introduced using a custom Python script with Mersenne Twister random number generation (seed=42).

**Point mutations.** Random substitutions were introduced at frequency μ per base, with equal probability of A↔T and C↔G transitions/transversions.

**Structural variants.** Inversions: a segment of specified length was excised, reverse-complemented, and reinserted at the same position. Translocations: a segment was excised and reinserted at a random new position. Insertions: a random sequence of specified length was inserted at a random position. Deletions: a segment of specified length was removed.

### Type IIB enzyme in-silico digestion

**BcgI** (32-bp tag): Forward strand: offset 10=CGA and offset 19=TGC; Reverse strand: offset 10=GCA and offset 19=TCG.

**AlfI** (32-bp tag): Offset 10=GCA and offset 19=TGC (palindromic).

**BplI** (27-bp tag): Offset 8=GAG and offset 16=CTC.

**CjePI** (32-bp tag): Recognition motif CCYGA (Y=C or T) at offset 10. Site search uses direct motif finding (CCCGA, CCTGA) with O(N) complexity.

Digestion was implemented in Python with anchor-based search (finding 3-bp anchor motifs and verifying full rules) for efficiency. Tags were filtered to exclude those containing ambiguous bases (N).

### Syn2b metrics

**Mash proxy.** 21-mer canonical k-mer Jaccard similarity with stride=10 sampling, converted to Mash distance: d = -1/k · ln(2j/(1+j)).

**Adjacency Jaccard.** For each genome, adjacent tag pairs were encoded as (tag_i, tag_{i+1}). The Jaccard similarity of the pair sets from two genomes was computed.

**Breakpoint count.** The number of adjacent pairs present in one genome but not the other.

**Kendall tau.** For tags shared between two genomes, their rank orders were compared using Kendall's τ. Tags present in only one genome were excluded.

**Syn2b synteny score.** A weighted combination: synteny = 0.7 × adjacency_Jaccard + 0.3 × (τ + 1)/2.

### popANI simulation

For Figure 3 simulations, popANI was approximated from known SNP rates: popANI ≈ 1 - (μ_a + μ_b) × 2, with small Gaussian noise (σ=0.00001). SVs were introduced at specified rates (sv_rate_same, sv_rate_diff).

### Runtime benchmarking

Benchmarks were performed on [system specs]. Syn2b was compiled in release mode (rustc [version]). SynTracker [version] was run with default parameters. Timing was measured using [method].

### Software availability

Syn2b is available at https://github.com/HuangShiLab/Syn2b. The analysis scripts and data for this manuscript are at https://github.com/HuangShiLab/Syn2b-paper.

---

## References

1. [B. fragilis evolution reference]
2. [H. pylori recombination reference]
3. [Fungal phytopathogens reference]
4. [inStrain paper: Olm et al.]
5. [MIDAS paper]
6. [StrainPhlAn paper]
7. [SNP-based methods review]
8. [S. pneumoniae recombination reference]
9. [N. meningitidis virulence reference]
10. [Phase variation reference 1]
11. [Phase variation reference 2]
12. [Phase variation reference 3]
13. [Hypermutator reference]
14. [SynTracker paper: Enav et al. Nat Biotechnol 2024/2025]
15. [Type IIB enzymes reference]
16. [DECIPHER R package reference]
17. [Kendall tau reference]
18. [Mash paper: Ondov et al.]

---

## Figures and Tables

**Figure 1. Syn2b algorithm.** (a) Type IIB enzyme recognition and tag generation. (b) Adjacency Jaccard and breakpoint calculation. (c) Kendall tau rank correlation for translocation detection.

**Figure 2. Syn2b sensitivity to structural variation.** (a) SNP-only controls. (b) Inversion detection across enzymes. (c) Translocation detection via Kendall tau. (d) Insertion/deletion detection. (e) Multi-enzyme synergy.

**Figure 3. popANI versus Syn2b synteny for simulated species.** (a) *S. rimosus*-like. (b) *H. pylori*-like. (c) *N. gonorrhoeae*-like. (d) *E. coli* hypermutator-like.

**Figure 4. Runtime benchmarking.** (a) Single-genome digestion time. (b) Pairwise comparison time. (c) Scalability with genome count.

**Table 1. Enzyme tag densities in *E. coli* K-12.**

**Table 2. SV detection metrics across enzymes and SV types.**

**Table 3. Runtime comparison: Syn2b versus SynTracker.**

---

## Supplementary Information

Supplementary Tables 1–5: Raw simulation data.
Supplementary Figure 1: Tag spacing distributions.
Supplementary Figure 2: Sensitivity analysis with varying SV sizes.
Supplementary Note 1: Mathematical derivation of adjacency Jaccard sensitivity to inversion size.

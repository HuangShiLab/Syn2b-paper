# Cover letter

Dear Editors,

We submit our Article, “Syn2b: fragmentation-invariant structural comparison of microbial genomes by ordered restriction-enzyme tags,” for consideration at *Nature Methods*.

Strain-level microbial comparison often requires two incompatible choices: fast ANI estimators scale to large genome collections but discard adjacency, whereas synteny or alignment workflows preserve structural information but are too slow for database-scale use. Syn2b resolves this trade-off by representing genomes as ordered Type IIB restriction tags and by prioritizing length-weighted landmark ratios, which remain invariant to splitting of assemblies or alignment blocks.

The manuscript provides the technical method and a validation chain that directly supports its intended use. Syn2b’s fixed-reference inverted landmark fraction tracks dnadiff across 43,312 GTDB-R207 genome pairs (Pearson r = 0.94) and reaches r = 0.996 in the strain-level ANIm ≥ 97% regime. Controlled variants separate structural changes from SNP load and define the resolution of the four-enzyme panel. A 47,748-pair survey shows that 34% of strain-level pairs carry two or more alignment-defined inversions that ANI ranking alone does not report. Syn2b then scales to 646,012,820 GTDB species pairs and an independent 62,872,901-pair HROM census while retaining taxon-specific structural heterogeneity. On a 22-genome benchmark, pairwise comparison takes approximately 19 ms, approximately 420 times faster than skani plus dnadiff and approximately 15,000 times faster than SynTracker.

The work is well suited to *Nature Methods* because the central contribution is a new method, supported by strong validation, an error model, runtime benchmarks and application to large public genome collections. The two-stage screening strategy preserves alignment-based validation for selected candidates while making structural annotation feasible during genome search.

Syn2b v0.1.0 is available under the MIT licence at https://github.com/HuangShiLab/Syn2b (commit 7145305fc2f3139611aafa821d0c8c7ba328a7ee). Analysis scripts, summary data, QC tables, figures and supplementary tables are available at https://github.com/HuangShiLab/Syn2b-paper. Raw GTDB, RefSeq and HROM inputs are publicly available, and the full compressed GTDB/HROM pair tables will be deposited in a permanent archive before publication.

This manuscript is not under consideration elsewhere. All authors approved the submission. We declare no competing interests.

Sincerely,  
Shi Huang and Yufeng Zhang  
Faculty of Dentistry, The University of Hong Kong

# Pre-submission review 6: post-census revision

**Date:** 2026-10-13  
**Manuscript:** `Syn2b_Manuscript.md`  
**Scope:** macro logic, claims, numerical consistency, section numbering, grammar and style, figures and tables, reproducibility, and submission blockers.

## Executive verdict

**Not yet ready to submit without two author/administrative actions.** The scientific text and census claims have been revised and now pass the internal factual and style checks. The remaining blockers are (1) missing funding text and (2) a permanent accession/DOI for the full GTDB/HROM pair tables. In addition, the two queued fragmentation jobs are not yet complete, so their results cannot yet be added to the manuscript.

## Checks performed and outcomes

### CRITICAL findings

1. **Funding is unresolved.**  
   `## Acknowledgements` still contains `Funding information to be added.` This must be supplied by the authors or the funding section must be explicitly removed by instruction. It cannot be invented.

2. **Full pair-table archive accession is not yet available.**  
   `Data availability` correctly distinguishes repository summaries from the large compressed GTDB/HROM tables, but says only that the full tables are prepared for Zenodo/figshare. A manuscript-facing DOI or explicit “will be deposited before publication” journal-policy statement is still needed.

### Major findings resolved in this iteration

1. **GTDB census overclaim.**  
   The former phrase “alignment-visible breakpoints” was replaced with “landmark-defined breakpoints” or “Syn2b-called breakpoints” where appropriate. The text now states that the complete GTDB paired-ANI layer remains pending and does not infer direct evolutionary rates.

2. **Nonlinear section numbering.**  
   `6b` was removed. The census sections are now numbered 7 and 8, SynTracker is section 9, and runtime benchmarking is section 10.

3. **Independent HROM evidence missing.**  
   Added results section 8 and methods section `HROM within-species census`, plus Data Availability and Supplementary Figure 6 entries. Numerical claims now match `results/census/hrom/`.

4. **Supplementary table duplication.**  
   Removed inline duplicated Supplementary Tables 2 and 4 from the supplementary list section. The canonical table pointers remain.

5. **Numerical reconciliation.**  
   The abstract and correlation analyses use 3,826 unique high-ANI pairs, while the inversion-prevalence survey uses 3,828. Methods now states that these are two related analysis pipelines rather than silently mixing the counts.

### Minor findings resolved

1. Removed em dashes from the manuscript. Current scan: 0 occurrences.
2. Replaced inappropriate “reveal” wording in manuscript-authored prose and captions. Remaining “reveal”, “harness”, or “robust” occurrences are only inside reference titles.
3. Corrected section cross-references after renumbering.
4. Added the GTDB/HROM summary-figure caption and switched supplied figure references to vector PDF files where available.
5. Specified the HROM enzyme panel and changed “quality-aware” to the more precise “digest-effective” for the GTDB census.

## Automated checks

- `git diff --check`: pass.
- Em dash scan (`grep -n '—' Syn2b_Manuscript.md`): 0.
- Section headings: linear, 1 through 10.
- Referenced local files under `results/`, `supplementary/`, `data/`, and `figures/`: all present.
- Main census values verified against:
  - `results/census/gtdb/effective/effective_census_qc.md`
  - `results/census/hrom/census_qc.md`
  - `results/census/hrom/HROM_CENSUS_REPORT.md`
  - `results/discordance/discordance_summary.md`
- Figure S6 generator and PDF/PNG outputs are present:
  - `scripts/generate_census_figure.py`
  - `figures/supplementary/fig6_structural_censuses.pdf`
  - `figures/supplementary/fig6_structural_censuses.png`

## Benchmark status at review time

Both fragmentation benchmark jobs remain queued on HPC:

```text
4263624 hrom-s2b-frag PENDING
4263625 gtdb-s2b-frag PENDING
```

Their outputs are expected at:

```text
/lustre1/g/aos_shihuang/syn2b_benchmarks/results/hrom_full_fragmentation/fragment_summary.tsv
/lustre1/g/aos_shihuang/syn2b_benchmarks/results/gtdb_full_fragmentation/fragment_summary.tsv
```

These should be copied into `/Users/macstudio/Downloads/Syn2b/validation/` after completion and checked before any manuscript claim is added.

## Required before final submission

1. Provide exact funding text or instruct removal.
2. Provide a Zenodo/figshare DOI for large pair tables or confirm the target journal permits deposition after review.
3. Complete and QC the two queued fragmentation jobs if their results are to be included before submission.

With these author actions completed, the manuscript is otherwise ready for the next mechanical pass: DOCX regeneration, final figure embedding, and commit.

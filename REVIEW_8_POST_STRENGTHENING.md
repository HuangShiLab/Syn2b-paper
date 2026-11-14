# Review 8: strengthened Nature Methods package audit

**Date:** 2026-11-14  
**Scope:** content strengthening after Review 7  
**Skills applied:** `supervisor-review`, `nature-reviewer`, `nature-statistics`, `nature-polishing` and `anti-defensive-writing`.

## Verdict

**Manuscript content is strengthened and initial-submission ready if the venue permits delayed deposition of the full pair tables.** The complete release remains blocked only by two external identifiers: the Syn2b Zenodo software DOI and the permanent GTDB/HROM pair-table DOI. The GTDB E. coli pair-count shortfall also remains a pre-deposition QC item.

## Content added

1. **HROM taxonomy mapping.** Representative species labels were joined to 367 clusters with at least 1,000 metric-complete pairs. The paper now names high- and low-prevalence taxa and cites Supplementary Table 12.
2. **HROM ANI stratification.** Structural prevalence is now reported at skani ANI >=97% and >=99%, not only in the full HROM census.
3. **Census sensitivity.** Added 10,000-replicate cluster bootstrap and leave-one-largest-cluster-out analyses for GTDB and HROM in Supplementary Table 13.
4. **Assembly-quality sensitivity.** Added the HROM high-quality species-mode sensitivity analysis and Supplementary Table 14.
5. **Metric-complete denominator.** Identified four HROM catalogue pairs with no shared landmarks and no finite structural metric. The catalogue retains 62,872,901 pairs; structural prevalence uses 62,872,897 metric-complete pairs.
6. **Recent method positioning.** Added ntSynt and KmerAperture citations and differentiated Syn2b's pair-level Type IIB landmark orientation metric.

## Release work completed

- Added GitHub Actions CI to the Syn2b repository.
- Verified the release build and 152 tests locally and in CI.
- Created and pushed Syn2b tag `v0.1.0`.
- CI commit: `248441eabb2c2c72c8b3a85da65d88577839b63e`.
- CI status on the tagged commit: success.

## Technical audit

- Abstract: 138 words.
- Results plus Discussion: approximately 2,195 words.
- Main figures: 6; main tables: 0.
- References: 26.
- Supplementary tables: 14.
- Em dashes in manuscript: 0.
- Duplicate SynTracker provenance sentence: removed.
- HROM taxonomy and robustness tables: present in the main and submission package.
- DOCX files were regenerated and page-render inspected.

## Remaining external actions

1. Link or create the Syn2b Zenodo software DOI for tag `v0.1.0`.
2. Resolve the GTDB *E. coli* 7,571,899-pair shortfall before archiving.
3. Deposit full GTDB/HROM pair tables and replace the deposition promise with a DOI.
4. Enter funding and ORCID information in the journal portal.

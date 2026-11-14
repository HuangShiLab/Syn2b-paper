# Nature Methods initial-submission checklist

## Article fit and format

- Article type: **Article** (new method with validation).
- Abstract: **138 words** (Nature Methods limit: up to 150; unreferenced).
- Main text excluding abstract, Online Methods, references, declarations and figure legends: **about 2,200 words** (Article target: 3,000; up to 5,000 with editorial discretion).
- Main display items: **6 figures, 0 main tables** (within the 6 figure/table limit).
- References: **26** (guideline: up to 50).
- Section order: Abstract; Introduction without heading; Results with topical subheadings; Discussion without subheadings; Online Methods; References; Data availability; Code availability; Figures; Supplementary Information; Author contributions; Competing interests.
- Citations are superscript Arabic numerals in first-appearance order.
- Figure citations occur in sequence as Fig. 1 through Fig. 6.
- References use Nature-style `Author, A. et al.` formatting.

## Reproducibility

- Method code: Syn2b v0.1.0, MIT licence, tag `v0.1.0`, commit `248441eabb2c2c72c8b3a85da65d88577839b63e`.
- Analysis repository: https://github.com/HuangShiLab/Syn2b-paper (tag `paper-v1`).
- Source inputs: GTDB-R207, NCBI RefSeq, published SynTracker cohorts and HROM.
- Per-pair validation table, runtime tables, validation summaries, census summaries, HROM taxonomy prevalence, census sensitivity, assembly-quality sensitivity, QC tables, notes and figure files are present.
- Data and code availability statements are separate and located after References.
- Compressed full GTDB/HROM pair tables are declared for permanent archival before publication. Versioned source is available at Syn2b tag `v0.1.0`; the Zenodo software DOI remains to be added.

## Submission files

- Manuscript: `submission/NatureMethods/Manuscript_NatureMethods.docx`
- Supplementary Information: `submission/NatureMethods/Supplementary_Information.docx`
- Cover letter: `submission/NatureMethods/Cover_Letter.md`
- Figure files: `submission/NatureMethods/Fig_1.pdf` through `Fig_6.pdf`
- Supplementary tables 1–14: `submission/NatureMethods/supplementary_tables/`
- Figure manifest: `submission/NatureMethods/figure_manifest.tsv`
- Package checksum manifest: `submission/NatureMethods/SHA256SUMS`

## Author actions at submission

1. Enter ORCID and affiliation details in the Nature submission portal.
2. Enter funding information in the grant/funding section of the portal before final submission or acceptance, as requested by the journal. Funding is omitted from the manuscript rather than represented by a placeholder.
3. Confirm the cover letter and author contributions.
4. Link the GitHub `v0.1.0` release to Zenodo or create the software archive manually, then replace this item with the DOI.
5. If the editor requests the Nature Portfolio reporting summary, complete it from the statistical and reproducibility information in Online Methods.
6. Deposit the full pair tables or provide the archive DOI when available; the manuscript commits to deposition before publication.

## External release blockers

1. Software archive DOI for Syn2b `v0.1.0`: pending Zenodo release linking.
2. Permanent GTDB/HROM pair-table DOI: pending deposit.

These blockers do not alter the manuscript's scientific claims, but the complete
computational release remains blocked until both DOIs resolve.

## Non-applicable items

- No new human participants or animals were studied.
- No new pathogen isolates were generated or released.
- No clinical trial, plant science or crystal-structure reporting summary applies.

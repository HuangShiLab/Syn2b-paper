# Review 7: Nature Methods submission audit

**Date:** 2026-10-27  
**Manuscript:** `Syn2b_Manuscript.md`  
**Target article type:** *Nature Methods* Article  
**Skills applied:** `nature-writing`, `nature-polishing`, `nature-reviewer`, `nature-shared` main-text discipline, `pre-submission-reviewer`, and `anti-defensive-writing-en`.

## Verdict

**Ready for initial submission after the authors enter grant/funding information in the Nature submission portal.** The manuscript now meets Nature Methods Article formatting and reproducibility requirements. The only substantive external item is a permanent archive DOI for the full pair tables, which may follow before publication according to the current Data availability statement.

## Findings resolved

### CRITICAL

1. **Abstract exceeded the 150-word limit.**  
   Reduced from 180 to 136 words and removed references.

2. **Main display items exceeded the six-item limit.**  
   Removed Tables 1–3 from the main manuscript and retained them as Supplementary Tables 7–9. The manuscript now has exactly six main figures and no main tables.

3. **Article structure did not match Nature Methods.**  
   The unheaded Introduction is now represented by an HTML production comment; Results has topical subheadings; Discussion has no subheadings; Methods is labelled Online Methods; Data and Code availability are separate.

4. **Funding placeholder prevented submission.**  
   Removed `Funding information to be added.` The manuscript no longer contains a fabricated or unresolved funding statement. Grant details should be entered in the Nature portal, as recorded in the submission checklist.

5. **Citation order was nonsequential.**  
   Reordered the reference list and updated citations so first appearance is 1 through 24.

### MAJOR

1. **Main text exceeded the standard Article target.**  
   Compressed the evidence chain from approximately 4,900 to approximately 2,400 words without removing the core validation: controlled variants, enzyme ablation, dnadiff validation, strain-level discordance, GTDB census, HROM census, cohort application and runtime benchmark.

2. **Reproducibility lacked exact method-code provenance.**  
   Code availability now gives Syn2b v0.1.0, MIT licence, URL and commit hash. The code commit has been pushed to GitHub.

3. **Supplementary tables were incomplete for the compressed main text.**  
   Added machine-readable Supplementary Tables 7–11 for enzyme density, SV detection, runtime scaling and GTDB validation bands.

4. **Human-data provenance was implicit.**  
   Data availability now states that SynTracker data came from published cohorts and that no new human participant data were generated.

### MINOR

1. Changed `Figure` to Nature-style `Fig.` in main citations and captions.
2. Changed `Supplementary Figure` to `Supplementary Fig.`.
3. Reformatted `et al.` to Nature style.
4. Removed defensive-only wording while retaining claim-defining limitations.
5. Expanded author contributions into separable contributions without changing scientific credit.

## Technical checks

- Abstract: 136 words.
- Main text: approximately 2,400 words.
- Main figures: 6.
- Main tables: 0.
- References: 24.
- Em dashes: 0.
- Funding placeholder: 0.
- Citation first-appearance order: 1–24.
- All local figure, results, supplementary and data paths named in the manuscript exist.
- All six main display files are vector PDFs.
- Code commit `7145305fc2f3139611aafa821d0c8c7ba328a7ee` is now present on `origin/main` of the public Syn2b repository.

## Remaining external actions

1. Enter funding information in the Nature Methods submission portal.
2. Provide ORCID and institutional submission metadata.
3. Deposit the compressed full pair tables before publication and replace the Data availability promise with the DOI.
4. Complete the journal reporting summary if requested by the editor.

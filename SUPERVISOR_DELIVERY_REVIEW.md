# Supervisor delivery review after repairs

## Verdict

`LOCAL PACKAGE READY - EXTERNAL DOI DEPOSIT BLOCKED`

## Gate results

| Gate | Result | Evidence |
|---|---|---|
| Ledger goal and artifact ownership | PASS | `DELIVERY_LEDGER.tsv` |
| Citation first-appearance order | PASS | `Syn2b_Manuscript.md`; ntSynt and KmerAperture are references 20-21 |
| Assembly-quality stability evidence | PASS | `results/census/hrom/hrom_quality_mode_stability.tsv`; `supplementary/Supplementary_Table_14_hrom_quality_species_modes.tsv` |
| Outgoing package hygiene | PASS | Review 7 and Review 8 removed from `submission/NatureMethods` |
| Package integrity | PASS | `submission/NatureMethods/SHA256SUMS`; all 29 listed deliverables verify |
| Software version | PASS | Syn2b `v0.1.0`, commit `248441eabb2c2c72c8b3a85da65d88577839b63e`, CI passed |
| Analysis version pin | PASS after commit | `paper-v1` tag points to this package commit |
| Data archive | BLOCKED | Permanent GTDB/HROM pair-table DOI requires external deposit |
| Software archive | BLOCKED | Zenodo DOI requires external release linking |

## Remaining actions

1. Create the Zenodo DOI for Syn2b `v0.1.0`.
2. Deposit and version the compressed GTDB/HROM pair tables.
3. Replace the pending-deposit wording in Data availability with resolving DOIs.

No scientific claim was expanded to compensate for the missing external DOI.

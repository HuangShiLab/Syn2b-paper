# MOCK REVIEW — Syn2b manuscript (Syn2b_Manuscript.md)

Review date: 2026-09-14. Manuscript state: commit `93410da` (post REVIEW_2 fixes).
Reviewers: three independent passes — (1) text/internal consistency, (2) figures &
tables, (3) claim-by-claim recomputation from committed data.

## Verdict

**Major revision before submission.** The quantitative core of the paper is sound:
the headline correlations (r = 0.9355 held-out; r = 0.9960, slope 1.006, SD 0.0135
at ANIm ≥ 97%, n = 3,826), the band tables, Table 1 tag counts, runtime/speedup
tables, and cohort statistics were all reproduced **exactly** from committed data.
The problems are concentrated in: two figure-level bugs that invert or mislabel
results (Figure 2 Mash panel, Figure 5a wrong enzyme panel), a revived REVIEW_2-B3
violation ("3,099 held-out pairs" is the predicted-ANI-selected set), several
stale/contradicted numbers in §3/§6/Discussion, and Methods completeness.

---

## Major findings

### M1. Figure 2 Mash panel: control is compared to itself; the panel contradicts its own caption

`scripts/generate_figure2_rust.py:200-221` computes the SNP-only control's Mash
distance as `mash_distance(control, control)` → exactly 0 (printed as "-0.00e+00"),
while SV variants are compared against the control. The intended claim — "Mash
reports no additional signal under any structural variant" (manuscript :564) — is
**inverted** by the figure: the translocation bar (6.05e-06) is 145× the control
bar (0.0), which a referee will read as "Mash detects translocations".

The correct display: compute the control's Mash distance against the *unmutated*
reference (≈ 0.01 for the 1% SNP background — the signal a k-mer method should
see), then show that SV conditions add nothing on top (all ≈ 1e-7, i.e. ~4 orders
of magnitude below the SNP baseline). Plot on a log scale that contains the SNP
baseline.

Related: Methods (:443-444) says "stride=10"; `canonical_kmers` enumerates every
position (stride 1). Align text and code, and cite the Mash reference (ref 18,
currently never cited).

### M2. Figure 5a times the legacy enzyme panel, not the production panel

`scripts/generate_manuscript_figures.py:379-392` times
`["BcgI","AlfI","BplI","CjePI","All 4"]` with "All 4" → "BcgI,AlfI,BplI,CjePI".
The production panel everywhere else in the paper is **BcgI+AlfI+AloI+FalI**;
AloI and FalI are never timed and BplI/CjePI are not in the panel. The "~45 ms"
in §1/Abstract and "30 ms (BplI) to 44 ms (All 4)" in §7 (:348-350) therefore
describe a measurement of the wrong panel. Also note the range claim is wrong
regardless of panel: BcgI alone measured 51.2 ms, above "All 4" at 44.0 ms.

Fix: re-time with `-e BcgI,AlfI,AloI,FalI`, regenerate Figure 5a, harmonize the
three text mentions (§1, §7, Abstract). Optionally archive the timing output so
the numbers are committed rather than measured live at figure-build time.

### M3. §6 misstates the E. coli hypermutator cohort's ANI regime

:328-329: "*E. coli* hypermutator pairs sit at lower ANI (94.5–96.5%) with low
breakpoint counts (median = 0)". The cohort skani data
(`data/syntracker_validation/skani/skani_Escherichia_coli_hypermutator.tsv`) show
**99.93–100.00% ANI** (median 99.99), median breakpoints = 0. The 94.5–96.5% range
is the *H. pylori* between-host cloud. As written, the sentence contradicts
Figure 4a (the E. coli points sit at ANI ≈ 100) and looks like a leftover from
the deleted simulation. Fix: near-clonal ANI ≥ 99.9% with median 0 breakpoints;
give the 94.6–96.5% range to H. pylori where it belongs.

### M4. "3,099 held-out pairs" is the predicted-ANI-selected set — REVIEW_2 B3 violated again

:292-295: "Among the 3,099 held-out pairs where dnadiff reports an inverted
fraction > 0.5, the fixed-reference version still agrees at r = 0.6826, whereas
the majority-frame version anticorrelates at r = −0.7438".

The numbers are verbatim from `results/metric_validation/high_ani_validation.tsv`
rows computed on `syn2b_inverted_fraction_high_ani_all.tsv` — the set **selected
on predicted ANI** that `results/gtdb50k/inverted_fraction_comparison_report.md:150-155`
explicitly says must not be quoted as a high-ANI result, and that REVIEW_2 B3
banned. The held-out set's saturated subset is 21,035 pairs, not 3,099. The
paragraph also silently switches sets (the r = 0.177 majority-frame number two
sentences earlier *is* from the held-out set) and contradicts :261-262 ("the
held-out set contains only two pairs at ≥97% ANIm").

Fix: recompute the saturated-subset analysis on `inverted_fraction_truth_four.tsv`
(held-out) and report those numbers, or name the actual set and justify it.

### M5. Discussion keeps the pre-correction "~9 ms per pair"

:388 "…orders of magnitude faster for large panels (~9 ms per pair)" contradicts
§7 and Table 4 (18.9 ms per *unique* pair, with the n² accounting explicitly
demoted). REVIEW_2 item E missed this instance. Fix: "~19 ms per unique pair".

### M6. §3 indel claims exceed the data; the resolution-limit citations point to nothing

- :185-189 "10-kb indels reliably recoverable, whereas BcgI alone … is less
  consistent": `results/supplementary_figure3_metrics.csv` shows insertions and
  deletions of 10/50/100 kb **all produce breakpoints = 0, scj = 0** — junction
  counts detect indels only when a landmark straddles the breakpoint. No BcgI-only
  indel data exist anywhere. The only measured indel signal is a drop in
  `shared_tags` (6,216 → 5,690 for 10 kb), which the text never names.
- The "~4-kb / ~8-kb resolution limits" are cited to Supplementary Note 2, which
  contains no such content; the numbers live only as a doc comment in the Rust
  source with no committed measurement.
Fix: rewrite to what the data show (indels detected via tag loss, not junction
counts), delete or data-support the BcgI comparison, and either commit the
resolution-limit measurement or cite the tool documentation instead of Note 2.
Same wrong pointer at :219 ("robust to the occasional missing or extra site
(Supplementary Note 2)").

### M7. §3 controlled claims are wider than the single measurement behind them

- :169-172 "Substitutions from 0.5% to 5% … produce zero structural junctions":
  only **1%** was ever run (`SNP_RATE = 0.01`, single replicate, seed 42).
- The control's zeros are **hardcoded, not measured**:
  `generate_figure2_rust.py:206-214` assigns `breakpoints: 0, scj_distance: 0`
  to the control instead of running `syn2b synteny` on the control pair. The
  paper's cleanest controlled claim rests on an asserted value.
Fix: run the control through the tool (and ideally a 0.5/1/2/5% sweep), report
what is measured.

### M8. §3 "400 kb" inversion does not exist

:175 "regardless of inversion size (400 kb or 500 kb)": the data
(`results/figure2_rust_metrics.csv`, Supplementary Table 3) contain 100 kb and
500 kb. Contradicts the Figure 2 caption (:559). Fix: "100 kb or 500 kb".

### M9. Methods: wrong genome length, and missing definitions for most reported metrics

- :423 "NC_000913.3, 4,543,028 bp" — the genome is 4,641,652 bp (§1, Table 1,
  Suppl. Table 3, and the bundled FASTA all agree). This value matches no K-12
  version.
- "Syn2b metrics" (:441-453) defines only Mash proxy, Adjacency Jaccard,
  breakpoint count, length-weighted inverted fraction. Results additionally report
  `scj_distance` (§3, §6, Figure 2 panel), `observable_fraction` (§6),
  `synteny_blocks` (§5), Kendall tau (§3), and the contig-term correction (:305-308)
  — none defined. "SCJ distance" is additionally described as "bp" (:336-338)
  though it is a count of operations.
- skani is used throughout §6/§7/Figure 4 with no Methods subsection, no version,
  no parameters, and no citation. Same for minimap2, MUMmer, BLAST.
- The high-ANI sample carrying the headline result is described only as "selected
  to cover the strain-level regime (ANIm ≥ 95%)" (:461) with no selection
  procedure. Given M4, this must be documented: source genomes, pair selection,
  deduplication, and why selection on predicted ANI is avoided.

### M10. Main-text Tables 2 and 3 are still title-only placeholders

:605-608. Table 2's content was withdrawn with the legacy prototype data and no
longer exists; Table 3 duplicates the §5 inline band table. A submission with
empty tables is not reviewable. Fix: populate (e.g., Table 2 = BcgI-only vs
4-enzyme × SV-type matrix from current Rust data; Table 3 = the §5 table) or
delete the entries.

### M11. Table 1 contains a false supporting sentence

:104 region: "The total is below the naive sum because recognition motifs occupy
overlapping sequence space" — the combined count 6,216 **equals** the naive sum
exactly (2,935+2,023+523+735); overlap only affects later landmark collapse
(`MIN_TAG_SEPARATION`), not digest counts. Delete or correct the sentence.

### M12. Figure 3b band identity is ambiguous vs the in-text ANI-band table

Figure 3b bands (80-85 n=12,152 … 95-100 n=404) are the TSV `band` column
(predicted-ANI bands); the Results table (:252-259) uses minimap2 ANIm bands
(80-85 n=1,850 …) on the same 43,312 pairs. Same labels, incompatible counts,
caption says only "identity band". Fix: state which identity measure defines each,
or use one banding throughout.

### M13. §5 closed-genome validation is glossed more strongly than the repo's own report allows

:297-303 quotes medians (4 matched inversions, 21.8 kb median distance, 31
dnadiff-only; all verified) but the underlying report
(`results/closed_inversions/JUNCTION_COORDINATE_REPORT.md`) adds: only 13.6% of
dnadiff boundaries have a Syn2b partner; only 40.4% of matched junctions are within
5 kb (p90 = 314.8 kb); greedy matching can pair different biological events. The
manuscript's "consistent with the sampling model" omits these caveats. Fix: report
the matched-fraction and within-5-kb percentages, state that 21.8 kb is the median
of per-pair medians, and temper the conclusion.

### M14. Figure 4 within-host claim has no statistical test, and the panel does not show it

§6 (:334-343): within-host H. pylori pairs have lower breakpoints/SCJ/higher
observable fraction (medians 0 vs 8; 11 vs 37; 0.995 vs 0.978 — verified). But no
test exists anywhere in `scripts/` (no Mann-Whitney/permutation). Figure 4d
overplots the 476 within-host pairs into a few dots at (≈100, 0) and the legend
sits on top of the between-host cloud; between-host pairs are colored by the
alphabetically-first isolate's host, which is arbitrary. Fix: add a permutation
test (report effect size + p), and rework panel d (within/between encoding,
jitter or inset).

---

## Minor findings

1. **Fig 4 legend overlaps** (b: legend over H. pylori cloud; d: legend on data) and
   Fig 4a: the E. coli median-0 cluster is invisible on a 0–30 axis.
2. **Fig 5b**: the "23.2" annotation sits on the curve line.
3. **Supplementary Figure 1 suptitle** contains literal `\*E. coli\*` asterisks
   (`generate_supplementary_figure1.py:87`).
4. **Supplementary Figure 1 caption arithmetic**: "6,215 intervals, median 1,533 bp,
   mean 2,980 bp" — if the intervals partitioned the 4.64 Mb genome, the mean must
   be 747 bp. The numbers reproduce from the TGT distances, so define what the
   intervals are (they are not a genome partition) or correct the caption.
5. **Fig 1 caption vs panels**: (b) caption says "order and orientation" but the
   panel shows order only; (a) panel title says "Type IIB/IIG", caption says "Type IIB".
6. **Fig 3c caption overclaims**: "whereas the fixed-reference raw_inverted_fraction
   removes the saturation" — panel c shows only the saturation; the fix is panel a.
7. **Fig 3b truncated y-axis** (0.85–1.0) exaggerates band differences; acceptable
   for correlations but note it.
8. **~2.3 GB RSS** (:359) has no committed source (no RSS column in any
   efficiency_v8 file); measure and archive or drop.
9. **Band table (:252-259)** rows sum to 43,310; the 2 pairs at ≥97% ANIm are not
   shown — add an "≥97 (n = 2, excluded)" footnote.
10. **"43,334 held-out pairs"** (:221, §4) vs 43,312 complete pairs and BcgI arm
    n = 41,485 — state complete-case n per arm.
11. **References**: refs 18 (Mash) and 19 (DECIPHER) are defined but never cited;
    skani, minimap2, MUMmer, BLAST are used with no reference. skani is the most
    serious (Figure 4 depends on it).
12. **Intro "the optimal landmark panel"** (:69) contradicts §4's own hedging
    (:236 "a high-performing … default rather than the absolute optimal"). Use
    "production panel".
13. **Suptitle style inconsistency**: Figure 2 and Suppl. Figs 1–3 embed
    "Figure N | …" titles in the image; Figures 1, 3, 4, 5 do not. Pick one.
14. **Discussion "library-scale surveys practical"** (:395) — largest benchmark is
    22 genomes; qualify.
15. **Placeholder text**: :702 "Funding information to be added."
16. **Suppl. Note 2 internal contradiction**: its per-pair SE table lists values
    below the stated 0.0205 floor (e.g., SE = 0.0153 at m = 5,000) — stale from an
    older fit.
17. **SCJ vs breakpoint_count ratios inconsistent on cohorts** (REVIEW_2 G2, still
    open): E. coli bp 0 / scj 24; H. pylori cross bp 8 / scj 37; S. rimosus bp 10 /
    scj 760. Explain or report one metric.
18. **Error-model internals**: Suppl. Note 2's application numbers (SD(z) = 1.006,
    95.3%, model SDs 0.0546/0.0874) do not reproduce from committed per-pair data
    (independent recomputation: 1.076/94.7%); the headline formula and coefficients
    do reproduce. Reconcile.
19. **"statistically indistinguishable"** (:233, fmh750 vs panel): CIs do not
    overlap (though the difference is negligible). Soften.
20. **Repo hygiene**: committed Word lock file `~$n2b_Manuscript.docx`; stale root-level
    `supplementary_note_1_fragmentation_principle.md` and
    `supplementary_table_3_h_pylori_simulation.md` (for the deleted simulation);
    `real_data_h_pylori.*` never renamed; Fig 3d depends on a sibling
    `../Syn2bANI-paper/data/gtdb_metadata` checkout; missing committed inputs
    (`high_ani_truth.tsv`, `data/gtdb50k_heldout_pairs.tsv`,
    `scripts/analyze_invfrac_error_model.py`) that the README/report cite for
    reproduction.
21. **Supplementary Tables 1–4 and Suppl. Figs 1–2 are never cited in the body**
    text (only listed in the SI section); Tables 2–4 likewise. Add citations.

## REVIEW_2 closure status

| Item | Status |
|---|---|
| A (Fig 4 on real cohorts) | Fixed — but new error M3 (E. coli ANI range) |
| B1 (mirroring) | Fixed, verified (r = 0.9355 / 0.1771 reproduce) |
| B2 (strain-range lead) | Fixed, verified exactly (n = 3,826, r = 0.9960) |
| B3 (don't quote predicted-ANI set) | **Violated again — M4** |
| B4 (drop n = 404 headline) | Fixed |
| B5 (error coefficients) | Fixed (formula + coefficients reproduce; note minor 18) |
| B6 (Fig 3d fabrication fallback) | Fixed (raises RuntimeError) |
| C (Rust Figure 2) | Partially — new errors M1, M7, M8 |
| D (panel justification → junction floor) | Fixed, verified (FracMinHash arms exact) |
| E (C(n,2) runtimes) | Partially — M2, M5 |
| F (real-data claim) | Fixed in text (files not renamed — minor 20) |
| G1 (banding sentence) | Fixed |
| G2 (SCJ vs breakpoints) | Not fixed — minor 17 |
| G3 (promised supplements) | Partially — files exist, M10 remains |
| G4 (figure quality) | Mostly fixed — M1, M2, minor 2/3 |

## Recommended fix order

1. M1 (Fig 2 Mash bug), M2 (Fig 5a panel), M7 (hardcoded control) — figure/data
   regenerations, 1 day, and they change what the paper claims.
2. M4, M3, M8, M9-genome-length — text corrections with recomputed numbers, hours.
3. M6, M7-sweep, M11 — decide what to measure vs soften.
4. M10, M9-definitions, M12, M13, M14 — Methods/table/figure polish.
5. Minor sweep in one pass.

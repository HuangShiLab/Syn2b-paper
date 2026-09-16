# REVIEW 5 — Nature Methods submission gate

Review date: 2026-09-17. Manuscript state: commit `7fa17d2` (post mock-review
closure, post retraction of Syn2bANI-derived statistics, post
landmark-fraction terminology fix).
Reviewer stance: Nature Methods editor + two referees (methods referee,
biology referee). Companion strategic plan:
[`APPLICATION_EXPANSION_PLAN.md`](APPLICATION_EXPANSION_PLAN.md).

---

## Verdict

**The methodological core is submission-grade. The application layer is not
yet.** The fragmentation principle is a genuine, general contribution; the
error model, dual-implementation cross-validation, and honest retraction
handling are strengths. What an NM referee will still object to is that the
flagship application (§6) surveys **47,748 pairs — 0.09% of GTDB R207's
within-species pair space (711,020,841 pairs)** — and that no analysis in the
paper demonstrates Syn2b recovering **known biology**. Both gaps are closable
in ~4–6 weeks with compute the lab already has (see the plan). With the census
and two biology vignettes added, the paper is a credible NM submission;
without them it reads as a strong Bioinformatics/AEM-level methods paper.

---

## What is already strong (keep, and lead with these)

1. **The fragmentation principle** (§2, Supplementary Note 1) — a general
   statement about a bias class in a field-standard tool (dnadiff/MUMmer
   breakpoint counts), with a formal derivation, three testable consequences,
   and empirical verification on 43,312 real pairs (dnadiff breakpoints vs
   contig count: ρ = 0.28, partial 0.34). This is the paper's NM card; the
   title already reflects it.
2. **Calibrated per-pair error model** — Var(err) = 1.504·p(1−p)/m + 0.0205²,
   fitted on held-out data, predicts the high-ANIm spread without refitting.
   Rare among comparable tools; make sure Discussion says why it matters
   (downstream users get a z-score per call).
3. **Dual independent implementations agreeing on controls** (Syn2b junction
   count and Syn2bANI v0.1.1 chain counter both: 0 self, 2 for EDL933 vs
   Sakai) — referees who check code will find this unusual and reassuring.
4. **The orientation-convention honesty** (§5 caveat, Limitations): the
   inverted fraction on drafts is a *tracking* metric; event counts are
   convention-free. Preempts the most damaging technical referee attack.
5. Runtime story: ~19 ms/pair vs 282 s/pair (SynTracker) and ~8 s/pair
   (dnadiff), with archived HPC measurements.

---

## Major vulnerabilities (referee-attack model)

### V1. Application scale — the decisive gap

§6's headline ("34% of ANIm ≥ 97% pairs carry ≥2 inversions") comes from a
47,748-pair sample that was **enriched for close relatives and then filtered
on measured ANIm**. A referee will ask: (a) what is the prevalence in the full
database; (b) what are the confidence intervals; (c) how much of the estimate
is driven by which taxa. The full within-species census of R207 is
**711,020,841 pairs across 317,542 genomes — computable in ≈3,700 core-hours**
(≈10 node-days on 16-core nodes; see plan). Doing it converts the paper's
central number from a sample statistic into a census statistic and adds the
headline "first structural-variation census of a genome database".

Note the concentration: **top-10 species = 95.6% of all pairs** (E. coli
50.7%, S. aureus 12.0%, Salmonella 10.6%, K. pneumoniae 9.0%...). The census
is, in effect, a structural census of the priority pathogens. This is a
feature — but the paper must then report per-species rates (a species ×
rearrangement-prevalence matrix), which is also the most biology-rich figure
possible (see V2).

### V2. No demonstration against known biology

Every real-data result is "agrees with dnadiff" or "reproduces SynTracker's
within-host signal". Nowhere does Syn2b recover a *documented* rearrangement.
Cheapest fixes, in order of value:

- **LTEE** (Tenaillon et al. 2016; 264 complete genomes, 12 populations ×
  50,000 generations): known population-specific inversions; complete
  genomes → no orientation artifact; a within-population structural timeline
  is one panel. Compute: 35k pairs, minutes.
- **Salmonella enterica / E. coli–Shigella complete genomes**: rDNA-mediated
  inversions are textbook; showing ANI-blind, junction-visible structure in
  the two biggest census species anchors the census in known mechanism.
- **MTB as the negative control** (23.4M pairs in the census, expected
  structurally hyper-conserved): a species that *doesn't* decouple makes the
  species matrix credible.

### V3. The inverted-fraction channel concedes too much to drafts

The paper itself states the fixed-reference fraction is a tracking metric on
drafts. The rebuttal the paper currently lacks: a **complete-genomes-only
analysis at scale** (within-species pairs among RefSeq/GenBank complete
genomes — tens of millions of pairs; see plan Phase 1c), where orientation is
real biology and dnadiff concordance is interpretable. One main-text panel.

### V4. Flagging performance is modest and the screening framing needs a costed two-stage analysis

AUC 0.80, 30% sensitivity at 4.6% FPR. Referees will call this weak. Two
defenses, both cheap:
- Improve the feature: junction count normalized by shared tags, or a
  junction + observable-fraction composite; recompute AUC from the census
  (n = tens of millions) with proper uncertainty.
- Cost the two-stage workflow explicitly: screening all within-species pairs
  with Syn2b (~X core-hours) + dnadiff confirmation only on flagged pairs
  (Y s × Z pairs) vs dnadiff on everything (711M × 8 s ≈ 6,300 node-days).
  The >100× end-to-end saving is the point; per-event sensitivity is not.

### V5. Comparator fairness leftovers

- No minimap2(+PAF parsing) timing as the fast alignment alternative
  (cross-repo review A8); add or argue scope.
- skani v0.1.0 now disclosed in Methods — fine, but the census should ship
  with skani ≥0.2.1 ANI columns (resolves the open item at the same time).
- SynTracker 16-core configuration is documented; keep the "workflows not
  identical quantities" table note.

### V6. Pre-submission text items

- Abstract: the 34% sentence needs the sample size and sampling caveat
  ("in a 47,748-pair survey"; becomes moot/much stronger after the census).
- §6: add cluster-bootstrap or per-species CIs; state that prevalence is
  pair-weighted and dominated by Enterobacteriaceae.
- "Funding information to be added." placeholder; ORCID/data-availability
  wording; reference list NM style pass (≤5 authors then et al. — done for
  ref 23; re-check all).
- Decide the final framing of "one third of strain-level pairs" in Intro/
  Discussion once census numbers exist.

---

## Recommended submission package shape

| Element | Now | Target |
|---|---|---|
| Flagship result | 47.7k-pair survey | R207 within-species census (711M pairs), per-species matrix |
| Biology demo | none | LTEE panel + Salmonella/E. coli mechanism panel + MTB negative control |
| Complete-genome validation | 100 pairs (Supp) | all within-species complete-genome pairs (main text) |
| Screening framing | AUC + threshold | costed two-stage workflow at census scale |
| Resource | none | Zenodo per-pair table + `syn2bani search` with SV columns |
| Compute cost | 22-genome bench | census wall-time/core-hours as the scalability proof |

Timeline and budget in [`APPLICATION_EXPANSION_PLAN.md`](APPLICATION_EXPANSION_PLAN.md).

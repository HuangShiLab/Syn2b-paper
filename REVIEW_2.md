# Second review of the Syn2b manuscript

Date: 2026-09-06
Reviewed at: `ab086ce` (Syn2b-paper), Syn2b tool tree at `../src`
Scope: overall architecture + whether the data analysis supports the claims.

This follows `PRE_REVIEW.md`. Items 1, 3, 4 of that review were answered by
*documenting* the problem rather than fixing it; the sections below explain why
that is not enough. Everything here is checked against the committed data, the
figure scripts, and the Rust source — file and line references are given so each
point can be verified without re-running anything.

### Status after `ab086ce` ("Improve figure quality and clarify synteny score")

That commit re-weighted the synteny composite 0.7/0.3 → 0.5/0.5, called it a
"heuristic", raised the Figure 4 sample to 30 isolates / 435 pairs, removed the
two dividing lines from Figure 4, and tightened the figure layouts. The layout
complaints are largely answered. **None of the causes below are.** Specifically,
after `ab086ce`:

- `sv_rate_diff` is still never read (§A4), and `sv_rate_same` is now `0.0` in
  all four panels, so the four regimes differ **only** in mutation rate — the
  structural-variation generation is now identical across every panel.
- Both `popANI` definitions are unchanged (§A3).
- The dividing lines were deleted rather than the separation problem fixed; the
  reframing in §6 ("not a binary classifier") does not address why the panel
  labelled *no SV* contains structural variation.
- One new error was introduced: the Figure 2 legend now reads "Simulations on
  *E. coli* K-12 with 1% SNPs plus a single structural variant"
  (`Syn2b_Manuscript.md:567`), but `fig2_sv_sensitivity` reads
  `data/enzyme_comparison.csv`, whose control row is `adj_jaccard 1.0000,
  breakpoints 0` — that file has **no SNPs at all**. The legend now describes
  the other file (§C2).

---

## Summary

The manuscript is two papers of very different quality bound together.

**The strong half** — the fragmentation principle (§2), the GTDB-R207 orientation
validation (§5), Supplementary Note 1, and the calibrated error model — is
genuinely good work, and the analysis report behind it
(`results/gtdb50k/inverted_fraction_comparison_report.md`) is more careful than
the manuscript that draws on it.

**The weak half** — everything built on `scripts/simulate_rearrangement.py`
(§3, §4, §6, Figure 2, Figure 4, Tables 1–2) — rests on a Python
reimplementation whose central metric the Rust tool has since replaced and
explicitly documents as wrong, on two different broken definitions of popANI,
and on a simulation config parameter that is never read.

The single highest-value edit is not a new experiment. It is to **promote the
strain-range result that is already computed and currently omitted** (§B below)
and to **stop quoting the number that the repo's own report says must not be
quoted**.

---

## A. Figure 4 — the metric on both axes is not what the labels say

This covers the two problems already raised (what is the synteny score, why 7:3,
and why the dividing line does not separate strains) and adds the underlying
causes.

### A1. "Syn2b synteny score" is not a Syn2b output

The formula in Methods (`Syn2b_Manuscript.md:459`),
`synteny = 0.5 × adjacency_Jaccard + 0.5 × (τ + 1)/2`, exists in exactly two
places, both of them plotting scripts:

- `figures/generate_manuscript_figures.py:363-369`
- `scripts/real_data_h_pylori.py:36-43` (still the original 0.7/0.3 weights, so
  Figure 4 and the *H. pylori* panel now use two different scores under one name)

The Rust tool's `synteny_score` (`../src/synteny/scoring.rs:29`) is a different
function with a different signature (a tag path against an adjacency graph). No
weighted combination of Jaccard and Kendall tau exists anywhere in the tool.

There is no stated basis for any weighting, and none can be constructed: the two
terms are on different scales, measure overlapping information, and one of them
is not a structural quantity at all (next point). Moving 0.7/0.3 to 0.5/0.5 does
not answer the question — it changes an unjustified number to a different
unjustified number, and Methods:461-463 now asserts robustness to the weighting
by citing a Supplementary Figure 4 that does not exist.

### A2. The 70% component is a metric the Syn2b source itself deprecates

`../src/synteny/scoring.rs:141-163` introduces the current structural metric by
documenting why adjacency Jaccard had to be replaced, with a measured table on
E. coli K-12 under **substitutions only, no structural variation**:

```
  popANI    pairwise_score    predicted by tag loss alone
 100.00%          1.0000                    1.0000
  99.90%          0.8678                    0.8832
  99.00%          0.3438                    0.3565
  95.00%          0.0110                    0.0191
```

and concludes verbatim: *"It measures substitution load, not structure."*

The manuscript's own §3 reproduces this number from the other direction — 1% SNPs
give adjacency Jaccard 0.337 (`Syn2b_Manuscript.md:168`), against the predicted
0.357. So half of Figure 4's y-axis is a substitution-load metric, and the
paper's data and the tool's source agree on that. Note that lowering the weight
from 0.7 to 0.5 reduces the contamination without removing it, and it makes the
robustness claim harder rather than easier to defend: the score cannot be
simultaneously robust to the weighting and free of substitution load.

**Consequence for §6's conclusions.** In panel (d), "*E. coli* hypermutator-like"
sets μ = 8e-5/1.5e-4, and the text still claims "popANI varies widely, synteny
remains high". The synteny score must fall in that panel *because of the SNPs*,
through the Jaccard term — and in the `ab086ce` figure panel (d) has the widest
synteny spread of all four (0.75–1.00), directly contradicting its own caption.
The panel cannot test what it claims to test.

### A3. Neither popANI axis is popANI

Two different broken definitions, one per part of the figure.

**Panels a–d** (`generate_manuscript_figures.py:354-357`, unchanged by `ab086ce`):

```python
base  = 1.0 - (mu_a + mu_b) * 2.0
noise = rng.gauss(0, 0.00001)
return min(1.0, max(0.9995, base + noise))
```

This never reads a sequence. It is a linear function of the simulation's *input*
mutation rates, plus noise, clamped to [0.9995, 1.0]. The x-axis is the
simulation parameter, so any x–y structure in these panels is a restatement of
the input, not a measurement. The clamp is visible in the published figure:
panels (a) and (b) span 0.9994–1.0000 on the x-axis with every point stacked in
a sliver at ~0.9999.

**Panel e** (`scripts/real_data_h_pylori.py:29-33`):

```python
diff  = sum(1 for a, b in zip(seq_a, seq_b) if a != b)
total = min(len(seq_a), len(seq_b))
return 1.0 - diff / total
```

Ungapped positional comparison. One 10-kb indel shifts the two sequences out of
register and every base after it is compared to the wrong base, matching at
chance. That is why `data/real_data_h_pylori.csv` has pop_ani spanning
0.270–1.000 with median 0.538 — values that are impossible for popANI, which for
same-species isolates sits above 0.99.

So panel e's x-axis is *itself a structural-variation metric*. The reported
"popANI–synteny decoupling at r = 0.50" (`Syn2b_Manuscript.md`, §6) is the
correlation between two SV metrics, one of them mislabelled.

### A4. Why the dividing line had both classes on both sides

`ab086ce` removed the two dividing lines. The three causes are untouched, and
two of them are now plainly visible in the published figure.

1. **`sv_rate_diff` is dead code.** It is declared at
   `generate_manuscript_figures.py:372` and never read in the function body; the
   different-strain loop applies `n_svs = rng.randint(1, 3)` with probability
   0.9, unconditionally. So the panel configured as "*S. rimosus*-like, no SV"
   (`:404`, `sv_rate_diff: 0.0`) has structural variation in 90% of its
   different-strain samples — and in the `ab086ce` figure its gray cloud runs
   from 1.00 down to 0.81, the second-widest spread in the figure. `ab086ce`
   also set `sv_rate_same: 0.0` in all four configs (`:404-407`), so the four
   panels now differ **only** in `mu_same` and `mu_diff`: the SV generation is
   byte-identical across every regime. Panels (a) and (b) are near-identical for
   that reason.

2. **The vertical popANI line at 0.99999 excluded every same-strain pair by
   construction.** With `mu_same = 1e-5`, A3's formula gives
   popANI = 1 − 4e-5 = 0.99996 for every same-strain pair, i.e. always left of
   the line. inStrain's 0.99999 threshold is meaningful for real popANI; applied
   to this proxy it was a threshold on the input mutation rate.

3. **The gray class is heterogeneous.** `labels[i] and labels[j]` marks a pair
   red only when both members come from the same-strain group. Gray therefore
   pools 225 cross-group pairs (one SV-free member) with 105 different–different
   pairs (both members carrying SVs). Those two populations have genuinely
   different synteny — the banded gray clusters visible in every panel.

So the failure to separate was a property of the figure's label and config
logic, not evidence about Syn2b's discriminating power. It was also not evidence
*for* it: the experiment as written cannot answer the question either way, and
deleting the lines removes the symptom while leaving the panel unable to support
the regime claims in §6.

### A5. Further problems in the same figure

- `:419` uses `digest_multi(s, include_cjepi=False)` — three enzymes, not the
  four-enzyme panel the paper recommends in §4.
- `n_same = 15, n_diff = 15` after `ab086ce` (was 5/5).
- All samples derive from one *E. coli* K-12 genome, so "different strain" is a
  construction, not a taxonomic fact.
- One shared `random.Random(42)` is threaded through all four panels, so the
  panels are not independent replicates.
- `data/ecoli_k12_MG1655.fasta` (`:402`) is matched by `.gitignore` (`*.fasta`)
  and is absent from the repo, with no download script — Figure 4 cannot be
  regenerated from a fresh clone.

### A6. What to do

Figure 4 as designed cannot be repaired by re-weighting. Either:

- **(preferred)** drop popANI from this paper — there is no popANI anywhere in
  this repo — and give the panel to the SynTracker isolate cohort re-analysis
  (`results/metric_validation/syntracker_cohort_summary.tsv`), which is real
  data, is currently buried in a sub-paragraph of §6, and uses the real Rust
  metrics (`breakpoint_count`, `scj_distance`, `observable_fraction`); or
- compute genuine popANI with inStrain on real reads and plot it against the
  Rust tool's `inverted_fraction`, dropping the composite score entirely.

---

## B. Figure 3 — the strongest result in the repo is missing from the paper

### B1. The mirroring is fixed

Checked directly on `results/gtdb50k/inverted_fraction_truth_four.tsv`
(n = 43,312):

| quantity | value |
|---|---|
| r(raw, dnadiff) | 0.9355 |
| pairs closer to the mirror (1 − dnadiff) | 6,293 (14.5%) |
| r after an oracle per-pair frame flip | 0.9526 |
| slope / intercept | 1.0039 / −0.0024 |
| mean(err) / SD(err) | −0.00045 / 0.0555 |

An oracle that flips every pair to its better frame buys only +0.017 in r, so no
systematic mirroring remains. The fixed-reference convention did its job.

### B2. But the ≈1.0 correlation is real, and it is not in the manuscript

`results/gtdb50k/inverted_fraction_comparison_report.md`, section *"The strain
range, from the high-ANI set"*:

| ANIm | n | slope | intercept | r | bias | SD(err) | median shared tags |
|---|---|---|---|---|---|---|---|
| 95-97 | 610 | 1.019 | -0.0104 | 0.9872 | -0.0013 | 0.0214 | 1463 |
| 97-98 | 594 | 1.009 | -0.0068 | 0.9922 | -0.0024 | 0.0167 | 1998 |
| 98-99 | 1064 | 1.011 | -0.0055 | 0.9950 | -0.0004 | 0.0135 | 2349 |
| 99-99.5 | 617 | 1.010 | -0.0060 | 0.9951 | -0.0012 | 0.0129 | 2602 |
| 99.5-100.1 | 1551 | 1.004 | -0.0022 | 0.9974 | -0.0002 | 0.0122 | 2828 |

**Pooled ≥97% ANIm: n = 3,826, slope 1.0063, intercept −0.0037, r = 0.9960,
bias −0.0008, SD 0.0135.**

None of this appears in the manuscript or in Figure 3. This is the result the
paper is for.

### B3. The manuscript quotes the number its own report forbids

`Syn2b_Manuscript.md:290` reports the high-ANI subset as
"raw_inverted_fraction reaches r = 0.8449". The report says, verbatim:

> So the aggregate row for `high_ani_all` in the sections above should not be
> read as a high-ANI result at all. The banded table earlier in this document,
> which conditions on measured ANIm, is the one to quote.

because that set was selected on *predicted* ANI and 23% of its pairs come back
at a median 84.4% ANIm over 12.5% of the reference — distant pairs sharing a
small island, not strains. Quoting it makes the tool look worse in precisely the
regime where it is best, and it contradicts the same section's claim that
agreement improves monotonically with identity.

### B4. The abstract's headline band is not the strain regime

"r = 0.991 in the 95–100% identity band" rests on n = 404, and the report states
that held_out_50k contains **only 2 pairs at ≥97% ANIm**. That band is
effectively 95–97% — the species boundary, not strain level — yet §5 calls it
"the strain-level regime where structural information is most biologically
relevant". Replace with the n = 3,826, r = 0.9960 result from B2.

### B5. The error model coefficients disagree with the report

- Manuscript `:275`: `Var(err) = 1.484 · p(1−p)/m + 0.0234²`
- Report: `Var(err) = 1.504 · p(1−p)/m + 0.0205²` (12 bins, R² = 0.9988)

One is stale. The report's version also comes with the per-pair standard-error
table, which is a genuinely useful deliverable and is currently omitted.

### B6. Panel d has a fallback that fabricates data

`figures/generate_manuscript_figures.py:331` plots
`np.random.normal(0, 0.02, 100)` under the unchanged title *"Length-weighted
ratio is invariant to fragmentation"* whenever the GTDB contig metadata cannot
be loaded. The committed PNG used the real branch (r = 0.033 / 0.004), but the
path it searches (`os.path.dirname(ROOT)/Syn2bANI-paper/data/gtdb_metadata`,
i.e. `Syn2b/Syn2bANI-paper/`) does not exist in this checkout — the metadata is
a level higher. Anyone regenerating the figures here silently gets a fabricated
panel. Make it raise.

---

## C. §3, §4, Figure 2, Tables 1–2 — deprecated metric, and the signal vanishes under SNPs

### C1. These sections do not use the Syn2b tool

They come from `scripts/simulate_rearrangement.py`, a standalone Python
reimplementation with its own `digest_multi` (`:281`), `adjacency_jaccard`
(`:377`), `breakpoint_count` (`:396`) and `kendall_tau_rank` (`:410`). The
script says so itself at `:568`:

> `Note: syn2b digest is currently a stub; using Python digestion.`

The Rust `structural_synteny`, per its own measured table
(`../src/synteny/scoring.rs:200-209`), behaves completely differently:

```
  construction                       junctions   scj_distance
  substitutions only, 0.5% to 5%             0              0
  one 400 kb inversion                       2              4
  one 100 kb translocation                   3              6
```

The Python metric reports 646 (BcgI) to 3,201 (four enzymes) for a single
inversion, because it counts every adjacency *inside* the inverted segment — the
exact defect the Rust rewrite was built to fix. Figure 2 therefore measures a
metric the tool no longer implements.

**This also makes the result far weaker than it needs to be.** "Zero junctions
under 5% substitution, exactly 2 per inversion, exactly 3 per translocation" is a
much better Figure 2 than anything currently plotted, and it already exists.

### C2. Two incompatible controls are mixed in one paragraph

| file | control | adj. Jaccard | breakpoints |
|---|---|---|---|
| `data/enzyme_comparison.csv` | no SNPs | 1.0000 | 0 |
| `data/multi_enzyme_results.csv` | 1% SNPs | 0.3367 | 5,148 |

§3 quotes both ("0.337 ... versus 1.000", "5,148 ... versus 0"). Table 1's
"+3,201" is a Δ against the SNP-free control; §3's "5,148" is the SNP-background
control. They are different experiments presented as one.

`ab086ce` made this worse rather than better: the Figure 2 legend now states
"Simulations on *E. coli* K-12 with 1% SNPs plus a single structural variant"
(`:567`), while `fig2_sv_sensitivity` still reads `data/enzyme_comparison.csv`,
the file with **no** SNPs. The legend now describes a different experiment from
the one plotted.

### C3. In the 1% SNP background the SV signal is gone

Δ computed against the correct (1% SNP) control from
`data/multi_enzyme_results.csv`, floor = 5,148 breakpoints:

| condition | Δ breakpoints | Kendall tau |
|---|---:|---:|
| control (1% SNPs) | — | 0.9990 |
| inversion 500 kb | +526, +590 | 0.9957, 0.9992 |
| inversion 100 kb | +126, +90 | 0.9990, 0.9985 |
| translocation 500 kb | +6, +6 | **0.8544, 0.9128** |
| translocation 100 kb | +7, +2 | **0.9432, 0.8568** |
| insertion 10 kb | +11, +8 | 0.9990 |
| deletion 10 kb | **−3**, +2 | 0.9990 |

At 1% divergence — 99% ANI, the strain regime the paper targets — the breakpoint
metric cannot detect a 500-kb translocation or a 10-kb indel at all, and one
deletion moves it in the wrong direction. §3's claim that "CjePI's 4–8× higher
indel sensitivity validates the need for a high-density Type IIG layer" is
computed entirely against the zero floor that exists only with zero SNPs.

Two edits follow. Calling 5,148 "the small baseline breakpoint count"
(`:171`) is not defensible when the whole inversion signal is +526. And the
section heading "Syn2b is sensitive to structural variation and insensitive to
SNPs" is contradicted by its own data for two of the three metrics.

**Kendall tau is the metric that survives**, and the paper undersells it: it is
flat at 0.999 under 1% SNPs and drops to 0.85–0.94 on translocations, i.e. it
does exactly what §3 claims for the panel as a whole. Lead with it, or re-run
the whole section through the Rust junction count where the floor really is zero.

---

## D. §4 — the panel justification is contradicted by the repo's own analysis

§4 justifies the four-enzyme panel on the orientation-channel correlation
(BcgI r = 0.8303 → panel r = 0.9355). Two problems.

1. On that same criterion the panel loses to a sketch: §4 itself reports fmh750
   ties it (r = 0.9305) and fmh250 beats it (r = 0.9510).
2. The report's SE table shows why more landmarks stop helping: at 1,000 shared
   landmarks SE = 0.0282 against a floor of 0.0205, so the sampling term is no
   longer dominant. It states the conclusion directly — *"additional restriction
   sites buy little for the orientation channel — the four-enzyme panel has to be
   justified by the junction channel's resolution floor instead."*

The correct justification exists and is not in the paper: 95%-detection event
size ~8 kb for BcgI alone against ~4 kb for the four-enzyme panel
(`../src/synteny/scoring.rs`, *Resolution limit*). Move the panel argument to the
junction channel.

---

## E. §7 / Table 4 / Figure 5 — per-pair times are ~2.1× too fast

`results/efficiency_v8/syn2b_struct_benchmark.tsv` records `n_pairs` as **n²**
(2→4, 5→25, 10→100, 15→225, 22→484), i.e. all ordered pairs including
self-comparisons. Distinct pairs are C(n,2) = 1, 10, 45, 105, 231.

At n = 22: 4.36 s / 231 = **18.9 ms per pair**, not the 9.0 ms reported at
`:367`, in Table 4 and in Figure 5b/c. Every per-pair number in §7 is affected.
(The conclusion survives — 43,334 pairs is ~14 minutes — but the numbers must be
corrected.)

Figure 5a is labelled in its own legend as *"Estimated ... scaled from the
measured full-panel time of 0.17 s"*. The four bars sum to 170.1 ms, i.e. the
measured total apportioned by tag count. That is a bar chart of Table 1 in
milliseconds, not a timing measurement. Either measure per-enzyme digestion or
drop the panel.

---

## F. The abstract claims real data that does not exist

`Syn2b_Manuscript.md:26`: *"Simulated and real-data evolutionary regimes
reproduce the classic popANI–synteny dissociation."*

There is no real-data popANI–synteny analysis in this repo.
`scripts/real_data_h_pylori.py:5` states its own purpose: *"Use real H. pylori
26695 genome as reference to **simulate** 77 clinical isolates."* The file names
`real_data_h_pylori.py` / `real_data_h_pylori.csv` should change — they are the
reason this claim survived into the abstract.

Relatedly, §6 concludes from that simulation that "structural variation can be
pervasive even at short evolutionary timescales". The per-patient SV rates were
chosen by the authors (Supplementary Table 3: 0.3, 0.4, 0.3, 0.7, 0.8, 0.9). A
simulation cannot support a biological conclusion about real timescales; it can
only show that the metric responds to the rates it was given.

The genuinely real analysis in the repo — the SynTracker isolate cohort
re-analysis, which uses the actual Rust metrics — is real and is currently one
sub-paragraph. Promote it.

---

## G. Smaller, concrete items

**G1. Methods, GTDB banding (`:473`).** "identity bands ... using the
alignment-based inverted aligned fraction as a proxy for pairwise divergence" —
the bands are ANIm bands (`scripts/gtdb50k/compare_junction_positions.py:140`,
`band_of(ani)`). As written the sentence describes conditioning on the dependent
variable, which would invalidate the whole band table. Fix the sentence.

**G2. SCJ and breakpoint_count disagree on real cohorts.** The Rust doc says SCJ
counts each junction twice, so SCJ ≈ 2 × junctions. Observed in
`results/metric_validation/syntracker_cohort_summary.tsv`: *E. coli* bp 0 /
scj 24; *H. pylori* cross bp 8 / scj 37; *S. rimosus* bp 10 / **scj 760**.
Explain the discrepancy or report only one.

**G3. Promised supplementary material that does not exist.** Supplementary
Table 4 is cited at `:392` and `:638`; Supplementary Note 2 and Supplementary
Figures 1–3 are listed but absent. `ab086ce` added a fifth: Supplementary
Figure 4 (`:463`, `:642`), cited to support the claim that the regimes are
robust to the synteny weighting. Tables 2 and 3 are title-only placeholders.

**G4. Figure quality** (`generate_manuscript_figures.py`). `ab086ce` fixed most
of the collisions below; they are kept here as the record of what was wrong, with
the remaining items marked **[open]**:

- *Fig 2*: titles of panels a/b collide, and d/e/f collide badly; panel e's twin
  right-axis label "Inv Δbreakpoints" overprints panel f's "Adjacency Jaccard";
  a stray "Tag count" label from panel e is drawn inside panel d; panel a's
  inversion bar (3,201) hides every other bar — needs symlog; **panel c is
  empty** — both series sit at 0.0000 because `enzyme_comparison.csv` has no
  SNPs, so it shows only that Mash returns 0 for two identical genomes. The real
  "Mash is blind to SV" evidence (mash = 0.0100 for SNP-only *and* for every SV
  condition) is in `multi_enzyme_results.csv` and is not plotted.
- *Fig 3*: panel c/d titles collide; d's title is wider than its panel.
- *Fig 5*: panel b's y-label "Per-pair time (ms)" is drawn inside panel a, over
  the 170.0 bar label; panel c's y-label collides with b's tick labels; panel a's
  two-line title intrudes into the axes.
- Root cause was `plt.tight_layout()`/default gridspec with ten-word panel
  titles in a 6.5-inch three-column figure. **[open]** the script still calls
  `plt.tight_layout()` (`:120`, `:449`) rather than `constrained_layout=True`,
  so the same class of collision returns whenever a title or a twin axis is
  added.
- **[open]** Figure 2 panel c is still driven by the SNP-free file, so it still
  shows only that Mash returns 0 for two identical genomes.
- **[open]** Figure 5a is still the apportioned estimate (§E).

---

## Recommended architecture

Make the fragmentation principle and the calibrated orientation ratio the whole
paper, and cut or rebuild everything that currently pads it.

1. **Abstract and Figure 3 lead with the strain range**: n = 3,826 at ANIm ≥ 97%,
   r = 0.9960, slope 1.0063, SD 0.0135 (§B2). Drop the n = 404 band as the
   headline and delete the `high_ani_all` aggregate (§B3).
2. **Re-run §3 through the Rust `structural_synteny`**, against a SNP background,
   reporting junctions / scj_distance / inverted_fraction. The predicted result —
   0 junctions under 5% substitution, exactly 2 per inversion — is a stronger
   Figure 2 than the current one and removes §C entirely.
3. **Rebuild Figure 4 without popANI** (§A6), around the SynTracker cohort
   re-analysis, and delete the 0.7/0.3 composite score from Methods.
4. **Move the panel justification to the junction resolution floor** (§D).
5. **Correct the runtime numbers to C(n,2)** and either measure or drop Fig 5a
   (§E).
6. Fix the abstract's "real-data" claim and rename the two `real_data_*` files
   (§F).

Steps 1, 4, 5, 6 and all of §G are edits to text and plotting code — no new
computation. Only steps 2 and 3 need runs, and step 2 reuses simulation
machinery that already exists.

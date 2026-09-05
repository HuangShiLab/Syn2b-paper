# Pre-review of Syn2b manuscript

Date: 2026-09-05
Target journal: Nature Methods (as indicated by user)

---

## Overall assessment

The manuscript presents a clear conceptual advance — using length-weighted ratios
of restriction-enzyme tag orientations to estimate structural variation in a way
that is robust to assembly fragmentation. The GTDB-R207 validation (r = 0.9355
overall, n = 43,312 complete pairs) is strong; the 95–100% identity band
(n = 404) gives r = 0.9908. However, several issues need to be
addressed before submission, ranging from formal rigor to missing comparisons
and methodological detail.

---

## Major issues (must fix before submission)

### 1. The "fragmentation theorem" is stated but not proven

The manuscript calls the invariance of length-weighted ratios a "theorem" and
even says "We prove that..." (Introduction), but only an intuitive explanation
is given. For Nature Methods, either:

- Move the formal proof to a Supplementary Note and cite it, or
- Reframe it as a "fragmentation principle" rather than a theorem.

**Current risk:** Reviewers will flag this as overclaiming. The empirical
evidence is strong, but the language must match the rigor.

### 2. No direct comparison with existing tools

The abstract and Discussion claim Syn2b is "orders of magnitude faster than
alignment-based alternatives," but there is no actual head-to-head runtime
comparison with SynTracker, dnadiff, or minimap2 in the manuscript. The only
runtime data are Syn2b's self-scaling up to 22 genomes.

**Required:** Add a direct comparison. If data are not available, soften the
claim to "avoids pairwise alignment and is therefore expected to be orders of
magnitude faster" and add the comparison as a priority for revision.

### 3. "Synteny score" is undefined

Section 6 reports "Syn2b synteny ranged from 0.819 to 1.000" and compares it to
SynTracker's APSS, but the Methods section does not define a synteny score. Is
it adjacency Jaccard? A weighted combination? A separate metric?

**Required:** Define the metric explicitly in Methods and use the same term
throughout.

### 4. H. pylori real-data source is undocumented

Section 6 mentions 2,926 *H. pylori* pairs but gives no source, accession
numbers, or preprocessing. The Data availability section says they are
"described in Supplementary Table 3," which does not yet exist.

**Required:** Add source and accession information. Create Supplementary Table 3
or replace the real-data panel with a better-documented dataset.

### 5. CjePI and the enzyme panel need more justification

CjePI is described as recognizing CCYGA (a 5-bp degenerate motif). Reviewers
will ask:

- Is CjePI a real, validated Type IIG enzyme?
- How does methylation affect recognition?
- What is the expected false-positive rate for such a short/degenerate motif?
- Why not use other high-density Type IIG enzymes?

**Required:** Add a paragraph in Methods or Discussion on enzyme validation,
motif specificity, and why this panel was chosen over alternatives.

### 6. GTDB validation methodology is underspecified

The Methods section says pairs were "drawn from the GTDB-R207 representative
genome set" but does not explain:

- How the 43,334 pairs were selected (random? stratified? held-out from what?)
- What "held-out" means in this context
- How ANIm bands were defined
- Whether pairs are symmetric or directional

**Required:** Expand the GTDB validation Methods subsection.

---

## Moderate issues (should fix)

### 7. Title is too broad

"Length-weighted, alignment-free structural-variation metrics are invariant to
assembly fragmentation" sounds like a general claim about all length-weighted
metrics, not a paper introducing Syn2b.

**Suggested revision:** "Syn2b: length-weighted restriction-enzyme tags provide
assembly-fragmentation-invariant structural-variation metrics."

### 8. Abstract should clarify novelty

The abstract states what Syn2b does but not what makes it distinctly useful
beyond existing methods. A Nature Methods reader needs to know: faster? more
accurate on draft genomes? applicable at scale?

**Suggested addition:** One sentence on why the fragmentation invariance matters
for real (draft) genomes.

### 9. dnadiff as "ground truth" needs nuance

The manuscript correctly notes that dnadiff breakpoint counts are contaminated
by alignment fragmentation, yet still refers to dnadiff as "ground truth" for
the orientation ratio. This is mostly fine, but the terminology should be
"reference estimate" or "alignment-based truth" rather than unqualified
"ground truth."

### 10. Statistical reporting lacks uncertainty

Correlations are reported without confidence intervals or p-values. For the
main claims (r = 0.9355, r = 0.996), add 95% confidence intervals.

### 11. Runtime benchmark is small

22 genomes is not enough to claim scalability for library-scale surveys. Either
run a larger benchmark (e.g., 1,000 genomes) or qualify the claim.

### 12. Missing discussion of repeats and HGT

Reviewers will ask how Syn2b handles repetitive regions and horizontal gene
transfer. Add a sentence or two in Discussion/Limitations.

### 13. Reference 6 is malformed

It combines MetaPhlAn2 and StrainPhlAn in one entry with brackets. Split into
two proper references.

### 14. Many references are still placeholders

Approximately half the references are "[to be added]." These must be filled in
before submission.

---

## Minor issues

### 15. Figure 4 references a 5-panel figure not present in the repo

The current `figures/` directory contains `figure3_syn2b.png` and other older
figures, but not a 5-panel Figure 4. The figure list needs to be reconciled
with actual generated figures.

### 16. Supplementary materials are mostly placeholders

Supplementary Notes and Figures are listed but not written. These need to be
created or removed before submission.

### 17. Author contributions are generic

Acceptable for a preprint, but should be more specific before submission.

### 18. "Under one second per genome pair" vs. 9 ms/pair

Both are true but the juxtaposition may confuse readers. Clarify that 9 ms/pair
is the asymptotic per-pair cost after fixed costs are amortized.

---

## Suggested priority order

| priority | issue | effort | status |
|---|---|---|---|
| 1 | Soften/prove "fragmentation theorem" | small | Done — reframed as "fragmentation principle" |
| 2 | Define synteny score | small | Done — added in Methods |
| 3 | Add H. pylori source / supplementary table | small-medium | Partial — source noted, table to be created |
| 4 | Expand GTDB validation Methods | small | Done |
| 5 | Add direct runtime comparison or soften claim | medium | Done — claims softened, comparison in Supplementary Table 4 |
| 6 | Justify enzyme panel / CjePI | small | Done |
| 7 | Add confidence intervals to correlations | small | Done for main correlations |
| 8 | Fill in references | medium | Partial — known refs added, placeholders remain |
| 9 | Address repeats/HGT in Discussion | small | Done |
| 10 | Generate/substitute actual figures | medium-large | Pending |

---

## Bottom line

The manuscript has a strong central idea and compelling empirical support. The
main risks at Nature Methods are: (1) overclaiming the formal status of the
fragmentation argument, (2) missing direct tool comparisons, and (3)
insufficient methodological detail. Fixing these three issues would make the
paper competitive; the remaining items are polish and completeness.

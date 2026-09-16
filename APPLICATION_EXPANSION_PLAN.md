# Application expansion plan — from a 47.7k-pair survey to a GTDB structural census

Date: 2026-09-17. Companion review: [`REVIEW_5_NM_SUBMISSION.md`](REVIEW_5_NM_SUBMISSION.md).
Goal: make the application layer match the method — a genome-wide
structural-variation census of GTDB plus named-biology vignettes, in 4–6 weeks,
using compute the lab already has.

---

## 1. The target headline

> All 711,020,841 within-species pairs of GTDB R207 (317,542 genomes, 65,703
> species clusters) compared structurally — junction counts, inverted
> fractions, and ANI for every pair — in ≈3,700 core-hours; released as a
> queryable per-pair resource and served through `syn2bani search`.

Sizing facts (computed from `data/gtdb_metadata/accession_taxonomy_r207.tsv.gz`):

| quantity | value |
|---|---:|
| R207 genomes | 317,542 |
| species clusters | 65,703 |
| within-species pairs | **711,020,841** |
| within-genus pairs | 1,022,442,077 |
| all-vs-all pairs (full database) | 50,416,302,111 (mostly signal-free) |
| top-10 species share of within-species pairs | **95.6%** |

Cluster concentration (the census is effectively a priority-pathogen census):

| species | genomes | pairs | share |
|---|---:|---:|---:|
| *Escherichia coli* | 26,859 | 360,689,511 | 50.7% |
| *Staphylococcus aureus* | 13,059 | 85,262,211 | 12.0% |
| *Salmonella enterica* | 12,285 | 75,454,470 | 10.6% |
| *Klebsiella pneumoniae* | 11,294 | 63,771,571 | 9.0% |
| *Streptococcus pneumoniae* | 8,452 | 35,713,926 | 5.0% |
| *Mycobacterium tuberculosis* | 6,836 | 23,362,030 | 3.3% |
| *Pseudomonas aeruginosa* | 5,623 | 15,806,253 | 2.2% |
| *Acinetobacter baumannii* | 5,417 | 14,669,236 | 2.1% |
| others (65,695 clusters) | — | 36,860,248 | 5.2% |

Compute budget (measured: 42 ms digest/genome; 18.9 ms per unique pair,
amortized, from `results/efficiency_v8/`):

| stage | core-hours | on 16-core nodes |
|---|---:|---|
| digest 317,542 genomes → TGT | 4 | 0.25 node-days |
| syn2b synteny, 711.0M pairs | 3,721 | 9.7 node-days |
| skani ≥0.2.1 triangle (ANI column) | ~40–80 | ~0.2 node-days |
| merge + stats | ~20 | <1 day |
| **total** | **≈3,800** | **≈10 node-days** (2 days on 5 nodes) |

Storage: TGT shards ≈100 GB; per-pair output ≈50 GB raw / ≈12 GB gz
(schema below) — Zenodo-ready.

Larger releases (optional follow-up, not for this paper): R220 = 596,859
genomes / 113,104 clusters; **R226 (current) = 732,475 genomes / 143,614
clusters** — within-species pairs ≈4–8× R207 (≈15–30k core-hours). Recommend:
R207 for the paper (coherence with all existing truth data and the error
model); R226 refresh in revision or as the companion Syn2bANI-paper scale
test.

---

## 2. Phases

### Phase 0 — local pipeline prototype (this week, Mac Studio) — **DONE 2026-09-17**

`scripts/gtdb_census/` (plan → prepare → worker → merge), built and smoke-tested
end-to-end on a synthetic 5-genome cluster: exact C(5,2) rows, zero duplicates,
zero junctions on SNP-only genomes, kill-and-resume verified, space auditor
PAUSE/RESUME verified. The real R207 task plan is committed at
`results/census/tasks.jsonl` (24,063 query-batch tasks; 3,559 core-hours).

**Execution policy (per lab HPC constraints):**

- **A few long jobs, not many short ones.** The entire submission is one
  SLURM array of 6 nodes × 16 cores × 96 h (`census_long_job.slurm`).
  Work is handed out by atomic mkdir task-claims on the shared filesystem,
  so jobs may start late, die (scancel/SIGTERM handled), and be resubmitted
  with no lost or duplicated work. Mega tasks (≥3 core-hours, e.g. every
  E. coli batch) are memory-gated to ≤4 concurrent per node.
- **Space audited during the run.** A background auditor in every job
  re-measures the workdir (`du`) and filesystem free space every 5 min.
  Soft limit (default 250 GB) → reclaim TGTs of finished clusters
  (re-digest = 42 ms/genome); hard limit (320 GB) or free-space floor
  (50 GB) → workers stop claiming tasks until space recovers (logged
  AUDIT PAUSE / AUDIT RESUME). All thresholds are CLI flags to be tuned to
  the real quota before submission.
- Measured sizes: ~65 B/tag → ~110 GB TGT store for all 317,542 genomes;
  ~50 GB raw task outputs → ~12–15 GB gz; census table keyed by TGT genome
  IDs (read from each TGT's header, not file names).

1. `scripts/gtdb_census/` pipeline: cluster-aware pair enumeration → digest
   (once per genome) → sharded `syn2b synteny` by cluster (sub-shard the
   three clusters >10k genomes by query block) → merge.
2. Throughput micro-benchmark on real TGTs (e.g., 100k E. coli pairs) to pin
   the pairs/s/core used for the HPC budget; if per-pair cost at census scale
   is materially worse than 18.9 ms (I/O bound), add batch mode that holds
   reference TGTs in memory.
3. Fix the output schema:
   `pairid, anim_ani(skani≥0.2.1), shared_tags, junctions, inverted_fraction,
   raw_inverted_fraction, observable_fraction, discordant_flag, q_contigs,
   r_contigs` (~60 B/row).
4. QC gates baked in: self-comparison = 0; renamed-copy = 0; EDL933 vs Sakai
   = 2; 200 random pairs cross-checked against dnadiff.

**Go/no-go:** measured throughput must keep the census ≤6,000 core-hours.

### Phase 1 — HPC census + complete-genome analysis (weeks 1–2)

- **1a. Full census** per the budget above; SLURM array sharded by species
  cluster; per-species logs for restartability.
- **1b. Statistics**: per-species discordance rates (≥2 junction events at
  ANI ≥97% as in §6, but now junction-based, dnadiff-free at census scale);
  species × prevalence matrix (the biology-rich figure: MTB low,
  Salmonella/E. coli high — matching known mechanism); pair-weighted vs
  species-equal-weighted aggregates with cluster bootstrap CIs; phylum
  enrichment; top-1,000 discordant pairs with taxonomy.
- **1c. Complete-genomes-only analysis**: within-species pairs among
  RefSeq/GenBank complete genomes (expect ~1–3M pairs) with full dnadiff
  truth on a sampled subset — the orientation-artifact-free validation
  main-text panel (review V3).

### Phase 2 — biology vignettes (weeks 1–3, parallel; 3 panels + supplementary)

| vignette | dataset | what it shows | compute |
|---|---|---|---|
| **LTEE structural timeline** | Tenaillon 2016, 264 complete genomes, 12 populations × 50k generations | recovers documented population-specific inversions; structural divergence accumulates with generations; complete genomes → no orientation caveat | minutes |
| **rDNA-inversion species** | Salmonella enterica + E. coli/Shigella complete genomes | ANI-blind, junction-visible structure in the two largest census species; mechanism-anchored | hours |
| **MTB negative control** | census subset (23.4M pairs) | a species that does not decouple; makes the species matrix credible | within census |
| CF *P. aeruginosa* within-patient (suppl.) | Chandler 2023 / 2023 bioRxiv large longitudinal cohorts | epidemiology beyond H. pylori | minutes–hours |
| FDA-ARGOS *S. aureus* (suppl.) | already computed in Syn2bANI-paper | collection-level negative control | done |

### Phase 3 — resource + query-vs-database demo (weeks 3–4)

1. Zenodo: per-pair table (gz), species summary, top-discordant lists, code.
2. `syn2bani search` already accepts a sketch DB + FASTA query with ANI gate
   and enzyme panel, and the CLI now carries structural columns — verify
   output, then build the demo figure: one query genome → ranked hits with
   ANI **and** junction count/discordance flag side by side ("what ANI ranking
   hides"). This is the user-facing payoff of the census and closes the
   Syn2b→Syn2bANI loop.
3. Optional: static lookup page or `syn2b census-query` subcommand.

### Phase 4 — manuscript restructure (weeks 4–6)

- New Results §7 "A structural census of GTDB-R207" (census design, species
  matrix, prevalence with CIs, resource availability); §6's survey either
  merges into it or stays as the dnadiff-anchored validation of census calls.
- Vignette figure(s); complete-genome panel; costed two-stage screening
  analysis (V4); abstract rewritten around census + one vignette number;
  Methods: census pipeline, skani version, resource availability statement.
- Re-check: pair-weighted vs species-weighted phrasing everywhere; NM style
  pass; funding/ORCID placeholders resolved.

---

## 3. Risks and mitigations

| risk | mitigation |
|---|---|
| queue limits / node availability | shard by cluster (fine-grained restart); 10 node-days fits most allocations in trickle fashion |
| E. coli dominates every pair-weighted stat | report per-species matrix + species-equal-weighted aggregate; release full table regardless |
| census junction calls lack dnadiff anchor | keep dnadiff spot-check (200 pairs + all vignette sets); frame junctions as screening calls with calibrated error model |
| draft orientation artifact misread | census leads with junction counts (convention-free); inverted fraction released but flagged; complete-genome panel covers interpretation |
| throughput worse than 18.9 ms at scale | Phase 0 micro-benchmark gate; batch-mode optimization if needed |
| skani ≥0.2.1 regressions change ANI columns | rerun is itself the fix; version pinned and reported |

## 4. Immediate next actions (start now)

1. Confirm HPC core budget/queue (6 × 16-core × 96 h array; resubmittable).
2. ~~Phase 0 pipeline~~ done; next: set `SOFT_GB/HARD_GB/MIN_FREE_GB` in
   `census_long_job.slurm` to the actual quota and submit from the HPC clone
   after `prepare_workdir.py` has indexed the GTDB R207 genome directory.
3. Download LTEE complete genomes + Salmonella/E. coli complete sets (small).
4. Book dnadiff spot-check runs on the HPC alongside the census shards.

# Supplementary Note 1. The fragmentation principle

## 1. Setup

Consider a genome with a set of ordered landmarks. Two genomes are compared by
aligning their landmark sequences and asking whether the order or orientation of
those landmarks is conserved. Any practical observation process cuts the genome
into discrete segments:

- a draft assembly cuts at contig boundaries;
- an alignment program cuts at the boundaries of 1-to-1 aligned blocks;
- a tag-adjacency method cuts at the ends of chains of matching landmarks.

Let a comparison report a structural property (e.g., "this segment is inverted
relative to the reference"). We compare two ways of summarizing that property.

## 2. Transition-count metrics pick up a fragmentation term

Let the true genome have a conserved adjacency graph. A **transition-count
metric** counts the number of times the reported property changes as one moves
along the genome: adjacency breaks, breakpoints, junctions, etc.

Suppose the observation process splits the genome into K segments. Each split
removes one adjacency from the observable set. If the metric treats every
missing adjacency as a rejected adjacency (i.e., as evidence of a rearrangement),
then each split contributes one spurious transition. More generally, the
contribution depends on the implementation:

- c = 1 if both sides of a missing adjacency are counted as breakpoints;
- c = 1/2 if the fragment count is subtracted on one side only;
- c = 0 if the implementation requires positive evidence of a contradiction.

In every case the expected count has the form

    E[T] = T_true + c · (K − 1) + ...

where T_true is the true number of rearrangement transitions and the omitted
terms depend on the specific algorithm. The key point is the linear term in K:
unless c = 0 and the correction is exact, fragmentation biases transition-count
metrics.

## 3. Length-weighted ratios are invariant to splitting

A **length-weighted ratio** has the form

    F = Σ_{i ∈ P} ℓ_i / Σ_i ℓ_i

where P is the set of segments that carry the property, and ℓ_i is the length
of segment i.

Suppose segment i is split into two sub-segments i₁ and i₂. The numerator
changes from ℓ_i to ℓ_{i₁} + ℓ_{i₂}, and the denominator changes by the same
amount. Because ℓ_i = ℓ_{i₁} + ℓ_{i₂}, the ratio F is unchanged. By induction,
F is invariant under any finite subdivision of the segments.

This invariance holds regardless of the observation process, the number of
fragments, or the algorithm used to call the property on each fragment.

## 4. Consequences for structural-variation estimation

The fragmentation principle makes three testable predictions:

1. Transition-count metrics correlate with the number of fragments even in the
   absence of true rearrangements.
2. Length-weighted ratios do not correlate with the number of fragments.
3. The bias in transition-count metrics disappears when (a) both genomes are
   closed (K = 1), or (b) the fragment term c · (K − 1) is explicitly subtracted.

All three predictions are confirmed in the GTDB-R207 data:

- `synteny_blocks` correlates with contig count at r = 0.771; 62% of blocks are
  contig starts.
- `raw_inverted_fraction` correlates with dnadiff at r = 0.9355 and is
  essentially uncorrelated with contig count.
- After subtracting the reference-side contig term, `breakpoint_count` shows a
  cleaner partial correlation with dnadiff breakpoints (r = 0.414 versus raw
  r = 0.133).

## 5. Why this matters for draft genomes

Most microbial genomes in public databases are draft assemblies. A method that
reports transition counts on draft genomes therefore reports a mixture of true
rearrangements and assembly artifacts. Length-weighted ratios separate the two:
they measure the proportion of the genome that carries a property, independent
of how the genome is chopped.

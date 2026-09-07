# Supplementary Note 1: Why length-weighted ratios are invariant to fragmentation

Consider a genome of total length \(L\) and a property \(P\) that is defined on
intervals of the genome (for example, "the interval is inverted relative to the
reference"). Let \(A\subseteq[0,L]\) be the set of positions that have property
\(P\), and define the population fraction

\[
f \;=\; \frac{|A|}{L} .
\]

Now observe the genome through a process that cuts it into \(K\) disjoint
fragments \(F_1,\dots,F_K\). The observation does not change the underlying
positions; it only changes which adjacencies are visible. If we recompute the
fraction from the fragments,

\[
\hat f \;=\; \frac{\sum_{k=1}^{K} |A\cap F_k|}{\sum_{k=1}^{K} |F_k|}
       \;=\; \frac{|A|}{L}
       \;=\; f ,
\]

because the numerator and denominator are both additive over the partition.
Therefore any statistic of the form

\[
\frac{\text{total length with property }P}{\text{total length}}
\]

is invariant to the number of fragments \(K\) and to where the cuts fall.

---

**Contrast with transition counts.** A transition-count statistic asks how many
times the state changes along a path through the fragments. If a true junction is
located inside a fragment, it is visible; if it falls on a cut boundary, it may be
absent from every fragment or, depending on implementation, counted as an extra
junction. Either way, the count acquires a term proportional to \(K\). In
Syn2b's case,

\[
\mathbb{E}[\text{observed transitions}]
\;=\;
T_{\text{true}} + c\,(K-1) + \cdots,
\]

where \(c\) depends on whether a missing adjacency is treated as a denied
adjacency. Syn2b's `raw_inverted_fraction` is a length-weighted ratio and
therefore has \(c=0\); `breakpoint_count` is a transition count and is
fragmentation-dependent.

This is not a special feature of restriction-enzyme tags. The same statement
holds for contigs, nucmer 1-to-1 blocks, or any other segmentation of the genome.

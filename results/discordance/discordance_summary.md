# ANI-synteny discordance across GTDB-R207

Rearrangement truth: dnadiff inversion event counts (.report), robust to the
assembly orientation-convention artifact that drives inverted *aligned-fraction*
metrics toward 0.5 on draft assemblies. Discordant = ANIm >= 97% & >= 2 inversions.

Combined pairs: 47,748 (held-out 43,312 after removing 0 overlap; ANIm-verified high-ANI sample 4,436)

## Inversion-event rates by ANIm

| ANIm range | n | median inversions | % >=1 | % >=2 | % >=5 |
|---|---:|---:|---:|---:|---:|
| 95+ | 5,090 | 1 | 54.2 | 43.2 | 23.5 |
| >=97 | 3,828 | 0 | 45.1 | 33.8 | 15.6 |
| >=99 | 2,169 | 0 | 29.1 | 18.4 | 7.7 |

Discordant pairs (ANIm >= 97, >= 2 inversions): 1,295

## Syn2b flagging (junction count, ANIm >= 97)
 threshold_inversion_events  n_pairs  n_discordant  syn2b_junctions_cutoff  sensitivity_pct  fpr_pct      auc
                          2     3828          1295                    12.0         29.72973  4.57955 0.799766

## Fragmentation vs dnadiff breakpoint counts (held-out, I3)
              metric  spearman_vs_max_contigs  partial_corr_controlling_anim     n
 dnadiff breakpoints                 0.279082                       0.335073 43312
Syn2b junction count                 0.022384                       0.014242 43312

Top discordant pairs: /Users/macstudio/Downloads/Syn2b-paper/results/discordance/discordant_pairs_top.tsv

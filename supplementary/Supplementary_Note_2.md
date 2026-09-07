# Supplementary Note 2: Error model and per-pair standard errors

## Model

For a pair with `raw_inverted_fraction` estimate \(\hat p\) and \(m\) shared
restriction-enzyme landmarks (`shared_tags`), the squared standard error is

\[
\mathrm{SE}(\hat p)^2 \;=\; \frac{1.504\,\hat p(1-\hat p)}{m} \;+\; 0.0205^2 .
\]

The model was fitted by binning the GTDB-R207 held-out set by \(m\) and
regressing the observed error variance against the binomial sampling term
\(\hat p(1-\hat p)/m\). The fit uses 12 bins and gives \(R^2 = 0.9988\).

## Interpretation of the coefficients

* **1.504** is a design effect. If landmarks were independently sampled, the
coefficient would be 1. The observed value near 1.5 reflects spatial clustering
of landmarks inside inverted segments.
* **0.0205** is a method floor. dnadiff averages the inverted aligned fraction
over aligned bases, whereas Syn2b averages over shared landmarks. The two
sets of positions are not identical, so even with infinite landmarks there is a
residual discrepancy of approximately 2 percentage points.

## Predictive performance

| set | \(n\) | observed SD(err) | model SD(err) | SD(\(z\)) | within \(\pm 2\) SE |
|---|---:|---:|---:|---:|---:|
| held_out_50k | 43,312 | 0.0555 | 0.0546 | 1.006 | 95.3% |
| high_ani_all | 6,922 | 0.0848 | 0.0874 | 0.763 | 97.9% |

The model was fitted only on held_out_50k, yet it predicts the aggregate spread
of the independent high-ANI sample without re-fitting.

## Per-pair SE table

| shared tags | SE at \(p=0.5\) | sampling share of variance |
|---|---:|---:|
| 50 | 0.0891 | 95% |
| 100 | 0.0647 | 90% |
| 250 | 0.0417 | 73% |
| 500 | 0.0300 | 51% |
| 1,000 | 0.0228 | 30% |
| 2,000 | 0.0182 | 14% |
| 5,000 | 0.0153 | 5% |

At high shared-tag counts the 0.0205 floor dominates; at low counts the sampling
term dominates. Because `shared_tags` is emitted for every pair, each Syn2b
output can carry its own SE without recomputing the model.

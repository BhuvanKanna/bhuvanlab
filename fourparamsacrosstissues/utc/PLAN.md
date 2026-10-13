# Upper-tail compression (UTC): analysis plan

**Status: locked 2026-10-13, committed before any UTC score was computed on real data.**
Any change after this commit is listed under "Deviations" at the bottom, with
the reason, so the evaluation stays honest.

## Hypothesis

Genes whose over-expression is deleterious have a compressed **upper** tail of
expression across healthy donors, beyond what their expression level predicts.
Only the upper tail is predicted to change.

## Why this metric

RTI, LTI, their σ-distances and newRTI all failed on expression-matched
controls. They depend on the single most extreme donor, on a Gaussian fit that
can go degenerate, and on the whole distribution, so donors dragged into the
**lower** tail by tissue composition or sample quality read as a capped upper
tail. Spearman(mean, skew) is about −0.8 in every tissue.

UTC uses only the upper half of the distribution, uses quantiles rather than a
fit, and is adjusted for expression level by construction.

Feasibility (2026-10-13, 10 tissues): pTriplo ≥ 0.94 genes reach a median
q97.5 − q50 of 0.71 log2 (1.5× = 0.585), so donor variation does reach
duplication dosage. An unmatched first look showed pTriplo genes narrower in
**both** tails, with the same up/down asymmetry as other genes. That look is
why the variability metric (S1) is pre-specified as a comparison.

## Data and eligibility (fixed)

- Matrices: `data/v11_log2_<tissue>.csv.gz`, all donors, converted to
  `x = log2(TPM + 1)`. **Raw values, no exclusion at −1.**
- Genes: GENCODE v47 `protein_coding` only (`docs/genetype/`).
- Tissues: those with **≥ 150 donors**.
- Gene × tissue eligible when median TPM ≥ 5 (so differences are close to
  log2 fold changes and the floor never touches the upper half).
- A gene gets a score when eligible in **≥ 5 tissues**.

## Metrics (fixed)

With per-gene quantiles `q` of `x` in one tissue:

| metric | definition | prediction for dosage-sensitive genes |
|---|---|---|
| **UTC** (primary) | `ln[(q97.5 − q50) / (q75 − q50)]` | lower |
| LTC (specificity) | `ln[(q50 − q2.5) / (q50 − q25)]` | unchanged for pTriplo; lower for pHaplo |
| VAR (secondary, S1) | `ln(q75 − q25)` | lower (general constraint, not truncation) |

A gene × tissue value is missing if a denominator is ≤ 0.

**Expression adjustment, within each tissue.** Sort eligible genes by `q50`.
Residual = value − rolling median of the same metric over the 201 genes nearest
in `q50` (centered window, truncated at the ends). z = residual / (1.4826 ×
MAD of all residuals in that tissue). Donor count is constant within a tissue,
so it needs no adjustment there.

**Aggregation.** Gene score = median z across its eligible tissues.

## Evaluation (fixed)

Expression for matching = mean of the gene's per-tissue `q50` over its eligible
tissues.

- **Positives**: `genelists/pos_pTriplo.txt` (pTriplo ≥ 0.94).
  **Negatives**: `genelists/neg_dup_tolerant.txt`. Protein-coding, scored genes
  only; any gene in both sets is dropped.
- **Matching**: 1:1 nearest neighbor on expression, caliper 0.1 log2, without
  replacement, positives processed in random order (seed 0).
- **AUC** = P(positive score < negative score) over matched pairs, ties count
  ½. 95% CI from 2,000 bootstrap resamples of pairs (seed 0).

| id | test | supports the hypothesis when |
|---|---|---|
| **P1 (primary)** | UTC AUC, matched pairs | 95% CI lower bound > 0.5 |
| P2 | UTC AUC − LTC AUC, same pairs, paired bootstrap | 95% CI excludes 0, positive |
| P3a | partial Spearman(UTC, pTriplo \| pHaplo, expression), all scored genes in Collins | < 0, 95% bootstrap CI excludes 0 |
| P3b | partial Spearman(LTC, pHaplo \| pTriplo, expression) | < 0 (the mirror) |
| P4 | P1 in each half of the gene set, split by parity of the ENSG number | point estimate > 0.5 in both |
| S1 | VAR AUC, same pairs | reported; not a test of truncation |

Partial Spearman: rank-transform all variables, regress the two of interest on
the rank covariates by OLS, correlate residuals.

**The hypothesis is supported only if P1 holds and P2 holds.** P1 alone could
be general variability; P2 is what makes it upper-specific.

## Validation, run before the evaluation (fixed)

Simulated `x` with n = 300 donors, 2,000 replicates per scenario, UTC computed
exactly as above (no residualization):

- V1: Gaussian, untruncated (reference).
- V2: Gaussian truncated above at μ + 1σ, + 1.5σ, + 2σ. UTC must fall below V1.
- V3: composition-like contamination, 85% N(0,1) + 15% N(−2,1). UTC must stay
  within ±0.05 of V1. This is the failure mode that sank RTI.

Known limitation, stated in advance: a distribution with an intrinsically short
right tail (e.g. a negatively skewed skew-normal) has a compressed upper tail
too, and UTC cannot tell it apart from soft truncation.

## Deferred

- Donor-covariate adjustment (sex, age, Hardy scale, RIN, ischemic time):
  the GTEx `SampleAttributesDS` / `SubjectPhenotypesDS` files are not in this
  repo. Run as a sensitivity analysis once they are.

## Deviations

None yet.

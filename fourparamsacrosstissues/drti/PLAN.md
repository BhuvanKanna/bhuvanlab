# dRTI: donor-level right truncation index — plan

**Status: locked 2026-10-13, after the viability check and before the
benchmark, metric selection or any comparison with dosage labels.** Later
changes go under "Deviations" with the reason.

## Goal

One per-gene score that measures **right-truncation shape**: the upper end of
a gene's expression distribution across people is a ceiling (density lost
abruptly) rather than a smooth tail, beyond what expression level, skew,
donor count, reliability and known confounders produce.

Two separate claims, tested separately:

1. **Measurement validity** (this plan, phases 1–5): the score responds to
   truncation and not to confounds. It is chosen and validated on
   semi-synthetic spike-ins in real genes, **without dosage labels**.
2. **Biological validity** (phase 6): whether dosage-sensitive genes score
   higher. Tested once, after the metric is frozen.

## What failed before, and why

| metric | failure |
|---|---|
| RTI / LTI (4-param fit) | single donor sets the ceiling; degenerate fits; tracks expression (ρ ≈ ±0.6) |
| newRTI (tail model) | reads left skew as truncation; below chance on matched controls |
| UTC (quantile ratio, per tissue) | at chance (AUC 0.501); robust to composition but measured per tissue |

**Root cause (viability, below):** truncation acts on a person's underlying
expression level, but each tissue measures that level plus a large
tissue-specific noise term (composition, sampling, sample quality). In one
tissue the ceiling is blurred to almost nothing. Averaging a donor's tissues
recovers it.

## Generative picture

For gene g, donor d, tissue t, residual expression after covariates and PCs:

    z[d,t] = λ_t · L_d + e[d,t]

L_d is the donor's shared level (genetics plus unmeasured systemic state).
e[d,t] is tissue-specific. Lethal over-expression removes donors with high
L_d, so the truncation lives in L, and the best measurement object is an
estimate of L_d.

Known limitation: a biochemical ceiling (e.g. a transcription rate that
saturates) produces the same shape. Shape alone cannot tell it from
selection. Expression-level calibration reduces, but does not remove, this.

## Viability (done 2026-10-13, 3,000 random protein-coding genes, seed 0)

`drti_prep.py` then `drti_viability.py`.

- 49 tissues with ≥ 100 donors. Median gene: 47 tissues, 934 donors with
  ≥ 5 tissues expressing it.
- ICC (donor-shared share of one tissue's variance): median 0.13 after
  covariates + 10 PCs (0.14 before).
- Spike-in: tissues split per gene into halves A and B; donors in the top 10%
  of their A-half mean removed (truncation on a latent level, measured on
  independent tissues). UTC drop in SDs of random removal of as many donors:

| measured on | median effect | detected (> 1.645) |
|---|---|---|
| one B tissue (median over tissues) | 0.18 | 2.1% |
| donor mean over B tissues | 1.61 | 49.5% |

## Phase 1: preprocessing (fixed)

`drti_prep.py`: per tissue with ≥ 100 donors, x = log2(TPM + 1) of all
donors, protein-coding genes with median TPM ≥ 1. Remove covariates (sex,
age, Hardy scale, RIN, ischemic time, exonic rate, center) and then the top
10 PCs of those residuals, then robust-standardize per gene
((r − median) / 1.4826·MAD). The excluded-table default does not apply here:
TPM ≥ 1 medians and the regression need every donor.

## Phase 2: candidates

**Objects** (from a gene's tissues; donors with ≥ 5 tissues):

| id | object |
|---|---|
| `mean` | donor mean of z |
| `blup` | mean × σL² / (σL² + σw²/k_d), variance components by one-way ANOVA |
| `factor` | weighted mean, w_t = λ_t / (1 − λ_t²), λ_t = corr(z_t, leave-one-out donor mean) clipped to [0.05, 0.95] |

**Statistics** (all oriented so higher = more truncated):

| id | statistic |
|---|---|
| `utc95`, `utc975`, `utc99` | −ln[(q_p − q50) / (q75 − q50)], p = 95, 97.5, 99 |
| `xi_mom` | −ξ, Dekkers–Einmahl–de Haan moment estimator on the top 10% (replaced by `xi_pwm`, see Deviations) |
| `xi_gpd` | −ξ, generalized Pareto MLE on exceedances over q90 |
| `asym` | −(UTC975 − LTC975), upper relative to lower |
| `edge_lr` | `mean` object only: signed √(2·ΔlogL) of a Gaussian latent truncated above at c plus known noise s_d² = σw²/k_d, against no truncation; also reports the missing fraction 1 − Φ((c − μ)/σ) |

**Calibrations:** `raw`, and `knn`: residual from the median of the
statistic among the 200 nearest genes in (expression, ICC, log donors)
standardized space, computed from **control** condition values and divided by
their 1.4826·MAD.

## Phase 3: benchmark (fixed)

Gene pool: protein-coding genes with ≥ 10 tissues, **excluding every gene in
any of the four published control sets** (pTriplo, DECIPHER, duplication-
tolerant, olfactory). Random sample of 4,000 (seed 1). Each gene's tissues
are split into halves A and B (seeded by gene id). Objects are built from B
only; truncation is defined from A only. Analysis donors are those with ≥ 3 A
tissues and ≥ 5 B tissues.

Conditions:

| id | what | should |
|---|---|---|
| `control` | as is | — |
| `t05`, `t10`, `t20` | remove donors in the top 5/10/20% of the A-half mean | fire |
| `hard05` | remove top 5% of the B `mean` object itself (upper bound) | fire |
| `comp_tissue` | in each B tissue, an independent random 15% of donors shifted −2 | not fire |
| `comp_systemic` | one random 15% of donors shifted −1 in every B tissue | not fire |
| `skew` | each B tissue mapped through g(z) = (1 − e^(−0.5z))/0.5, re-standardized | not fire |
| `outliers` | 1% of donors +4 in one random B tissue each | not fire |
| `fewer` | random 30% of donors removed | not fire |

Truncation conditions are compared with removing the same number of donors at
random (`control` uses that random removal when paired with a truncation).

**Scoring a candidate** (object × statistic × calibration):

- **Power**: genes split 50/50 at random (seed 2) into condition vs control;
  AUC of the score separating them **across genes** (the realistic use).
  Main power = mean of t05 and t10.
- **False fire**: the same AUC for each confound; FF = max |AUC − 0.5|.
- **Confound correlation**: on control values, CC = max |Spearman| with
  expression, ICC and donor count.

## Phase 4: selection rule (fixed)

1. Admissible: FF ≤ 0.05 and CC ≤ 0.10.
2. Best single = highest main power among admissible candidates.
3. Stacked: logistic regression on the admissible candidates' scores, trained
   on genes with even ENSG number (positives t05 + t10, negatives control +
   all five confounds), evaluated on odd. Adopted only if it is admissible
   on the odd half and beats the best single's main power there by ≥ 0.02.
4. If nothing is admissible, report that, relax nothing silently, and record
   the closest candidates.

## Phase 5: final score

Run the chosen recipe on every eligible gene using all of its tissues (no A/B
split). Output per gene: `drti` (calibrated, higher = more truncated),
`missing_frac` from `edge_lr`, a 200-replicate donor-bootstrap SE, `icc`,
`n_donors`, `n_tissues`. Genes with < 10 tissues get no dRTI (tissue-level
power is ~chance, so a score would mislead).

## Phase 6: biological test (fixed now, run once after phase 5)

Same design as `utc/PLAN.md`, on dRTI:

- P1: pTriplo ≥ 0.94 vs duplication-tolerant, 1:1 matched on expression
  (caliper 0.1) **and ICC (caliper 0.03)**; AUC = P(positive more truncated);
  supported if the 95% bootstrap CI lower bound > 0.5.
- P2: partial Spearman(dRTI, pTriplo | pHaplo, expression, ICC) > 0, CI
  excludes 0.
- P3: P1 on each ENSG-parity half, point estimate > 0.5 in both.
- S1: the same test for donor-level variability (log IQR of `mean`),
  reported for comparison with `utc` S1.

**dRTI is supported as a dosage marker only if P1 and P2 hold.** Measurement
validity (phases 3–4) stands on its own either way.

## Deviations

1. **`xi_mom` → `xi_pwm`** (2026-10-13, before the benchmark ran). The
   Dekkers–Einmahl–de Haan moment estimator works on log(X), so it depends on
   where zero is and cannot be applied to standardized residuals that go
   negative. Replaced with the probability-weighted-moments GPD estimator
   (Hosking & Wallis 1987) on exceedances over q90, which is location- and
   scale-invariant: −ξ = a0/(a0 − 2a1) − 2.
2. **Realistic-strength confounds, added as a second benchmark**
   (2026-10-13, after a 60-gene dry run and **before** the full benchmark
   finished). The dry run showed every candidate firing on `skew` and
   `comp_systemic`. Checked against 2,290 real genes (viability sample),
   both confounds as locked are far outside the data: `skew` (a = 0.5) gives
   a median donor-level Bowley skewness of −0.17, more left-skewed than 98.7%
   of real genes (real median +0.004); `comp_systemic` (shift −1) gives a
   median lower-tail weight (q50 − q2.5)/(q75 − q25) of 1.89, heavier than 97%
   of real genes. A realistic-strength variant sets each so that the median
   gene sits at the real 95th percentile: `skew` a = 0.25 (Bowley −0.083;
   real 5th percentile −0.085) and `comp_systemic` shift −0.9 (weight ≈ 1.82;
   real 95th percentile 1.822). Everything else is unchanged.
   **Rule, fixed now:** the locked benchmark is reported in full and stays
   primary. If it leaves no admissible candidate, selection uses the same rule
   on the realistic-strength benchmark, and the write-up states that the
   chosen metric is not robust to the locked (extreme) strengths.

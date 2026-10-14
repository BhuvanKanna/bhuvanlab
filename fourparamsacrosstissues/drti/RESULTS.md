# dRTI: results

**Outcome (2026-10-14): no candidate is admissible, on either the locked or
the realistic-strength benchmark. Under step 4 of the selection rule, no dRTI
is chosen, so phase 5 (final scores) and phase 6 (the biology test) do not
run.** Nothing was relaxed after seeing results. See `PLAN.md` for the locked
design and its two deviations.

## What was run

`drti_bench.py` on the 49-tissue prep (covariates + 10 PCs, robust z).
4,000 genes were sampled (seed 1) from the 13,301 protein-coding genes with
≥ 10 tissues, after excluding the 2,477 genes in the four control sets.
3,897 had ≥ 100 analysis donors (≥ 3 A-half and ≥ 5 B-half tissues) and were
benchmarked. Median gene: 867 analysis donors, ICC 0.137, expression 3.48.

38 candidates: 3 objects × 6 statistics × 2 calibrations, plus `edge_lr` × 2.
Power, false fire (FF) and confound correlation (CC) as in the plan.
Admissible = FF ≤ 0.05 and CC ≤ 0.10.

| file | contents |
|---|---|
| `bench_raw.csv.gz` | every statistic, every gene × 14 conditions (locked strengths) |
| `bench_raw_realistic.csv.gz` | the `skew` / `comp_systemic` rows at realistic strength |
| `bench_scores.csv`, `bench_scores_realistic.csv` | per-candidate AUCs, FF, CC |
| `selection.json`, `selection_realistic.json` | `n_admissible: 0`, `chosen: null` |

## Locked benchmark

Top by power, and the most confound-resistant (AUC, 0.5 = chance):

| candidate | t05 | t10 | hard05 | power | comp_tissue | comp_systemic | skew | outliers | fewer | FF | CC |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `raw:mean.edge_lr` | 0.685 | 0.708 | 0.957 | **0.696** | 0.551 | 0.776 | 0.989 | 0.498 | 0.542 | 0.489 | 0.198 |
| `knn:mean.edge_lr` | 0.681 | 0.696 | 0.919 | 0.689 | 0.543 | 0.747 | 0.950 | 0.495 | 0.549 | 0.450 | 0.016 |
| `knn:blup.utc99` | 0.667 | 0.691 | 0.972 | 0.679 | 0.597 | 0.644 | 0.946 | 0.490 | 0.503 | 0.446 | 0.009 |
| `knn:mean.utc99` | 0.665 | 0.682 | 0.986 | 0.674 | 0.584 | 0.640 | 0.957 | 0.493 | 0.500 | 0.457 | 0.012 |
| `knn:blup.xi_gpd` | 0.639 | 0.633 | 0.947 | 0.636 | 0.529 | 0.529 | 0.730 | 0.503 | 0.508 | 0.230 | 0.028 |
| `knn:blup.xi_pwm` | 0.622 | 0.619 | 0.912 | 0.621 | 0.524 | 0.520 | 0.714 | 0.490 | 0.481 | **0.214** | 0.019 |

Every candidate fires on `skew` (AUC 0.71–0.99). No candidate is admissible.

## Realistic-strength benchmark (deviation 2)

`skew` a = 0.25 and `comp_systemic` shift −0.9, each set so the median gene
sits at the real 95th percentile. All other conditions are identical.

| candidate | power | comp_systemic | skew | FF | CC |
|---|---|---|---|---|---|
| `knn:mean.edge_lr` | 0.689 | 0.714 | 0.863 | 0.363 | 0.016 |
| `knn:blup.utc99` | 0.679 | 0.633 | 0.804 | 0.304 | 0.009 |
| `knn:blup.xi_gpd` | 0.636 | 0.528 | 0.638 | 0.138 | 0.028 |
| `knn:blup.xi_pwm` | 0.621 | 0.517 | 0.624 | **0.124** | 0.019 |

Still no admissible candidate. The lowest FF is 0.122 (`raw:blup.xi_pwm`),
more than twice the 0.05 limit, and it comes entirely from `skew`. Its other
four confounds are all within 0.03 of chance.

## What this shows

1. **Donor-level truncation is detectable.** A 5–10% cut on the latent level,
   measured on independent tissues, is picked up at AUC 0.62–0.70 across
   genes, and a hard 5% ceiling at 0.91–0.99. Single-tissue metrics were at
   chance (viability: 2.1% of genes detected from one tissue vs 49.5% from the
   donor mean).
2. **Calibration solves the confounding by level and reliability.** kNN
   calibration takes CC to ≤ 0.03 for every top candidate (raw `edge_lr`:
   0.198 → 0.016).
3. **Composition splits the statistics.** The quantile-ratio statistics
   (`utc*`, `asym`) and `edge_lr` fire on a 15% systemic lower shift (AUC
   0.62–0.78), because it moves the median and the IQR they are scaled by. The
   pure upper-tail shape estimators (`xi_pwm`, `xi_gpd` on exceedances over
   q90) do not (0.52–0.53).
4. **Skew is the barrier, and it looks fundamental rather than a matter of
   tuning.** A gene with natural left skew at the real 5th percentile is read
   as truncated by every statistic (AUC 0.62–0.90), and by the best of them
   about as strongly as a real 5–10% cut is detected (0.62 vs 0.62). The
   reason is the measurement itself: the donor mean is the latent level plus
   noise of variance σw²/k. Convolving a truncated latent with that noise turns
   the edge into a smoothly compressed upper tail, which is exactly what a
   natural left skew is. At the median ICC (0.14), nothing in the observed
   shape tells the two apart. This is the known limitation stated in the
   plan (selection vs saturation), measured.

So the best statement the data supports is: **a right-truncation shape score
exists, but it is a skew-confounded one.** The closest candidate,
`blup.xi_pwm` (−ξ of a GPD fit by probability-weighted moments to the BLUP
donor level above its 90th percentile), is robust to expression, ICC, donor
count, composition, outliers and donor loss. It cannot tell a ceiling from a
naturally short upper tail.

## Not run

- Phase 5 (`drti_final.py`) and phase 6 (`drti_evaluate.py`) are written and
  tested (`tests/test_drti.py`) but were not run, because the rule chose
  nothing. Running phase 6 on a non-admissible candidate would not be the
  pre-registered test. If it is ever run, it must be labeled exploratory, and a
  pTriplo difference could not be attributed to truncation over skew.

## Operational notes

- The first full run (10 workers) was stopped by the OS memory guard at 2,500
  genes, with nothing saved. `drti_bench.py` now checkpoints every 500 genes in
  `bench_parts/` (not committed; `bench_raw.csv.gz` is their concatenation) and
  defaults to 4 workers. The second run was also stopped once, after 6 of 8
  batches, and resumed from them. The checkpoints contain identical rows
  whichever run wrote them, because every random draw is seeded by gene id.

# Upper-tail compression (UTC): results

Run 2026-10-13, exactly as locked in `PLAN.md` (commit `4f83a6d`, made before
any score was computed). No deviations.

```bash
cd fourparam
python simulate_utc.py
python utc.py --jobs 4
python evaluate_utc.py --collins <collins2022-TableS7-pHaplo-pTriplo.tsv>
```

## Verdict: not supported

The pre-specified rule (P1 and P2) fails at P1. Expression-matched pTriplo
genes do not have a more compressed upper tail than duplication-tolerant genes.
The signal that does exist is **lower overall variability** (S1), which is not
specific to the upper tail.

## Validation (before evaluation)

n = 300, 2,000 replicates (`simulation.csv`).

| scenario | UTC mean | LTC mean |
|---|---|---|
| Gaussian | 1.056 | 1.059 |
| truncated at +1σ | 0.735 | 1.102 |
| truncated at +1.5σ | 0.848 | 1.075 |
| truncated at +2σ | 0.960 | 1.061 |
| 15% of donors shifted −2σ (composition-like) | 1.024 | 1.197 |

V2 pass (truncation lowers UTC), V3 pass (contamination moves UTC −0.032).
Per gene × tissue, a +2σ cut is only about 1 SD of noise; power comes from
aggregating across tissues.

## Data

48 tissues with ≥ 150 donors, 467,276 eligible protein-coding gene × tissue
values, 13,487 genes scored in ≥ 5 tissues. Positives: 1,317 scored; negatives:
174 scored (most of the duplication-tolerant list is not protein-coding or not
expressed at ≥ 5 TPM); 3 genes in both lists dropped. 170 matched pairs,
expression AUC 0.500.

## Tests

AUC = P(positive scores lower than its matched negative).

| test | result | locked criterion |
|---|---|---|
| **P1** UTC AUC | **0.501 [0.439–0.560]** | **FAIL** (lower bound > 0.5) |
| P2 LTC AUC | 0.375 [0.318–0.434] | |
| P2 UTC − LTC | +0.126 [+0.037, +0.214] | pass, but see below |
| P3a partial ρ(UTC, pTriplo \| pHaplo, expr), n = 12,813 | −0.053 [−0.071, −0.036] | pass |
| P3b partial ρ(LTC, pHaplo \| pTriplo, expr), n = 12,813 | −0.029 [−0.045, −0.011] | pass |
| P4 UTC AUC by ENSG parity | even 0.487 (95 pairs), odd 0.478 (75 pairs) | FAIL |
| S1 VAR AUC | **0.665 [0.613–0.718]** | reported |

## Reading it

- **P2 passes for the wrong reason.** The gap comes from LTC (0.375): pTriplo
  genes have a relatively *heavier* lower tail, not a compressed upper one.
  Across genes LTC and VAR correlate at −0.56, so this is mostly what a
  narrower core looks like: when the IQR is tight, the same few donors
  dragged down by composition or sample quality sit relatively further out.
- **P3a is statistically clear but very small.** ρ = −0.05 explains about
  0.3% of rank variance, in the predicted direction. With 12,813 genes almost
  any residual confound reaches significance, and the matched comparison (P1)
  and both halves (P4) show nothing. It is not evidence of truncation on its
  own.
- **S1 is the real signal.** pTriplo genes are less variable across donors
  than expression-matched tolerant genes (AUC 0.665). That agrees with the
  feasibility check and with published work on dosage-sensitive genes. It
  says expression is held tighter in both directions, which is a different
  claim from the truncation hypothesis.
- **Power caveat.** Only 170 matched pairs, so P1's CI is ±0.06. A real UTC
  AUC of about 0.55 could have been missed. The continuous P3a, which uses
  every gene, bounds the effect as small if it exists.

## Note

Gene-level UTC still correlates with mean expression (Spearman +0.20) even
though each tissue was residualized, because aggregation across tissues
reintroduces it. P1 is matched (expression AUC 0.500) and P3 includes
expression as a covariate, so neither is affected, but any new use of
`gene_scores.csv` should adjust for `expr`.

## Files

- `gene_tissue.csv.gz`: q50, raw UTC/LTC/VAR and their z per eligible gene × tissue.
- `gene_scores.csv`: per gene, `n_tissues`, `expr`, median-z `utc`, `ltc`, `var`.
- `evaluation.json`: every number above.

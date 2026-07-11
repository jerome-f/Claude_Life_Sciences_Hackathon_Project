# DBNascent Concordance Test — Does Signed Co-transcription Predict Causal Direction?

**Date:** 2026-07-11
**The gate:** the entire signed-E2G premise rests on one question — does the *sign* of DBNascent nascent co-transcription (pcc) agree with the *sign* of a measured causal perturbation? This test answers it against the ENCODE-rE2G CRISPRi benchmark (Gschwind et al. 2023), the field-standard causal truth set.

## Design (and why the naive test is wrong)

**Two independently-derived signs on the same (enhancer, gene) pairs:**
- **DBNascent:** `sign(pcc)` — positive = enhancer co-transcribed with gene (activating).
- **CRISPRi:** the assay *represses* the element, so a true activating enhancer shows expression **decrease** on knockdown → **enhancer sign = −sign(EffectSize)**.

If the label tracks causation, `pcc` and `EffectSize` should be **negatively correlated**.

**The base-rate trap (handled):** both distributions are ~80–85% activating (85% of DBNascent pairs positive; 71% of CRISPRi effects activating). A raw sign-match rate is therefore ~64% *by chance alone* and is meaningless. The load-bearing statistics are the **continuous Spearman ρ** (uses magnitude, immune to sign base rate) and the **cell-type-matched odds ratio**, not concordance %.

**Data:** 666 significant matched pairs (pValueAdjusted<0.05), from 14,734 valid CRISPRi pairs joined to DBNascent bidirectional→gene pairs by genomic overlap + gene symbol, across K562/HCT116/WTC11 and others.

## Result: QUALIFIED PASS — the label tracks *direction*, weakly but robustly

### Primary — continuous concordance
- **Spearman ρ = -0.224** (95% CI [-0.292, -0.150], p = 4.8e-09). Negative, i.e. **the predicted direction**: higher co-transcription → larger expression drop when the element is silenced.
- **Permutation test** (5,000 label shuffles): p = 2.0e-04. The correlation is not a base-rate artifact.
- **Survives distance control:** partial ρ | distance-to-TSS = **-0.128** (p = 9.6e-04). Not driven by the fact that closer pairs are both more co-transcribed and more perturbable.

### Cell-type matching sharpens it (hallmark of real biology)
- **DBNascent blood × K562 CRISPRi: ρ = -0.286** (n=326), Fisher **OR = 5.4** (p = 2.2e-04). Matching the tissue to the perturbed cell line strengthens the signal — an artifact would not respect cell-type.
- Negative in every cell type tested: K562 ρ=−0.16 (n=500), HCT116 ρ=−0.47 (n=30), WTC11 ρ=−0.12 (n=113).

### Robustness controls (all negative ρ)
| Filter | n | ρ | 
|---|---|---|
| Pooled, significant | 666 | -0.224 |
| Well-powered (Power@ES25>0.8) | 469 | -0.218 |
| H3K27ac enhancers only | 591 | -0.215 |
| Direct effects only (>0.8) | 473 | -0.074 |

## The three honest caveats (these bound what can be built)

1. **The effect is modest (|ρ| ≈ 0.13–0.22).** DBNascent gets direction right *on average*, but the correlation is weak. This matches the negatives-test preview (1.43× enrichment) — the label carries real but noisy causal signal.

2. **Sign only, not magnitude.** Among *regulated* pairs, |pcc| does **not** track |EffectSize| (ρ = 0.070, p = 0.13, n = 466). The label licenses a **signed** edge, not a quantitative effect-size predictor.

3. **The repression arm is underpowered.** Only 69 significant *repressive* CRISPRi pairs and 102 DBNascent-negative pairs exist in the overlap; the both-repressive cell of the 2×2 has just 12. We can confirm the label captures **activating** direction; we **cannot** yet confirm it captures **repressive** direction. The confident-sign Fisher OR (1.9, p=0.13) is not significant precisely because of this.

## Verdict for the project

**The premise survives — with a firm ceiling.** Signed nascent co-transcription predicts causal *direction* at genome scale, significantly, robustly, and in the correct sign, and cell-type matching improves it. This **licenses a signed E2G edge** — but the design must respect three limits the data just set:
- Predict **sign**, not magnitude (an abstention/confidence class, not a regression on effect size).
- Expect a **weak ceiling** (~ρ 0.2 against causal truth); do not over-claim calibration.
- Treat the **repressive class as unvalidated** until a repression-enriched perturbation set is added.

This is a **publishable finding on its own**: "genome-scale nascent co-transcription sign agrees with CRISPRi causal direction (ρ≈−0.22, cell-type-matched OR up to 5.4), but carries no effect-magnitude information and is unvalidated for repression." That result stands whether or not the full model is ever built — exactly the deliverable framing recommended earlier.

## Outputs
| File | Description |
|---|---|
| `concordance_matched_pairs.parquet` | 666 significant CRISPRi×DBNascent matched pairs (the test set) |
| `concordance_results.json` | Primary + cell-type + all confound-control statistics |
| `concordance_extra.json` | Permutation, partial-correlation, magnitude, per-cell-type |
| `dbnascent_concordance_test.png` | 3-panel figure (scatter, robustness forest, base-rate caveat) |

CRISPRi benchmark staged to `gs://claude_hackathon/crispri_benchmark/20260711/`; results under `gs://claude_hackathon/dbnascent/20260711/`.


---

## Addendum — Magnitude concordance in commensurate (z-score) terms

**Motivation (user):** the sign test is rank-based (unit-free), but the *magnitude* test in the primary analysis compared a raw correlation coefficient (`pcc`) against a raw log-fold-change — incommensurate scales, and `pcc` ignores per-pair measurement precision. DBNascent carries a precision-aware standardized statistic — the signed **t-statistic** (`t`), which upweights high-nObs pairs. Re-tested with `t`, cell-type-matched (blood<->K562), and restricted to well-measured pairs.

**Result: the poor magnitude concordance was largely a measurement-noise artifact.**

| Comparison | DBNascent metric | n (regulated) | Magnitude rho (abs metric vs abs effect) | p |
|---|---|---|---|---|
| Pooled, raw (primary test) | abs(pcc) | 466 | 0.070 | n.s. |
| Blood x K562 matched | abs(pcc) | 291 | 0.118 | 0.044 |
| Blood x K562 matched | abs(signed t) | 291 | 0.148 | 0.012 |
| **Blood x K562, well-measured both sides** | **abs(pcc)** | 159 | **0.257** | 0.0011 |
| **Blood x K562, well-measured both sides** | **abs(signed t)** | 159 | **0.324** | 3.1e-05 |

"Well-measured both sides" = DBNascent nObs>=50 AND CRISPRi PowerAtEffectSize25>0.8. Magnitude concordance climbs from **0.07 (n.s.) -> 0.32 (p=3e-5)** once the readouts are put on comparable footing.

**Sign vs magnitude — a metric division of labor.** For the *sign* call, raw `pcc` slightly outperforms the t-statistic (blood x K562: raw Spearman -0.286 vs signed-t -0.189; precision subset: -0.378 vs -0.318) — the t-stat inflates with nObs regardless of directional cleanliness, adding precision noise to a pure sign call. So: **use raw `pcc` for sign, signed `t` on well-measured pairs for magnitude.**

**Revised design implication.** Caveat (2) in the main report is softened: DBNascent *does* carry quantitative magnitude information (rho~0.32), but only on the precision-filtered, cell-type-matched subset. A signed edge may therefore reasonably carry a **calibrated magnitude head**, not merely sign + abstention — trained on well-measured pairs, with confidence tied to measurement depth (nObs). Direction remains the primary, most robust signal.

Output: `concordance_zscore.json` (all sign + magnitude statistics under raw-pcc, signed-t, and Stouffer-combined-z), `zscore_matched_bk.parquet`.

# DBNascent E2G — Processing & QC Report

**Date:** 2026-07-11
**Source:** DBNascent (Sigauke et al., BMC Genomics 2025). Zenodo record 14519113, `DBNascent_data.tar.gz` (4.6 GB), member `bidirectional_gene_pairs/dbnascent_pairs.txt.gz`.
**License:** CC-BY-NC-ND 4.0 — the NoDerivatives clause governs any redistribution of transformed tables; these processed products are staged privately for training/validation use, not public re-release.
**Compute:** GCP genomics-vm (e2-standard-16, 64 GB), duckdb.
**Staged to:** `gs://claude_hackathon/dbnascent/20260711/` (raw / processed / qc).

## 1. Input

The published pairs file is comma-separated, 27 columns, bed12 layout (cols 1–6 gene, 7–12 bidirectional) plus correlation summary statistics. The signed label is `pcc` (Pearson correlation of nascent transcription between a gene and a bidirectional/enhancer within 1 Mb, per tissue).

| Metric | Value |
|---|---|
| Per-tissue observations (rows) | 12,697,055 |
| Unique gene–enhancer pairs | 6,700,460 |
| Unique genes | 27,822 |
| Unique bidirectionals (enhancers) | 589,339 |
| Tissues | 11 |
| Null values in pcc / adj_p / nObs / distance / gene_id | 0 |

## 2. Key QC findings

**The file is pre-filtered to FDR < 0.01.** Max observed `adj_p_BH` = 0.0099998, max `pval` = 0.00163. DBNascent publishes **only significant correlations** — there are no tested-but-uncorrelated pairs in the distributed file. This directly affects the negatives design (§4).

**Sign balance:** 82.3% of observations positive (mean pcc +0.645), 17.7% negative (mean pcc −0.697). Negative-correlation pairs are real and substantial — consistent with true repressive/antisense relationships and with the design premise that E2G sign must be modelled, not assumed positive.

**FDR stringency tiers:** 8,315,182 observations at FDR < 0.001 (the stringent tier the concordance test calls for); the remaining 4,381,873 fall in [0.001, 0.01).

**3D support flags:** 270,279 observations overlap BOTH a GTEx eQTL pair and a PolII ChIA-PET loop; 7,070,242 overlap ChIA-PET only; 338,868 overlap GTEx (any). These define the "3D-supported" subset.

**Distance to TSS** (|distance|, kb): 5th=18, 25th=162, median=409, 75th=693, 95th=937. Pairs span the full ±1 Mb window; median link is ~400 kb from the TSS — far beyond a promoter-proximal window, underscoring the need for the distance-decay prior.

**Genomic position:** 79% of observations are intragenic (enhancer within a gene body), split roughly evenly upstream/downstream of the TSS.

## 3. Cross-tissue sign consistency (the load-bearing statistic)

Aggregating the 6,700,460 unique pairs across tissues:

| Sign class | Pairs | % |
|---|---|---|
| Consistent positive | 4,806,830 | 71.7% |
| Consistent negative | 1,182,741 | 17.7% |
| **Sign-flip (both signs across tissues)** | **710,889** | **10.6%** |

**10.6% of pairs are significant with opposite signs in different tissues.** This is direct empirical confirmation of the design-doc concern: if these were mean-aggregated, ~711k pairs would collapse toward zero and masquerade as non-links. The trainable subset therefore preserves, per pair: `n_tissues`, `n_pos`, `n_neg`, `mean/median/min/max_pcc`, `sd_pcc`, a `sign_concordance` statistic (fraction of tissues agreeing with the majority sign), the ordered `pcc_by_tissue` list, and a `sign_class` label. Among multi-tissue pairs, median sign_concordance = 1.0 (most pairs are unanimous), but the 5th percentile is 0.5 — the conflicted tail the later cell-type layer needs.

## 4. Trainable subset

**Positives** (`dbnascent_e2g_positives.parquet`, 6,700,460 pairs): one row per unique gene–enhancer pair with the cross-tissue consistency block above, `distance_tss`, `position`, and a `stringent_3d` flag (best_adj_p < 0.001 AND (GTEx or ChIA-PET support)) — 2,729,458 pairs (40.7%) qualify as stringent 3D-supported.

**Per-tissue signed table** (`dbnascent_e2g_pertissue_signed.parquet`, 12,697,055 rows): the full cleaned per-tissue observations, retained so the deferred cell-type layer can recover conflicted pairs.

**Negatives** (`dbnascent_e2g_negatives.parquet`, 9,977,500 pairs): because the published file contains no tested-but-uncorrelated pairs, negatives are constructed as **distance-matched hard negatives** — (gene, transcribed-enhancer) pairs within ±1 Mb that are NOT among the significant partners of that gene. Both the gene and the enhancer are transcribed in the atlas, so these are plausible-but-unlinked pairs, not random genomic background.

> **Caveat carried forward:** the within-window non-significant candidate pool (9,977,500) is smaller than a 3:1 target relative to positives (~20M), so all implied negatives were kept — the realised positive:negative ratio is ~1.5:1 and distance-matched down-sampling was not triggered. A true "tested-but-uncorrelated" negative set would require recomputing correlations from the raw count matrices (the 2.4 GB `normalized_counts` / per-sample ZIPs), which was out of scope here.

## 5. Outputs

| File | Rows | Description |
|---|---|---|
| `dbnascent_e2g_positives.parquet` | 6,700,460 | Unique signed pairs + cross-tissue consistency block + stringent_3d flag |
| `dbnascent_e2g_pertissue_signed.parquet` | 12,697,055 | Full per-tissue signed observations |
| `dbnascent_e2g_negatives.parquet` | 9,977,500 | Distance-matched hard negatives (within-window, non-significant) |

All staged under `gs://claude_hackathon/dbnascent/20260711/processed/`, with QC tables/figures under `.../qc/` and raw inputs under `.../raw/`.

# Multi-Biobank Credible-Set Inventory — QC Report

**Date:** 2026-07-11
**Build:** hg38 (all sources; FinnGen `23`→`X` normalized; 18 metabolomics rows with missing hg38 positions dropped)
**Combined table:** `multibiobank_credible_sets_combined.parquet` — **10,072,543 credible-set rows**, 5,080,965 distinct variants, **48,411 at PIP > 0.9**

This inventory assembles in-sample-LD fine-mapped credible sets across biobanks and molecular-QTL studies, unified to a common hg38 schema for the V2G2T program. Every source is tagged with `assay_basis` (wgs/imputation), `ancestry`, `build`, and `data_layer` (GWAS / molQTL / metaboliteQTL) so provenance is never ambiguous.

## Per-source summary

| Source | Layer | Assay | Ancestry | Rows | Traits | Cred. sets | PIP>0.9 | Has r² | Has gene |
|---|---|---|---|---|---|---|---|---|---|
| AoU | GWAS | wgs | LAT+EUR/AFR/EUR/AFR+LAT+EUR/AFR+EUR/LAT | 276,767 | 14 | 154 | 1,158 | 0 | 0 |
| UKB_WGS_WGS | GWAS | wgs | EUR | 471,513 | 34 | 30,059 | 10,993 | 0 | 0 |
| UKB_WGS_TOPMed | GWAS | imputation | EUR | 442,027 | 34 | 29,092 | 9,728 | 0 | 0 |
| UKB_WGS_UK10K | GWAS | imputation | EUR | 458,647 | 34 | 27,568 | 8,857 | 0 | 0 |
| FinnGen | GWAS | imputation | FIN | 542,605 | 1,424 | 19,947 | 3,257 | 542,605 | 421,623 |
| BBJ_JCTF | molQTL | imputation | EAS | 7,837,680 | 19,685 | 21,565 | 5,725 | 0 | 7,837,680 |
| UKB_metabolomics | metaboliteQTL | imputation | multi_ancestry_UKB | 43,304 | 249 | 27,746 | 8,693 | 43,298 | 42,670 |

## Source provenance & scientific notes

**AoU MultiSuSiE** (Nat Genet 2026, Zenodo 14458399) — WGS, hg38, genuinely **multi-ancestry** (African 47,041 + Latino/Admixed-American 36,378 + European 115,620). MultiSuSiE fine-mapping on summary stats + in-sample LD (LDStore2). 14 quantitative traits. Per-ancestry effect sizes retained (`beta` = primary-ancestry effect). This is the only source that breaks the European monoculture with in-sample fine-mapping.

**UKB-WGS** (Nature s41586-025-09720-6) — hg38, EUR, 34 traits, SuSiE in-sample. Split by fine-mapping panel via the paper's `Set` column into three provenance sub-sources: **WGS** (471,513 rows — the genuine whole-genome set), **TOPMed** and **UK10K** (imputation panels). Only the `UKB_WGS_WGS` sub-source is `assay_basis=wgs`.

**FinnGen R13** (public bucket `finngen-public-data-r13`) — hg38, Finnish, 1,424 endpoints, SuSiE in-sample. **Newer than the R12 credible sets in Open Targets**, and uniquely carries **r²-to-lead** (`cs_avg_r2`, from the `.cred.summary` files — 100% coverage) and gene annotation (`gene_most_severe`, 78% coverage) — fields Open Targets strips. Imputation-based (SISu panel), not WGS.

**BBJ / JCTF molQTL** (Wang et al. Nat Genet 2024, NBDC hum0343) — hg38, **East-Asian**, blood cis-eQTL (19,026 genes) + cis-pQTL (2,522 proteins), SuSiE with in-sample covariate-adjusted LD. The molQTL layer; "trait" = the regulated gene/protein. 7.84M rows is the full PIP>0.001 tail (thresholdable via the carried PIP).

**UKB metabolomics** (Nat Genet 2025, s41588-025-02355-3) — hg38, multi-ancestry UKB, **249 NMR metabolite traits**, SuSiE fine-mapped (42,498 SuSiE / 824 Wakefield). Carries r²-to-lead and lead-gene annotation. Upgrades the earlier keyword-matched metabolite-GWAS proxy to real fine-mapped credible sets.

## Layer & assay composition

- **By data layer:** GWAS 2,191,559 · molQTL 7,837,680 · metaboliteQTL 43,304
- **By assay basis:** WGS 748,280 · imputation 9,324,263

## QC checks passed

- **Build:** 100% hg38 across all rows. Chromosomes uniform `1-22, X` (FinnGen/metabolomics `23` normalized to `X`).
- **Nulls:** zero nulls in `variant_hg38`, `chromosome`, `PIP`, `trait`, `source`. 18 metabolomics rows with blank hg38 position dropped (unusable for coordinate join).
- **Cross-biobank GWAS-variant overlap** (biological sanity): all-3-way (AoU∩FinnGen∩UKB-WGS) = 9,617 variants; FinnGen∩UKB-WGS = 63,801 (both European-heavy, expected highest); AoU∩FinnGen = 18,008, AoU∩UKB-WGS = 24,802 (AoU multi-ancestry, expected lower overlap). Overlaps are non-trivial but far from identical — sources are complementary, not redundant.

## Companion G2T burden resources (significant hits only — NOT part of V2G2T)

- **AoU 392k All-by-All burden** (`aou_burden_g2t_sig.parquet`) — 1,580 gene-phenotype pLoF burden hits, 235 genes, 327 phenotypes, hg38 (medRxiv 2026, Supplementary Data 3/5/7: 216 pLoF + 44 cross-domain + 193 novel).
- **FinnGen R13 LoF burden** (`finngen_lof_burden_sig.parquet`) — 13,968 gene-trait burden hits, 3,984 genes, 304 genome-wide significant (mlogp>6).

These are held as trait/disease-association (G2T) validation resources per the acquisition scope, kept distinct from the variant-level credible-set (V2G2T) layer.

## Excluded (with reason)

- **Multi-cohort pQTL** (Cell 2026, 10.1016/j.cell.2026.03.049, 24,738 CS) — 38-study meta-analysis using reference-panel LD, not in-sample. Fails the in-sample gate; meta fine-mapping is miscalibrated at single-variant resolution.
- **Kanai 2021 cross-biobank fine-mapping** — legacy imputation-panel-based; dropped per the user's own analysis (noisy 94-trait fine-mapping).
- **FinnGen native pairwise LD matrix** — only the SISu v4.2 reference panel (n=3,775) is published (FinnGen says not for fine-mapping); the r²-to-lead summary we DO carry rides in the credible-set download.

## Staging

All outputs at `gs://claude_hackathon/multibiobank/20260711/`:
- `multibiobank_credible_sets_combined.parquet` (combined)
- `persource/` (7 per-source parquets)
- `aou_burden_g2t_sig.parquet`, `finngen_lof_burden_sig.parquet` (G2T resources)
- `qc/multibiobank_qc.json`, `qc/qc_pipdist.parquet`

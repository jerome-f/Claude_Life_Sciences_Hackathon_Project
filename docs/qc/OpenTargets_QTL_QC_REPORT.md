# Open Targets Molecular-QTL — Extraction & QC Report

**Date:** 2026-07-11
**Source:** Open Targets Platform release **26.06**, `output/credible_set/` (200 parquet parts, 3,518,871 total credible sets) and `output/study/` (requester-pays bucket, billed to project `flash-hour-452305-m7`).
**Compute:** GCP genomics-vm (e2-standard-16, 64 GB), duckdb + gcsfs.
**Staged to:** `gs://claude_hackathon/opentargets/20260711/` (processed / qc).

## 1. Definition of "high confidence"

The `confidence` field in credible_set is not a High/Med/Low label — it encodes the fine-mapping method and LD source. Observed vocabulary and quality order:

1. **SuSiE fine-mapped credible set with in-sample LD** ← selected (highest quality; matches the project's stated in-sample-LD preference)
2. SuSiE fine-mapped credible set with out-of-sample LD
3. PICS fine-mapped credible set extracted from summary statistics
4. PICS fine-mapped credible set based on reported top hit

**Filter applied:** `studyType ∈ {eqtl, pqtl, sqtl, tuqtl, sceqtl}` AND `confidence = "SuSiE fine-mapped credible set with in-sample LD"` (strict, per user decision).

## 2. Molecular-QTL result

**Total high-confidence molecular-QTL credible sets: 2,012,086.**

| studyType | Credible sets | Studies | Lead variants |
|---|---|---|---|
| eqtl | 1,349,418 | 1,210,080 | 560,428 |
| tuqtl (transcript usage) | 384,849 | 353,257 | 175,884 |
| sqtl (splicing) | 223,500 | 209,582 | 114,601 |
| sceqtl (single-cell eQTL) | 52,738 | 48,717 | 32,865 |
| pqtl (protein) | 1,581 | 789 | 1,528 |

### QC facts
- **All 2,012,086 sets are cis** (`isTransQtl = False` for every row) — in-sample SuSiE molecular-QTL credible sets in 26.06 are cis-only.
- **Zero quality-control flags** on any retained set (`qualityControls` empty for all).
- **`purityMeanR2` / `purityMinR2` are entirely null** for QTL credible sets (populated only for GWAS-type sets) — so purity-based filtering is not available for this layer; carried as a known gap.
- **Credible-set size** (variants): 5th=1, 25th=3, median=9, 75th=27, 95th=93.
- **Max lead-variant PIP** (posterior probability): 5th=0.036, 25th=0.121, median=0.301, 75th=0.698, 95th=0.950. Roughly a quarter of sets have a lead PIP > 0.70 — the confidently fine-mapped tail.

### The pQTL caveat (explicit)
Strict in-sample SuSiE keeps only **1,581 pQTL credible sets** — 95% of pQTL (32,137 out-of-sample SuSiE sets) is excluded by the uniform in-sample standard. Per user decision, the uniform standard was held rather than loosened for the protein channel. **The protein layer is effectively minimal in this extract.** If pQTL coverage is needed later, re-run allowing out-of-sample SuSiE for `studyType='pqtl'` only.

## 3. Metabolite GWAS (proxy — NOT molecular QTL)

**Finding: there is no metabolite-QTL studyType in OT credible_set 26.06.** Molecular QTLs are e/p/s/tu/sceqtl only. Metabolite signals exist as GWAS traits, not molecular QTLs. As requested, metabolite-level GWAS were hunted as a proxy and are clearly labelled as such.

**Method:** among `studyType='gwas'` studies (from `output/study/`), trait-keyword match on `traitFromSource` (metabolite / lipid / fatty acid / amino acid / cholesterol / lipoprotein / creatinine / urate / concentration-of / etc.), then linked to their SuSiE credible sets.

| Metric | Value |
|---|---|
| Metabolite-GWAS studies (keyword-matched) | 33,141 |
| Metabolite-GWAS SuSiE credible sets | 144,418 |
| — of which in-sample SuSiE | 238 |
| — of which out-of-sample SuSiE | 144,180 |

Top traits: metabolite levels, triglycerides, total/HDL/LDL cholesterol, fatty-acid ratios, creatinine, glycine, glucose, tyrosine, glutamine, apolipoprotein B. These are metabolite *concentrations measured as GWAS phenotypes* — a legitimate metabolite-genetics resource, but distinct from molecular QTLs (which measure a molecular readout in cis). **Do not merge these into the molecular-QTL layer.**

## 4. Outputs

| File | Rows | Description |
|---|---|---|
| `ot_molecular_qtl_highconf.parquet` | 2,012,086 | High-confidence molecular-QTL credible sets (e/p/s/tu/sceqtl, in-sample SuSiE) + derived credible_set_size, max_pip, n_in_95/99 |
| `ot_metabolite_gwas_studies.parquet` | 33,141 | Metabolite-GWAS study records (trait-keyword matched) |
| `ot_metabolite_gwas_crediblesets.parquet` | 144,418 | Their SuSiE credible sets (proxy set; clearly NOT molecular QTL) |

All staged under `gs://claude_hackathon/opentargets/20260711/processed/`, QC under `.../qc/`.

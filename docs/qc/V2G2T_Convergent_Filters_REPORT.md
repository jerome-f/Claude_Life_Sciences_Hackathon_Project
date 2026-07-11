# GWAS Convergent-Evidence V2G2T Gold-Standard Filters

**Source:** Open Targets 26.06 (`gs://open-targets-data-releases/26.06/output`), in-sample SuSiE credible sets only.
**Purpose:** extract high-confidence variant->gene->trait (V2G2T) triplets where the causal gene is nailed by *orthogonal coding evidence*, so the non-coding GWAS variant inherits an assay-independent gene label. These are the strongest positives for training/validating a V2G model — they do not depend on the eQTL/correlation signal the model itself will use.

Both filters are restricted to `confidence = "SuSiE fine-mapped credible set with in-sample LD"` (the strict in-sample tier agreed for this project).

---

## Filter 1 — LD-linked non-coding + pLoF credible sets

**Intended logic (user):** a GWAS locus with two credible sets, one prioritising a non-coding variant and one a pLoF variant in the proximal gene, both variants in high LD -> the non-coding variant->gene->trait is a valid V2G2T.

**What OT 26.06 actually supports — an important constraint.** OT ships **no variant-level LD r2 for in-sample SuSiE** credible sets (`ldSet` is empty, `locus.r2Overall` is null — expected, since in-sample LD uses the study's own genotypes, which OT cannot redistribute). The strict two-credible-set-with-r2 construction therefore returns **0** — and this is mechanistically correct, not a bug: SuSiE splits two credible sets *precisely because* it resolved them as distinct signals, i.e. NOT in high LD. Requiring two separate CS to be in high mutual LD fights the fine-mapper.

**The in-sample-faithful equivalent that IS populated:** a single credible set whose non-coding **lead** and a **pLoF variant** are *co-members of the same 95% credible set* — meaning the fine-mapper could not statistically separate them in-sample (maximal LD by construction), while the pLoF still assigns the causal gene. LD proxy = 95%-CS co-membership, **bounded by physical distance** (a 95% set spanning megabases is loose fine-mapping, not real LD).

**Result:**

| Tier (lead<->pLoF distance) | Pairs | Distinct (variant, gene) |
|---|---|---|
| **<=250 kb (high-confidence)** | **86** | **53** |
| 250-500 kb (medium) | 39 | (74 cumulative) |
| >500 kb (wide, flagged) | 27 | — |

**High-confidence examples (distance 1.5-12 kb, textbook causal assignments):**

| Trait | Non-coding lead | pLoF variant | Gene | Dist (bp) |
|---|---|---|---|---|
| Leiomyoma of uterus | 12_915604_G_A | 12_914052_G_T | **RAD52** | 1,552 |
| Venous thromboembolism | 9_133252218_C_G | 9_133255669_CG_C | **ABO** | 3,451 |
| Chronic diseases of tonsils and adenoids | 19_10279738_T_C | 19_10283576_C_CG | **ICAM1** | 3,838 |
| Von Willebrand disease | 9_133261662_G_A | 9_133257521_T_TC | **ABO** | 4,141 |
| Phlebitis and thrombophlebitis of lower extr | 9_133261662_G_A | 9_133257521_T_TC | **ABO** | 4,141 |
| Benign neoplasm: Peripheral nerves and auton | 22_20989430_G_A | 22_20993977_G_A | **LZTR1** | 4,547 |
| Benign neoplasm of other and unspecified sit | 22_20989430_G_A | 22_20993977_G_A | **LZTR1** | 4,547 |
| Autoimmune diseases excluding thyroid diseas | 9_133266092_G_A | 9_133257521_T_TC | **ABO** | 8,571 |

Recognisable hits: **ABO** (venous thromboembolism), **CHEK2** / **LZTR1** (cancers), **ICAM1**, **RAD52**, **HTRA1** (macular degeneration).

**Output:** `v2g2t_filter1_ld_noncoding_plof.parquet` (152 rows, all <=1 Mb; `ld_proxy_tier` column for filtering). pLoF calls use LOFTEE high-confidence (`lofteePrediction='HC'`) plus SO pLoF terms (stop-gained, frameshift, splice donor/acceptor, start/stop-lost).

---

## Filter 2 — burden + GWAS + molecular colocalization convergence

**Logic (user):** a burden test identifies a pLoF gene G for trait D; a GWAS credible set for D is proximal to G; and molecular colocalization shows an eQTL for G at the locus -> strong likelihood the GWAS variant points to the burden gene.

**Construction (three orthogonal lines):**
1. **Burden:** `evidence_gene_burden`, `directionOnTarget='LoF'` -> 8,454 gene-disease pairs (2,091 genes; AstraZeneca PheWAS + Genebass/UKB 450k).
2. **Proximity:** an in-sample SuSiE GWAS credible set for a matching disease (shared EFO/MONDO ID) within +/-500 kb of the gene body -> 1,221 candidate rows / 306 genes.
3. **Colocalization:** that GWAS CS colocalises (**H4 >= 0.8**) with an in-sample e/p/sQTL credible set for the **same** gene.

**Result: 166 convergent triplets, 133 distinct (variant, gene), 84 genes.** QTL support: 146 eQTL, 56 sQTL, 11 pQTL (many triplets carry multiple QTL types; 46 have >=2). Every triplet carries a colocalization-derived **effect direction** (`beta_sign`): 84 concordant / 82 discordant between QTL and GWAS effect — an even split, which is expected biology (higher expression is protective at some genes, risk-increasing at others), so the sign is an informative feature, not a bias.

**Top examples (H4=1.0):**

| Gene | Trait | GWAS lead | QTL types | beta sign |
|---|---|---|---|---|
| **FGFR3** | Height, inverse-rank normalized | 4_1804792_C_T | eqtl,sqtl | + |
| **APOB** | Disorders of lipoprotein metabolism and  | 2_21043902_GGCAGCGCCA_G | eqtl | - |
| **FES** | Hypertension | 15_90884462_G_A | eqtl | - |
| **FES** | Hypertension, essential | 15_90884462_G_A | eqtl | - |
| **GRN** | Dementia | 17_44352876_C_T | pqtl,eqtl | - |
| **TP53** | Malignant neoplasm, excluding all cancer | 17_7668434_T_G | eqtl,sqtl | + |
| **SCMH1** | Height, inverse-rank normalized | 1_41078607_G_T | sqtl,eqtl | + |
| **PTOV1** | Body-mass index, inverse-rank normalized | 19_49860328_CAG_C | eqtl | - |
| **CFI** | Age-related macular degeneration (whethe | 4_109740713_T_A | eqtl | - |
| **TP53** | Myeloproliferative diseases (CML exclude | 17_7668434_T_G | sqtl,eqtl | + |

Recognisable hits: **GRN** (dementia; pQTL+eQTL), **APOB**/**LDLR** (lipids), **TP53** (cancers), **CFI**/**CFH** (macular degeneration), **SLC2A9**/**ABCG2** (urate), **TERT** (leiomyoma), **PROC** (coagulation).

**Output:** `v2g2t_filter2_burden_gwas_coloc.parquet` (166 rows; `best_h4`, `qtl_types`, `beta_sign`, `burden_score`, `n_coloc_qtl_cs`).

---

## Why these matter for the model

Both filters yield **assay-independent positives**: the causal gene is fixed by coding/burden evidence, not by the eQTL-correlation signal a V2G model trains on. Per the held-out design (partition by assay class to avoid eQTL circularity), these are ideal for a *coding-anchored* validation slice — a non-coding variant with a gold gene label that no correlation-based method can trivially recover. Filter 1 gives 53 high-confidence non-coding->gene pairs; Filter 2 gives 133 convergent triplets with direction.

**Caveats:**
- Filter 1 LD is a **co-membership + distance proxy**, not r2 (OT limitation for in-sample SuSiE). The <=250 kb tier is defensible; treat 250-500 kb as medium and >500 kb as unreliable.
- Filter 2 disease matching is on shared EFO/MONDO IDs; trait granularity differs between burden (often broad) and GWAS (often specific), so some matches are at the disease-family level.
- Both are European-biased (UK Biobank / FinnGen dominate the in-sample SuSiE GWAS layer).

**Staged:** `gs://claude_hackathon/opentargets/20260711/v2g2t_goldstandard/` (parquets + `qc/` summaries). Scripts: `filter1.py`, `filter2.py`, `dedup2.py`.

# Non-European molecular-omics QTLs + Validated Variant–Gene–Trait Triplets

**Date:** 2026-07-11
**Purpose:** (1) expand the QTL layer beyond UKB-PPP/European sources with non-European expression/proteomics/metabolomics credible sets; (2) inventory experimentally *solved* variant→gene→trait triplets to serve as a curated truth set beyond statistical fine-mapping.

---

## Part 1 — Non-European / non-UKBB molecular-QTL credible sets

### Structural finding (important for interpretation)
**Open Targets 26.06 does not populate ancestry for molecular QTLs.** Every eQTL/sQTL/tuQTL/sceQTL study has `ancestry = null` in the study table; only pQTL carries a label (European). So ancestry cannot be filtered inside OT — it must be curated at the project level from the source publications. Two consequences:
- The catalogue is overwhelmingly European. GTEx alone contributes **1,156,194 of 2,012,086** high-confidence QTL credible sets (57%), and GTEx is ~85% European.
- Non-European signal in OT comes from a handful of **pooled** immune/response cohorts. Their credible sets were fine-mapped on the *pooled* sample, so OT cannot resolve any individual credible set as ancestry-specific — they are "includes non-European samples", not "non-European credible sets".

### Extracted from OT (data in hand): 158,208 credible sets from non-European-inclusive projects
`ot_noneur_qtl_highconf.parquet` — high-confidence (in-sample SuSiE) credible sets from the 7 OT projects with established non-European content, ancestry-annotated at project level:

| Project | Ancestry | High-conf CS |
|---|---|---|
| Quach_2016 | African + European (monocytes) | 83,107 |
| GEUVADIS | European + African YRI | 30,098 |
| Nedelec_2016 | African + European (macrophages) | 24,079 |
| Nathan_2022 | Peruvian / Indigenous American (T cells) | 10,879 |
| iPSCORE | Multi-ethnic Californian iPSC | 6,776 |
| Perez_2022 | East Asian + European (SLE) | 1,835 |
| Randolph_2021 | African + European (flu) | 1,434 |

Nathan_2022 is the only substantially **Indigenous-American** QTL resource in OT; Perez_2022 the main **East Asian**. **No dedicated African-American, Japanese, or metabolite-QTL resource exists in OT 26.06.**

### Major non-European resources NOT in OT (inventory — sourcing targets)
The real non-European omics data lives outside OT. Full table in `noneuropean_omics_resource_inventory.csv`; the high-value, openly-sourceable ones:

| Resource | Ancestry | Omic | Access |
|---|---|---|---|
| **ImmuNexUT** (Ota 2021) | Japanese | eQTL, 28 immune cell types | open sumstats |
| **Japan Omics Browser** | Japanese | eQTL + pQTL fine-mapped + MPRA | open portal |
| **BioBank Japan** (Sakaue 2021) | Japanese | metabolite + biomarker GWAS, SuSiE | open sumstats (pheweb.jp) |
| **ARIC** (Zhang 2022) | African American + European American | pQTL, SomaScan 4,877 proteins | open sumstats |
| **Jackson Heart Study** | African American | pQTL, Olink ~3,000 | portal (indiv gated) |
| **deCODE** (Ferkingstad 2021) | Icelandic (distinct LD) | pQTL, SomaScan 4,907 | open sumstats |
| **Multi-ancestry brain pQTL** (2025) | AA + Hispanic + NHW | brain pQTL, MESuSiE fine-mapped | open (Nat Genet ST) |

**Metabolomics specifically:** there is no metabolite *QTL* in OT (confirmed earlier). The non-European metabolomics path is **BioBank Japan** (metabolite GWAS, SuSiE-fine-mapped, hg19, public) — the East-Asian analogue of the metabolite-GWAS layer already extracted for European cohorts. ARIC and deCODE also carry metabolomics arms.

**Recommendation:** the two highest-yield, fully-open sourcing targets that add a genuinely new ancestry axis are **ImmuNexUT (Japanese eQTL)** and **ARIC AA pQTL** — both parallel resources already in the pipeline (immune eQTL; plasma pQTL) but in East-Asian and African-American backgrounds. Each needs the same harmonization ST16 needed (build liftover + variant-ID normalization).

---

## Part 2 — Experimentally solved variant–gene / variant–gene–trait triplets

Statistical fine-mapping (SuSiE credible sets) gives *probable* causal variants; it does not establish the gene or the mechanism. These resources carry **experimental** validation (CRISPRi perturbation, MPRA allelic activity) — the truth set for auditing V2G direction and calibration. Full table in `validated_triplet_resource_inventory.csv`.

### Extracted (data in hand): the CRISPRi-trained E2G layer
OT ingests **ENCODE-rE2G** (Gschwind et al. 2023) — the field-standard enhancer→gene model trained on a gold-standard of **10,411 CRISPR-tested element–gene pairs (1,075 validated positives)** in K562.

`ot_encode_re2g_distal.parquet` — **10,370,321 distal enhancer→gene links** (18,575 genes, 322 biosamples):
- **QC caveat:** the OT E2G `score` is ≈1.0 for nearly all links (median 0.9997), so a score threshold is *not* selective. The meaningful filter is `intervalType`: of 37.1M links ≥0.8, **24.4M are promoter self-links** (distanceToTss≈0). The **10.37M distal** links (intergenic+genic, ≥1 kb from TSS; distance median 7.3 kb) are the enhancer signal comparable to DBNascent's bidirectional→gene pairs.
- This is the direct, CRISPRi-anchored comparator for the DBNascent signed E2G edge — and the substrate for the signed-direction audit.

### Inventory of validated triplet resources (sourcing targets)

| Resource | Layer | Evidence | N | Access |
|---|---|---|---|---|
| **ENCODE-rE2G CRISPRi benchmark** | E2G | CRISPRi element-gene | 10,411 tested / 1,075 pos | EngreitzLab/CRISPR_comparison GitHub repo (open); exact Synapse mirror ID unverified |
| **cS2G** (Gazal 2022) | V2G2T | calibrated SNP-gene-disease | 7,111 triplets, mean PIP 0.80 | open |
| **Psychiatric MPRA** (Cell 2025) | V2G2T | MPRA+eQTL+HiC+CRISPR | 8 disorders | open |
| **MPRAu 3'UTR** (Griesemer 2021) | V2G | MPRA 3'UTR causal variants | thousands, 2 CRISPR-confirmed | open |
| **SCZ neurogenesis MPRA** (McAfee 2022) | V2G | MPRA regulatory variants | 5,173 tested / 439 pos | open |
| **OT L2G gold standard** | V2G2T | curated causal genes | ~445 GSP | open (in OT) |

**Why this matters for the project:** the ENCODE-rE2G CRISPRi set is exactly test 5 flagged in the negatives-validation work — the orthogonal, *causal* (not co-expression, not statistical) truth set. It closes the loop: DBNascent gives the signed co-transcription label, the CRISPRi benchmark says whether that sign matches a real perturbation effect. cS2G and the MPRA sets extend the audit to the full V→G→T chain and to non-enhancer mechanisms (3'UTR/RNA-stability) the enhancer-only view misses.

---

## Outputs

| File | Rows | Description |
|---|---|---|
| `ot_noneur_qtl_highconf.parquet` | 158,208 | Non-European-inclusive OT QTL credible sets, project-level ancestry annotation |
| `ot_encode_re2g_distal.parquet` | 10,370,321 | CRISPRi-trained distal enhancer→gene links (validated E2G layer) |
| `ot_encode_re2g_highconf.parquet` | 37,110,868 | All ENCODE-rE2G links ≥0.8 (incl. promoter self-links) |
| `noneuropean_omics_resource_inventory.csv` | 15 | Non-European eQTL/pQTL/metabolomics resources + access status |
| `validated_triplet_resource_inventory.csv` | 8 | Experimentally validated V2G/V2G2T truth sets + access status |

Staged under `gs://claude_hackathon/opentargets/20260711/` (`noneuropean/`, `validated_e2g/`, `qc/`).

**Caveats.** (1) OT ancestry annotation for QTLs is literature-curated at project level, not per-credible-set; pooled cohorts are not ancestry-resolvable. (2) The inventoried non-OT resources are *sourcing targets*, not yet extracted — each needs build/ID harmonization. (3) ENCODE-rE2G score is non-selective; always filter by intervalType/distance, not score.

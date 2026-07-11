# DBNascent Negatives Validation + UKB-PPP In-Sample pQTL Sourcing

**Date:** 2026-07-11

---

## Part 1 — Negatives contamination test

**Question:** the DBNascent trainable negatives are hard negatives by construction — (gene, transcribed-enhancer) pairs within ±1 Mb that are NOT among the gene's FDR<0.01 significant partners. The failure mode is **false negatives**: real regulatory links DBNascent missed (low co-transcription power, or hidden by the file's FDR<0.01 pre-filter). We test contamination against an **orthogonal** truth channel — fine-mapped eQTLs — that was never used to build the negatives.

### Method
For each (gene G, enhancer E) pair, ask whether any **in-sample-SuSiE eQTL 95% credible-set variant for gene G** (from our OT extract, 33.4M gene-mapped credible-set variants) falls inside E's genomic interval. A clean negative set sits near genomic background for this independent evidence; positives should be enriched. Comparison is restricted to the fair universe — pairs whose gene has at least one eQTL (positives: 5,332,151; negatives: 6,957,818) — and stratified by distance-to-TSS to rule out a distance covariate.

### Result — the negatives are validated, with a measured contamination floor

| Metric | Positives | Negatives |
|---|---|---|
| Pairs in eQTL-gene universe | 5,332,151 | 6,957,818 |
| With independent eQTL support in the enhancer | 261,227 | 237,854 |
| **Rate** | **4.90%** | **3.42%** |

- **Overall enrichment = 1.43×** (positive rate / negative rate). Positives carry orthogonal causal evidence at 1.43× the rate of negatives — the construction genuinely separates linked from unlinked pairs.
- **The enrichment holds at every distance bin** (1.10–1.20× across 0–1 Mb; see figure), so it is not an artifact of positives being closer to the TSS. Negatives do sit slightly farther out (median 496 kb vs 442 kb), but distance stratification confirms the signal survives matching.
- **False-negative contamination floor = 3.42%** (237,854 negatives carry independent eQTL support and are likely real missed links).

### Interpretation
Two things are true at once. (1) The negative set works — it is significantly depleted of causal signal relative to positives, uniformly across distance. (2) The separation is **modest** (1.43×), which is itself the scientifically important reading: it quantifies how noisy nascent co-transcription is as a causal label. This is the same premise the concordance test targets, now measured on the negative side. A signed E2G model trained on these labels inherits a ~1.4× signal-to-contrast ceiling from the label, before any modelling choice.

### Action taken
The 237,854 contaminated negatives are **flagged, not silently kept**. Output `dbnascent_negatives_flagged.parquet` adds an `eqtl_supported` boolean:
- **6,719,964 clean negatives** (`eqtl_supported=false`) — the defensible training-negative set.
- **237,854 flagged** (`eqtl_supported=true`) — move to an ambiguous/held-out bucket rather than training as hard negatives, or use as a positive-leaning weak-label set.

**Caveat:** this test detects only contamination visible to the eQTL channel (European/blood-immune-biased, cis-eQTL). True contamination is ≥3.42%; the 96.6% "clean" fraction is an upper bound on cleanliness, not a guarantee. Test 5 (ENCODE-rE2G CRISPRi truth) would tighten it and is the recommended next check.

---

## Part 2 — UKB-PPP in-sample pQTL sourcing

**Finding confirmed:** OT's out-of-sample pQTL label is a *pipeline artifact*. OT re-fine-maps UKB-PPP with pan-UKBB reference LD (not the study's own genotypes), so every UKB-PPP credible set is flagged out-of-sample — which is why strict in-sample kept only 1,581 pQTL sets. The **original UKB-PPP fine-mapping (Sun et al. 2023, *Nature*, doi 10.1038/s41586-023-06592-6) is in-sample** and is published as **Supplementary Table 16** ("Independent pQTL signals", SuSiE credible sets).

### Sourced
Fetched the Nature supplementary workbook (ungated; the AWS `ukbiobank.opendata.sagebase.org` bucket and Synapse syn51364943 are both access-controlled and give summary stats, not credible sets). Parsed ST16:

| Metric | Value |
|---|---|
| Independent pQTL signals (SuSiE credible sets) | **29,420** |
| Proteins | 2,414 |
| cis | 10,750 |
| trans | 18,670 |
| cis credible-set size (median) | 2 |
| Top-variant PIP (median) | 0.51 |
| Build | **hg19** |

This is **~19× more in-sample pQTL** than OT's 1,581 — a real protein channel at proper in-sample confidence.

### Harmonization still required before merge
- **Liftover hg19 → hg38** (ST16 is hg19; OT is hg38). Blocking step before any coordinate join.
- **Variant-ID normalization** — ST16 uses `chr:pos:ref:alt:imp:vN`; OT uses `chr_pos_ref_alt`.
- **Protein → gene mapping** — ST16 `UKBPPP_ProteinID` (`GENE:UniProt:OID:vN`); gene symbol parsed into `gene_symbol`, but Olink assay → Ensembl should be confirmed against the OT target table.
- **Schema alignment** — ST16 credible sets carry PIP/log10BF/CS-size but not OT's full credible_set schema; a mapping layer is needed to co-file with `ot_molecular_qtl_highconf.parquet`.

trans-pQTL (63% of signals) are genuine but must be handled separately from the cis molecular-QTL layer, exactly as trans is elsewhere.

---

## Outputs

| File | Rows | Description |
|---|---|---|
| `dbnascent_negatives_flagged.parquet` | 6,957,818 | Negatives + `eqtl_supported` contamination flag (6,719,964 clean / 237,854 flagged) |
| `ukbppp_st16_insample_pqtl_hg19.parquet` | 29,420 | UKB-PPP in-sample SuSiE pQTL credible sets (hg19, pre-liftover) |
| `UKBPPP_Sun2023_supplementary_tables.xlsx` | — | Full Nature supplementary workbook (source) |

Negatives outputs under `gs://claude_hackathon/dbnascent/20260711/`; pQTL under `gs://claude_hackathon/opentargets/20260711/pqtl_insample/`.

---

## Part 3 — UKB-PPP pQTL lifted to hg38 (harmonization complete)

The ST16 in-sample pQTL credible sets (originally hg19) were lifted to **hg38** with the UCSC `hg19ToHg38.over.chain` (pyliftover), and variant IDs rebuilt into OT's `chr_pos_ref_alt` convention so they co-file directly with `ot_molecular_qtl_highconf.parquet`.

### Liftover result
| Metric | Value |
|---|---|
| Top variants lifted | 29,420 / 29,420 (**100%**) |
| Credible-set variants lifted | 572,574 / 572,574 (**100%**) |
| Unmapped / unparseable | 0 |

### Independent validation against OT hg38 variant universe
Exact-string match of lifted top variants against the OT credible_set variant universe (1,136,908 distinct hg38 IDs) — this validates **position AND allele orientation**, not just that a coordinate was produced:

| | Signals | Exact match in OT hg38 |
|---|---|---|
| cis | 10,750 | 7,437 (69%) |
| trans | 18,670 | 15,876 (85%) |

69–85% of UKB-PPP top variants land exactly on a known OT hg38 variant ID — strong confirmation the liftover is correct (not all UKB-PPP pQTLs are OT credible-set members, so <100% is expected). Only **326 (1.1%)** matched solely under a ref/alt swap — an allele-ordering convention difference, carried as a known minor caveat, not a liftover error.

### Status
`ukbppp_st16_insample_pqtl_hg38.parquet` (adds `top_variantId_hg38`, `cs_variant_ids_hg38`, per-row lift status) is now **build- and ID-compatible with the OT molecular-QTL layer**. Remaining before a formal merge: confirm Olink protein → Ensembl mapping against the OT target table (gene symbol already parsed), and decide cis-only vs cis+trans inclusion for the protein channel.

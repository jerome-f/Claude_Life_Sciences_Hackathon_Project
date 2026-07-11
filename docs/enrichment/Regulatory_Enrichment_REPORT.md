# Regulatory-Class Enrichment of Fine-Mapped GWAS Credible-Set Variants

**Question.** Before committing to element-class modeling heads (promoter / enhancer / silencer,
with CTCF handled by physics and splice handled by SpliceAI), how is each regulatory class
represented among fine-mapped GWAS variants, and which classes are actually *enriched* for
causal signal?

## Method

- **Substrate:** the GWAS layer of the multi-biobank credible-set set (AoU, UKB-WGS, FinnGen R13),
  2,191,559 credible-set member rows on hg38.
- **Sparse credible sets:** the 51,658 credible sets (of 106,547) that
  contain at least one member at **PIP ≥ 0.5** — i.e. sets where fine-mapping actually resolved a
  likely-causal variant.
- **Foreground:** 27,040 distinct variants at PIP ≥ 0.5 within those sets.
- **Background:** 182,739 distinct variants at PIP < 0.5 within the *same* sets.
  Using the low-PIP members of the same loci as background controls for locus ascertainment and LD
  structure without an external random-variant panel — the enrichment then reads as "as fine-mapping
  concentrates posterior on a variant, which regulatory classes does it concentrate *into*."
- **Annotation:** ENCODE SCREEN Registry of cCREs, v3 / ENCODE4 (hg38), 2,345,453
  elements, pulled via the UCSC REST API (`cCREregistry` track). Each variant is assigned **all**
  overlapping cCRE classes (multi-label: a variant in a promoter that also overlaps a CTCF element
  gets both). Fold-enrichment = (foreground fraction) / (background fraction) per class; significance
  by Fisher exact test; 95% CIs by the Katz binomial log-ratio approximation.

## Headline

- **38.7%** of foreground variants overlap a cCRE vs **26.9%**
  of background — fine-mapped variants land in annotated regulatory DNA markedly more often.

## Per-class enrichment

| cCRE class | element head | FG % | BG % | fold | 95% CI | Fisher p |
|---|---|---|---|---|---|---|
| Promoter | promoter | 2.53 | 0.84 | 3.002 | 2.746–3.282 | 1.2e-108 |
| Proximal enhancer | enhancer/silencer | 6.97 | 3.61 | 1.928 | 1.835–2.026 | 4.7e-129 |
| CA-H3K4me3 | ambiguous OC | 0.92 | 0.65 | 1.405 | 1.226–1.61 | 2.5e-06 |
| Distal enhancer | enhancer/silencer | 24.04 | 17.27 | 1.392 | 1.36–1.425 | 2.4e-150 |
| CA-TF | ambiguous OC | 0.31 | 0.26 | 1.175 | 0.933–1.481 | 1.7e-01 |
| CA-CTCF | loop anchor (CTCF) | 0.94 | 0.94 | 1.003 | 0.88–1.143 | 9.7e-01 |
| TF | ambiguous OC | 1.3 | 1.35 | 0.962 | 0.861–1.075 | 5.1e-01 |
| CA | ambiguous OC | 1.72 | 1.94 | 0.884 | 0.803–0.973 | 1.2e-02 |

## Reading the result

1. **Promoter — 3.0× (p≈1e-108).** The strongest enrichment, as expected: promoter-proximal causal
   variants are the easiest class and validate the design of a dedicated promoter head.
2. **Enhancer-like (proximal 1.9× + distal 1.4×, both p<1e-120).** The largest bucket by far —
   24% of foreground variants overlap a distal enhancer. This is the class the enhancer/silencer head
   targets, and it is where the modeling effort should concentrate. **Note:** cCREs cannot separate
   activating enhancers from silencers (both are "enhancer-like" by chromatin signature); that
   activating-vs-repressive sign is precisely the signed edge the model is being built to learn, so it
   is not annotatable upstream and is reported here as one combined bucket.
3. **CA-CTCF — 1.00× (ns).** Loop-anchor/insulator elements are **not** enriched at fine-mapped
   variants relative to background. This supports handling CTCF separately (a physics/loop-extrusion
   module) rather than folding it into the expression head — the causal-variant density simply is not
   concentrated there.
4. **Bare open-chromatin (CA, TF, CA-TF) — flat to depleted.** Accessibility without a promoter/enhancer
   signature carries little causal enrichment; these are the ambiguous residue, not a modeling target.

## Caveats

- **Cross-cell-type union.** The cCRE registry is a union across biosamples; it answers "is this a
  regulatory element in some context," not "in the trait-relevant tissue." Tissue-resolved ChromHMM
  would sharpen the enhancer call but needs trait→tissue mapping.
- **Silencers are invisible to annotation** (point 2) — the enhancer-like fold-enrichment bounds the
  *combined* activating+repressive bucket.
- **Scope = GWAS layer only.** Molecular-QTL (BBJ) and metabolite-QTL members were excluded; their lead
  variants concentrate in different (often coding/promoter-proximal) classes and would answer a
  different question.

## Files
- `ccre_enrichment_by_class.csv` / `.parquet` — the enrichment table above.
- `credible_set_variant_ccre_labels.parquet` — every foreground/background variant with its multi-label cCRE assignment.
- `ccre_enrichment_credible_sets.png` (+ `ccre_enrichment_legend.png`) — the figure.
- `enrichment_summary.json` — headline counts.
- `ccre_v3_hg38.parquet` — the cCRE registry snapshot used.

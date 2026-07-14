# GWAS Effector-Gene Benchmark — Summary

_Assembled 2026-07-13. Purpose: support the question "at GWAS loci where cis-eQTL
colocalization fails to nominate the effector gene, does a sequence-based signed E2G model
recover the known causal gene?" This is a **data-assembly** deliverable — no model was run._

## Headline numbers

| Set | n loci | COLOC_POSITIVE | COLOC_WRONG | COLOC_SILENT | Winnable (SILENT+WRONG) |
|-----|-------:|---------------:|------------:|-------------:|------------------------:|
| **All resolved** | **1,160** | **24.5%** (284) | **30.9%** (358) | **44.7%** (518) | **75.6%** |
| High-confidence effectors only | 611 | 34.4% (210) | 19.1% (117) | 46.5% (284) | 65.6% |

- **The coloc-silent fraction — the headline winnable gap — is 44.7% of resolved loci (45–47% robust across confidence tiers).** Adding coloc-wrong loci (where cis-eQTL colocalises with a *different* candidate gene than the curated effector), the subset where cis-eQTL coloc does NOT correctly name the effector is **75.6%** of all resolved loci, or **65.6%** restricting to high-confidence effectors.
- This sits squarely inside the field-wide eQTL–GWAS coloc miss rate (~60–80%) and independently corroborates the project thesis: even a mature, multi-tissue coloc resource pinpoints the curated effector at only ~1 in 4 loci (~1 in 3 for high-confidence loci). The remaining ~3/4 is the space a sequence-based E2G model is meant to recover.
- 128 curated (locus, gene) rows are **UNKNOWN** — the sentinel variant is not present in the Open Targets variant index (mostly older-curation SNVs not observed in any OT GWAS/QTL credible set). These are reported as UNKNOWN, never guessed.

## What a row is

One row per **(curated locus sentinel variant, curated effector gene)** pair. 1,189 unique sentinel
variants, 526 unique effector genes, 1,288 rows total (a handful of variants carry >1 curated
effector; a handful of effectors recur across traits). Both the effector gene and all candidate
genes are given as **HGNC symbols** (the model index key), with Ensembl gene IDs alongside.

## Coloc-status definitions (per locus)

Classification uses **Open Targets Platform release 26.06** colocalisation (`COLOC_PIP_ECAVIAR`,
h4 ≥ 0.8) — the same release family as this project's harvested FinnGen GWAS×QTL coloc table.
For each sentinel variant we retrieved every **cis-QTL credible set that contains the variant**
(eQTL/sQTL/tuQTL/edQTL = expression-QTL family; pQTL recorded separately), resolved its target
gene, and asked whether that QTL credible set colocalises (h4 ≥ 0.8) with **any GWAS credible set**.

- **COLOC_POSITIVE** — an **expression-QTL** for the *curated effector gene* colocalises with a GWAS signal at the locus (h4 ≥ 0.8). Column `effector_eqtl_coloc_h4` gives the value.
- **COLOC_WRONG** — no expression-QTL coloc for the effector, but an expression-QTL for a **different candidate gene** (within ±500 kb) does colocalise. `coloc_nominated_gene` lists the gene(s) coloc points to instead.
- **COLOC_SILENT** — no expression-QTL colocalisation (h4 ≥ 0.8) for *any* candidate gene at the locus.
- **UNKNOWN** — the sentinel variant is absent from the OT variant index; coloc status cannot be determined and is not guessed.

`effector_pqtl_coloc_h4` is reported separately: some coloc-silent-for-eQTL loci (e.g. APOE) do
colocalise at the **protein** level, which is directly relevant to the project's cis-pQTL bridge.

## Candidate genes per locus

`candidate_genes` = all **protein-coding** genes whose gene body falls within **±500 kb** of the
sentinel (GRCh38), from Ensembl REST `overlap/region` (`biotype=protein_coding`). Median 12
candidates/locus (range 0–71). The curated effector falls inside this window for **98.9%** of
(locus, gene) pairs; the 1.1% that fall outside are flagged `effector_outside_500kb` in `notes`
(known long-range cases, e.g. the FTO-intronic obesity variant acting on IRX3/IRX5 ~1 Mb away).

## Sources

1. **Effector-gene gold standard** — Open Targets Genetics L2G gold-standard training set
   (Mountjoy et al. 2021, *Nat Genet* 53:1527–1533, doi 10.1038/s41588-021-00945-5), file
   `gwas_gold_standards.191108.tsv` from github.com/opentargets/genetics-gold-standards.
   2,435 curated rows; sentinel variant (GRCh38 pos + alleles + rsID where available), Ensembl
   effector gene, confidence tier, trait. Provenance within it: ChEMBL drug-target evidence,
   T2D Knowledge Portal effector genes, ProGeM (metabolite/mQTL), and expert curation
   (Fauman/Ghoussaini/Mountjoy).
2. **Colocalisation** — Open Targets Platform GraphQL API, release **26.06**
   (`api.platform.opentargets.org`), `COLOC_PIP_ECAVIAR` h4 ≥ 0.8, cis-QTL credible sets from
   GTEx v8, EBI eQTL Catalogue, UKB-PPP/plasma pQTL, etc. (OT's QTL ingest).
3. **Candidate genes** — Ensembl REST (GRCh38), protein-coding gene overlap ±500 kb.
4. Gene-symbol resolution — Ensembl REST `lookup/id` (524/526 gold ENSGs resolved to current
   symbols; 2 retired IDs kept as ENSG).

## Caveats (read before using)

- **"Coloc-silent" means no PP.H4 ≥ 0.8 expression-QTL hit in OT 26.06 — not proof of no eQTL.**
  It can reflect incomplete QTL tissue/cell-type coverage, the causal cell type not being assayed,
  or the effector acting in a context GTEx/eQTL-Catalogue does not capture. The silent fraction is
  the *opportunity* for E2G, bounded by current QTL panels, not a claim of biological absence.
- **Coloc source is OT 26.06, not this project's harvested FinnGen table.** The harvested FinnGen
  GWAS×QTL parquet is keyed by opaque studyLocusId hashes (no gene/variant columns) and is
  FinnGen-ancestry-specific, whereas the gold-standard traits are predominantly European
  UKB/GWAS-Catalog. OT 26.06 (same release family, gene-resolved, pan-resource) is the defensible
  instrument for a per-locus gene-level coloc verdict. Both draw on the same underlying QTL panels.
- **Effector-gene curation is not infallible.** The L2G gold standard is enriched for drug-target
  loci (ChEMBL, ~55% of rows, mostly Low/Medium confidence) and metabolite/T2D loci; it is not a
  random sample of the genome. Confidence tiers are provided (`effector_confidence`); the
  high-confidence subset (n=611) is the most defensible core.
- **Ancestry** — gold standard and OT GWAS credible sets are predominantly European.
- **80 loci are `cs_truncated_500`** — extremely pleiotropic variants where the credible-set scan
  was capped at 500 sets. QTL credible sets sort early, so coloc capture is complete in spot
  checks (e.g. SORT1 correctly POSITIVE), but the flag is retained for transparency.
- **UNKNOWN loci (128)** must be excluded from any accuracy denominator.

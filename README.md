# Signed Enhancer→Gene Mapping for GWAS Effector-Gene Identification

Translational human genetics project (Claude Life Sciences Hackathon, 2026).
Maps GWAS-associated non-coding variants to their causal **effector gene** and
the **direction** of the regulatory effect, from DNA sequence.

## The problem

Genome-wide association studies flag thousands of trait-linked loci, but the
causal *gene* usually stays hidden — the variant sits in non-coding regulatory
DNA, often far from the gene it controls. The field's default tool, cis-eQTL
colocalization, resolves this poorly, and no method reports the *direction*
(activation vs. repression) of a regulatory effect.

## What's here

A two-phase project:

1. **Data + fine-mapping phase** (`scripts/data_processing`, `scripts/multibiobank`,
   `scripts/v2g2t`, `docs/qc`, `docs/inventory`, `docs/enrichment`): harvesting and
   QC of DBNascent nascent-transcription E2G pairs, Open Targets / multi-biobank
   fine-mapped credible sets, colocalization results, and a convergent-evidence
   gold standard.
2. **Modeling phase** (`scripts/model`, `scripts/scoring`, `results`, `docs/model`):
   a signed enhancer→gene (E2G) model (**CATE** — Convolution-Attention Two-tower
   Encoder) over frozen Nucleotide-Transformer-v2 embeddings, and a locus-level
   utility benchmark against cis-eQTL colocalization.

## Key findings

On a benchmark of **1,160 GWAS loci** with curated effector genes
(Open Targets L2G × colocalization h4≥0.8):

- cis-eQTL colocalization names the correct effector at only **24.5%** of loci;
  it points to the **wrong** gene at 30.9% and is **silent** at 44.7% —
  **75.6% not correctly resolved**.
- The E2G model does not beat a nearest-gene distance baseline outright, but is
  **complementary**: it recovers **39 effector genes distance misses** (18 at
  coloc-silent loci), lifting union recovery from 47.7% → 51.1%.
- Uniquely, it assigns a **direction of effect** (held-out repressive-sign
  AUROC ≈ 0.61) that neither colocalization nor distance provides.
- Held-out enhancer→gene link recovery (7,666 genes): p@1 0.729 (model) vs
  0.719 (distance); LINK AUROC 0.611.

See `docs/model/project_description.md` for the narrative summary and
`docs/model/e2g_showcase.html` for the results webpage.

## Repository layout

```
scripts/
  data_processing/   DBNascent + Open Targets extraction & QC
  multibiobank/      AoU / UKB-WGS / FinnGen / BBJ credible-set pulls
  v2g2t/             convergent-evidence gold-standard filters
  model/             CATE model, training, per-token embedding, v5 two-tower
  scoring/           held-out scoring + locus-level effector nomination
  analysis/          concordance & replication tests
results/             benchmark CSV, locus nominations, utility metrics (JSON)
docs/
  model/             project description + results webpage
  figures/           architecture diagram + results figures
  qc/  inventory/  enrichment/    phase-1 reports
```

## License

Released under the [MIT License](LICENSE). Copyright (c) 2026 Jerome Irudayanathan.

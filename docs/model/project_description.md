# A signed sequence-based enhancer→gene model for GWAS effector-gene identification

**What I built and investigated.** GWAS have catalogued tens of thousands of
trait-associated loci, but naming the causal *effector gene* — usually acting
through a non-coding regulatory variant far from the gene — remains unsolved at
scale. The field's default evidence, cis-QTL colocalization, has a ceiling that
is larger than commonly stated. I set out to (i) measure that gap rigorously and
(ii) build a sequence-based line of evidence that adds what colocalization
cannot: a *directional* enhancer→gene link available where colocalization is
silent. Over the past week this ran end to end: I assembled a harmonized
five-source, 10M-row multi-biobank credible-set substrate (All-of-Us, UKB-WGS,
FinnGen R13, BBJ molecular-QTL, UKB metabolomics); harvested published
colocalization results from UKB-PPP, FinnGen R13 and BBJ; derived 6.7M
FDR-controlled *signed* enhancer→gene pairs from the DBNascent nascent-RNA atlas
(validated to Pearson = 1.0 against published pairs, plus a clean-room mouse
reimplementation); embedded ~1.2M regulatory windows through a frozen Nucleotide
Transformer v2; and trained a family of signed-E2G models culminating in CATE, a
convolution-attention two-tower encoder.

**What I found.** Three results. (1) The gap is real and large: on my own
substrate only 1–2% of confident QTL signals share a causal lead with a GWAS
signal despite ~90% being feasibly close, and on a 1,160-locus benchmark
cis-eQTL colocalization names the curated effector at only 24.5% of loci
(44.7% silent, 30.9% wrong — 75.6% unresolved), corroborated by three flagship
resources. (2) A DBNascent-derived sign predicts causal *direction* (Spearman
ρ = −0.224, p = 5×10⁻⁹; blood×K562 odds ratio 5.4) — weak but robust, and
unique among E2G methods. (3) The model is *complementary* to distance: it
recovers 39 effector genes distance misses (18 at coloc-silent loci, including
GLP1R, MYC, SNCA, VDR), and a methodological finding — link AUROC plateaus at
~0.61 across gradient-boosted trees, a pooled MLP, and cross-attention — localizes
the ceiling to pooling away sequence positions, motivating the per-token CATE model.

**Why it matters.** Effector-gene assignment is the bottleneck between a GWAS hit
and a drug target: get the gene wrong and every downstream decision is wrong too.
Quantifying exactly where the standard tool fails, and adding a sequence-only,
directional signal available at the silent majority of loci, gives a more honest
and more complete locus-to-gene map.

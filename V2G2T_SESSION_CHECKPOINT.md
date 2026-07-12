# V2G2T Project — Session Checkpoint
**Generated:** 2026-07-11 | **Project:** proj_557c1e6c8c4d | **Resume from:** any new session (read this file first)

This is a full-state handoff. Read top-to-bottom; every artifact ID is verbatim and current.

---

## 1. PROJECT SCOPE

Building a **V2G2T** (variant → gene → trait) pipeline that assigns a **calibrated probability at
every hop** and treats **signed direction** as first-class — not just "variant regulates gene" but
"variant *up*-regulates gene, which *raises* the trait." The distinguishing bet vs. Open Targets /
ENCODE-rE2G / AlphaGenome / cS2G is the **signed edge** plus a multiplicative **V2E × E2G × G2T**
sign-chain with explicit per-hop uncertainty and abstention.

**5-step pipeline:** V2E (variant→element activity, signed) → E2G (enhancer→gene, signed) →
V2G = V2E×E2G → G2T (gene→trait, burden+GWAS) → V2G2T (full signed chain).

**Novelty thesis (verified mid-2026):** the *unsigned* E2G edge is already occupied (ENCODE-rE2G
>13M links, in Open Targets as L2G features; scE2G). Direction is largely held by AlphaGenome.
Calibrated composition pre-exists as cS2G. The genuinely unoccupied slice = **nascent
co-transcription as a SUPERVISED label to predict the SIGNED E2G edge from sequence**, plus the
explicit sign-chain with abstention.

**User:** human-genetics scientist (startups + consulting), translational focus. Expert co-designer —
consistently catches methodological subtleties (z-score magnitude, in-sample LD requirements). Being
deliberately exhaustive on data to prove the V2G2T strategy is worth investing in (vs. a prior 2-year
manual-curation timeline). Prefers polars/duckdb over pandas (house skill). House plot style =
inclusive-plot-style skill (separate legend images, never inline ax.legend()).

---

## 2. WHERE WE ARE — STATUS

| Phase | Status |
|---|---|
| Data acquisition (credible sets + QTL + burden) | **COMPLETE** |
| Concordance design-gate (does DBNascent signed label track truth?) | **COMPLETE — qualified pass** |
| Two GWAS convergent V2G2T filters (gold loci) | **COMPLETE** |
| Deduplicated data inventory PDF | **COMPLETE** |
| Regulatory-class enrichment (cCRE, pre-modeling diagnostic) | **COMPLETE (this session)** |
| genomics-vm rebuild (fixed stale-mount blocker) | **COMPLETE (this session)** |
| **Model architecture — build the signed E2G model** | **NEXT — design agreed, not yet built** |

---

## 3. THE CONCORDANCE GATE (passed, qualified)

Does DBNascent's nascent-co-transcription signed E2G label track a real regulatory truth?
- **Primary: Spearman rho = -0.224 (p=4.8e-9, 666 matched pairs)** vs ENCODE-rE2G CRISPRi benchmark.
  Sign convention: enhancer regulatory sign = -sign(CRISPRi EffectSize). Robust across all cell types
  + confound controls (H3K27ac, well-powered, direct-effect). Blood×K562 rho=-0.286 (n=326).
- **Magnitude:** poor in raw units (different readouts) but the user's z-score insight rescued it —
  well-measured blood×K562 |signed-t| vs |EffectSize| rho=0.324 (p=3.1e-5). Poor magnitude was
  largely a measurement-noise artifact.
- **DC-TAP replication:** directionally consistent (all cell-matched rho<0) but n.s. — DC-TAP is
  deeply-powered-but-narrow, only ~28 usable extra pairs; full pool n=694 rho=-0.229 (unchanged).
- **VERDICT: qualified pass — the label tracks DIRECTION, weakly but robustly.** Three caveats bound
  the model design (see §6).

---

## 4. DATA SUBSTRATE (all hg38, complete)

### Combined credible-set table — THE substrate
`multibiobank_credible_sets_combined.parquet` — **10,072,543 rows, 5,080,965 distinct variants,
48,411 at PIP>0.9.** 17-col unified schema: variant_hg38, chromosome, position, ref, alt, trait,
study, PIP, cs_id, cs_membership, r2_to_lead, ancestry, assay_basis, build, source, beta, gene,
data_layer.
- artifact_id=`887360b8-75e4-453a-90dc-3c3d4bfc50a7` version=`3e94a838-d003-408a-b3e6-39fe941b54ed` (CHECKPOINT)
- Layers: molQTL 7,837,680 / GWAS 2,191,559 / metaboliteQTL 43,322.
- GWAS sources: FinnGen 542,605 / UKB_WGS_WGS 471,513 / UKB_WGS_UK10K 458,647 / UKB_WGS_TOPMed 442,027 / AoU 276,767.
- Ancestry of 10.07M rows: East Asian 7,837,680 (78%, BBJ molQTL artifact — flag in any ancestry plot) / European 1,444,626 / Finnish 542,605 / multi/non-Eur 247,632.

### Per-source parquets (slices of combined — do NOT sum with combined)
| file | rows | ancestry/assay | artifact / version |
|---|---|---|---|
| aou_credible_sets.parquet | 276,767 | multi-anc WGS (MultiSuSiE) | `46a4be96-e84d-441f-855a-32db19d72338` / `1a2579d8-55c5-47fc-9dde-d9cc48b1036c` |
| ukbwgs_credible_sets.parquet | 1,372,187 | EUR WGS+imputation | `33f5de49-689e-490c-bff8-1ff8e830f1a7` / `caeb64c2-4a7f-4d0e-a4a7-59c7c3ff8162` |
| finngen_credible_sets.parquet | 542,605 | FIN, has r2_to_lead+gene | `f079f164-3bdf-4d75-9747-c97adddc202a` / `97fad35c-1b49-4687-9128-6de330d91628` |
| bbj_molqtl_credible_sets.parquet | 7,837,680 | EAS in-sample molQTL | `8edaeba8-e1d8-4c9f-a66d-74e186633596` / `bef34302-f2f6-4f5e-8357-b869d93ab95a` |
| ukb_metabolomics_credible_sets.parquet | 43,322 | 249 NMR traits | `bca8223d-e973-4af6-9c3f-9c118b434d82` / `ece7dcbc-f55c-4228-b83d-ea547300313f` |

### G2T burden companions (significant hits only, NOT part of V2G2T substrate)
- aou_burden_g2t_sig.parquet (1,580 hits) — `e5377bfe-c7d7-4434-8a2b-ac7882bcbd88` / `67eb9025-4734-453a-9195-9e123a34a71e`
- finngen_lof_burden_sig.parquet (13,968 hits) — `19f11f70-f5b5-4f8f-b9eb-6ca8b3265fe1` / `d026c880-9f26-4ab9-81b7-9ea42a6425a5`

### DROPPED / EXCLUDED (do not re-add)
- **Kanai** UKB 94-trait fine-mapping — legacy imputation panel, user analyzed it, too noisy.
- **Cell meta-pQTL (2026)** — 38-study meta-analysis on reference LD, fails in-sample gate.
- OT metabolite-keyword proxy — superseded by real UKB NMR SuSiE sets.

### QC
`Multibiobank_Credible_Sets_QC_REPORT.md` = `980aa648-8d3e-47c2-be82-876f0978d099` / `b82997c5-754d-4c68-9903-a3c61598e077`;
overview png `448c70f3-6288-4f88-8f3e-d12329482646` / `2fa95b89-e1cb-4650-aab3-848002841d03`; mb_qc.json `ffe5bfe4-...`.

---

## 5. GOLD V2G2T LOCI (two convergent filters, complete)

- **Filter 1** (LD non-coding + pLoF co-membership): 86 pairs ≤250kb / 125 ≤500kb. Genes ABO/CHEK2/LZTR1/RAD52/ICAM1.
  `v2g2t_filter1_ld_noncoding_plof.parquet` = `eb209077-e63a-46fd-88cf-a0145106f59e` / `9f4686d7-47fd-43de-98c5-a71f58b0c4bb`
- **Filter 2** (burden + GWAS + eQTL coloc): 166 triplets, 133 distinct variant-gene, 84 genes. GRN/APOB/LDLR/TP53/CFH.
  `v2g2t_filter2_burden_gwas_coloc.parquet` = `c5b4bd68-5123-4111-8ce0-f9b5c7b2e934` / `5662671e-4c35-47c8-b9bc-4d0c7856b009`
- Report `V2G2T_Convergent_Filters_REPORT.md` = `048e7164-...` / `f116765b-...`; fig `91d10f8e-...` / `9f30caab-...`.
- NOTE: OT ships no r2 for in-sample SuSiE, so Filter 1 LD proxy = 95%-CS co-membership bounded by physical distance.

---

## 6. REGULATORY ENRICHMENT (this session — the pre-modeling diagnostic)

**Question:** before committing to element-class modeling heads, how is each regulatory class
represented among fine-mapped GWAS variants, and which carry causal signal?

**Method:** GWAS layer only (2,191,559 rows). Sparse sets = 51,658 sets with >=1 member PIP>=0.5.
Foreground = 27,040 distinct PIP>=0.5 variants; background = 182,739 PIP<0.5 members of SAME sets
(controls locus/LD). Annotation = ENCODE SCREEN cCRE v3 / ENCODE4 (hg38, 2,345,453 elements, pulled
via UCSC REST API track `cCREregistry` — api.genome.ucsc.edu is allowlisted; wenglab/UCSC-download/
ENCODE-portal are NOT). Multi-label interval overlap; Fisher exact + Katz 95% CI.

**Result:** 38.7% of FG vs 26.9% of BG variants overlap a cCRE.

| cCRE class | element head | fold | Fisher p |
|---|---|---|---|
| Promoter | promoter | 3.00× | 1e-108 |
| Proximal enhancer | enhancer/silencer | 1.93× | 5e-129 |
| Distal enhancer | enhancer/silencer | 1.39× | 2e-150 |
| CA-H3K4me3 | ambiguous OC | 1.41× | 3e-6 |
| CA-TF / TF / CA | ambiguous OC | ~1.0 / depleted | ns / 0.01 |
| **CA-CTCF** | loop anchor | **1.00×** | **ns (0.97)** |

**Three findings that shape the architecture:**
1. Promoter + enhancer-like classes hold the causal signal; enhancer/silencer head is the big target
   (24% of FG hit distal enhancers — largest bucket).
2. **CTCF NOT enriched (1.00×, ns)** → data supports handling loop anchors via a **physics module**
   (loop-extrusion, Sabaté params), NOT the expression head.
3. **Enhancer-vs-silencer is NOT annotatable** — cCREs call both "enhancer-like" by chromatin
   signature; the activating-vs-repressive SIGN is exactly the signed edge the model must learn.
   Enrichment sizes the combined bucket; silencer fraction is non-partitionable upstream.

**Artifacts:** report `989a98eb-fdc1-4caf-97a1-603415cd2546` / `56033326-001d-4279-ad16-28f4a41e1b92`;
table csv `4a2f8ee5-...` / `be75fc5b-...`; labeled variants parquet `34f18bc0-...` / `27844991-...`;
figure `103c8758-...` / `5360cf57-...` + legend `b7ab213a-...` / `f5fa0a12-...`; script
`64033fc0-...` / `fb191f8d-...`; cCRE registry checkpoint `7b827956-3fd9-4aab-9d9d-d14fe064c180` /
`911569fc-8a9c-471b-b1c5-ffc71954ef58`.

---

## 7. THE NEXT STEP — MODEL ARCHITECTURE (agreed design, not yet built)

**User decision:** building own lean **sequence-based signed E2G model** (NOT Borzoi/Borzoi-Prime),
using an **S2F (sequence-to-feature) model** to embed a **1kb element window + 1kb gene-TSS window**.
Goal = learn ONLY the E2G relationship + consensus direction, not an expression track.

**Agreed architecture (from design discussion, in v2g_model_design_doc.md = `0d03c3f6-...` / `1edcac24-...`):**
- **Asymmetric two-tower encoder + interaction head** (EPI lineage: TargetFinder/SPEID/EPIANN).
  Target = EDGE classifier `(element_1kb, tss_1kb, signed_distance, strand, tissue) → {P(link), sign|link}`,
  NOT seq-to-seq (candidate generation is external/distance-based, so no large receptive field needed).
- Two towers (promoter vs enhancer grammar differ), interaction features `[e_E, e_G, e_E−e_G, e_E⊙e_G]`
  (difference term antisymmetric = directional signal). Multi-task heads: link (binary, pos+hard-neg),
  **sign (binary, masked to positives)**, optional magnitude (regress precision-weighted signed-t).
- **KEY REFRAME:** using S2F for features collapses "sequence-based vs gradient-boost" into ONE
  experiment — same embedding matrix, MLP head vs GBDT head, decided by held-out repressive metric.
  Real fork is freeze-vs-train backbone. **FREEZE FIRST** (no GPU on VM, weak sign signal).
- **Element routing (5 consequences, each gets a probability):** splice→frozen SpliceAI; CTCF→physics
  module (not learned — confirmed by enrichment §6); promoter/enhancer/silencer→the E2G expression
  head (the active ideation target).
- **Encoder ranked:** (1) frozen DNA-LM embeddings (Nucleotide Transformer / Caduceus / HyenaDNA) —
  START HERE, CPU-feasible; (2) frozen Enformer/Borzoi trunk embedding; (3) from-scratch dilated CNN
  (needs GPU, do last).

**THREE reviewer failure modes to bake in from the start:**
1. **Distance shortcut** — report sequence's MARGINAL AUPRC over distance (partial-corr logic). Novelty claim lives/dies here.
2. **Sign class imbalance** — 85% activating; eval sign with AUPRC on the REPRESSIVE class, balanced sampling/focal loss. Repression arm is underpowered (concordance couldn't nail it, DC-TAP didn't rescue).
3. **Cross-cell-type / cross-chromosome leakage** — split by CHROMOSOME + hold out whole tissues.

**Recommended first move:** extract frozen DNA-LM embeddings for ~1.3M unique element+TSS windows
(one-time cache), train interaction head TWICE (MLP + GBDT) on identical features with distance
ablation + chromosome split. Answers 3 questions at once: is signed E2G learnable from sequence, does
neural beat GBDT, does sequence beat distance. GPU (Modal or GPU instance) only if that's encouraging.

**Two things to pin before building:** (a) which DNA-LM to embed with (Nucleotide Transformer = safe
default); (b) where GPU comes from (VM has none).

### DBNascent training corpus (the labels)
6.7M FDR<0.01 signed pairs (4.85M FDR<0.001). Base unit = 6,700,460 pairs; also per-tissue obs, negatives,
3D-stringent views. doi 10.1186/s12864-025-11568-z, Zenodo 10.5281/zenodo.14519113. Local tarball at
/Users/jerome/hackathon/pub_resources/DBNascent_data.tar.gz. Trainable subset: 16,677,960 candidate
within-window pairs, 9,977,500 implied negatives. Negatives have 3.42% false-negative contamination
(flagged not dropped) — 1.43× signal-to-contrast ceiling on the label.

### V2G benchmark (validation truth, separate from training)
2,164 positives / 13,613 negatives from 8 sources (kanai_multiome 852, cortex_mpra 494, schnitzler_cad
307, tardaguila_haemvar 208, ghatan_panten 180, fair_hdels 55, sting_seq 46, autoimmune_mpra 22).
positives `964f4b90-...`/`9e13bf36-1fc1-42b3-9432-9aa8d618d2ac`; negatives `e3136dae-...`/`6eceacdb-...`.
NOTE benchmark's kanai_multiome ≠ the dropped Kanai fine-mapping — different source.

---

## 8. INFRASTRUCTURE

**genomics-vm (REBUILT this session, currently TERMINATED):** e2-standard-16 (16 vCPU/64GB/300GB),
project flash-hour-452305-m7, zone us-central1-a, python /opt/micromamba/envs/py313. Agent verified
working (rc=0). Data stack (polars 1.42/pyarrow 25/scipy 1.18/numpy/pandas) installed post-provision —
**will need reinstall if VM rebuilt again** (not baked into provision.sh yet — offered to add).
Start with `start_vm()`, stop with `stop_vm()`, submit via `submit_job()`/`wait_job()`/`run_job()`
(fire-and-forget submit_job then wait_job is most robust; run_job's polling can exceed the 600s
python-tool foreground cap — use background or explicit wait_for_notification for long jobs).

**Bucket:** `gs://claude-hackathon` (HYPHEN — the old underscore `claude_hackathon` was deleted;
underscore is invalid DNS host so virtual-hosted S3 addressing fails). All outputs staged to
`gs://claude-hackathon/multibiobank/20260711/`. FUSE-mounted at `/mnt/claude-hackathon` (now
best-effort/bounded-retry so it can never block the agent — that was the rebuild's root-cause fix).
**Credential metadata still lists buckets=['claude_hackathon']** — user may want to update in
Customize→Credentials.

**Local sandbox:** 48GB RAM handles the full 10M-row combined table via polars lazy scan — the
enrichment ran entirely local. NO HTML→PDF tooling (weasyprint/wkhtmltopdf/pandoc/chrome absent);
use reportlab for PDFs. Blocked domains: oauth2.googleapis.com (can't mint SA token locally),
downloads.wenglab.org, hgdownload.soe.ucsc.edu, encodeproject.org. Allowlisted: api.genome.ucsc.edu,
compute.googleapis.com, standard science APIs + package managers.

**GitHub:** jerome-f/Claude_Life_Sciences_Hackathon_Project (branch master). Latest commits: daab254
(inventory), 2a7992e (enrichment), 7c858d5 (enrichment report fix). git identity jerome-f /
jerome-f@users.noreply.github.com, GIT_CONFIG_GLOBAL=/dev/null, auth x-access-token:${GITHUB_TOKEN}.
Enrichment lives in docs/enrichment/; inventory in docs/inventory/.

**Skills (personal, in github.com/jerome-f/my_claude_skills_repo):** gcp-genomics-vm (updated this
session — bucket rename + non-blocking mount, published skill_01MCAyVaYDoXAFhLS8B1pF8T),
inclusive-plot-style (separate-legend rule tightened), python-data-style (polars/duckdb mandate).

---

## 9. OPEN THREADS / DECISIONS PENDING

1. **Build the signed E2G model** — design agreed (§7), not started. Next concrete action = frozen
   DNA-LM embedding learnability probe. Pin DNA-LM choice + GPU source first.
2. **Bake data-stack install into provision.sh** — offered, user hasn't said yes. Would make VM
   rebuilds fully self-serve.
3. **Update GCP credential bucket metadata** claude_hackathon→claude-hackathon (user-side).
4. **Standing reframe recommendation (carried, never explicitly accepted):** reframe the deliverable
   to (a) concordance-test result + (b) signed-direction audit of SOTA (rE2G/AlphaGenome/L2G) rather
   than building another E2G model. User has instead chosen to build their own model — treat this
   reframe as superseded unless user revives it.
5. **Repression arm underpowered** — needs a repression-enriched perturbation/silencer/CTCF screen if
   the sign head's repressive class is to be validated well.
6. **V2E link is unmeasured** — the variant→element-activity sign (MPRA allelic-skew, base-editing) is
   a SEPARATE factor DBNascent doesn't provide; V2G direction = product V2E×E2G, each with own noise.
   Skipping V2E produces undiagnosable sign inconsistencies (flagged in memory).

---

## 10. KEY REFERENCE ARTIFACTS (read these to go deep)
- **v2g_model_design_doc.md** `0d03c3f6-f1ae-4117-be56-6d1108b7608e` / `1edcac24-331b-474b-a27a-71ad9b466df6` — the consolidated model design.
- **DBNascent_Concordance_Test_REPORT.md** `d00fddd8-...` / `0f5c34aa-777f-45a7-af91-d637e6c20027` — the gate result.
- **Regulatory_Enrichment_REPORT.md** `989a98eb-...` / `56033326-...` — this session's diagnostic.
- **V2G2T_Data_Inventory.pdf** `6de4154a-50b9-47ca-b13e-b91a227b9a0c` / `501a13e1-84c9-4bbf-91dd-23aa12366541` — the full dedup inventory.
- **DATA_DICTIONARY.md** `1d41bf62-...` / `3d6e5374-...` — schema definitions.
- Papers on disk: AoU MultiSuSiE s41588-025-02450-5.pdf (`5bc19828-...`), DBNascent s12864-025-11568-z.pdf (`62ae987c-...`).

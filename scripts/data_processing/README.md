# Data processing — DBNascent E2G & Open Targets molecular QTLs

Processing scripts for the two training/validation data sources of the
E2G→V2G→V2G2T project. All heavy processing runs on the GCP genomics-vm
(e2-standard-16, 64 GB) with duckdb; outputs are staged to
`gs://claude_hackathon/` (dated prefixes).

## DBNascent (signed E2G training label)

| Script | Purpose |
|---|---|
| `dbnascent_parse_qc.py` | Parse the published `dbnascent_pairs.txt.gz` (Zenodo 14519113), QC the 12.7M signed per-tissue observations (sign balance, FDR tiers, 3D-support flags, distance-to-TSS, nulls). |
| `dbnascent_trainable_subset.py` | Build the trainable subset: 6.7M unique pairs with a **cross-tissue consistency statistic** (sign_concordance, sign_class, per-tissue pcc list), a `stringent_3d` flag, and distance-matched hard negatives. |

Key finding: the published file is pre-filtered to FDR<0.01 (no tested-but-uncorrelated
pairs), and 10.6% of pairs are sign-flip pairs (both signs across tissues) — hence
the consistency statistic is preserved rather than mean-collapsed.

## Open Targets (causal variant universe)

| Script | Purpose |
|---|---|
| `ot_extract_qtl.py` | Extract high-confidence molecular-QTL credible sets from OT 26.06 `credible_set` — studyType ∈ {eqtl,pqtl,sqtl,tuqtl,sceqtl} AND confidence = in-sample SuSiE. 2.01M credible sets. |
| `ot_hunt_metabolite.py` | Hunt metabolite-level GWAS (trait-keyword matched among studyType=gwas) as a metabolite proxy — clearly NOT molecular QTL (no metabolite-QTL type exists in this release). |

Note: the OT bucket is requester-pays; reads are billed to the user project via
`gsutil -u <project>` / gcsfs `requester_pays=True`.

## Outputs (staged to GCS, not committed here)

- `gs://claude_hackathon/dbnascent/<date>/processed/` — positives, per-tissue signed, negatives (parquet)
- `gs://claude_hackathon/opentargets/<date>/processed/` — molecular QTLs, metabolite GWAS (parquet)

QC reports & figures: `docs/qc/`.

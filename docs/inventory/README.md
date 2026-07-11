# Data Inventory

Deduplicated catalogue of all datasets acquired for the V2G2T (variant→gene→trait) program.

## Contents
- **V2G2T_Data_Inventory.pdf** — full report: 17 canonical datasets across 6 tiers, deduplication ledger, ancestry composition, provenance & exclusions.
- **data_inventory_dedup.csv** — the canonical dataset table (tier, dataset, layer, ancestry, build, N, key detail, role).
- **dedup_ledger.csv** — every redundancy resolved (combined-vs-per-source, byte-identical copies, dual-role CRISPRi, DBNascent views, etc.) so headline counts are not double-counted.
- **noneuropean_omics_resource_inventory.csv**, **validated_triplet_resource_inventory.csv** — supporting sub-inventories.
- **figures/** — dataset-scale-by-tier (with separate legend), ancestry representation.

## Headline (deduplicated)
- Credible-set substrate: **10,072,543 rows** (5,080,965 distinct variants; 48,411 at PIP>0.9), all hg38, counted once via the combined table.
- Layers: 7,837,680 molecular-QTL · 2,191,559 GWAS · 43,304 metabolite-QTL.
- Single source of truth: **gs://claude-hackathon/multibiobank/20260711/**.

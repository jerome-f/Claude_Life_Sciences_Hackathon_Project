#!/usr/bin/env python3
"""Batch-download + parse FinnGen R13 SuSiE credible sets (finemap/summary/) and LoF burden.
Bucket: finngen-public-data-r13 (public, anonymous listing). hg38.
snp.filter.tsv -> per-variant PIP (cs_specific_prob); cred.summary.tsv -> r2/purity (cs_avg_r2).
Download step uses GNU parallel over the ~1,488 endpoint files (see shell wrapper).
"""
import polars as pl, glob

# ---- credible sets ----
snp = pl.concat([pl.read_csv(f, separator="\t", infer_schema_length=5000,
                 schema_overrides={"chromosome":pl.Utf8,"cs":pl.Utf8}) for f in glob.glob("snp/*.tsv")],
                how="vertical_relaxed")
cred = pl.concat([pl.read_csv(f, separator="\t", infer_schema_length=5000,
                  schema_overrides={"cs":pl.Utf8}) for f in glob.glob("cred/*.tsv")], how="vertical_relaxed")
cred_key = cred.select([pl.col("trait"),pl.col("region"),pl.col("cs").cast(pl.Utf8),
    pl.col("cs_avg_r2").cast(pl.Float64),pl.col("cs_min_r2").cast(pl.Float64),
    pl.col("low_purity").cast(pl.Utf8),pl.col("cs_size").cast(pl.Int64),
    pl.col("good_cs").cast(pl.Utf8),pl.col("cs_log10bf").cast(pl.Float64)]).unique(subset=["trait","region","cs"])
u = snp.join(cred_key, on=["trait","region","cs"], how="left").with_columns([
    pl.col("v").str.replace_all(":","_").alias("variant_hg38"),
    pl.col("chromosome").cast(pl.Utf8).str.replace(r"^chr","").alias("chromosome"),
    pl.col("position").cast(pl.Int64).alias("position"),
    pl.col("allele1").alias("ref"), pl.col("allele2").alias("alt"),
    pl.lit("FinnGen_R13").alias("study"), pl.col("cs_specific_prob").cast(pl.Float64).alias("PIP"),
    (pl.col("trait")+"|"+pl.col("region")+"|"+pl.col("cs")).alias("cs_id"),
    (pl.col("cs").cast(pl.Int64,strict=False)>0).alias("cs_membership"),
    pl.col("cs_avg_r2").alias("r2_to_lead"), pl.lit("FIN").alias("ancestry"),
    pl.lit("imputation").alias("assay_basis"), pl.lit("hg38").alias("build"), pl.lit("FinnGen").alias("source"),
    pl.col("beta").cast(pl.Float64).alias("beta"), pl.col("gene_most_severe").alias("gene"),
])
SCHEMA=["variant_hg38","chromosome","position","ref","alt","trait","study","PIP","cs_id","cs_membership","r2_to_lead","ancestry","assay_basis","build","source","beta","gene"]
u.select(SCHEMA).write_parquet("finngen_credible_sets.parquet")

# ---- LoF burden (G2T resource, significant hits) ----
d = pl.read_csv("finngen_R13_lof_sig.txt", separator="\t", infer_schema_length=20000, null_values=["NA"])
d.select([pl.col("GENE").alias("gene"), pl.col("PHENO").alias("trait"),
    pl.col("CHROM").cast(pl.Utf8).alias("chromosome"),
    pl.col("GENE_BETA").cast(pl.Float64).alias("gene_beta"),
    pl.col("GENE_MLOGP").cast(pl.Float64).alias("gene_mlogp"),
    pl.col("TOP_VAR_MLOGP").cast(pl.Float64,strict=False).alias("top_var_mlogp"),
    pl.lit("FinnGen_R13_LoF").alias("source"), pl.lit("hg38").alias("build"),
    pl.lit("burden_LoF").alias("evidence_type")]).write_parquet("finngen_lof_burden_sig.parquet")

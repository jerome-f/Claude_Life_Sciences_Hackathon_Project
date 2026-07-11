#!/usr/bin/env python3
"""Download + parse AoU MultiSuSiE and UKB-WGS credible sets to unified hg38 schema.
AoU: Zenodo 14458399 pips.tsv (WGS, multi-ancestry). UKB-WGS: Nature s41586-025-09720-6 MOESM4.gz (WGS+imputation panels).
Unified schema: variant_hg38, chromosome, position, ref, alt, trait, study, PIP, cs_id,
cs_membership, r2_to_lead, ancestry, assay_basis, build, source, beta, gene.
"""
import polars as pl

def norm_variant(col):
    return pl.col(col).str.replace(r"^chr","").str.replace_all(":","_")

# ---- AoU ----
aou = pl.read_csv("aou_pips.tsv", separator="\t", infer_schema_length=10000)
def anc(c):
    parts=[p.upper() for p in ("afr","lat","eur") if p in c]
    return "+".join(parts) if parts else "NA"
aou = aou.with_columns([
    norm_variant("rsid").alias("variant_hg38"),
    pl.col("chromosome").cast(pl.Utf8).str.replace(r"^chr","").alias("chromosome"),
    pl.col("bp").cast(pl.Int64).alias("position"),
    pl.lit("AoU_MultiSuSiE").alias("study"),
    pl.col("PIP").cast(pl.Float64).alias("PIP"),
    (pl.col("trait")+"|"+pl.col("cohort")).alias("cs_id"),
    pl.col("cohort").map_elements(anc, return_dtype=pl.Utf8).alias("ancestry"),
    pl.lit("wgs").alias("assay_basis"), pl.lit("hg38").alias("build"), pl.lit("AoU").alias("source"),
    pl.col("beta_0").cast(pl.Float64).alias("beta"),
    pl.lit(None,dtype=pl.Float64).alias("r2_to_lead"), pl.lit(None,dtype=pl.Utf8).alias("gene"),
    pl.lit(True).alias("cs_membership"),
]).with_columns([
    pl.col("variant_hg38").str.split("_").list.get(2).alias("ref"),
    pl.col("variant_hg38").str.split("_").list.get(3).alias("alt"),
])
SCHEMA=["variant_hg38","chromosome","position","ref","alt","trait","study","PIP","cs_id","cs_membership","r2_to_lead","ancestry","assay_basis","build","source","beta","gene"]
aou.select(SCHEMA).write_parquet("aou_credible_sets.parquet")

# ---- UKB-WGS ----
ukb = pl.read_csv("ukbwgs_M4.gz", separator="\t", infer_schema_length=20000, null_values=["NA"]).filter(pl.col("PIP").is_not_null())
ukb = ukb.with_columns([
    norm_variant("CS_SNP").alias("variant_hg38"),
    pl.col("CS_SNP").str.split(":").list.get(0).alias("chromosome"),
    pl.col("CS_SNP").str.split(":").list.get(1).cast(pl.Int64).alias("position"),
    pl.lit("UKB_WGS_h2").alias("study"), pl.col("PIP").cast(pl.Float64).alias("PIP"),
    (pl.col("trait")+"|"+pl.col("loci")+"|"+pl.col("Component")).alias("cs_id"),
    pl.lit("EUR").alias("ancestry"),
    pl.when(pl.col("Set")=="WGS").then(pl.lit("wgs")).otherwise(pl.lit("imputation")).alias("assay_basis"),
    pl.lit("hg38").alias("build"), (pl.lit("UKB_WGS_")+pl.col("Set")).alias("source"),
    pl.lit(None,dtype=pl.Float64).alias("beta"), pl.lit(None,dtype=pl.Float64).alias("r2_to_lead"),
    pl.lit(None,dtype=pl.Utf8).alias("gene"), pl.lit(True).alias("cs_membership"),
]).with_columns([
    pl.col("variant_hg38").str.split("_").list.get(2).alias("ref"),
    pl.col("variant_hg38").str.split("_").list.get(3).alias("alt"),
])
ukb.select(SCHEMA).write_parquet("ukbwgs_credible_sets.parquet")

#!/usr/bin/env python3
"""Extract + parse Japan COVID-19 Task Force (BBJ) blood cis-eQTL/pQTL credible sets.
Source: NBDC hum0343 v3 QTL zip (hum0343.v3.qtl.v1.zip). hg38 native. East-Asian, in-sample LD.
Keeps fine-mapped members (pip_susie > 0.001). "trait" = regulated gene/protein.
"""
import zipfile, polars as pl
z=zipfile.ZipFile("qtl.zip")
for member in [
  "hum0343.v3.qtl.v1/hum0343.v3.eqtl.v1/eqtl_sumstats_wang_qs_et_al_2024_japan_covid19_taskforce.tsv.gz",
  "hum0343.v3.qtl.v1/hum0343.v3.pqtl.v1/pqtl_sumstats_wang_qs_et_al_2024_japan_covid19_taskforce.tsv.gz"]:
    tgt=member.split("/")[-1]
    with z.open(member) as s, open(tgt,"wb") as o:
        while (b:=s.read(1<<20)): o.write(b)

def parse_qtl(path, qtl_type):
    df=pl.read_csv(path, separator="\t", infer_schema_length=10000,
                   schema_overrides={"variant_id_hg38":pl.Utf8,"variant_id_hg19":pl.Utf8})
    df=df.filter(pl.col("pip_susie").cast(pl.Float64)>0.001)
    v=pl.col("variant_id_hg38").str.replace(r"^chr","").str.replace_all(":","_")
    out=df.with_columns([v.alias("variant_hg38"),
        pl.col("variant_id_hg38").str.replace(r"^chr","").str.split(":").list.get(0).alias("chromosome"),
        pl.col("variant_id_hg38").str.replace(r"^chr","").str.split(":").list.get(1).cast(pl.Int64,strict=False).alias("position"),
        pl.lit(qtl_type).alias("study"), pl.col("gene_name").alias("trait"),
        pl.col("pip_susie").cast(pl.Float64).alias("PIP"),
        (pl.lit(qtl_type)+"|"+pl.col("gene_id")).alias("cs_id"),
        (pl.col("pip_susie").cast(pl.Float64)>=0.001).alias("cs_membership"),
        pl.lit(None,dtype=pl.Float64).alias("r2_to_lead"), pl.lit("EAS").alias("ancestry"),
        pl.lit("imputation").alias("assay_basis"), pl.lit("hg38").alias("build"), pl.lit("BBJ_JCTF").alias("source"),
        pl.col("slope").cast(pl.Float64).alias("beta"), pl.col("gene_name").alias("gene")]).with_columns([
        pl.col("variant_hg38").str.split("_").list.get(2).alias("ref"),
        pl.col("variant_hg38").str.split("_").list.get(3).alias("alt")])
    SCHEMA=["variant_hg38","chromosome","position","ref","alt","trait","study","PIP","cs_id","cs_membership","r2_to_lead","ancestry","assay_basis","build","source","beta","gene"]
    return out.select(SCHEMA)

pl.concat([parse_qtl("eqtl_sumstats_wang_qs_et_al_2024_japan_covid19_taskforce.tsv.gz","eqtl"),
           parse_qtl("pqtl_sumstats_wang_qs_et_al_2024_japan_covid19_taskforce.tsv.gz","pqtl")]).write_parquet("bbj_molqtl_credible_sets.parquet")

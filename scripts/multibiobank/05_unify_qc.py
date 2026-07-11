#!/usr/bin/env python3
"""Unify 5 credible-set sources to one hg38 table + cross-source QC.
Normalizes chrom 23->X, drops rows with null position, tags data_layer.
"""
import polars as pl, json
SCHEMA=["variant_hg38","chromosome","position","ref","alt","trait","study","PIP","cs_id","cs_membership","r2_to_lead","ancestry","assay_basis","build","source","beta","gene"]
srcs=["aou_credible_sets.parquet","ukbwgs_credible_sets.parquet","finngen_credible_sets.parquet","bbj_molqtl_credible_sets.parquet","ukb_metabolomics_credible_sets.parquet"]
c=pl.concat([pl.read_parquet(f).select(SCHEMA) for f in srcs], how="vertical_relaxed")
c=c.with_columns(pl.when(pl.col("study").is_in(["eqtl","pqtl"])).then(pl.lit("molQTL"))
    .when(pl.col("source")=="UKB_metabolomics").then(pl.lit("metaboliteQTL")).otherwise(pl.lit("GWAS")).alias("data_layer"))
c=c.with_columns(pl.when(pl.col("chromosome")=="23").then(pl.lit("X")).otherwise(pl.col("chromosome")).alias("chromosome"))
c=c.filter(pl.col("position").is_not_null())
c.write_parquet("multibiobank_credible_sets_combined.parquet")

rows=[]
for src in c["source"].unique().to_list():
    s=c.filter(pl.col("source")==src)
    rows.append({"source":src,"n_rows":s.height,"n_variants":s["variant_hg38"].n_unique(),
      "n_traits":s["trait"].n_unique(),"n_cs":s["cs_id"].n_unique(),
      "pip_gt05":s.filter(pl.col("PIP")>0.5).height,"pip_gt09":s.filter(pl.col("PIP")>0.9).height,
      "pip_median":float(s["PIP"].median() or 0),
      "has_r2":int(s.filter(pl.col("r2_to_lead").is_not_null()).height),
      "has_gene":int(s.filter(pl.col("gene").is_not_null()&(pl.col("gene")!="NA")&(pl.col("gene")!="")).height),
      "ancestry":s["ancestry"].unique().to_list(),"assay":s["assay_basis"].unique().to_list()[0],
      "layer":s["data_layer"].unique().to_list()[0]})
json.dump({"per_source":rows,"total_rows":c.height,"by_layer":c["data_layer"].value_counts().to_dicts(),
    "by_assay":c["assay_basis"].value_counts().to_dicts(),"pip_gt09_total":c.filter(pl.col("PIP")>0.9).height,
    "distinct_variants":c["variant_hg38"].n_unique()}, open("multibiobank_qc.json","w"), indent=2)

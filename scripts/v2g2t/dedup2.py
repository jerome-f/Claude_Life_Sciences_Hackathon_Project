
import duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/v2g2t")
con=duckdb.connect(); con.execute("PRAGMA threads=16")
con.execute(f"CREATE TABLE f2 AS SELECT * FROM read_parquet('{OUT}/v2g2t_filter2_burden_gwas_coloc.parquet')")
con.execute("""CREATE TABLE f2d AS
SELECT gene_id, gene_symbol, any_value(trait) trait, any_value(burden_disease) burden_disease,
  max(burden_score) burden_score, any_value(burden_method) burden_method,
  gwas_study, gwas_lead, any_value(gene_chrom) gene_chrom,
  max(h4) best_h4, max(clpp) best_clpp,
  string_agg(DISTINCT qtl_type, ',') qtl_types, count(DISTINCT qtl_type) n_qtl_types,
  any_value(betaRatioSignAverage) beta_sign, count(DISTINCT qtl_sl) n_coloc_qtl_cs
FROM f2 GROUP BY gene_id, gene_symbol, gwas_study, gwas_lead""")
n=con.execute("SELECT count(*) FROM f2d").fetchone()[0]
ndg=con.execute("SELECT count(DISTINCT gwas_lead||'|'||gene_id) FROM f2d").fetchone()[0]
con.execute(f"COPY f2d TO '{OUT}/v2g2t_filter2_burden_gwas_coloc.parquet' (FORMAT parquet)")
withdir=con.execute("SELECT count(*) FROM f2d WHERE beta_sign IS NOT NULL").fetchone()[0]
multiqtl=con.execute("SELECT count(*) FROM f2d WHERE n_qtl_types>=2").fetchone()[0]
json.dump(dict(deduped_rows=int(n), distinct_variant_gene=int(ndg), with_direction=int(withdir),
  multi_qtl_type_support=int(multiqtl)), open(OUT/"filter2_dedup.json","w"), indent=2)
print(json.dumps(dict(deduped_rows=int(n), distinct_variant_gene=int(ndg), with_direction=int(withdir), multi_qtl_type_support=int(multiqtl)), indent=2))
print(con.execute("SELECT gene_symbol, trait, gwas_lead, qtl_types, round(best_h4,3) h4, round(beta_sign,2) beta_sign, round(burden_score,2) burden FROM f2d ORDER BY best_h4 DESC, burden_score DESC LIMIT 18").df().to_string())

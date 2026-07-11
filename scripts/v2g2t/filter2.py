
import gcsfs, duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/v2g2t")
fs=gcsfs.GCSFileSystem(project="flash-hour-452305-m7",requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='40GB'"); con.register_filesystem(fs)
ROOT="open-targets-data-releases/26.06/output"
def P(ds): return "['"+"','".join("gcs://"+x for x in fs.ls(f"{ROOT}/{ds}/") if x.endswith(".parquet"))+"']"
INSAMPLE="SuSiE fine-mapped credible set with in-sample LD"

# 1) burden LoF gene-disease (score>=0.6 default OT cutoff already; keep all, tag)
con.execute(f"""CREATE TABLE burden AS
SELECT targetId AS gene_id, diseaseId AS burden_disease, max(score) AS burden_score,
  min(pValueExponent) AS burden_pexp, any_value(statisticalMethod) AS method, any_value(projectId) AS burden_project
FROM read_parquet({P('evidence_gene_burden')}) WHERE directionOnTarget='LoF'
GROUP BY 1,2""")
nb=con.execute("SELECT count(*) FROM burden").fetchone()[0]
log.info(f"burden LoF gene-disease pairs: {nb}")

# 2) gene coords
con.execute(f"""CREATE TABLE genes AS
SELECT id AS gene_id, approvedSymbol AS gene_symbol,
  genomicLocation.chromosome AS chrom, genomicLocation.start AS gstart, genomicLocation.end AS gend
FROM read_parquet({P('target')})""")

# 3) GWAS in-sample credible sets + disease ids (explode diseaseIds)
con.execute(f"""CREATE TABLE gwas AS
SELECT c.studyLocusId, c.studyId, c.chromosome, c.position, c.variantId AS lead, d AS disease_id
FROM read_parquet({P('credible_set')}) c
LEFT JOIN read_parquet({P('study')}) s ON c.studyId=s.studyId, unnest(s.diseaseIds) t(d)
WHERE c.studyType='gwas' AND c.confidence='{INSAMPLE}'""")
ng=con.execute("SELECT count(*) FROM gwas").fetchone()[0]
log.info(f"GWAS in-sample CS-disease rows: {ng}")

# 4) burden gene x matching-disease GWAS CS, proximal (<=500kb of gene body)
con.execute("""CREATE TABLE bg AS
SELECT b.gene_id, b.burden_disease, b.burden_score, b.method AS burden_method,
  g.gene_symbol, g.chrom AS gene_chrom, g.gstart, g.gend,
  w.studyLocusId AS gwas_sl, w.studyId AS gwas_study, w.lead AS gwas_lead, w.position AS gwas_pos, w.disease_id
FROM burden b JOIN genes g ON b.gene_id=g.gene_id
JOIN gwas w ON w.disease_id=b.burden_disease AND w.chromosome=g.chrom
WHERE w.position BETWEEN g.gstart-500000 AND g.gend+500000""")
nbg=con.execute("SELECT count(*) FROM bg").fetchone()[0]
nbg_gene=con.execute("SELECT count(DISTINCT gene_id) FROM bg").fetchone()[0]
log.info(f"burden-gene x disease-matched proximal GWAS CS: {nbg} rows, {nbg_gene} genes")

# 5) colocalisation: GWAS CS (left) <-> eQTL CS (right) for SAME gene, H4>=0.8
# need eQTL CS -> gene mapping via study.geneId
con.execute(f"""CREATE TABLE qtl AS
SELECT c.studyLocusId, s.geneId AS qtl_gene, c.studyType AS qtl_type
FROM read_parquet({P('credible_set')}) c JOIN read_parquet({P('study')}) s ON c.studyId=s.studyId
WHERE c.studyType IN ('eqtl','pqtl','sqtl') AND c.confidence='{INSAMPLE}' AND s.geneId IS NOT NULL""")
con.execute(f"""CREATE TABLE coloc AS
SELECT leftStudyLocusId, rightStudyLocusId, h4, clpp, rightStudyType, betaRatioSignAverage
FROM read_parquet({P('colocalisation')}) WHERE h4>=0.8""")
nc=con.execute("SELECT count(*) FROM coloc").fetchone()[0]
log.info(f"coloc pairs H4>=0.8: {nc}")

# 6) join: GWAS CS colocalizes (either side) with an eQTL CS for the burden gene
con.execute("""CREATE TABLE final2 AS
SELECT DISTINCT bg.*, q.qtl_type, co.h4, co.clpp, co.betaRatioSignAverage, q.studyLocusId AS qtl_sl
FROM bg
JOIN coloc co ON bg.gwas_sl=co.leftStudyLocusId
JOIN qtl q ON co.rightStudyLocusId=q.studyLocusId AND q.qtl_gene=bg.gene_id
UNION
SELECT DISTINCT bg.*, q.qtl_type, co.h4, co.clpp, co.betaRatioSignAverage, q.studyLocusId
FROM bg
JOIN coloc co ON bg.gwas_sl=co.rightStudyLocusId
JOIN qtl q ON co.leftStudyLocusId=q.studyLocusId AND q.qtl_gene=bg.gene_id""")
nf=con.execute("SELECT count(*) FROM final2").fetchone()[0]
ndg=con.execute("SELECT count(DISTINCT gwas_lead||'|'||gene_id) FROM final2").fetchone()[0]
ndgene=con.execute("SELECT count(DISTINCT gene_id) FROM final2").fetchone()[0]
log.info(f"FILTER2 convergent (burden+GWAS+eQTL coloc): {nf} rows, {ndg} distinct (variant,gene), {ndgene} genes")

con.execute(f"""CREATE TABLE final2t AS SELECT f.*, s.traitFromSource trait
  FROM final2 f LEFT JOIN read_parquet({P('study')}) s ON f.gwas_study=s.studyId""")
con.execute(f"""COPY (SELECT gene_id,gene_symbol,burden_disease,burden_score,burden_method,trait,
  gwas_study,gwas_lead,gwas_pos,gene_chrom,gstart,gend,disease_id,qtl_type,qtl_sl,h4,clpp,betaRatioSignAverage
  FROM final2t) TO '{OUT}/v2g2t_filter2_burden_gwas_coloc.parquet' (FORMAT parquet)""")
topg=con.execute("SELECT gene_symbol, count(DISTINCT gwas_lead) n FROM final2t GROUP BY 1 ORDER BY n DESC LIMIT 12").df()
summ=dict(burden_lof_gene_disease=int(nb), gwas_insample_cs_disease_rows=int(ng),
  burden_x_proximal_gwas=int(nbg), coloc_h4_ge08=int(nc),
  filter2_convergent_rows=int(nf), distinct_variant_gene=int(ndg), distinct_genes=int(ndgene),
  top_genes={r['gene_symbol']:int(r['n']) for _,r in topg.iterrows()},
  qtl_type_breakdown={r['qtl_type']:int(r['n']) for _,r in con.execute("SELECT qtl_type, count(*) n FROM final2 GROUP BY 1").df().iterrows()},
  definition="Convergent V2G2T: (1) burden test implicates gene G via LoF for disease D (evidence_gene_burden), (2) an in-sample SuSiE GWAS credible set for D lies within 500kb of G, (3) that GWAS CS colocalises (H4>=0.8) with an in-sample eQTL/pQTL/sQTL credible set for the SAME gene G. Three orthogonal lines (coding burden + GWAS fine-map + molecular coloc) all point to G -> strong causal V2G2T.")
json.dump(summ, open(OUT/"filter2_summary.json","w"), indent=2)
print(json.dumps(summ, indent=2))
print("\nsample:")
print(con.execute("SELECT gene_symbol, trait, gwas_lead, qtl_type, round(h4,3) h4, round(burden_score,2) burden, round(betaRatioSignAverage,2) beta_sign FROM final2t ORDER BY h4 DESC LIMIT 15").df().to_string())

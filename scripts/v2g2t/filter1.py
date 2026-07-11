
import gcsfs, duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/v2g2t")
fs=gcsfs.GCSFileSystem(project="flash-hour-452305-m7",requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='40GB'"); con.register_filesystem(fs)
ROOT="open-targets-data-releases/26.06/output"
def P(ds): return "['"+"','".join("gcs://"+x for x in fs.ls(f"{ROOT}/{ds}/") if x.endswith(".parquet"))+"']"
INSAMPLE="SuSiE fine-mapped credible set with in-sample LD"
PLOF_SO=('SO_0001587','SO_0001589','SO_0001574','SO_0001575','SO_0002012','SO_0001578')
CODING_SO=PLOF_SO+('SO_0001583','SO_0001819','SO_0001821','SO_0001822','SO_0001567','SO_0001818','SO_0001580','SO_0001626','SO_0001650')
con.execute(f"""CREATE TABLE gcs AS
SELECT studyLocusId, studyId, chromosome, position AS lead_pos, region, variantId AS lead, locus
FROM read_parquet({P('credible_set')}) WHERE studyType='gwas' AND confidence='{INSAMPLE}'""")
# members with parsed position from variantId (chr_pos_ref_alt)
con.execute("""CREATE TABLE mem AS
SELECT studyLocusId, lead, lead_pos, x.variantId AS member, x.posteriorProbability AS pip,
  CAST(split_part(x.variantId,'_',2) AS BIGINT) AS member_pos
FROM gcs, unnest(locus) t(x) WHERE x.is95CredibleSet""")
con.execute("CREATE TABLE allv AS SELECT DISTINCT member AS v FROM mem UNION SELECT DISTINCT lead FROM gcs")
con.execute(f"""CREATE TABLE vann AS
SELECT v.variantId,
  (SELECT tc.approvedSymbol FROM unnest(v.transcriptConsequences) t(tc) WHERE tc.lofteePrediction='HC' LIMIT 1) AS lof_symbol,
  (SELECT tc.targetId FROM unnest(v.transcriptConsequences) t(tc) WHERE tc.lofteePrediction='HC' LIMIT 1) AS lof_gene,
  (EXISTS(SELECT 1 FROM unnest(v.transcriptConsequences) t(tc) WHERE tc.lofteePrediction='HC') OR v.mostSevereConsequenceId IN {PLOF_SO}) AS is_plof,
  v.mostSevereConsequenceId NOT IN {CODING_SO} AS is_noncoding, v.mostSevereConsequenceId AS conseq
FROM read_parquet({P('variant')}) v WHERE v.variantId IN (SELECT v FROM allv)""")
# lead consequence
con.execute("""CREATE TABLE leadflag AS SELECT g.studyLocusId, g.studyId, g.chromosome, g.lead, g.lead_pos,
  la.is_noncoding AS lead_noncoding, la.conseq AS lead_conseq FROM gcs g LEFT JOIN vann la ON g.lead=la.variantId""")
# pLoF members with distance + pip
con.execute("""CREATE TABLE plof_mem AS
SELECT m.studyLocusId, m.member AS plof_variant, m.pip AS plof_pip, va.lof_symbol AS plof_gene_symbol, va.lof_gene AS plof_gene_id,
  va.conseq AS plof_conseq, abs(m.member_pos - m.lead_pos) AS dist_to_lead
FROM mem m JOIN vann va ON m.member=va.variantId WHERE va.is_plof AND va.lof_gene IS NOT NULL""")
# join: non-coding lead CS that contains a pLoF member
con.execute("""CREATE TABLE joined AS
SELECT lf.*, pm.plof_variant, pm.plof_pip, pm.plof_gene_symbol, pm.plof_gene_id, pm.plof_conseq, pm.dist_to_lead
FROM leadflag lf JOIN plof_mem pm ON lf.studyLocusId=pm.studyLocusId
WHERE lf.lead_noncoding AND lf.lead<>pm.plof_variant""")
tot=con.execute("SELECT count(*) FROM joined").fetchone()[0]
distq=con.execute("SELECT quantile_cont(dist_to_lead,[0.25,0.5,0.75,0.9,1.0]) FROM joined").fetchone()[0]
log.info(f"raw co-member pLoF joins: {tot}; dist quantiles {distq}")
# tiers by physical distance (LD proxy bound)
for lab,cap in [("<=250kb",250000),("<=500kb",500000),("<=1Mb",1000000)]:
    c=con.execute(f"SELECT count(*), count(DISTINCT lead||'|'||plof_gene_id) FROM joined WHERE dist_to_lead<={cap}").fetchone()
    log.info(f"  dist {lab}: {c[0]} rows, {c[1]} distinct (nc,gene)")
# FINAL: keep <=500kb (defensible LD range) as primary, tag tier
con.execute("""CREATE TABLE final1 AS SELECT *,
  CASE WHEN dist_to_lead<=250000 THEN 'high_conf_<=250kb' WHEN dist_to_lead<=500000 THEN 'medium_<=500kb' ELSE 'wide_>500kb' END AS ld_proxy_tier
  FROM joined""")
con.execute(f"""CREATE TABLE final1t AS SELECT f.*, s.traitFromSource trait, s.traitFromSourceMappedIds, s.diseaseIds, s.pubmedId
  FROM final1 f LEFT JOIN read_parquet({P('study')}) s ON f.studyId=s.studyId""")
con.execute(f"""COPY (SELECT studyLocusId,studyId,trait,traitFromSourceMappedIds,diseaseIds,chromosome,
  lead AS nc_variant,lead_conseq,plof_variant,plof_pip,plof_gene_symbol,plof_gene_id,plof_conseq,dist_to_lead,ld_proxy_tier,pubmedId
  FROM final1t WHERE dist_to_lead<=1000000) TO '{OUT}/v2g2t_filter1_ld_noncoding_plof.parquet' (FORMAT parquet)""")
n250=con.execute("SELECT count(*), count(DISTINCT lead||'|'||plof_gene_id) FROM final1t WHERE dist_to_lead<=250000").fetchone()
n500=con.execute("SELECT count(*), count(DISTINCT lead||'|'||plof_gene_id) FROM final1t WHERE dist_to_lead<=500000").fetchone()
summ=dict(gwas_insample_cs=int(con.execute('SELECT count(*) FROM gcs').fetchone()[0]),
  raw_co_member_joins=int(tot), dist_quantiles=[int(x) for x in distq],
  high_conf_le250kb=dict(rows=int(n250[0]), distinct_variant_gene=int(n250[1])),
  medium_le500kb=dict(rows=int(n500[0]), distinct_variant_gene=int(n500[1])),
  definition="Single in-sample SuSiE GWAS credible set: non-coding LEAD + a LOFTEE-HC pLoF variant as a 95% co-member. LD proxy = 95%-CS co-membership BOUNDED by physical distance (OT ships no r2 for in-sample SuSiE). Tier by distance: <=250kb high-conf, <=500kb medium. Distances >500kb flagged wide (loose fine-mapping, weak LD inference).")
json.dump(summ, open(OUT/"filter1_summary.json","w"), indent=2)
print(json.dumps(summ, indent=2))
print("\nHIGH-CONF sample (<=250kb):")
print(con.execute("SELECT trait, lead AS nc_variant, lead_conseq, plof_variant, round(plof_pip,3) plof_pip, plof_gene_symbol, dist_to_lead FROM final1t WHERE dist_to_lead<=250000 ORDER BY dist_to_lead LIMIT 15").df().to_string())

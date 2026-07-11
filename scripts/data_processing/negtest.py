
import gcsfs, duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log=logging.getLogger()
PROJ="flash-hour-452305-m7"
OUT=Path("/data/negtest"); OUT.mkdir(exist_ok=True,parents=True)
fs=gcsfs.GCSFileSystem(project=PROJ, requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'")
con.register_filesystem(fs)
DBN="gs://claude_hackathon/dbnascent/20260711/processed"

# reuse persisted eqtl_g if present, else build
import os
if os.path.exists("/data/negtest/eqtl_g.parquet"):
    log.info("loading cached eqtl_g")
    con.execute("CREATE TABLE eqtl_g AS SELECT * FROM read_parquet('/data/negtest/eqtl_g.parquet')")
else:
    tbase="open-targets-data-releases/26.06/output/target/"
    tp="['"+"','".join("gcs://"+p for p in fs.ls(tbase) if p.endswith('.parquet'))+"']"
    con.execute(f"CREATE TABLE sym2ens AS SELECT approvedSymbol AS symbol, id AS ensg FROM read_parquet({tp}) WHERE approvedSymbol IS NOT NULL")
    sbase="open-targets-data-releases/26.06/output/study/"
    sp="['"+"','".join("gcs://"+p for p in fs.ls(sbase) if p.endswith('.parquet'))+"']"
    con.execute(f"CREATE TABLE study AS SELECT studyId, geneId FROM read_parquet({sp}) WHERE geneId IS NOT NULL")
    cbase="open-targets-data-releases/26.06/output/credible_set/"
    cp="['"+"','".join("gcs://"+p for p in fs.ls(cbase) if p.endswith('.parquet'))+"']"
    INSAMPLE="SuSiE fine-mapped credible set with in-sample LD"
    log.info("expanding eQTL 95%-CS variants (chrom normalized, no chr prefix)")
    con.execute(f"""
    CREATE TABLE eqtl_g AS
    SELECT m.symbol, c.chromosome AS chrom,
           CAST(regexp_extract(u.variantId,'_([0-9]+)_',1) AS BIGINT) AS pos,
           u.posteriorProbability AS pip
    FROM read_parquet({cp}) c
    JOIN study s ON c.studyId=s.studyId
    JOIN sym2ens m ON s.geneId=m.ensg
    CROSS JOIN UNNEST(c.locus) AS t(u)
    WHERE c.studyType='eqtl' AND c.confidence='{INSAMPLE}' AND u.is95CredibleSet=true
      AND m.symbol IS NOT NULL
    """)
    con.execute("COPY eqtl_g TO '/data/negtest/eqtl_g.parquet' (FORMAT parquet)")
nv=con.execute("SELECT count(*) FROM eqtl_g").fetchone()[0]
log.info(f"eqtl_g rows (symbol-mapped)={nv}")
con.execute("CREATE TABLE eqtl_genes AS SELECT DISTINCT symbol FROM eqtl_g")

# pos/neg enhancer tables, chrom normalized to NO prefix to match OT
con.execute(f"""CREATE TABLE pos AS
  SELECT gene_id AS symbol, replace(chrom,'chr','') AS chrom,
    CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) AS enh_start,
    CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) AS enh_stop,
    abs(CAST(distance_tss AS BIGINT)) AS abs_dist
  FROM read_parquet('{DBN}/dbnascent_e2g_positives.parquet')""")
con.execute(f"""CREATE TABLE neg AS
  SELECT gene_id AS symbol, replace(chrom,'chr','') AS chrom,
    CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) AS enh_start,
    CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) AS enh_stop,
    CAST(dist AS BIGINT) AS abs_dist
  FROM read_parquet('{DBN}/dbnascent_e2g_negatives.parquet')""")

def test1(tbl):
    con.execute(f"CREATE OR REPLACE TABLE {tbl}_u AS SELECT p.* FROM {tbl} p JOIN eqtl_genes g ON p.symbol=g.symbol")
    n_univ=con.execute(f"SELECT count(*) FROM {tbl}_u").fetchone()[0]
    hit=con.execute(f"""SELECT count(*) FROM {tbl}_u p WHERE EXISTS (
      SELECT 1 FROM eqtl_g e WHERE e.symbol=p.symbol AND e.chrom=p.chrom
        AND e.pos BETWEEN p.enh_start AND p.enh_stop)""").fetchone()[0]
    return n_univ, hit
log.info("TEST1 pos"); pu,ph=test1("pos")
log.info("TEST1 neg"); nu,nh=test1("neg")
t1=dict(pos_universe=pu,pos_hits=ph,pos_rate=ph/pu,neg_universe=nu,neg_hits=nh,neg_rate=nh/nu,
        enrichment=(ph/pu)/(nh/nu) if nh else None,
        neg_contaminated=nh, neg_contam_pct=100*nh/nu)
log.info(f"TEST1 FIXED: {json.dumps(t1)}")

t4=con.execute("""SELECT 'pos' g, quantile_cont(abs_dist,[0.25,0.5,0.75]) q FROM pos
                  UNION ALL SELECT 'neg', quantile_cont(abs_dist,[0.25,0.5,0.75]) FROM neg""").df()
t4d={r['g']:[float(x) for x in r['q']] for _,r in t4.iterrows()}
json.dump(dict(test1=t1,test4_dist=t4d), open(OUT/"negtest_results.json","w"), indent=2, default=float)
print("RESULTS:",json.dumps(dict(test1=t1,test4_dist=t4d),indent=2,default=float))
log.info("DONE")


import duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log=logging.getLogger()
OUT=Path("/data/negtest")
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'")
DBN="gs://claude_hackathon/dbnascent/20260711/processed"
import gcsfs
fs=gcsfs.GCSFileSystem(project="flash-hour-452305-m7",requester_pays=True); con.register_filesystem(fs)

con.execute("CREATE TABLE eqtl_g AS SELECT * FROM read_parquet('/data/negtest/eqtl_g.parquet')")
con.execute("CREATE TABLE eqtl_genes AS SELECT DISTINCT symbol FROM eqtl_g")

con.execute(f"""CREATE TABLE pos AS
  SELECT gene_id AS symbol, replace(chrom,'chr','') AS chrom,
    CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) es,
    CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) ee,
    abs(CAST(distance_tss AS BIGINT)) d FROM read_parquet('{DBN}/dbnascent_e2g_positives.parquet')""")
con.execute(f"""CREATE TABLE neg AS
  SELECT gene_id AS symbol, replace(chrom,'chr','') AS chrom, bidirectional_id,
    CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) es,
    CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) ee,
    CAST(dist AS BIGINT) d FROM read_parquet('{DBN}/dbnascent_e2g_negatives.parquet')""")

# distance-stratified enrichment: bins of 100kb, restricted to eqtl-gene universe
def strat_rate(tbl):
    con.execute(f"CREATE OR REPLACE TABLE {tbl}_u AS SELECT p.*, CAST(p.d/100000 AS INT) bin FROM {tbl} p JOIN eqtl_genes g ON p.symbol=g.symbol")
    con.execute(f"""CREATE OR REPLACE TABLE {tbl}_hit AS
      SELECT p.*, EXISTS(SELECT 1 FROM eqtl_g e WHERE e.symbol=p.symbol AND e.chrom=p.chrom
        AND e.pos BETWEEN p.es AND p.ee) AS hit FROM {tbl}_u p""")
    return con.execute(f"SELECT bin, count(*) n, sum(hit::int) hits, avg(hit::int) rate FROM {tbl}_hit GROUP BY 1 ORDER BY 1").df()
pr=strat_rate("pos"); nr=strat_rate("neg")
m=pr.merge(nr, on="bin", suffixes=("_pos","_neg"))
m=m[m["bin"]<=10]
m["enrichment"]=m["rate_pos"]/m["rate_neg"]
m["dist_kb"]=m["bin"]*100
print("DISTANCE-STRATIFIED ENRICHMENT (100kb bins):")
print(m[["dist_kb","n_pos","rate_pos","n_neg","rate_neg","enrichment"]].to_string(index=False))
m.to_csv(OUT/"negtest_stratified.csv", index=False)

# emit cleaned negatives with contamination flag (the neg_hit table has hit=true for contaminated)
con.execute(f"""COPY (
  SELECT symbol AS gene_id, 'chr'||chrom AS chrom, bidirectional_id, d AS dist, hit AS eqtl_supported
  FROM neg_hit
) TO '/data/negtest/dbnascent_negatives_flagged.parquet' (FORMAT parquet, COMPRESSION zstd)""")
# note: neg_hit lost bidirectional_id? it has it (from neg). check
nclean=con.execute("SELECT count(*) FROM neg_hit WHERE hit=false").fetchone()[0]
ncontam=con.execute("SELECT count(*) FROM neg_hit WHERE hit=true").fetchone()[0]
log.info(f"clean_negatives={nclean} contaminated_flagged={ncontam}")
json.dump(dict(clean=int(nclean), contaminated=int(ncontam),
               strat=m[["dist_kb","rate_pos","rate_neg","enrichment"]].to_dict(orient="records")),
          open(OUT/"negtest_summary2.json","w"), indent=2, default=float)
print("clean_neg",nclean,"contam",ncontam)

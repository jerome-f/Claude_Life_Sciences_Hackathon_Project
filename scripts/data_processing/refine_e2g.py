import duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/e2g")
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'")
con.execute(f"CREATE TABLE e2g AS SELECT * FROM read_parquet('{OUT}/ot_encode_re2g_highconf.parquet')")
# distal = intergenic + genic, exclude promoter self-links; require real distance
con.execute("""CREATE TABLE distal AS SELECT * FROM e2g
  WHERE intervalType IN ('intergenic','genic') AND abs(distanceToTss)>=1000""")
n=con.execute("SELECT count(*) FROM distal").fetchone()[0]
ng=con.execute("SELECT count(DISTINCT geneId) FROM distal").fetchone()[0]
nb=con.execute("SELECT count(DISTINCT biosampleId) FROM distal").fetchone()[0]
dq=con.execute("SELECT quantile_cont(abs(distanceToTss),[0.25,0.5,0.75,0.95]) FROM distal").fetchone()[0]
con.execute(f"COPY distal TO '{OUT}/ot_encode_re2g_distal.parquet' (FORMAT parquet, COMPRESSION zstd)")
log.info(f"distal E2G links={n} genes={ng} biosamples={nb}")
print("distal_links", n, "genes", ng, "biosamples", nb)
print("distance quantiles [25,50,75,95]:", [int(x) for x in dq])
s=json.load(open(OUT/"ot_encode_re2g_summary.json"))
s["distal_enhancer_links"]=dict(n=int(n), genes=int(ng), biosamples=int(nb),
    distance_quantiles=[int(x) for x in dq],
    note="intervalType in (intergenic,genic) & |distanceToTss|>=1kb; excludes 24.4M promoter self-links. This is the distal-enhancer E2G signal comparable to DBNascent bidirectional->gene pairs.")
json.dump(s, open(OUT/"ot_encode_re2g_summary.json","w"), indent=2)

import gcsfs, duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/e2g"); OUT.mkdir(exist_ok=True,parents=True)
fs=gcsfs.GCSFileSystem(project="flash-hour-452305-m7",requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'"); con.register_filesystem(fs)
base="open-targets-data-releases/26.06/output/enhancer_to_gene/"
cp="['"+"','".join("gcs://"+p for p in fs.ls(base) if p.endswith('.parquet'))+"']"

# high-confidence E2G (score>=0.8), flatten DNase/HiC
log.info("extracting high-confidence ENCODE-rE2G links (score>=0.8)")
con.execute(f"""CREATE TABLE e2g AS
SELECT geneId, chromosome, start, "end", score, distanceToTss,
  biosampleId, biosampleName, datasourceId, intervalType, studyId, intervalId,
  list_extract(list_transform(list_filter(resourceScore, x->x.name='DNase'), x->x.value),1) AS dnase,
  list_extract(list_transform(list_filter(resourceScore, x->x.name='HiC_contacts'), x->x.value),1) AS hic
FROM read_parquet({cp}) WHERE score>=0.8""")
n=con.execute("SELECT count(*) FROM e2g").fetchone()[0]
log.info(f"high-conf E2G links (>=0.8)={n}")
con.execute(f"COPY e2g TO '{OUT}/ot_encode_re2g_highconf.parquet' (FORMAT parquet, COMPRESSION zstd)")

# QC
sc=con.execute(f"SELECT quantile_cont(score,[0.5,0.8,0.9,0.99]) FROM read_parquet({cp})").fetchone()[0]
ng=con.execute("SELECT count(DISTINCT geneId) FROM e2g").fetchone()[0]
nb=con.execute("SELECT count(DISTINCT biosampleId) FROM e2g").fetchone()[0]
bytype=con.execute("SELECT intervalType, count(*) n FROM e2g GROUP BY 1 ORDER BY n DESC").df()
dqc=con.execute("SELECT quantile_cont(abs(distanceToTss),[0.25,0.5,0.75]) FROM e2g").fetchone()[0]
summary=dict(total_links_all=48810390, highconf_links_ge08=int(n), genes=int(ng), biosamples=int(nb),
   score_quantiles_all=[float(x) for x in sc], dist_quantiles_highconf=[float(x) for x in dqc],
   by_intervalType={r['intervalType']:int(r['n']) for _,r in bytype.iterrows()},
   source="Open Targets 26.06 enhancer_to_gene (ENCODE-rE2G, Gschwind et al 2023; CRISPRi-trained; OT pre-filters at 0.6)")
json.dump(summary, open(OUT/"ot_encode_re2g_summary.json","w"), indent=2)
print("SUMMARY:", json.dumps(summary, indent=2))

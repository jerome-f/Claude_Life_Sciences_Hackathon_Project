
import gcsfs, duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log=logging.getLogger()
PROJ="flash-hour-452305-m7"
OUT=Path("/data/ot_out"); OUT.mkdir(exist_ok=True,parents=True)

fs=gcsfs.GCSFileSystem(project=PROJ, requester_pays=True)
base="open-targets-data-releases/26.06/output/credible_set/"
parts=[p for p in fs.ls(base) if p.endswith(".parquet")]
paths="['" + "','".join("gcs://"+p for p in parts) + "']"
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'")
con.register_filesystem(fs)

QTL_TYPES=('eqtl','pqtl','sqtl','tuqtl','sceqtl')
INSAMPLE="SuSiE fine-mapped credible set with in-sample LD"

log.info("extracting high-confidence molecular QTL credible sets")
con.execute(f"""
CREATE TABLE qtl AS
SELECT
  studyLocusId, studyId, variantId, chromosome, "position",
  studyType, confidence, finemappingMethod, credibleSetIndex,
  beta, zScore, pValueMantissa, pValueExponent, standardError,
  effectAlleleFrequencyFromSource, credibleSetlog10BF,
  purityMeanR2, purityMinR2, locusStart, locusEnd, sampleSize,
  isTransQtl, region, subStudyDescription,
  len(locus) AS credible_set_size,
  list_max(list_transform(locus, x -> x.posteriorProbability)) AS max_pip,
  list_sum(list_transform(locus, x -> CASE WHEN x.is95CredibleSet THEN 1 ELSE 0 END)) AS n_in_95,
  list_sum(list_transform(locus, x -> CASE WHEN x.is99CredibleSet THEN 1 ELSE 0 END)) AS n_in_99,
  len(qualityControls) AS n_qc_flags
FROM read_parquet({paths})
WHERE studyType IN {QTL_TYPES} AND confidence = '{INSAMPLE}'
""")
n=con.execute("SELECT count(*) FROM qtl").fetchone()[0]
log.info(f"high-confidence molecular QTL credible sets = {n}")

# QC tables
by_type=con.execute("SELECT studyType, count(*) n, count(DISTINCT studyId) studies, count(DISTINCT variantId) lead_variants FROM qtl GROUP BY 1 ORDER BY n DESC").df()
by_type.to_csv(OUT/"qc_ot_by_studytype.csv", index=False)
print("BY TYPE:\n", by_type.to_string(index=False))

trans=con.execute("SELECT studyType, isTransQtl, count(*) n FROM qtl GROUP BY 1,2 ORDER BY 1,2").df()
trans.to_csv(OUT/"qc_ot_trans.csv", index=False)
print("TRANS:\n", trans.to_string(index=False))

qc_flags=con.execute("SELECT n_qc_flags, count(*) n FROM qtl GROUP BY 1 ORDER BY 1").df()
qc_flags.to_csv(OUT/"qc_ot_qcflags.csv", index=False)
print("QC FLAG COUNTS:\n", qc_flags.to_string(index=False))

qd={}
for col,expr in [("cs_size_q","credible_set_size"),("max_pip_q","max_pip"),
                 ("purity_q","purityMeanR2"),("purity_min_q","purityMinR2")]:
    nn=con.execute(f"SELECT count(*) FROM qtl WHERE {expr} IS NOT NULL").fetchone()[0]
    if nn>0:
        vals=con.execute(f"SELECT quantile_cont({expr},[0.05,0.25,0.5,0.75,0.95]) FROM qtl WHERE {expr} IS NOT NULL").fetchone()[0]
        qd[col]=[float(x) for x in vals]
    else:
        qd[col]=None
    qd[col+"_nonnull"]=int(nn)
json.dump(qd, open(OUT/"qc_ot_quantiles.json","w"), indent=2)
print("QUANTILES:", json.dumps(qd, indent=2))

chrom=con.execute("SELECT chromosome, count(*) n FROM qtl GROUP BY 1 ORDER BY try_cast(chromosome AS INT) NULLS LAST").df()
chrom.to_csv(OUT/"qc_ot_chrom.csv", index=False)

# write parquet
con.execute(f"COPY (SELECT * FROM qtl) TO '{OUT}/ot_molecular_qtl_highconf.parquet' (FORMAT parquet, COMPRESSION zstd)")
n_null_pip=con.execute("SELECT sum(CASE WHEN max_pip IS NULL THEN 1 ELSE 0 END) FROM qtl").fetchone()[0]
summary=dict(n_credible_sets=int(n), n_null_maxpip=int(n_null_pip),
             by_studytype=by_type.set_index('studyType')['n'].to_dict())
json.dump(summary, open(OUT/"ot_qtl_summary.json","w"), indent=2, default=int)
print("SUMMARY:", json.dumps(summary, indent=2, default=int))
log.info("DONE")

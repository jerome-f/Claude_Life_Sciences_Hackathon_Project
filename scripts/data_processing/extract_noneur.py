import gcsfs, duckdb, pandas as pd, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/noneur"); OUT.mkdir(exist_ok=True,parents=True)
PROJ="flash-hour-452305-m7"
fs=gcsfs.GCSFileSystem(project=PROJ,requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'"); con.register_filesystem(fs)

# Literature-curated ancestry annotation for OT QTL projects (OT carries NO ancestry for QTLs).
# Only projects with ESTABLISHED non-European content are flagged non_european=True.
ANNOT = {
 "GEUVADIS":      ("EUR (CEU/FIN/GBR/TSI) + African (YRI)", True,  "African"),
 "Quach_2016":    ("African + European monocytes (Belgium)", True, "African"),
 "Nedelec_2016":  ("African + European macrophages", True,        "African"),
 "Randolph_2021": ("African + European (flu response)", True,     "African"),
 "Perez_2022":    ("East Asian + European (SLE)", True,           "East Asian"),
 "Nathan_2022":   ("Peruvian TB cohort (admixed/Indigenous American)", True, "Admixed American"),
 "iPSCORE":       ("Multi-ethnic Californian iPSC (EUR/EAS/Hispanic/AA)", True, "Multi"),
}
# projects that POOL an unlabelled non-European subset but OT can't resolve it
POOLED = {"GTEx":"~85% EUR + ~12% African American (pooled, unresolved)",
          "BrainSeq":"EUR + African American (pooled, unresolved)"}

sbase="open-targets-data-releases/26.06/output/study/"
sp="['"+"','".join("gcs://"+p for p in fs.ls(sbase) if p.endswith('.parquet'))+"']"
cbase="open-targets-data-releases/26.06/output/credible_set/"
cp="['"+"','".join("gcs://"+p for p in fs.ls(cbase) if p.endswith('.parquet'))+"']"
INSAMPLE="SuSiE fine-mapped credible set with in-sample LD"
con.execute(f"CREATE TABLE st AS SELECT studyId, projectId, studyType, geneId, biosampleFromSourceId FROM read_parquet({sp})")

projs=list(ANNOT.keys())
plist="('"+"','".join(projs)+"')"
log.info(f"extracting high-conf credible sets for {len(projs)} non-European projects")
con.execute(f"""CREATE TABLE ne AS
SELECT c.studyLocusId, c.studyId, s.projectId, c.studyType, s.geneId,
       c.chromosome, c.position AS lead_pos, c.variantId AS lead_variant,
       c.beta, c.credibleSetlog10BF, c.finemappingMethod,
       len(c.locus) AS cs_size,
       (SELECT max(u.posteriorProbability) FROM UNNEST(c.locus) AS t(u)) AS max_pip
FROM read_parquet({cp}) c JOIN st s ON c.studyId=s.studyId
WHERE s.projectId IN {plist} AND c.confidence='{INSAMPLE}'
  AND c.studyType IN ('eqtl','pqtl','sqtl','tuqtl','sceqtl')""")
n=con.execute("SELECT count(*) FROM ne").fetchone()[0]
log.info(f"non-European high-conf credible sets={n}")

# annotate ancestry
adf=pd.DataFrame([(k,v[0],v[2]) for k,v in ANNOT.items()], columns=["projectId","ancestry_desc","non_eur_ancestry"])
con.execute("CREATE TABLE anc AS SELECT * FROM adf")
con.execute("""CREATE TABLE ne_annot AS
  SELECT n.*, a.ancestry_desc, a.non_eur_ancestry FROM ne n LEFT JOIN anc a ON n.projectId=a.projectId""")
con.execute(f"COPY ne_annot TO '{str(OUT)}/ot_noneur_qtl_highconf.parquet' (FORMAT parquet, COMPRESSION zstd)")

# summary
byp=con.execute("SELECT projectId, studyType, count(*) n, count(DISTINCT geneId) genes FROM ne_annot GROUP BY 1,2 ORDER BY n DESC").df()
print("=== non-European OT QTL credible sets ===")
print(byp.to_string(index=False))
print("\nTOTAL:", n)
summary=dict(total_credible_sets=int(n),
   by_project=con.execute("SELECT projectId, count(*) n FROM ne_annot GROUP BY 1 ORDER BY n DESC").df().set_index("projectId")["n"].to_dict(),
   by_type=con.execute("SELECT studyType, count(*) n FROM ne_annot GROUP BY 1 ORDER BY n DESC").df().set_index("studyType")["n"].to_dict(),
   annotation=ANNOT, pooled_unresolved=POOLED,
   note="OT 26.06 carries NO ancestry for molecular QTLs; ancestry is literature-curated at project level. Only projects with established non-European content included. GTEx/BrainSeq pool unlabelled non-European subsets OT cannot resolve.")
json.dump(summary, open(OUT/"ot_noneur_qtl_summary.json","w"), indent=2, default=str)
log.info("DONE")

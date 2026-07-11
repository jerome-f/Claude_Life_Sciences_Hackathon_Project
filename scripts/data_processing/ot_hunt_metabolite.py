
import gcsfs, duckdb, json, logging, re
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log=logging.getLogger()
PROJ="flash-hour-452305-m7"
OUT=Path("/data/ot_out"); OUT.mkdir(exist_ok=True,parents=True)
fs=gcsfs.GCSFileSystem(project=PROJ, requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'")
con.register_filesystem(fs)

sbase="open-targets-data-releases/26.06/output/study/"
sparts=[p for p in fs.ls(sbase) if p.endswith(".parquet")]
spaths="['" + "','".join("gcs://"+p for p in sparts) + "']"

# flatten study with mapped EFO ids
con.execute(f"""
CREATE TABLE study AS
SELECT studyId, studyType, traitFromSource,
  array_to_string(list_transform(traitFromSourceMappedIds, x->x), ',') AS efo_ids,
  publicationJournal, publicationFirstAuthor, nSamples
FROM read_parquet({spaths})
""")
nstudy=con.execute("SELECT count(*) FROM study").fetchone()[0]
log.info(f"total studies={nstudy}")

# metabolite detection: EFO/metabolite concentration terms + trait keyword patterns
# EFO metabolite measurement subtree markers + common metabolomics keywords
kw = ["metabolite","metabolomic","lipid","fatty acid","amino acid","lipoprotein",
      "cholesterol","triglyceride","phospholipid","sphingomyelin","acylcarnitine",
      "creatinine","urate","glucose","cortisol","bilirubin","ketone","citrate",
      "phosphatidyl","steroid","bile acid","carnitine","choline","betaine",
      "glutamine","glycine","alanine","valine","leucine","tyrosine","serum metabolite",
      "plasma metabolite","concentration of"]
kwpat = "|".join(re.escape(k) for k in kw)

con.execute(f"""
CREATE TABLE metab_studies AS
SELECT *,
  CASE
    WHEN regexp_matches(lower(coalesce(traitFromSource,'')), '{kwpat}') THEN 'trait_keyword'
    ELSE 'other' END AS metab_reason
FROM study
WHERE studyType='gwas'
  AND regexp_matches(lower(coalesce(traitFromSource,'')), '{kwpat}')
""")
nmetab=con.execute("SELECT count(*) FROM metab_studies").fetchone()[0]
log.info(f"candidate metabolite GWAS studies={nmetab}")
top=con.execute("SELECT traitFromSource, count(*) n FROM metab_studies GROUP BY 1 ORDER BY n DESC LIMIT 30").df()
print("TOP METAB TRAITS:\n", top.to_string(index=False))

# link to gwas credible sets (in-sample OR out-of-sample SuSiE, not PICS) for these studies
cbase="open-targets-data-releases/26.06/output/credible_set/"
cparts=[p for p in fs.ls(cbase) if p.endswith(".parquet")]
cpaths="['" + "','".join("gcs://"+p for p in cparts) + "']"
con.execute(f"""
CREATE TABLE metab_cs AS
SELECT c.studyLocusId, c.studyId, c.variantId, c.chromosome, c."position",
  c.studyType, c.confidence, c.finemappingMethod, c.beta, c.zScore,
  c.pValueMantissa, c.pValueExponent, c.credibleSetlog10BF,
  c.purityMeanR2, c.isTransQtl,
  len(c.locus) AS credible_set_size,
  list_max(list_transform(c.locus, x -> x.posteriorProbability)) AS max_pip,
  m.traitFromSource, m.efo_ids, m.metab_reason, m.publicationFirstAuthor
FROM read_parquet({cpaths}) c
JOIN metab_studies m ON c.studyId=m.studyId
WHERE c.confidence LIKE 'SuSiE%'
""")
ncs=con.execute("SELECT count(*) FROM metab_cs").fetchone()[0]
log.info(f"metabolite-GWAS SuSiE credible sets={ncs}")
byconf=con.execute("SELECT confidence, count(*) n FROM metab_cs GROUP BY 1 ORDER BY n DESC").df()
print("METAB CS BY CONFIDENCE:\n", byconf.to_string(index=False))

con.execute(f"COPY (SELECT * FROM metab_studies) TO '{OUT}/ot_metabolite_gwas_studies.parquet' (FORMAT parquet, COMPRESSION zstd)")
con.execute(f"COPY (SELECT * FROM metab_cs) TO '{OUT}/ot_metabolite_gwas_crediblesets.parquet' (FORMAT parquet, COMPRESSION zstd)")
top.to_csv(OUT/"qc_metab_top_traits.csv", index=False)
summary=dict(n_metab_studies=int(nmetab), n_metab_crediblesets=int(ncs),
             note="Metabolite GWAS (trait-keyword matched among studyType=gwas). NOT molecular QTLs; supplied as proxy per user request.")
json.dump(summary, open(OUT/"ot_metabolite_summary.json","w"), indent=2, default=int)
print("SUMMARY:", json.dumps(summary, indent=2, default=int))
log.info("DONE")

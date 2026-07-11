import gcsfs, duckdb, pandas as pd, numpy as np, json, logging
from pathlib import Path
from scipy import stats
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/concord"); OUT.mkdir(exist_ok=True,parents=True); rng=np.random.default_rng(11)
fs=gcsfs.GCSFileSystem(project="flash-hour-452305-m7",requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='40GB'"); con.register_filesystem(fs)
DBN="gs://claude_hackathon/dbnascent/20260711/processed"

# CRISPRi
tr=pd.read_csv("/data/crispri/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz", sep="\t")
ho=pd.read_csv("/data/crispri/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz", sep="\t")
cr=pd.concat([tr,ho],ignore_index=True); cr["chrom"]=cr["chrom"].astype("string").astype(str)
cr=cr[cr["ValidConnection"]==True].copy()
con.register("cr", cr)

# DBNascent per-tissue with signed t-stat + nObs (the precision-aware z-score)
con.execute(f"""CREATE TABLE dbn AS
SELECT gene_id, chrom,
  CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) es,
  CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) ee,
  pcc, t AS tstat, nObs, adj_p_BH, tissue
FROM read_parquet('{DBN}/dbnascent_e2g_pertissue_signed.parquet')""")

# cell-type map: which DBNascent tissue matches each CRISPRi cell line
ctmap={"K562":"blood","HCT116":"intestine","Jurkat":"blood","THP1":"blood","GM12878":"blood","NCIH1299":None,"WTC11":None,"MCF7":"breast","A549":None}
def matched(cell, tissue):
    q=f"""SELECT c.EffectSize, c.pValueAdjusted, c.Significant, c.Regulated, c.PowerAtEffectSize25, c.CellType,
      d.pcc, d.tstat, d.nObs
      FROM cr c JOIN dbn d
       ON c.chrom=d.chrom AND c.measuredGeneSymbol=d.gene_id AND d.tissue='{tissue}'
       AND c.chromStart<=d.ee AND c.chromEnd>=d.es
      WHERE c.CellType='{cell}'"""
    return con.execute(q).df()

# primary: blood x K562 (largest matched cell)
m=matched("K562","blood")
log.info(f"blood x K562 matched (per-tissue, with t): {len(m)}")

def analyze(df, label, sig=True):
    d=df.dropna(subset=["pcc","tstat","EffectSize"]).copy(); d=d[d["pcc"]!=0]
    if sig: d=d[d["pValueAdjusted"]<0.05]
    n=len(d)
    if n<10: return dict(label=label,n=n,note="insufficient")
    out=dict(label=label,n=int(n))
    # SIGN+MAGNITUDE combined, three readouts of the DBNascent side
    for nm,col in [("raw_pcc","pcc"),("signed_t","tstat")]:
        rho,p=stats.spearmanr(d[col], d["EffectSize"])
        pear,pp=stats.pearsonr(d[col], d["EffectSize"])
        out[f"spearman_{nm}_vs_effect"]=float(rho); out[f"spearman_{nm}_p"]=float(p)
        out[f"pearson_{nm}_vs_effect"]=float(pear)
    # MAGNITUDE only, among regulated (the test that was 0.07 with |pcc|)
    reg=d[d["Regulated"]==True]
    if len(reg)>=15:
        for nm,col in [("abs_pcc","pcc"),("abs_t","tstat")]:
            rho,p=stats.spearmanr(reg[col].abs(), reg["EffectSize"].abs())
            out[f"mag_{nm}_vs_absEffect_rho"]=float(rho); out[f"mag_{nm}_p"]=float(p)
        out["n_regulated"]=int(len(reg))
    return out, d

res={}
res["blood_K562"],d_bk = analyze(m,"blood x K562 (per-tissue signed-t)")
# precision-weighted: well-measured both sides (nObs>=50 AND well-powered CRISPRi)
mp=m[(m["nObs"]>=50)&(m["PowerAtEffectSize25"]>0.8)]
r2=analyze(mp,"blood x K562, well-measured both sides (nObs>=50 & power>0.8)")
res["blood_K562_precision"]=r2[0] if isinstance(r2,tuple) else r2

# pooled Stouffer combined z across tissues per pair, vs pooled CRISPRi
con.execute("""CREATE TABLE dbn_z AS
SELECT gene_id, chrom, es, ee,
  sum(sign(tstat)*abs(tstat))/sqrt(count(*)) AS stouffer_z,
  avg(pcc) AS mean_pcc, sum(nObs) AS tot_obs, count(*) AS k
FROM dbn GROUP BY 1,2,3,4""")
mm=con.execute("""SELECT c.EffectSize,c.pValueAdjusted,c.Significant,c.Regulated,c.CellType,
   z.stouffer_z, z.mean_pcc, z.tot_obs, z.k
   FROM cr c JOIN dbn_z z ON c.chrom=z.chrom AND c.measuredGeneSymbol=z.gene_id
    AND c.chromStart<=z.ee AND c.chromEnd>=z.es""").df()
mm=mm.dropna(subset=["stouffer_z","EffectSize"]); mm=mm[mm["pValueAdjusted"]<0.05]
rho_sz,p_sz=stats.spearmanr(mm["stouffer_z"],mm["EffectSize"])
rho_mp,p_mp=stats.spearmanr(mm["mean_pcc"],mm["EffectSize"])
reg=mm[mm["Regulated"]==True]
mag_sz=stats.spearmanr(reg["stouffer_z"].abs(),reg["EffectSize"].abs())
mag_mp=stats.spearmanr(reg["mean_pcc"].abs(),reg["EffectSize"].abs())
res["pooled_stouffer"]=dict(n=int(len(mm)), n_regulated=int(len(reg)),
  spearman_stouffer_z_vs_effect=float(rho_sz), p=float(p_sz),
  spearman_mean_pcc_vs_effect=float(rho_mp),
  mag_stouffer_z_vs_absEffect=float(mag_sz.statistic), mag_stouffer_p=float(mag_sz.pvalue),
  mag_mean_pcc_vs_absEffect=float(mag_mp.statistic))

json.dump(res, open(OUT/"concordance_zscore.json","w"), indent=2, default=float)
d_bk.to_parquet(OUT/"zscore_matched_bk.parquet", index=False)
print(json.dumps(res, indent=2, default=float))

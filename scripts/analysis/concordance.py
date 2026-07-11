import gcsfs, duckdb, pandas as pd, numpy as np, json, logging
from pathlib import Path
from scipy import stats
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/concord"); OUT.mkdir(exist_ok=True,parents=True)
rng=np.random.default_rng(42)
PROJ="flash-hour-452305-m7"
fs=gcsfs.GCSFileSystem(project=PROJ,requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='40GB'"); con.register_filesystem(fs)
DBN="gs://claude_hackathon/dbnascent/20260711/processed"

# ---- load CRISPRi (training K562 + heldout), combine ----
tr=pd.read_csv("/data/crispri/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz", sep="\t")
ho=pd.read_csv("/data/crispri/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz", sep="\t")
cr=pd.concat([tr,ho], ignore_index=True)
cr["chrom"]=cr["chrom"].astype("string").astype(str)
cr=cr[cr["ValidConnection"]==True].copy()
log.info(f"CRISPRi valid pairs: {len(cr)} (K562={len(tr)}, heldout={len(ho)})")
con.register("cr", cr)

# ---- DBNascent per-tissue (blood) + pooled ----
# pooled per-pair sign
con.execute(f"""CREATE TABLE dbn_pool AS
SELECT gene_id, chrom,
  CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) es,
  CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) ee,
  mean_pcc, consensus_sign, sign_class, n_tissues, abs(distance_tss) abs_dist
FROM read_parquet('{DBN}/dbnascent_e2g_positives.parquet')""")
# blood per-tissue
con.execute(f"""CREATE TABLE dbn_blood AS
SELECT gene_id, chrom,
  CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) es,
  CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) ee,
  pcc AS pcc_blood, adj_p_BH
FROM read_parquet('{DBN}/dbnascent_e2g_pertissue_signed.parquet') WHERE tissue='blood'""")

# ---- overlap join: enhancer interval overlaps CRISPRi element, same gene symbol, same chrom ----
def match(dbtbl, pcc_col, extra=""):
    q=f"""
    SELECT c.*, d.{pcc_col} AS dbn_pcc {extra}
    FROM cr c JOIN {dbtbl} d
      ON c.chrom=d.chrom AND c.measuredGeneSymbol=d.gene_id
     AND c.chromStart <= d.ee AND c.chromEnd >= d.es
    """
    return con.execute(q).df()

m_pool=match("dbn_pool","mean_pcc", ", d.consensus_sign, d.sign_class, d.abs_dist, d.n_tissues")
m_blood=match("dbn_blood","pcc_blood", ", d.adj_p_BH AS dbn_adjp")
log.info(f"matched pairs: pooled={len(m_pool)}, blood={len(m_blood)}")

def concord_stats(df, pcc_col="dbn_pcc", sig_only=True, label=""):
    d=df.dropna(subset=[pcc_col,"EffectSize"]).copy()
    d=d[d[pcc_col]!=0]
    if sig_only:
        d=d[d["pValueAdjusted"]<0.05]
    n=len(d)
    if n<10: return dict(label=label, n=n, note="insufficient N")
    # continuous: Spearman pcc vs EffectSize (expect NEGATIVE)
    rho,p_rho=stats.spearmanr(d[pcc_col], d["EffectSize"])
    # bootstrap CI on rho
    boots=[]
    idx=np.arange(n)
    for _ in range(2000):
        bi=rng.choice(idx,n,replace=True)
        r,_=stats.spearmanr(d[pcc_col].values[bi], d["EffectSize"].values[bi])
        boots.append(r)
    rho_ci=[float(np.percentile(boots,2.5)), float(np.percentile(boots,97.5))]
    # sign concordance
    dbn_sign=np.sign(d[pcc_col].values)
    cr_reg_sign=-np.sign(d["EffectSize"].values)  # activating(+) if EffectSize<0
    concordant=(dbn_sign==cr_reg_sign)
    obs_conc=concordant.mean()
    # base-rate expected concordance under independence
    p_dbn_pos=(dbn_sign>0).mean(); p_cr_pos=(cr_reg_sign>0).mean()
    exp_conc=p_dbn_pos*p_cr_pos + (1-p_dbn_pos)*(1-p_cr_pos)
    # 2x2 Fisher
    a=int(((dbn_sign>0)&(cr_reg_sign>0)).sum())  # both activating
    b=int(((dbn_sign>0)&(cr_reg_sign<0)).sum())
    c=int(((dbn_sign<0)&(cr_reg_sign>0)).sum())
    e=int(((dbn_sign<0)&(cr_reg_sign<0)).sum())  # both repressive
    OR,p_fish=stats.fisher_exact([[a,b],[c,e]])
    # negative-arm power: significant repressive CRISPRi (EffectSize>0 & sig)
    n_cr_repress=int((d["EffectSize"]>0).sum())
    n_dbn_neg=int((dbn_sign<0).sum())
    return dict(label=label, n=n,
        spearman_rho=float(rho), spearman_p=float(p_rho), rho_ci95=rho_ci,
        obs_concordance=float(obs_conc), expected_concordance_baserate=float(exp_conc),
        concordance_excess=float(obs_conc-exp_conc),
        fisher_OR=float(OR), fisher_p=float(p_fish),
        table_2x2=dict(dbn_pos_cr_activ=a, dbn_pos_cr_repress=b, dbn_neg_cr_activ=c, dbn_neg_cr_repress=e),
        n_cr_significant_repressive=n_cr_repress, n_dbn_negative=n_dbn_neg,
        pct_dbn_positive=float(p_dbn_pos), pct_cr_activating=float(p_cr_pos))

results={}
results["primary_pooled_sig"]=concord_stats(m_pool, "dbn_pcc", True, "Pooled DBNascent x all-CRISPRi, significant only")
results["blood_K562_sig"]=concord_stats(m_blood[m_blood["CellType"]=="K562"], "dbn_pcc", True, "DBNascent blood x K562-CRISPRi, significant (cell-type matched)")
results["pooled_all_pairs"]=concord_stats(m_pool, "dbn_pcc", False, "Pooled, all tested pairs (incl non-sig) — attenuated")

# confound controls on primary pooled significant set
base=m_pool.dropna(subset=["dbn_pcc","EffectSize"]); base=base[(base["dbn_pcc"]!=0)&(base["pValueAdjusted"]<0.05)]
results["ctrl_well_powered"]=concord_stats(base[base["PowerAtEffectSize25"]>0.8], "dbn_pcc", True, "Well-powered (PowerAtES25>0.8)")
results["ctrl_direct_effects"]=concord_stats(base[base["direct_vs_indirect_negative"]>0.8], "dbn_pcc", True, "Direct effects only (>0.8)")
results["ctrl_h3k27ac"]=concord_stats(base[base["H3K27ac_peak_overlap"]==1], "dbn_pcc", True, "H3K27ac enhancers only")
# distance-stratified Spearman
base2=base.copy(); base2["dbin"]=(base2["abs_dist"]//100000).clip(upper=10)
dstrat=[]
for b,g in base2.groupby("dbin"):
    if len(g)>=20:
        r,p=stats.spearmanr(g["dbn_pcc"],g["EffectSize"])
        dstrat.append(dict(dist_kb=int(b*100), n=len(g), rho=float(r), p=float(p)))
results["distance_stratified"]=dstrat

json.dump(results, open(OUT/"concordance_results.json","w"), indent=2, default=float)
base.to_parquet(OUT/"concordance_matched_pairs.parquet", index=False)
# print key
print("=== PRIMARY (pooled, significant) ===")
r=results["primary_pooled_sig"]
for k in ["n","spearman_rho","spearman_p","rho_ci95","obs_concordance","expected_concordance_baserate","concordance_excess","fisher_OR","fisher_p","table_2x2","n_cr_significant_repressive","n_dbn_negative"]:
    print(f"  {k}: {r.get(k)}")
print("\n=== confound controls (rho, n) ===")
for k in ["ctrl_well_powered","ctrl_direct_effects","ctrl_h3k27ac","blood_K562_sig"]:
    rr=results[k]; print(f"  {k}: n={rr.get('n')} rho={rr.get('spearman_rho')} OR={rr.get('fisher_OR')} p={rr.get('spearman_p')}")
log.info("DONE")

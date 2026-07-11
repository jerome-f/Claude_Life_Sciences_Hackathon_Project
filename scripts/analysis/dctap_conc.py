import gcsfs, duckdb, pandas as pd, numpy as np, json, logging
from pathlib import Path
from scipy import stats
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/concord"); rng=np.random.default_rng(23)
fs=gcsfs.GCSFileSystem(project="flash-hour-452305-m7",requester_pays=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='40GB'"); con.register_filesystem(fs)
DBN="gs://claude_hackathon/dbnascent/20260711/processed"

def load_dctap(path, cell):
    df=pd.read_csv(path, sep="\t")
    vc=df["ValidConnection"].astype(str)
    bad=vc.str.contains("promoter|TSS targeting|gene body|exon", case=False, regex=True)
    df=df[~bad].copy()  # clean distal only
    df["chrom"]=df["chrom"].astype("string").astype(str)
    df["CellType"]=cell; df["Dataset"]="DC_TAP"
    df["Regulated"]=(df["pValueAdjusted"]<0.05)&(df["EffectSize"]<0)
    return df[["chrom","chromStart","chromEnd","EffectSize","pValueAdjusted","Significant","Regulated","CellType","Dataset","measuredGeneSymbol","PowerAtEffectSize25"]]

k=load_dctap("/data/dctap/ENCODE_K562_DC_TAP_Seq_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz","K562")
w=load_dctap("/data/dctap/ENCODE_WTC11_DC_TAP_Seq_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz","WTC11")

# EPCrisprBenchmark (already used) for pooling
tr=pd.read_csv("/data/crispri/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz", sep="\t")
ho=pd.read_csv("/data/crispri/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz", sep="\t")
ep=pd.concat([tr,ho],ignore_index=True); ep=ep[ep["ValidConnection"]==True].copy()
ep["chrom"]=ep["chrom"].astype("string").astype(str); ep["Dataset"]="EPCrispr"
ep=ep[["chrom","chromStart","chromEnd","EffectSize","pValueAdjusted","Significant","Regulated","CellType","Dataset","measuredGeneSymbol","PowerAtEffectSize25"]]

# DBNascent per-tissue with t
con.execute(f"""CREATE TABLE dbn AS
SELECT gene_id, chrom,
  CAST(regexp_extract(bidirectional_id,':([0-9]+)-',1) AS BIGINT) es,
  CAST(regexp_extract(bidirectional_id,'-([0-9]+)$',1) AS BIGINT) ee,
  pcc, t AS tstat, nObs, tissue
FROM read_parquet('{DBN}/dbnascent_e2g_pertissue_signed.parquet')""")
ctmap={"K562":"blood","WTC11":None}  # WTC11=iPSC, no DBNascent match; use pooled all-tissue for it

def match(cr_df, tissue=None):
    con.register("cr", cr_df)
    if tissue:
        q=f"""SELECT c.*, d.pcc, d.tstat, d.nObs FROM cr c JOIN dbn d
          ON c.chrom=d.chrom AND c.measuredGeneSymbol=d.gene_id AND d.tissue='{tissue}'
          AND c.chromStart<=d.ee AND c.chromEnd>=d.es"""
    else:
        # pooled: mean pcc / stouffer t across tissues
        q="""SELECT c.*, z.pcc, z.tstat, z.nObs FROM cr c JOIN
          (SELECT gene_id,chrom,es,ee, avg(pcc) pcc, sum(sign(tstat)*abs(tstat))/sqrt(count(*)) tstat, sum(nObs) nObs
           FROM dbn GROUP BY 1,2,3,4) z
          ON c.chrom=z.chrom AND c.measuredGeneSymbol=z.gene_id AND c.chromStart<=z.ee AND c.chromEnd>=z.es"""
    return con.execute(q).df()

def stat_block(d, name):
    d=d.dropna(subset=["pcc","EffectSize"]); d=d[d["pcc"]!=0]; d=d[d["pValueAdjusted"]<0.05]
    n=len(d)
    if n<8: return dict(label=name, n=int(n), note="insufficient")
    rho,p=stats.spearmanr(d["pcc"], d["EffectSize"])
    ds=np.sign(d["pcc"].values); cs=-np.sign(d["EffectSize"].values)
    a=int(((ds>0)&(cs>0)).sum()); b=int(((ds>0)&(cs<0)).sum()); c=int(((ds<0)&(cs>0)).sum()); e=int(((ds<0)&(cs<0)).sum())
    OR,pf=stats.fisher_exact([[a,b],[c,e]]) if min(a,b,c,e)>=0 and (b+c)>0 else (np.nan,np.nan)
    out=dict(label=name, n=int(n), spearman_rho=float(rho), spearman_p=float(p),
        n_sig_repressive=int((d["EffectSize"]>0).sum()), n_dbn_neg=int((ds<0).sum()),
        fisher_OR=float(OR), fisher_p=float(pf), table=dict(pp=a,pn=b,np=c,nn=e))
    reg=d[d["Regulated"]==True]
    if len(reg)>=10:
        mr,mp=stats.spearmanr(reg["pcc"].abs(), reg["EffectSize"].abs())
        mt,mtp=stats.spearmanr(reg["tstat"].abs(), reg["EffectSize"].abs())
        out.update(mag_abs_pcc_rho=float(mr), mag_abs_pcc_p=float(mp), mag_abs_t_rho=float(mt), mag_abs_t_p=float(mtp), n_regulated=int(len(reg)))
    return out

res={}
# 1) DC-TAP K562 x DBNascent blood (independent replication, cell-type matched)
res["dctap_k562_blood"]=stat_block(match(k,"blood"), "DC-TAP K562 x DBNascent blood (independent replication)")
# 2) DC-TAP WTC11 x DBNascent pooled (iPSC, no tissue match)
res["dctap_wtc11_pooled"]=stat_block(match(w,None), "DC-TAP WTC11 x DBNascent pooled")
# 3) POOLED EPCrispr + DC-TAP, K562-relevant x blood
pooled=pd.concat([ep[ep["CellType"]=="K562"], k], ignore_index=True)
res["pooled_ep_dctap_k562_blood"]=stat_block(match(pooled,"blood"), "EPCrispr+DC-TAP K562 x blood (pooled)")
# 4) FULL pool: all EPCrispr + both DC-TAP, x pooled DBNascent
fullpool=pd.concat([ep,k,w], ignore_index=True)
res["full_pool"]=stat_block(match(fullpool,None), "All CRISPRi (EPCrispr+DC-TAP) x DBNascent pooled")

json.dump(res, open(OUT/"concordance_dctap.json","w"), indent=2, default=float)
print(json.dumps(res, indent=2, default=float))

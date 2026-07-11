import pandas as pd, numpy as np, json, logging
from pathlib import Path
from scipy import stats
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s"); log=logging.getLogger()
OUT=Path("/data/concord"); rng=np.random.default_rng(7)
m=pd.read_parquet(OUT/"concordance_matched_pairs.parquet")  # pooled significant set
log.info(f"loaded matched significant pairs: {len(m)}")

d=m.dropna(subset=["dbn_pcc","EffectSize"]).copy(); d=d[d["dbn_pcc"]!=0]

# 1) PERMUTATION test: is Spearman rho more extreme than under label shuffle?
obs_rho,_=stats.spearmanr(d["dbn_pcc"],d["EffectSize"])
perm=np.empty(5000)
ev=d["EffectSize"].values.copy()
for i in range(5000):
    perm[i],_=stats.spearmanr(d["dbn_pcc"].values, rng.permutation(ev))
p_perm=(np.sum(perm<=obs_rho)+1)/(5001)  # one-sided (expect negative)
log.info(f"permutation p (rho<=obs): {p_perm}")

# 2) PARTIAL correlation controlling for |distance| (Spearman partial via rank residuals)
def partial_spearman(x,y,z):
    rx=stats.rankdata(x); ry=stats.rankdata(y); rz=stats.rankdata(z)
    # residualize rx and ry on rz
    def resid(a,b):
        b=np.c_[np.ones_like(b),b]; beta,_,_,_=np.linalg.lstsq(b,a,rcond=None); return a-b@beta
    return stats.pearsonr(resid(rx,rz),resid(ry,rz))
pr,pp=partial_spearman(d["dbn_pcc"].values,d["EffectSize"].values,d["abs_dist"].values)
log.info(f"partial rho|distance = {pr:.4f} (p={pp:.2e})")

# 3) magnitude-magnitude: does |pcc| track |EffectSize| among regulated?
reg=d[d["Regulated"]==True]
if len(reg)>20:
    rho_reg,p_reg=stats.spearmanr(reg["dbn_pcc"].abs(), reg["EffectSize"].abs())
else:
    rho_reg,p_reg=None,None

# 4) sign concordance restricted to pairs where BOTH signals are confident
#    DBNascent |pcc|>0.3, CRISPRi significant & regulated
conf=d[(d["dbn_pcc"].abs()>0.3)]
ds=np.sign(conf["dbn_pcc"].values); cs=-np.sign(conf["EffectSize"].values)
a=int(((ds>0)&(cs>0)).sum()); b=int(((ds>0)&(cs<0)).sum())
c=int(((ds<0)&(cs>0)).sum()); e=int(((ds<0)&(cs<0)).sum())
OR_conf,p_conf=stats.fisher_exact([[a,b],[c,e]]) if (b+c)>0 else (np.inf,0)

# 5) per-cell-type Spearman (all matched, sig)
allm=pd.read_parquet(OUT/"concordance_matched_pairs.parquet")
byct=[]
for ct,g in d.groupby("CellType"):
    if len(g)>=20:
        r,p=stats.spearmanr(g["dbn_pcc"],g["EffectSize"])
        byct.append(dict(cell_type=str(ct), n=int(len(g)), rho=float(r), p=float(p)))

extra=dict(
  permutation_p_onesided=float(p_perm), obs_rho=float(obs_rho),
  partial_rho_given_distance=float(pr), partial_p=float(pp),
  magnitude_rho_regulated=(float(rho_reg) if rho_reg is not None else None),
  magnitude_p_regulated=(float(p_reg) if p_reg is not None else None), n_regulated=int(len(reg)),
  confident_signs_2x2=dict(n=int(len(conf)), dbn_pos_cr_activ=a, dbn_pos_cr_repress=b, dbn_neg_cr_activ=c, dbn_neg_cr_repress=e, fisher_OR=float(OR_conf), fisher_p=float(p_conf)),
  per_cell_type=byct)
json.dump(extra, open(OUT/"concordance_extra.json","w"), indent=2, default=float)
print(json.dumps(extra, indent=2, default=float))

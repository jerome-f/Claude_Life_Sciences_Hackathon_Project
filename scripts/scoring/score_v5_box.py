#!/usr/bin/env python3
"""v5 utility scoring on box. Fetches raw 1024-d embeddings, projects to PCA-613 (v4 basis),
loads v5, computes LINK metrics + gene-centric prioritization + SIGN. Prints RESULT_JSON."""
import os, json, time, subprocess, numpy as np, torch
import polars as pl
D="."
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)

# --- fetch raw embeddings via signed URLs ---
urls=json.load(open("get2.json"))
for local,u in urls.items():
    if os.path.exists(local) and os.path.getsize(local)>1000: 
        log(local,"cached"); continue
    r=subprocess.run(["curl","-s","-o",local,"-w","%{http_code}","--max-time","500",u],capture_output=True,text=True)
    log("fetch",local,r.stdout.strip(),os.path.getsize(local) if os.path.exists(local) else 0)

# --- reconstruct raw, project to PCA-613 (v4 basis) ---
comp=np.load("pca_components_v4.npy")   # (613,1024)
pmean=np.load("pca_mean_v4.npy")        # (1024,)
hmean=np.load("human_mean.npy")         # (1024,) human common-mode
Eh_c=np.load("E_centered.npy", mmap_mode="r")   # human centered (N,1024) f4
Em_r=np.load("E_mouse_raw.npy", mmap_mode="r")  # mouse raw (M,1024) f2
log(f"Eh_c {Eh_c.shape} Em_r {Em_r.shape} comp {comp.shape}")
def project(E, is_human):
    E=np.asarray(E, dtype=np.float32)
    if is_human: E=E + hmean.astype(np.float32)   # restore raw human
    return ((E - pmean.astype(np.float32)) @ comp.T.astype(np.float32)).astype(np.float32)  # (N,613)
Eh=project(Eh_c, True); Em=project(Em_r, False)
log(f"projected PCA-613: human {Eh.shape} mouse {Em.shape}")
ih=json.load(open("id_index_human.json")); im=json.load(open("id_index_mouse.json"))

# --- HELD-OUT test pairs: read the v4 TEST parquet ONLY (split_v4=='test' by construction).
# NOTE: chr3/chr13 are v5 TRAIN chromosomes; the held-out test split is
# human {chr16,chr2,chr22,chr9} / mouse {chr2,chr9,chrX}. Using the test parquet avoids leakage.
rows=[]
for sp,idx,pfx in (("human",ih,"h:"),("mouse",im,"m:")):
    df=pl.read_parquet(f"v4_pairs_{sp}_test.parquet")
    for r in df.iter_rows(named=True):
        eid=r["elem_id"]; gid=r["gene_id"]
        if eid in idx and gid in idx:
            rows.append((pfx, idx[eid], idx[gid], eid, gid,
                         float(r.get("abs_distance") or 0), int(r.get("consensus_sign") or 0),
                         int(r["link"]), r.get("label_route","")))
log(f"test pairs resolved: {len(rows)}")

# --- features (match v5 build_pair_features) ---
CONTACT_D0=1000.0; CONTACT_GAMMA=1.0; N=len(rows)
E_arr=np.zeros((N,613),np.float32); G_arr=np.zeros((N,613),np.float32); S_arr=np.zeros((N,3),np.float32)
for i,(pfx,ei,gi,eid,gid,dist,csign,link,route) in enumerate(rows):
    E=Eh if pfx=="h:" else Em
    e=E[ei]; g=E[gi]
    contact=1.0/(1.0+(dist/CONTACT_D0)**CONTACT_GAMMA); abc=contact*float(np.linalg.norm(e))
    E_arr[i]=e; G_arr[i]=g; S_arr[i]=[dist/1000.0, contact, abc]

import v5_model as v5m
net=v5m.TwoTowerE2G(emb_dim=613,n_scalar=3)
ck=torch.load("v5_twotower.pt",map_location="cpu"); sd=ck.get("state",ck)
mi,un=net.load_state_dict(sd,strict=False); net.eval()
DEV="cuda" if torch.cuda.is_available() else "cpu"; net.to(DEV)
log(f"v5 loaded on {DEV} missing={len(mi)} unexpected={len(un)}")
link_p=np.zeros(N); sign_rep_p=np.zeros(N); BS=8192
with torch.no_grad():
    for b in range(0,N,BS):
        e=torch.from_numpy(E_arr[b:b+BS]).to(DEV); g=torch.from_numpy(G_arr[b:b+BS]).to(DEV); s=torch.from_numpy(S_arr[b:b+BS]).to(DEV)
        o=net(e,g,s)
        link_p[b:b+BS]=torch.sigmoid(o["link"]).cpu().numpy()
        sign_rep_p[b:b+BS]=torch.softmax(o["sign"],dim=1).cpu().numpy()[:,1]
log("scored")

from sklearn.metrics import average_precision_score, roc_auc_score
y=np.array([r[7] for r in rows]); res={"n_test":int(N),"link_prev":round(float(y.mean()),4)}
res["link"]={"AUPRC":round(float(average_precision_score(y,link_p)),4),"AUROC":round(float(roc_auc_score(y,link_p)),4)}
mask=np.array([(r[7]==1 and r[6]!=0 and r[8]=="confident") for r in rows])
if mask.sum()>10:
    yr=np.array([1 if r[6]==-1 else 0 for r in rows])[mask]
    if len(set(yr.tolist()))>1:
        res["sign_repressive"]={"AUPRC":round(float(average_precision_score(yr,sign_rep_p[mask])),4),
            "AUROC":round(float(roc_auc_score(yr,sign_rep_p[mask])),4),"n":int(mask.sum()),"rep_prev":round(float(yr.mean()),4)}
from collections import defaultdict
byg=defaultdict(list)
for i,r in enumerate(rows): byg[(r[0],r[4])].append(i)
p1=[]; rr=[]; dp1=[]
for k,idxs in byg.items():
    if len(idxs)<2: continue
    ys=np.array([rows[i][7] for i in idxs])
    if ys.sum()<1: continue
    ps=link_p[idxs]; ds=np.array([rows[i][5] for i in idxs]); order=np.argsort(-ps)
    p1.append(1.0 if rows[idxs[order[0]]][7]==1 else 0.0)
    ranks=np.where(ys[order]==1)[0]; rr.append(1.0/(ranks[0]+1) if len(ranks) else 0.0)
    dp1.append(1.0 if rows[idxs[int(np.argmin(ds))]][7]==1 else 0.0)
res["gene_prioritization"]={"n_genes":len(p1),"model_precision_at_1":round(float(np.mean(p1)),4),
    "model_MRR":round(float(np.mean(rr)),4),"distance_precision_at_1":round(float(np.mean(dp1)),4)}
print("RESULT_JSON="+json.dumps(res))

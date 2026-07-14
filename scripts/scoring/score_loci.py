#!/usr/bin/env python3
"""End-to-end locus demonstration: at each benchmark GWAS locus, use v5 signed-E2G to nominate
the effector gene, and compare to the curated effector, stratified by cis-eQTL coloc status.

Anchor: nearest DBNascent enhancer element to the sentinel variant (same chrom). For each candidate
gene (HGNC symbol) present in our TSS index, score LINK(anchor_enhancer -> gene_TSS). Top-ranked
gene = E2G nomination. distance baseline = candidate gene whose TSS is nearest the anchor.
Reports, per coloc_status, the fraction of loci where the E2G top gene == curated effector."""
import os, json, time, re, numpy as np, torch
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]",*a,flush=True)

comp=np.load("pca_components_v4.npy"); pmean=np.load("pca_mean_v4.npy"); hmean=np.load("human_mean.npy")
Eh_c=np.load("E_centered.npy", mmap_mode="r")
def project(E): 
    E=np.asarray(E,dtype=np.float32)+hmean.astype(np.float32)
    return ((E-pmean.astype(np.float32))@comp.T.astype(np.float32)).astype(np.float32)
ih=json.load(open("id_index_human.json"))
tss_coords=json.load(open("tss_coords.json"))   # symbol -> {chrom, mid, strand}
log("loaded index", len(ih), "tss_coords", len(tss_coords))

# parse element coords (chr:start-end) and gene symbols
elem_re=re.compile(r'(chr[\w]+):(\d+)-(\d+)')
elems_by_chr={}   # chrom -> list of (mid, key)
gene_rows={}      # symbol -> row idx
for k,row in ih.items():
    m=elem_re.match(k)
    if m:
        c=m.group(1); mid=(int(m.group(2))+int(m.group(3)))//2
        elems_by_chr.setdefault(c,[]).append((mid,k,row))
    else:
        gene_rows[k]=row
for c in elems_by_chr: elems_by_chr[c].sort()
log("elements chroms", len(elems_by_chr), "genes", len(gene_rows))

import v5_model as v5m
net=v5m.TwoTowerE2G(emb_dim=613,n_scalar=3)
ck=torch.load("v5_twotower.pt",map_location="cpu"); net.load_state_dict(ck.get("state",ck),strict=False); net.eval()
log("v5 loaded")

loci=json.load(open("benchmark_loci.json"))
import bisect
CONTACT_D0=1000.0
MAXA=int(os.environ.get("MAX_ANCHOR_KB","50"))*1000   # require anchor within this of sentinel
results=[]; skipped_noanchor=0; skipped_nocand=0
# cache projected embeddings for rows we touch
row_cache={}
def emb(row):
    if row not in row_cache: row_cache[row]=project(Eh_c[row:row+1])[0]
    return row_cache[row]

for L in loci:
    chrom="chr"+str(L["chr"]); pos=int(L["pos_hg38"])
    cands=[g for g in str(L["candidate_genes"]).split(";") if g in gene_rows]
    if not cands: skipped_nocand+=1; continue
    arr=elems_by_chr.get(chrom)
    if not arr: skipped_noanchor+=1; continue
    mids=[a[0] for a in arr]; j=bisect.bisect_left(mids,pos)
    best=None
    for jj in (j-1,j,j+1):
        if 0<=jj<len(arr):
            d=abs(arr[jj][0]-pos)
            if best is None or d<best[0]: best=(d,arr[jj])
    if best is None or best[0]>MAXA: skipped_noanchor+=1; continue
    anch_mid,anch_key,anch_row=best[1]
    e=emb(anch_row)
    # score anchor -> each candidate TSS
    scores=[]; gdists=[]
    E_batch=[]; G_batch=[]; S_batch=[]
    for g in cands:
        gv=emb(gene_rows[g])
        # real enhancer(anchor)->gene TSS distance, matching v5 training features
        tc=tss_coords.get(g)
        if tc and tc["chrom"]==chrom: dist=abs(anch_mid - tc["mid"])
        else: dist=abs(anch_mid - pos)   # fallback: anchor->sentinel
        gdists.append(dist)
        dist_kb=dist/1000.0; contact=1.0/(1.0+(dist/CONTACT_D0)**1.0); abc=contact*float(np.linalg.norm(e))
        E_batch.append(e); G_batch.append(gv); S_batch.append([dist_kb,contact,abc])
    with torch.no_grad():
        o=net(torch.from_numpy(np.array(E_batch,dtype=np.float32)),
              torch.from_numpy(np.array(G_batch,dtype=np.float32)),
              torch.from_numpy(np.array(S_batch,dtype=np.float32)))
        scores=torch.sigmoid(o["link"]).numpy().tolist()
    top=cands[int(np.argmax(scores))]
    # distance-only baseline: nearest-TSS candidate gene
    dist_top=cands[int(np.argmin(gdists))]
    eff=str(L["known_effector_gene"])
    results.append({"locus":L["locus_id"],"coloc":L["coloc_status"],"eff":eff,"conf":L.get("effector_confidence"),
                    "top":top,"hit":int(top==eff),"dist_top":dist_top,"dist_hit":int(dist_top==eff),
                    "n_cand":len(cands),"anchor_kb":round(best[0]/1000,1),"eff_in_cand":int(eff in cands)})
log(f"scored {len(results)} loci; skipped noanchor={skipped_noanchor} nocand={skipped_nocand}")

# metrics stratified by coloc status (only loci where the effector is among scored candidates -> fair p@1)
from collections import defaultdict
def rate(rows): 
    rr=[r for r in rows if r["eff_in_cand"]==1]
    if not rr: return (0,None,None)
    return (len(rr), round(np.mean([r["hit"] for r in rr]),4), round(np.mean([r["dist_hit"] for r in rr]),4))
out={"n_scored":len(results),"skipped_noanchor":skipped_noanchor,"skipped_nocand":skipped_nocand,
     "max_anchor_kb":MAXA//1000,"by_coloc":{}}
for st in ("SILENT","WRONG","POSITIVE"):
    sub=[r for r in results if r["coloc"]==st]
    n,p,dp=rate(sub); out["by_coloc"][st]={"n_loci":len(sub),"n_eff_in_cand":n,"e2g_precision_at_1":p,"distance_precision_at_1":dp}
n_all,p_all,dp_all=rate(results); out["overall"]={"n_eff_in_cand":n_all,"e2g_precision_at_1":p_all,"distance_precision_at_1":dp_all}
# random baseline = 1/median candidates
med=int(np.median([r["n_cand"] for r in results])) if results else 0
out["random_baseline_p1"]=round(1.0/med,4) if med else None
out["median_candidates"]=med
json.dump(results, open("loci_results.json","w"))
print("LOCI_JSON="+json.dumps(out))

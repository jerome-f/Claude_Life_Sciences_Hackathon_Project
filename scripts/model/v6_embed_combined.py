#!/usr/bin/env python3
"""
v6 COMBINED per-token embed job (human hg38 + mouse mm10), Option A substrate.

Two-pass, writes ONLY the reduced codec (never the 414 GB full array):
  PASS 1  fit a token-level PCA 1024->CODEC_DIM on a random sample of tokens
          (frozen NT last-hidden-states), pooled across a window sample.
  PASS 2  embed ALL combined windows, project each token 1024->CODEC_DIM on the GPU,
          write (N, L_FIX, CODEC_DIM) fp16 memmap + token-length vector + id_index.

Conventions locked (match v4/v5 so tokens share the embedding manifold):
  * FORWARD reference strand for ALL windows (NO reverse-complement) — preserves the
    shared human/mouse space and is correct for a DIRECTIONAL task (RC-consistency hurts
    orientation-dependent tasks; strand enters as a FEATURE in the trainer, not here).
  * fp32 NT weights + bf16 autocast (never hard-cast NT weights: breaks rotary buffers).

CODEC_DIM is a STORAGE codec only; the model's representation is the learned projection
on top of it (v6_train). id keys are species-prefixed ('h:'/'m:') to stay unique.
"""
import os, sys, json, time, gzip, urllib.request
import numpy as np

MODEL_ID  = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
L_FIX     = 168
H         = 1024
CODEC_DIM = int(os.environ.get("V6_CODEC_DIM", "256"))
HB = "embed_hb.json"

GENOME_URL = {"human": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
              "mouse": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"}
GENOME_FA  = {"human": "hg38.fa", "mouse": "mm10.fa"}

def hb(phase, frac=0.0, **extra):
    d = {"phase": phase, "frac": round(frac, 4), "t": time.time()}; d.update(extra)
    json.dump(d, open(HB, "w")); print(f"[HB] {phase} {frac:.3f} {extra}", flush=True)

def ensure_fasta(species):
    fa = GENOME_FA[species]
    if os.path.exists(fa) and os.path.getsize(fa) > 2_000_000_000:
        print(f"{fa} present", flush=True); return fa
    print(f"downloading {species} genome ...", flush=True); t0 = time.time()
    urllib.request.urlretrieve(GENOME_URL[species], fa + ".gz")
    with gzip.open(fa + ".gz", "rt") as fi, open(fa, "w") as fo:
        for line in fi: fo.write(line)
    print(f"  {fa} ready in {time.time()-t0:.0f}s", flush=True); return fa

def load_genome(fa, keep):
    keep = set(keep); genome = {}; name = None; buf = []
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None and name in keep: genome[name] = "".join(buf)
                name = line[1:].strip().split()[0]; buf = []
            elif name in keep:
                buf.append(line.strip())
    if name is not None and name in keep: genome[name] = "".join(buf)
    return {c: genome[c].upper() for c in genome}

def extract_seq(genome, chrom, start, end):
    s = genome[chrom][start:end]
    if len(s) < (end - start): s = s + "N" * ((end - start) - len(s))
    return s

def load_model():
    import torch
    from transformers import AutoTokenizer, AutoModelForMaskedLM
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    tok = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(
        MODEL_ID, output_hidden_states=True, trust_remote_code=True).to("cuda").eval()
    return tok, model

def embed_tokens(tok, model, seqs):
    """-> (b, L_FIX, H) fp32 hidden states + (b, L_FIX) attention mask."""
    import torch
    enc = tok(seqs, return_tensors="pt", padding="max_length", truncation=True,
              max_length=L_FIX).to("cuda")
    with torch.inference_mode(), torch.autocast("cuda", dtype=torch.bfloat16):
        o = model(**enc)
    return o.hidden_states[-1].float(), enc["attention_mask"]

def main():
    import torch, polars as pl
    from sklearn.decomposition import PCA
    win = pl.read_parquet("v6_windows_combined.parquet")
    genomes = {}
    for sp in ("human", "mouse"):
        chroms = sorted(win.filter(pl.col("species") == sp)["chrom"].unique().to_list())
        hb(f"genome_{sp}"); genomes[sp] = load_genome(ensure_fasta(sp), chroms)
    tok, model = load_model()
    rows = win.to_dicts()
    N = len(rows)
    print(f"combined windows={N} codec_dim={CODEC_DIM}", flush=True)

    # ---- PASS 1: fit token-level PCA on a window sample ----
    hb("pca_sample")
    rng = np.random.RandomState(0)
    samp = rng.choice(N, size=min(4000, N), replace=False)
    tok_acc = []
    BS = 256
    for s in range(0, len(samp), BS):
        chunk = [rows[i] for i in samp[s:s+BS]]
        seqs = [extract_seq(genomes[r["species"]], r["chrom"], r["win_start"], r["win_end"]) for r in chunk]
        h, am = embed_tokens(tok, model, seqs)          # (b,L,H)
        m = am.bool()
        # take valid tokens only, subsample to keep memory bounded
        for j in range(h.shape[0]):
            v = h[j][m[j]].cpu().numpy()
            if len(v) > 40: v = v[rng.choice(len(v), 40, replace=False)]
            tok_acc.append(v)
    X = np.concatenate(tok_acc, 0).astype(np.float32)
    print(f"PCA fit on {X.shape[0]} tokens x {X.shape[1]}", flush=True)
    pca = PCA(n_components=CODEC_DIM, svd_solver="randomized", random_state=0).fit(X)
    evr = float(pca.explained_variance_ratio_.sum())
    np.save("v6_codec_mean.npy", pca.mean_.astype(np.float32))
    np.save("v6_codec_components.npy", pca.components_.astype(np.float32))  # (CODEC_DIM, H)
    print(f"codec EVR@{CODEC_DIM}={evr:.4f}", flush=True)
    Wc = torch.from_numpy(pca.components_.astype(np.float32)).cuda()          # (C,H)
    mu = torch.from_numpy(pca.mean_.astype(np.float32)).cuda()                # (H,)

    # ---- PASS 2: embed all, project on GPU, write reduced memmap ----
    E = np.lib.format.open_memmap("E_combined_pertoken.npy", mode="w+",
                                  dtype=np.float16, shape=(N, L_FIX, CODEC_DIM))
    lens = np.zeros(N, dtype=np.int16)
    t0 = time.time()
    for i in range(0, N, BS):
        chunk = rows[i:i+BS]
        seqs = [extract_seq(genomes[r["species"]], r["chrom"], r["win_start"], r["win_end"]) for r in chunk]
        h, am = embed_tokens(tok, model, seqs)          # (b,L,H)
        red = torch.matmul(h - mu, Wc.T)                # (b,L,C)
        E[i:i+len(chunk)] = red.to(torch.float16).cpu().numpy()
        lens[i:i+len(chunk)] = am.sum(1).to(torch.int16).cpu().numpy()
        if (i // BS) % 40 == 0:
            wps = (i + len(chunk)) / max(time.time() - t0, 1e-9)
            hb("embed", (i + len(chunk)) / N, wps=round(wps, 1),
               eta_min=round((N - i) / max(wps, 1e-9) / 60, 1))
    E.flush()
    np.save("E_combined_pertoken_lens.npy", lens)
    json.dump({r["uid"]: k for k, r in enumerate(rows)}, open("id_index_combined.json", "w"))
    perf = {"n": N, "L_fix": L_FIX, "codec_dim": CODEC_DIM, "codec_evr": round(evr, 4),
            "wall_s": round(time.time()-t0, 1), "gb": round(E.nbytes/1e9, 1)}
    json.dump(perf, open("embed_perf.json", "w"))
    hb("done", 1.0, **perf)
    print("COMBINED PERTOKEN DONE", perf, flush=True)

if __name__ == "__main__":
    main()

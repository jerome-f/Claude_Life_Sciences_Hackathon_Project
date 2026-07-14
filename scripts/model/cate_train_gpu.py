#!/usr/bin/env python3
"""
v6 trainer (Option A per-token cross-encoder). Rough-and-dirty subset run:
train on 12 chroms, val on 2, test on 2 (all leakage-safe subsets of the canonical
v4/v5 partitions). bf16 everywhere + torch.compile. Loads the per-token memmap
(N,168,1024) fp16, builds padding masks from token lengths, gathers per batch.

Reuses v5's label logic (build_pair_features), interaction/heads (via v6_model
CrossEncoderE2G), and masked multitask_loss verbatim so v6 vs v5 is controlled.
"""
import os, json, time
import numpy as np
import torch
import polars as pl
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.isotonic import IsotonicRegression
from cate_model import CATE                       # CNN -> cross-attention hybrid tower
from v6_model import multitask_loss, build_pair_features, SIGN_REP

D = os.environ.get("E2G_DIR", ".")
DEV = "cuda" if torch.cuda.is_available() else "cpu"
L_FIX = 168
CODEC_DIM = int(os.environ.get("V6_CODEC_DIM", "64"))   # per-token storage codec width
def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ---- best-effort GCS checkpoint PUT (signed URL from v6_urls.json) ----
_PUT_URLS = None
def _put_gcs(local_path, key):
    """PUT a file to GCS via a signed URL for `key`. Never raises — resilience only.
    Reads cate_put_urls.json (fresh virtual-hosted-style V4 signed PUTs bound to the
    HYPHEN bucket claude-hackathon). The old v6_urls.json used the underscore-bucket
    path and 404-failed on the epoch-1 checkpoint; this file is the fix."""
    global _PUT_URLS
    try:
        import urllib.request
        if _PUT_URLS is None:
            _PUT_URLS = json.load(open(f"{D}/cate_put_urls.json"))
        url = _PUT_URLS.get(key)
        if not url:
            log(f"[ckpt] no PUT url for {key} (skipping GCS, local save intact)")
            return
        with open(local_path, "rb") as f:
            req = urllib.request.Request(url, data=f.read(), method="PUT",
                                         headers={"Content-Type": "application/octet-stream"})
            r = urllib.request.urlopen(req, timeout=180)
        log(f"[ckpt] PUT {key} -> GCS {r.status}")
    except Exception as e:
        log(f"[ckpt] GCS PUT failed for {key}: {e} (continuing; local save intact)")

# ---- combined per-token codec (memmap) + token lengths + id_index (species-prefixed uids) ----
def load_pertoken():
    # GPU-RESIDENT codec (the bottleneck fix). At CODEC_DIM=128 the array is
    # (N,168,128) fp16 ~= 52 GB and FITS a single 80 GB A100 alongside the model
    # + activations (measured ~61 GB total at bs=8192, ~19 GB headroom). Loading the
    # codec ONTO the GPU turns the per-batch gather into a GPU index (microseconds),
    # eliminating the single-threaded cache-cold CPU gather that was 90% of every batch.
    # No prefetch loader, no host->device copy of embeddings, no fork/COW.
    # Load via a CPU staging array (streamed, never a full RAM+GPU double), then to GPU.
    E_np = np.load(f"{D}/E_combined_pertoken.npy", mmap_mode="r")     # header + lazy pages
    gb = E_np.nbytes / 1e9
    log(f"combined codec E {E_np.shape} {E_np.dtype} ({gb:.1f} GB) -> loading GPU-RESIDENT on {DEV}")
    # move to GPU in row-chunks so we never hold a second full copy in host RAM
    E_gpu = torch.empty(E_np.shape, dtype=torch.float16, device=DEV)
    CH = 100000
    for s in range(0, E_np.shape[0], CH):
        e = min(E_np.shape[0], s + CH)
        E_gpu[s:e] = torch.from_numpy(np.ascontiguousarray(E_np[s:e]))
    del E_np
    log(f"codec resident on {DEV}: {E_gpu.element_size()*E_gpu.nelement()/1e9:.1f} GB, "
        f"cuda_alloc={torch.cuda.memory_allocated()/1e9:.1f} GB")
    lens = np.load(f"{D}/E_combined_pertoken_lens.npy")               # (N,) int16  (stays CPU, tiny)
    idx = json.load(open(f"{D}/id_index_combined.json"))             # uid ('h:'/'m:' prefixed) -> row
    return E_gpu, lens, idx

# ---- split membership from human+mouse v4 pair tables, by chromosome ----
def build_split(split, idx, manifest):
    chroms = set(manifest[f"{split}_chroms"])
    frames = []
    for sp, pfx, fdr in (("human", "h:", "min_fdr"), ("mouse", "m:", "min_adj_p")):
        for s in ("train", "val", "test"):
            t = pl.read_parquet(f"{D}/v4_pairs_{sp}_{s}.parquet").filter(pl.col("chrom").is_in(list(chroms)))
            if len(t) == 0:
                continue
            # species-prefix the window ids so they key into the combined index
            t = t.with_columns([(pl.lit(pfx) + pl.col("elem_id")).alias("elem_id"),
                                 (pl.lit(pfx) + pl.col("gene_id")).alias("gene_id")])
            # normalize shared-column dtypes: human/mouse v4 tables disagree on widths
            # (consensus_sign Int32 vs Int64, abs_distance Float64 vs Int64), which breaks
            # diagonal concat. Cast to canonical types before stacking.
            casts = []
            if "consensus_sign" in t.columns:
                casts.append(pl.col("consensus_sign").cast(pl.Int64))
            if "abs_distance" in t.columns:
                casts.append(pl.col("abs_distance").cast(pl.Float64))
            if casts:
                t = t.with_columns(casts)
            frames.append(t)
    tbl = pl.concat(frames, how="diagonal")
    F = build_pair_features(tbl, np.zeros((len(idx), 1), np.float32), idx)
    # strand as a 4th scalar (directional task: orientation is a FEATURE, not a symmetry) —
    # +1 for '+', -1 for '-', 0 for '.' (bidirectional elements). Aligned to resolved rows.
    strand_map = {"+": 1.0, "-": -1.0, ".": 0.0}
    strand_all = np.array([strand_map.get(s, 0.0) for s in tbl["strand"].to_list()], np.float32)
    F["strand"] = strand_all[F["keep_idx"]] if "keep_idx" in F else strand_all[:len(F["e_rows"])]
    log(f"split={split} chroms={sorted(chroms)} pairs={len(F['e_rows'])} miss={F['n_missing']} link_prev={F['y_link'].mean():.3f}")
    return F

# Finalize the scalar vector to [dist_kb, contact_freq, strand] (n_scalar=3).
# ABC (contact_freq * element-L2-activity) is DROPPED: the v4 contact-feature ablation
# showed it adds ~0 lift (human GBDT +0.0001, mouse +0.001) because DBNascent's empirical
# distance-decay is near-flat (gamma~0.05), so ABC carries nothing beyond dist_kb. Removing
# it also removes the expensive per-token activity gather. build_pair_features still emits a
# 3-wide [dist_kb, contact, abc]; we replace column 2 (abc) with strand.
def fix_activity_scalar(F, E, lens):
    strand = F.get("strand")
    if strand is None:
        strand = np.zeros(len(F["e_rows"]), np.float32)
    F["scalars"] = np.stack([F["scalars"][:, 0],          # dist_kb
                             F["scalars"][:, 1],          # contact_freq
                             strand.astype(np.float32)],  # strand (directional feature)
                            axis=1).astype(np.float32)
    return F

def make_pad(lens_batch):
    # True = padding
    ar = np.arange(L_FIX)[None, :]
    return torch.from_numpy(ar >= lens_batch[:, None])

def gather(F, E, lens, idx_rows):
    ei = F["e_rows"][idx_rows]; gi = F["g_rows"][idx_rows]
    # E is GPU-RESIDENT (torch fp16 on DEV). The gather is now a GPU index_select:
    # microseconds, no CPU work, no H2D copy of embeddings. Keep fp16 here; the model's
    # autocast(bf16) handles precision inside. This is the whole point of the redesign.
    ei_t = torch.as_tensor(ei, dtype=torch.long, device=DEV)
    gi_t = torch.as_tensor(gi, dtype=torch.long, device=DEV)
    e = E.index_select(0, ei_t)                                 # (b,168,C) fp16 on GPU
    g = E.index_select(0, gi_t)
    return {
        "e_tok": e, "g_tok": g,
        "e_pad": make_pad(lens[ei]), "g_pad": make_pad(lens[gi]),
        "scalars": torch.from_numpy(F["scalars"][idx_rows]),
        "y_link": torch.from_numpy(F["y_link"][idx_rows]),
        "y_sign": torch.from_numpy(F["y_sign"][idx_rows]),
        "mask_sign": torch.from_numpy(F["mask_sign"][idx_rows]),
        "y_elem": torch.from_numpy(F["y_elem"][idx_rows]),
        "mask_elem": torch.from_numpy(F["mask_elem"][idx_rows]),
        "y_splice": torch.from_numpy(F["y_splice"][idx_rows]),
        "mask_splice": torch.from_numpy(F["mask_splice"][idx_rows]),
    }

def to_dev(b):
    return {k: (v.to(DEV, non_blocking=True) if torch.is_tensor(v) else v) for k, v in b.items()}

@torch.no_grad()
def evaluate(net, F, E, lens, bs=2048):
    net.eval(); n = len(F["e_rows"]); yl, pl_, ys, ps = [], [], [], []
    for s in range(0, n, bs):
        idx = np.arange(s, min(s+bs, n))
        b = to_dev(gather(F, E, lens, idx))
        with torch.autocast("cuda", dtype=torch.bfloat16, enabled=(DEV=="cuda")):
            o = net(b["e_tok"], b["g_tok"], b["e_pad"], b["g_pad"], b["scalars"])
        yl.append(b["y_link"].cpu().numpy()); pl_.append(torch.sigmoid(o["link"]).float().cpu().numpy())
        ms = b["mask_sign"]
        if ms.any():
            sl = o["sign"][ms][:, [0, 1]]; p_rep = torch.softmax(sl.float(), -1)[:, 1]
            ys.append((b["y_sign"][ms] == SIGN_REP).float().cpu().numpy()); ps.append(p_rep.cpu().numpy())
    yl = np.concatenate(yl); pl_ = np.concatenate(pl_)
    res = {"link": {"AUPRC": round(float(average_precision_score(yl, pl_)), 4),
                    "AUROC": round(float(roc_auc_score(yl, pl_)), 4)}}
    if ys:
        ys = np.concatenate(ys); ps = np.concatenate(ps)
        if ys.sum() > 5 and len(set(ys.tolist())) > 1:
            res["sign_repressive"] = {"AUPRC": round(float(average_precision_score(ys, ps)), 4),
                                      "AUROC": round(float(roc_auc_score(ys, ps)), 4),
                                      "prev": round(float(ys.mean()), 4)}
    return res

def run(epochs=12, bs=1024, lr=3e-4, wd=1e-2, d_model=512, n_cross=2, out_dir="."):
    manifest = json.load(open(f"{D}/v6_split_manifest_combined.json"))
    E, lens, idx = load_pertoken()
    trF = fix_activity_scalar(build_split("train", idx, manifest), E, lens)
    vaF = build_split("val", idx, manifest)   # activity scalar defaults; fix for consistency
    vaF = fix_activity_scalar(vaF, E, lens)

    ms = trF["mask_sign"]; ys = trF["y_sign"][ms]
    rep_frac = float((ys == SIGN_REP).mean()) if ms.any() else 0.5
    scw = torch.tensor([1.0, (1-rep_frac)/max(rep_frac, 1e-6), 1.0], dtype=torch.float32, device=DEV)
    n_tr = len(trF["e_rows"])
    log(f"train={n_tr} val={len(vaF['e_rows'])} bs={bs} batches/ep={n_tr//bs+1} rep_frac={rep_frac:.4f}")

    net = CATE(emb_dim=CODEC_DIM, d_model=d_model, n_scalar=3, nhead=8,
                    conv_kernel=int(os.environ.get("CATE_KERNEL", "5")),
                    conv_stride=int(os.environ.get("CATE_STRIDE", "2")),
                    extra_conv=(os.environ.get("CATE_EXTRA_CONV", "1") == "1"),
                    pool=int(os.environ.get("CATE_POOL", "2")),
                    n_self=1, n_cross=n_cross, trunk_dim=d_model).to(DEV)
    log(f"CATE params={sum(p.numel() for p in net.parameters())/1e6:.2f}M "
        f"(CODEC_DIM={CODEC_DIM}, d_model={d_model}, kernel={os.environ.get('CATE_KERNEL','5')}, stride={os.environ.get('CATE_STRIDE','2')})")
    # torch.compile: OFF by default. Its fork-based inductor worker pool (--kind=fork --workers=16)
    # DEADLOCKS when forking this 80GB-resident process — hung 40 min at 0% before first batch.
    # Eager bf16 on a 20M-param model trains fine; correctness/liveness > the compile speedup.
    # Set V6_COMPILE=1 to re-enable (only worth it once the fork-pool issue is resolved).
    if DEV == "cuda" and os.environ.get("V6_COMPILE", "0") == "1":
        net = torch.compile(net)
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
    rng = np.random.RandomState(0); best, best_state = -1, None
    n_batches = n_tr // bs + 1
    HB_EVERY = int(os.environ.get("CATE_HB_EVERY", "50"))
    for ep in range(epochs):
        net.train(); perm = rng.permutation(n_tr); t0 = time.time()
        for bi, s in enumerate(range(0, n_tr, bs)):
            idxb = perm[s:s+bs]; b = to_dev(gather(trF, E, lens, idxb)); opt.zero_grad()
            with torch.autocast("cuda", dtype=torch.bfloat16, enabled=(DEV=="cuda")):
                o = net(b["e_tok"], b["g_tok"], b["e_pad"], b["g_pad"], b["scalars"])
                loss, parts = multitask_loss(o, b, sign_class_weight=scw)
            loss.backward(); opt.step()
            if bi % HB_EVERY == 0:
                el = time.time() - t0; bps = (bi + 1) / max(el, 1e-9)
                json.dump({"phase": "train", "epoch": ep+1, "epochs": epochs, "batch": bi,
                           "n_batches": n_batches, "batches_per_s": round(bps, 3),
                           "epoch_eta_min": round((n_batches - bi) / max(bps, 1e-9) / 60, 1),
                           "run_loss": round(float(loss.item()), 4), "t": time.time(),
                           "gpu_mem_gb": round(torch.cuda.max_memory_allocated()/1e9, 1) if DEV=="cuda" else 0},
                          open(f"{out_dir}/progress.json", "w"))
                log(f"ep{ep+1} batch {bi}/{n_batches} {bps:.2f} b/s "
                    f"eta_ep {(n_batches-bi)/max(bps,1e-9)/60:.1f}min loss {loss.item():.4f}")
        res = evaluate(net, vaF, E, lens); a = res["link"]["AUPRC"]
        if a > best:
            best = a
            sd = (net._orig_mod if hasattr(net, "_orig_mod") else net).state_dict()
            best_state = {k: v.detach().cpu().clone() for k, v in sd.items()}
            # crash-resilience: persist the best checkpoint to disk + GCS the moment it improves,
            # so a box death mid-run loses at most the epochs since the last best (not everything).
            torch.save({"state": best_state, "epoch": ep+1, "best_val_link_AUPRC": best,
                        "d_model": d_model, "n_cross": n_cross}, f"{out_dir}/cate_model.pt")
            _put_gcs(f"{out_dir}/cate_model.pt", "cate_model.pt")
        json.dump({"epoch": ep+1, "epochs": epochs, "val_link_AUPRC": a, "best": best, "full": res,
                   "sec_ep": round(time.time()-t0, 1), "t": time.time(),
                   "gpu_mem_gb": round(torch.cuda.max_memory_allocated()/1e9, 1) if DEV=="cuda" else 0},
                  open(f"{out_dir}/progress.json", "w"))
        _put_gcs(f"{out_dir}/progress.json", "progress.json")
        log(f"ep{ep+1}/{epochs} val_link_AUPRC={a} full={res} ({time.time()-t0:.0f}s) "
            f"peak_gpu={torch.cuda.max_memory_allocated()/1e9:.1f}GB")

    core = (net._orig_mod if hasattr(net, "_orig_mod") else net)
    core.load_state_dict(best_state)
    iso = IsotonicRegression(out_of_bounds="clip", increasing=False)
    iso.fit(trF["scalars"][:, 0], trF["y_link"])
    torch.save({"state": best_state, "epoch": epochs, "best_val_link_AUPRC": best,
                "d_model": d_model, "n_cross": n_cross}, f"{out_dir}/cate_model.pt")
    _put_gcs(f"{out_dir}/cate_model.pt", "cate_model.pt")
    json.dump({"d_model": d_model, "n_cross": n_cross, "best_val_link_AUPRC": best,
               "codec_dim": CODEC_DIM,
               "params_M": round(sum(p.numel() for p in core.parameters())/1e6, 2)},
              open(f"{out_dir}/cate_meta.json", "w"))
    _put_gcs(f"{out_dir}/cate_meta.json", "cate_meta.json")

    teF = fix_activity_scalar(build_split("test", idx, manifest), E, lens)
    test_res = evaluate(core.to(DEV), teF, E, lens)
    json.dump(test_res, open(f"{out_dir}/cate_test_results.json", "w"), indent=1)
    _put_gcs(f"{out_dir}/cate_test_results.json", "cate_test_results.json")
    log("TEST", test_res)
    return test_res

if __name__ == "__main__":
    run(epochs=int(os.environ.get("V6_EPOCHS", "10")),
        bs=int(os.environ.get("V6_BS", "8192")),
        d_model=int(os.environ.get("V6_DMODEL", "512")),
        n_cross=int(os.environ.get("V6_NCROSS", "2")),
        out_dir=os.environ.get("V6_OUT", "."))

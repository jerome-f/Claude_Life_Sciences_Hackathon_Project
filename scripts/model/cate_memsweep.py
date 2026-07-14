#!/usr/bin/env python3
"""
CATE GPU-memory smoke: load the 128d codec RESIDENT on GPU, then run a real
forward+backward at a range of batch sizes, measuring torch.cuda.max_memory_allocated().
Picks the largest bs that leaves >HEADROOM_GB free. This is the make-or-break
pre-launch check for the GPU-resident design — do NOT assume bs=8192 fits when
the codec already occupies 52 GB of the 80 GB card.

Runs on ONE GPU (cuda:0). Uses a SUBSAMPLE of real train rows so the gather path,
pad masks, and loss are all exercised exactly as in the full run.
"""
import os, json, time
import numpy as np, torch
os.environ.setdefault("V6_CODEC_DIM", "128")
CODEC_DIM = int(os.environ["V6_CODEC_DIM"])
D = os.environ.get("E2G_DIR", ".")
DEV = "cuda"
HEADROOM_GB = float(os.environ.get("CATE_HEADROOM_GB", "5"))
import sys; sys.path.insert(0, D)
from cate_model import CATE
from v6_model import multitask_loss

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

torch.cuda.set_device(0)
# 1) codec resident on GPU (streamed, chunked)
E_np = np.load(f"{D}/E_combined_pertoken.npy", mmap_mode="r")
N = E_np.shape[0]; gb = E_np.nbytes/1e9
log(f"codec {E_np.shape} {gb:.1f} GB -> GPU")
E = torch.empty(E_np.shape, dtype=torch.float16, device=DEV)
CH = 100000
for s in range(0, N, CH):
    e = min(N, s+CH); E[s:e] = torch.from_numpy(np.ascontiguousarray(E_np[s:e]))
del E_np
codec_gb = torch.cuda.memory_allocated()/1e9
log(f"codec resident: {codec_gb:.1f} GB allocated")

# 2) model
d_model = int(os.environ.get("V6_DMODEL", "256"))
net = CATE(emb_dim=CODEC_DIM, d_model=d_model, n_scalar=3, nhead=8,
           conv_kernel=5, conv_stride=2, n_self=1, n_cross=2, trunk_dim=d_model).to(DEV)
opt = torch.optim.AdamW(net.parameters(), lr=3e-4)
scw = torch.tensor([1.0, 3.0, 1.0], device=DEV)
log(f"model params={sum(p.numel() for p in net.parameters())/1e6:.2f}M d_model={d_model}")

L = 168
def synth_batch(bs):
    # index real codec rows (exercises index_select gather); synth labels/masks/pad
    ridx = torch.randint(0, N, (bs,), device=DEV)
    e = E.index_select(0, ridx); g = E.index_select(0, torch.randint(0, N, (bs,), device=DEV))
    e_pad = torch.zeros(bs, L, dtype=torch.bool, device=DEV)
    g_pad = torch.zeros(bs, L, dtype=torch.bool, device=DEV)
    b = {"e_tok": e, "g_tok": g, "e_pad": e_pad, "g_pad": g_pad,
         "scalars": torch.randn(bs, 3, device=DEV),
         "y_link": torch.randint(0, 2, (bs,), device=DEV).float(),
         "y_sign": torch.randint(0, 3, (bs,), device=DEV),
         "mask_sign": torch.ones(bs, dtype=torch.bool, device=DEV),
         "y_elem": torch.zeros(bs, 4, device=DEV).scatter_(1, torch.randint(0, 4, (bs, 1), device=DEV), 1.0),
         "mask_elem": torch.ones(bs, dtype=torch.bool, device=DEV),
         "y_splice": torch.randint(0, 2, (bs,), device=DEV).float(),
         "mask_splice": torch.ones(bs, dtype=torch.bool, device=DEV)}
    return b

results = {}
best_bs = None
for bs in [512, 1024, 2048, 4096, 8192]:
    torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
    try:
        # 2 steps to capture optimizer state allocation
        for _ in range(2):
            b = synth_batch(bs); opt.zero_grad()
            with torch.autocast("cuda", dtype=torch.bfloat16):
                o = net(b["e_tok"], b["g_tok"], b["e_pad"], b["g_pad"], b["scalars"])
                loss, _ = multitask_loss(o, b, sign_class_weight=scw)
            loss.backward(); opt.step()
        peak = torch.cuda.max_memory_allocated()/1e9
        free = 80 - peak
        ok = free > HEADROOM_GB
        results[bs] = {"peak_gb": round(peak, 1), "free_gb": round(free, 1), "ok": ok}
        log(f"bs={bs}: peak={peak:.1f}GB free={free:.1f}GB {'OK' if ok else 'TOO TIGHT'}")
        if ok: best_bs = bs
    except torch.cuda.OutOfMemoryError:
        results[bs] = {"peak_gb": None, "oom": True}
        log(f"bs={bs}: OOM")
        break
    except RuntimeError as ex:
        if "out of memory" in str(ex).lower():
            results[bs] = {"peak_gb": None, "oom": True}; log(f"bs={bs}: OOM (RuntimeError)"); break
        raise

json.dump({"codec_gb": round(codec_gb, 1), "results": results, "recommended_bs": best_bs},
          open(f"{D}/memsweep.json", "w"))
log(f"RECOMMENDED bs={best_bs}")
print("MEMSWEEP DONE", json.dumps(results), flush=True)

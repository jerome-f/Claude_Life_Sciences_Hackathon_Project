#!/usr/bin/env python3
"""
v5 training driver — pooled two-tower multi-task signed-E2G (Option B).

Runs identically local (CPU smoke) and on tnr-0 A100. Reads combined-basis PCA
embeddings E_{human,mouse}_v4pca.npy + id_index_{human,mouse}.json and the v4
pair tables (same inputs as v4_train.py). Trains the structured multi-task head
end-to-end, early-stops on val link AUPRC, reports LINK + SIGN metrics head-to-head
against the v4 XGBoost champion, and saves the model + a calibrated distance prior
(compose-time, isotonic) for the CRISPRi correction.

Env: E2G_DIR (workdir with inputs). Device auto (cuda if available).
"""
import os, json, time, gc
import numpy as np, polars as pl, torch
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.isotonic import IsotonicRegression
from v5_model import TwoTowerE2G, multitask_loss, build_pair_features, SIGN_ACT, SIGN_REP

D = os.environ.get("E2G_DIR", ".")
DEV = "cuda" if torch.cuda.is_available() else "cpu"
t0 = time.time()
def log(*a): print(f"[{time.strftime('%H:%M:%S')} +{time.time()-t0:.0f}s]", *a, flush=True)


class PairDS(Dataset):
    """Holds resolved row-indices + labels; gathers embeddings per-item from a shared
    memmap E (kept on CPU; batches moved to device in the loop)."""
    def __init__(self, feat, E):
        self.f = feat; self.E = E
    def __len__(self): return len(self.f["e_rows"])
    def __getitem__(self, i):
        return {
            "e_emb": torch.from_numpy(np.asarray(self.E[self.f["e_rows"][i]], dtype=np.float32)),
            "g_emb": torch.from_numpy(np.asarray(self.E[self.f["g_rows"][i]], dtype=np.float32)),
            "scalars": torch.from_numpy(self.f["scalars"][i]),
            "y_link": torch.tensor(self.f["y_link"][i]),
            "y_sign": torch.tensor(self.f["y_sign"][i]),
            "mask_sign": torch.tensor(self.f["mask_sign"][i]),
            "y_elem": torch.from_numpy(self.f["y_elem"][i]),
            "mask_elem": torch.tensor(self.f["mask_elem"][i]),
            "y_splice": torch.tensor(self.f["y_splice"][i]),
            "mask_splice": torch.tensor(self.f["mask_splice"][i]),
        }

def collate(items):
    out = {}
    for k in items[0]:
        out[k] = torch.stack([it[k] for it in items])
    return out


def to_dev(batch, dev):
    return {k: v.to(dev, non_blocking=True) for k, v in batch.items()}


def evaluate(net, loader, dev):
    net.eval()
    yl, pl_, ys, ps = [], [], [], []
    with torch.no_grad():
        for b in loader:
            b = to_dev(b, dev)
            o = net(b["e_emb"], b["g_emb"], b["scalars"])
            yl.append(b["y_link"].cpu().numpy())
            pl_.append(torch.sigmoid(o["link"]).float().cpu().numpy())
            ms = b["mask_sign"]
            if ms.any():
                # P(repressive | linked) among {act,rep} — 2-class softmax on the two sign logits
                sl = o["sign"][ms][:, [SIGN_ACT, SIGN_REP]]
                p_rep = torch.softmax(sl, dim=-1)[:, 1]
                ys.append((b["y_sign"][ms] == SIGN_REP).float().cpu().numpy())
                ps.append(p_rep.float().cpu().numpy())
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


def load_emb(species):
    E = np.load(f"{D}/E_{species}_v4pca.npy", mmap_mode="r")
    idx = json.load(open(f"{D}/id_index_{species}.json"))
    return E, idx


def build_split(split, species_list=("human", "mouse")):
    feats = []
    Es = {}
    for sp in species_list:
        if sp not in Es:
            Es[sp], idx = load_emb(sp)
        p = pl.read_parquet(f"{D}/v4_pairs_{sp}_{split}.parquet")
        f = build_pair_features(p, Es[sp], idx)
        # offset row indices into a per-species stacked E later; here train each species
        # against its own E via a species tag
        f["species"] = sp
        feats.append(f)
        log(f"{sp} {split}: n={len(f['e_rows'])} miss={f['n_missing']} "
            f"link_prev={f['y_link'].mean():.3f} sign_conf={f['mask_sign'].sum()}")
    return feats, Es


def concat_feats_single_species(feats):
    """For a clean baseline we train per-species then pool metrics; but to train ONE joint
    model we stack embeddings into a single matrix with remapped indices."""
    # remap: build a combined E by stacking, offset g/e rows
    return feats


def train(epochs=15, bs=8192, lr=1e-3, wd=1e-4, emb_dim=None, out_dir=None):
    out_dir = out_dir or f"{D}/out_v5"; os.makedirs(out_dir, exist_ok=True)
    tr_feats, Es = build_split("train")
    va_feats, _ = build_split("val")

    # Build a single stacked embedding matrix across species so one model trains jointly.
    # Offset each species' rows into a combined index space.
    def stack(feats):
        mats, offset, per = [], 0, []
        # collect unique species order
        return feats
    # Simplest correct approach: concatenate per-species E and shift row indices.
    order = [f["species"] for f in tr_feats]
    E_parts, offsets, cur = [], {}, 0
    for sp in dict.fromkeys(order):
        offsets[sp] = cur; E_parts.append(np.asarray(Es[sp], dtype=np.float32)); cur += Es[sp].shape[0]
    E_all = np.concatenate(E_parts, axis=0)
    K = E_all.shape[1]
    log(f"combined E {E_all.shape}")

    def shift(feats):
        merged = {k: [] for k in ("e_rows","g_rows","scalars","y_link","y_sign","mask_sign",
                                  "y_elem","mask_elem","y_splice","mask_splice")}
        for f in feats:
            off = offsets[f["species"]]
            merged["e_rows"].append(f["e_rows"] + off)
            merged["g_rows"].append(f["g_rows"] + off)
            for k in ("scalars","y_link","y_sign","mask_sign","y_elem","mask_elem","y_splice","mask_splice"):
                merged[k].append(f[k])
        out = {k: (np.concatenate(v, axis=0) if v[0].ndim else np.concatenate(v)) for k, v in merged.items()}
        return out
    trF, vaF = shift(tr_feats), shift(va_feats)

    # sign class weight: upweight repressive (rarer)
    ms = trF["mask_sign"]; ys = trF["y_sign"][ms]
    rep_frac = float((ys == SIGN_REP).mean()) if ms.any() else 0.5
    scw = torch.tensor([1.0, (1 - rep_frac) / max(rep_frac, 1e-6), 1.0], dtype=torch.float32, device=DEV)
    log(f"sign confident train={int(ms.sum())} rep_frac={rep_frac:.4f}")

    tr = DataLoader(PairDS(trF, E_all), batch_size=bs, shuffle=True, collate_fn=collate,
                    num_workers=min(8, os.cpu_count() or 2), drop_last=False)
    va = DataLoader(PairDS(vaF, E_all), batch_size=32768, shuffle=False, collate_fn=collate,
                    num_workers=min(8, os.cpu_count() or 2))

    net = TwoTowerE2G(emb_dim=K, n_scalar=3).to(DEV)
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
    use_amp = (DEV == "cuda")
    best, best_state = -1, None
    for ep in range(epochs):
        net.train()
        for b in tr:
            b = to_dev(b, DEV); opt.zero_grad()
            with torch.autocast("cuda", dtype=torch.bfloat16, enabled=use_amp):
                o = net(b["e_emb"], b["g_emb"], b["scalars"])
                loss, parts = multitask_loss(o, b, sign_class_weight=scw)
            loss.backward(); opt.step()
        res = evaluate(net, va, DEV)
        auprc = res["link"]["AUPRC"]
        if auprc > best:
            best = auprc; best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
        log(f"ep{ep+1}/{epochs} val_link_AUPRC={auprc} full={res}")
    net.load_state_dict(best_state)

    # distance calibrated prior (compose-time CRISPRi correction, isotonic)
    iso = IsotonicRegression(out_of_bounds="clip", increasing=False)
    iso.fit(trF["scalars"][:, 0], trF["y_link"])

    torch.save(net.state_dict(), f"{out_dir}/v5_twotower.pt")
    json.dump({"K": int(K), "best_val_link_AUPRC": best}, open(f"{out_dir}/v5_meta.json", "w"))

    # test eval
    te_feats, _ = build_split("test")
    teF = shift(te_feats)
    te = DataLoader(PairDS(teF, E_all), batch_size=32768, shuffle=False, collate_fn=collate,
                    num_workers=min(8, os.cpu_count() or 2))
    test_res = evaluate(net, te, DEV)
    json.dump(test_res, open(f"{out_dir}/v5_test_results.json", "w"), indent=1)
    log("TEST", test_res)
    return net, test_res


if __name__ == "__main__":
    train()

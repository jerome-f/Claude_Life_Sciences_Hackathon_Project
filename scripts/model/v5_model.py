#!/usr/bin/env python3
"""
v5 signed-E2G — end-to-end differentiable POOLED two-tower multi-task model (Option B).

The honest realization of "gene on one side, enhancer on the other": a separate
enhancer tower and gene/TSS tower over the FROZEN combined-basis PCA-K pooled
embeddings, joined by an explicit interaction block, feeding a shared trunk that
fans out to a STRUCTURED multi-task head:

    link      : P(link)                       binary        (label: `link`)
    sign      : {activating, repressive, none} 3-way, link-conditioned
                                               (label: `consensus_sign` on confident positives)
    elem_class: {enhancer, promoter, CTCF, silencer} multi-label  (masked; annotation join TODO)
    splice    : P(near splice site)           binary, masked   (annotation join TODO)

Design commitments (settled with user, see V2G2T_PROJECT_CHECKPOINT / SESSION_CHECKPOINT):
  * Backbone stays FROZEN — embeddings are precomputed. The freeze-vs-train fork is
    a property of Option A (per-token cross-encoder), not this pooled baseline.
  * Structured head, NOT 8 flat sigmoids. Sign is a 3-way head conditioned on link,
    trained only on sign-confident linked positives. Element-class is MULTI-LABEL
    (enhancer/promoter/CTCF/silencer are not mutually exclusive).
  * Masked multi-task loss: each head contributes loss only on rows where its label
    exists, so aux heads (elem_class, splice) sit in the graph now and switch on when
    the annotation layers are joined — no rearchitecting.
  * Distance scalars stay as INPUT features here (in-domain DBNascent), matching v4.
    The CRISPRi-specific distance prior is a COMPOSE-TIME correction, not encoder input,
    and is applied at scoring, not here.

This module is backbone-agnostic about *where* it runs: build_pair_features() mirrors
v4_train.resolve_pairs/fill_into exactly (same id_index lookup, same feature order),
so the same code smoke-tests locally on synthetic tensors and trains on the A100
against the real GCS-staged E_{human,mouse}_v4pca.npy embeddings.
"""
from __future__ import annotations
import json, math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# ----- sign label encoding (3-way head) -----
SIGN_ACT, SIGN_REP, SIGN_NONE = 0, 1, 2   # class indices
ELEM_CLASSES = ["enhancer", "promoter", "CTCF", "silencer"]


# =====================================================================================
# Model
# =====================================================================================
class TwoTowerE2G(nn.Module):
    """Pooled two-tower multi-task signed-E2G head over frozen PCA embeddings.

    Args
    ----
    emb_dim   : K, dimension of the pooled PCA embedding (v4 combined basis = 613)
    n_scalar  : number of scalar pair features appended to the interaction vector
                (v4: dist_kb, contact_freq, abc = 3)
    tower_dim : width of each tower's output representation
    trunk_dim : width of the shared trunk
    tie_towers: share weights between enhancer and gene towers (default False —
                element vs promoter grammar differ)
    n_elem    : number of multi-label element classes
    """
    def __init__(self, emb_dim=613, n_scalar=3, tower_dim=256, trunk_dim=256,
                 dropout=0.3, tie_towers=False, n_elem=len(ELEM_CLASSES)):
        super().__init__()
        self.emb_dim, self.n_scalar, self.n_elem = emb_dim, n_scalar, n_elem

        def tower():
            return nn.Sequential(
                nn.Linear(emb_dim, tower_dim), nn.LayerNorm(tower_dim),
                nn.GELU(), nn.Dropout(dropout),
                nn.Linear(tower_dim, tower_dim), nn.LayerNorm(tower_dim), nn.GELU(),
            )
        self.enh_tower = tower()
        self.gene_tower = self.enh_tower if tie_towers else tower()

        # interaction: [t_e, t_g, |t_e - t_g|, t_e * t_g] + scalars
        inter_dim = 4 * tower_dim + n_scalar
        self.trunk = nn.Sequential(
            nn.Linear(inter_dim, trunk_dim), nn.LayerNorm(trunk_dim),
            nn.GELU(), nn.Dropout(dropout),
            nn.Linear(trunk_dim, trunk_dim), nn.LayerNorm(trunk_dim), nn.GELU(),
        )
        # structured heads
        self.head_link   = nn.Linear(trunk_dim, 1)          # P(link)
        self.head_sign   = nn.Linear(trunk_dim, 3)          # {act, rep, none}
        self.head_elem   = nn.Linear(trunk_dim, n_elem)     # multi-label
        self.head_splice = nn.Linear(trunk_dim, 1)          # masked binary

    def forward(self, e_emb, g_emb, scalars):
        t_e = self.enh_tower(e_emb)
        t_g = self.gene_tower(g_emb)
        inter = torch.cat([t_e, t_g, (t_e - t_g).abs(), t_e * t_g, scalars], dim=-1)
        z = self.trunk(inter)
        return {
            "link":   self.head_link(z).squeeze(-1),   # logit
            "sign":   self.head_sign(z),               # 3 logits
            "elem":   self.head_elem(z),               # n_elem logits
            "splice": self.head_splice(z).squeeze(-1), # logit
            "trunk":  z,
        }


# =====================================================================================
# Masked multi-task loss
# =====================================================================================
def multitask_loss(out, batch, weights=None, sign_class_weight=None):
    """Sum of per-head losses, each masked to rows where that head's label exists.

    batch keys (all tensors, same first dim):
      y_link   : {0,1} float, always present
      y_sign   : long class in {0=act,1=rep,2=none}; mask_sign selects confident linked positives
      mask_sign: bool
      y_elem   : float [N, n_elem] multi-label; mask_elem bool [N]
      mask_elem: bool
      y_splice : {0,1} float; mask_splice bool
      mask_splice: bool
    """
    w = {"link": 1.0, "sign": 1.0, "elem": 0.5, "splice": 0.5}
    if weights:
        w.update(weights)
    losses, parts = 0.0, {}

    # LINK — always on
    l_link = F.binary_cross_entropy_with_logits(out["link"], batch["y_link"])
    parts["link"] = float(l_link.detach()); losses = losses + w["link"] * l_link

    # SIGN — 3-way CE on confident linked positives only
    ms = batch["mask_sign"]
    if ms.any():
        l_sign = F.cross_entropy(out["sign"][ms], batch["y_sign"][ms],
                                 weight=sign_class_weight)
        parts["sign"] = float(l_sign.detach()); losses = losses + w["sign"] * l_sign

    # ELEM — multi-label BCE, masked
    me = batch["mask_elem"]
    if me.any():
        l_elem = F.binary_cross_entropy_with_logits(out["elem"][me], batch["y_elem"][me])
        parts["elem"] = float(l_elem.detach()); losses = losses + w["elem"] * l_elem

    # SPLICE — masked binary
    mp = batch["mask_splice"]
    if mp.any():
        l_sp = F.binary_cross_entropy_with_logits(out["splice"][mp], batch["y_splice"][mp])
        parts["splice"] = float(l_sp.detach()); losses = losses + w["splice"] * l_sp

    return losses, parts


# =====================================================================================
# Feature / label builder — mirrors v4_train.resolve_pairs + fill_into EXACTLY
# =====================================================================================
CONTACT_GAMMA, CONTACT_D0 = 1.0, 1000.0

def build_pair_features(pairs_pl, E, id_index):
    """Resolve a polars pair table into (e_rows, g_rows, scalars, labels/masks) using the
    frozen embedding matrix E (rows indexed by id_index: win_id -> row).

    Returns a dict of numpy arrays. Embedding gather is deferred to the training loop
    (index into E per-batch) so we never materialize the 4K-wide feature matrix v4 built —
    the two-tower model consumes E_e, E_g separately.
    """
    cols = pairs_pl.columns
    elem = pairs_pl["elem_id"].to_list()
    gene = pairs_pl["gene_id"].to_list()
    link = np.asarray(pairs_pl["link"].to_list(), dtype=np.float32)
    csign = pairs_pl["consensus_sign"].to_list() if "consensus_sign" in cols else [None]*len(elem)
    lroute = pairs_pl["label_route"].to_list() if "label_route" in cols else [None]*len(elem)
    adist = (pairs_pl["abs_distance"].to_list() if "abs_distance" in cols
             else pairs_pl["dist_bin_1kb"].to_list())

    ei, gi, yl, ysign, mask_sign, dd, keep = [], [], [], [], [], [], []
    miss = 0
    for i in range(len(elem)):
        je, jg = id_index.get(elem[i]), id_index.get(gene[i])
        if je is None or jg is None:
            miss += 1; continue
        keep.append(i)
        ei.append(je); gi.append(jg); yl.append(link[i])
        dd.append(adist[i] if adist[i] is not None else 0.0)
        # sign: confident linked positive with a defined ±1 sign
        s = csign[i]; conf = (lroute[i] == "confident") if lroute[i] is not None else False
        if link[i] == 1 and conf and s in (1, -1):
            ysign.append(SIGN_ACT if s == 1 else SIGN_REP); mask_sign.append(True)
        else:
            ysign.append(SIGN_NONE); mask_sign.append(False)

    ei = np.asarray(ei); gi = np.asarray(gi)
    dd = np.asarray(dd, dtype=np.float32)
    dist_kb = dd / 1000.0
    contact_freq = np.power(np.abs(dd) + CONTACT_D0, -CONTACT_GAMMA).astype(np.float32)
    contact_freq /= (contact_freq.max() + 1e-12)
    activity = np.linalg.norm(np.asarray(E[ei], dtype=np.float32), axis=1)
    abc = (contact_freq * activity).astype(np.float32); abc /= (abc.max() + 1e-12)
    scalars = np.stack([dist_kb, contact_freq, abc], axis=1).astype(np.float32)

    n = len(ei)
    return dict(
        e_rows=ei, g_rows=gi, scalars=scalars,
        y_link=np.asarray(yl, dtype=np.float32),
        y_sign=np.asarray(ysign, dtype=np.int64),
        mask_sign=np.asarray(mask_sign, dtype=bool),
        # aux heads: no labels yet -> all-masked-off placeholders (annotation join TODO)
        y_elem=np.zeros((n, len(ELEM_CLASSES)), dtype=np.float32),
        mask_elem=np.zeros(n, dtype=bool),
        y_splice=np.zeros(n, dtype=np.float32),
        mask_splice=np.zeros(n, dtype=bool),
        n_missing=miss,
        keep_idx=np.asarray(keep, dtype=np.int64),   # input-row indices that resolved (for aligning extra cols e.g. strand)
    )

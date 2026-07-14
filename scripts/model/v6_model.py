#!/usr/bin/env python3
"""
v6 signed-E2G — Option A: per-token TWO-TOWER CROSS-ENCODER over FROZEN NT-v2-500M
per-token embeddings.

Difference from v5 (Option B, pooled):
  * v5 consumed ONE pooled vector per window. v6 consumes the full per-token stream
    (~168 tokens x 1024) so a cross-attention layer can let enhancer tokens attend to
    gene tokens and back — the piece pooling threw away, and the reason A (not B) can
    later do V2E (a variant-token delta propagates through cross-attention instead of
    being averaged into a ~0.6% pooled-vector nudge).
  * PCA is GONE. A LEARNED linear projection (1024 -> d_model), TIED across the two
    towers, replaces it. PCA is variance-optimal (keeps GC/composition axes); the
    learned projection is task-optimal and cross-species-tying is the transfer mechanism
    that the shared v4 PCA basis used to provide.

Everything downstream of the pooled context vectors — interaction block, shared trunk,
the four structured masked heads, and multitask_loss — is REUSED verbatim from v5_model
so v5/v6 are a controlled comparison (only the encoder differs).

Design commitments unchanged: backbone frozen (embeddings precomputed & cached), masked
multi-task loss, distance stays a compose-time prior (NOT an encoder input). bf16 + compile.
"""
from __future__ import annotations
import math
import torch
import torch.nn as nn
import torch.nn.functional as F

# reuse the exact heads' loss + label encoding from v5 (controlled comparison)
from v5_model import (SIGN_ACT, SIGN_REP, SIGN_NONE, ELEM_CLASSES,
                      multitask_loss, build_pair_features)


class PositionalEncoding(nn.Module):
    """Standard sinusoidal PE, added to projected tokens (max ~256 tok >> our ~168)."""
    def __init__(self, d_model, max_len=512):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        pos = torch.arange(0, max_len).unsqueeze(1).float()
        div = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(pos * div)
        pe[:, 1::2] = torch.cos(pos * div)
        self.register_buffer("pe", pe.unsqueeze(0))  # (1, max_len, d_model)

    def forward(self, x):  # x: (B, L, d)
        return x + self.pe[:, : x.size(1)]


class CrossBlock(nn.Module):
    """One bidirectional cross-attention block: each stream attends to the OTHER,
    then a per-stream FFN. Pre-norm, residual. key_padding_mask handles variable length."""
    def __init__(self, d_model, nhead, ff_mult=4, dropout=0.1):
        super().__init__()
        self.ln_e1 = nn.LayerNorm(d_model); self.ln_g1 = nn.LayerNorm(d_model)
        self.attn_e = nn.MultiheadAttention(d_model, nhead, dropout=dropout, batch_first=True)
        self.attn_g = nn.MultiheadAttention(d_model, nhead, dropout=dropout, batch_first=True)
        self.ln_e2 = nn.LayerNorm(d_model); self.ln_g2 = nn.LayerNorm(d_model)
        def ff():
            return nn.Sequential(nn.Linear(d_model, ff_mult * d_model), nn.GELU(),
                                 nn.Dropout(dropout), nn.Linear(ff_mult * d_model, d_model))
        self.ff_e = ff(); self.ff_g = ff()

    def forward(self, e, g, e_pad, g_pad):
        # enhancer attends to gene; gene attends to enhancer (queries from self, kv from other)
        en, gn = self.ln_e1(e), self.ln_g1(g)
        e = e + self.attn_e(en, gn, gn, key_padding_mask=g_pad, need_weights=False)[0]
        g = g + self.attn_g(gn, en, en, key_padding_mask=e_pad, need_weights=False)[0]
        e = e + self.ff_e(self.ln_e2(e))
        g = g + self.ff_g(self.ln_g2(g))
        return e, g


class CrossEncoderE2G(nn.Module):
    """Per-token two-tower cross-encoder. Encoder differs from v5; heads identical.

    emb_dim  : NT per-token hidden width (1024) OR the storage-codec width if tokens are
               cached reduced (then set emb_dim to that width).
    d_model  : learned projection / attention width (user choice: 512).
    """
    def __init__(self, emb_dim=1024, d_model=512, n_scalar=3, nhead=8,
                 n_self=1, n_cross=2, trunk_dim=512, dropout=0.1,
                 n_elem=len(ELEM_CLASSES)):
        super().__init__()
        self.emb_dim, self.d_model = emb_dim, d_model
        # LEARNED projection, TIED across towers (replaces PCA; cross-species transfer)
        self.proj = nn.Linear(emb_dim, d_model)
        self.pe = PositionalEncoding(d_model)
        # per-tower self-attention (separate weights: enhancer vs promoter grammar differ)
        enc = lambda: nn.TransformerEncoderLayer(d_model, nhead, dim_feedforward=4 * d_model,
                                                 dropout=dropout, batch_first=True, norm_first=True)
        self.self_e = nn.TransformerEncoder(enc(), n_self)
        self.self_g = nn.TransformerEncoder(enc(), n_self)
        # bidirectional cross-attention stack
        self.cross = nn.ModuleList([CrossBlock(d_model, nhead, dropout=dropout) for _ in range(n_cross)])
        self.ln_out = nn.LayerNorm(d_model)

        # ---- downstream identical to v5 (interaction + trunk + heads) ----
        inter_dim = 4 * d_model + n_scalar
        self.trunk = nn.Sequential(
            nn.Linear(inter_dim, trunk_dim), nn.LayerNorm(trunk_dim), nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(trunk_dim, trunk_dim), nn.LayerNorm(trunk_dim), nn.GELU(),
        )
        self.head_link   = nn.Linear(trunk_dim, 1)
        self.head_sign   = nn.Linear(trunk_dim, 3)
        self.head_elem   = nn.Linear(trunk_dim, n_elem)
        self.head_splice = nn.Linear(trunk_dim, 1)

    @staticmethod
    def _masked_mean(x, pad):  # x:(B,L,d) pad:(B,L) True=pad
        keep = (~pad).unsqueeze(-1).float()
        return (x * keep).sum(1) / keep.sum(1).clamp(min=1.0)

    def forward(self, e_tok, g_tok, e_pad, g_pad, scalars):
        # e_tok/g_tok: (B, L, emb_dim) float; *_pad: (B, L) bool True=padding
        e = self.pe(self.proj(e_tok)); g = self.pe(self.proj(g_tok))
        e = self.self_e(e, src_key_padding_mask=e_pad)
        g = self.self_g(g, src_key_padding_mask=g_pad)
        for blk in self.cross:
            e, g = blk(e, g, e_pad, g_pad)
        t_e = self._masked_mean(self.ln_out(e), e_pad)
        t_g = self._masked_mean(self.ln_out(g), g_pad)
        inter = torch.cat([t_e, t_g, (t_e - t_g).abs(), t_e * t_g, scalars], dim=-1)
        z = self.trunk(inter)
        return {
            "link":   self.head_link(z).squeeze(-1),
            "sign":   self.head_sign(z),
            "elem":   self.head_elem(z),
            "splice": self.head_splice(z).squeeze(-1),
            "trunk":  z,
        }

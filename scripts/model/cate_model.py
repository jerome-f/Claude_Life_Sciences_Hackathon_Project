#!/usr/bin/env python3
"""
v6h signed-E2G — HYBRID CNN -> cross-attention over per-token PCA-64 embeddings.

Motivation (user choice): the pooled ceiling (~0.61 AUROC across v4/v5/v5d) is the pooling,
not the head. The only axis where CNN/attention are legitimate is the 168 REAL sequence
positions. So we go back to per-token (PCA-64 codec, 24 GB resident) and put a hybrid tower
on the token axis:

  proj(64 -> d_model) + positional encoding
  -> Conv1d stack over the 168 positions   (LOCAL MOTIF detection; the CNN inductive bias
     that PROcapNet/Puffin/Danko element models use), stride reduces length 168 -> ~L'
  -> per-tower self-attention              (within-window long-range)
  -> bidirectional cross-attention E<->G   (the E<->G interaction pooling threw away)
  -> masked-mean pool -> [t_e,t_g,|t_e-t_g|,t_e*t_g,scalars] interaction
  -> shared trunk -> v5 masked multitask heads (link/sign/elem/splice)

REUSES v6_model.CrossBlock + PositionalEncoding and v5_model heads/loss verbatim
(controlled comparison; only the CNN front-end is new).
"""
from __future__ import annotations
import os
import torch
import torch.nn as nn
from v6_model import PositionalEncoding, CrossBlock
from v5_model import SIGN_ACT, SIGN_REP, SIGN_NONE, ELEM_CLASSES, multitask_loss


class ConvStack(nn.Module):
    """Conv1d motif front-end over the token axis. Two conv blocks; the second strides by
    `stride` to downsample the 168 positions (cheaper attention, wider effective receptive
    field). Pre-norm-ish: Conv -> GELU -> LayerNorm (over channels). Returns (B, L', d) and
    the downsampled padding mask."""
    def __init__(self, d_model, kernel=5, stride=2, dropout=0.1, extra_conv=True, pool=2):
        super().__init__()
        pad = kernel // 2
        self.conv1 = nn.Conv1d(d_model, d_model, kernel, padding=pad)
        self.conv2 = nn.Conv1d(d_model, d_model, kernel, stride=stride, padding=pad)
        self.ln1 = nn.LayerNorm(d_model); self.ln2 = nn.LayerNorm(d_model)
        self.act = nn.GELU(); self.drop = nn.Dropout(dropout); self.stride = stride
        # optional 3rd conv + maxpool: deeper motif-syntax hierarchy + local translation
        # invariance (a motif shifted a few bp gives the same pooled response — the standard
        # Basset/DeepSEA/Basenji conv-pool stack). Also halves the token count into attention
        # (84 -> 42), which cuts the O(L^2) attention cost and activation memory.
        self.extra_conv = extra_conv; self.pool = pool
        if extra_conv:
            self.conv3 = nn.Conv1d(d_model, d_model, kernel, padding=pad)
            self.ln3 = nn.LayerNorm(d_model)
            self.maxpool = nn.MaxPool1d(pool)

    def forward(self, x, pad):
        # x: (B, L, d), pad: (B, L) True=padding. Conv1d wants (B, d, L).
        # zero padded positions so they don't leak into the conv
        x = x.masked_fill(pad.unsqueeze(-1), 0.0)
        h = x.transpose(1, 2)                          # (B, d, L)
        h = self.act(self.conv1(h))
        h = self.ln1(h.transpose(1, 2)).transpose(1, 2)
        h = self.drop(h)
        h = self.act(self.conv2(h))                    # (B, d, L') strided
        h = self.ln2(h.transpose(1, 2)).transpose(1, 2)  # (B, d, L')  keep channels-first for conv3
        # downsample the padding mask to match L' (a position is 'pad' if its stride window was all pad)
        pad_ds = pad[:, ::self.stride][:, : h.size(2)]
        if pad_ds.size(1) < h.size(2):                 # pad-right if conv kept one extra
            pad_ds = torch.cat([pad_ds, pad_ds[:, -1:].expand(-1, h.size(2) - pad_ds.size(1))], 1)
        if self.extra_conv:
            # zero padded positions before conv3 so they don't leak, then conv3 + maxpool
            h = h.masked_fill(pad_ds.unsqueeze(1), 0.0)
            h = self.act(self.conv3(h))                # (B, d, L')
            h = self.ln3(h.transpose(1, 2)).transpose(1, 2)
            h = self.maxpool(h)                        # (B, d, L'//pool)
            # downsample pad mask by pool: a pooled position is 'pad' only if ALL its window was pad
            Lp = h.size(2)
            pad_pool = pad_ds[:, : Lp * self.pool].reshape(pad_ds.size(0), Lp, self.pool).all(-1)
            pad_ds = pad_pool
        h = h.transpose(1, 2)                          # (B, L_out, d)
        return h, pad_ds


class CATE(nn.Module):
    def __init__(self, emb_dim=64, d_model=256, n_scalar=3, nhead=8,
                 conv_kernel=5, conv_stride=2, n_self=1, n_cross=2,
                 trunk_dim=256, dropout=0.1, n_elem=len(ELEM_CLASSES),
                 extra_conv=True, pool=2):
        super().__init__()
        self.emb_dim, self.d_model = emb_dim, d_model
        self.proj = nn.Linear(emb_dim, d_model)        # PCA-64 -> d_model, tied across towers
        self.pe = PositionalEncoding(d_model)
        # separate conv stacks (enhancer vs promoter grammar differ)
        self.conv_e = ConvStack(d_model, conv_kernel, conv_stride, dropout, extra_conv, pool)
        self.conv_g = ConvStack(d_model, conv_kernel, conv_stride, dropout, extra_conv, pool)
        enc = lambda: nn.TransformerEncoderLayer(d_model, nhead, dim_feedforward=4 * d_model,
                                                 dropout=dropout, batch_first=True, norm_first=True)
        self.self_e = nn.TransformerEncoder(enc(), n_self)
        self.self_g = nn.TransformerEncoder(enc(), n_self)
        self.cross = nn.ModuleList([CrossBlock(d_model, nhead, dropout=dropout) for _ in range(n_cross)])
        self.ln_out = nn.LayerNorm(d_model)
        # FiLM: condition each tower's pooled vector on the scalars (distance/contact/strand)
        # BEFORE the trunk. Distance is the dominant E-G prior; FiLM lets it gate the sequence
        # signal per-channel (scale+shift) so "trust this motif" depends on how far apart E and G
        # are. A shared MLP -> 2*d_model (gamma||beta), applied t = gamma*t + beta to t_e and t_g.
        self.use_film = os.environ.get("CATE_FILM", "1") == "1"
        if self.use_film:
            self.film = nn.Sequential(nn.Linear(n_scalar, d_model), nn.GELU(),
                                      nn.Linear(d_model, 2 * d_model))
            # init last layer to identity-ish (gamma~1, beta~0) so FiLM starts as a no-op
            nn.init.zeros_(self.film[-1].weight); nn.init.zeros_(self.film[-1].bias)
        inter_dim = 4 * d_model + n_scalar
        self.trunk = nn.Sequential(
            nn.Linear(inter_dim, trunk_dim), nn.LayerNorm(trunk_dim), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(trunk_dim, trunk_dim), nn.LayerNorm(trunk_dim), nn.GELU(),
        )
        self.head_link   = nn.Linear(trunk_dim, 1)
        self.head_sign   = nn.Linear(trunk_dim, 3)
        self.head_elem   = nn.Linear(trunk_dim, n_elem)
        self.head_splice = nn.Linear(trunk_dim, 1)

    @staticmethod
    def _masked_mean(x, pad):
        keep = (~pad).unsqueeze(-1).float()
        return (x * keep).sum(1) / keep.sum(1).clamp(min=1.0)

    def forward(self, e_tok, g_tok, e_pad, g_pad, scalars):
        e = self.pe(self.proj(e_tok)); g = self.pe(self.proj(g_tok))
        e, e_pad = self.conv_e(e, e_pad)               # CNN motif front-end + downsample
        g, g_pad = self.conv_g(g, g_pad)
        e = self.self_e(e, src_key_padding_mask=e_pad)
        g = self.self_g(g, src_key_padding_mask=g_pad)
        for blk in self.cross:
            e, g = blk(e, g, e_pad, g_pad)
        t_e = self._masked_mean(self.ln_out(e), e_pad)
        t_g = self._masked_mean(self.ln_out(g), g_pad)
        if self.use_film:
            gamma, beta = self.film(scalars).chunk(2, dim=-1)     # each (B, d_model)
            t_e = (1.0 + gamma) * t_e + beta                       # FiLM: scale+shift, identity-init
            t_g = (1.0 + gamma) * t_g + beta
        inter = torch.cat([t_e, t_g, (t_e - t_g).abs(), t_e * t_g, scalars], dim=-1)
        z = self.trunk(inter)
        return {
            "link":   self.head_link(z).squeeze(-1),
            "sign":   self.head_sign(z),
            "elem":   self.head_elem(z),
            "splice": self.head_splice(z).squeeze(-1),
            "trunk":  z,
        }

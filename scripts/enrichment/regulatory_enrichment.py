#!/usr/bin/env python
"""
Regulatory-class enrichment of fine-mapped GWAS credible-set variants.

Foreground = credible-set members at PIP>=0.5 (in sets with >=1 such member).
Background = PIP<0.5 members of the SAME sets (controls for locus/LD).
Annotation = ENCODE SCREEN cCRE registry v3 (hg38), multi-label overlap.
Fold-enrichment (FG frac / BG frac) per class + Fisher exact + Katz 95% CI.
"""
import json, urllib.request
import numpy as np
import polars as pl
from scipy import stats

COMBINED = "multibiobank_credible_sets_combined.parquet"   # gs://claude-hackathon/multibiobank/20260711/
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]

def fetch_ccre():
    """Pull the ENCODE4 cCRE registry (hg38) from the UCSC REST API, per chromosome."""
    rows = []
    for ch in CHROMS:
        u = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cCREregistry;chrom={ch}"
        req = urllib.request.Request(u, headers={"User-Agent": "python-urllib"})
        with urllib.request.urlopen(req, timeout=180) as r:
            d = json.loads(r.read())
        for f in d.get("cCREregistry", []):
            rows.append((f["chrom"].replace("chr", ""), f["chromStart"], f["chromEnd"],
                         f["name"], f["cCRE_class"]))
    return pl.DataFrame(rows, schema=["chrom", "start", "end", "name", "cCRE_class"], orient="row")

def overlap_labels(variants, ccre):
    """variant -> ';'-joined set of overlapping cCRE classes (multi-label containment)."""
    out = []
    for (ch,), sub in variants.group_by("chromosome"):
        cc = ccre.filter(pl.col("chrom") == ch)
        if cc.height == 0:
            out += [(v, "none") for v in sub["variant_hg38"].to_list()]; continue
        s = cc["start"].to_numpy(); e = cc["end"].to_numpy(); cl = cc["cCRE_class"].to_list()
        o = np.argsort(s); s = s[o]; e = e[o]; cl = [cl[i] for i in o]
        maxlen = int((e - s).max())
        p0 = sub["position"].to_numpy() - 1                     # 1-based -> 0-based
        idx = np.searchsorted(s, p0, side="right")
        for k, (pp, vid) in enumerate(zip(p0, sub["variant_hg38"].to_list())):
            hit = set(); j = idx[k] - 1
            while j >= 0 and s[j] > pp - maxlen - 1:
                if s[j] <= pp < e[j]: hit.add(cl[j])
                j -= 1
            out.append((vid, ";".join(sorted(hit)) if hit else "none"))
    return pl.DataFrame(out, schema=["variant_hg38", "labels"], orient="row")

def main():
    g = (pl.scan_parquet(COMBINED)
           .filter(pl.col("data_layer") == "GWAS")
           .select(["chromosome", "position", "PIP", "cs_id", "variant_hg38"]).collect())
    sparse = g.filter(pl.col("PIP") >= 0.5).select("cs_id").unique()
    gs = g.join(sparse, on="cs_id", how="inner")
    fg_v = gs.filter(pl.col("PIP") >= 0.5).select(["chromosome","position","variant_hg38"]).unique()
    bg_v = (gs.filter(pl.col("PIP") < 0.5).select(["chromosome","position","variant_hg38"]).unique()
              .join(fg_v.select("variant_hg38"), on="variant_hg38", how="anti"))

    ccre = fetch_ccre()
    fg = overlap_labels(fg_v, ccre); bg = overlap_labels(bg_v, ccre)
    n_fg, n_bg = fg.height, bg.height

    def expl(l): return (l.filter(pl.col("labels") != "none")
                          .with_columns(pl.col("labels").str.split(";").alias("cls")).explode("cls"))
    fe, be = expl(fg), expl(bg)
    rows = []
    for cls in ccre["cCRE_class"].unique().to_list():
        a = fe.filter(pl.col("cls") == cls).height; b = be.filter(pl.col("cls") == cls).height
        ff, bf = a / n_fg, b / n_bg
        fold = ff / bf if bf else np.nan
        _, p = stats.fisher_exact([[a, n_fg - a], [b, n_bg - b]])
        se = np.sqrt(1/a - 1/n_fg + 1/b - 1/n_bg) if a and b else np.nan
        rows.append((cls, a, b, round(100*ff,2), round(100*bf,2), round(fold,3),
                     round(np.exp(np.log(fold)-1.96*se),3), round(np.exp(np.log(fold)+1.96*se),3), p))
    enr = pl.DataFrame(rows, schema=["cCRE_class","fg_n","bg_n","fg_pct","bg_pct","fold","ci_lo","ci_hi","fisher_p"],
                       orient="row").sort("fold", descending=True)
    enr.write_csv("ccre_enrichment_by_class.csv")
    pl.concat([fg.with_columns(pl.lit("fg").alias("group")),
               bg.with_columns(pl.lit("bg").alias("group"))]).write_parquet("credible_set_variant_ccre_labels.parquet")
    print(enr)

if __name__ == "__main__":
    main()

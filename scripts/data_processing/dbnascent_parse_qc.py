
import duckdb, json, gzip
import pandas as pd
from pathlib import Path

SRC = "/data/dbnascent/DBNascent_data/bidirectional_gene_pairs/dbnascent_pairs.txt.gz"
OUT = Path("/data/out"); OUT.mkdir(exist_ok=True, parents=True)

con = duckdb.connect()
con.execute("PRAGMA threads=16")
# read as CSV (comma), header present
rel = f"read_csv_auto('{SRC}', delim=',', header=true, sample_size=200000)"

# --- global counts ---
n_rows = con.execute(f"SELECT count(*) FROM {rel}").fetchone()[0]
n_pairs = con.execute(f"SELECT count(DISTINCT pair_id) FROM {rel}").fetchone()[0]
n_genes = con.execute(f"SELECT count(DISTINCT gene_id) FROM {rel}").fetchone()[0]
n_bidir = con.execute(f"SELECT count(DISTINCT transcript_2) FROM {rel}").fetchone()[0]
n_tissue = con.execute(f"SELECT count(DISTINCT tissue) FROM {rel}").fetchone()[0]
print(f"rows={n_rows} pairs={n_pairs} genes={n_genes} bidir={n_bidir} tissues={n_tissue}")

# --- tissue breakdown ---
tissue_tbl = con.execute(f"SELECT tissue, count(*) n, count(DISTINCT pair_id) pairs FROM {rel} GROUP BY 1 ORDER BY n DESC").df()
tissue_tbl.to_csv(OUT/"qc_by_tissue.csv", index=False)
print("TISSUES:\n", tissue_tbl.to_string(index=False))

# --- pcc sign distribution overall + by significance ---
sign_tbl = con.execute(f"""
SELECT
  CASE WHEN pcc>0 THEN 'positive' WHEN pcc<0 THEN 'negative' ELSE 'zero' END sign,
  count(*) n,
  round(avg(pcc),4) mean_pcc,
  round(avg(nObs),1) mean_nObs
FROM {rel} GROUP BY 1 ORDER BY n DESC
""").df()
sign_tbl.to_csv(OUT/"qc_pcc_sign.csv", index=False)
print("SIGN:\n", sign_tbl.to_string(index=False))

# --- FDR stringency tiers ---
fdr_tbl = con.execute(f"""
SELECT
  CASE WHEN adj_p_BH < 0.001 THEN 'FDR<0.001'
       WHEN adj_p_BH < 0.01 THEN 'FDR<0.01'
       WHEN adj_p_BH < 0.05 THEN 'FDR<0.05'
       ELSE 'FDR>=0.05' END tier,
  count(*) n, count(DISTINCT pair_id) pairs
FROM {rel} GROUP BY 1 ORDER BY tier
""").df()
fdr_tbl.to_csv(OUT/"qc_fdr_tiers.csv", index=False)
print("FDR:\n", fdr_tbl.to_string(index=False))

# --- 3D support flags ---
flag_tbl = con.execute(f"""
SELECT gtex, chiapet, count(*) n FROM {rel} GROUP BY 1,2 ORDER BY n DESC
""").df()
flag_tbl.to_csv(OUT/"qc_3d_flags.csv", index=False)
print("3D FLAGS:\n", flag_tbl.to_string(index=False))

# --- distance / nObs / percent_transcribed distributions (quantiles) ---
dist_tbl = con.execute(f"""
SELECT
  quantile_cont(abs(distance_tss),[0.05,0.25,0.5,0.75,0.95]) abs_dist_tss_q,
  quantile_cont(nObs,[0.05,0.25,0.5,0.75,0.95]) nObs_q,
  quantile_cont(percent_transcribed_both,[0.05,0.25,0.5,0.75,0.95]) pct_trans_q,
  quantile_cont(pcc,[0.05,0.25,0.5,0.75,0.95]) pcc_q
FROM {rel}
""").df()
qd = {c: [float(x) for x in dist_tbl[c].iloc[0]] for c in dist_tbl.columns}
json.dump(qd, open(OUT/"qc_quantiles.json","w"), indent=2)
print("QUANTILES:", json.dumps(qd, indent=2))

# --- position / bidirectional_location breakdown ---
pos_tbl = con.execute(f"SELECT position, bidirectional_location, count(*) n FROM {rel} GROUP BY 1,2 ORDER BY n DESC").df()
pos_tbl.to_csv(OUT/"qc_position_location.csv", index=False)
print("POSITION/LOC:\n", pos_tbl.to_string(index=False))

# --- null/coverage check on key cols ---
nulls = con.execute(f"""
SELECT
  sum(CASE WHEN pcc IS NULL THEN 1 ELSE 0 END) pcc_null,
  sum(CASE WHEN adj_p_BH IS NULL THEN 1 ELSE 0 END) adjp_null,
  sum(CASE WHEN nObs IS NULL THEN 1 ELSE 0 END) nObs_null,
  sum(CASE WHEN distance_tss IS NULL THEN 1 ELSE 0 END) dist_null,
  sum(CASE WHEN gene_id IS NULL THEN 1 ELSE 0 END) gene_null
FROM {rel}
""").df()
nulls.to_csv(OUT/"qc_nulls.csv", index=False)
print("NULLS:\n", nulls.to_string(index=False))

summary = dict(n_rows=int(n_rows), n_pairs=int(n_pairs), n_genes=int(n_genes),
               n_bidir=int(n_bidir), n_tissues=int(n_tissue))
json.dump(summary, open(OUT/"qc_summary_dbnascent.json","w"), indent=2)
print("DONE")

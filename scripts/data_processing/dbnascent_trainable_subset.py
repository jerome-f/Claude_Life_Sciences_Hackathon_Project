
import duckdb, json, logging
from pathlib import Path
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log=logging.getLogger()

SRC="/data/dbnascent/DBNascent_data/bidirectional_gene_pairs/dbnascent_pairs.txt.gz"
OUT=Path("/data/out"); OUT.mkdir(exist_ok=True,parents=True)
con=duckdb.connect(); con.execute("PRAGMA threads=16"); con.execute("PRAGMA memory_limit='48GB'")
rel=f"read_csv_auto('{SRC}',delim=',',header=true,sample_size=200000)"

# ---- 1. clean per-tissue signed observations (the full significant table) ----
log.info("materializing cleaned per-tissue table")
con.execute(f"""
CREATE TABLE obs AS
SELECT
  pair_id, gene_id, transcript_1, transcript_2 AS bidirectional_id,
  transcript1_chrom AS chrom, transcript1_strand AS gene_strand,
  transcript1_start AS gene_start, transcript1_stop AS gene_stop,
  transcript2_start AS bidir_start, transcript2_stop AS bidir_stop,
  CASE WHEN transcript1_strand='+' THEN transcript1_start ELSE transcript1_stop END AS tss,
  (transcript2_start+transcript2_stop)/2.0 AS bidir_center,
  pcc, pval, adj_p_BH, nObs, t,
  distance_tss, distance_tes, "position", bidirectional_location,
  tissue, percent_transcribed_both, gtex, chiapet
FROM {rel}
""")
n_obs=con.execute("SELECT count(*) FROM obs").fetchone()[0]
log.info(f"obs rows={n_obs}")

# ---- 2. per-pair aggregate with cross-tissue consistency statistic ----
log.info("aggregating per-pair cross-tissue stats")
con.execute("""
CREATE TABLE pos AS
SELECT
  pair_id, gene_id, any_value(transcript_1) transcript_1,
  any_value(bidirectional_id) bidirectional_id, any_value(chrom) chrom,
  any_value(gene_strand) gene_strand, any_value(tss) tss,
  any_value(bidir_center) bidir_center,
  any_value(distance_tss) distance_tss, any_value("position") "position",
  any_value(bidirectional_location) bidirectional_location,
  count(DISTINCT tissue) n_tissues,
  sum(CASE WHEN pcc>0 THEN 1 ELSE 0 END) n_pos,
  sum(CASE WHEN pcc<0 THEN 1 ELSE 0 END) n_neg,
  avg(pcc) mean_pcc, median(pcc) median_pcc, min(pcc) min_pcc, max(pcc) max_pcc,
  stddev_samp(pcc) sd_pcc,
  min(adj_p_BH) best_adj_p, max(nObs) max_nObs, avg(percent_transcribed_both) mean_pct_trans,
  max(gtex) gtex_any, max(chiapet) chiapet_any,
  list(tissue) tissues, list(round(pcc,4)) pcc_by_tissue
FROM obs GROUP BY pair_id, gene_id
""")
# derived consistency fields
con.execute("""
CREATE TABLE positives AS
SELECT *,
  greatest(n_pos,n_neg)*1.0/n_tissues AS sign_concordance,
  CASE WHEN n_pos>0 AND n_neg>0 THEN 'sign_flip'
       WHEN n_pos>0 THEN 'consistent_positive'
       ELSE 'consistent_negative' END AS sign_class,
  CASE WHEN mean_pcc>0 THEN 1 ELSE -1 END AS consensus_sign,
  CASE WHEN best_adj_p<0.001 AND (gtex_any=1 OR chiapet_any=1) THEN true ELSE false END AS stringent_3d
FROM pos
""")
n_pos=con.execute("SELECT count(*) FROM positives").fetchone()[0]
log.info(f"positives (unique pairs)={n_pos}")

# QC on positives subset
sc=con.execute("SELECT sign_class,count(*) n FROM positives GROUP BY 1 ORDER BY n DESC").df()
st=con.execute("SELECT stringent_3d,count(*) n FROM positives GROUP BY 1").df()
print("SIGN_CLASS:\n",sc.to_string(index=False))
print("STRINGENT_3D:\n",st.to_string(index=False))
# concordance distribution among multi-tissue pairs
conc=con.execute("""SELECT
  quantile_cont(sign_concordance,[0.05,0.25,0.5,0.75,0.95]) conc_q
FROM positives WHERE n_tissues>1""").df()
print("CONCORDANCE_Q(multi-tissue):", [float(x) for x in conc['conc_q'].iloc[0]])

# ---- 3. distance-matched hard negatives ----
# universe of transcribed bidirectionals (appear in atlas), with coords
log.info("building bidirectional universe")
con.execute("""
CREATE TABLE bidirs AS
SELECT DISTINCT bidirectional_id, chrom, bidir_center FROM obs
""")
nb=con.execute("SELECT count(*) FROM bidirs").fetchone()[0]
log.info(f"transcribed bidirs={nb}")
# gene TSS universe (one representative TSS per gene_id)
con.execute("""
CREATE TABLE genes AS
SELECT gene_id, any_value(chrom) chrom, any_value(gene_strand) strand,
       min(tss) tss_min, max(tss) tss_max, any_value(tss) tss
FROM obs GROUP BY gene_id
""")
ng=con.execute("SELECT count(*) FROM genes").fetchone()[0]
log.info(f"genes={ng}")

# candidate within-window pairs (|tss - bidir_center| <= 1Mb), NOT significant partners
log.info("enumerating within-window candidates + anti-join significant")
con.execute("""
CREATE TABLE cand AS
SELECT g.gene_id, g.chrom, g.tss, b.bidirectional_id, b.bidir_center,
       abs(g.tss - b.bidir_center) AS dist
FROM genes g JOIN bidirs b
  ON g.chrom=b.chrom AND abs(g.tss-b.bidir_center)<=1000000
""")
nc=con.execute("SELECT count(*) FROM cand").fetchone()[0]
log.info(f"within-window candidate pairs={nc}")
con.execute("""
CREATE TABLE negatives_all AS
SELECT c.gene_id, c.chrom, c.tss, c.bidirectional_id, c.bidir_center, c.dist,
       c.gene_id || '~' || c.bidirectional_id AS neg_key
FROM cand c
LEFT JOIN (SELECT DISTINCT gene_id, bidirectional_id FROM obs) p
  ON c.gene_id=p.gene_id AND c.bidirectional_id=p.bidirectional_id
WHERE p.gene_id IS NULL
""")
nn=con.execute("SELECT count(*) FROM negatives_all").fetchone()[0]
log.info(f"implied negatives (within-window, not significant)={nn}")

# distance-match: sample negatives to match positive distance distribution, ratio ~3:1
# bucket by 50kb distance bins; sample negatives proportional to positive counts * ratio
RATIO=3
con.execute(f"""
CREATE TABLE pos_bins AS
SELECT CAST(abs(distance_tss)/50000 AS INT) bin, count(*) n_pos FROM positives GROUP BY 1
""")
con.execute(f"""
CREATE TABLE neg_binned AS
SELECT *, CAST(dist/50000 AS INT) bin,
       row_number() OVER (PARTITION BY CAST(dist/50000 AS INT) ORDER BY random()) rn
FROM negatives_all
""")
con.execute(f"""
CREATE TABLE negatives AS
SELECT n.gene_id, n.chrom, n.tss, n.bidirectional_id, n.bidir_center, n.dist, n.bin
FROM neg_binned n JOIN pos_bins p ON n.bin=p.bin
WHERE n.rn <= p.n_pos*{RATIO}
""")
nneg=con.execute("SELECT count(*) FROM negatives").fetchone()[0]
log.info(f"distance-matched negatives (ratio {RATIO}:1 target)={nneg}")

# ---- 4. write parquet outputs ----
log.info("writing parquet")
con.execute(f"COPY (SELECT * FROM positives) TO '{OUT}/dbnascent_e2g_positives.parquet' (FORMAT parquet, COMPRESSION zstd)")
con.execute(f"COPY (SELECT * FROM obs) TO '{OUT}/dbnascent_e2g_pertissue_signed.parquet' (FORMAT parquet, COMPRESSION zstd)")
con.execute(f"COPY (SELECT *, false AS is_positive FROM negatives) TO '{OUT}/dbnascent_e2g_negatives.parquet' (FORMAT parquet, COMPRESSION zstd)")

summary=dict(
  n_pertissue_obs=int(n_obs), n_positives_pairs=int(n_pos),
  n_transcribed_bidirs=int(nb), n_genes=int(ng),
  n_candidate_within_window=int(nc), n_implied_negatives=int(nn),
  n_matched_negatives=int(nneg), neg_ratio_target=RATIO,
  sign_class=sc.set_index('sign_class')['n'].to_dict(),
  stringent_3d=st.set_index('stringent_3d')['n'].to_dict(),
)
json.dump(summary, open(OUT/"trainable_subset_summary.json","w"), indent=2, default=int)
print("SUMMARY:", json.dumps(summary, indent=2, default=int))
log.info("DONE")

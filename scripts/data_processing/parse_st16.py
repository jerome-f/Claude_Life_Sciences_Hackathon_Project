
import openpyxl, json, re
import pandas as pd
from pathlib import Path
OUT=Path("/data/ukbppp"); 
wb=openpyxl.load_workbook("/data/ukbppp/MOESM3.xlsx", read_only=True)
s=wb["ST16"]
rows=list(s.iter_rows(values_only=True))
# header is row index 4 (0-based); data starts row 5
hdr=[str(x) if x is not None else f"col{i}" for i,x in enumerate(rows[4])]
# first column header is on row 3 for protein id
hdr[0]="UKBPPP_ProteinID"
data=[r for r in rows[5:] if any(v is not None for v in r)]
df=pd.DataFrame(data, columns=hdr[:len(data[0])])
print("rows", len(df))
print("cols", list(df.columns))
# clean names
df=df.rename(columns={
  'Variant ID':'top_variantId_hg19','rsID':'top_rsid','PIP':'top_pip',
  'log10(BF)':'log10BF','Cis/trans':'cis_trans','Test region hg19':'test_region_hg19',
  'CS size':'cs_size','CS variant IDs':'cs_variant_ids_hg19','CS rsIDs':'cs_rsids'})
# protein id -> gene symbol (format: GENE:UniProt:OID:vX)
df['gene_symbol']=df['UKBPPP_ProteinID'].astype(str).str.split(':').str[0]
print(df[['UKBPPP_ProteinID','gene_symbol','top_variantId_hg19','top_pip','cis_trans','cs_size']].head(6).to_string(index=False))
# counts
print("n_signals", len(df))
print("n_proteins", df['gene_symbol'].nunique())
print("cis/trans:\n", df['cis_trans'].value_counts().to_string())
print("cs_size present:", df['cs_size'].notna().sum())
df.to_parquet(OUT/"ukbppp_st16_insample_pqtl_hg19.parquet", index=False)
json.dump(dict(n_signals=int(len(df)), n_proteins=int(df['gene_symbol'].nunique()),
               cis_trans=df['cis_trans'].value_counts().to_dict(),
               build="hg19", note="UKB-PPP Sun2023 Nature ST16 in-sample SuSiE credible sets; coords hg19, need liftover to hg38"),
          open(OUT/"ukbppp_st16_summary.json","w"), indent=2, default=int)
print("DONE")

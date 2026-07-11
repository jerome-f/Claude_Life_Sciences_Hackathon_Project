import pandas as pd, json, re, logging
from pathlib import Path
from pyliftover import LiftOver
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log=logging.getLogger()
OUT=Path("/data/ukbppp")
lo=LiftOver('/data/ukbppp/hg19ToHg38.over.chain.gz')

df=pd.read_parquet(OUT/"ukbppp_st16_insample_pqtl_hg19.parquet")
log.info(f"input signals={len(df)}")

def parse_vid(v):
    # format chr:pos:ref:alt:imp:vN  (chr may be numeric or with 'chr'); some may be missing
    if v is None or str(v)=='nan': return None
    p=str(v).split(':')
    if len(p)<4: return None
    c=p[0]; pos=p[1]; ref=p[2]; alt=p[3]
    return c,int(pos),ref,alt

def lift_one(c,pos):
    chrom = c if str(c).startswith('chr') else 'chr'+str(c)
    r=lo.convert_coordinate(chrom, pos-1)  # pyliftover is 0-based
    if not r: return None
    nchrom, npos0, strand, _ = r[0]
    return nchrom.replace('chr',''), npos0+1  # back to 1-based, strip chr

# lift top variant
def lift_vid(v):
    pv=parse_vid(v)
    if not pv: return None, 'unparseable'
    c,pos,ref,alt=pv
    l=lift_one(c,pos)
    if not l: return None, 'unmapped'
    nc,npos=l
    return f"{nc}_{npos}_{ref}_{alt}", 'ok'

res=df['top_variantId_hg19'].map(lift_vid)
df['top_variantId_hg38']=[r[0] for r in res]
df['top_lift_status']=[r[1] for r in res]

# lift each CS variant list (comma or semicolon separated)
def lift_cs(cell):
    if cell is None or str(cell)=='nan': return None, 0, 0
    ids=re.split(r'[;,]\s*', str(cell))
    ids=[i for i in ids if i.strip()]
    out=[]; nmap=0
    for v in ids:
        vid,st=lift_vid(v)
        if vid: out.append(vid); nmap+=1
    return ('|'.join(out) if out else None), len(ids), nmap
cs=df['cs_variant_ids_hg19'].map(lift_cs)
df['cs_variant_ids_hg38']=[c[0] for c in cs]
df['cs_n_input']=[c[1] for c in cs]
df['cs_n_mapped']=[c[2] for c in cs]

# QC
top_ok=(df['top_lift_status']=='ok').sum()
top_unmapped=(df['top_lift_status']=='unmapped').sum()
top_unparse=(df['top_lift_status']=='unparseable').sum()
cs_map_rate = df['cs_n_mapped'].sum()/max(df['cs_n_input'].sum(),1)
log.info(f"top variant: ok={top_ok} unmapped={top_unmapped} unparseable={top_unparse}")
log.info(f"CS variants: {df['cs_n_mapped'].sum()}/{df['cs_n_input'].sum()} mapped ({100*cs_map_rate:.2f}%)")

df.to_parquet(OUT/"ukbppp_st16_insample_pqtl_hg38.parquet", index=False)
summary=dict(n_signals=int(len(df)), top_ok=int(top_ok), top_unmapped=int(top_unmapped),
   top_unparseable=int(top_unparse), top_map_rate=float(top_ok/len(df)),
   cs_variants_input=int(df['cs_n_input'].sum()), cs_variants_mapped=int(df['cs_n_mapped'].sum()),
   cs_map_rate=float(cs_map_rate), build_from="hg19", build_to="hg38",
   variantId_format="chr_pos_ref_alt (OT-compatible, no chr prefix)")
json.dump(summary, open(OUT/"ukbppp_st16_hg38_liftover_qc.json","w"), indent=2)
print("SUMMARY:", json.dumps(summary, indent=2))
print(df[['gene_symbol','top_variantId_hg19','top_variantId_hg38','top_lift_status','cs_n_input','cs_n_mapped']].head(6).to_string(index=False))

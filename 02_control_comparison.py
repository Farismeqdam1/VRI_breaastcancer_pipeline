

import os, sys, gzip, subprocess
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


WORK_DIR    = Path(".")
CONTROL_DIR = WORK_DIR / "healthy_control"
CONTROL_DIR.mkdir(exist_ok=True)

MASTER_XLSX = WORK_DIR / "output" / "methylation" / "Methylation_Analysis_Master.xlsx"

GSE_GZ  = CONTROL_DIR / "GSE40279_series_matrix.txt.gz"
GPL_FILE = CONTROL_DIR / "GPL13534.annot.txt"
PROBE_STATS_CACHE = CONTROL_DIR / "GSE40279_probe_stats.tsv"
PROBE_GENE_CACHE  = CONTROL_DIR / "probe_gene_mapping.tsv"

OUTPUT_DIR = WORK_DIR / "output" / "methylation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GSE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz"

DELTA_THRESH = 0.10
DELTA_STRONG = 0.20
FDR_THRESH   = 0.05
BC_N         = 90     # Number of BC samples with methylation data
HEALTHY_N    = 656    # GSE40279 healthy samples

# Tissue-specific gene filter (removes false positives)
TISSUE_GENES = {
    "neural": ["GFAP","MBP","OLIG1","OLIG2","NEUROD1","NEUROD2","SOX1","SOX10","PAX6"],
    "liver": ["ALB","AFP","APOA1","APOB","APOC3","SERPINA1","FGA","FGB","FGG"],
    "kidney": ["AQP1","AQP2","SLC12A1","UMOD","NPHS1","NPHS2"],
    "muscle": ["MYH6","MYH7","TNNI3","TNNT2","ACTA1","DES"],
    "pancreas": ["INS","GCG","SST","PPY","PDX1"],
    "reproductive": ["DAZL","DDX4","SYCP1","SYCP3","SPO11"],
    "eye": ["RHO","OPN1SW","CRX","NRL"],
    "skin": ["KRT1","KRT5","KRT10","KRT14","LOR","IVL","FLG"],
    "lung": ["SFTPA1","SFTPB","SFTPC","SFTPD","SCGB1A1"],
}

CANCER_PROTECT = {"BRCA1","RASSF1","APC","CDH1","GSTP1","PTEN","TP53","DAPK1",
                  "CDKN2A","MLH1","MSH2","MSH6","MGMT"}

ALL_TISSUE = set()
GENE2TISSUE = {}
for t, gl in TISSUE_GENES.items():
    for g in gl:
        ALL_TISSUE.add(g)
        GENE2TISSUE.setdefault(g, []).append(t)
ALL_TISSUE -= CANCER_PROTECT


def hdr(s):
    print(f"\n{'='*70}\n  {s}\n{'='*70}")


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction."""
    pv = np.array(pvals, dtype=float)
    v = ~np.isnan(pv)
    n = v.sum()
    if n == 0: return pv.copy()
    r = np.full_like(pv, np.nan)
    idx_valid = np.where(v)[0]
    pv_valid = pv[idx_valid]
    order = np.argsort(pv_valid)
    adj = pv_valid[order] * n / (np.arange(1, n+1))
    adj = np.clip(adj, 0, 1)
    for i in range(len(adj)-2, -1, -1):
        if adj[i] > adj[i+1]: adj[i] = adj[i+1]
    out = np.empty(n)
    out[order] = adj
    r[idx_valid] = out
    return r


def download(url, dest, label):
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"  [cached] {dest.name}")
        return True
    print(f"  Downloading {label}...")
    try:
        subprocess.run(["wget", "-q", "-O", str(dest), url], check=True, timeout=600)
        return True
    except Exception:
        try:
            subprocess.run(["curl", "-L", "-o", str(dest), url], check=True, timeout=600)
            return True
        except Exception as e:
            print(f"  FAILED: {e}")
            return False


def parse_gpl(fp):
    """Parse GPL13534 annotation to get probe→gene mapping."""
    print(f"  Parsing {fp.name}...")
    opener = gzip.open if str(fp).endswith('.gz') else open
    with opener(fp, 'rt') as f:
        hl = None
        for i, line in enumerate(f):
            if line.startswith('#') or line.startswith('!') or line.startswith('^'):
                continue
            if 'ID' in line and '\t' in line:
                hl = i; break
    df = pd.read_csv(fp, sep='\t', skiprows=hl, low_memory=False)

    ic = gc = rc = None
    for c in df.columns:
        cl = c.lower().strip()
        if cl == 'id': ic = c
        if 'refgene_name' in cl: gc = c
        if 'refgene_group' in cl: rc = c
    if ic is None: ic = df.columns[0]

    recs = []
    for _, row in df.iterrows():
        gs = str(row.get(gc, '')).strip()
        rs = str(row.get(rc, '')).strip()
        if gs in ('', 'nan', 'NA', '.'): continue
        for j, g in enumerate(gs.split(';')):
            rl = rs.split(';')
            r = rl[j].strip() if j < len(rl) else 'Unknown'
            recs.append({'probe_id': str(row[ic]).strip(), 'gene': g.strip(), 'region': r})
    result = pd.DataFrame(recs)
    print(f"    {len(result):,} probe-gene pairs, {result['gene'].nunique():,} genes")
    return result


def main():
    hdr("CONTROL COMPARISON PIPELINE")

    # Load BC data
    hdr("STEP 1: Load BC cohort")
    gs = pd.read_excel(MASTER_XLSX, sheet_name='Gene_Level_Summary')
    print(f"  {len(gs):,} genes loaded")

    # Download healthy reference
    hdr("STEP 2: Download GSE40279")
    if not download(GSE_URL, GSE_GZ, "GSE40279"):
        sys.exit(1)
    if not GPL_FILE.exists():
        print(f"  ERROR: Place GPL13534.annot.txt in {CONTROL_DIR}/")
        sys.exit(1)

    # Probe→gene mapping
    hdr("STEP 3: Probe-gene mapping")
    if PROBE_GENE_CACHE.exists():
        pg = pd.read_csv(PROBE_GENE_CACHE, sep='\t')
    else:
        pg = parse_gpl(GPL_FILE)
        pg.to_csv(PROBE_GENE_CACHE, sep='\t', index=False)

    pp = pg[pg['region'].isin(['TSS200','TSS1500'])]
    bp = pg[pg['region'] == 'Body']

    # Probe stats
    hdr("STEP 4: Healthy probe statistics")
    if PROBE_STATS_CACHE.exists():
        ps = pd.read_csv(PROBE_STATS_CACHE, sep='\t', index_col=0)
    else:
        print("  Parsing series matrix (this takes a few minutes)...")
        opener = gzip.open if str(GSE_GZ).endswith('.gz') else open
        with opener(GSE_GZ, 'rt') as f:
            for line in f:
                if line.startswith('"ID_REF"'): break
        bm = pd.read_csv(GSE_GZ, sep='\t', skiprows=lambda x: x < 0,
                          index_col=0, na_values=['null','NA'], low_memory=False)
        # Re-read properly
        import io
        lines = []
        with opener(GSE_GZ, 'rt') as f:
            found = False
            for line in f:
                if line.startswith('"ID_REF"'):
                    found = True
                if found:
                    if line.startswith('!series_matrix_table_end'):
                        break
                    lines.append(line)
        bm = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', index_col=0,
                          na_values=['null','NA'], low_memory=False)
        ps = pd.DataFrame({'healthy_mean': bm.mean(1), 'healthy_sd': bm.std(1)})
        ps.to_csv(PROBE_STATS_CACHE, sep='\t')
    print(f"  {len(ps):,} probes")

    # Gene-level healthy reference
    hdr("STEP 5: Gene-level healthy reference")
    pps = pp.merge(ps[['healthy_mean','healthy_sd']], left_on='probe_id', right_index=True, how='inner')
    hp = pps.groupby('gene').agg(
        h_prom_mean=('healthy_mean','mean'),
        h_prom_sd=('healthy_sd','mean'),
        n_probes=('probe_id','nunique')).reset_index()
    print(f"  Promoter reference: {len(hp):,} genes")

    # Merge + Δβ
    hdr("STEP 6: Compute Δβ")
    c = gs.merge(hp, on='gene', how='left')
    c['delta_beta'] = c['mean_beta'] - c['h_prom_mean']
    nw = c['h_prom_mean'].notna().sum()
    print(f"  Coverage: {nw:,}/{len(c):,} ({nw/len(c)*100:.1f}%)")

    # Tissue filter
    c['tissue'] = c['gene'].apply(lambda g: '; '.join(GENE2TISSUE[g]) if g in ALL_TISSUE else '')
    c['is_ts'] = c['tissue'] != ''
    c['is_cp'] = c['gene'].isin(CANCER_PROTECT)

    # Welch's t-test
    hdr("STEP 7: Welch's t-test + BH-FDR")
    ta = np.full(len(c), np.nan)
    pa = np.full(len(c), np.nan)
    for i, r in c.iterrows():
        m1, s1 = r.get('mean_beta', np.nan), r.get('sd_beta', np.nan)
        m2, s2 = r.get('h_prom_mean', np.nan), r.get('h_prom_sd', np.nan)
        if any(pd.isna(x) for x in [m1, s1, m2, s2]): continue
        se = np.sqrt(s1**2/BC_N + s2**2/HEALTHY_N)
        if se == 0: continue
        t = (m1 - m2) / se
        nu = (s1**2/BC_N + s2**2/HEALTHY_N)**2
        de = (s1**2/BC_N)**2/(BC_N-1) + (s2**2/HEALTHY_N)**2/(HEALTHY_N-1)
        if de == 0: continue
        ta[i] = t
        pa[i] = 2 * stats.t.sf(abs(t), nu/de)
    c['t_stat'] = ta
    c['p_value'] = pa
    c['p_fdr'] = bh_fdr(pa)

    # Classify
    hdr("STEP 8: Classify")
    def classify(r):
        has_ctrl = pd.notna(r.get('h_prom_mean'))
        ts, cp = r.get('is_ts', False), r.get('is_cp', False)
        db = r.get('delta_beta', np.nan)
        fdr = r.get('p_fdr', np.nan)
        b = r.get('mean_beta', 0)
        if cp:
            return f"KNOWN_CANCER_GENE (Δβ={db:+.3f})" if has_ctrl and pd.notna(db) else "KNOWN_CANCER_GENE"
        if has_ctrl and pd.notna(db):
            sig = pd.notna(fdr) and fdr < FDR_THRESH
            if abs(db) < 0.05: return "NORMAL_BLOOD"
            if db > DELTA_STRONG and sig: return "TISSUE_SPECIFIC" if ts else "CANCER_SPECIFIC_STRONG"
            if db > DELTA_THRESH and sig: return "TISSUE_SPECIFIC" if ts else "CANCER_SPECIFIC"
            if db > DELTA_THRESH: return "ELEVATED_NS"
            if db < -DELTA_THRESH: return "HYPOMETHYLATED_VS_HEALTHY"
            return "MARGINAL"
        if b >= 0.5: return "POSSIBLY_CANCER (no ctrl)"
        if b >= 0.2: return "ELEVATED_NO_CTRL"
        return "LOW_METHYLATION"

    c['classification'] = c.apply(classify, axis=1)
    print(c['classification'].value_counts().to_string())

    # Save
    hdr("STEP 9: Save")
    oc = [x for x in ['gene','is_bc_gene','primary_promoter','promoter_class','chr',
                       'mean_beta','sd_beta','cv','h_prom_mean','h_prom_sd','n_probes',
                       'delta_beta','t_stat','p_value','p_fdr','is_ts','tissue','is_cp',
                       'classification','interpretation'] if x in c.columns]

    cancer = c[c['classification'].str.contains('CANCER_SPECIFIC|KNOWN_CANCER|POSSIBLY_CANCER', na=False)]
    normal = c[c['classification'].str.contains('NORMAL_BLOOD|TISSUE_SPECIFIC', na=False)]

    out1 = OUTPUT_DIR / "Control_Comparison_Results.xlsx"
    with pd.ExcelWriter(out1, engine='openpyxl') as w:
        cancer.sort_values('delta_beta', ascending=False).to_excel(w, 'Cancer_Specific_Hits', index=False)
        normal.to_excel(w, 'Normal_Blood_Filtered', index=False)
        c[oc].sort_values('delta_beta', ascending=False).to_excel(w, 'Full_Comparison', index=False)
    print(f"  Saved: {out1}")

    hdr("DONE")


if __name__ == '__main__':
    main()

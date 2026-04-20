
import os, sys, gzip, urllib.request
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


WORK_DIR    = Path(".")
OUTPUT_DIR  = WORK_DIR / "output" / "methylation"
TCGA_DIR    = WORK_DIR / "tcga_cache"
TCGA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input files
PRIORITY_GENES = WORK_DIR / "output" / "methylation" / "priority_genes.txt"  # One gene per line
PROBE_CACHE    = WORK_DIR / "healthy_control" / "probe_gene_mapping.tsv"
BLOOD_RESULTS  = WORK_DIR / "output" / "methylation" / "Control_Comparison_Results.xlsx"

OUTPUT = OUTPUT_DIR / "TCGA_Validation_Results.xlsx"

# TCGA datasets (UCSC Xena)
XENA_BASE = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download"
TCGA_DATASETS = {
    'BRCA': f"{XENA_BASE}/TCGA.BRCA.sampleMap%2FHumanMethylation450.gz",
    'LUAD': f"{XENA_BASE}/TCGA.LUAD.sampleMap%2FHumanMethylation450.gz",
    'COAD': f"{XENA_BASE}/TCGA.COAD.sampleMap%2FHumanMethylation450.gz",
    'KIRC': f"{XENA_BASE}/TCGA.KIRC.sampleMap%2FHumanMethylation450.gz",
    'LIHC': f"{XENA_BASE}/TCGA.LIHC.sampleMap%2FHumanMethylation450.gz",
    'LUSC': f"{XENA_BASE}/TCGA.LUSC.sampleMap%2FHumanMethylation450.gz",
    'PRAD': f"{XENA_BASE}/TCGA.PRAD.sampleMap%2FHumanMethylation450.gz",
    'UCEC': f"{XENA_BASE}/TCGA.UCEC.sampleMap%2FHumanMethylation450.gz",
}

PAN_CANCER_MIN = 3
HYPER_THRESH   = 0.3
DELTA_TISSUE   = 0.10
FDR_THRESH     = 0.05


def hdr(msg):
    print(f"\n{'='*70}\n  {msg}\n{'='*70}")


def classify_sample(barcode):
    parts = barcode.replace('.', '-').split('-')
    if len(parts) >= 4:
        try:
            st = int(parts[3][:2])
            if 1 <= st <= 9: return 'tumor'
            if 10 <= st <= 19: return 'normal'
        except ValueError: pass
    return 'unknown'


def bh_fdr(pvals):
    pv = np.array(pvals, dtype=float)
    valid = ~np.isnan(pv)
    result = np.full(len(pv), np.nan)
    if valid.sum() == 0: return result
    idx = np.where(valid)[0]
    pv_valid = pv[idx]
    order = np.argsort(pv_valid)
    rank = np.empty_like(order); rank[order] = np.arange(1, len(pv_valid)+1)
    adjusted = pv_valid * len(pv_valid) / rank
    adjusted = np.minimum.accumulate(adjusted[np.argsort(-rank)])[np.argsort(np.argsort(-rank))]
    result[idx] = np.clip(adjusted, 0, 1)
    return result


def main():
    hdr("TCGA CROSS-CANCER VALIDATION")

    # Load genes
    hdr("STEP 1: Load priority genes")
    if PRIORITY_GENES.exists():
        all_genes = sorted(set(
            g.strip() for g in open(PRIORITY_GENES) if g.strip()))
    else:
        print(f"  WARNING: {PRIORITY_GENES} not found")
        print("  Provide a gene list file (one gene per line)")
        sys.exit(1)
    print(f"  {len(all_genes)} genes to validate")

    # Load probe mapping
    hdr("STEP 2: Load probe mapping")
    pg = pd.read_csv(PROBE_CACHE, sep='\t')
    prom = pg[pg['region'].isin(['TSS200', 'TSS1500'])]
    target = prom[prom['gene'].isin(all_genes)]
    print(f"  {len(target):,} promoter probes → {target['gene'].nunique()} genes")

    probe2genes = {}
    for _, r in target.iterrows():
        probe2genes.setdefault(r['probe_id'], set()).add(r['gene'])

    # Download + process TCGA datasets
    hdr("STEP 3: Process TCGA datasets")
    results = {}

    for cancer, url in TCGA_DATASETS.items():
        fname = TCGA_DIR / f"TCGA_{cancer}_450K.gz"

        # Download if needed
        if not fname.exists() or fname.stat().st_size < 1000:
            print(f"\n  Downloading {cancer}...")
            try:
                urllib.request.urlretrieve(url, str(fname))
                print(f"  Saved: {fname.stat().st_size/1024**2:.0f} MB")
            except Exception as e:
                print(f"  FAILED: {e}")
                continue

        print(f"\n  Processing {cancer} ({fname.stat().st_size/1024**2:.0f} MB)...")
        opener = gzip.open if str(fname).endswith('.gz') else open

        with opener(fname, 'rt') as f:
            header = f.readline().strip().split('\t')
        samples = header[1:]
        stypes = {s: classify_sample(s) for s in samples}
        tumor_idx = [i for i, s in enumerate(samples) if stypes[s] == 'tumor']
        normal_idx = [i for i, s in enumerate(samples) if stypes[s] == 'normal']
        print(f"  Samples: {len(tumor_idx)} tumor, {len(normal_idx)} normal")

        gene_tumor = {g: [] for g in all_genes}
        gene_normal = {g: [] for g in all_genes}

        with opener(fname, 'rt') as f:
            f.readline()
            for line in f:
                parts = line.strip().split('\t')
                probe = parts[0]
                if probe not in probe2genes: continue
                vals = np.array([float(v) if v not in ('', 'NA', 'null', 'nan') else np.nan
                                 for v in parts[1:]])
                for g in probe2genes[probe]:
                    tv = vals[tumor_idx]; tv = tv[~np.isnan(tv)]
                    gene_tumor[g].extend(tv.tolist())
                    if normal_idx:
                        nv = vals[normal_idx]; nv = nv[~np.isnan(nv)]
                        gene_normal[g].extend(nv.tolist())

        # Check for M-values
        sample_vals = []
        for g in all_genes:
            sample_vals.extend(gene_tumor[g][:50])
            if len(sample_vals) > 500: break
        if sample_vals and np.any(np.array(sample_vals) < -0.01):
            print(f"  M-values detected → converting to β")
            for g in all_genes:
                gene_tumor[g]  = [(2**v)/(1+2**v) for v in gene_tumor[g]]
                gene_normal[g] = [(2**v)/(1+2**v) for v in gene_normal[g]]

        # Per-gene stats
        cancer_results = {}
        for g in all_genes:
            t, n = gene_tumor[g], gene_normal[g]
            entry = {
                'tumor_mean': np.mean(t) if t else np.nan,
                'normal_mean': np.mean(n) if n else np.nan,
                'tumor_n': len(t), 'normal_n': len(n),
            }
            if t and n and len(t) > 1 and len(n) > 1:
                entry['delta_tissue'] = entry['tumor_mean'] - entry['normal_mean']
                try: _, pval = stats.ttest_ind(t, n, equal_var=False); entry['p_value'] = pval
                except: entry['p_value'] = np.nan
            else:
                entry['delta_tissue'] = entry['p_value'] = np.nan
            cancer_results[g] = entry
        results[cancer] = cancer_results

    # Build cross-reference table
    hdr("STEP 4: Build validation table")

    # Load blood data
    blood_data = {}
    if BLOOD_RESULTS.exists():
        comp = pd.read_excel(BLOOD_RESULTS, sheet_name='Full_Comparison')
        for _, r in comp.iterrows():
            g = r.get('gene')
            if pd.notna(g):
                blood_data[g] = {
                    'mean_beta': r.get('mean_beta', np.nan),
                    'h_ref_mean': r.get('h_prom_mean', np.nan),
                    'delta_beta': r.get('delta_beta', np.nan),
                }

    rows = []
    for gene in all_genes:
        row = {'gene': gene}
        if gene in blood_data:
            row['blood_bc_beta'] = blood_data[gene]['mean_beta']
            row['blood_healthy_beta'] = blood_data[gene]['h_ref_mean']
            row['blood_delta'] = blood_data[gene]['delta_beta']

        if 'BRCA' in results and gene in results['BRCA']:
            brca = results['BRCA'][gene]
            row['tcga_tumor_beta'] = brca['tumor_mean']
            row['tcga_normal_beta'] = brca['normal_mean']
            row['tcga_delta_tissue'] = brca['delta_tissue']
            row['tcga_pvalue'] = brca['p_value']

        hyper_cancers = []
        for cancer in TCGA_DATASETS:
            if cancer == 'BRCA': continue
            if cancer in results and gene in results[cancer]:
                cr = results[cancer][gene]
                if (pd.notna(cr['tumor_mean']) and cr['tumor_mean'] >= HYPER_THRESH
                    and pd.notna(cr['delta_tissue']) and cr['delta_tissue'] >= DELTA_TISSUE):
                    hyper_cancers.append(cancer)
                row[f'{cancer}_tumor'] = cr['tumor_mean']
                row[f'{cancer}_delta'] = cr['delta_tissue']

        row['pan_cancer_count'] = len(hyper_cancers)
        row['pan_cancer_types'] = '; '.join(hyper_cancers)
        rows.append(row)

    df = pd.DataFrame(rows)
    if 'tcga_pvalue' in df.columns:
        df['tcga_fdr'] = bh_fdr(df['tcga_pvalue'].values)

    # Classify
    hdr("STEP 5: Classify validation")
    def classify(row):
        blood_hyper = pd.notna(row.get('blood_delta')) and row['blood_delta'] > 0.10
        tcga_hyper = (pd.notna(row.get('tcga_delta_tissue')) and row['tcga_delta_tissue'] > DELTA_TISSUE
                      and pd.notna(row.get('tcga_fdr')) and row['tcga_fdr'] < FDR_THRESH)
        tcga_high = pd.notna(row.get('tcga_tumor_beta')) and row['tcga_tumor_beta'] >= HYPER_THRESH
        pan = row.get('pan_cancer_count', 0) >= PAN_CANCER_MIN
        no_tcga = pd.isna(row.get('tcga_tumor_beta'))
        if no_tcga: return 'NO_TCGA_DATA'
        if blood_hyper and tcga_hyper and pan: return 'PAN_CANCER_CONFIRMED'
        if blood_hyper and tcga_hyper and not pan: return 'BREAST_SPECIFIC_CONFIRMED'
        if blood_hyper and tcga_high and not tcga_hyper: return 'TISSUE_ELEVATED'
        if blood_hyper and not tcga_high: return 'BLOOD_ONLY_SIGNAL'
        if not blood_hyper and tcga_hyper: return 'TISSUE_ONLY'
        return 'NOT_VALIDATED'

    df['validation'] = df.apply(classify, axis=1)
    print(df['validation'].value_counts().to_string())

    # Save
    hdr("STEP 6: Save")
    with pd.ExcelWriter(str(OUTPUT), engine='openpyxl') as w:
        df.to_excel(w, 'All_Validation', index=False)
        for cat in ['BREAST_SPECIFIC_CONFIRMED','PAN_CANCER_CONFIRMED','BLOOD_ONLY_SIGNAL']:
            sub = df[df['validation']==cat]
            if not sub.empty:
                sheet = cat.replace('_',' ')[:31]
                sub.to_excel(w, sheet, index=False)
        summary = pd.DataFrame({
            'Category': df['validation'].value_counts().index,
            'Count': df['validation'].value_counts().values
        })
        summary.to_excel(w, 'Summary', index=False)

    print(f"  Saved: {OUTPUT}")
    hdr("DONE")


if __name__ == '__main__':
    main()

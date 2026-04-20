

import pandas as pd
import numpy as np
import os
import argparse
import warnings
warnings.filterwarnings('ignore')


WORK_DIR     = "."
METH_MATRIX  = os.path.join(WORK_DIR, "promoter_unweighted_meth_matrix.tsv")
PROMO_ANNO   = os.path.join(WORK_DIR, "promoter_annotation.tsv")
CLASS_FILE   = os.path.join(WORK_DIR, "promoter_classification.tsv")
BC_GENES_FILE = None  # Optional: path to BC gene list (one per line) or Excel

OUTPUT_DIR   = os.path.join(WORK_DIR, "output", "methylation")
MIN_CPGS     = 10  # Minimum CpGs per promoter window

os.makedirs(OUTPUT_DIR, exist_ok=True)

# BC susceptibility genes (default panel — customize as needed)
DEFAULT_BC_GENES = [
    "BRCA1","BRCA2","TP53","PTEN","CDH1","STK11","PALB2","CHEK2","ATM","NBN",
    "NF1","BARD1","BRIP1","RAD51C","RAD51D","MLH1","MSH2","MSH6","PMS2","EPCAM",
    "PIK3CA","AKT1","ESR1","GATA3","FOXA1","MAP3K1","CDK12","RB1","ARID1A","CBFB",
    "RUNX1","ERBB2","MYC","CCND1","FGFR1","FGFR2","NOTCH1","NOTCH2"
]


def hdr(msg):
    print(f"\n{'='*70}\n  {msg}\n{'='*70}")


def interpret(cls, beta):
    """Interpret methylation level based on promoter classification."""
    if pd.isna(cls) or cls == 'unclassified':
        return "Unclassified — needs manual verification"
    elif cls == 'gene_body':
        if beta >= 0.8: return "Gene body: high methylation (active transcription)"
        elif beta >= 0.4: return "Gene body: intermediate methylation"
        else: return "Gene body: low methylation"
    elif cls == 'alternative':
        if beta >= 0.5: return "Alternative promoter: high methylation (likely silenced)"
        else: return "Alternative promoter: low/intermediate"
    elif cls == 'canonical':
        if beta >= 0.8: return "HYPERMETHYLATED — potential silencing"
        elif beta >= 0.5: return "Substantial methylation at canonical promoter"
        elif beta >= 0.2: return "Partial methylation at canonical promoter"
        elif beta >= 0.05: return "Low methylation (normal range)"
        else: return "Unmethylated (active gene)"
    return ""


def main():
    hdr("METHYLATION MASTER PIPELINE")

    # --- Load data ---
    hdr("STEP 1: Load data")
    meth = pd.read_csv(METH_MATRIX, sep='\t')
    if 'gene' not in meth.columns:
        meth.insert(1, 'gene', meth['promoter_name'].str.rsplit('_', n=1).str[0])

    sample_cols = [c for c in meth.columns if c not in ['promoter_name', 'gene']]
    blood_cols  = [s for s in sample_cols if not s.endswith('s')]
    saliva_cols = [s for s in sample_cols if s.endswith('s')]
    print(f"  {meth.shape[0]:,} promoters, {len(sample_cols)} samples "
          f"({len(blood_cols)} blood, {len(saliva_cols)} saliva)")

    anno = pd.read_csv(PROMO_ANNO, sep='\t')
    classif = pd.read_csv(CLASS_FILE, sep='\t')

    # Detect CpG count column
    cpg_col = None
    for c in ['n_cpgs', 'mean_n_cpgs', 'min_n_cpgs', 'n_probes', 'num_cpgs']:
        if c in anno.columns:
            cpg_col = c; break
    if cpg_col:
        print(f"  CpG column: '{cpg_col}', median={anno[cpg_col].median():.0f}")

    # BC genes
    bc_genes = DEFAULT_BC_GENES
    if BC_GENES_FILE and os.path.exists(BC_GENES_FILE):
        bc_genes = [g.strip() for g in open(BC_GENES_FILE) if g.strip()]

    # --- Merge ---
    hdr("STEP 2: Merge and filter")
    merge_cols = [c for c in ['promoter_name','chr','start','end','tss',
                               'has_cpg_island','cpg_island_overlap_bp',
                               'position_vs_gene','classification',
                               'distance_to_gene_tss'] if c in classif.columns]
    merged = meth.merge(classif[merge_cols], on='promoter_name', how='left')

    if cpg_col:
        merged = merged.merge(
            anno[['promoter_name', cpg_col]].rename(columns={cpg_col: 'n_cpgs'}),
            on='promoter_name', how='left')
        n_before = len(merged)
        merged = merged[merged['n_cpgs'] >= MIN_CPGS].copy()
        print(f"  CpG filter (>={MIN_CPGS}): {n_before:,} → {len(merged):,}")
    else:
        merged['n_cpgs'] = np.nan

    # --- Statistics ---
    hdr("STEP 3: Calculate statistics")
    merged['mean_beta']   = merged[sample_cols].mean(axis=1)
    merged['median_beta'] = merged[sample_cols].median(axis=1)
    merged['sd_beta']     = merged[sample_cols].std(axis=1)
    merged['min_beta']    = merged[sample_cols].min(axis=1)
    merged['max_beta']    = merged[sample_cols].max(axis=1)
    merged['range_beta']  = merged['max_beta'] - merged['min_beta']
    merged['n_samples']   = merged[sample_cols].notna().sum(axis=1)
    merged['cv'] = merged['sd_beta'] / merged['mean_beta'].replace(0, np.nan)

    merged['n_patients_beta_ge_0.2'] = (merged[sample_cols] >= 0.2).sum(axis=1)
    merged['n_patients_beta_ge_0.5'] = (merged[sample_cols] >= 0.5).sum(axis=1)
    merged['pct_patients_beta_ge_0.2'] = (
        merged['n_patients_beta_ge_0.2'] / merged['n_samples'] * 100).round(1)
    merged['pct_patients_beta_ge_0.5'] = (
        merged['n_patients_beta_ge_0.5'] / merged['n_samples'] * 100).round(1)

    if blood_cols:
        merged['mean_beta_blood'] = merged[blood_cols].mean(axis=1)
    if saliva_cols:
        merged['mean_beta_saliva'] = merged[saliva_cols].mean(axis=1)

    merged['is_bc_gene'] = merged['gene'].isin(bc_genes) if 'gene' in merged.columns else False
    merged['interpretation'] = merged.apply(
        lambda r: interpret(r.get('classification'), r['mean_beta']), axis=1)

    # --- Build sheets ---
    hdr("STEP 4: Build output sheets")

    meta = [c for c in ['promoter_name','gene','chr','start','end','classification',
                         'position_vs_gene','has_cpg_island','distance_to_gene_tss',
                         'n_cpgs','is_bc_gene'] if c in merged.columns]
    stats = [c for c in ['mean_beta','median_beta','sd_beta','min_beta','max_beta',
                          'range_beta','cv','n_samples','n_patients_beta_ge_0.2',
                          'pct_patients_beta_ge_0.2','n_patients_beta_ge_0.5',
                          'pct_patients_beta_ge_0.5','mean_beta_blood',
                          'mean_beta_saliva','interpretation'] if c in merged.columns]

    def make_sheet(df):
        cols = [c for c in meta + stats if c in df.columns]
        return df[cols].sort_values('mean_beta', ascending=False).reset_index(drop=True)

    canonical = merged[(merged['classification']=='canonical') & (merged['mean_beta']>=0.2)]
    gene_body = merged[merged['classification']=='gene_body']
    hyper     = merged[(merged['classification']=='canonical') & (merged['mean_beta']>=0.5)]
    high_cv   = merged[(merged['mean_beta']>=0.1) & (merged['cv']>=0.15)]

    bc_all    = merged[merged['is_bc_gene']]
    bc_canon  = bc_all[bc_all['classification']=='canonical']
    bc_gb     = bc_all[bc_all['classification']=='gene_body']

    print(f"  Canonical hits (β≥0.2): {len(canonical):,}")
    print(f"  Gene body: {len(gene_body):,}")
    print(f"  Hypermethylated (β≥0.5): {len(hyper):,}")
    print(f"  High CV: {len(high_cv):,}")
    print(f"  BC canonical: {len(bc_canon)}, BC gene body: {len(bc_gb)}")

    # Gene-level summary
    gene_rows = []
    for gene in sorted(merged['gene'].unique()):
        gp = merged[merged['gene']==gene]
        canon = gp[gp['classification']=='canonical']
        if len(canon) > 0:
            best = canon.iloc[canon['mean_beta'].argmax()]
            src = 'canonical'
        else:
            best = gp.iloc[0]
            src = best.get('classification','unknown')
        gene_rows.append({
            'gene': gene, 'is_bc_gene': gene in bc_genes,
            'primary_promoter': best['promoter_name'], 'promoter_class': src,
            'chr': best.get('chr',''), 'mean_beta': best['mean_beta'],
            'sd_beta': best['sd_beta'], 'cv': best.get('cv', np.nan),
            'n_total_promoters': len(gp),
            'n_canonical': len(canon),
            'interpretation': best.get('interpretation','')
        })
    gene_summary = pd.DataFrame(gene_rows).sort_values('mean_beta', ascending=False)

    # QC summary
    qc = pd.DataFrame([
        {'Metric': 'Total promoters (pre-filter)', 'Value': len(meth)},
        {'Metric': f'After CpG filter (>={MIN_CPGS})', 'Value': len(merged)},
        {'Metric': 'Canonical promoters', 'Value': (merged['classification']=='canonical').sum()},
        {'Metric': 'Gene body promoters', 'Value': (merged['classification']=='gene_body').sum()},
        {'Metric': 'Unique genes', 'Value': merged['gene'].nunique()},
        {'Metric': 'Total samples', 'Value': len(sample_cols)},
        {'Metric': 'Blood samples', 'Value': len(blood_cols)},
        {'Metric': 'Saliva samples', 'Value': len(saliva_cols)},
        {'Metric': 'Canonical with β≥0.2', 'Value': len(canonical)},
        {'Metric': 'Canonical with β≥0.5', 'Value': len(hyper)},
    ])

    # --- Save ---
    hdr("STEP 5: Save Excel files")

    genome_file = os.path.join(OUTPUT_DIR, "Methylation_Analysis_Master.xlsx")
    with pd.ExcelWriter(genome_file, engine='openpyxl') as w:
        make_sheet(canonical).to_excel(w, 'Canonical_Hits', index=False)
        make_sheet(gene_body).to_excel(w, 'Gene_Body_Track', index=False)
        make_sheet(hyper).to_excel(w, 'Hypermethylated_Canonical', index=False)
        make_sheet(high_cv).to_excel(w, 'High_Variability', index=False)
        gene_summary.to_excel(w, 'Gene_Level_Summary', index=False)
        qc.to_excel(w, 'QC_Summary', index=False)
    print(f"  Saved: {genome_file}")

    bc_file = os.path.join(OUTPUT_DIR, "BC_Genes_Methylation_Master.xlsx")
    with pd.ExcelWriter(bc_file, engine='openpyxl') as w:
        make_sheet(bc_canon).to_excel(w, 'BC_Canonical', index=False)
        make_sheet(bc_gb).to_excel(w, 'BC_Gene_Body', index=False)
        make_sheet(bc_all).to_excel(w, 'BC_All_Promoters', index=False)
        gene_summary[gene_summary['is_bc_gene']].to_excel(w, 'BC_Gene_Summary', index=False)
    print(f"  Saved: {bc_file}")

    hdr("PIPELINE COMPLETE")


if __name__ == '__main__':
    main()

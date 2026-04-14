# Genomic and Epigenomic Characterization of Breast Cancer in Peripheral DNA

Multi-omics analysis combining germline whole-exome sequencing (WES) and genome-wide promoter methylation profiling in a UAE breast cancer cohort.

## Study Overview

- **Cohort**: 94 breast cancer patients, 415 cancer-free controls (United Arab Emirates)
- **WES**: Germline variant discovery using dual-pipeline approach (Case-Only + Background Subtraction)
- **Methylation**: Promoter-level DNA methylation (Illumina EPIC 850K) with healthy blood reference comparison and TCGA tissue validation

## Repository Structure

```
├── wes/                          # Whole-Exome Sequencing Pipelines
│   ├── 01_pipeline_a_case_only.R           # Pipeline A: Case-only rare deleterious model
│   ├── 02_pipeline_b_background.R          # Pipeline B: Control-derived background subtraction
│   └── 03_tier_classification.R            # Five-tier variant classification + Excel output
│
├── methylation/                  # Methylation Analysis Pipeline
│   ├── 01_methylation_master.py            # EPD promoter annotation, classification, QC
│   ├── 02_control_comparison.py            # Healthy blood (GSE40279) comparison + Δβ
│   ├── 03_tcga_validation.py               # TCGA cross-cancer tissue validation
│   └── 04_priority_tiers.py                # Priority investigation tier classification
│
├── figures/                      # Publication Figure Scripts
│   ├── fig_bc300_gene_panel.py             # BC300 gene panel 4-panel figure
│   └── fig_pca_blood_saliva.py             # PCA blood vs saliva (Supp. Fig S1)
│
└── data/                         # Example input file descriptions
    └── INPUT_FILES.md                      # Required input file formats
```

## Prerequisites

### R packages
```r
install.packages(c("tidyverse", "data.table", "openxlsx"))
```

### Python packages
```bash
pip install pandas numpy scipy openpyxl scikit-learn matplotlib
```

### External Tools & Databases
| Tool | Version | Reference |
|------|---------|-----------|
| ANNOVAR | 2025Mar02 | Wang et al., 2010 (PMID: 20601685) |
| BWA-MEM | v0.7.17 | Li, 2013 (arXiv:1303.3997) |
| GATK HaplotypeCaller | v4.2 | McKenna et al., 2010 (PMID: 20644199) |
| gnomAD | v4.1 | Karczewski et al., 2020 (PMID: 32461654) |
| dbNSFP | v4.7a | Liu et al., 2020 (PMID: 33261662) |
| ClinVar | clinvar_20250721 | Landrum et al., 2018 (PMID: 29165669) |
| REVEL | — | Ioannidis et al., 2016 (PMID: 27666373) |
| CADD | — | Rentzsch et al., 2019 (PMID: 30371827) |
| SIFT4G | — | Vaser et al., 2016 (PMID: 26633127) |
| PolyPhen-2 | — | Adzhubei et al., 2010 (PMID: 20354512) |
| DANN | — | Quang et al., 2015 (PMID: 25338716) |

## Usage

### 1. WES Analysis
```bash
# Step 1: Run Pipeline A (case-only filtering)
Rscript wes/01_pipeline_a_case_only.R

# Step 2: Run Pipeline B (background subtraction)
Rscript wes/02_pipeline_b_background.R

# Step 3: Combine and classify into tiers
Rscript wes/03_tier_classification.R
```

### 2. Methylation Analysis
```bash
# Step 1: Promoter annotation and master pipeline
python methylation/01_methylation_master.py

# Step 2: Healthy blood control comparison
python methylation/02_control_comparison.py

# Step 3: TCGA tissue validation
python methylation/03_tcga_validation.py

# Step 4: Priority tier classification
python methylation/04_priority_tiers.py
```

### 3. Generate Figures
```bash
python figures/fig_bc300_gene_panel.py
python figures/fig_pca_blood_saliva.py
```

## Configuration

All scripts use a `CONFIG` section at the top. Modify paths and parameters before running:

```python
# Python scripts
WORK_DIR = "/path/to/your/working/directory"
```

```r
# R scripts
setwd("/path/to/your/working/directory")
```

## Citation

If you use this code, please cite:

> Alnaqbi H, Zayed N, Al Meqdam F, et al. Genomic and epigenomic characterization of breast cancer in peripheral DNA. *[Journal]*, 2026.

## License

This project is licensed under the MIT License.

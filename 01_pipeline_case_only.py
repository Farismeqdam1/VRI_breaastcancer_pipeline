library(tidyverse)
library(data.table)
library(openxlsx)

# =============================================================================
# CONFIGURATION — MODIFY THESE PATHS
# =============================================================================

WORK_DIR <- "."                    # Working directory
setwd(WORK_DIR)

ANNOVAR_FILE <- "VRI_II_annotated.hg38_multianno.txt"   # ANNOVAR output
CASE_SAMPLE_LIST <- "case_samples.txt"                   # One sample ID per line
ALL_SAMPLE_LIST  <- "sample_list.txt"                    # All samples in VCF

OUTPUT_DIR <- "output/pipeline_a"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Filtering thresholds
GNOMAD_AF_THRESHOLD <- 0.01       # gnomAD global MAF < 1%
REVEL_DELETERIOUS   <- 0.5        # REVEL >= 0.5 for "predicted deleterious"
REVEL_HIGH          <- 0.75       # REVEL >= 0.75 for high-confidence VUS
CADD_THRESHOLD      <- 20         # CADD PHRED >= 20
DANN_THRESHOLD      <- 0.95       # DANN >= 0.95
SIFT_THRESHOLD      <- 0.05       # SIFT4G <= 0.05 (deleterious)
PP2_THRESHOLD       <- 0.85       # PolyPhen-2 >= 0.85 (probably damaging)
MIN_PREDICTORS      <- 3          # Minimum concordant deleterious predictions

# Breast cancer susceptibility gene panel (customizable)
BC_GENES <- c(
  "BRCA1", "BRCA2", "TP53", "PTEN", "CDH1", "STK11",
  "PALB2", "CHEK2", "ATM", "NBN", "NF1", "BARD1", "BRIP1", "RAD51C", "RAD51D",
  "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM",
  "PIK3CA", "AKT1", "ESR1", "GATA3", "FOXA1", "MAP3K1", "CDK12", "RB1",
  "ARID1A", "CBFB", "RUNX1", "ERBB2", "MYC", "CCND1", "FGFR1", "FGFR2",
  "NOTCH1", "NOTCH2"
)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("=============================================================\n")
cat("   PIPELINE A: Case-Only Rare Deleterious Variant Model      \n")
cat("=============================================================\n\n")

cat("Loading ANNOVAR annotated file...\n")
dt <- fread(ANNOVAR_FILE, sep = "\t", header = TRUE, na.strings = c(".", "NA", ""))
df <- as.data.frame(dt); rm(dt); gc()
cat("Loaded", nrow(df), "variants\n\n")

# Load sample lists
case_samples <- readLines(CASE_SAMPLE_LIST)
all_samples  <- readLines(ALL_SAMPLE_LIST)
cat("Case samples:", length(case_samples), "\n\n")

# Map sample names to ANNOVAR Otherinfo columns
sample_to_col <- setNames(paste0("Otherinfo", 12 + seq_along(all_samples)), all_samples)
case_cols <- sample_to_col[case_samples]
case_cols <- case_cols[!is.na(case_cols)]

# =============================================================================
# STEP 1: FILTER VARIANTS PRESENT IN CASES
# =============================================================================

cat("=== Step 1: Filter to variants present in cases ===\n")

bc_present <- rep(FALSE, nrow(df))
for (col in case_cols) {
  if (col %in% names(df)) {
    gt <- df[[col]]
    bc_present <- bc_present | (!is.na(gt) & grepl("1", gt))
  }
}
df_bc <- df[bc_present, ]
cat("Variants in cases:", nrow(df_bc), "\n\n")
rm(df); gc()

# =============================================================================
# STEP 2: EXONIC / SPLICING ONLY
# =============================================================================

cat("=== Step 2: Filter to exonic/splicing ===\n")
df_exonic <- df_bc %>%
  filter(Func.refGene %in% c("exonic", "exonic;splicing", "splicing"))
cat("Exonic/splicing:", nrow(df_exonic), "\n\n")

# =============================================================================
# STEP 3: REMOVE SYNONYMOUS
# =============================================================================

cat("=== Step 3: Remove synonymous ===\n")
df_nonsyn <- df_exonic %>%
  filter(ExonicFunc.refGene != "synonymous SNV" | is.na(ExonicFunc.refGene))
cat("Non-synonymous:", nrow(df_nonsyn), "\n\n")

# =============================================================================
# STEP 4: ALLELE FREQUENCY FILTER (gnomAD < 1%)
# =============================================================================

cat("=== Step 4: gnomAD AF filter (<", GNOMAD_AF_THRESHOLD * 100, "%) ===\n")

# Auto-detect gnomAD AF column
gnomad_col <- NULL
for (candidate in c("gnomad41_genome_AF", "AF", "gnomAD_genome_ALL")) {
  if (candidate %in% names(df_nonsyn)) {
    gnomad_col <- candidate
    break
  }
}

if (!is.null(gnomad_col)) {
  df_nonsyn[[gnomad_col]] <- as.numeric(as.character(df_nonsyn[[gnomad_col]]))
  df_rare <- df_nonsyn %>%
    filter(is.na(.data[[gnomad_col]]) | .data[[gnomad_col]] <= GNOMAD_AF_THRESHOLD)
  cat("Using column:", gnomad_col, "\n")
} else {
  warning("No gnomAD AF column found — skipping AF filter")
  df_rare <- df_nonsyn
}
cat("Rare variants:", nrow(df_rare), "\n\n")

# =============================================================================
# STEP 5: PATHOGENICITY PREDICTION
# =============================================================================

cat("=== Step 5: Multi-tool pathogenicity scoring ===\n")

score_cols <- c("REVEL_score", "CADD_phred", "SIFT4G_score",
                "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "DANN_score")
for (col in score_cols) {
  if (col %in% names(df_rare)) {
    df_rare[[col]] <- as.numeric(as.character(df_rare[[col]]))
  }
}

df_rare <- df_rare %>%
  mutate(
    is_REVEL  = !is.na(REVEL_score) & REVEL_score >= REVEL_DELETERIOUS,
    is_CADD   = !is.na(CADD_phred) & CADD_phred >= CADD_THRESHOLD,
    is_DANN   = !is.na(DANN_score) & DANN_score >= DANN_THRESHOLD,
    is_SIFT   = !is.na(SIFT4G_score) & SIFT4G_score <= SIFT_THRESHOLD,
    is_PP2_HD = !is.na(Polyphen2_HDIV_score) & Polyphen2_HDIV_score >= PP2_THRESHOLD,
    is_PP2_HV = !is.na(Polyphen2_HVAR_score) & Polyphen2_HVAR_score >= PP2_THRESHOLD
  ) %>%
  mutate(
    pathogenic_count = rowSums(select(., starts_with("is_")) == TRUE, na.rm = TRUE)
  )

# Keep if REVEL >= 0.5 OR >= 3 concordant predictors
df_pathogenic <- df_rare %>%
  filter(pathogenic_count >= MIN_PREDICTORS | is_REVEL == TRUE)

cat("Pathogenic (>=3 tools OR REVEL>=0.5):", nrow(df_pathogenic), "\n\n")

# =============================================================================
# STEP 6: CLINVAR PATHOGENIC
# =============================================================================

cat("=== Step 6: ClinVar pathogenic variants ===\n")

clinvar_terms <- c("Pathogenic", "Pathogenic/Likely_pathogenic",
                   "Likely_pathogenic", "Pathogenic|risk_factor")
df_clinvar <- df_rare %>% filter(CLNSIG %in% clinvar_terms)
cat("ClinVar P/LP:", nrow(df_clinvar), "\n\n")

# =============================================================================
# STEP 7: BC GENE PANEL
# =============================================================================

cat("=== Step 7: BC gene panel variants ===\n")

df_bc_genes <- df_pathogenic %>%
  separate_rows(Gene.refGene, sep = ";") %>%
  filter(Gene.refGene %in% BC_GENES)
cat("Variants in BC genes:", nrow(df_bc_genes), "\n\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("=== Saving outputs ===\n")

write_tsv(df_rare,       file.path(OUTPUT_DIR, "filtered_rare_variants.tsv"))
write_tsv(df_pathogenic, file.path(OUTPUT_DIR, "pathogenic_variants.tsv"))
write_tsv(df_clinvar,    file.path(OUTPUT_DIR, "clinvar_pathogenic.tsv"))

# Gene lists
genes_all <- df_rare %>% separate_rows(Gene.refGene, sep=";") %>%
  pull(Gene.refGene) %>% unique() %>% na.omit()
genes_path <- df_pathogenic %>% separate_rows(Gene.refGene, sep=";") %>%
  pull(Gene.refGene) %>% unique() %>% na.omit()

writeLines(genes_all,  file.path(OUTPUT_DIR, "genes_all_filtered.txt"))
writeLines(genes_path, file.path(OUTPUT_DIR, "genes_pathogenic.txt"))

saveRDS(df_rare,       file.path(OUTPUT_DIR, "filtered_rare.RDS"))
saveRDS(df_pathogenic, file.path(OUTPUT_DIR, "pathogenic.RDS"))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=======================================================\n")
cat("                  PIPELINE A SUMMARY                   \n")
cat("=======================================================\n")
cat("Variants in cases:          ", nrow(df_bc), "\n")
cat("Exonic/splicing:            ", nrow(df_exonic), "\n")
cat("Non-synonymous:             ", nrow(df_nonsyn), "\n")
cat("Rare (AF <=", GNOMAD_AF_THRESHOLD, "):     ", nrow(df_rare), "\n")
cat("Pathogenic (multi-tool):    ", nrow(df_pathogenic), "\n")
cat("ClinVar P/LP:               ", nrow(df_clinvar), "\n")
cat("BC gene panel variants:     ", nrow(df_bc_genes), "\n")
cat("Unique genes (all):         ", length(genes_all), "\n")
cat("Unique genes (pathogenic):  ", length(genes_path), "\n")
cat("=======================================================\n")
cat("Done!\n")

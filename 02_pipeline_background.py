library(tidyverse)
library(data.table)

# CONFIGURATION — MODIFY THESE PATHS

WORK_DIR <- "."
setwd(WORK_DIR)

ANNOVAR_FILE         <- "VRI_II_annotated.hg38_multianno.txt"
CASE_SAMPLE_LIST     <- "case_samples.txt"
CONTROL_SAMPLE_LIST  <- "control_samples.txt"
ALL_SAMPLE_LIST      <- "sample_list.txt"

OUTPUT_DIR <- "output/pipeline_b"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Thresholds
MENA_AF_THRESHOLD        <- 0.005    # Keep if MENA AF >= 0.5%
REVEL_STRONG_THRESHOLD   <- 0.75     # "Strongly damaging"
CONTROL_COMMON_THRESHOLD <- 0.10     # >= 10% of controls = "common"

# Cohort sizes (auto-detected from sample lists)
N_CASES    <- NULL  # Set to NULL for auto-detection
N_CONTROLS <- NULL

# HELPER FUNCTIONS

has_variant <- function(gt) {
  if (is.na(gt)) return(FALSE)
  grepl("1", gt)
}

count_carriers <- function(df_subset, sample_cols) {
  apply(df_subset[, ..sample_cols], 1, function(row) {
    sum(sapply(row, has_variant), na.rm = TRUE)
  })
}


# LOAD DATA

cat("Loading data...\n")
df <- fread(ANNOVAR_FILE, sep = "\t", header = TRUE, na.strings = c(".", "", "NA"))

case_samples    <- readLines(CASE_SAMPLE_LIST)
control_samples <- readLines(CONTROL_SAMPLE_LIST)
all_samples     <- readLines(ALL_SAMPLE_LIST)

if (is.null(N_CASES))    N_CASES    <- length(case_samples)
if (is.null(N_CONTROLS)) N_CONTROLS <- length(control_samples)

MIN_CONTROL_CARRIERS <- ceiling(N_CONTROLS * CONTROL_COMMON_THRESHOLD)

cat("Loaded", nrow(df), "total variants\n")
cat("Cases:", N_CASES, "| Controls:", N_CONTROLS, "\n\n")

# Map sample columns
case_cols    <- paste0("Otherinfo", 12 + match(case_samples, all_samples))
control_cols <- paste0("Otherinfo", 12 + match(control_samples, all_samples))

# STEP 0: COUNT CARRIERS

cat("Step 0: Counting carriers (this may take several minutes)...\n")
df$n_case_carriers    <- count_carriers(df, case_cols)
df$n_control_carriers <- count_carriers(df, control_cols)
df$case_carrier_freq    <- df$n_case_carriers / N_CASES
df$control_carrier_freq <- df$n_control_carriers / N_CONTROLS
cat("Done.\n\n")

# STEP 1: MENA ALLELE FREQUENCY FILTER

cat("=== Step 1: MENA AF filter (>= ", MENA_AF_THRESHOLD * 100, "%) ===\n")

df_controls <- df %>% filter(n_control_carriers > 0)

# Auto-detect AF column
af_col <- NULL
for (candidate in c("AF", "gnomad41_genome_AF", "gnomAD_genome_ALL", "AF_popmax")) {
  if (candidate %in% names(df_controls)) { af_col <- candidate; break }
}

if (!is.null(af_col)) {
  df_controls$MENA_AF <- as.numeric(df_controls[[af_col]])
  cat("Using AF column:", af_col, "\n")
} else {
  df_controls$MENA_AF <- df_controls$control_carrier_freq
  cat("WARNING: No AF column found — using control frequency as proxy\n")
}
df_controls$MENA_AF[is.na(df_controls$MENA_AF)] <- 0

df_step1 <- df_controls %>%
  filter(MENA_AF >= MENA_AF_THRESHOLD | control_carrier_freq >= CONTROL_COMMON_THRESHOLD)

cat("Retained:", nrow(df_step1), "| Removed:", nrow(df_controls) - nrow(df_step1), "\n\n")

# STEP 2: CLINVAR FILTER

cat("=== Step 2: ClinVar classification filter ===\n")

df_step1$clinvar_class <- case_when(
  grepl("Pathogenic", df_step1$CLNSIG, ignore.case = TRUE) &
    !grepl("Conflicting|Benign", df_step1$CLNSIG, ignore.case = TRUE) ~ "P_LP",
  grepl("Likely_pathogenic", df_step1$CLNSIG, ignore.case = TRUE) &
    !grepl("Conflicting|Benign", df_step1$CLNSIG, ignore.case = TRUE) ~ "P_LP",
  grepl("Benign", df_step1$CLNSIG, ignore.case = TRUE) &
    !grepl("Pathogenic|Conflicting", df_step1$CLNSIG, ignore.case = TRUE) ~ "B_LB",
  grepl("Likely_benign", df_step1$CLNSIG, ignore.case = TRUE) &
    !grepl("Pathogenic|Conflicting", df_step1$CLNSIG, ignore.case = TRUE) ~ "B_LB",
  TRUE ~ "VUS_other"
)

df_benign     <- df_step1 %>% filter(clinvar_class == "B_LB")
df_plp_removed <- df_step1 %>% filter(clinvar_class == "P_LP")
df_vus        <- df_step1 %>% filter(clinvar_class == "VUS_other")

cat("Benign/LB (kept):", nrow(df_benign), "\n")
cat("P/LP (removed from background):", nrow(df_plp_removed), "\n")
cat("VUS (to Step 3):", nrow(df_vus), "\n\n")


# STEP 3: REVEL + CONTROL FREQUENCY FOR VUS

cat("=== Step 3: VUS filtering (REVEL + control frequency) ===\n")

df_vus$REVEL_score <- as.numeric(df_vus$REVEL_score)
df_vus$REVEL_score[is.na(df_vus$REVEL_score)] <- 0

df_vus <- df_vus %>%
  mutate(vus_decision = case_when(
    REVEL_score >= REVEL_STRONG_THRESHOLD &
      control_carrier_freq >= CONTROL_COMMON_THRESHOLD ~ "KEEP_common",
    REVEL_score >= REVEL_STRONG_THRESHOLD &
      control_carrier_freq < CONTROL_COMMON_THRESHOLD ~ "REMOVE_rare_damaging",
    TRUE ~ "KEEP_low_REVEL"
  ))

df_vus_keep    <- df_vus %>% filter(vus_decision != "REMOVE_rare_damaging")
df_vus_removed <- df_vus %>% filter(vus_decision == "REMOVE_rare_damaging")

cat("VUS kept:", nrow(df_vus_keep), "| Removed:", nrow(df_vus_removed), "\n\n")

# BUILD HEALTHY BACKGROUND

cat("=== Building Healthy Background Set ===\n")

healthy_background <- bind_rows(
  df_benign  %>% mutate(source = "ClinVar_Benign"),
  df_vus_keep %>% mutate(source = "VUS_passed")
)
healthy_background$var_id <- paste(
  healthy_background$Chr, healthy_background$Start,
  healthy_background$Ref, healthy_background$Alt, sep = "_"
)

cat("Total healthy background:", nrow(healthy_background), "\n\n")

# Save
fwrite(
  healthy_background %>%
    select(Chr, Start, End, Ref, Alt, Gene.refGene, REVEL_score, CLNSIG,
           MENA_AF, n_control_carriers, control_carrier_freq, source, var_id),
  file.path(OUTPUT_DIR, "healthy_background_variants.tsv"), sep = "\t"
)

# SUBTRACT BACKGROUND FROM CASES

cat("=== Subtracting background from cases ===\n")

df_cases <- df %>% filter(n_case_carriers > 0)
df_cases$var_id <- paste(df_cases$Chr, df_cases$Start, df_cases$Ref, df_cases$Alt, sep = "_")

df_candidates <- df_cases %>% filter(!var_id %in% healthy_background$var_id)

cat("Case variants before:", nrow(df_cases), "\n")
cat("After background subtraction:", nrow(df_candidates), "\n")
cat("Removed:", nrow(df_cases) - nrow(df_candidates), "\n\n")

# Save
fwrite(
  df_candidates %>%
    select(Chr, Start, End, Ref, Alt, Gene.refGene, AAChange.refGene,
           REVEL_score, CLNSIG, n_case_carriers, n_control_carriers,
           case_carrier_freq, control_carrier_freq, var_id),
  file.path(OUTPUT_DIR, "case_candidates.tsv"), sep = "\t"
)


# GENE-LEVEL PATIENT COUNTS

cat("=== Gene-level patient counts ===\n")

df$var_id <- paste(df$Chr, df$Start, df$Ref, df$Alt, sep = "_")
df_counting <- df %>% filter(var_id %in% df_candidates$var_id)

gene_counts <- df_counting %>%
  select(Gene.refGene, var_id, all_of(case_cols)) %>%
  pivot_longer(cols = all_of(case_cols), names_to = "sample_col", values_to = "genotype") %>%
  mutate(
    sample_idx = as.numeric(gsub("Otherinfo", "", sample_col)) - 12,
    sample_name = case_samples[sample_idx],
    has_var = sapply(genotype, has_variant)
  ) %>%
  filter(has_var) %>%
  group_by(Gene.refGene) %>%
  summarise(
    n_variants = n_distinct(var_id),
    n_patients = n_distinct(sample_name),
    pct = round(n_patients / N_CASES * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients))

fwrite(gene_counts, file.path(OUTPUT_DIR, "gene_patient_counts.tsv"), sep = "\t")

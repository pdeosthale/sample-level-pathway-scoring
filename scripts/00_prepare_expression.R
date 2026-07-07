# ==============================================================================
# 00_prepare_expression.R
# Download and prepare normalized expression matrix for GSVA/ssGSEA
#
# Part of: sample-level-pathway-scoring
# Dataset: GSE161731 — Blood transcriptomics in respiratory infections
#          McClain et al. 2021, Nat Commun
#
# Input:
#   GSE161731_xpr_nlcpm.csv.gz  — normalized log-CPM (authors' pre-computed)
#   GSE161731_key.csv.gz        — sample metadata
#
# Output:
#   data/processed/expr_matrix.rds    — genes x samples log-CPM matrix
#   data/processed/metadata.rds       — sample metadata, aligned to matrix
#
# Notes:
#   - nlcpm file is samples x genes — transpose to genes x samples for GSVA
#   - Sample IDs use "-" in metadata, need to match column names
#   - One sample per subject (earliest timepoint) to avoid pseudoreplication
# ==============================================================================

library(tidyverse)
library(GEOquery)

setwd("~/Desktop/github_repos/sample-level-pathway-scoring")
dir.create("data/raw",       showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results",        showWarnings = FALSE)
dir.create("figures",        showWarnings = FALSE)


# ==============================================================================
# 1. DOWNLOAD FILES
# ==============================================================================

cat("=== Step 1: Downloading files ===\n")

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161731/suppl/"

files <- list(
  nlcpm = "GSE161731_xpr_nlcpm.csv.gz",
  key   = "GSE161731_key.csv.gz"
)

for (fname in files) {
  dest <- file.path("data/raw", fname)
  if (!file.exists(dest)) {
    cat("Downloading", fname, "...\n")
    options(timeout = 300)
    download.file(paste0(base_url, fname), dest, mode = "wb")
    cat("  Done\n")
  } else {
    cat(" ", fname, "already exists\n")
  }
}


# ==============================================================================
# 2. LOAD EXPRESSION MATRIX
# ==============================================================================

cat("\n=== Step 2: Loading expression matrix ===\n")

# File is samples x genes — transpose to genes x samples
expr_raw <- read.csv("data/raw/GSE161731_xpr_nlcpm.csv.gz",
                     row.names = 1,
                     check.names = FALSE)

cat("Raw dimensions (samples x genes):", dim(expr_raw), "\n")

# Transpose to genes x samples
expr_mat <- t(as.matrix(expr_raw))

cat("After transpose (genes x samples):", dim(expr_mat), "\n")
cat("Example gene names:", head(rownames(expr_mat), 3), "\n")
cat("Example sample IDs:", head(colnames(expr_mat), 3), "\n")


# ==============================================================================
# 3. LOAD AND CLEAN METADATA
# ==============================================================================

cat("\n=== Step 3: Loading metadata ===\n")

key <- read.csv("data/raw/GSE161731_key.csv.gz")
cat("Metadata columns:", colnames(key), "\n")
cat("Cohorts:", unique(key$cohort), "\n")
cat("Total samples in key:", nrow(key), "\n")


# ==============================================================================
# 4. USE DE METADATA AS SAMPLE FILTER
# ==============================================================================
#
# Rather than re-deduplicating from scratch, we use the QC-cleaned metadata
# from the infectious-disease-transcriptomics DE analysis as the authoritative
# sample list. This ensures:
#   - Same outlier removal as the DE analysis
#   - Same one-sample-per-subject deduplication (earliest timepoint)
#   - Consistent sample set across both repositories
#
# Path to DE metadata — adjust if your repo is in a different location
DE_META_PATH <- "~/Desktop/github_repos/infectious-disease-transcriptomics/results/GSE161731_metadata_clean.csv"

cat("\n=== Step 4: Loading DE metadata as sample filter ===\n")

de_meta <- read.csv(DE_META_PATH)
cat("DE metadata samples (before deduplication):", nrow(de_meta), "\n")

# Keep one sample per subject — earliest timepoint
# Matches the approach used in the DE analysis
time_order <- c("early", "middle", "late", "", NA)
de_meta <- de_meta %>%
  mutate(time_order = match(time_since_onset,
                            time_order,
                            nomatch = length(time_order))) %>%
  arrange(subject_id, time_order) %>%
  distinct(subject_id, .keep_all = TRUE)

cat("DE metadata samples (after deduplication):", nrow(de_meta), "\n")
cat("Cohort breakdown in DE metadata:\n")
print(table(de_meta$cohort))

# Find common samples between expression matrix and DE metadata
common_samples <- intersect(colnames(expr_mat), de_meta$rna_id)
cat("\nSamples in expression matrix:", ncol(expr_mat), "\n")
cat("Samples in DE metadata:", nrow(de_meta), "\n")
cat("Common samples:", length(common_samples), "\n")

if (length(common_samples) < nrow(de_meta)) {
  missing <- setdiff(de_meta$rna_id, colnames(expr_mat))
  cat("Samples in DE metadata but not in nlcpm matrix:", length(missing), "\n")
  cat("Missing:", missing, "\n")
}

# Subset expression matrix and align metadata
expr_mat <- expr_mat[, common_samples]
meta <- de_meta %>%
  filter(rna_id %in% common_samples) %>%
  arrange(match(rna_id, common_samples))

stopifnot(all(colnames(expr_mat) == meta$rna_id))

cat("\nFinal sample set:", nrow(meta), "samples\n")
cat("Cohort breakdown:\n")
print(table(meta$cohort))


# ==============================================================================
# 6. FILTER LOW-VARIANCE GENES
# ==============================================================================
#
# GSVA performance improves when very low-variance genes are removed.
# Keep top 10,000 genes by variance (standard practice).

cat("\n=== Step 6: Filtering to top 10,000 genes by variance ===\n")

gene_var <- apply(expr_mat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:10000]
expr_filtered <- expr_mat[top_genes, ]

cat("Genes retained:", nrow(expr_filtered), "\n")


# ==============================================================================
# 7. SAVE
# ==============================================================================

cat("\n=== Step 7: Saving ===\n")

saveRDS(expr_filtered, "data/processed/expr_matrix.rds")
saveRDS(meta,          "data/processed/metadata.rds")

cat("Saved:\n")
cat("  data/processed/expr_matrix.rds —",
    nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")
cat("  data/processed/metadata.rds    —", nrow(meta), "samples\n")
cat("\nCohort breakdown:\n")
print(table(meta$cohort))
cat("\nNext step: run 01_gsva.R\n")
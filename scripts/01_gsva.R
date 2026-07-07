# ==============================================================================
# 01_gsva.R
# GSVA pathway activity scoring — GSE161731
#
# Part of: sample-level-pathway-scoring
# Purpose: Assign per-sample pathway activity scores using GSVA.
#          Scores reflect relative pathway activity within each sample
#          on a scale of approximately -1 to +1.
#
# Resources: MSigDB Hallmark (50 gene sets), Reactome (~1,500 pathways)
#
# Inputs:
#   data/processed/expr_matrix.rds   — from 00_prepare_expression.R
#   data/processed/metadata.rds      — from 00_prepare_expression.R
#
# Outputs:
#   results/gsva_scores_hallmark.rds  — samples x Hallmark pathway scores
#   results/gsva_scores_reactome.rds  — samples x Reactome pathway scores
#   figures/gsva_heatmap_hallmark.png — heatmap annotated by cohort
#   figures/gsva_heatmap_reactome.png — heatmap top variable Reactome pathways
#   figures/gsva_boxplot_hallmark.png — top discriminating pathways per cohort
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(GSVA)           # BiocManager::install("GSVA")
library(msigdbr)
library(GSEABase)
library(ComplexHeatmap) # BiocManager::install("ComplexHeatmap")
library(circlize)
library(ggpubr)

set.seed(42)
setwd("~/Desktop/github_repos/sample-level-pathway-scoring")


# ==============================================================================
# 2. LOAD INPUTS
# ==============================================================================

cat("=== Loading inputs ===\n")

expr <- readRDS("data/processed/expr_matrix.rds")
meta <- readRDS("data/processed/metadata.rds")

cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Cohort breakdown:\n")
print(table(meta$cohort))

# Define cohort colors for all plots
cohort_colors <- c(
  "COVID-19"  = "#E64B35",
  "Influenza" = "#4DBBD5",
  "Bacterial" = "#F39B7F",
  "CoV other" = "#7E6148",
  "healthy"   = "#52B14A"
)


# ==============================================================================
# 3. LOAD GENE SETS
# ==============================================================================

cat("\n=== Step 1: Loading gene sets ===\n")

# Hallmark
h_sets    <- msigdbr(species = "Homo sapiens", collection = "H")
h_list    <- split(h_sets$gene_symbol, h_sets$gs_name)
names(h_list) <- gsub("HALLMARK_", "", names(h_list))
names(h_list) <- gsub("_", " ", names(h_list))
cat("Hallmark gene sets:", length(h_list), "\n")

# Reactome
r_sets    <- msigdbr(species = "Homo sapiens",
                     collection = "C2", subcollection = "CP:REACTOME")
r_list    <- split(r_sets$gene_symbol, r_sets$gs_name)
names(r_list) <- gsub("REACTOME_", "", names(r_list))
names(r_list) <- gsub("_", " ", names(r_list))
cat("Reactome gene sets:", length(r_list), "\n")


# ==============================================================================
# 4. RUN GSVA
# ==============================================================================

cat("\n=== Step 2: Running GSVA ===\n")

run_gsva <- function(expr_mat, gene_sets, label) {
  cat("Running GSVA:", label, "...\n")
  param <- gsvaParam(
    exprData  = expr_mat,
    geneSets  = gene_sets,
    minSize   = 15,
    maxSize   = 500
  )
  scores <- gsva(param, verbose = FALSE)
  cat("  Done —", nrow(scores), "pathways x", ncol(scores), "samples\n")
  scores
}

gsva_hallmark  <- run_gsva(expr, h_list, "Hallmark")
gsva_reactome  <- run_gsva(expr, r_list, "Reactome")

saveRDS(gsva_hallmark, "results/gsva_scores_hallmark.rds")
saveRDS(gsva_reactome, "results/gsva_scores_reactome.rds")
cat("Saved GSVA score matrices\n")


# ==============================================================================
# 5. HEATMAP — HALLMARK
# ==============================================================================

cat("\n=== Step 3: Hallmark heatmap ===\n")

# Column annotation
col_anno <- HeatmapAnnotation(
  Cohort = meta$cohort,
  Batch  = factor(meta$batch),
  col    = list(
    Cohort = cohort_colors,
    Batch  = c("1" = "#999999", "2" = "#DDDDDD")
  ),
  annotation_name_side = "left"
)

# Order samples by cohort
sample_order <- order(meta$cohort)

png("figures/gsva_heatmap_hallmark.png",
    width = 18, height = 9, units = "in", res = 300)

top25_hallmark <- names(sort(apply(gsva_hallmark, 1, var),
                             decreasing = TRUE))[1:25]

draw(Heatmap(
  gsva_hallmark[top25_hallmark, sample_order],
  name             = "GSVA\nscore",
  col              = colorRamp2(c(-0.6, 0, 0.6),
                                c("#2166AC", "white", "#D6604D")),
  top_annotation   = col_anno[sample_order],
  show_column_names = FALSE,
  rect_gp          = gpar(col = "black", lwd = 0.5),
  row_names_gp     = gpar(fontsize = 8),
  column_split          = meta$cohort[match(colnames(gsva_hallmark[top25_hallmark, sample_order]), meta$rna_id)],
  cluster_column_slices = FALSE,
  cluster_columns       = TRUE,
  column_title     = "GSE161731 — GSVA Hallmark Pathway Scores",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  clustering_distance_rows    = "pearson",
  clustering_distance_columns = "pearson",
  heatmap_legend_param = list(direction = "horizontal"),
  heatmap_width = unit(0.65, "npc")
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()
cat("Saved: figures/gsva_heatmap_hallmark.png\n")


# ==============================================================================
# 6. HEATMAP — REACTOME (TOP 50 VARIABLE PATHWAYS)
# ==============================================================================

cat("\n=== Step 4: Reactome heatmap (top 50 variable pathways) ===\n")

pathway_var <- apply(gsva_reactome, 1, var)
top25_reactome <- names(sort(pathway_var, decreasing = TRUE))[1:25]

# Clean long Reactome names for display
clean_names <- top25_reactome

reactome_subset <- gsva_reactome[top25_reactome, sample_order]
rownames(reactome_subset) <- clean_names
library(stringr)
clean_names_wrapped <- str_wrap(rownames(reactome_subset), width = 33)
rownames(reactome_subset) <- clean_names_wrapped

png("figures/gsva_heatmap_reactome.png",
    width = 22, height = 14, units = "in", res = 300)

draw(Heatmap(
  reactome_subset,
  name             = "GSVA\nscore",
  col              = colorRamp2(c(-0.6, 0, 0.6),
                                c("#2166AC", "white", "#D6604D")),
  top_annotation   = col_anno[sample_order],
  show_column_names = FALSE,
  rect_gp          = gpar(col = "black", lwd = 0.5),
  row_names_gp     = gpar(fontsize = 8),
  height           = unit(35 * 0.8, "cm"),
  column_split          = meta$cohort[match(colnames(reactome_subset), meta$rna_id)],
  cluster_column_slices = FALSE,
  cluster_columns       = TRUE,
  column_title     = "GSVA Reactome Scores",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  clustering_distance_rows    = "pearson",
  clustering_distance_columns = "pearson",
  heatmap_legend_param = list(direction = "horizontal"),
  heatmap_width = unit(0.58, "npc")
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()
cat("Saved: figures/gsva_heatmap_reactome.png\n")


# ==============================================================================
# 7. BOXPLOT — TOP DISCRIMINATING HALLMARK PATHWAYS
# ==============================================================================

cat("\n=== Step 5: Boxplot — top discriminating Hallmark pathways ===\n")

# Select top 9 pathways by variance across samples
top9 <- names(sort(apply(gsva_hallmark, 1, var),
                   decreasing = TRUE))[1:9]

plot_df <- as.data.frame(t(gsva_hallmark[top9, ])) %>%
  rownames_to_column("rna_id") %>%
  left_join(meta %>% dplyr::select(rna_id, cohort), by = "rna_id") %>%
  pivot_longer(-c(rna_id, cohort),
               names_to  = "pathway",
               values_to = "gsva_score") %>%
  mutate(cohort = factor(cohort,
                         levels = c("healthy","COVID-19","Influenza",
                                    "Bacterial","CoV other")))

ggplot(plot_df, aes(x = cohort, y = gsva_score, fill = cohort)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = cohort_colors, guide = "none") +
  facet_wrap(~ pathway, ncol = 3, scales = "free_y") +
  labs(
    title    = "GSVA Hallmark pathway scores by cohort",
    subtitle = "Top 9 pathways by cross-sample variance — GSE161731",
    x        = NULL, y = "GSVA score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
    strip.text   = element_text(size = 8, face = "bold"),
    plot.title   = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9)
  )

ggsave("figures/gsva_boxplot_hallmark.png",
       width = 12, height = 10, dpi = 300)
cat("Saved: figures/gsva_boxplot_hallmark.png\n")

cat("\n=== GSVA complete ===\n")
cat("Next step: run 02_ssgsea.R\n")
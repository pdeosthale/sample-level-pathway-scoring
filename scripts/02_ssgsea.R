# ==============================================================================
# 02_ssgsea.R
# ssGSEA pathway activity scoring — GSE161731
#
# Part of: sample-level-pathway-scoring
# Purpose: Assign per-sample pathway activity scores using ssGSEA.
#          Unlike GSVA, ssGSEA scores are not bounded — they reflect
#          absolute enrichment of a gene set in each sample's ranked list.
#
# Resources: MSigDB Hallmark (50 gene sets), Reactome (~1,500 pathways)
#
# Inputs:
#   data/processed/expr_matrix.rds   — from 00_prepare_expression.R
#   data/processed/metadata.rds      — from 00_prepare_expression.R
#
# Outputs:
#   results/ssgsea_scores_hallmark.rds
#   results/ssgsea_scores_reactome.rds
#   figures/ssgsea_heatmap_hallmark.png
#   figures/ssgsea_heatmap_reactome.png
#   figures/ssgsea_boxplot_hallmark.png
# ==============================================================================

library(tidyverse)
library(GSVA)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)

set.seed(42)
setwd("~/Desktop/github_repos/sample-level-pathway-scoring")


# ==============================================================================
# 1. LOAD INPUTS
# ==============================================================================

cat("=== Loading inputs ===\n")

expr <- readRDS("data/processed/expr_matrix.rds")
meta <- readRDS("data/processed/metadata.rds")

cohort_colors <- c(
  "COVID-19"  = "#E64B35",
  "Influenza" = "#4DBBD5",
  "Bacterial" = "#F39B7F",
  "CoV other" = "#7E6148",
  "healthy"   = "#52B14A"
)

# Load gene sets (same as 01_gsva.R)
h_sets <- msigdbr(species = "Homo sapiens", collection = "H")
h_list <- split(h_sets$gene_symbol, h_sets$gs_name)
names(h_list) <- gsub("HALLMARK_", "", names(h_list))
names(h_list) <- gsub("_", " ", names(h_list))

r_sets <- msigdbr(species = "Homo sapiens",
                  collection = "C2", subcollection = "CP:REACTOME")
r_list <- split(r_sets$gene_symbol, r_sets$gs_name)
names(r_list) <- gsub("REACTOME_", "", names(r_list))
names(r_list) <- gsub("_", " ", names(r_list))

cat("Expression:", nrow(expr), "genes x", ncol(expr), "samples\n")


# ==============================================================================
# 2. RUN ssGSEA
# ==============================================================================

cat("\n=== Step 1: Running ssGSEA ===\n")

run_ssgsea <- function(expr_mat, gene_sets, label) {
  cat("Running ssGSEA:", label, "...\n")
  param <- ssgseaParam(
    exprData = expr_mat,
    geneSets = gene_sets,
    minSize  = 15,
    maxSize  = 500,
    normalize = TRUE   # normalize scores to [-1, 1] for comparability with GSVA
  )
  scores <- gsva(param, verbose = FALSE)
  cat("  Done —", nrow(scores), "pathways x", ncol(scores), "samples\n")
  scores
}

ssgsea_hallmark <- run_ssgsea(expr, h_list, "Hallmark")
ssgsea_reactome <- run_ssgsea(expr, r_list, "Reactome")

saveRDS(ssgsea_hallmark, "results/ssgsea_scores_hallmark.rds")
saveRDS(ssgsea_reactome, "results/ssgsea_scores_reactome.rds")
cat("Saved ssGSEA score matrices\n")


# ==============================================================================
# 3. HEATMAP — HALLMARK
# ==============================================================================

cat("\n=== Step 2: Hallmark heatmap ===\n")

col_anno <- HeatmapAnnotation(
  Cohort = meta$cohort,
  Batch  = factor(meta$batch),
  col    = list(
    Cohort = cohort_colors,
    Batch  = c("1" = "#999999", "2" = "#DDDDDD")
  ),
  annotation_name_side = "left"
)

sample_order <- order(meta$cohort)

png("figures/ssgsea_heatmap_hallmark.png",
    width = 18, height = 9, units = "in", res = 300)

top25_hallmark <- names(sort(apply(ssgsea_hallmark, 1, var),
                             decreasing = TRUE))[1:25]

draw(Heatmap(
  ssgsea_hallmark[top25_hallmark, sample_order],
  name             = "ssGSEA\nscore",
  col              = colorRamp2(c(-1, 0, 1),
                                c("#2166AC", "white", "#D6604D")),
  top_annotation   = col_anno[sample_order],
  show_column_names = FALSE,
  rect_gp          = gpar(col = "black", lwd = 0.5),
  row_names_gp     = gpar(fontsize = 8),
  column_split          = meta$cohort[match(colnames(ssgsea_hallmark[top25_hallmark, sample_order]), meta$rna_id)],
  cluster_column_slices = FALSE,
  cluster_columns       = TRUE,
  column_title     = "ssGSEA Hallmark Pathway Scores",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  clustering_distance_rows    = "pearson",
  clustering_distance_columns = "pearson",
  heatmap_legend_param = list(direction = "horizontal"),
  heatmap_width = unit(0.65, "npc")
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()
cat("Saved: figures/ssgsea_heatmap_hallmark.png\n")


# ==============================================================================
# 4. HEATMAP — REACTOME TOP 50
# ==============================================================================

cat("\n=== Step 3: Reactome heatmap ===\n")

pathway_var   <- apply(ssgsea_reactome, 1, var)
top25_reactome <- names(sort(pathway_var, decreasing = TRUE))[1:25]
clean_names   <- strtrim(top25_reactome, 55)

reactome_subset <- ssgsea_reactome[top25_reactome, sample_order]
rownames(reactome_subset) <- clean_names
library(stringr)
clean_names_wrapped <- str_wrap(rownames(reactome_subset), width = 33)
rownames(reactome_subset) <- clean_names_wrapped

png("figures/ssgsea_heatmap_reactome.png",
    width = 22, height = 14, units = "in", res = 300)

draw(Heatmap(
  reactome_subset,
  name             = "ssGSEA\nscore",
  col              = colorRamp2(c(-1, 0, 1),
                                c("#2166AC", "white", "#D6604D")),
  top_annotation   = col_anno[sample_order],
  show_column_names = FALSE,
  rect_gp          = gpar(col = "black", lwd = 0.5),
  row_names_gp     = gpar(fontsize = 8),
  height           = unit(35 * 0.8, "cm"),
  column_split          = meta$cohort[match(colnames(reactome_subset), meta$rna_id)],
  cluster_column_slices = FALSE,
  cluster_columns       = TRUE,
  column_title     = "ssGSEA Reactome Scores",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  clustering_distance_rows    = "pearson",
  clustering_distance_columns = "pearson",
  heatmap_legend_param = list(direction = "horizontal"),
  heatmap_width = unit(0.58, "npc")
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()
cat("Saved: figures/ssgsea_heatmap_reactome.png\n")


# ==============================================================================
# 5. BOXPLOT — TOP DISCRIMINATING HALLMARK PATHWAYS
# ==============================================================================

cat("\n=== Step 4: Boxplot ===\n")

top9 <- names(sort(apply(ssgsea_hallmark, 1, var),
                   decreasing = TRUE))[1:9]

plot_df <- as.data.frame(t(ssgsea_hallmark[top9, ])) %>%
  rownames_to_column("rna_id") %>%
  left_join(meta %>% dplyr::select(rna_id, cohort), by = "rna_id") %>%
  pivot_longer(-c(rna_id, cohort),
               names_to  = "pathway",
               values_to = "score") %>%
  mutate(cohort = factor(cohort,
                         levels = c("healthy","COVID-19","Influenza",
                                    "Bacterial","CoV other")))

ggplot(plot_df, aes(x = cohort, y = score, fill = cohort)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = cohort_colors, guide = "none") +
  facet_wrap(~ pathway, ncol = 3, scales = "free_y") +
  labs(
    title    = "ssGSEA Hallmark pathway scores by cohort",
    subtitle = "Top 9 pathways by cross-sample variance",
    x        = NULL, y = "ssGSEA score (normalized)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x   = element_text(angle = 35, hjust = 1, size = 8),
    strip.text    = element_text(size = 8, face = "bold"),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9)
  )

ggsave("figures/ssgsea_boxplot_hallmark.png",
       width = 12, height = 10, dpi = 300)
cat("Saved: figures/ssgsea_boxplot_hallmark.png\n")

cat("\n=== ssGSEA complete ===\n")
cat("Next step: run 03_compare_methods.R\n")
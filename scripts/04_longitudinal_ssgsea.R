# ==============================================================================
# 04_longitudinal_ssgsea.R

# ssGSEA longitudinal trajectory analysis — COVID-19 patients
#
# Part of: sample-level-pathway-scoring
# Purpose: Demonstrate ssGSEA's strength for longitudinal single-patient
#          data. Uses COVID-19 samples (early/middle/late timepoints)
#          from 14 subjects excluded from the cross-cohort analysis.
#          Raw metadata contains 45 rows for these subjects; after
#          deduplication to one sample per subject per timepoint (keeping
#          first rna_id alphabetically when visit replicates exist), 31
#          samples are retained for analysis.
#
# Why ssGSEA for longitudinal data:
#   ssGSEA scores are computed independently per sample — adding or removing
#   timepoints does not change existing scores. This makes them stable for
#   longitudinal tracking. GSVA scores depend on all samples in the matrix
#   and shift when the sample set changes, making it less appropriate for
#   temporal trajectory analysis.
#
# Design:
#   14 COVID-19 subjects × up to 3 timepoints (early/middle/late)
#   Unbalanced: not all subjects have all three timepoints
#   Post-deduplication: 31 samples (early: n=7, middle: n=13, late: n=11)
#
# Gene sets: MSigDB Hallmark (50 gene sets)
#
# Inputs:
#   data/raw/GSE161731_xpr_nlcpm.csv.gz   — full nlcpm matrix
#   data/raw/GSE161731_key.csv.gz         — sample metadata
#
# Outputs:
#   data/processed/expr_longitudinal.rds    — expression matrix, 31 samples
#   data/processed/meta_longitudinal.rds    — metadata, 31 samples
#   results/ssgsea_longitudinal_hallmark.rds — ssGSEA scores (50 x 31)
#   figures/longitudinal_trajectories.png   — per-pathway trajectory plots
#                                             (top 9 by temporal variance)
#   figures/longitudinal_heatmap.png        — z-scored heatmap ordered by
#                                             timepoint
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(GSVA)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

set.seed(42)
setwd("~/Desktop/github_repos/sample-level-pathway-scoring")
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results",        showWarnings = FALSE)
dir.create("figures",        showWarnings = FALSE)


# ==============================================================================
#  2. LOAD LONGITUDINAL SAMPLES
# ==============================================================================

cat("=== Step 1: Loading longitudinal COVID-19 samples ===\n")

key <- read.csv("data/raw/GSE161731_key.csv.gz")

# Extract subjects with multiple timepoints
covid_longitudinal <- key %>%
  filter(cohort == "COVID-19") %>%
  group_by(subject_id) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  mutate(
    timepoint = factor(time_since_onset,
                       levels = c("early", "middle", "late"))
  ) %>%
  # Keep one sample per subject per timepoint (first rna_id alphabetically)
  # Multiple samples with the same timepoint label are replicates from
  # the same visit — not additional longitudinal timepoints
  arrange(subject_id, timepoint, rna_id) %>%
  distinct(subject_id, timepoint, .keep_all = TRUE)

cat("Longitudinal samples after deduplication:", nrow(covid_longitudinal), "\n")
cat("Subjects:", n_distinct(covid_longitudinal$subject_id), "\n")
cat("\nTimepoints per subject:\n")
covid_longitudinal %>%
  group_by(subject_id) %>%
  summarise(
    timepoints = paste(as.character(timepoint), collapse = " → "),
    n          = n(),
    .groups    = "drop"
  ) %>%
  print()


# ==============================================================================
# 3. BUILD EXPRESSION MATRIX
# ==============================================================================

cat("\n=== Step 2: Building expression matrix ===\n")

# Load full nlcpm matrix (samples x genes) — already downloaded
expr_raw <- read.csv("data/raw/GSE161731_xpr_nlcpm.csv.gz",
                     row.names = 1, check.names = FALSE)

# Subset to longitudinal samples and transpose to genes x samples
expr_long <- t(as.matrix(expr_raw[covid_longitudinal$rna_id, ]))

cat("Expression matrix:", nrow(expr_long), "genes x",
    ncol(expr_long), "samples\n")
stopifnot(all(colnames(expr_long) == covid_longitudinal$rna_id))

# Filter to top 10,000 genes by variance
gene_var    <- apply(expr_long, 1, var)
top_genes   <- names(sort(gene_var, decreasing = TRUE))[1:10000]
expr_long   <- expr_long[top_genes, ]

cat("After variance filter:", nrow(expr_long), "genes retained\n")

saveRDS(expr_long,          "data/processed/expr_longitudinal.rds")
saveRDS(covid_longitudinal, "data/processed/meta_longitudinal.rds")


# ==============================================================================
# 4. RUN ssGSEA
# ==============================================================================

cat("\n=== Step 3: Running ssGSEA ===\n")

h_sets <- msigdbr(species = "Homo sapiens", collection = "H")
h_list <- split(h_sets$gene_symbol, h_sets$gs_name)
names(h_list) <- gsub("HALLMARK_", "", names(h_list))
names(h_list) <- gsub("_", " ", names(h_list))

param  <- ssgseaParam(
  exprData  = expr_long,
  geneSets  = h_list,
  minSize   = 15,
  maxSize   = 500,
  normalize = TRUE
)

ssgsea_long <- gsva(param, verbose = FALSE)
cat("ssGSEA scores:", nrow(ssgsea_long), "pathways x",
    ncol(ssgsea_long), "samples\n")

saveRDS(ssgsea_long, "results/ssgsea_longitudinal_hallmark.rds")
cat("Saved: results/ssgsea_longitudinal_hallmark.rds\n")


# ==============================================================================
# 5. TRAJECTORY PLOTS — TOP DISCRIMINATING PATHWAYS
# ==============================================================================
#
# For each pathway, plot the score trajectory per subject across timepoints.
# Lines connect the same subject across early → middle → late.
#
# Pathway selection: variance of per-timepoint means rather than
# cross-sample variance. Cross-sample variance captures individual
# heterogeneity but not temporal signal — a pathway with high between-
# subject noise and a flat mean trajectory would rank highly. Variance of
# the timepoint means specifically captures early → middle → late change,
# selecting pathways that actually move over the disease course.

cat("\n=== Step 4: Trajectory plots ===\n")

# Build full tidy data frame across all pathways (needed for selection)
traj_df_all <- as.data.frame(t(ssgsea_long)) %>%
  rownames_to_column("rna_id") %>%
  pivot_longer(-rna_id, names_to = "pathway", values_to = "score") %>%
  left_join(covid_longitudinal %>%
              dplyr::select(rna_id, subject_id, timepoint),
            by = "rna_id") %>%
  mutate(timepoint = factor(timepoint, levels = c("early", "middle", "late")))

# Select top 9 pathways by variance of per-timepoint means
# (captures genuine temporal signal, not just inter-subject noise)
tp_means <- traj_df_all %>%
  group_by(pathway, timepoint) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = timepoint, values_from = mean_score)

pathway_tp_var <- apply(tp_means[, c("early", "middle", "late")], 1, var)
top9           <- tp_means$pathway[order(pathway_tp_var, decreasing = TRUE)][1:9]

cat("Top 9 pathways by temporal variance:\n")
print(top9)

# Subset tidy data frame to top 9
traj_df <- traj_df_all %>% filter(pathway %in% top9)

# (traj_df_all no longer needed)
rm(traj_df_all)

# Summarise: mean + SE per timepoint for background ribbon
summary_df <- traj_df %>%
  group_by(pathway, timepoint) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    se_score   = sd(score, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

ggplot() +
  # Individual subject lines
  geom_line(data = traj_df,
            aes(x = timepoint, y = score,
                group = subject_id),
            color = "grey70", linewidth = 0.4, alpha = 0.7) +
  geom_point(data = traj_df,
             aes(x = timepoint, y = score),
             color = "grey60", size = 1, alpha = 0.7) +
  # Mean trajectory with ribbon
  geom_ribbon(data = summary_df,
              aes(x = timepoint,
                  ymin = mean_score - se_score,
                  ymax = mean_score + se_score,
                  group = 1),
              fill = "#E64B35", alpha = 0.2) +
  geom_line(data = summary_df,
            aes(x = timepoint, y = mean_score, group = 1),
            color = "#E64B35", linewidth = 1.2) +
  geom_point(data = summary_df,
             aes(x = timepoint, y = mean_score),
             color = "#E64B35", size = 3) +
  facet_wrap(~ pathway, ncol = 3, scales = "free_y",
             labeller = labeller(pathway = label_wrap_gen(25))) +
  labs(
    title    = "COVID-19 longitudinal pathway trajectories (ssGSEA)",
    subtitle = "Grey lines = individual subjects | Red = mean ± SE | n=14 subjects",
    x        = "Timepoint", y = "ssGSEA score (normalized)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text    = element_text(size = 8, face = "bold"),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9),
    axis.text.x   = element_text(size = 9)
  )

ggsave("figures/longitudinal_trajectories.png",
       width = 12, height = 10, dpi = 300)
cat("Saved: figures/longitudinal_trajectories.png\n")


# ==============================================================================
# 6. HEATMAP — SAMPLES ORDERED BY TIMEPOINT
# ==============================================================================

cat("\n=== Step 5: Heatmap ordered by timepoint ===\n")

# Row-scale ssGSEA scores to z-scores (scale per pathway across samples).
# Raw ssGSEA scores cluster between ~0.3–0.7 for most Hallmark pathways,
# compressing the color range and washing out cross-timepoint patterns.
# Z-scoring each pathway independently highlights relative changes across
# the disease course while preserving row clustering structure.
ssgsea_z <- t(scale(t(ssgsea_long)))

# Subset to top 25 pathways by variance of z-scores across samples.
# Reduces visual clutter while retaining the most variable pathways.
top25_rows <- names(sort(apply(ssgsea_z, 1, var), decreasing = TRUE))[1:25]
ssgsea_z   <- ssgsea_z[top25_rows, ]

# Order samples: early → middle → late
timepoint_order <- order(covid_longitudinal$timepoint)

col_anno <- HeatmapAnnotation(
  Timepoint = covid_longitudinal$timepoint[timepoint_order],
  Subject   = covid_longitudinal$subject_id[timepoint_order],
  col = list(
    Timepoint = c("early"  = "#FED976",
                  "middle" = "#FC4E2A",
                  "late"   = "#800026"),
    Subject   = setNames(
      colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(
        n_distinct(covid_longitudinal$subject_id)),
      unique(covid_longitudinal$subject_id[timepoint_order])
    )
  ),
  annotation_name_side = "left"
)

png("figures/longitudinal_heatmap.png",
    width = 14, height = 10, units = "in", res = 300)

draw(Heatmap(
  ssgsea_z[, timepoint_order],
  name             = "ssGSEA\n(z-score)",
  col              = colorRamp2(c(-2, 0, 2),
                                c("#2166AC", "white", "#D6604D")),
  top_annotation   = col_anno,
  show_column_names = FALSE,
  rect_gp          = gpar(col = NA),              # no cell borders (31 cols × 50 rows)
  row_names_gp     = gpar(fontsize = 8),
  column_title     = "COVID-19 Longitudinal ssGSEA Hallmark Scores (z-scored per pathway)",
  column_title_gp  = gpar(fontsize = 13, fontface = "bold"),
  cluster_columns  = FALSE,                       # preserve timepoint order
  clustering_distance_rows = "pearson",
  heatmap_legend_param = list(direction = "horizontal"),
  heatmap_width    = unit(0.75, "npc")
), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()
cat("Saved: figures/longitudinal_heatmap.png\n")

cat("\n=== Longitudinal ssGSEA complete ===\n")
cat("Samples analysed: 31 (14 subjects × up to 3 timepoints, post-dedup)\n")
cat("figures/longitudinal_trajectories.png — top 9 pathways by temporal variance\n")
cat("figures/longitudinal_heatmap.png      — all 50 Hallmark pathways, z-scored\n")
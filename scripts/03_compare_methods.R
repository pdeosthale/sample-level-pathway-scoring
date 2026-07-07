# ==============================================================================
# 03_compare_methods.R
# Compare GSVA vs ssGSEA scores — GSE161731
#
# Part of: sample-level-pathway-scoring
# Purpose: Evaluate concordance between GSVA and ssGSEA and assess
#          whether either method separates infection cohorts in
#          an unsupervised setting.
#
# Comparisons:
#   1. Score correlation: are GSVA and ssGSEA scores concordant per pathway?
#   2. Cohort separation: does unsupervised clustering recover infection groups?
#   3. Top discriminating pathways: which pathways best distinguish cohorts?
#   4. Score distribution: GSVA vs ssGSEA range and spread per pathway
#
# Inputs:
#   results/gsva_scores_hallmark.rds
#   results/ssgsea_scores_hallmark.rds
#   data/processed/metadata.rds
#
# Outputs:
#   figures/compare_top_pathways.png          — top discriminating pathways
# ==============================================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

set.seed(42)
setwd("~/Desktop/github_repos/sample-level-pathway-scoring")


# ==============================================================================
gsva_h   <- readRDS("results/gsva_scores_hallmark.rds")
ssgsea_h <- readRDS("results/ssgsea_scores_hallmark.rds")
meta     <- readRDS("data/processed/metadata.rds")

cohort_colors <- c(
  "COVID-19"  = "#E64B35",
  "Influenza" = "#4DBBD5",
  "Bacterial" = "#F39B7F",
  "CoV other" = "#7E6148",
  "healthy"   = "#52B14A"
)

# Align: both matrices should have same samples and pathways
common_samples  <- intersect(colnames(gsva_h), colnames(ssgsea_h))
common_pathways <- intersect(rownames(gsva_h), rownames(ssgsea_h))

gsva_h   <- gsva_h[common_pathways, common_samples]
ssgsea_h <- ssgsea_h[common_pathways, common_samples]
meta     <- meta %>% filter(rna_id %in% common_samples) %>%
  arrange(match(rna_id, common_samples))

cat("Aligned:", nrow(gsva_h), "pathways x", ncol(gsva_h), "samples\n")


# ==============================================================================
# 2. TOP DISCRIMINATING PATHWAYS — KRUSKAL-WALLIS
# ==============================================================================

cat("\n=== Step 2: Top discriminating pathways ===\n")

get_top_discriminating <- function(score_mat, meta, n = 15) {
  map_dfr(rownames(score_mat), function(pw) {
    scores <- score_mat[pw, ]
    groups <- meta$cohort[match(names(scores), meta$rna_id)]
    kw     <- kruskal.test(scores ~ groups)
    tibble(pathway = pw, kw_pvalue = kw$p.value,
           kw_statistic = kw$statistic)
  }) %>%
    mutate(kw_fdr = p.adjust(kw_pvalue, method = "BH")) %>%
    arrange(kw_fdr) %>%
    head(n)
}

top_gsva   <- get_top_discriminating(gsva_h,   meta)
top_ssgsea <- get_top_discriminating(ssgsea_h, meta)

cat("Top discriminating pathways — GSVA (Kruskal-Wallis FDR < 0.05):\n")
print(top_gsva %>% dplyr::select(pathway, kw_fdr) %>%
        mutate(kw_fdr = signif(kw_fdr, 3)) %>% as.data.frame())

cat("\nTop discriminating pathways — ssGSEA:\n")
print(top_ssgsea %>% dplyr::select(pathway, kw_fdr) %>%
        mutate(kw_fdr = signif(kw_fdr, 3)) %>% as.data.frame())

cat("\nPathways significant in both methods:\n")
print(intersect(top_gsva$pathway, top_ssgsea$pathway))

# Boxplot of top 6 discriminating pathways, both methods side by side
top6 <- union(head(top_gsva$pathway, 6), head(top_ssgsea$pathway, 6)) %>%
  head(6)

plot_df <- bind_rows(
  as.data.frame(t(gsva_h[top6, ])) %>%
    rownames_to_column("rna_id") %>%
    pivot_longer(-rna_id, names_to = "pathway", values_to = "score") %>%
    mutate(method = "GSVA"),
  as.data.frame(t(ssgsea_h[top6, ])) %>%
    rownames_to_column("rna_id") %>%
    pivot_longer(-rna_id, names_to = "pathway", values_to = "score") %>%
    mutate(method = "ssGSEA")
) %>%
  left_join(meta %>% dplyr::select(rna_id, cohort), by = "rna_id") %>%
  mutate(cohort = factor(cohort,
                         levels = c("healthy","COVID-19","Influenza",
                                    "Bacterial","CoV other")))

ggplot(plot_df, aes(x = cohort, y = score, fill = cohort)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = cohort_colors, guide = "none") +
  facet_grid(method ~ pathway, scales = "free_y",
             labeller = labeller(pathway = label_wrap_gen(20))) +
  labs(
    title    = "Top discriminating Hallmark pathways: GSVA vs ssGSEA",
    subtitle = "Selected by Kruskal-Wallis test across cohorts",
    x        = NULL, y = "Score"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 7),
    strip.text   = element_text(size = 7, face = "bold"),
    plot.title   = element_text(size = 11, face = "bold")
  )

ggsave("figures/compare_top_pathways.png",
       width = 14, height = 8, dpi = 300)
cat("Saved: figures/compare_top_pathways.png\n")

cat("\n=== Method comparison complete ===\n")
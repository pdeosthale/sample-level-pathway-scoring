# Sample-Level Pathway Scoring

Demonstration of GSVA and ssGSEA for per-sample pathway activity scoring using blood transcriptomics from respiratory infection patients (GSE161731, McClain et al. 2021, *Nat Commun*).

---

## Biological Question

GSVA and ssGSEA are both sample-level pathway scoring methods but answer fundamentally different questions:

- **GSVA** asks: *which pathways are relatively more active in this sample compared to other samples?*
- **ssGSEA** asks: *which pathways are absolutely enriched within this individual sample, independent of any reference group?*

This repository demonstrates both methods on the same dataset using study designs matched to each method's strengths — cross-cohort comparison for GSVA, longitudinal patient trajectories for ssGSEA.

---

## When to Use Each Method

**Use GSVA when** your primary goal is to compare pathway activity across
biological groups. GSVA estimates relative pathway enrichment by considering
variation across all samples in the matrix, which often provides greater
dynamic range and improved sensitivity for detecting between-group differences.

**Use ssGSEA when** you need pathway scores that are computed independently
for each sample — for example, in single-patient analyses, longitudinal
studies, or when ranking pathways within individual samples. Each sample's
enrichment is computed independently of the others. With `normalize = FALSE`,
adding or removing samples does not alter existing scores; with the default
`normalize = TRUE`, scores are rescaled by the min–max range across the
dataset, so their absolute values depend on which samples are included.

**For between-group comparisons,** ssGSEA can still be used with standard
statistical tests (Wilcoxon, limma, etc.), but its scores may exhibit a
narrower dynamic range than GSVA in some datasets and implementations,
potentially reducing statistical power and making visual separation of groups
less pronounced. This repository demonstrates both use cases on the same
dataset — see `compare_top_pathways.png` for a direct comparison of score
distributions on GSE161731.

| Scenario | Method |
|---|---|
| Which pathways differ between infection groups? | GSVA |
| Unsupervised clustering / heatmap across many samples | GSVA |
| How does a patient's pathway activity change over time? | ssGSEA |
| Which pathways are most active in this one sample? | ssGSEA |
| Adding new samples without recalculating existing scores | ssGSEA with `normalize = FALSE` |

---

## Dataset

| Accession | Description | Design |
|---|---|---|
| [GSE161731](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161731) | McClain et al. 2021, *Nat Commun* | Blood transcriptomics: COVID-19, CoV other, Bacterial, Influenza, Healthy |

**Analysis 1 — Cross-cohort (GSVA):** 152 samples, one per subject (earliest timepoint). COVID-19 n=46, CoV other n=48, Bacterial n=23, Influenza n=17, Healthy n=18.

**Analysis 2 — Longitudinal COVID-19 (ssGSEA):** 31 COVID-19 samples with multiple timepoints (early/middle/late) from 14 subjects. These samples were excluded from Analysis 1 to avoid pseudoreplication — here they are the primary data.

---

## Repository Structure

```
sample-level-pathway-scoring/
│
├── scripts/
│   ├── 00_prepare_expression.R        # Download nlcpm, align to DE metadata
│   ├── 01_gsva.R                      # GSVA: cross-cohort analysis (152 samples)
│   ├── 02_ssgsea.R                    # ssGSEA: cross-cohort (for comparison)
│   ├── 03_compare_methods.R           # Top discriminating pathways, GSVA vs ssGSEA
│   └── 04_longitudinal_ssgsea.R       # ssGSEA: COVID-19 longitudinal trajectories
│
├── data/
│   ├── raw/                           # Downloaded files (not tracked)
│   └── processed/                     # Expression matrices, metadata (not tracked)
│
├── results/                           # Score matrices (not tracked)
│
└── figures/
    ├── gsva_heatmap_hallmark.png
    ├── gsva_boxplot_hallmark.png
    ├── ssgsea_heatmap_hallmark.png
    ├── ssgsea_boxplot_hallmark.png
    ├── compare_top_pathways.png
    ├── longitudinal_trajectories.png
    └── longitudinal_heatmap.png
```

---

## Methods

### Expression Data

Authors' pre-computed normalized log-CPM values (`GSE161731_xpr_nlcpm.csv.gz`) downloaded from GEO. Transposed to genes × samples. Filtered to top 10,000 genes by variance.

### Gene-Set Resources

MSigDB Hallmark (50 gene sets) and Reactome (~1,839 pathways, filtered to 15–500 genes = 625 retained). Loaded via `msigdbr`.

### GSVA (Analysis 1)

`GSVA::gsvaParam()` with `minSize = 15`, `maxSize = 500`. Run on 152 cross-cohort samples. Scores reflect relative pathway activity across the sample matrix. In this dataset, scores ranged approximately −0.5 to +0.5 with most pathways distributed around zero, though the exact range is dataset- and parameter-dependent.

### ssGSEA (Analysis 2)

`GSVA::ssgseaParam()` with `normalize = TRUE`. Run on 31 COVID-19 longitudinal samples. Each sample's enrichment is computed independently, but with `normalize = TRUE` the scores are rescaled by the dataset-wide min–max range, so adding new timepoints can shift existing scores (raw, un-normalized scores would not change). Scores are therefore comparable within this 31-sample run, not across runs with different sample sets.

### Discriminating Pathways

Kruskal-Wallis test across 5 cohorts per pathway, BH-corrected FDR. Top 6 pathways shown for both methods side by side to illustrate score distribution differences.

---

## Key Results

### Analysis 1: Cross-cohort GSVA

GSVA Hallmark scores clearly separate infection cohorts in unsupervised clustering. Bacterial infection shows elevated angiogenesis, IL-6/JAK-STAT3, and TNFα/NF-κB signaling, and suppressed MYC targets, oxidative phosphorylation, and G2M checkpoint. Influenza shows the strongest interferon gamma and interferon alpha response signals. COVID-19 and CoV other cluster together with intermediate interferon scores.

Running the same GSVA analysis with Reactome (625 pathways after size filtering) illustrates an important property of gene set choice: Reactome's finer granularity surfaces pathogen-specific signatures that Hallmark's curated gene sets miss entirely. The top variable Reactome pathways include SARS-CoV-2- and influenza-specific translation machinery pathways (e.g. "SARS-CoV-2 modulates host translation machinery", "Influenza viral RNA transcription and replication"), which cluster cleanly by pathogen in unsupervised analysis — only possible because Reactome contains pathogen-specific gene sets that Hallmark does not. B cell receptor signaling and complement activation pathways also rank highly, consistent with the differential expression findings in `infectious-disease-transcriptomics`. Reactome heatmaps are not included as figures in this repository — the Hallmark figures carry the GSVA vs ssGSEA comparison more cleanly — but the scripts produce them and the observation stands: when pathogen identity rather than host response category is the question, a more granular gene set database is the appropriate choice.

### Cross-cohort ssGSEA vs GSVA

Running ssGSEA on the same 152-sample cross-cohort dataset illustrates a practical limitation for this use case. In this dataset, ssGSEA scores are all positive and span a narrower range than GSVA — the same biological differences are present but with smaller absolute differences between cohorts, reducing statistical power and visual interpretability.

The figure below shows the top discriminating Hallmark pathways for both methods side by side. GSVA (top row) shows cohorts clearly separated for each pathway. ssGSEA (bottom row) shows the same ordering of cohorts but compressed into a narrower range — illustrating that ssGSEA may be less sensitive than GSVA for cross-cohort comparison in this dataset.

### Analysis 2: Longitudinal COVID-19 ssGSEA

ssGSEA was applied to 31 COVID-19 samples from 14 subjects with multiple timepoints (early/middle/late), excluded from Analysis 1 to avoid pseudoreplication. Pathway selection for trajectory plots used variance of per-timepoint means rather than cross-sample variance, to select pathways with genuine temporal signal rather than high inter-subject noise.

**Innate immune contraction:** Interferon alpha and interferon gamma response both show the strongest early → late decline of any pathways in the panel (IFN alpha: ~0.82 → 0.50; IFN gamma: ~0.72 → 0.63). The parallel trajectories of both interferon arms — with tight SE ribbons and consistent subject-level directionality — suggest coordinated contraction of the acute antiviral response over the disease course rather than noise.

**Proliferative and metabolic suppression:** E2F targets, G2M checkpoint, MTORC1 signaling, and MYC targets V1 all decline monotonically from early to late. This pattern is consistent with progressive lymphopenia and T cell exhaustion reported in severe COVID-19, and with metabolic reprogramming away from anabolic growth programs during sustained infection.

**Erythroid/myeloid expansion:** Heme metabolism rises steadily from early to late, potentially reflecting erythropoietic recovery or expanding neutrophil/monocyte populations as infection resolves.

**Note:** Hedgehog signaling rises late but shows wide inter-subject variance. Hallmark HEDGEHOG is a small gene set and this trajectory should be interpreted cautiously.

The heatmap (z-scored per pathway) recapitulates the same structure across all 50 Hallmark pathways. A coherent early-high, late-low cluster encompasses the interferon, proliferative (E2F, G2M, MYC, MTORC1), and metabolic (oxidative phosphorylation, cholesterol homeostasis) pathways. A second cluster of pathways (inflammatory response, complement, coagulation, IL6-JAK-STAT3) shows more heterogeneous subject-level patterns without a clear directional trajectory, reflecting the known inter-individual variability in innate immune activation in COVID-19.

---

## Tools and Languages

**Language:** R

**Sample-level scoring:** GSVA v2.4.9

**Gene sets:** msigdbr (Hallmark, Reactome)

**Visualization:** ComplexHeatmap, ggplot2, tidyverse

---

## Relationship to Other Repositories

- `infectious-disease-transcriptomics` — group-level DE and pathway enrichment on the same GSE161731 dataset
- `Pathway-enrichment-comparison` — ORA/GSEA/CAMERA benchmarking on GSE193336 macrophage LPS data
- `sample-level-pathway-scoring` (this repo) — per-sample pathway activity scoring, demonstrating GSVA for cross-cohort and ssGSEA for longitudinal use cases

---

## Reference

McClain MT, Constantine FJ, Henao R, et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nat Commun.* 2021;12:1079. https://doi.org/10.1038/s41467-021-21289-y

## License

MIT License — see [LICENSE](LICENSE) for details.

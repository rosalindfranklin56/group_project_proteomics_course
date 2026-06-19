# group_project_proteomics_course
Part of "Comparative Proteomics of Human Biofluids: Tears, Cerebrospinal Fluid and Plasma". DIA-MS comparative analysis of tears vs depleted plasma proteomes. Evaluated protein overlap, concordance, abundance distributions, and completeness, while identifying unique and enriched proteins.

# Project 3 — Comparative Proteomic Analysis of Tears and Depleted Plasma

Comparative proteomics project evaluating protein overlap, abundance concordance, and proteome completeness between **tear fluid** and **depleted plasma** using DIA mass spectrometry datasets and R.

---

## Overview

This project compares two peripheral human biofluids — **tears** and **depleted plasma** — at the protein level using DIA proteomics data.

The analysis focuses on:

* protein identification overlap
* shared and biofluid-specific proteins
* abundance concordance of shared proteins
* abundance distribution patterns
* protein completeness across replicates
* visualisation of highly differential shared proteins

The goal is to determine whether tears primarily reflect systemic plasma proteins or contain a substantial tear-enriched proteome.

---

## Repository contents

```text
Project3_Tears_vs_Plasma.R      # Full analysis script
README.md                       # This file
```

All figures and tables are generated automatically after running the script and are stored in the output directory:

```text
Project3_Tears_vs_Plasma_results/
```

---

## Data

Two input files are required.

Place both files in the working directory before running the script.

| File                               | Biofluid        | Format                                   |
| ---------------------------------- | --------------- | ---------------------------------------- |
| Report_training_data.pg_matrix.tsv | Tears           | Wide protein × sample matrix             |
| Depleted_plasma_quant.tsv          | Depleted plasma | Long-format protein quantification table |

---

## Dataset description

### Tears dataset

The tears dataset originates from DIA-NN protein-group output.

Only pooled technical replicates are used:

```text
Pool_DIA_R1
Pool_DIA_R2
Pool_DIA_R3
Pool_DIA_R4
Pool_DIA_R5
Pool_DIA_R6
Pool_DIA_R7
Pool_DIA_R8
```

These represent repeated measurements of the same pooled tear sample.

Consequently, variation within tears primarily reflects technical rather than biological variability.

---

### Plasma dataset

The depleted plasma dataset contains protein quantification values in long format.

Protein quantities are aggregated at the protein level and converted into a protein × sample matrix.

---

## Analysis structure

The script is divided into six sequential sections.

---

### Part A — Data loading and preprocessing

Both datasets are imported into R.

Protein identifiers are standardised:

* UniProt accession extraction
* removal of isoform suffixes
* removal of duplicated protein-group representations

Tear and plasma datasets are converted into comparable protein-level matrices.

---

### Part B — Protein overlap analysis

Protein sets are generated for both biofluids.

The following categories are identified:

* proteins detected in Tears
* proteins detected in Plasma
* proteins shared between both biofluids
* proteins unique to Tears
* proteins unique to Plasma

All protein lists are exported as CSV tables.

---

### Part C — Venn diagram

A custom geometry-based Venn diagram is generated to visualise:

* proteins unique to Tears
* proteins unique to Plasma
* proteins shared by both biofluids

The diagram provides a direct overview of proteome overlap between the two fluids.

---

### Part D — Abundance comparison of shared proteins

For every shared protein:

* mean abundance is calculated in Tears
* mean abundance is calculated in Plasma

A concordance plot is then generated.

Pearson correlation is calculated to quantify the agreement between protein abundances across biofluids.

Proteins with the strongest abundance differences are labelled.

Proteins are classified as:

* tear-enriched
* plasma-enriched
* shared/balanced

using log2 fold-change thresholds.

---

### Part E — Mean abundance distribution

Protein abundance distributions are compared between Tears and Plasma.

Kernel density curves are generated using:

```r
log2(mean abundance)
```

to visualise overall proteome composition.

This analysis highlights whether one biofluid is dominated by highly abundant proteins relative to the other.

---

### Part F — Heatmap of top shared proteins

Shared proteins are ranked by:

```r
|log2 fold change|
```

The proteins showing the largest abundance differences between Tears and Plasma are selected.

A heatmap is generated using:

* row-wise Z-score normalisation
* fixed sample ordering
* biofluid annotation colours

Missing values are imputed using row means solely for visualisation purposes.

No imputation is used for statistical summaries.

---

### Part G — Protein completeness analysis

Protein completeness is defined as:

```text
Number of replicates in which a protein is detected
```

Separate completeness distributions are generated for:

* Tears
* Depleted Plasma

This analysis evaluates detection consistency across samples and provides a measure of proteome robustness.

---

## Biological questions addressed

The project investigates several biological questions:

### 1. How similar are Tears and Plasma?

Similarity is assessed through:

* overlap analysis
* abundance concordance
* shared protein composition

---

### 2. Which proteins are tear-enriched?

Proteins showing substantially higher abundance in tears may represent:

* ocular surface proteins
* lacrimal gland products
* tear-specific defence proteins

---

### 3. Which proteins are plasma-enriched?

Proteins showing substantially higher abundance in plasma likely represent:

* systemic proteins
* transport proteins
* coagulation-related proteins
* complement proteins

---

### 4. Do tears largely reflect systemic background?

The abundance concordance analysis helps determine whether tear composition is primarily driven by systemic proteins or whether tears exhibit a distinct proteomic signature.

---

## Outputs

All outputs are written to:

```text
Project3_Tears_vs_Plasma_results/
```

| File                                    | Description                                   |
| --------------------------------------- | --------------------------------------------- |
| overlap_summary.csv                     | Summary counts of overlap and unique proteins |
| shared_proteins.csv                     | Shared protein list                           |
| tear_unique_proteins.csv                | Tears-only proteins                           |
| plasma_unique_proteins.csv              | Plasma-only proteins                          |
| shared_proteins_summary.csv             | Shared proteins with abundance statistics     |
| shared_proteins_classified.csv          | Shared proteins classified by enrichment      |
| Figure1_Venn_Tears_vs_Plasma.pdf        | Protein overlap Venn diagram                  |
| Figure2_Concordance_SharedProteins.pdf  | Shared protein concordance plot               |
| Figure3_Mean_Abundance_Distribution.pdf | Abundance distribution comparison             |
| Figure4_TopSharedProteins_Heatmap.pdf   | Heatmap of top differential shared proteins   |
| Figure5A_Completeness_Tears.pdf         | Tears completeness distribution               |
| Figure5B_Completeness_Plasma.pdf        | Plasma completeness distribution              |

---

## Dependencies

### CRAN packages

```r
install.packages(
  c(
    "tidyverse",
    "patchwork",
    "scales",
    "ggrepel",
    "pheatmap"
  )
)
```

### Additional package

```r
install.packages("RColorBrewer")
```

used for heatmap colour palettes.

---

## How to run

1. Download the datasets.
2. Place both input files in the project folder.
3. Open the R script in RStudio.
4. Set the working directory:

```r
setwd("C:/proteomics/group_project")
```

5. Run the entire script.

All figures and tables will be generated automatically.

---

## Interpretation notes

Tears and depleted plasma represent biologically distinct fluids but share a substantial fraction of proteins.

Shared proteins are expected to reflect:

* systemic circulation
* inflammatory pathways
* transport proteins
* complement proteins

Proteins enriched in tears may indicate:

* ocular surface biology
* lacrimal gland secretion
* local immune defence mechanisms

Therefore, overlap alone should not be interpreted as biological equivalence between the two fluids.

The abundance comparison and enrichment patterns provide more informative evidence regarding biofluid specificity.

---

## Reproducibility

To save the R environment information:

```r
sink("session_info.txt")
sessionInfo()
sink()
```

This records the exact R version and package versions used for the analysis.

### Session info used for this analysis

```text
R version 4.5.3 (2026-03-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United Kingdom.utf8 
[2] LC_CTYPE=English_United Kingdom.utf8   
[3] LC_MONETARY=English_United Kingdom.utf8
[4] LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.utf8    

time zone: Europe/Kiev
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pheatmap_1.0.13 ggrepel_0.9.8   scales_1.4.0    patchwork_1.3.2
 [5] lubridate_1.9.5 forcats_1.0.1   stringr_1.6.0   dplyr_1.2.1    
 [9] purrr_1.2.2     readr_2.2.0     tidyr_1.3.2     tibble_3.3.1   
[13] ggplot2_4.0.3   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] Matrix_1.7-4       bit_4.6.0          gtable_0.3.6       crayon_1.5.3      
 [5] compiler_4.5.3     tidyselect_1.2.1   Rcpp_1.1.1-1.1     parallel_4.5.3    
 [9] splines_4.5.3      textshaping_1.0.5  systemfonts_1.3.2  lattice_0.22-9    
[13] R6_2.6.1           labeling_0.4.3     generics_0.1.4     pillar_1.11.1     
[17] RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.2.0        stringi_1.8.7     
[21] S7_0.2.2           bit64_4.8.0        timechange_0.4.0   cli_3.6.6         
[25] mgcv_1.9-4         withr_3.0.2        magrittr_2.0.5     grid_4.5.3        
[29] vroom_1.7.1        rstudioapi_0.18.0  hms_1.1.4          nlme_3.1-168      
[33] lifecycle_1.0.5    vctrs_0.7.3        glue_1.8.1         farver_2.1.2      
[37] ragg_1.5.2         tools_4.5.3        pkgconfig_2.0.3   
```

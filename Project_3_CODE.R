# =============================================================================
# PROJECT 3 — Tears_DIA vs Depleted_plasma_DIA:
#   - protein-level overlap
#   - Venn diagram
#   - concordance plot
#   - mean abundance distribution
#   - heatmap of top shared proteins
#   - completeness plot
# =============================================================================

setwd("C:/proteomics/group_project")
# 0) Packages -----------------------------------------------------------------
packages <- c("tidyverse", "patchwork", "scales", "ggrepel", "pheatmap")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(patchwork)
library(scales)
library(ggrepel)
library(pheatmap)

# 1) Working directory and output folder --------------------------------------
setwd("C:/proteomics/group_project")

out_dir <- file.path(
  "C:/proteomics/group_project",
  "Project3_Tears_vs_Plasma_results"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Clear any open graphics devices
while (dev.cur() != 1) dev.off()

# 2) Colours and theme ---------------------------------------------------------
col_tears  <- "#2F6DB3"
col_plasma <- "#E76F6A"
col_dark   <- "#1F2937"

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text          = element_text(color = col_dark),
      axis.text     = element_text(color = col_dark),
      axis.title    = element_text(face = "bold"),
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.85)),
      legend.title  = element_text(face = "bold"),
      strip.text    = element_text(face = "bold"),
      plot.margin   = margin(8, 8, 8, 8)
    )
}

# 3) Helper functions ----------------------------------------------------------

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

clean_protein_id <- function(x) {
  x <- as.character(x)
  
  # If the ID is in a "sp|P02788|TRFL_HUMAN" format, keep only the accession.
  pipe_id <- stringr::str_match(x, "\\|([^|]+)\\|")[, 2]
  x <- ifelse(!is.na(pipe_id), pipe_id, x)
  
  # If several IDs are grouped, keep only the first one.
  x <- stringr::str_split_fixed(x, "[;,]", 2)[, 1]
  
  # Remove isoform suffixes like P01876-1
  x <- stringr::str_remove(x, "-\\d+$")
  x <- stringr::str_trim(x)
  
  x
}

detect_sample_cols <- function(df, pattern) {
  cols <- names(df)[stringr::str_detect(names(df), pattern)]
  if (length(cols) == 0) {
    stop("Could not find sample columns using pattern: ", pattern)
  }
  cols
}

impute_and_clean <- function(mat) {
  # Row-mean imputation is only for visualisation.
  row_means <- rowMeans(mat, na.rm = TRUE)
  
  for (i in seq_len(nrow(mat))) {
    na_idx <- is.na(mat[i, ]) | is.nan(mat[i, ]) | is.infinite(mat[i, ])
    if (any(na_idx)) mat[i, na_idx] <- row_means[i]
  }
  
  # Drop rows with zero variance after imputation.
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  bad_rows <- is.na(row_sd) | row_sd == 0
  
  if (any(bad_rows)) {
    mat <- mat[!bad_rows, , drop = FALSE]
  }
  
  mat
}

summarise_biofluid <- function(wide_df, sample_cols) {
  mat <- wide_df %>% select(all_of(sample_cols)) %>% as.matrix()
  det <- !is.na(mat) & mat > 0
  
  tibble(
    protein = wide_df$protein,
    n_present = rowSums(det),
    mean_abundance = ifelse(
      rowSums(det) > 0,
      rowMeans(ifelse(det, mat, NA_real_), na.rm = TRUE),
      NA_real_
    )
  )
}

# 4) Read data ----------------------------------------------------------------
# Tears data: use only the pooled technical replicates Pool_DIA_R1-R8
# This matches the design described in Project 5 and keeps the comparison clean.
tears_raw <- readr::read_tsv(
  file.path("Tears_DIA", "Report_training_data.pg_matrix.tsv"),
  show_col_types = FALSE
)

# Plasma data: long format, one row per protein/precursor/sample record.
plasma_raw <- readr::read_tsv(
  "Depleted_plasma_quant.tsv",
  show_col_types = FALSE
)

# 5) Build Tears protein matrix ------------------------------------------------
tear_sample_cols <- detect_sample_cols(tears_raw, "Pool_DIA_R[1-8]")
tear_sample_names <- stringr::str_extract(tear_sample_cols, "Pool_DIA_R[1-8]")
names(tear_sample_names) <- tear_sample_cols

tears_matrix <- tears_raw %>%
  select(protein_raw = Protein.Group, all_of(tear_sample_cols)) %>%
  mutate(protein = clean_protein_id(protein_raw)) %>%
  select(-protein_raw) %>%
  mutate(across(all_of(tear_sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  group_by(protein) %>%
  summarise(across(all_of(tear_sample_cols), safe_mean), .groups = "drop")

# Rename the long file-path sample columns to short names: Pool_DIA_R1 ... R8
names(tears_matrix)[-1] <- tear_sample_names

tear_sample_cols <- names(tears_matrix)[-1]

# 6) Build Plasma protein matrix ----------------------------------------------
plasma_sample_order <- plasma_raw %>%
  distinct(R.FileName) %>%
  pull(R.FileName)

plasma_matrix <- plasma_raw %>%
  mutate(
    protein  = clean_protein_id(PG.ProteinAccessions),
    sample   = R.FileName,
    quantity = suppressWarnings(as.numeric(PG.Quantity))
  ) %>%
  filter(!is.na(protein), protein != "", !is.na(sample), sample != "") %>%
  group_by(protein, sample) %>%
  summarise(quantity = safe_mean(quantity), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = quantity)

# Make sure sample columns are in the original sample order
missing_plasma_cols <- setdiff(plasma_sample_order, names(plasma_matrix))
for (s in missing_plasma_cols) plasma_matrix[[s]] <- NA_real_
plasma_matrix <- plasma_matrix %>% select(protein, all_of(plasma_sample_order))

plasma_sample_cols <- names(plasma_matrix)[-1]

# 7) Create protein sets -------------------------------------------------------
tear_set <- tears_matrix %>%
  filter(if_any(all_of(tear_sample_cols), ~ !is.na(.x) & .x > 0)) %>%
  pull(protein)

plasma_set <- plasma_matrix %>%
  filter(if_any(all_of(plasma_sample_cols), ~ !is.na(.x) & .x > 0)) %>%
  pull(protein)

shared_proteins <- intersect(tear_set, plasma_set)
tear_only       <- setdiff(tear_set, plasma_set)
plasma_only     <- setdiff(plasma_set, tear_set)

# 8) Save summary tables -------------------------------------------------------
overlap_summary <- tibble(
  tears_total   = length(tear_set),
  plasma_total  = length(plasma_set),
  shared        = length(shared_proteins),
  tears_unique  = length(tear_only),
  plasma_unique = length(plasma_only)
)

write_csv(overlap_summary, file.path(out_dir, "overlap_summary.csv"))
write_csv(tibble(protein = shared_proteins), file.path(out_dir, "shared_proteins.csv"))
write_csv(tibble(protein = tear_only),       file.path(out_dir, "tear_unique_proteins.csv"))
write_csv(tibble(protein = plasma_only),     file.path(out_dir, "plasma_unique_proteins.csv"))

# 9) Venn diagram --------------------------------------------------------------
# Geometry version, similar to the style you asked for.
ellipse_df <- function(x0, y0, a, b, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(
    x = x0 + a * cos(t),
    y = y0 + b * sin(t)
  )
}

plot_venn_project3 <- function(set1, set2, labels = c("Tears DIA", "Plasma DIA")) {
  n1_only <- length(setdiff(set1, set2))
  n2_only <- length(setdiff(set2, set1))
  common  <- length(intersect(set1, set2))
  
  left_ellipse <- ellipse_df(0, 0, 1.4, 0.9) %>% mutate(fluid = labels[1])
  right_ellipse <- ellipse_df(1.2, 0, 1.4, 0.9) %>% mutate(fluid = labels[2])
  ellipses <- bind_rows(left_ellipse, right_ellipse)
  
  ggplot() +
    geom_polygon(
      data = ellipses,
      aes(x = x, y = y, group = fluid, fill = fluid),
      alpha = 0.5,
      color = "black",
      linewidth = 0.5
    ) +
    annotate("text", x = -0.5, y = 0,    label = comma(n1_only), fontface = "bold", size = 4.5) +
    annotate("text", x =  0.6, y = 0,    label = comma(common),  fontface = "bold", size = 4.5) +
    annotate("text", x =  1.7, y = 0,    label = comma(n2_only), fontface = "bold", size = 4.5) +
    annotate("text", x = -0.8, y = 1.05, label = labels[1], fontface = "bold", size = 5, color = col_tears) +
    annotate("text", x =  2.0, y = 1.05, label = labels[2], fontface = "bold", size = 5, color = col_plasma) +
    scale_fill_manual(values = setNames(c(col_tears, col_plasma), labels)) +
    coord_equal() +
    labs(title = "Protein Identification Overlap") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
    )
}

fig_venn <- plot_venn_project3(tear_set, plasma_set, c("Tears DIA", "Plasma DIA"))

ggsave(file.path(out_dir, "Figure1_Venn_Tears_vs_Plasma.pdf"),
       fig_venn, width = 6.5, height = 5.2)

# 10) Mean abundance summaries -------------------------------------------------
tear_summary   <- summarise_biofluid(tears_matrix, tear_sample_cols)
plasma_summary <- summarise_biofluid(plasma_matrix, plasma_sample_cols)

write_csv(tear_summary,   file.path(out_dir, "tear_mean_abundance_summary.csv"))
write_csv(plasma_summary, file.path(out_dir, "plasma_mean_abundance_summary.csv"))

# 11) Concordance plot for shared proteins ------------------------------------
shared_summary <- tear_summary %>%
  rename(
    mean_abundance_tears = mean_abundance,
    n_present_tears      = n_present
  ) %>%
  inner_join(
    plasma_summary %>%
      rename(
        mean_abundance_plasma = mean_abundance,
        n_present_plasma      = n_present
      ),
    by = "protein"
  ) %>%
  filter(n_present_tears > 0, n_present_plasma > 0) %>%
  mutate(
    log2FC = log2(mean_abundance_tears + 1) - log2(mean_abundance_plasma + 1),
    class = case_when(
      log2FC >  1 ~ "tear-enriched",
      log2FC < -1 ~ "plasma-enriched / systemic-like",
      TRUE        ~ "shared / balanced"
    )
  )

# Label the strongest discordant proteins, to keep the plot readable.
top_labels <- shared_summary %>%
  arrange(desc(abs(log2FC))) %>%
  slice_head(n = 15)

cor_test <- cor.test(
  log2(shared_summary$mean_abundance_tears + 1),
  log2(shared_summary$mean_abundance_plasma + 1),
  method = "pearson"
)

fig_concordance <- ggplot(
  shared_summary,
  aes(x = log2(mean_abundance_tears + 1),
      y = log2(mean_abundance_plasma + 1))
) +
  geom_point(alpha = 0.65, size = 1.8, colour = col_plasma) +
  geom_smooth(method = "lm", se = FALSE, colour = col_tears, linewidth = 0.8) +
  geom_text_repel(
    data = top_labels,
    aes(label = protein),
    size = 3,
    max.overlaps = 20
  ) +
  annotate(
    "text",
    x = min(log2(shared_summary$mean_abundance_tears + 1), na.rm = TRUE),
    y = max(log2(shared_summary$mean_abundance_plasma + 1), na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    label = paste0(
      "r = ", round(unname(cor_test$estimate), 3),
      "\np = ", signif(cor_test$p.value, 3),
      "\nn = ", nrow(shared_summary)
    ),
    fontface = "bold",
    size = 3.6
  ) +
  labs(
    title = "Concordance of shared proteins",
    x = "log2 mean abundance in Tears + 1",
    y = "log2 mean abundance in Plasma + 1"
  ) +
  theme_pub()

ggsave(file.path(out_dir, "Figure2_Concordance_SharedProteins.pdf"),
       fig_concordance, width = 7.2, height = 6.2)

write_csv(shared_summary, file.path(out_dir, "shared_proteins_concordance_summary.csv"))

# 12) Mean abundance distribution ---------------------------------------------
density_df <- bind_rows(
  tear_summary   %>% transmute(protein, biofluid = "Tears DIA",   mean_abundance),
  plasma_summary %>% transmute(protein, biofluid = "Plasma DIA",  mean_abundance)
) %>%
  filter(!is.na(mean_abundance), mean_abundance > 0)

fig_distribution <- ggplot(
  density_df,
  aes(x = log2(mean_abundance), fill = biofluid)
) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Tears DIA" = col_tears, "Plasma DIA" = col_plasma)) +
  labs(
    title = "Mean abundance distribution",
    x = "log2(mean abundance)",
    y = "Density",
    fill = NULL
  ) +
  theme_pub()

ggsave(file.path(out_dir, "Figure3_Mean_Abundance_Distribution.pdf"),
       fig_distribution, width = 7.2, height = 5.8)

# 13) Readable heatmap of top shared proteins ---------------------------------
#   - choose only top proteins
#   - order rows by effect size
#   - impute only for plotting
#   - do NOT cluster the sample columns

make_shared_heatmap <- function(shared_summary, tears_matrix, plasma_matrix, top_n = 30) {
  
  top_prots <- shared_summary %>%
    arrange(desc(abs(log2FC))) %>%
    slice_head(n = top_n) %>%
    pull(protein)
  
  if (length(top_prots) < 2) {
    message("Not enough shared proteins for heatmap.")
    return(invisible(NULL))
  }
  
  get_mat <- function(wide_df, proteins, sample_cols) {
    m <- wide_df %>%
      filter(protein %in% proteins) %>%
      arrange(match(protein, proteins)) %>%
      select(protein, all_of(sample_cols)) %>%
      column_to_rownames("protein") %>%
      as.matrix()
    
    m
  }
  
  m_tears  <- get_mat(tears_matrix,  top_prots, tear_sample_cols)
  m_plasma <- get_mat(plasma_matrix, top_prots, plasma_sample_cols)
  
  # Align the same proteins in the same row order.
  common <- intersect(rownames(m_tears), rownames(m_plasma))
  m_tears  <- m_tears[common, , drop = FALSE]
  m_plasma <- m_plasma[common, , drop = FALSE]
  
  combined <- cbind(m_tears, m_plasma)
  combined <- impute_and_clean(combined)
  
  # Re-order rows by absolute effect size, so the strongest differences are on top.
  row_order <- shared_summary %>%
    filter(protein %in% rownames(combined)) %>%
    arrange(desc(abs(log2FC))) %>%
    pull(protein)
  
  combined <- combined[row_order, , drop = FALSE]
  
  # Row scaling, just like the “sample heatmap” style in Project 5.
  z <- t(scale(t(combined)))
  rownames(z) <- rownames(combined)
  
  ann_col <- data.frame(
    Biofluid = factor(
      c(rep("Tears DIA", ncol(m_tears)), rep("Plasma DIA", ncol(m_plasma))),
      levels = c("Tears DIA", "Plasma DIA")
    )
  )
  rownames(ann_col) <- colnames(z)
  
  ann_colors <- list(
    Biofluid = setNames(c(col_tears, col_plasma), c("Tears DIA", "Plasma DIA"))
  )
  
  pdf(
    file.path(out_dir, "Figure4_TopSharedProteins_Heatmap.pdf"),
    width = 12,
    height = max(8, 0.28 * nrow(z))
  )
  
  pheatmap(
    z,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    fontsize_row = 6,
    fontsize_col = 7,
    border_color = NA,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    main = paste0("Top ", nrow(z), " shared proteins\nordered by |log2FC|")
  )
  
  dev.off()
  
  invisible(z)
}

heatmap_mat <- make_shared_heatmap(shared_summary, tears_matrix, plasma_matrix, top_n = 30)

# 14) Completeness -------------------------------------------------------------
tears_completeness <- tears_matrix %>%
  mutate(
    detected_replicates = rowSums(
      across(all_of(tear_sample_cols), ~ !is.na(.x) & .x > 0)
    )
  ) %>%
  filter(detected_replicates > 0) %>%
  select(protein, detected_replicates)

plasma_completeness <- plasma_matrix %>%
  mutate(
    detected_replicates = rowSums(
      across(all_of(plasma_sample_cols), ~ !is.na(.x) & .x > 0)
    )
  ) %>%
  select(protein, detected_replicates)

fig_completeness_tears <- ggplot(
  tears_completeness,
  aes(x = factor(detected_replicates, levels = 0:length(tear_sample_cols)))
) +
  geom_bar(
    fill = col_tears,
    color = "black",
    linewidth = 0.25
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Protein Completeness Distribution",
    subtitle = "Tears DIA",
    x = "Detected replicates",
    y = "Number of proteins"
  ) +
  theme_pub()

fig_completeness_plasma <- ggplot(
  plasma_completeness,
  aes(x = detected_replicates)
) +
  geom_histogram(
    bins = length(plasma_sample_cols),
    fill = col_plasma,
    color = "black",
    linewidth = 0.25
  ) +
  labs(
    title = "Protein Completeness Distribution",
    subtitle = "Depleted Plasma DIA",
    x = "Detected replicates",
    y = "Number of proteins"
  ) +
  theme_pub()

ggsave(
  file.path(out_dir, "Figure5A_Completeness_Tears.pdf"),
  fig_completeness_tears,
  width = 7,
  height = 5
)

ggsave(
  file.path(out_dir, "Figure5B_Completeness_Plasma.pdf"),
  fig_completeness_plasma,
  width = 7,
  height = 5
)
# 15) Export the main shared-annotation table ---------------------------------
write_csv(shared_summary, file.path(out_dir, "shared_proteins_summary.csv"))

# 16) Optional quick interpretation table -------------------------------------
# This is a simple heuristic, not a biological truth.
# Positive log2FC = higher in tears.
# Negative log2FC = higher in plasma, often more “systemic-like”.
shared_summary_classified <- shared_summary %>%
  mutate(
    class_simple = case_when(
      log2FC >  1 ~ "tear-enriched",
      log2FC < -1 ~ "plasma-enriched / systemic-like",
      TRUE        ~ "shared / balanced"
    )
  )

write_csv(shared_summary_classified, file.path(out_dir, "shared_proteins_classified.csv"))

cat("\nAll results saved to:\n", normalizePath(out_dir), "\n")
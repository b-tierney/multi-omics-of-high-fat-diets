#!/usr/bin/env Rscript

# =============================================================================
# scRNA-seq DEG calling and plotting
# =============================================================================
# Purpose
# - Reproduce single-cell differential expression analyses and figure panels:
#   Fig 7a, 7b, 7c, 7d and Supplementary Fig 6a, 6b.
#
# Starting object
# - scrnaseq_seurat_obj.RDS (Seurat object, see required metadata below)
#
# Required Seurat metadata columns
# - mice_code, diet_time, sex, celltype, orig.ident
#
# Key analysis choices
# - Male-only DE analyses.
# - Palm_12mo excluded as outlier.
# - DESeq2 uses RNA raw counts (slot = "counts") in pseudobulk.
# - Wilcoxon DE uses Libra::run_de with per-diet vs Control comparisons.
#
# Main outputs written by this script
# - DESeq2 tables:
#   DESeq2_pseudobulk_males_rev_final.rds
#   DESeq2_pseudobulk_males_rev_final.csv
# - Wilcoxon tables:
#   DE_male_mice_singlecell_wilcox_final.rds
#   DE_<diet>_seuratDEG2_wilcox_<diet>_final.(rds/txt)
# - Figure panels (PDF):
#   plots/fig7a_single_cell_umap_final.pdf
#   plots/fig7b_deseq2_stacked_bar_final.pdf
#   plots/fig7c_deseq2_individual_degs_final.pdf
#   plots/fig7d_mhc2_deg_dotplot_faceted_final.pdf
#   plots/supp6a_wilcox_stacked_bar_final.pdf
#   plots/supp6b_wilcox_individual_degs_final.pdf
#
# Notes for reproducibility
# - Script is written to be non-interactive (Rscript compatible).
# - Some memory cleanup steps (rm/gc) are included to reduce peak RAM usage.
# - Optional dependency: cowplot (used only for Fig 7d theme fallback).
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(DESeq2)
  library(Matrix)
  library(Libra)
  library(future)
})

set.seed(123)
options(stringsAsFactors = FALSE)
plan(sequential)

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 0) Load Seurat object
# -------------------------------
seurat_obj <- readRDS("scrnaseq_seurat_obj.RDS")
stopifnot(all(c("mice_code", "diet_time", "sex", "celltype", "orig.ident") %in% colnames(seurat_obj@meta.data)))

# -------------------------------
# 1) Fig 7a: UMAP of Seurat cells
# -------------------------------
cells_non_palm <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$diet_time != "Palm_12mo"]
p_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "celltype",
  cells = cells_non_palm,
  label = TRUE,
  repel = TRUE
) + ggtitle("UMAP Colored by Cell Type")

ggsave("plots/fig7a_single_cell_umap_final.pdf", p_umap, width = 10, height = 7)

# Drop the full object after UMAP to reduce peak memory before DE.
rm(p_umap, seurat_obj, cells_non_palm)
gc()

# -------------------------------
# 2) DESeq2 pseudobulk DEGs (male)
# -------------------------------
# Build pseudobulk matrices by summing RNA counts across cells grouped by
# (celltype, mouse, diet). Returns a sparse count matrix and sample metadata.
make_pseudobulk <- function(seu, cells_use = NULL) {
  meta <- seu@meta.data
  if (!is.null(cells_use)) {
    meta <- meta[cells_use, , drop = FALSE]
  }
  counts <- GetAssayData(seu, assay = "RNA", slot = "counts")[, rownames(meta), drop = FALSE]

  meta$group_id <- paste(meta$celltype, meta$mice_code, meta$diet_time, sep = "_")
  groups <- unique(meta$group_id)

  pb_list <- vector("list", length(groups))
  names(pb_list) <- groups
  for (g in groups) {
    cells <- rownames(meta)[meta$group_id == g]
    if (length(cells) == 0) next
    pb_list[[g]] <- Matrix::rowSums(counts[, cells, drop = FALSE])
  }

  pb_mat <- do.call(cbind, pb_list)
  pb_mat <- as(pb_mat, "dgCMatrix")

  pb_meta <- meta %>%
    dplyr::select(group_id, mice_code, diet_time, sex, celltype) %>%
    distinct() %>%
    filter(group_id %in% colnames(pb_mat))

  pb_mat <- pb_mat[, pb_meta$group_id]
  list(counts = pb_mat, meta = pb_meta)
}

run_deseq2 <- function(pb_counts, pb_meta, celltype_use, diet_use) {
  meta_sub <- pb_meta %>%
    filter(celltype == celltype_use, diet_time %in% c("Control", diet_use))

  if (length(unique(meta_sub$mice_code[meta_sub$diet_time == "Control"])) < 2 ||
      length(unique(meta_sub$mice_code[meta_sub$diet_time == diet_use])) < 2) {
    message("Skipping (not enough mice): ", celltype_use, " ", diet_use)
    return(NULL)
  }

  counts_sub <- pb_counts[, meta_sub$group_id, drop = FALSE]
  rownames(meta_sub) <- meta_sub$group_id
  stopifnot(all(colnames(counts_sub) == rownames(meta_sub)))

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData = meta_sub,
    design = ~ diet_time
  )

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("diet_time", diet_use, "Control"))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res$celltype <- celltype_use
  res$diet <- diet_use
  res
}

# Reload object for DE steps to avoid holding both full and subsetted copies.
seurat_de <- readRDS("scrnaseq_seurat_obj.RDS")
stopifnot(all(c("mice_code", "diet_time", "sex", "celltype", "orig.ident") %in% colnames(seurat_de@meta.data)))

DefaultAssay(seurat_de) <- "RNA"
gc()

male_cells <- rownames(seurat_de@meta.data)[
  seurat_de@meta.data$sex == "Male" &
    seurat_de@meta.data$diet_time != "Palm_12mo"
]
pb_male <- make_pseudobulk(seurat_de, cells_use = male_cells)

diets_use <- c(
  "Coco_12mo",
  "Fish_12mo", "Fish_4mr", "Fish_9mr",
  "Keto_12mo", "Keto_4mr", "Keto_9mr",
  "Lard_12mo", "Lard_4mr", "Lard_9mr",
  "Olive_12mo"
)

male_results_list <- list()
for (ct in sort(unique(pb_male$meta$celltype))) {
  message("DESeq2 celltype: ", ct)
  for (diet in diets_use) {
    res <- run_deseq2(pb_male$counts, pb_male$meta, ct, diet)
    if (!is.null(res)) {
      male_results_list[[paste(ct, diet, sep = "_")]] <- res
    }
  }
}

male_results <- bind_rows(male_results_list)
saveRDS(pb_male, "pb_male_rev_final.rds")
saveRDS(male_results, "DESeq2_pseudobulk_males_rev_final.rds")
write.csv(male_results, "DESeq2_pseudobulk_males_rev_final.csv", row.names = FALSE)

# Free DESeq2 pseudobulk objects before Wilcoxon.
rm(pb_male, male_results, male_results_list)
gc()

# -------------------------------
# 3) Wilcoxon single-cell DEGs (male)
# -------------------------------
# Run Wilcoxon DE for one diet versus Control using Libra::run_de.
# Output schema follows Libra output columns (e.g., cell_type, gene, avg_logFC).
male_meta <- seurat_de@meta.data[male_cells, , drop = FALSE]
dietTypes <- setdiff(unique(male_meta$diet_time), "Control")

singlecell_wilcox <- function(cur_diet, sc_obj, male_cells_use) {
  message("Wilcoxon diet: ", cur_diet)

  meta_sub <- sc_obj@meta.data[male_cells_use, , drop = FALSE]
  cells_use <- rownames(meta_sub)[meta_sub$diet_time %in% c(cur_diet, "Control")]
  if (length(cells_use) == 0) return(NULL)

  cur_sc <- subset(sc_obj, cells = cells_use)
  cur_sc$diet_time <- ifelse(cur_sc$diet_time == "Control", "Z_Control", cur_sc$diet_time)
  cur_sc$diet_time <- factor(cur_sc$diet_time)

  if (length(unique(cur_sc$diet_time)) < 2) return(NULL)

  ct_tab <- table(cur_sc$celltype)
  keep_ct <- names(ct_tab[ct_tab >= 50])
  cur_sc <- subset(cur_sc, subset = celltype %in% keep_ct)

  de <- run_de(
    cur_sc,
    cell_type_col = "celltype",
    label_col = "diet_time",
    replicate_col = "orig.ident",
    de_family = "singlecell",
    de_method = "wilcox"
  )

  out_prefix <- paste0("DE_", cur_diet)
  saveRDS(de, paste0(out_prefix, "_seuratDEG2_wilcox_", cur_diet, "_final.rds"))
  write.table(
    de,
    file = paste0(out_prefix, "_seuratDEG2_wilcox_", cur_diet, "_final.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  de
}

wilcox_results <- lapply(dietTypes, function(d) singlecell_wilcox(d, seurat_de, male_cells))
names(wilcox_results) <- dietTypes
saveRDS(wilcox_results, "DE_male_mice_singlecell_wilcox_final.rds")

# Free Seurat object prior to plotting from saved result tables.
rm(seurat_de, male_cells, male_meta)
gc()

# -------------------------------
# 4) Plot helpers used for Fig 7 / Supp 6
# -------------------------------
colorsDiet <- c(
  Coco_12mo = "#a65628",
  Fish_12mo = "#f781bf",
  Keto_12mo = "#ff7f00",
  Lard_12mo = "#e41a1c",
  Olive_12mo = "#984ea3",
  Fish_4mr = "#f781bf",
  Fish_9mr = "#f781bf",
  Keto_4mr = "#ff7f00",
  Keto_9mr = "#ff7f00",
  Lard_4mr = "#e41a1c",
  Lard_9mr = "#e41a1c"
)

# Helper to reproduce "individual DEG" panel:
# - restrict to 12-month diets
# - choose top cell types and top genes per cell type
# - point size by adjusted p-value, color by diet
plot_deg_individual <- function(diffdata, out_file, top_celltype_mode = c("top_n", "slice_max")) {
  top_celltype_mode <- match.arg(top_celltype_mode)
  diffdata_12mo <- diffdata %>%
    filter(str_ends(diet, "_12mo")) %>%
    mutate(diet = str_remove(diet, "_12mo"))

  colorsDiet12 <- c(
    Coco = "#a65628",
    Fish = "#f781bf",
    Keto = "#ff7f00",
    Lard = "#e41a1c",
    Olive = "#984ea3"
  )

  top_celltype_tbl <- diffdata_12mo %>%
    filter(!is.na(padj), padj < 0.05) %>%
    dplyr::count(celltype)
  top_celltypes <- if (top_celltype_mode == "top_n") {
    top_celltype_tbl %>% top_n(6, n) %>% pull(celltype)
  } else {
    top_celltype_tbl %>% slice_max(order_by = n, n = 6, with_ties = TRUE) %>% pull(celltype)
  }

  top_genes <- diffdata_12mo %>%
    filter(!is.na(padj), padj < 0.05, celltype %in% top_celltypes) %>%
    group_by(celltype) %>%
    slice_max(order_by = abs(log2FoldChange), n = 5, with_ties = FALSE) %>%
    pull(gene) %>%
    unique()

  plot_data <- diffdata_12mo %>%
    filter(celltype %in% top_celltypes, gene %in% top_genes) %>%
    mutate(padj_plot = pmax(padj, .Machine$double.xmin))

  gene_levels <- plot_data %>%
    group_by(gene) %>%
    summarize(mean_logFC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_logFC) %>%
    pull(gene)

  plot_data <- plot_data %>%
    mutate(
      gene = factor(gene, levels = gene_levels),
      celltype = factor(celltype, levels = top_celltypes)
    )

  p <- ggplot(plot_data, aes(
    x = log2FoldChange,
    y = gene,
    color = diet,
    size = -log10(padj_plot),
    shape = celltype
  )) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dotted") +
    facet_wrap(~celltype, nrow = 1, scales = "free_x", strip.position = "top") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.spacing = grid::unit(0.5, "lines"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(0.25, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_blank()
    ) +
    xlab("log2 Fold Change") +
    scale_color_manual(values = colorsDiet12) +
    scale_size_continuous(name = "-log10(padj)", range = c(1, 6)) +
    scale_shape_manual(values = c(16, 17, 15, 18, 19, 25)) +
    coord_cartesian(clip = "off")

  ggsave(out_file, p, width = 16, height = 6)
}

# -------------------------------
# 5) Fig 7b + Fig 7c from DESeq2
# -------------------------------
# Fig 7b (stacked bar):
# - computed from DESeq2 results with Palm removed
# - dedicated 12mo-only output for manuscript panel
# Fig 7c (individual DEGs):
# - derived from same DESeq2 result table
diffdata_deseq2 <- read.csv("DESeq2_pseudobulk_males_rev_final.csv", stringsAsFactors = FALSE)
diffdata_deseq2 <- diffdata_deseq2 %>% filter(diet != "Palm_12mo")

bar_data_deseq2 <- diffdata_deseq2 %>%
  filter(!is.na(padj), padj < 0.05) %>%
  group_by(celltype, diet, Dir = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  dplyr::count() %>%
  mutate(n = ifelse(Dir == "Increased", n, -n))

total_counts_deseq2 <- bar_data_deseq2 %>%
  group_by(celltype, Dir) %>%
  summarize(total_n = sum(abs(n)), .groups = "drop") %>%
  mutate(
    total_n_posneg = ifelse(Dir == "Increased", total_n, -total_n),
    vjust_label = ifelse(Dir == "Increased", -0.5, 1.5)
  )

celltypes_keep_deseq2 <- bar_data_deseq2 %>%
  group_by(celltype) %>%
  summarize(total_DEGs = sum(abs(n)), .groups = "drop") %>%
  filter(total_DEGs >= 100) %>%
  pull(celltype)

bar_data_deseq2_filt <- bar_data_deseq2 %>% filter(celltype %in% celltypes_keep_deseq2)
total_counts_deseq2_filt <- total_counts_deseq2 %>% filter(celltype %in% celltypes_keep_deseq2)
celltype_order_deseq2 <- bar_data_deseq2_filt %>%
  group_by(celltype) %>%
  summarize(total_DEGs = sum(abs(n)), .groups = "drop") %>%
  arrange(desc(total_DEGs)) %>%
  pull(celltype)

p_deseq2_bar <- ggplot(
  bar_data_deseq2_filt,
  aes(x = factor(celltype, levels = celltype_order_deseq2), y = n, fill = diet)
) +
  geom_bar(stat = "identity") +
  geom_text(
    data = total_counts_deseq2_filt,
    aes(
      x = factor(celltype, levels = celltype_order_deseq2),
      y = total_n_posneg,
      label = total_n,
      vjust = vjust_label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  ylab("Count of DEGs") +
  xlab("") +
  scale_fill_manual(values = colorsDiet) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  scale_y_continuous(
    breaks = pretty_breaks(),
    labels = function(x) abs(x)
  ) +
  ggtitle("Male 12-month DEGs by cell type (>100 DEGs) and diet")

ggsave("plots/degs_singlecell_summary_male_final.pdf", p_deseq2_bar, width = 5.5, height = 6)

# Fig 7b should be 12-month diets only.
bar_data_deseq2_12mo <- bar_data_deseq2 %>% filter(str_ends(diet, "_12mo"))
total_counts_deseq2_12mo <- bar_data_deseq2_12mo %>%
  group_by(celltype, Dir) %>%
  summarize(total_n = sum(abs(n)), .groups = "drop") %>%
  mutate(
    total_n_posneg = ifelse(Dir == "Increased", total_n, -total_n),
    vjust_label = ifelse(Dir == "Increased", -0.5, 1.5)
  )
celltypes_keep_deseq2_12mo <- bar_data_deseq2_12mo %>%
  group_by(celltype) %>%
  summarize(total_DEGs = sum(abs(n)), .groups = "drop") %>%
  filter(total_DEGs >= 100) %>%
  pull(celltype)
bar_data_deseq2_12mo_filt <- bar_data_deseq2_12mo %>% filter(celltype %in% celltypes_keep_deseq2_12mo)
total_counts_deseq2_12mo_filt <- total_counts_deseq2_12mo %>% filter(celltype %in% celltypes_keep_deseq2_12mo)
celltype_order_deseq2_12mo <- bar_data_deseq2_12mo_filt %>%
  group_by(celltype) %>%
  summarize(total_DEGs = sum(abs(n)), .groups = "drop") %>%
  arrange(desc(total_DEGs)) %>%
  pull(celltype)

p_fig7b_12mo <- ggplot(
  bar_data_deseq2_12mo_filt,
  aes(x = factor(celltype, levels = celltype_order_deseq2_12mo), y = n, fill = diet)
) +
  geom_bar(stat = "identity") +
  geom_text(
    data = total_counts_deseq2_12mo_filt,
    aes(
      x = factor(celltype, levels = celltype_order_deseq2_12mo),
      y = total_n_posneg,
      label = total_n,
      vjust = vjust_label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  ylab("Count of DEGs") +
  xlab("") +
  scale_fill_manual(values = colorsDiet) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  scale_y_continuous(
    breaks = pretty_breaks(),
    labels = function(x) abs(x)
  ) +
  ggtitle("Fig 7b: DESeq2 DEGs by cell type and diet (12mo only)")

ggsave("plots/fig7b_deseq2_stacked_bar_final.pdf", p_fig7b_12mo, width = 5.5, height = 6)

plot_deg_individual(
  diffdata = diffdata_deseq2,
  out_file = "plots/degs_singlecell_genes_male_final.pdf",
  top_celltype_mode = "top_n"
)
plot_deg_individual(
  diffdata = diffdata_deseq2,
  out_file = "plots/fig7c_deseq2_individual_degs_final.pdf",
  top_celltype_mode = "top_n"
)

# -------------------------------
# 6) Supp 6a + Supp 6b from Wilcoxon
# -------------------------------
# Supplementary Fig 6 is generated from Wilcoxon outputs transformed into a
# DESeq2-like plotting table (columns: log2FoldChange, padj, gene, celltype, diet).
wilcox_list <- readRDS("DE_male_mice_singlecell_wilcox_final.rds")
required_cols <- c("cell_type", "gene", "avg_logFC", "p_val", "p_val_adj")

diffdata_wilcox <- lapply(names(wilcox_list), function(diet_name) {
  df <- wilcox_list[[diet_name]]
  if (!is.data.frame(df)) return(NULL)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Wilcoxon entry '%s' missing columns: %s", diet_name, paste(missing_cols, collapse = ", ")))
  }
  df %>%
    transmute(
      baseMean = NA_real_,
      log2FoldChange = as.numeric(avg_logFC),
      lfcSE = NA_real_,
      stat = NA_real_,
      pvalue = as.numeric(p_val),
      padj = as.numeric(p_val_adj),
      gene = as.character(gene),
      celltype = as.character(cell_type),
      diet = diet_name
    )
}) %>% bind_rows()

diffdata_wilcox_12mo <- diffdata_wilcox %>%
  filter(str_ends(diet, "_12mo"))

bar_data_wilcox <- diffdata_wilcox_12mo %>%
  filter(!is.na(padj), padj < 0.05) %>%
  group_by(celltype, diet, Dir = ifelse(log2FoldChange > 0, "Increased", "Decreased")) %>%
  dplyr::count() %>%
  mutate(n = ifelse(Dir == "Increased", n, -n))

total_counts_wilcox <- bar_data_wilcox %>%
  group_by(celltype, Dir) %>%
  summarize(total_n = sum(abs(n)), .groups = "drop") %>%
  mutate(
    total_n_posneg = ifelse(Dir == "Increased", total_n, -total_n),
    vjust_label = ifelse(Dir == "Increased", -0.5, 1.5)
  )

celltypes_keep_wilcox <- bar_data_wilcox %>%
  group_by(celltype) %>%
  summarize(total_DEGs = sum(abs(n)), .groups = "drop") %>%
  filter(total_DEGs >= 100) %>%
  pull(celltype)

bar_data_wilcox_filt <- bar_data_wilcox %>% filter(celltype %in% celltypes_keep_wilcox)
total_counts_wilcox_filt <- total_counts_wilcox %>% filter(celltype %in% celltypes_keep_wilcox)
celltype_order_wilcox <- bar_data_wilcox_filt %>%
  group_by(celltype) %>%
  summarize(total_DEGs = sum(abs(n)), .groups = "drop") %>%
  arrange(desc(total_DEGs)) %>%
  pull(celltype)

p_wilcox_bar <- ggplot(
  bar_data_wilcox_filt,
  aes(x = factor(celltype, levels = celltype_order_wilcox), y = n, fill = diet)
) +
  geom_bar(stat = "identity") +
  geom_text(
    data = total_counts_wilcox_filt,
    aes(
      x = factor(celltype, levels = celltype_order_wilcox),
      y = total_n_posneg,
      label = total_n,
      vjust = vjust_label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  ylab("Count of DEGs") +
  xlab("") +
  scale_fill_manual(values = colorsDiet) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  scale_y_continuous(
    breaks = pretty_breaks(),
    labels = function(x) abs(x)
  ) +
  ggtitle("Male 12-month DEGs by cell type (>100 DEGs) and diet (Wilcoxon)")

ggsave("plots/supp6a_wilcox_stacked_bar_final.pdf", p_wilcox_bar, width = 10, height = 6)

plot_deg_individual(
  diffdata = diffdata_wilcox,
  out_file = "plots/degs_singlecell_genes_male_wilcox_final.pdf",
  top_celltype_mode = "slice_max"
)
plot_deg_individual(
  diffdata = diffdata_wilcox,
  out_file = "plots/supp6b_wilcox_individual_degs_final.pdf",
  top_celltype_mode = "slice_max"
)

# -------------------------------
# 7) Fig 7d (Wilcoxon MHC-II faceted plot)
# -------------------------------

wilcox_for_fig7d <- readRDS("DE_male_mice_singlecell_wilcox_final.rds")

mhc2_genes_fig7d <- c("H2-Aa", "H2-Ab1", "Cd74")
colors_diet_fig7d <- c(
  Lard = "#e41a1c",
  Keto = "#ff7f00",
  Fish = "#f781bf",
  Coco = "#a65628",
  Olive = "#984ea3",
  Palm = "#4daf4a"
)

combined_fig7d <- bind_rows(
  lapply(names(wilcox_for_fig7d), function(nm) {
    df <- wilcox_for_fig7d[[nm]]
    df$comparison <- nm
    df
  })
)

combined_fig7d <- combined_fig7d %>%
  mutate(
    diet = str_extract(comparison, "^[^_]+"),
    category = case_when(
      str_detect(comparison, "_12mo$") ~ "Diet",
      str_detect(comparison, "_4mr$") ~ "4MR",
      str_detect(comparison, "_9mr$") ~ "9MR",
      TRUE ~ NA_character_
    )
  )

mhc2_data_fig7d <- combined_fig7d %>%
  filter(gene %in% mhc2_genes_fig7d) %>%
  mutate(category = factor(category, levels = c("Diet", "4MR", "9MR")))

mean_bars_facet_fig7d <- mhc2_data_fig7d %>%
  group_by(cell_type, category) %>%
  summarise(mean_lfc = mean(avg_logFC, na.rm = TRUE), .groups = "drop")

ct_order_fig7d <- mhc2_data_fig7d %>%
  group_by(cell_type) %>%
  summarise(overall = mean(avg_logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(overall) %>%
  pull(cell_type)

mhc2_data_fig7d$cell_type <- factor(mhc2_data_fig7d$cell_type, levels = ct_order_fig7d)
mean_bars_facet_fig7d$cell_type <- factor(mean_bars_facet_fig7d$cell_type, levels = ct_order_fig7d)

base_theme_fig7d <- if (requireNamespace("cowplot", quietly = TRUE)) {
  cowplot::theme_cowplot()
} else {
  theme_minimal(base_size = 11)
}

p_fig7d <- ggplot() +
  geom_col(
    data = mean_bars_facet_fig7d,
    aes(x = mean_lfc, y = cell_type),
    fill = "grey70",
    width = 0.6,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  geom_point(
    data = mhc2_data_fig7d %>% filter(diet != "Palm"),
    aes(x = avg_logFC, y = cell_type, color = diet, shape = gene),
    size = 2.5,
    alpha = 0.8,
    position = position_jitter(width = 0, height = 0.2)
  ) +
  scale_shape_discrete(name = "MHC-II Gene") +
  scale_color_manual(name = "Diet", values = colors_diet_fig7d) +
  facet_wrap(~ category, ncol = 3) +
  labs(
    x = "Log2 Fold Change",
    y = NULL
  ) +
  base_theme_fig7d +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 11)
  )

ggsave("plots/fig7d_mhc2_deg_dotplot_faceted_final.pdf", p_fig7d, width = 14, height = 5)

message("Done: DESeq2 + Wilcoxon DEG generation and Fig7/Supp6/Fig7d panel outputs.")

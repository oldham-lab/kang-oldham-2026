library(qs)
library(data.table)
library(tidyverse)

version <- "v2"

# ============================================================
# PLACEHOLDERS
# ============================================================

# --- DE result files (output from full_DE_pipeline_ADvsCon.R) ---
MIT_EDGER_PATH  <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH  <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH  <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH  <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_DESeq2_ADvsCon_by_celltype.qs")

# --- External gene list ---
# CSV/TSV with at minimum two columns: one for celltype, one for gene name
GENE_LIST_PATH     <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v2/Allct_genes_allComparisons.csv")
GENE_LIST_CT_COL   <- "Celltype"    # column name for celltype in gene list
GENE_LIST_GENE_COL <- "Gene"        # column name for gene in gene list
GENE_LIST_COMP_COL <- "Comparison"  # column name for comparison grouping in gene list

# --- DE thresholds ---
FDR_THRESH <- 0.05
LFC_THRESH <- 0      # applied as |logFC| > LFC_THRESH

# --- Output ---
SAVE_DIR <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "DE_overlap_analysis/")

# ============================================================
# HELPERS
# ============================================================

# Standardise edgeR result to common column names
tidy_edger <- function(res_list) {
  lapply(res_list, function(df) {
    df |>
      dplyr::rename(
        gene   = any_of(c("genes", "gene", "Gene")),
        logFC  = any_of(c("logFC")),
        FDR    = any_of(c("FDR", "adj.P.Val", "padj"))
      ) |>
      dplyr::select(gene, logFC, FDR)
  })
}

# Standardise DESeq2 result to common column names
tidy_deseq <- function(res_list) {
  lapply(res_list, function(df) {
    df |>
      dplyr::rename(
        gene   = any_of(c("genes", "gene", "Gene")),
        logFC  = any_of(c("logFC", "log2FoldChange")),
        FDR    = any_of(c("FDR", "padj"))
      ) |>
      dplyr::select(gene, logFC, FDR)
  })
}

# Apply FDR + LFC thresholds and return DE gene names
get_de_genes <- function(df, fdr_thresh, lfc_thresh) {
  df |>
    dplyr::filter(!is.na(FDR), !is.na(logFC),
                  FDR < fdr_thresh,
                  abs(logFC) > lfc_thresh) |>
    dplyr::pull(gene) |>
    unique()
}

# Fisher's exact test for overlap enrichment.
# Tests whether the overlap between DE genes and the reference gene list
# is greater than expected by chance given the background (all tested genes).
#   a = in DE and in list
#   b = in DE but not in list
#   c = not DE but in list
#   d = not DE and not in list
fisher_overlap <- function(de_genes, list_genes, background_genes) {
  a <- length(intersect(de_genes,   list_genes))
  b <- length(setdiff(de_genes,     list_genes))
  c <- length(intersect(setdiff(background_genes, de_genes), list_genes))
  d <- length(setdiff(background_genes, union(de_genes, list_genes)))
  mat <- matrix(c(a, b, c, d), nrow = 2,
                dimnames = list(c("DE", "not_DE"), c("in_list", "not_in_list")))
  ft  <- fisher.test(mat, alternative = "greater")
  data.frame(
    n_de           = length(de_genes),
    n_list         = length(list_genes),
    n_background   = length(background_genes),
    n_overlap      = a,
    pct_de_in_list = ifelse(length(de_genes) > 0, round(100 * a / length(de_genes), 1), NA),
    pct_list_in_de = ifelse(length(list_genes) > 0, round(100 * a / length(list_genes), 1), NA),
    odds_ratio     = round(ft$estimate, 3),
    p_value        = ft$p.value,
    overlap_genes  = paste(sort(intersect(de_genes, list_genes)), collapse = ";")
  )
}

# Run overlap analysis across all celltypes for one DE result set
compare_to_list <- function(de_list,         # named list of per-celltype DE data frames
                             gene_list_by_ct, # named list of gene vectors (from external list)
                             fdr_thresh,
                             lfc_thresh,
                             label) {         # e.g. "MIT_edgeR"
  # Union of celltypes present in both
  cts <- union(names(de_list), names(gene_list_by_ct))

  results <- lapply(cts, function(ct) {
    de_df      <- de_list[[ct]]
    list_genes <- gene_list_by_ct[[ct]]

    # If celltype missing from DE results, skip
    if (is.null(de_df)) {
      message("  ", label, " | ", ct, ": no DE results — skipping")
      return(NULL)
    }

    background <- unique(de_df$gene)

    # If celltype not in external list, report NA row
    if (is.null(list_genes) || length(list_genes) == 0) {
      message("  ", label, " | ", ct, ": not in external gene list — skipping")
      return(NULL)
    }

    de_genes <- get_de_genes(de_df, fdr_thresh, lfc_thresh)
    row      <- fisher_overlap(de_genes, list_genes, background)
    row$celltype <- ct
    row$dataset  <- label
    row
  })

  dplyr::bind_rows(results)
}

# ============================================================
# LOAD DATA
# ============================================================
message("Loading DE results ...")
mit_edger <- tidy_edger(readRDS(MIT_EDGER_PATH))
mit_deseq <- tidy_deseq(qread(MIT_DESEQ_PATH))
sea_edger <- tidy_edger(readRDS(SEA_EDGER_PATH))
sea_deseq <- tidy_deseq(qread(SEA_DESEQ_PATH))

message("Loading external gene list ...")
gene_list_raw <- fread(GENE_LIST_PATH, data.table = FALSE)

for (col in c(GENE_LIST_CT_COL, GENE_LIST_GENE_COL, GENE_LIST_COMP_COL)) {
  if (!col %in% colnames(gene_list_raw))
    stop("Column '", col, "' not found in gene list. ",
         "Available columns: ", paste(colnames(gene_list_raw), collapse = ", "))
}

# Convert to nested named list: comparison -> celltype -> character vector of genes
gene_list_nested <- gene_list_raw |>
  dplyr::group_by(.data[[GENE_LIST_COMP_COL]], .data[[GENE_LIST_CT_COL]]) |>
  dplyr::summarise(genes = list(unique(.data[[GENE_LIST_GENE_COL]])), .groups = "drop")

# Outer list keyed by comparison, inner list keyed by celltype
comparisons <- unique(gene_list_raw[[GENE_LIST_COMP_COL]])
gene_list_nested <- setNames(lapply(comparisons, function(comp) {
  sub <- gene_list_nested[gene_list_nested[[GENE_LIST_COMP_COL]] == comp, ]
  setNames(sub$genes, sub[[GENE_LIST_CT_COL]])
}), comparisons)

message("  External list: ", length(comparisons), " comparisons, ",
        length(unique(gene_list_raw[[GENE_LIST_CT_COL]])), " celltypes, ",
        nrow(gene_list_raw), " total gene entries")

# ============================================================
# RUN OVERLAP ANALYSIS
# ============================================================
message("Running overlap analysis ...")
if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)

de_sets <- list(
  MIT_edgeR  = mit_edger,
  MIT_DESeq2 = mit_deseq,
  SEA_edgeR  = sea_edger,
  SEA_DESeq2 = sea_deseq
)

# Loop over every combination of DE dataset and comparison group
results_list <- lapply(names(de_sets), function(ds_name) {
  lapply(names(gene_list_nested), function(comp) {
    res <- compare_to_list(de_sets[[ds_name]],
                           gene_list_nested[[comp]],
                           FDR_THRESH, LFC_THRESH,
                           ds_name)
    if (!is.null(res)) res$comparison <- comp
    res
  })
}) |> unlist(recursive = FALSE)

# Combined table across all datasets and comparisons
results_combined <- dplyr::bind_rows(results_list) |>
  dplyr::select(dataset, comparison, celltype, n_de, n_list, n_background,
                n_overlap, pct_de_in_list, pct_list_in_de,
                odds_ratio, p_value, overlap_genes) |>
  dplyr::arrange(dataset, comparison, p_value)

# FDR-correct p-values within each dataset x comparison group
results_combined <- results_combined |>
  dplyr::group_by(dataset, comparison) |>
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::relocate(FDR, .after = p_value)

# ============================================================
# SAVE OUTPUTS
# ============================================================

# Full results table
fwrite(results_combined,
       file = file.path(SAVE_DIR, "overlap_results_all.csv"))
message("Saved: overlap_results_all.csv")

# Per-dataset CSVs (without the long overlap_genes column, for readability)
for (ds in unique(results_combined$dataset)) {
  out <- results_combined |>
    dplyr::filter(dataset == ds) |>
    dplyr::select(-overlap_genes)
  fwrite(out, file = file.path(SAVE_DIR, paste0("overlap_", ds, ".csv")))
  message("Saved: overlap_", ds, ".csv")
}

# Per-comparison CSVs
for (comp in unique(results_combined$comparison)) {
  safe_comp <- gsub("[^A-Za-z0-9_-]", "_", comp)  # sanitise for filename
  out <- results_combined |>
    dplyr::filter(comparison == comp) |>
    dplyr::select(-overlap_genes)
  fwrite(out, file = file.path(SAVE_DIR, paste0("overlap_comparison_", safe_comp, ".csv")))
  message("Saved: overlap_comparison_", safe_comp, ".csv")
}

# Separate file with just the overlap gene lists (long format, easy to parse)
overlap_genes_long <- results_combined |>
  dplyr::select(dataset, comparison, celltype, overlap_genes) |>
  dplyr::filter(nchar(overlap_genes) > 0) |>
  tidyr::separate_rows(overlap_genes, sep = ";") |>
  dplyr::rename(gene = overlap_genes)
fwrite(overlap_genes_long,
       file = file.path(SAVE_DIR, "overlap_genes_long.csv"))
message("Saved: overlap_genes_long.csv")

# ============================================================
# SUMMARY PRINT
# ============================================================
message("\n--- Summary (FDR < 0.05) ---")
results_combined |>
  dplyr::filter(FDR < 0.05) |>
  dplyr::select(dataset, comparison, celltype, n_overlap, pct_list_in_de, odds_ratio, FDR) |>
  as.data.frame() |>
  print()

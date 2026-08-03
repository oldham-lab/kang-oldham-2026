library(qs)
library(data.table)
library(tidyverse)

# ============================================================
# PLACEHOLDERS
# ============================================================

MIT_EDGER_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_DESeq2_ADvsCon_by_celltype.qs")

DCOPA_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v2/dcopa_overlap.csv")
SAVE_DIR   <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "dcopa_overlap/")

# DE thresholds
FDR_THRESH <- 0.05
LFC_THRESH <- 1      # |logFC| > 1

# ============================================================
# LOAD & PREPARE DCOPA GENE LIST
# ============================================================
message("Loading dcopa gene list ...")

dcopa_raw <- fread(DCOPA_PATH, data.table = FALSE)

# Filter to AllADVsCon_DFC and celltypes with "all" suffix only
dcopa <- dcopa_raw |>
  dplyr::filter(
    Comparison == "AllADVsCon_DFC",
    grepl(" all$", Celltype)
  ) |>
  dplyr::mutate(
    # Strip " all" suffix so celltype names match DE result keys
    celltype_clean = sub(" all$", "", Celltype),
    direction      = dplyr::case_when(
      Direction == "Higher in more severe" ~ "higher",
      Direction == "Lower in more severe"  ~ "lower",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(direction))

message("  dcopa: ", length(unique(dcopa$celltype_clean)), " celltypes after filtering")

# Split into per-direction, per-celltype gene vectors
dcopa_higher <- split(dcopa$Gene[dcopa$direction == "higher"],
                      dcopa$celltype_clean[dcopa$direction == "higher"])
dcopa_lower  <- split(dcopa$Gene[dcopa$direction == "lower"],
                      dcopa$celltype_clean[dcopa$direction == "lower"])
dcopa_higher <- lapply(dcopa_higher, unique)
dcopa_lower  <- lapply(dcopa_lower,  unique)

# ============================================================
# LOAD & STANDARDISE DE RESULTS
# ============================================================
message("Loading DE results ...")

# Standardise column names to: gene, logFC, FDR
tidy_edger <- function(res_list) {
  lapply(res_list, function(df) {
    df |>
      dplyr::rename(
        gene  = any_of(c("genes", "gene", "Gene")),
        logFC = any_of(c("logFC")),
        FDR   = any_of(c("FDR", "adj.P.Val"))
      ) |>
      dplyr::select(gene, logFC, FDR)
  })
}

tidy_deseq <- function(res_list) {
  lapply(res_list, function(df) {
    df |>
      dplyr::rename(
        gene  = any_of(c("genes", "gene", "Gene")),
        logFC = any_of(c("logFC", "log2FoldChange")),
        FDR   = any_of(c("FDR", "padj"))
      ) |>
      dplyr::select(gene, logFC, FDR)
  })
}

mit_edger <- tidy_edger(readRDS(MIT_EDGER_PATH))
mit_deseq <- tidy_deseq(qread(MIT_DESEQ_PATH))
sea_edger <- tidy_edger(readRDS(SEA_EDGER_PATH))
sea_deseq <- tidy_deseq(qread(SEA_DESEQ_PATH))

# Diagnostic: print celltype names from DE results and dcopa so mismatches are visible
message("\n--- Celltype name check ---")
message("dcopa celltypes (after stripping ' all'): ",
        paste(sort(unique(dcopa$celltype_clean)), collapse = ", "))
message("MIT edgeR celltypes: ",  paste(sort(names(mit_edger)), collapse = ", "))
message("SEA edgeR celltypes: ",  paste(sort(names(sea_edger)), collapse = ", "))

# Check for any matches at all — warn loudly if none
dcopa_cts <- unique(dcopa$celltype_clean)
de_cts_all <- unique(c(names(mit_edger), names(mit_deseq),
                        names(sea_edger), names(sea_deseq)))
matched <- intersect(dcopa_cts, de_cts_all)
if (length(matched) == 0) {
  stop(
    "No celltype names match between dcopa and DE results.\n",
    "dcopa example names: ", paste(head(dcopa_cts, 5), collapse = ", "), "\n",
    "DE result example names: ", paste(head(de_cts_all, 5), collapse = ", "), "\n",
    "Check capitalisation, spacing, or whether the DE results use different celltype labels."
  )
} else {
  message("Matched ", length(matched), " celltypes: ", paste(sort(matched), collapse = ", "))
  unmatched_dcopa <- setdiff(dcopa_cts, de_cts_all)
  if (length(unmatched_dcopa) > 0)
    message("dcopa celltypes with no DE results: ",
            paste(sort(unmatched_dcopa), collapse = ", "))
}
message("---\n")

# ============================================================
# HELPER: extract directional DE gene sets from one DE result
# ============================================================
get_de_directional <- function(de_df, fdr_thresh, lfc_thresh) {
  sig <- de_df |>
    dplyr::filter(!is.na(FDR), !is.na(logFC),
                  FDR < fdr_thresh,
                  abs(logFC) > lfc_thresh)
  list(
    higher = unique(sig$gene[sig$logFC > 0]),
    lower  = unique(sig$gene[sig$logFC < 0])
  )
}

# ============================================================
# MAIN ANALYSIS
# ============================================================
# All celltypes present across dcopa and any DE result
all_cts <- sort(unique(c(
  names(dcopa_higher), names(dcopa_lower),
  names(mit_edger), names(mit_deseq),
  names(sea_edger), names(sea_deseq)
)))

de_sets <- list(
  MIT_edgeR  = mit_edger,
  MIT_DESeq2 = mit_deseq,
  SEA_edgeR  = sea_edger,
  SEA_DESeq2 = sea_deseq
)

# ── TABLE 1: Summary counts and % overlap per celltype ───────────────────────
message("Building summary table ...")

summary_rows <- lapply(all_cts, function(ct) {

  # dcopa gene counts
  n_dcopa_higher <- length(dcopa_higher[[ct]])
  n_dcopa_lower  <- length(dcopa_lower[[ct]])

  row <- data.frame(celltype = ct,
                    n_dcopa_higher = n_dcopa_higher,
                    n_dcopa_lower  = n_dcopa_lower)

  for (ds_name in names(de_sets)) {
    de_df <- de_sets[[ds_name]][[ct]]

    if (is.null(de_df)) {
      # Celltype not tested in this dataset
      row[[paste0("n_DE_higher_", ds_name)]] <- NA_integer_
      row[[paste0("n_DE_lower_",  ds_name)]] <- NA_integer_
      row[[paste0("pct_overlap_higher_", ds_name)]] <- NA_real_
      row[[paste0("pct_overlap_lower_",  ds_name)]] <- NA_real_
    } else {
      directional  <- get_de_directional(de_df, FDR_THRESH, LFC_THRESH)
      de_higher    <- directional$higher
      de_lower     <- directional$lower
      dcopa_h      <- dcopa_higher[[ct]]
      dcopa_l      <- dcopa_lower[[ct]]

      n_overlap_higher <- length(intersect(de_higher, dcopa_h))
      n_overlap_lower  <- length(intersect(de_lower,  dcopa_l))

      # % of dcopa genes recovered in DE (what fraction of the reference does DE capture)
      pct_h <- if (length(dcopa_h) > 0) round(100 * n_overlap_higher / length(dcopa_h), 1) else NA_real_
      pct_l <- if (length(dcopa_l) > 0) round(100 * n_overlap_lower  / length(dcopa_l), 1) else NA_real_

      row[[paste0("n_DE_higher_", ds_name)]] <- length(de_higher)
      row[[paste0("n_DE_lower_",  ds_name)]] <- length(de_lower)
      row[[paste0("pct_overlap_higher_", ds_name)]] <- pct_h
      row[[paste0("pct_overlap_lower_",  ds_name)]] <- pct_l
    }
  }
  row
})

summary_table <- dplyr::bind_rows(summary_rows)

# Reorder columns for readability:
# celltype | DE counts (all datasets, higher then lower) |
#          | dcopa counts | % overlaps
count_cols   <- c(
  outer(c("n_DE_higher_", "n_DE_lower_"), names(de_sets), paste0) |> as.vector()
)
dcopa_cols   <- c("n_dcopa_higher", "n_dcopa_lower")
overlap_cols <- c(
  outer(c("pct_overlap_higher_", "pct_overlap_lower_"), names(de_sets), paste0) |> as.vector()
)

summary_table <- summary_table |>
  dplyr::select(celltype, all_of(count_cols), all_of(dcopa_cols), all_of(overlap_cols))

# ── TABLE 2: Overlapping gene names, long format ─────────────────────────────
message("Building overlap gene table ...")

overlap_rows <- lapply(all_cts, function(ct) {
  lapply(names(de_sets), function(ds_name) {
    de_df <- de_sets[[ds_name]][[ct]]
    if (is.null(de_df)) return(NULL)

    directional <- get_de_directional(de_df, FDR_THRESH, LFC_THRESH)

    higher_genes <- intersect(directional$higher, dcopa_higher[[ct]])
    lower_genes  <- intersect(directional$lower,  dcopa_lower[[ct]])

    rows <- list()
    if (length(higher_genes) > 0)
      rows[[1]] <- data.frame(celltype  = ct,
                              dataset   = ds_name,
                              direction = "higher_in_AD",
                              gene      = higher_genes)
    if (length(lower_genes) > 0)
      rows[[2]] <- data.frame(celltype  = ct,
                              dataset   = ds_name,
                              direction = "lower_in_AD",
                              gene      = lower_genes)
    dplyr::bind_rows(rows)
  }) |> dplyr::bind_rows()
}) |> dplyr::bind_rows()

if (nrow(overlap_rows) == 0) {
  warning("No overlapping genes found between DE results and dcopa list. ",
          "Check that celltype names match (see diagnostic output above).")
  overlap_table <- data.frame(celltype  = character(),
                               dataset   = character(),
                               direction = character(),
                               gene      = character())
} else {
  overlap_table <- overlap_rows |>
    dplyr::arrange(celltype, direction, dataset, gene)
}

# ============================================================
# SAVE OUTPUTS
# ============================================================
if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)

fwrite(summary_table,
       file = file.path(SAVE_DIR, "dcopa_overlap_summary.csv"))
message("Saved: dcopa_overlap_summary.csv")

fwrite(overlap_table,
       file = file.path(SAVE_DIR, "dcopa_overlap_genes.csv"))
message("Saved: dcopa_overlap_genes.csv")

# ============================================================
# PRINT PREVIEW
# ============================================================
message("\n--- Summary table (", nrow(summary_table), " celltypes) ---")
print(summary_table, n = Inf)

message("\n--- Overlap genes (", nrow(overlap_table), " entries across all datasets/celltypes) ---")
print(head(overlap_table, 30))

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

FDR_THRESH <- 0.05
LFC_THRESH <- 0      # applied as |logFC| > LFC_THRESH

# ============================================================
# HELPERS
# ============================================================

tidy_edger <- function(res_list) {
  lapply(res_list, function(df) {
    df |>
      dplyr::rename(
        gene  = any_of(c("genes", "gene", "Gene")),
        logFC = any_of(c("logFC")),
        FDR   = any_of(c("FDR", "adj.P.Val", "padj"))
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

# ============================================================
# LOAD DATA
# ============================================================

message("Loading DE results ...")
datasets <- list(
  MIT_edgeR  = tidy_edger(readRDS(MIT_EDGER_PATH)),
  MIT_DESeq2 = tidy_deseq(qread(MIT_DESEQ_PATH)),
  SEA_edgeR  = tidy_edger(readRDS(SEA_EDGER_PATH)),
  SEA_DESeq2 = tidy_deseq(qread(SEA_DESEQ_PATH))
)

# ============================================================
# DISPLAY DE GENES (count + top 5 by FDR)
# ============================================================

for (ds_name in names(datasets)) {
  cat("\n===", ds_name, "===\n")
  for (ct in names(datasets[[ds_name]])) {
    de <- datasets[[ds_name]][[ct]] |>
      dplyr::filter(!is.na(FDR), !is.na(logFC),
                    FDR < FDR_THRESH,
                    abs(logFC) > LFC_THRESH) |>
      dplyr::arrange(FDR)
    n    <- nrow(de)
    top5 <- paste(head(de$gene, 5), collapse = ", ")
    cat(sprintf("  %s: %d DE genes  |  top 5 (by FDR): %s\n", ct, n, top5))
  }
}

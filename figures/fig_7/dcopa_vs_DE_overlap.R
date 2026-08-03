# Overlap: dCoPA (gene, celltype, direction) tuples vs significant DE tuples.
# For each (dataset, region) combination, reports what % of dCoPA tuples are also
# significant DE tuples, separately for edgeR and DESeq2.
# Direction is matched: dCoPA "higher"→"up", "lower"→"down" (aligned to DE logFC sign).
# AD_modules (ROSMAP) dCoPA comparisons are excluded.

library(qs)
library(data.table)

# ============================================================
# PATHS
# ============================================================

MIT_EDGER_PATH_DFC  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/MIT_DE_results_DFC/mit_dfc_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH_DFC  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/MIT_DE_results_DFC/mit_dfc_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH_DFC  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/SEA_DE_results_DFC/sea_dfc_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH_DFC  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/SEA_DE_results_DFC/sea_dfc_DESeq2_ADvsCon_by_celltype.qs")
MIT_EDGER_PATH_MTG  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/MIT_DE_results_MTG/mit_mtg_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH_MTG  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/MIT_DE_results_MTG/mit_mtg_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH_MTG  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/SEA_DE_results_MTG/sea_mtg_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH_MTG  <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/DE/SEA_DE_results_MTG/sea_mtg_DESeq2_ADvsCon_by_celltype.qs")

DCOPA_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.1/panel_B_dcopa_genelist.csv")

SAVE_DIR   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7")

FDR_THRESH <- 0.05
LFC_THRESH <- 0

# ============================================================
# LOAD DE DATA (panel A)
# ============================================================

load_edger <- function(path) {
  res <- readRDS(path)
  lapply(res, function(df) data.frame(gene = df$genes, logFC = df$logFC, FDR = df$FDR,
                                      stringsAsFactors = FALSE))
}

load_deseq <- function(path) {
  res <- qread(path)
  lapply(res, function(df) data.frame(gene = df$genes, logFC = df$logFC, FDR = df$FDR,
                                      stringsAsFactors = FALSE))
}

sig_gene_ct <- function(de_list) {
  pairs <- lapply(names(de_list), function(ct) {
    df <- de_list[[ct]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    keep <- !is.na(df$FDR) & !is.na(df$logFC) &
            df$FDR < FDR_THRESH & abs(df$logFC) > LFC_THRESH
    df <- df[keep, ]
    if (nrow(df) == 0) return(NULL)
    direction <- ifelse(df$logFC > 0, "up", "down")
    paste(df$gene, ct, direction, sep = "|||")
  })
  unique(unlist(Filter(Negate(is.null), pairs)))
}

message("Loading DE data ...")
de_sets <- list(
  Gabitto_DFC_edgeR  = sig_gene_ct(load_edger(SEA_EDGER_PATH_DFC)),
  Gabitto_DFC_DESeq2 = sig_gene_ct(load_deseq(SEA_DESEQ_PATH_DFC)),
  Liu_DFC_edgeR      = sig_gene_ct(load_edger(MIT_EDGER_PATH_DFC)),
  Liu_DFC_DESeq2     = sig_gene_ct(load_deseq(MIT_DESEQ_PATH_DFC)),
  Gabitto_MTG_edgeR  = sig_gene_ct(load_edger(SEA_EDGER_PATH_MTG)),
  Gabitto_MTG_DESeq2 = sig_gene_ct(load_deseq(SEA_DESEQ_PATH_MTG)),
  Liu_MTG_edgeR      = sig_gene_ct(load_edger(MIT_EDGER_PATH_MTG)),
  Liu_MTG_DESeq2     = sig_gene_ct(load_deseq(MIT_DESEQ_PATH_MTG))
)

# ============================================================
# LOAD dCoPA DATA (panel B) — non-ROSMAP comparisons only
# ============================================================

message("Loading dCoPA data ...")
dcopa_raw <- fread(DCOPA_PATH, data.table = FALSE)

DCOPA_COMPS <- c(
  "Gabitto_AllADVsCon_DFC",
  "Liu_AllADVsCon_DFC",
  "Gabitto_AllADVsCon_MTG",
  "Liu_AllADVsCon_MTG"
)

dcopa_sets <- setNames(lapply(DCOPA_COMPS, function(comp) {
  rows <- dcopa_raw[dcopa_raw$Comparison == comp, ]
  direction <- ifelse(rows$Direction == "Higher in more severe", "up", "down")
  unique(paste(rows$Gene, rows$Celltype, direction, sep = "|||"))
}), DCOPA_COMPS)

# ============================================================
# COMPUTE OVERLAPS
# ============================================================

results <- do.call(rbind, lapply(DCOPA_COMPS, function(comp) {
  dataset <- sub("_AllADVsCon_.*", "", comp)
  region  <- sub(".*_AllADVsCon_", "", comp)

  dcopa_pairs <- dcopa_sets[[comp]]
  n_dcopa <- length(dcopa_pairs)

  do.call(rbind, lapply(c("edgeR", "DESeq2"), function(method) {
    de_key    <- paste(dataset, region, method, sep = "_")
    n_overlap <- sum(dcopa_pairs %in% de_sets[[de_key]])
    data.frame(
      dataset             = dataset,
      region              = region,
      de_method           = method,
      n_dcopa_tuples      = n_dcopa,
      n_de_tuples         = length(de_sets[[de_key]]),
      n_overlap           = n_overlap,
      pct_dcopa_in_de     = round(100 * n_overlap / n_dcopa, 1),
      stringsAsFactors    = FALSE
    )
  }))
}))

print(results, row.names = FALSE)

out_path <- file.path(SAVE_DIR, "dcopa_vs_DE_overlap.csv")
write.csv(results, out_path, row.names = FALSE)
message("Saved: ", out_path)

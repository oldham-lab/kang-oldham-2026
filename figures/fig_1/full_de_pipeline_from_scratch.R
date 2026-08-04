# slightly different prompt where I asked to rewrite this script:
# /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_7/DE_old/full_DE_pipeline_ADvsCon_nodonorblocking_v3.R
# except for celltype vs all, with donor blocking, control samples only
# probably don't need to use this

library(qs)
library(data.table)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(future.apply)
library(BiocParallel)

# ============================================================
# PARALLELISM
# ============================================================
N_WORKERS      <- min(parallelly::availableCores() - 1L, 8L)
plan(multisession, workers = N_WORKERS)
BPPARAM_SERIAL <- BiocParallel::SerialParam()


# ============================================================
# DATASET DEFINITIONS
# ============================================================
# Required fields per dataset:
#   expr_path     — CSV/CSV.gz, genes x samples, first column = "Gene"
#   anno_path     — CSV, one row per sample
#   celltype_col  — celltype column name in anno
#   donor_col     — donor column name in anno
#   de_save_dir   — output directory
#
# Optional (disease filtering):
#   dx_col        — diagnosis column. NULL = use all samples
#   dx_control    — value indicating controls. NULL = use all samples

DATASETS <- list(

  MIT_Multiome = list(
    expr_path    = file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), ".../mit_pfc_pseudobulk.csv.gz"),
    anno_path    = file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), ".../mit_pfc_pseudobulk_annotations.csv"),
    celltype_col = "Cell_Type",
    donor_col    = "Donor_ID",
    dx_col       = "Diagnosis",
    dx_control   = "Control",
    de_save_dir  = file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_celltype_vs_all/")
  ),

  SEA_AD = list(
    expr_path    = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/cell_expression_by_donor_sum_subclass.csv"),
    anno_path    = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/annotations_by_donor_sum_subclass.csv"),
    celltype_col = "Subclass",
    donor_col    = "Donor",
    dx_col       = "Dx",
    dx_control   = "Control",
    de_save_dir  = file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_celltype_vs_all/")
  )

  # No-disease dataset example:
  # My_Dataset = list(
  #   expr_path    = "/path/to/pseudobulk.csv.gz",
  #   anno_path    = "/path/to/annotations.csv",
  #   celltype_col = "celltype",
  #   donor_col    = "donor",
  #   dx_col       = NULL,
  #   dx_control   = NULL,
  #   de_save_dir  = "/path/to/results/"
  # )
)

MIN_SAMPLES_PER_GROUP <- 3L


# ============================================================
# HELPERS
# ============================================================

load_expr <- function(path) {
  message("  Loading: ", path)
  mat <- fread(path, data.table = FALSE) |>
    tibble::column_to_rownames("Gene") |>
    as.matrix() |>
    round()
  storage.mode(mat) <- "integer"
  mat
}

filter_to_controls <- function(mat, anno, dx_col, dx_control) {
  if (is.null(dx_col) || is.null(dx_control)) {
    message("  No disease column — using all ", ncol(mat), " samples")
    return(list(mat = mat, anno = anno))
  }
  if (!dx_col %in% colnames(anno)) {
    message("  WARNING: '", dx_col, "' not found in annotation — using all samples")
    return(list(mat = mat, anno = anno))
  }
  keep <- anno[[dx_col]] == dx_control
  if (!any(keep))
    stop("No samples match dx_control = '", dx_control, "' in '", dx_col, "'")
  if (all(keep)) {
    message("  All samples are controls — no filtering needed")
    return(list(mat = mat, anno = anno))
  }
  anno <- anno[keep, , drop = FALSE]
  mat  <- mat[, keep, drop = FALSE]
  message("  Filtered to controls: ", ncol(mat), " samples retained")
  list(mat = mat, anno = anno)
}

check_feasibility <- function(is_target, donor_vec, ct) {
  n_target <- sum(is_target)
  n_other  <- sum(!is_target)
  if (n_target < MIN_SAMPLES_PER_GROUP)
    return(paste0(ct, ": only ", n_target, " target samples"))
  if (n_other < MIN_SAMPLES_PER_GROUP)
    return(paste0(ct, ": only ", n_other, " other samples"))
  shared <- intersect(unique(donor_vec[is_target]), unique(donor_vec[!is_target]))
  if (length(shared) == 0)
    return(paste0(ct, ": no donor in both groups — blocking not possible"))
  NULL
}


# ============================================================
# EDGER: celltype-vs-all
# ============================================================
run_edger_ct_vs_all <- function(counts_mat, anno, celltype_col,
                                 donor_col, de_save_dir, save_prefix = "") {
  if (!dir.exists(de_save_dir)) dir.create(de_save_dir, recursive = TRUE)
  celltypes <- unique(anno[[celltype_col]])
  message("  edgeR: ", length(celltypes), " celltypes")

  result_list <- future_lapply(celltypes, function(ct) {
    is_target <- anno[[celltype_col]] == ct
    donor_f   <- factor(anno[[donor_col]])
    target_f  <- factor(ifelse(is_target, "target", "other"),
                        levels = c("other", "target"))

    reason <- check_feasibility(is_target, anno[[donor_col]], ct)
    if (!is.null(reason)) { message("    Skipping: ", reason); return(NULL) }

    keep       <- filterByExpr(counts_mat, group = target_f,
                               min.count = 5, min.total.count = 10)
    sub_counts <- counts_mat[keep, ]
    design     <- model.matrix(~ donor_f + target_f)

    if (qr(design)$rank < ncol(design)) {
      message("    Skipping ", ct, ": design not full rank"); return(NULL)
    }

    y   <- DGEList(counts = sub_counts, genes = rownames(sub_counts))
    y   <- calcNormFactors(y, method = "TMM")
    y   <- estimateDisp(y, design, robust = TRUE)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = ncol(design))

    res          <- as.data.frame(topTags(lrt, n = nrow(sub_counts)))
    res$celltype <- ct
    res
  }, future.seed = TRUE)

  names(result_list) <- celltypes
  result_list        <- Filter(Negate(is.null), result_list)
  out <- file.path(de_save_dir, paste0(save_prefix, "edgeR_ct_vs_all.RDS"))
  saveRDS(result_list, file = out)
  message("  Saved: ", out)
  invisible(result_list)
}


# ============================================================
# DESEQ2: celltype-vs-all
# ============================================================
run_deseq2_ct_vs_all <- function(counts_mat, anno, celltype_col,
                                  donor_col, de_save_dir, save_prefix = "") {
  if (!dir.exists(de_save_dir)) dir.create(de_save_dir, recursive = TRUE)
  celltypes <- unique(anno[[celltype_col]])
  message("  DESeq2: ", length(celltypes), " celltypes")

  result_list <- future_lapply(celltypes, function(ct) {
    is_target <- anno[[celltype_col]] == ct

    reason <- check_feasibility(is_target, anno[[donor_col]], ct)
    if (!is.null(reason)) { message("    Skipping: ", reason); return(NULL) }

    col_data <- data.frame(
      donor  = factor(anno[[donor_col]]),
      target = factor(ifelse(is_target, "target", "other"),
                      levels = c("other", "target"))
    )
    keep       <- rowSums(counts_mat >= 5) >= max(2L, floor(ncol(counts_mat) * 0.2))
    sub_counts <- counts_mat[keep, ]

    dds  <- DESeqDataSetFromMatrix(sub_counts, col_data, design = ~ donor + target)
    dds1 <- DESeq(dds, test = "LRT", reduced = ~ donor,
                  BPPARAM = BPPARAM_SERIAL, quiet = TRUE)

    res   <- results(dds1, name = "target_target_vs_other")
    resdf <- as.data.frame(res@listData) |>
      `rownames<-`(res@rownames) |>
      arrange(padj)
    resdf$celltype <- ct
    resdf$gene     <- rownames(resdf)
    colnames(resdf)[c(2, 5, 6)] <- c("logFC", "PValue", "FDR")
    resdf
  }, future.seed = TRUE)

  names(result_list) <- celltypes
  result_list        <- Filter(Negate(is.null), result_list)
  out <- file.path(de_save_dir, paste0(save_prefix, "DESeq2_ct_vs_all.qs"))
  qsave(result_list, file = out)
  message("  Saved: ", out)
  invisible(result_list)
}


# ============================================================
# MAIN
# ============================================================
all_results <- lapply(names(DATASETS), function(ds_name) {
  cfg <- DATASETS[[ds_name]]
  message("\n====== Dataset: ", ds_name, " ======")

  mat  <- load_expr(cfg$expr_path)
  anno <- fread(cfg$anno_path, data.table = FALSE)

  # Align annotation to matrix column order
  if (!identical(colnames(mat), anno[[cfg$donor_col]])) {
    idx <- match(colnames(mat), anno[[cfg$donor_col]])
    if (anyNA(idx))
      stop(ds_name, ": ", sum(is.na(idx)),
           " matrix columns unmatched in annotation via '", cfg$donor_col, "'")
    anno <- anno[idx, , drop = FALSE]
  }

  filtered <- filter_to_controls(mat, anno, cfg$dx_col, cfg$dx_control)
  mat  <- filtered$mat
  anno <- filtered$anno

  message("  Samples: ", ncol(mat), " | Celltypes: ",
          length(unique(anno[[cfg$celltype_col]])), " | Genes: ", nrow(mat))

  edger_res <- run_edger_ct_vs_all(mat, anno, cfg$celltype_col,
                                    cfg$donor_col, cfg$de_save_dir,
                                    paste0(ds_name, "_"))
  deseq_res <- run_deseq2_ct_vs_all(mat, anno, cfg$celltype_col,
                                     cfg$donor_col, cfg$de_save_dir,
                                     paste0(ds_name, "_"))
  list(edgeR = edger_res, DESeq2 = deseq_res)
})

names(all_results) <- names(DATASETS)
message("\n====== All datasets complete ======")
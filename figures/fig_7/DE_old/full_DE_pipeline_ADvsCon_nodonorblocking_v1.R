library(qs)
library(data.table)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(future.apply)
library(BiocParallel)
options(future.globals.maxSize = 3000*1024^2)

# ============================================================
# PLACEHOLDERS — replace all caps variables before running
# ============================================================

# --- MIT_Multiome file paths ---
# These are the outputs from running prepare_mit_multiome.py WITHOUT --no-pseudobulk.
# Pseudobulking is done in Python where 64-bit sparse ops are fully supported;
# R receives a small dense genes x samples matrix directly.
#
# To generate:
#   python prepare_mit_multiome.py \
#       --input  /path/to/mit_multiome.h5ad \
#       --outdir /home/gugene/bdata/@shared/.../PFC/pfc_data/
#
MIT_PB_EXPR_PATH <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/PFC/pfc_data_pseudobulk/mit_pfc_pseudobulk.csv.gz"
MIT_PB_ANNO_PATH <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/PFC/pfc_data_pseudobulk/mit_pfc_pseudobulk_annotations.csv"

# --- MIT_Multiome column names in annotation file ---
MIT_CELLTYPE_COL <- "RNA.Subclass"    # celltype label column
MIT_DX_COL       <- "Diagnosis2"    # diagnosis column
MIT_DX_CONTROL   <- "Con"      # value in diagnosis column for controls
MIT_DX_AD        <- "AD"           # value in diagnosis column for AD cases

# --- SEA-AD file paths (already pseudobulked) ---
SEA_EXPR_PATH  <- file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/cell_expression_by_donor_sum_subclass.csv")
SEA_ANNO_PATH  <- file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/annotations_by_donor_sum_subclass.csv")

# --- SEA-AD column names in annotation file (pseudobulk-level) ---
SEA_CELLTYPE_COL <- "Subclass"
SEA_DX_COL       <- "Dx"
SEA_DX_CONTROL   <- "Control"
SEA_DX_AD        <- "Alzheimers"

# --- Output directories ---
MIT_DATA_DIR   <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_pseudobulk_output")
MIT_DE_DIR     <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_pseudobulk_output")
SEA_DE_DIR     <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/")

# ============================================================
# Parallelism setup
# ============================================================
N_WORKERS <- min(parallelly::availableCores() - 1L, 8L)
plan(multisession, workers = N_WORKERS)
BPPARAM_SERIAL <- BiocParallel::SerialParam()


# ============================================================
# STEP 1: Load MIT_Multiome pseudobulked data
# ============================================================
# Pseudobulking was done in Python (prepare_mit_multiome.py) where 64-bit
# sparse operations are fully supported. At 2.6B+ non-zeros, R cannot load
# the raw sparse matrix regardless of format — the triplet vectors exceed
# R's internal vector length limit. Pseudobulking in Python reduces the
# data to a small genes x samples dense matrix before R ever touches it.

message("--- Loading MIT_Multiome pseudobulked data ---")

mit_rawcounts <- fread(MIT_PB_EXPR_PATH, data.table = FALSE) |>
  tibble::column_to_rownames(var = "Gene")

mit_anno_pb <- fread(MIT_PB_ANNO_PATH, data.table = FALSE) |>
  dplyr::mutate(
    Dx       = factor(.data[[MIT_DX_COL]], levels = c(MIT_DX_CONTROL, MIT_DX_AD)),
    celltype = factor(.data[[MIT_CELLTYPE_COL]])
  )

mit_rawcounts <- data.frame(Gene = rownames(mit_rawcounts), mit_rawcounts)

message("MIT_Multiome loaded: ",
        nrow(mit_anno_pb), " pseudobulk samples, ",
        length(unique(mit_anno_pb$celltype)), " celltypes")


# ============================================================
# STEP 2: Load SEA-AD pseudobulked data
# ============================================================
message("--- Loading SEA-AD ---")

sea_rawcounts <- fread(SEA_EXPR_PATH, data.table = FALSE) |>
  tibble::column_to_rownames(var = "Gene")

sea_anno_pb <- fread(SEA_ANNO_PATH, data.table = FALSE) |>
  mutate(
    Dx       = factor(.data[[SEA_DX_COL]], levels = c(SEA_DX_CONTROL, SEA_DX_AD)),
    celltype = factor(.data[[SEA_CELLTYPE_COL]])
  )


# ============================================================
# STEP 3: Align gene universes, process celltypes independently
# ============================================================
# Genes are intersected so both datasets use the same gene set,
# but celltypes are processed separately per dataset.
common_genes      <- intersect(mit_rawcounts$Gene, rownames(sea_rawcounts))
mit_rawcounts_mat <- mit_rawcounts[mit_rawcounts$Gene %in% common_genes, ]
sea_rawcounts_mat <- data.frame(Gene = rownames(sea_rawcounts),
                                sea_rawcounts)[rownames(sea_rawcounts) %in% common_genes, ]

message(length(common_genes), " genes in common between datasets")


# ============================================================
# STEP 4: edgeR AD vs Control DE — one function, two datasets
# ============================================================
# For AD vs Control we run DE *within* each celltype using a simple ~ Dx
# design (no donor blocking — each donor is exclusively AD or Control,
# so donor and Dx are the same variable).

run_edger_dx_de <- function(rawcounts_df,   # data.frame: col1 = Gene, rest = samples
                             anno,           # data.frame: one row per sample
                             celltype_col,   # column name for celltype in anno
                             dx_col,         # column name for diagnosis in anno
                             dx_control,     # value for control group in dx_col
                             dx_ad,          # value for AD group in dx_col
                             de_save_dir,
                             save_prefix = "") {

  if (!dir.exists(de_save_dir)) dir.create(de_save_dir, recursive = TRUE)

  genes_vec  <- rawcounts_df[, 1]
  counts_mat <- as.matrix(rawcounts_df[, -1])

  cts <- unique(anno[[celltype_col]])
  result_list <- future_lapply(cts, function(ct) {

    these      <- which(anno[[celltype_col]] == ct)
    sub_counts <- counts_mat[, these]
    sub_anno   <- anno[these, ]
    dx_fac     <- factor(sub_anno[[dx_col]], levels = c(dx_control, dx_ad))

    # Require at least 3 samples per group
    if (sum(dx_fac == dx_control) < 3 || sum(dx_fac == dx_ad) < 3) {
      message("  Skipping ", ct, ": fewer than 3 samples in one group")
      return(NULL)
    }

    # Filter lowly expressed genes within this celltype
    keep       <- filterByExpr(sub_counts, group = dx_fac,
                               min.count = 5, min.total.count = 10)
    sub_counts <- sub_counts[keep, ]
    sub_genes  <- genes_vec[keep]

    design <- model.matrix(~ dx_fac)

    y   <- DGEList(counts = sub_counts, genes = sub_genes)
    y   <- calcNormFactors(y, method = "TMM")
    y   <- estimateDisp(y, design, robust = TRUE)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)   # coef 2 = dx_facAD

    res <- as.data.frame(topTags(lrt, n = nrow(sub_counts)))
    res$celltype <- ct
    res

  }, future.seed = TRUE)

  names(result_list) <- cts
  result_list <- Filter(Negate(is.null), result_list)

  saveRDS(result_list,
          file = file.path(de_save_dir, paste0(save_prefix, "edgeR_ADvsCon_by_celltype.RDS")))
  message(save_prefix, " edgeR DE done")
  invisible(result_list)
}


# ============================================================
# STEP 5: DESeq2 AD vs Control DE — one function, two datasets
# ============================================================

run_deseq2_dx_de <- function(counts_mat,     # matrix or data.frame (genes x samples)
                              anno,           # data.frame: one row per sample
                              celltype_col,
                              dx_col,
                              dx_control,     # value for control group in dx_col
                              dx_ad,          # value for AD group in dx_col
                              de_save_dir,
                              save_prefix = "") {

  if (!dir.exists(de_save_dir)) dir.create(de_save_dir, recursive = TRUE)

  if (is.data.frame(counts_mat)) {
    counts_mat <- as.matrix(counts_mat[, -1]) |>
      `rownames<-`(counts_mat[, 1])
  }

  cts <- unique(anno[[celltype_col]])
  result_list <- future_lapply(cts, function(ct) {

    these      <- which(anno[[celltype_col]] == ct)
    sub_counts <- counts_mat[, these]
    sub_anno   <- anno[these, ]
    sub_anno$Dx <- factor(sub_anno[[dx_col]], levels = c(dx_control, dx_ad))

    # Require at least 3 samples per group
    if (sum(sub_anno$Dx == dx_control) < 3 || sum(sub_anno$Dx == dx_ad) < 3) {
      message("  Skipping ", ct, ": fewer than 3 samples in one group")
      return(NULL)
    }

    # Basic low-count pre-filter
    keep       <- rowSums(sub_counts >= 5) >= max(2, floor(ncol(sub_counts) * 0.2))
    sub_counts <- sub_counts[keep, ]

    dds  <- DESeqDataSetFromMatrix(
      countData = sub_counts,
      colData   = sub_anno,
      design    = ~ Dx
    )
    dds1 <- DESeq(dds,
                  test    = "LRT",
                  reduced = ~ 1,
                  BPPARAM = BPPARAM_SERIAL,
                  quiet   = TRUE)

    res_name <- paste0("Dx_", dx_ad, "_vs_", dx_control)
    res   <- results(dds1, name = res_name)
    resdf <- as.data.frame(res@listData) |>
      `rownames<-`(res@rownames) |>
      arrange(padj)
    resdf$celltype <- ct
    resdf$genes    <- rownames(resdf)
    colnames(resdf)[c(2, 5, 6)] <- c("logFC", "PValue", "FDR")
    resdf

  }, future.seed = TRUE)

  names(result_list) <- cts
  result_list <- Filter(Negate(is.null), result_list)

  qsave(result_list,
        file = file.path(de_save_dir, paste0(save_prefix, "DESeq2_ADvsCon_by_celltype.qs")))
  message(save_prefix, " DESeq2 DE done")
  invisible(result_list)
}


# ============================================================
# STEP 6: Run DE on both datasets
# ============================================================

message("--- Running edgeR: MIT_Multiome ---")
mit_edger <- run_edger_dx_de(
  rawcounts_df  = mit_rawcounts_mat,
  anno          = mit_anno_pb,
  celltype_col  = "celltype",
  dx_col        = "Dx",
  dx_control    = MIT_DX_CONTROL,
  dx_ad         = MIT_DX_AD,
  de_save_dir   = MIT_DE_DIR,
  save_prefix   = "mit_"
)

message("--- Running edgeR: SEA-AD ---")
sea_edger <- run_edger_dx_de(
  rawcounts_df  = sea_rawcounts_mat,
  anno          = sea_anno_pb,
  celltype_col  = "celltype",
  dx_col        = "Dx",
  dx_control    = SEA_DX_CONTROL,
  dx_ad         = SEA_DX_AD,
  de_save_dir   = SEA_DE_DIR,
  save_prefix   = "sea_"
)

message("--- Running DESeq2: MIT_Multiome ---")
mit_deseq2 <- run_deseq2_dx_de(
  counts_mat   = mit_rawcounts_mat,
  anno         = mit_anno_pb,
  celltype_col = "celltype",
  dx_col       = "Dx",
  dx_control   = MIT_DX_CONTROL,
  dx_ad        = MIT_DX_AD,
  de_save_dir  = MIT_DE_DIR,
  save_prefix  = "mit_"
)

message("--- Running DESeq2: SEA-AD ---")
sea_deseq2 <- run_deseq2_dx_de(
  counts_mat   = sea_rawcounts_mat,
  anno         = sea_anno_pb,
  celltype_col = "celltype",
  dx_col       = "Dx",
  dx_control   = SEA_DX_CONTROL,
  dx_ad        = SEA_DX_AD,
  de_save_dir  = SEA_DE_DIR,
  save_prefix  = "sea_"
)

message("--- All DE complete ---")

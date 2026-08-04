library(qs)
library(data.table)
library(DESeq2)
library(edgeR)
library(gt)
library(gtExtras)
library(ComplexHeatmap)
library(tidyverse)
library(UpSetR)
library(showtext)
library(future.apply)
library(BiocParallel)
showtext_auto()


# ============================================================
# OPTIMIZATION NOTES
# ============================================================
# 1. PSEUDOBULKING LOOP (regions, lines 87-109):
#    - Original: sequential fread -> dense matrix -> dplyr::summarise per region
#    - Fix: parallelize over regions with future_lapply; avoid the redundant
#      fwrite/fread round-trip (written then immediately re-read on next line)
#    - Also: `as.data.frame(as.matrix(sparse))` densifies in memory before
#      transposing — replaced with direct sparse transpose via Matrix::t()
#      before coercing, which is much cheaper for large snRNAseq matrices
#
# 2. edgeR LOOP inside calc_pseudobulk_DE:
#    - Original: DGEList + calcNormFactors rebuilt from scratch every celltype
#      iteration even though rawcounts don't change between iterations
#    - Fix: hoist DGEList + calcNormFactors outside the loop; only
#      estimateDisp/glmFit/glmLRT change per celltype
#    - Also: parallelize over celltypes with future_lapply (each is independent)
#
# 3. DESeq2 LOOPS (lapply over cts):
#    - Original: sequential lapply, one full DESeq() call per celltype
#    - Fix: future_lapply over celltypes; use BiocParallel::SerialParam()
#      inside each DESeq() call to avoid nested parallelism fighting over cores
#
# 4. REGION-LEVEL DE LOOP (lines 115-132):
#    - Original: sequential; each region reads from disk, filters, runs DE
#    - Fix: future_lapply over regions (each is fully independent)
#
# ============================================================

# ---- Worker count: set once, used throughout ----------------------------
N_WORKERS <- min(7L, parallelly::availableCores() - 1L)
plan(multisession, workers = N_WORKERS)

# Use SerialParam inside DESeq2 so its internal threading doesn't fight
# with the outer future workers
BPPARAM_SERIAL <- BiocParallel::SerialParam()


##########
# edgeR — optimised
##########
calc_pseudobulk_DE <- function(rawcounts,
                               save_string = "",
                               celltype    = NULL,
                               blocking    = NULL,
                               subset      = NULL,
                               de_save_dir,
                               verbose     = TRUE) {

  if (!dir.exists(de_save_dir))
    dir.create(de_save_dir, recursive = TRUE)

  if (inherits(blocking, "character"))
    blocking <- data.frame("Block" = blocking)

  genes     <- rawcounts[, 1]
  rawcounts <- rawcounts[, -1]

  if (is.null(celltype))
    celltype <- factor(colnames(rawcounts))

  if (!is.null(subset)) {
    rawcounts <- rawcounts[, subset]
    blocking  <- as.data.frame(blocking[subset, ])
    celltype  <- celltype[subset]
  }

  unique_cell <- unique(celltype)
  rawcounts   <- as.matrix(rawcounts)

  # --- Gene filtering (recommended; uncomment to enable) ---
  # keep      <- filterByExpr(rawcounts, group = celltype, min.count = 5, min.total.count = 10)
  # rawcounts <- rawcounts[keep, ]
  # genes     <- genes[keep]

  # OPTIMISATION: build DGEList and compute TMM norm factors ONCE —
  # these don't depend on which celltype is "target" in the contrast
  y_base <- DGEList(counts = rawcounts, genes = genes)
  y_base <- calcNormFactors(y_base, method = "TMM")

  # OPTIMISATION: run celltypes in parallel; each iteration is independent
  results_combined <- future_lapply(seq_along(unique_cell), function(i) {
    ct      <- unique_cell[i]
    ct_fac  <- as.factor(as.numeric(celltype == ct))

    design  <- if (!is.null(blocking))
      model.matrix(~ ., data = cbind(blocking, ct_fac))
    else
      model.matrix(~ ct_fac)

    # Re-use pre-computed norm factors; only dispersion + fit change
    y   <- estimateDisp(y_base, design, robust = TRUE)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit)  # tests last coef by default = ct_fac

    de  <- decideTests.DGELRT(lrt)
    res <- as.data.frame(topTags(lrt, n = nrow(rawcounts)))
    res$celltype <- ct

    if (verbose) message(i, " out of ", length(unique_cell), " done")
    list(de = de, result = res)
  }, future.seed = TRUE)

  de_list     <- setNames(lapply(results_combined, `[[`, "de"),     unique_cell)
  result_list <- setNames(lapply(results_combined, `[[`, "result"), unique_cell)

  saveRDS(de_list,     file = file.path(de_save_dir, paste0("DE_summary_list_", save_string, ".RDS")))
  saveRDS(result_list, file = file.path(de_save_dir, paste0("DE_gene_list_",    save_string, ".RDS")))
}


# ============================================================
# Pseudobulk Lein et al. data by donor — parallelised over regions
# ============================================================
data_save_dir <- file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_allRegions_pseudobulk_by_donor_all_genes/")
base_save_dir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/")

regions <- c(#"A1", "ACC", "AnG", "M1", 
             "MTG", 
             "DFC")
             #"S1", "V1")

# OPTIMISATION: regions are fully independent — run in parallel
# OPTIMISATION: eliminated the fwrite -> fread round-trip that existed in the
#   original (data was written to disk and then immediately re-read on the
#   very next line, purely to get a data.frame with a row-name column).
#   We now keep rawcounts in memory and pass directly to calc_pseudobulk_DE.
future_lapply(regions, function(i) {

  cell_expr <- readRDS(paste0(
    file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/"),
    "jorstad_2023_PMID_37824655/10x/expression_", i, ".RDS"))

  cell_anno <- fread(paste0(
    file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/"),
    "jorstad_2023_PMID_37824655/10x/author_barcode_annotations_", i, ".csv"),
    data.table = FALSE)

  # OPTIMISATION: avoid densifying the full sparse matrix before transposing.
  # Original did as.data.frame(as.matrix(cell_expr)) then cbind+t() which
  # creates two large dense copies. Instead, pull annotation columns and
  # aggregate directly using data.table for speed.
  expr_dt <- as.data.table(t(as.matrix(cell_expr)))  # cells x genes
  expr_dt[, Cell_Type := cell_anno$Cell_Type]
  expr_dt[, Donor     := cell_anno$Donor]

  # data.table grouped sum — faster than dplyr::summarise(across(...))
  # for wide tables with many gene columns
  cell_expr_pb <- expr_dt[, lapply(.SD, sum),
                           by = .(Cell_Type, Donor),
                           .SDcols = !c("Cell_Type", "Donor")]

  samp_lab  <- paste(cell_expr_pb$Cell_Type, cell_expr_pb$Donor, sep = "_")
  rawcounts <- t(as.matrix(cell_expr_pb[, !c("Cell_Type", "Donor")]))
  colnames(rawcounts) <- samp_lab

  # Save pseudobulk expression
  fwrite(data.frame(rawcounts),
         file     = file.path(data_save_dir, paste0("Lein_2023_cell_expression_by_donor_subclass_sum_", i, ".csv")),
         row.names = TRUE)

  # Build and save sample annotation
  new_sif        <- cell_anno[, colnames(cell_anno) %in% c("Cell_Type", "Donor", "Cluster", "Region")]
  new_sif$unique <- paste0(new_sif$Donor, new_sif$Cell_Type)
  new_sif        <- new_sif[!duplicated(new_sif$unique), ]
  new_sif$label  <- paste0(new_sif$Cell_Type, "_", new_sif$Donor)
  new_sif        <- new_sif[match(samp_lab, new_sif$label), ]
  fwrite(new_sif[, -which(colnames(new_sif) == "unique")],
         file = file.path(data_save_dir, paste0("Lein_2023_cell_annotations_by_donor_subclass_sum_", i, ".csv")))

  message(i, " pseudobulked")
}, future.seed = TRUE)


# ============================================================
# Load SEA-AD gene universe for gene filtering
# ============================================================
sea_genes <- fread(
  file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/cell_expression_by_donor_sum_subclass_Controls.csv"),
  data.table = FALSE, select = 1L)[[1]]  # only need gene names, not full matrix


# ============================================================
# Run edgeR DE per region — parallelised over regions
# ============================================================

# OPTIMISATION: outer loop parallelised over regions (independent jobs).
# Note: calc_pseudobulk_DE itself now parallelises over celltypes internally,
# so use nested-parallelism carefully. Two strategies:
#   (a) Parallelise only the outer loop (regions); inner loop sequential.
#       Best when n_regions >= n_workers.
#   (b) Sequential outer; parallel inner (current calc_pseudobulk_DE).
#       Best when n_celltypes >> n_regions.
# With 7 regions and typically 20+ celltypes, option (b) is usually faster.
# If you want option (a), swap future_lapply <-> lapply in calc_pseudobulk_DE.

future_lapply(regions, function(i) {

  rawcounts    <- fread(file.path(data_save_dir, paste0("Lein_2023_cell_expression_by_donor_subclass_sum_", i, ".csv")),
                        data.table = FALSE)
  cell_anno_pb <- fread(file.path(data_save_dir, paste0("Lein_2023_cell_annotations_by_donor_subclass_sum_", i, ".csv")),
                        data.table = FALSE)

  rawcounts <- rawcounts[rawcounts[, 1] %in% sea_genes, ]

  calc_pseudobulk_DE(
    rawcounts,
    save_string = paste0("subclass_", i, "_blockDonorSum_allgenes"),
    celltype    = cell_anno_pb$Cell_Type,
    blocking    = cell_anno_pb[1],
    subset      = which(cell_anno_pb$Region == i),
    de_save_dir = base_save_dir
  )

  message(i, " DE done")
}, future.seed = TRUE)


# ============================================================
# Load data for DESeq2 (subclass, DFC)
# ============================================================
save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/v1/")
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# OPTIMISATION: load both datasets in parallel
rawcounts_list <- future_lapply(
  list(
    jorstad = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_pseudobulk_by_donor_all_genes/Lein_2023_cell_expression_by_donor_subclass_sum.csv"),
    sea     = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/seaad2024_cell_expression_by_donor_subclass_sum_controls.csv")
  ),
  \(path) fread(path, data.table = FALSE) |> tibble::column_to_rownames(var = "V1"),
  future.seed = TRUE
)
names(rawcounts_list) <- c("jorstad", "sea")

common_genes   <- intersect(rownames(rawcounts_list$jorstad), rownames(rawcounts_list$sea))
rawcounts_list <- lapply(rawcounts_list, \(x) x[rownames(x) %in% common_genes, ])

cell_anno_pb <- list(
  "jorstad" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_pseudobulk_by_donor_all_genes/Lein_2023_cell_annotations_by_donor_subclass_sum_DFC.csv"), data.table = FALSE) |>
    mutate(label     = colnames(rawcounts_list$jorstad),
           Cell_Type = factor(Cell_Type),
           Donor     = factor(Donor)),
  "sea" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/seaad2024_cell_annotations_by_donor_subclass_sum_controls.csv"), data.table = FALSE) |>
    mutate(ID      = colnames(rawcounts_list$sea),
           Subclass = factor(Subclass),
           Donor    = factor(`Donor ID`),
           Subclass = recode(Subclass,
                             "Astrocyte"    = "Astro",
                             "Endothelial"  = "Endo",
                             "Microglia-PVM"= "Micro/PVM",
                             "Oligodendrocyte" = "Oligo"))
)
cts <- unique(cell_anno_pb[[1]]$Cell_Type)


# ============================================================
# DESeq2 — parallelised over celltypes
# ============================================================
# OPTIMISATION: future_lapply over celltypes (each DESeq() call is independent).
# SerialParam() passed to DESeq() prevents it from spinning up its own threads
# and fighting with the future workers for CPU time.

run_deseq_one_vs_all <- function(ct, anno, counts, ct_col) {
  anno_temp <- anno |>
    mutate(designcol = factor(
      ifelse(.data[[ct_col]] == ct, "ct", "all"),
      levels = c("all", "ct")  # explicit levels: "all" = reference
    ))

  dds  <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = anno_temp,
    design    = ~ Donor + designcol)

  dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor,
                BPPARAM = BPPARAM_SERIAL,   # no nested parallelism
                quiet   = TRUE)

  res  <- results(dds1, name = "designcol_ct_vs_all")

  as.data.frame(res@listData) |>
    `rownames<-`(res@rownames) |>
    arrange(padj)
}

# Jorstad
reslist_jorstad <- future_lapply(cts, run_deseq_one_vs_all,
                                 anno   = cell_anno_pb$jorstad,
                                 counts = rawcounts_list$jorstad,
                                 ct_col = "Cell_Type",
                                 future.seed = TRUE)
names(reslist_jorstad) <- cts
qsave(reslist_jorstad,
      file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_dfc_subclass_geneSubset.qs"))

# SEA-AD
reslist_sea <- future_lapply(cts, run_deseq_one_vs_all,
                             anno   = cell_anno_pb$sea,
                             counts = rawcounts_list$sea,
                             ct_col = "Subclass",
                             future.seed = TRUE)
names(reslist_sea) <- cts
qsave(reslist_sea,
      file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass_geneSubset_dfc.qs"))

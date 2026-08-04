# 2.2 did not run DESeq2 MTG - 2.3 fixes this
# 2.4: edgeR gene matching now uses the PER-REGION common gene universe
#      (Jorstad ∩ Gabitto within each region), matching the DESeq2 logic.
#      DFC uses Jorstad-DFC ∩ SEAAD-DFC; MTG uses Jorstad-MTG ∩ SEAAD-MTG.
#      Regions are kept separate — no matching across all four.

library(qs)
library(data.table)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(future.apply)
library(BiocParallel)


# ============================================================
# OPTIMIZATION NOTES
# ============================================================
# 1. PSEUDOBULKING LOOP (regions):
#    - Parallelised over regions with future_lapply; grouped sum via data.table
#    - CAVEAT: still densifies via as.matrix() before the grouped sum. A sparse
#      sum-by-group (cell->group indicator matrix-multiply, or
#      Matrix.utils::aggregate.Matrix) would avoid this but is not implemented.
#
# 2. edgeR LOOP inside calc_pseudobulk_DE:
#    - Original: DGEList + calcNormFactors rebuilt from scratch every celltype
#      iteration even though rawcounts don't change between iterations
#    - Fix: hoist DGEList + calcNormFactors outside the loop; only
#      estimateDisp/glmFit/glmLRT change per celltype
#    - Also: parallelize over celltypes with future_lapply (each is independent)
#
# 3. DESeq2 / edgeR per-celltype loops:
#    - future_lapply over celltypes (each fit is independent); SerialParam()
#      inside each DESeq() call so its internal threads don't fight the workers
#
# 4. REGION-LEVEL loops (STEP 3 & STEP 5):
#    - Regions run SEQUENTIALLY (plain lapply) on purpose: future does not nest
#      parallelism, so keeping the outer loop serial lets the inner per-celltype
#      future_lapply use all workers (24 celltypes >> 2 regions).
#      STEP 1 (pseudobulk) stays future_lapply over regions — it has no inner
#      futures, so the 2 regions pseudobulk concurrently.
#
# ============================================================


# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================

# --- Jorstad (Lein 2023) --------------------------------------------------
JORSTAD_EXPR_DIR <- file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x")
JORSTAD_PB_DIR   <- file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_allRegions_pseudobulk_by_donor_all_genes")

# --- SEAAD 2024 (DFC) ----------------------------------------------------
SEAAD_EXPR_DFC   <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/gabitto_2024/expr_UMI_notADsamples.qs")
SEAAD_ANNO_DFC   <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")
SEAAD_PB_DIR     <- file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes")

# --- SEAAD 2024 (MTG) ----------------------------------------------------
SEAAD_EXPR_MTG   <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/gabitto_2024/expr_UMI_notADsamples_mtg.qs")
SEAAD_ANNO_MTG   <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")

# --- HGNC symbol harmonization -------------------------------------------
# Jorstad uses current HGNC symbols, SEAAD uses an older vintage. Genes renamed
# between vintages (AARS<->AARS1, MARCH1<->MARCHF1, histones, tRNA synthetases,
# ...) otherwise fail the symbol-string intersection below and get dropped from
# DE in BOTH datasets. We canonicalize both matrices to current HGNC before
# intersecting. Alias map: prev_symbol -> current_symbol (from HGNC complete set).
ALIAS_MAP_PATH <- Sys.getenv("HGNC_MAP", "/home/gugene/code/git/COPA/dataset_processing_unified/hgnc_prev_to_current.tsv")

# --- Regions to process ---------------------------------------------------
REGIONS <- c("MTG", "DFC")

# --- Output root (all results saved hierarchically under here) ------------
OUT_DIR <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/DE_v2.4")

# Subdirectories (derived from OUT_DIR — do not edit individually)
PB_DIR     <- file.path(OUT_DIR, "pseudobulk")
EDGER_DIR  <- file.path(OUT_DIR, "edgeR")
DESEQ2_DIR <- file.path(OUT_DIR, "DESeq2")

# Output files: naming convention [method]_[dataset]_[region]_subclass_matchedGenes.qs
# edgeR — Jorstad
EDGER_JORSTAD_DFC  <- file.path(EDGER_DIR,  "jorstad", "edgeR_jorstad_DFC_subclass_matchedGenes.qs")
EDGER_JORSTAD_MTG  <- file.path(EDGER_DIR,  "jorstad", "edgeR_jorstad_MTG_subclass_matchedGenes.qs")
# edgeR — SEAAD
EDGER_SEAAD_DFC    <- file.path(EDGER_DIR,  "seaad",   "edgeR_seaad_DFC_subclass_matchedGenes.qs")
EDGER_SEAAD_MTG    <- file.path(EDGER_DIR,  "seaad",   "edgeR_seaad_MTG_subclass_matchedGenes.qs")
# DESeq2 — Jorstad
DESEQ2_JORSTAD_DFC <- file.path(DESEQ2_DIR, "jorstad", "DESeq2_jorstad_DFC_subclass_matchedGenes.qs")
DESEQ2_JORSTAD_MTG <- file.path(DESEQ2_DIR, "jorstad", "DESeq2_jorstad_MTG_subclass_matchedGenes.qs")
# DESeq2 — SEAAD
DESEQ2_SEAAD_DFC   <- file.path(DESEQ2_DIR, "seaad",   "DESeq2_seaad_DFC_subclass_matchedGenes.qs")
DESEQ2_SEAAD_MTG   <- file.path(DESEQ2_DIR, "seaad",   "DESeq2_seaad_MTG_subclass_matchedGenes.qs")


# ============================================================
# PARALLELISM SETUP
# ============================================================
N_WORKERS <- min(7L, parallelly::availableCores() - 1L)
# multicore (fork) on Linux: workers share memory copy-on-write, so the large
# count matrices / ALIAS_VEC are not re-serialised to each worker and packages
# need not be reloaded. Forked workers also inherit the parent's attached
# packages. Falls back to sequential where fork is unavailable (e.g. Windows).
plan(multicore, workers = N_WORKERS)

# Use SerialParam inside DESeq2 so its internal threading doesn't fight
# with the outer future workers
BPPARAM_SERIAL <- BiocParallel::SerialParam()


# ============================================================
# HGNC symbol harmonization helpers
# ============================================================
# prev_symbol -> current_symbol lookup (unique prev symbols, so unambiguous)
.alias_tbl <- fread(ALIAS_MAP_PATH, data.table = FALSE)
ALIAS_VEC  <- setNames(.alias_tbl$current_symbol, .alias_tbl$prev_symbol)

# Map a vector of gene symbols to current HGNC; leave unmapped symbols as-is.
canonicalize_symbols <- function(x) {
  hit <- x %in% names(ALIAS_VEC)
  x[hit] <- ALIAS_VEC[x[hit]]
  x
}

# Read a pseudobulk CSV (col 1 "V1" = gene symbols, remaining cols = per-sample
# counts), canonicalize the gene symbols to current HGNC, and collapse any rows
# that collide after renaming by SUMMING their counts (few collisions: <=110 per
# dataset). Returns a data.frame with unique canonical symbols in column "V1" —
# uniqueness is required by DESeq2/edgeR rownames. Column order (samples) is
# preserved; row order may change after the group-sum.
read_pb_canonical <- function(path) {
  df    <- fread(path, data.table = FALSE)
  genes <- canonicalize_symbols(df[[1]])
  counts <- df[, -1, drop = FALSE]
  if (anyDuplicated(genes)) {
    dt <- as.data.table(counts)
    dt[, `__gene__` := genes]
    agg <- dt[, lapply(.SD, sum), by = `__gene__`,
              .SDcols = setdiff(names(dt), "__gene__")]
    out <- data.frame(V1 = agg[["__gene__"]],
                      agg[, !"__gene__"], check.names = FALSE)
  } else {
    out <- data.frame(V1 = genes, counts, check.names = FALSE)
  }
  out
}


# ============================================================
# edgeR DE function — optimised
# ============================================================
calc_pseudobulk_DE <- function(rawcounts,
                               save_path,        # full path to output .qs file
                               celltype  = NULL,
                               blocking  = NULL,
                               subset    = NULL,
                               verbose   = TRUE) {

  out_dir <- dirname(save_path)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

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
    ct     <- unique_cell[i]
    ct_fac <- as.factor(as.numeric(celltype == ct))

    design <- if (!is.null(blocking))
      model.matrix(~ ., data = cbind(blocking, ct_fac))
    else
      model.matrix(~ ct_fac)

    # Re-use pre-computed norm factors; only dispersion + fit change
    y   <- estimateDisp(y_base, design, robust = TRUE)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit)  # tests last coef by default = ct_fac

    # de  <- decideTests.DGELRT(lrt)   # DE summary (up/down/not-DE) — commented out
    res <- as.data.frame(topTags(lrt, n = nrow(rawcounts)))
    res$celltype <- ct

    if (verbose) message(i, " out of ", length(unique_cell), " done")
    res
  }, future.seed = TRUE)

  result_list <- setNames(results_combined, unique_cell)

  # de_list <- setNames(lapply(results_combined, `[[`, "de"), unique_cell)  # commented out
  # qsave(de_list, file = sub(".qs", "_summary.qs", save_path))             # commented out

  qsave(result_list, file = save_path)
}


# ============================================================
# STEP 1: Pseudobulk Jorstad (Lein 2023) data by donor
#         Parallelised over regions
# ============================================================
if (!dir.exists(PB_DIR)) dir.create(PB_DIR, recursive = TRUE)

future_lapply(REGIONS, function(i) {

  cell_expr <- readRDS(file.path(JORSTAD_EXPR_DIR, paste0("expression_", i, ".RDS")))

  cell_anno <- fread(
    file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x"),
              paste0("author_barcode_annotations_", i, ".csv")),
    data.table = FALSE)

  # OPTIMISATION: data.table grouped sum — faster than dplyr::summarise(across(...))
  expr_dt <- as.data.table(t(as.matrix(cell_expr)))  # cells x genes
  expr_dt[, Cell_Type := cell_anno$Cell_Type]
  expr_dt[, Donor     := cell_anno$Donor]

  cell_expr_pb <- expr_dt[, lapply(.SD, sum),
                           by = .(Cell_Type, Donor),
                           .SDcols = !c("Cell_Type", "Donor")]

  samp_lab  <- paste(cell_expr_pb$Cell_Type, cell_expr_pb$Donor, sep = "_")
  rawcounts <- t(as.matrix(cell_expr_pb[, !c("Cell_Type", "Donor")]))
  colnames(rawcounts) <- samp_lab

  fwrite(data.frame(rawcounts),
         file      = file.path(PB_DIR, paste0("jorstad_cell_expression_by_donor_subclass_sum_", i, ".csv")),
         row.names = TRUE)

  new_sif        <- cell_anno[, colnames(cell_anno) %in% c("Cell_Type", "Donor", "Cluster", "Region")]
  new_sif$unique <- paste0(new_sif$Donor, new_sif$Cell_Type)
  new_sif        <- new_sif[!duplicated(new_sif$unique), ]
  new_sif$label  <- paste0(new_sif$Cell_Type, "_", new_sif$Donor)
  new_sif        <- new_sif[match(samp_lab, new_sif$label), ]
  fwrite(new_sif[, -which(colnames(new_sif) == "unique")],
         file = file.path(PB_DIR, paste0("jorstad_cell_annotations_by_donor_subclass_sum_", i, ".csv")))

  message("Jorstad ", i, " pseudobulked")
}, future.seed = TRUE)


# ============================================================
# STEP 2: Pseudobulk SEAAD 2024 data by donor
#         Applied to both DFC and MTG using the same logic
# ============================================================
pseudobulk_seaad <- function(expr_path, anno_path, region) {

  cell_expr <- qread(expr_path)

  cell_anno <- fread(anno_path, data.table = FALSE) |>
    dplyr::filter(`Overall AD neuropathological Change` == "Not AD")
  cell_anno <- cell_anno[match(colnames(cell_expr), cell_anno$sample_id), ]

  # OPTIMISATION: data.table grouped sum (same approach as Jorstad pseudobulking)
  expr_dt <- as.data.table(t(as.matrix(cell_expr)))
  expr_dt[, Subclass   := cell_anno$Subclass]
  expr_dt[, `Donor ID` := cell_anno$`Donor ID`]

  cell_expr_pb <- expr_dt[, lapply(.SD, sum),
                           by = .(Subclass, `Donor ID`),
                           .SDcols = !c("Subclass", "Donor ID")]

  samp_lab  <- paste(cell_expr_pb$Subclass, cell_expr_pb$`Donor ID`, sep = "_")
  rawcounts <- t(as.matrix(cell_expr_pb[, !c("Subclass", "Donor ID")]))
  colnames(rawcounts) <- samp_lab

  fwrite(data.frame(rawcounts),
         file      = file.path(PB_DIR, paste0("seaad_cell_expression_by_donor_subclass_sum_", region, ".csv")),
         row.names = TRUE)

  new_sif        <- cell_anno[, colnames(cell_anno) %in% c("Subclass", "Donor ID", "Supertype")]
  new_sif$unique <- paste0(new_sif$`Donor ID`, new_sif$Subclass)
  new_sif        <- new_sif[!duplicated(new_sif$unique), ]
  new_sif$label  <- paste0(new_sif$Subclass, "_", new_sif$`Donor ID`)
  new_sif        <- new_sif[match(samp_lab, new_sif$label), ]
  fwrite(new_sif[, -which(colnames(new_sif) == "unique")],
         file = file.path(PB_DIR, paste0("seaad_cell_annotations_by_donor_subclass_sum_", region, ".csv")))

  message("SEAAD ", region, " pseudobulked")
}

pseudobulk_seaad(SEAAD_EXPR_DFC, SEAAD_ANNO_DFC, region = "DFC")
pseudobulk_seaad(SEAAD_EXPR_MTG, SEAAD_ANNO_MTG, region = "MTG")


# ============================================================
# STEP 3: edgeR DE — Jorstad, matched to the PER-REGION common gene
#         universe (Jorstad ∩ Gabitto within each region)
#         Parallelised over regions
# ============================================================

# Per-region common gene universe (gene names only — no need to read full matrices).
# DFC and MTG are intersected independently; no matching across regions.
common_genes_by_region <- setNames(lapply(REGIONS, function(i) {
  jg <- fread(file.path(PB_DIR, paste0("jorstad_cell_expression_by_donor_subclass_sum_", i, ".csv")),
              data.table = FALSE, select = 1L)[[1]]
  sg <- fread(file.path(PB_DIR, paste0("seaad_cell_expression_by_donor_subclass_sum_", i, ".csv")),
              data.table = FALSE, select = 1L)[[1]]
  # Canonicalize to current HGNC before intersecting (see helpers above)
  intersect(unique(canonicalize_symbols(jg)), unique(canonicalize_symbols(sg)))
}), REGIONS)

# NOTE: regions run sequentially so the per-celltype future_lapply inside
# calc_pseudobulk_DE can use all workers (future does not nest parallelism).
lapply(REGIONS, function(i) {

  rawcounts    <- read_pb_canonical(
    file.path(PB_DIR, paste0("jorstad_cell_expression_by_donor_subclass_sum_", i, ".csv")))
  cell_anno_pb <- fread(
    file.path(PB_DIR, paste0("jorstad_cell_annotations_by_donor_subclass_sum_", i, ".csv")),
    data.table = FALSE)

  # Restrict to the per-region common genes
  rawcounts <- rawcounts[rawcounts[, 1] %in% common_genes_by_region[[i]], ]

  save_path <- if (i == "DFC") EDGER_JORSTAD_DFC else EDGER_JORSTAD_MTG

  calc_pseudobulk_DE(
    rawcounts,
    save_path = save_path,
    celltype  = cell_anno_pb$Cell_Type,
    blocking  = cell_anno_pb[1],
    subset    = which(cell_anno_pb$Region == i)
  )

  message("Jorstad ", i, " edgeR DE done")
})


# ============================================================
# STEP 4: edgeR DE — SEAAD 2024, matched to the PER-REGION common gene
#         universe (same set used for Jorstad in each region)
#         Run for DFC and MTG
# ============================================================

# DFC
seaad_rc_dfc   <- read_pb_canonical(file.path(PB_DIR, "seaad_cell_expression_by_donor_subclass_sum_DFC.csv"))
seaad_anno_dfc <- fread(file.path(PB_DIR, "seaad_cell_annotations_by_donor_subclass_sum_DFC.csv"), data.table = FALSE)
seaad_rc_dfc   <- seaad_rc_dfc[seaad_rc_dfc[, 1] %in% common_genes_by_region[["DFC"]], ]

calc_pseudobulk_DE(
  seaad_rc_dfc,
  save_path = EDGER_SEAAD_DFC,
  celltype  = seaad_anno_dfc$Subclass,
  blocking  = seaad_anno_dfc[1]  # block by Donor ID
)
message("SEAAD DFC edgeR DE done")

# MTG
seaad_rc_mtg   <- read_pb_canonical(file.path(PB_DIR, "seaad_cell_expression_by_donor_subclass_sum_MTG.csv"))
seaad_anno_mtg <- fread(file.path(PB_DIR, "seaad_cell_annotations_by_donor_subclass_sum_MTG.csv"), data.table = FALSE)
seaad_rc_mtg   <- seaad_rc_mtg[seaad_rc_mtg[, 1] %in% common_genes_by_region[["MTG"]], ]

calc_pseudobulk_DE(
  seaad_rc_mtg,
  save_path = EDGER_SEAAD_MTG,
  celltype  = seaad_anno_mtg$Subclass,
  blocking  = seaad_anno_mtg[1]  # block by Donor ID
)
message("SEAAD MTG edgeR DE done")


# ============================================================
# STEP 5: DESeq2 — Jorstad vs SEAAD 2024, both DFC and MTG
#         Load matched pseudobulk data, find common genes, run DESeq2
#         Parallelised over regions; celltypes parallelised within each region
# ============================================================

# DESeq2 helper — one-vs-all per celltype
run_deseq_one_vs_all <- function(ct, anno, counts, ct_col) {
  anno_temp <- anno |>
    mutate(designcol = factor(
      ifelse(.data[[ct_col]] == ct, "ct", "all"),
      levels = c("all", "ct")  # "all" = reference
    ))

  dds  <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = anno_temp,
    design    = ~ Donor + designcol)

  dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor,
                BPPARAM = BPPARAM_SERIAL,  # no nested parallelism
                quiet   = TRUE)

  res <- results(dds1, name = "designcol_ct_vs_all")

  as.data.frame(res@listData) |>
    `rownames<-`(res@rownames) |>
    arrange(padj)
}

# Output path lookup — indexed by region
DESEQ2_OUT <- list(
  jorstad = list(DFC = DESEQ2_JORSTAD_DFC, MTG = DESEQ2_JORSTAD_MTG),
  sea     = list(DFC = DESEQ2_SEAAD_DFC,   MTG = DESEQ2_SEAAD_MTG)
)

# NOTE: regions run sequentially so the per-celltype future_lapply below can
# use all workers (future does not nest parallelism).
lapply(REGIONS, function(region) {

  # Load both datasets for this region in parallel
  rawcounts_list <- future_lapply(
    list(
      jorstad = file.path(PB_DIR, paste0("jorstad_cell_expression_by_donor_subclass_sum_", region, ".csv")),
      sea     = file.path(PB_DIR, paste0("seaad_cell_expression_by_donor_subclass_sum_",   region, ".csv"))
    ),
    \(path) read_pb_canonical(path) |> tibble::column_to_rownames(var = "V1"),
    future.seed    = TRUE,
    future.globals = c("read_pb_canonical", "canonicalize_symbols", "ALIAS_VEC")
  )
  names(rawcounts_list) <- c("jorstad", "sea")

  # Match genes between the two datasets
  common_genes   <- intersect(rownames(rawcounts_list$jorstad), rownames(rawcounts_list$sea))
  rawcounts_list <- lapply(rawcounts_list, \(x) x[rownames(x) %in% common_genes, ])

  cell_anno_pb <- list(
    "jorstad" = fread(
      file.path(PB_DIR, paste0("jorstad_cell_annotations_by_donor_subclass_sum_", region, ".csv")),
      data.table = FALSE) |>
      mutate(label     = colnames(rawcounts_list$jorstad),
             Cell_Type = factor(Cell_Type),
             Donor     = factor(Donor)),
    "sea" = fread(
      file.path(PB_DIR, paste0("seaad_cell_annotations_by_donor_subclass_sum_", region, ".csv")),
      data.table = FALSE) |>
      mutate(ID       = colnames(rawcounts_list$sea),
             Subclass = factor(Subclass),
             Donor    = factor(`Donor ID`),
             Subclass = recode(Subclass,
                               "Astrocyte"       = "Astro",
                               "Endothelial"     = "Endo",
                               "Microglia-PVM"   = "Micro/PVM",
                               "Oligodendrocyte" = "Oligo"))
  )

  cts <- unique(cell_anno_pb[["jorstad"]]$Cell_Type)

  # Jorstad
  reslist_jorstad <- future_lapply(cts, run_deseq_one_vs_all,
                                   anno   = cell_anno_pb$jorstad,
                                   counts = rawcounts_list$jorstad,
                                   ct_col = "Cell_Type",
                                   future.seed = TRUE)
  names(reslist_jorstad) <- cts
  out_path <- DESEQ2_OUT$jorstad[[region]]
  if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
  qsave(reslist_jorstad, file = out_path)
  message("Jorstad ", region, " DESeq2 done")

  # SEAAD
  reslist_sea <- future_lapply(cts, run_deseq_one_vs_all,
                               anno   = cell_anno_pb$sea,
                               counts = rawcounts_list$sea,
                               ct_col = "Subclass",
                               future.seed = TRUE)
  names(reslist_sea) <- cts
  out_path <- DESEQ2_OUT$sea[[region]]
  if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
  qsave(reslist_sea, file = out_path)
  message("SEAAD ", region, " DESeq2 done")
})

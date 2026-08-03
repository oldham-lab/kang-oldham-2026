library(qs)
library(data.table)
library(DESeq2)
library(tidyverse)
library(future.apply)
library(BiocParallel)

# ============================================================
# INPUT / OUTPUT PATHS  —  edit these before running
# ============================================================

# Pseudobulk CSVs (already generated)
PB_DIR <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/DE_v2.4/pseudobulk")

# Output files
DESEQ2_JORSTAD_MTG <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/DE_v2.4/DESeq2/jorstad/DESeq2_jorstad_MTG_subclass_matchedGenes.qs")
DESEQ2_SEAAD_MTG   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/DE_v2.4/DESeq2/seaad/DESeq2_seaad_MTG_subclass_matchedGenes.qs")


# ============================================================
# PARALLELISM SETUP
# ============================================================
N_WORKERS      <- min(7L, parallelly::availableCores() - 1L)
plan(multisession, workers = N_WORKERS)
BPPARAM_SERIAL <- BiocParallel::SerialParam()


# ============================================================
# Load MTG pseudobulk data and match genes
# ============================================================
rawcounts_list <- future_lapply(
  list(
    jorstad = file.path(PB_DIR, "jorstad_cell_expression_by_donor_subclass_sum_MTG.csv"),
    sea     = file.path(PB_DIR, "seaad_cell_expression_by_donor_subclass_sum_MTG.csv")
  ),
  \(path) fread(path, data.table = FALSE) |> tibble::column_to_rownames(var = "V1"),
  future.seed = TRUE
)
names(rawcounts_list) <- c("jorstad", "sea")

common_genes   <- intersect(rownames(rawcounts_list$jorstad), rownames(rawcounts_list$sea))
rawcounts_list <- lapply(rawcounts_list, \(x) x[rownames(x) %in% common_genes, ])

cell_anno_pb <- list(
  "jorstad" = fread(
    file.path(PB_DIR, "jorstad_cell_annotations_by_donor_subclass_sum_MTG.csv"),
    data.table = FALSE) |>
    mutate(label     = colnames(rawcounts_list$jorstad),
           Cell_Type = factor(Cell_Type),
           Donor     = factor(Donor)),
  "sea" = fread(
    file.path(PB_DIR, "seaad_cell_annotations_by_donor_subclass_sum_MTG.csv"),
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


# ============================================================
# DESeq2 helper — one-vs-all per celltype
# ============================================================
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
                BPPARAM = BPPARAM_SERIAL,
                quiet   = TRUE)

  res <- results(dds1, name = "designcol_ct_vs_all")

  as.data.frame(res@listData) |>
    `rownames<-`(res@rownames) |>
    arrange(padj)
}


# ============================================================
# Run DESeq2 — Jorstad MTG
# ============================================================
reslist_jorstad <- future_lapply(cts, run_deseq_one_vs_all,
                                 anno   = cell_anno_pb$jorstad,
                                 counts = rawcounts_list$jorstad,
                                 ct_col = "Cell_Type",
                                 future.seed = TRUE)
names(reslist_jorstad) <- cts
if (!dir.exists(dirname(DESEQ2_JORSTAD_MTG))) dir.create(dirname(DESEQ2_JORSTAD_MTG), recursive = TRUE)
qsave(reslist_jorstad, file = DESEQ2_JORSTAD_MTG)
message("Jorstad MTG DESeq2 done")


# ============================================================
# Run DESeq2 — SEAAD MTG
# ============================================================
reslist_sea <- future_lapply(cts, run_deseq_one_vs_all,
                             anno   = cell_anno_pb$sea,
                             counts = rawcounts_list$sea,
                             ct_col = "Subclass",
                             future.seed = TRUE)
names(reslist_sea) <- cts
if (!dir.exists(dirname(DESEQ2_SEAAD_MTG))) dir.create(dirname(DESEQ2_SEAAD_MTG), recursive = TRUE)
qsave(reslist_sea, file = DESEQ2_SEAAD_MTG)
message("SEAAD MTG DESeq2 done")

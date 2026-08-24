# UpSet plot using UpSetR package
#
# Unit of analysis: (Gene, Celltype, Direction) tuples.
# A tuple belongs to a set only if that exact combination appears in that dataset.
# A tuple is considered shared between sets only if the same Gene, Celltype,
# and Direction all match.
#
# One plot with 8 groups: 4 datasets × 2 regions (DFC, MTG).

library(qs)
library(data.table)
library(tidyverse)
library(UpSetR)
library(svglite)
library(showtext)
showtext_auto()

# ============================================================
# PLACEHOLDERS
# ============================================================

# --- DFC ---
MIT_EDGER_PATH_DFC <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/MIT_DE_results_DFC/mit_dfc_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH_DFC <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/MIT_DE_results_DFC/mit_dfc_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH_DFC <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/SEA_DE_results_DFC/sea_dfc_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH_DFC <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/SEA_DE_results_DFC/sea_dfc_DESeq2_ADvsCon_by_celltype.qs")

# --- MTG ---
MIT_EDGER_PATH_MTG <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/MIT_DE_results_MTG/mit_mtg_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH_MTG <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/MIT_DE_results_MTG/mit_mtg_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH_MTG <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/SEA_DE_results_MTG/sea_mtg_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH_MTG <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/DE_ADvsCon/SEA_DE_results_MTG/sea_mtg_DESeq2_ADvsCon_by_celltype.qs")

#SAVE_DIR   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4")
SAVE_DIR <- Sys.getenv("SCRATCH_DIR", "~/test")

FDR_THRESH <- 0.05
LFC_THRESH <- 0

PLOT_W <- 6
PLOT_H <- 3.5

# ============================================================
# LOAD FILES
# ============================================================
message("Loading DE results ...")

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

de_dfc <- list(
  MIT_edgeR  = load_edger(MIT_EDGER_PATH_DFC),
  MIT_DESeq2 = load_deseq(MIT_DESEQ_PATH_DFC),
  SEA_edgeR  = load_edger(SEA_EDGER_PATH_DFC),
  SEA_DESeq2 = load_deseq(SEA_DESEQ_PATH_DFC)
)

de_mtg <- list(
  MIT_edgeR  = load_edger(MIT_EDGER_PATH_MTG),
  MIT_DESeq2 = load_deseq(MIT_DESEQ_PATH_MTG),
  SEA_edgeR  = load_edger(SEA_EDGER_PATH_MTG),
  SEA_DESeq2 = load_deseq(SEA_DESEQ_PATH_MTG)
)

# ============================================================
# EXTRACT SIGNIFICANT (GENE, CELLTYPE, DIRECTION) TUPLES PER DATASET
# A tuple is included if the gene passes thresholds in that celltype.
# ============================================================
get_sig_tuples <- function(de_list) {
  tuples <- lapply(names(de_list), function(ct) {
    df <- de_list[[ct]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    keep <- !is.na(df$FDR) & !is.na(df$logFC) &
            df$FDR < FDR_THRESH & abs(df$logFC) > LFC_THRESH
    df <- df[keep, ]
    if (nrow(df) == 0) return(NULL)
    direction <- ifelse(df$logFC > 0, "up", "down")
    paste(df$gene, ct, direction, sep = "|||")
  })
  unique(unlist(Filter(Negate(is.null), tuples)))
}

# ============================================================
# BUILD COMBINED 8-GROUP MEMBERSHIP MATRIX
# ============================================================
SET_KEYS <- c(
  "MIT_DESeq2_DFC",  "MIT_edgeR_DFC",  "SEA_DESeq2_DFC",  "SEA_edgeR_DFC",
  "MIT_DESeq2_MTG",  "MIT_edgeR_MTG",  "SEA_DESeq2_MTG",  "SEA_edgeR_MTG"
)

message("\n--- Building combined 8-group membership matrix ---")

gene_sets <- list(
  MIT_DESeq2_DFC = get_sig_tuples(de_dfc$MIT_DESeq2),
  MIT_edgeR_DFC  = get_sig_tuples(de_dfc$MIT_edgeR),
  SEA_DESeq2_DFC = get_sig_tuples(de_dfc$SEA_DESeq2),
  SEA_edgeR_DFC  = get_sig_tuples(de_dfc$SEA_edgeR),
  MIT_DESeq2_MTG = get_sig_tuples(de_mtg$MIT_DESeq2),
  MIT_edgeR_MTG  = get_sig_tuples(de_mtg$MIT_edgeR),
  SEA_DESeq2_MTG = get_sig_tuples(de_mtg$SEA_DESeq2),
  SEA_edgeR_MTG  = get_sig_tuples(de_mtg$SEA_edgeR)
)

all_ids <- unique(unlist(gene_sets))
message("Total unique (gene, celltype, direction) combinations across all datasets: ", length(all_ids))

# Unique gene count per set (gene is the first element of each tuple key)
# Used to label the set-size bars in the UpSet plot.
set_gene_counts_raw <- sapply(SET_KEYS, function(k) {
  length(unique(sub("\\|\\|\\|.*", "", gene_sets[[k]])))
})

mat <- as.data.frame(
  setNames(
    lapply(gene_sets, function(tuples) as.integer(all_ids %in% tuples)),
    SET_KEYS
  )
)

n_tuples <- length(all_ids)
n_shared <- sum(rowSums(mat) > 1)
message(n_tuples, " unique (gene, celltype, direction) combinations; ", n_shared, " shared in >=2 datasets")
message("Column sums:"); print(colSums(mat))
message("Row sum distribution:"); print(table(rowSums(mat)))

stopifnot(nrow(mat) > 0, any(mat > 0))

# Rename columns for display: SEA -> Gabitto, MIT -> Liu
# Reorder parts to {Dataset}_{Region}_{Method}
plot_keys <- gsub("^SEA", "Gabitto", gsub("^MIT", "Liu", SET_KEYS))
plot_keys <- gsub("^(Gabitto|Liu)_(edgeR|DESeq2)_(DFC|MTG)$", "\\1_\\3_\\2", plot_keys)
names(mat) <- plot_keys

# Gene counts in plot_keys order for use inside do_upset()
names(set_gene_counts_raw) <- plot_keys
set_gene_counts <- set_gene_counts_raw[plot_keys]

# ============================================================
# PLOT
# ============================================================
do_upset <- function() {
  # Patch 1: fix bottom_margin overflow and fix bar label alignment for
  # vertical text (angle=90): vjust=0.5 centres on bar, hjust=-0.25 floats
  # the label just above the bar top.
  orig_fn  <- getFromNamespace("Make_main_bar", "UpSetR")
  fn_src   <- paste(deparse(body(orig_fn)), collapse = "\n")
  fn_src   <- sub("bottom_margin <- \\(-1\\) \\* 0\\.65", "bottom_margin <- 0", fn_src)
  fn_src   <- gsub("vjust = -1,", "vjust = 0.5, hjust = -0.25,", fn_src)
  patched_fn        <- orig_fn
  body(patched_fn)  <- parse(text = fn_src)[[1]]
  assignInNamespace("Make_main_bar", patched_fn, "UpSetR")
  on.exit(assignInNamespace("Make_main_bar", orig_fn, "UpSetR"), add = TRUE)

  # Patch 2: replace set-size bars with unique gene counts per set.
  orig_size_fn <- getFromNamespace("Make_size_plot", "UpSetR")
  gene_counts  <- set_gene_counts
  patched_size_fn <- function(Set_size_data, ...) {
    Set_size_data$y <- unname(gene_counts)[Set_size_data$x]
    orig_size_fn(Set_size_data, ...)
  }
  assignInNamespace("Make_size_plot", patched_size_fn, "UpSetR")
  on.exit(assignInNamespace("Make_size_plot", orig_size_fn, "UpSetR"), add = TRUE)

  # Patch 3: prevent UpSetR's internal grid.newpage() from advancing the PDF
  # device to a new page, which would leave page 1 blank.
  orig_newpage <- grid::grid.newpage
  unlockBinding("grid.newpage", asNamespace("grid"))
  assign("grid.newpage", function(...) invisible(NULL), envir = asNamespace("grid"))
  on.exit({
    assign("grid.newpage", orig_newpage, envir = asNamespace("grid"))
    lockBinding("grid.newpage", asNamespace("grid"))
  }, add = TRUE)

  UpSetR::upset(
    mat,
    sets            = plot_keys,
    keep.order      = TRUE,
    order.by        = "freq",
    decreasing      = TRUE,
    mb.ratio        = c(0.6, 0.4),
    text.scale      = c(1.3, 1.2, 1, 1, 1.2, 1),
    point.size      = 2.5,
    line.size       = 0.8,
    mainbar.y.label = "# of shared genes\n(by celltype and direction)",
    sets.x.label    = "Total unique DE genes",
    number.angles   = 90
  )
}

if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)
stem <- file.path(SAVE_DIR, "panel_A_upset_DE_celltype_intersection_combined")

cairo_pdf(paste0(stem, ".pdf"), width = PLOT_W, height = PLOT_H)
print(do_upset())
dev.off()
message("Saved: ", basename(stem), ".pdf")

svglite::svglite(paste0(stem, ".svg"), width = PLOT_W, height = PLOT_H)
print(do_upset())
dev.off()
message("Saved: ", basename(stem), ".svg")

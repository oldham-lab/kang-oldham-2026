# UpSet plot using UpSetR package
#
# Membership: a gene is in dataset X's set if it is significantly DE in at least
#             one celltype in dataset X (either direction).
# Intersection: two datasets share a gene if that gene is DE in either dataset,
#               regardless of celltype or direction.
#
# Unit of analysis: individual genes (counted once per dataset).
# One plot with 8 groups: 4 datasets × 2 regions (DFC, MTG).

library(qs)
library(data.table)
library(tidyverse)
library(UpSetR)
library(svglite)  # install.packages("svglite")
library(showtext)
showtext_auto()

# ============================================================
# PLACEHOLDERS
# ============================================================

# --- DFC ---
MIT_EDGER_PATH_DFC <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH_DFC <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH_DFC <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH_DFC <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_DESeq2_ADvsCon_by_celltype.qs")

# --- MTG ---
MIT_EDGER_PATH_MTG <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results_MTG/mit_mtg_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH_MTG <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results_MTG/mit_mtg_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH_MTG <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results_MTG/sea_mtg_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH_MTG <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results_MTG/sea_mtg_DESeq2_ADvsCon_by_celltype.qs")

SAVE_DIR   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3")

FDR_THRESH <- 0.05
LFC_THRESH <- 0

PLOT_W <- 6
PLOT_H <- 3.5

# ============================================================
# CELLTYPE MAPPING
# ============================================================
map_list <- list(
  c("Endothelial",  "SMC", "VLMC", "End", "Per"),
  c("L4 IT",        "Exc L4-5 IT-2", "Exc L3-4 IT", "Exc L4-5 IT-1"),
  c("L5 ET",        "Exc L5 ET"),
  c("L5 IT",        "Exc L4-5 IT-2", "Exc L4-5 IT-1", "Exc L3-5 IT", "Exc L5-6 IT"),
  c("L5/6 NP",      "Exc L5/6 NP"),
  c("Lamp5",        "Inh LAMP5"),
  c("Pvalb",        "Inh PVALB"),
  c("Sst",          "Inh SST"),
  c("L6 IT",        "Exc L5-6 IT"),
  c("L6 IT Car3",   "Exc L5/6 IT Car3"),
  c("L6 CT",        "Exc L6 CT"),
  c("Pax6",         "Inh PAX6"),
  c("Astrocyte",    "Ast"),
  c("OPC",          "OPC"),
  c("Vip",          "Inh VIP"),
  c("L6b",          "Exc L6b"),
  c("L2/3 IT",      "Exc L2-3 IT")
)

sea_names  <- sapply(map_list, `[`, 1)
sea_to_mit <- setNames(lapply(map_list, `[`, -1), sea_names)

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

raw_dfc <- list(
  MIT_edgeR  = load_edger(MIT_EDGER_PATH_DFC),
  MIT_DESeq2 = load_deseq(MIT_DESEQ_PATH_DFC),
  SEA_edgeR  = load_edger(SEA_EDGER_PATH_DFC),
  SEA_DESeq2 = load_deseq(SEA_DESEQ_PATH_DFC)
)

raw_mtg <- list(
  MIT_edgeR  = load_edger(MIT_EDGER_PATH_MTG),
  MIT_DESeq2 = load_deseq(MIT_DESEQ_PATH_MTG),
  SEA_edgeR  = load_edger(SEA_EDGER_PATH_MTG),
  SEA_DESeq2 = load_deseq(SEA_DESEQ_PATH_MTG)
)

# ============================================================
# TRANSLATE MIT TO SEA CANONICAL CELLTYPE NAMES
# ============================================================
pool_mit <- function(de_list) {
  setNames(lapply(sea_names, function(sea_ct) {
    parts <- Filter(Negate(is.null), lapply(sea_to_mit[[sea_ct]], function(m) de_list[[m]]))
    if (length(parts) == 0) return(NULL)
    dt <- rbindlist(parts)
    as.data.frame(dt[dt[, .I[which.min(FDR)], by = gene]$V1])
  }), sea_names)
}

make_de <- function(raw) {
  list(
    MIT_edgeR  = pool_mit(raw$MIT_edgeR),
    MIT_DESeq2 = pool_mit(raw$MIT_DESeq2),
    SEA_edgeR  = raw$SEA_edgeR,
    SEA_DESeq2 = raw$SEA_DESeq2
  )
}

de_dfc <- make_de(raw_dfc)
de_mtg <- make_de(raw_mtg)

# ============================================================
# EXTRACT SIGNIFICANT GENES PER DATASET
# A gene is included if it passes thresholds in at least one celltype.
# ============================================================
get_sig_genes <- function(de_list) {
  genes <- lapply(names(de_list), function(ct) {
    df <- de_list[[ct]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    keep <- !is.na(df$FDR) & !is.na(df$logFC) &
            df$FDR < FDR_THRESH & abs(df$logFC) > LFC_THRESH
    df$gene[keep]
  })
  unique(unlist(Filter(Negate(is.null), genes)))
}

# ============================================================
# BUILD COMBINED 8-GROUP MEMBERSHIP MATRIX
# ============================================================
SET_KEYS <- c(
  "MIT_edgeR_DFC",  "MIT_DESeq2_DFC",  "SEA_edgeR_DFC",  "SEA_DESeq2_DFC",
  "MIT_edgeR_MTG",  "MIT_DESeq2_MTG",  "SEA_edgeR_MTG",  "SEA_DESeq2_MTG"
)

message("\n--- Building combined 8-group membership matrix ---")

gene_sets <- list(
  MIT_edgeR_DFC  = get_sig_genes(de_dfc$MIT_edgeR),
  MIT_DESeq2_DFC = get_sig_genes(de_dfc$MIT_DESeq2),
  SEA_edgeR_DFC  = get_sig_genes(de_dfc$SEA_edgeR),
  SEA_DESeq2_DFC = get_sig_genes(de_dfc$SEA_DESeq2),
  MIT_edgeR_MTG  = get_sig_genes(de_mtg$MIT_edgeR),
  MIT_DESeq2_MTG = get_sig_genes(de_mtg$MIT_DESeq2),
  SEA_edgeR_MTG  = get_sig_genes(de_mtg$SEA_edgeR),
  SEA_DESeq2_MTG = get_sig_genes(de_mtg$SEA_DESeq2)
)

all_ids <- unique(unlist(gene_sets))
message("Total unique DE genes across all datasets: ", length(all_ids))

mat <- as.data.frame(
  setNames(
    lapply(gene_sets, function(genes) as.integer(all_ids %in% genes)),
    SET_KEYS
  )
)

n_genes  <- length(all_ids)
n_shared <- sum(rowSums(mat) > 1)
message(n_genes, " unique genes; ", n_shared, " shared in >=2 datasets")
message("Column sums:"); print(colSums(mat))
message("Row sum distribution:"); print(table(rowSums(mat)))

stopifnot(nrow(mat) > 0, any(mat > 0))

# Rename columns for display: SEA -> Gabitto, MIT -> Liu
# Rename SEA->Gabitto, MIT->Liu, then reorder parts to {Dataset}_{Region}_{Method}
plot_keys <- gsub("^SEA", "Gabitto", gsub("^MIT", "Liu", SET_KEYS))
plot_keys <- gsub("^(Gabitto|Liu)_(edgeR|DESeq2)_(DFC|MTG)$", "\\1_\\3_\\2", plot_keys)
names(mat) <- plot_keys

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

  # Patch 2: prevent UpSetR's internal grid.newpage() from advancing the PDF
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
    mainbar.y.label = "Intersection size",
    sets.x.label    = "Set size",
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

# UpSet plot using UpSetR package — dCoPA gene list input
#
# Unit of analysis: (Gene, Celltype, Direction) tuples.
# A tuple belongs to a set only if that exact combination appears in that comparison.
# Sharing between sets requires matching Gene, Celltype, AND Direction.
#
# Sets (8): one per Comparison value in the input CSV.

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

#DCOPA_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/panel_B_dcopa_genelist.csv")
DCOPA_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.1/panel_B_dcopa_genelist.csv")
SAVE_DIR   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4")
#SAVE_DIR <- "/home/gugene/test/"

PLOT_W <- 8
PLOT_H <- 3.5

# ============================================================
# SET ORDER
# ============================================================
SET_KEYS <- c(
  "Gabitto_AllADVsCon_DFC_ROSMAP",
  "Liu_AllADVsCon_DFC_ROSMAP",
  "Gabitto_AllADVsCon_MTG_ROSMAP",
  "Liu_AllADVsCon_MTG_ROSMAP",
  "Gabitto_AllADVsCon_DFC",
  "Liu_AllADVsCon_DFC",
  "Gabitto_AllADVsCon_MTG",
  "Liu_AllADVsCon_MTG"
)

# ============================================================
# LOAD & NORMALISE
# ============================================================
message("Loading dCoPA gene list ...")
dcopa_raw <- fread(DCOPA_PATH, data.table = FALSE)

stopifnot(
  all(c("Comparison", "Dataset", "Direction", "Celltype", "Gene") %in% names(dcopa_raw))
)

missing_comp <- setdiff(SET_KEYS, unique(dcopa_raw$Comparison))
if (length(missing_comp) > 0)
  stop("Missing Comparisons in input: ", paste(missing_comp, collapse = ", "))

dcopa <- dcopa_raw |>
  dplyr::mutate(
    direction = dplyr::case_when(
      Direction == "Higher in more severe" ~ "higher",
      Direction == "Lower in more severe"  ~ "lower",
      TRUE                                 ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(direction))

n_excluded <- nrow(dcopa_raw) - nrow(dcopa)
message("  Excluded ", n_excluded, " rows (unknown direction)")
message("  Retained ", nrow(dcopa), " rows across ",
        length(unique(dcopa$Comparison)), " comparisons")

# ============================================================
# BUILD (GENE, CELLTYPE, DIRECTION) TUPLE SETS PER COMPARISON
# ============================================================
gene_sets <- setNames(lapply(SET_KEYS, function(comp) {
  rows <- dcopa[dcopa$Comparison == comp, ]
  unique(paste(rows$Gene, rows$Celltype, rows$direction, sep = "|||"))
}), SET_KEYS)

all_ids <- unique(unlist(gene_sets))
message("Total unique (gene, celltype, direction) combinations across all comparisons: ", length(all_ids))

# Unique gene counts per set (for set-size bar display)
set_gene_counts_raw <- sapply(SET_KEYS, function(k) {
  length(unique(sub("\\|\\|\\|.*", "", gene_sets[[k]])))
})

# ============================================================
# BUILD MEMBERSHIP MATRIX
# ============================================================
mat <- as.data.frame(
  setNames(
    lapply(gene_sets, function(tuples) as.integer(all_ids %in% tuples)),
    SET_KEYS
  )
)

n_tuples <- length(all_ids)
n_shared <- sum(rowSums(mat) > 1)
message(n_tuples, " unique tuples; ", n_shared, " shared in >=2 comparisons")
message("Set sizes (tuples):"); print(colSums(mat))
message("Unique gene counts:"); print(set_gene_counts_raw)
message("Row sum distribution:"); print(table(rowSums(mat)))

stopifnot(nrow(mat) > 0, any(mat > 0))

# Rename columns for display: remove AllADVsCon_, replace ROSMAP with AD_modules
plot_keys <- gsub("_AllADVsCon_", "_", SET_KEYS)
plot_keys <- gsub("_ROSMAP$", "_AD_modules", plot_keys)
names(mat) <- plot_keys

plot_order <- c(
  "Liu_DFC_AD_modules",     "Gabitto_DFC_AD_modules",
  "Liu_MTG_AD_modules",     "Gabitto_MTG_AD_modules",
  "Liu_DFC",                "Gabitto_DFC",
  "Liu_MTG",                "Gabitto_MTG"
)
mat <- mat[, plot_order]

names(set_gene_counts_raw) <- plot_keys
set_gene_counts <- set_gene_counts_raw[plot_order]

# ============================================================
# PLOT
# ============================================================
do_upset <- function() {
  # Patch 1: fix bottom_margin overflow and bar label alignment for vertical text.
  orig_fn  <- getFromNamespace("Make_main_bar", "UpSetR")
  fn_src   <- paste(deparse(body(orig_fn)), collapse = "\n")
  fn_src   <- sub("bottom_margin <- \\(-1\\) \\* 0\\.65", "bottom_margin <- 0", fn_src)
  fn_src   <- gsub("vjust = -1,", "vjust = 0.5, hjust = -0.1,", fn_src)
  fn_src   <- sub("scale_y_continuous\\(", "scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.8)), ", fn_src)
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

  # Patch 3: prevent UpSetR's internal grid.newpage() from creating a blank PDF page.
  orig_newpage <- grid::grid.newpage
  unlockBinding("grid.newpage", asNamespace("grid"))
  assign("grid.newpage", function(...) invisible(NULL), envir = asNamespace("grid"))
  on.exit({
    assign("grid.newpage", orig_newpage, envir = asNamespace("grid"))
    lockBinding("grid.newpage", asNamespace("grid"))
  }, add = TRUE)

  UpSetR::upset(
    mat,
    sets            = plot_order,
    keep.order      = TRUE,
    order.by        = "freq",
    decreasing      = TRUE,
    mb.ratio        = c(0.6, 0.4),
    text.scale      = c(1.3, 1.2, 1, 1, 1.2, 1),
    point.size      = 2.5,
    line.size       = 0.8,
    show.numbers    = "yes",
    mainbar.y.label = "# of shared genes\n(by celltype and direction)",
    sets.x.label    = "Total unique dCoPA genes",
    number.angles   = 90
  )
}

if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)
stem <- file.path(SAVE_DIR, "panel_B_upset_dcopa_celltype_intersection_combined")

cairo_pdf(paste0(stem, ".pdf"), width = PLOT_W, height = PLOT_H)
print(do_upset())
dev.off()
message("Saved: ", basename(stem), ".pdf")

svglite::svglite(paste0(stem, ".svg"), width = PLOT_W, height = PLOT_H)
print(do_upset())
dev.off()
message("Saved: ", basename(stem), ".svg")

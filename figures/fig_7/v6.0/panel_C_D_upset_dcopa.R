# UpSet plot using UpSetR package — dCoPA gene list input
#
# Produces two separate plots: one for DFC, one for MTG.
# Each plot contains 4 sets (Gabitto + Liu × cohort vs AD-modules).
#
# Unit of analysis: individual genes.
# A gene belongs to a set if it appears in ANY celltype or direction in that comparison.

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

DCOPA_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v4/panel_B_dcopa_genelist.csv")
SAVE_DIR   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6")

PLOT_W <- 5.821
PLOT_H <- 2.94

# ============================================================
# REGION-SPECIFIC SET KEYS
# ============================================================
SET_KEYS_DFC <- c(
  "Gabitto_AllADVsCon_DFC_ROSMAP",
  "Liu_AllADVsCon_DFC_ROSMAP",
  "Gabitto_AllADVsCon_DFC",
  "Liu_AllADVsCon_DFC"
)

SET_KEYS_MTG <- c(
  "Gabitto_AllADVsCon_MTG_ROSMAP",
  "Liu_AllADVsCon_MTG_ROSMAP",
  "Gabitto_AllADVsCon_MTG",
  "Liu_AllADVsCon_MTG"
)

# ============================================================
# LOAD & NORMALISE
# ============================================================
message("Loading dCoPA gene list ...")
dcopa_raw <- fread(DCOPA_PATH, data.table = FALSE)

stopifnot(
  all(c("Comparison", "Gene") %in% names(dcopa_raw))
)

# Check all expected comparisons are present
all_keys <- c(SET_KEYS_DFC, SET_KEYS_MTG)
missing_comp <- setdiff(all_keys, unique(dcopa_raw$Comparison))
if (length(missing_comp) > 0)
  stop("Missing Comparisons in input: ", paste(missing_comp, collapse = ", "))

dcopa <- dcopa_raw
message("  Loaded ", nrow(dcopa), " rows across ",
        length(unique(dcopa$Comparison)), " comparisons")

# ============================================================
# HELPER: build membership matrix for a given set of keys
# ============================================================
build_mat <- function(set_keys, label) {
  gene_sets <- setNames(lapply(set_keys, function(comp) {
    unique(dcopa$Gene[dcopa$Comparison == comp])
  }), set_keys)

  all_ids <- unique(unlist(gene_sets))
  message(label, ": ", length(all_ids), " unique genes")

  mat <- as.data.frame(
    setNames(
      lapply(gene_sets, function(genes) as.integer(all_ids %in% genes)),
      set_keys
    )
  )

  # Rename for display
  label_map <- c(
    "Gabitto_AllADVsCon_DFC"        = "CTRL modules | Gabitto SN",
    "Liu_AllADVsCon_DFC"            = "CTRL modules | Liu SN",
    "Gabitto_AllADVsCon_DFC_ROSMAP" = "AD modules | Gabitto SN",
    "Liu_AllADVsCon_DFC_ROSMAP"     = "AD modules | Liu SN",
    "Gabitto_AllADVsCon_MTG"        = "CTRL modules | Gabitto SN",
    "Liu_AllADVsCon_MTG"            = "CTRL modules | Liu SN",
    "Gabitto_AllADVsCon_MTG_ROSMAP" = "AD modules | Gabitto SN",
    "Liu_AllADVsCon_MTG_ROSMAP"     = "AD modules | Liu SN"
  )
  names(mat) <- label_map[set_keys]
  mat
}

# ============================================================
# BUILD MATRICES
# ============================================================
mat_dfc <- build_mat(SET_KEYS_DFC, "DFC")
mat_mtg <- build_mat(SET_KEYS_MTG, "MTG")

# Display order within each plot (Liu before Gabitto so Gabitto renders on top)
plot_order_dfc <- c(
  "AD modules | Liu SN", "AD modules | Gabitto SN",
  "CTRL modules | Liu SN", "CTRL modules | Gabitto SN"
)

plot_order_mtg <- c(
  "AD modules | Liu SN", "AD modules | Gabitto SN",
  "CTRL modules | Liu SN", "CTRL modules | Gabitto SN"
)

mat_dfc <- mat_dfc[, plot_order_dfc]
mat_mtg <- mat_mtg[, plot_order_mtg]

# ============================================================
# PLOT FUNCTION
# ============================================================
do_upset <- function(mat, plot_order) {
  # Patch 1: fix bottom_margin overflow and bar label alignment for vertical text
  orig_fn  <- getFromNamespace("Make_main_bar", "UpSetR")
  fn_src   <- paste(deparse(body(orig_fn)), collapse = "\n")
  fn_src   <- sub("bottom_margin <- \\(-1\\) \\* 0\\.65", "bottom_margin <- 0", fn_src)
  fn_src   <- gsub("vjust = -1,", "vjust = 0.5, hjust = -0.25,", fn_src)
  fn_src   <- sub("top_margin_main <- 3", "top_margin_main <- 6", fn_src)
  patched_fn        <- orig_fn
  body(patched_fn)  <- parse(text = fn_src)[[1]]
  assignInNamespace("Make_main_bar", patched_fn, "UpSetR")
  on.exit(assignInNamespace("Make_main_bar", orig_fn, "UpSetR"), add = TRUE)

  # Patch 2: prevent UpSetR's internal grid.newpage() from advancing the device
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
    mb.ratio        = c(0.7, 0.3),
    text.scale      = c(1.3, 1.2, 1, 1, 1.2, 1),
    point.size      = 2.5,
    line.size       = 0.8,
    show.numbers    = "yes",
    mainbar.y.label = "# of overlapping genes",
    sets.x.label    = "# of unique genes",
    number.angles   = 90
  )
}

# ============================================================
# SAVE PLOTS
# ============================================================
if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)

for (region in c("DFC", "MTG")) {
  mat        <- if (region == "DFC") mat_dfc else mat_mtg
  plot_order <- if (region == "DFC") plot_order_dfc else plot_order_mtg
  stem       <- file.path(SAVE_DIR,
                          paste0("panel_C_D_upset_dcopa_", tolower(region)))

  cairo_pdf(paste0(stem, ".pdf"), width = PLOT_W, height = PLOT_H)
  print(do_upset(mat, plot_order))
  dev.off()
  message("Saved: ", basename(stem), ".pdf")

  svglite::svglite(paste0(stem, ".svg"), width = PLOT_W, height = PLOT_H)
  print(do_upset(mat, plot_order))
  dev.off()
  message("Saved: ", basename(stem), ".svg")
}

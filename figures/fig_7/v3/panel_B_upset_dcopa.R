# UpSet plot using UpSetR package — dCoPA gene list input
#
# Analogous to panel_A_upset_DE.R but using dCoPA module genes as input
# instead of pseudobulk DE results.
#
# Unit of analysis: individual genes.
# A gene belongs to a set if it appears in ANY celltype or direction in that comparison.
#
# Celltype normalization (mirrors panel_A_upset_DE.R):
#   Gabitto rows → treated as SEA-AD; celltype names kept as-is (already canonical).
#   Liu rows     → treated as MIT; celltype names translated to SEA canonical via map_list.
#   Liu celltypes with no map entry are excluded.
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

DCOPA_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3/panel_B_dcopa_genelist.csv")
#SAVE_DIR   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3")
SAVE_DIR <- "/home/gugene/test/"

PLOT_W <- 8
PLOT_H <- 3.5

# ============================================================
# SET ORDER
# ============================================================
SET_KEYS <- c(
  # ROSMAP first → renders at bottom of matrix (UpSetR reverses set order)
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
# CELLTYPE MAPPING  (Liu MIT names → SEA canonical names)
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

# Flat lookup: MIT name -> SEA canonical name
mit_to_sea <- unlist(lapply(map_list, function(x) setNames(rep(x[1], length(x) - 1), x[-1])))

# ============================================================
# LOAD & NORMALISE
# ============================================================
message("Loading dCoPA gene list ...")
dcopa_raw <- fread(DCOPA_PATH, data.table = FALSE)

stopifnot(
  all(c("Comparison", "Dataset", "Direction", "Celltype", "Gene") %in% names(dcopa_raw))
)

# Check all expected comparisons are present
missing_comp <- setdiff(SET_KEYS, unique(dcopa_raw$Comparison))
if (length(missing_comp) > 0)
  stop("Missing Comparisons in input: ", paste(missing_comp, collapse = ", "))

dcopa <- dcopa_raw |>
  dplyr::mutate(
    canonical_ct = dplyr::case_when(
      Dataset == "Gabitto"                                    ~ Celltype,
      Dataset == "Liu" & Celltype %in% names(mit_to_sea)     ~ mit_to_sea[Celltype],
      TRUE                                                    ~ NA_character_
    ),
    direction = dplyr::case_when(
      Direction == "Higher in more severe" ~ "higher",
      Direction == "Lower in more severe"  ~ "lower",
      TRUE                                 ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(canonical_ct), !is.na(direction))

n_excluded <- nrow(dcopa_raw) - nrow(dcopa)
message("  Excluded ", n_excluded, " rows (unmapped Liu celltypes or unknown direction)")
message("  Retained ", nrow(dcopa), " rows across ",
        length(unique(dcopa$Comparison)), " comparisons")

# ============================================================
# BUILD GENE SETS PER COMPARISON
# A gene is counted once regardless of how many celltypes/directions it appears in.
# ============================================================
gene_sets <- setNames(lapply(SET_KEYS, function(comp) {
  unique(dcopa$Gene[dcopa$Comparison == comp])
}), SET_KEYS)

all_ids <- unique(unlist(gene_sets))
message("Total unique genes across all comparisons: ", length(all_ids))

# ============================================================
# BUILD MEMBERSHIP MATRIX
# ============================================================
mat <- as.data.frame(
  setNames(
    lapply(gene_sets, function(genes) as.integer(all_ids %in% genes)),
    SET_KEYS
  )
)

n_genes  <- length(all_ids)
n_shared <- sum(rowSums(mat) > 1)
message(n_genes, " unique genes; ", n_shared, " shared in >=2 comparisons")
message("Set sizes:"); print(colSums(mat))
message("Row sum distribution:"); print(table(rowSums(mat)))

stopifnot(nrow(mat) > 0, any(mat > 0))

# Rename columns for display: remove AllADVsCon_, replace ROSMAP with AD_modules
plot_keys <- gsub("_AllADVsCon_", "_", SET_KEYS)
plot_keys <- gsub("_ROSMAP$", "_AD_modules", plot_keys)
names(mat) <- plot_keys

# Reorder so Gabitto is above Liu (UpSetR renders sets bottom-to-top,
# so Liu must come before Gabitto in the vector)
plot_order <- c(
  "Liu_DFC_AD_modules",     "Gabitto_DFC_AD_modules",
  "Liu_MTG_AD_modules",     "Gabitto_MTG_AD_modules",
  "Liu_DFC",                "Gabitto_DFC",
  "Liu_MTG",                "Gabitto_MTG"
)
mat <- mat[, plot_order]

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
    sets            = plot_order,
    keep.order      = TRUE,
    order.by        = "freq",
    decreasing      = TRUE,
    mb.ratio        = c(0.6, 0.4),
    text.scale      = c(1.3, 1.2, 1, 1, 1.2, 1),
    point.size      = 2.5,
    line.size       = 0.8,
    show.numbers    = "yes",
    mainbar.y.label = "Intersection size",
    sets.x.label    = "Set size",
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

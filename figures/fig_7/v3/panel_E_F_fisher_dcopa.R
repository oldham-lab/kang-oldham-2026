# Pairwise Fisher's exact tests for dCoPA gene set overlap
#
# For every pair of the 8 groups × every celltype × direction (higher/lower),
# tests whether the gene overlap is greater than expected by chance.
#
# Universe: set UNIVERSE_PATH to a plain-text file (one gene per line) to use a
# fixed background (e.g. all expressed genes). If NULL, falls back to the union
# of all genes present across the 8 groups for each (celltype, direction).
#
# FDR correction: Benjamini-Hochberg across all tests.
#
# Outputs:
#   panel_B_fisher_dcopa_full.csv    — one row per pair × celltype × direction
#   panel_B_fisher_dcopa_summary.csv — one row per pair, summarised across celltypes

library(data.table)
library(tidyverse)

# ============================================================
# PLACEHOLDERS
# ============================================================
DCOPA_PATH    <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3/panel_B_dcopa_genelist.csv")
SAVE_DIR      <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3")

# Optional: path to a plain-text file of background genes (one gene per line).
# Set to NULL to use the union of all genes across the 8 groups per celltype×direction.
UNIVERSE_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v3/all_bulk_genes.txt")

# Minimum genes required in BOTH groups for a test to run
MIN_GENES     <- 5

# ============================================================
# SET ORDER  (must match panel_B_upset_dcopa.R)
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
  c("L6 IT",        "Exc L5-6 IT", "Exc L6 IT"),
  c("L6 IT Car3",   "Exc L5/6 IT Car3"),
  c("L6 CT",        "Exc L6 CT"),
  c("Pax6",         "Inh PAX6"),
  c("Astrocyte",    "Ast"),
  c("OPC",          "OPC"),
  c("Vip",          "Inh VIP"),
  c("L6b",          "Exc L6b"),
  c("L2/3 IT",      "Exc L2-3 IT"),
  c("Microglia-PVM", "Mic"),
  c("Oligodendrocyte", "Oli")
)

# Build order-independent canonical lookup after data is loaded (see below).
# For each synonym group, whichever member appears in the Gabitto data becomes
# the canonical name; all other members map to it.
build_canonical_map <- function(map_list, gabitto_celltypes) {
  lookup <- c()
  for (grp in map_list) {
    in_gab    <- grp[grp %in% gabitto_celltypes]
    canonical <- if (length(in_gab) >= 1) in_gab[1] else grp[1]
    lookup    <- c(lookup, setNames(rep(canonical, length(grp)), grp))
  }
  lookup
}

# ============================================================
# LOAD & NORMALISE
# ============================================================
message("Loading dCoPA gene list ...")
dcopa_raw <- fread(DCOPA_PATH, data.table = FALSE)

canonical_map <- build_canonical_map(
  map_list,
  unique(dcopa_raw$Celltype[dcopa_raw$Dataset == "Gabitto"])
)

dcopa <- dcopa_raw |>
  dplyr::mutate(
    canonical_ct = dplyr::case_when(
      Celltype %in% names(canonical_map) ~ unname(canonical_map[Celltype]),
      Dataset == "Gabitto"               ~ Celltype,   # unmapped Gabitto names kept
      TRUE                               ~ NA_character_
    ),
    direction = dplyr::case_when(
      Direction == "Higher in more severe" ~ "higher",
      Direction == "Lower in more severe"  ~ "lower",
      TRUE                                 ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(canonical_ct), !is.na(direction), Comparison %in% SET_KEYS)

message("  ", nrow(dcopa), " rows retained after normalisation")

# ============================================================
# BUILD GENE SETS
# gene_sets[[comparison]][[celltype]][[direction]] = character vector
# ============================================================
all_cts  <- sort(unique(dcopa$canonical_ct))
all_dirs <- c("higher", "lower")

gene_sets <- setNames(lapply(SET_KEYS, function(comp) {
  sub <- dcopa[dcopa$Comparison == comp, ]
  setNames(lapply(all_cts, function(ct) {
    setNames(lapply(all_dirs, function(dir) {
      unique(sub$Gene[sub$canonical_ct == ct & sub$direction == dir])
    }), all_dirs)
  }), all_cts)
}), SET_KEYS)

# Universe: either a user-supplied gene list or the union across all 8 groups.
# NOTE: using the pairwise union would make d = 0 always, causing p = 1.
if (!is.null(UNIVERSE_PATH)) {
  user_universe <- unique(readLines(UNIVERSE_PATH))
  user_universe <- user_universe[nzchar(user_universe)]
  message("  User-supplied universe: ", length(user_universe), " genes from ", UNIVERSE_PATH)
  # Same vector used for every (celltype, direction) test
  global_universe <- setNames(lapply(all_cts, function(ct) {
    setNames(lapply(all_dirs, function(dir) user_universe), all_dirs)
  }), all_cts)
} else {
  message("  No UNIVERSE_PATH set — using union of all 8 groups per celltype×direction")
  global_universe <- setNames(lapply(all_cts, function(ct) {
    setNames(lapply(all_dirs, function(dir) {
      unique(unlist(lapply(SET_KEYS, function(comp) gene_sets[[comp]][[ct]][[dir]])))
    }), all_dirs)
  }), all_cts)
}

# ============================================================
# PAIRWISE FISHER'S EXACT TESTS
# ============================================================
message("Running pairwise Fisher's tests ...")

pairs <- combn(SET_KEYS, 2, simplify = FALSE)
message("  ", length(pairs), " pairs × ", length(all_cts),
        " celltypes × ", length(all_dirs), " directions = ",
        length(pairs) * length(all_cts) * length(all_dirs), " possible tests")

run_fisher <- function(genes_A, genes_B, universe) {
  n_univ <- length(universe)
  if (n_univ == 0) return(NULL)
  if (length(genes_A) < MIN_GENES || length(genes_B) < MIN_GENES) return(NULL)

  a <- length(intersect(genes_A, genes_B))
  b <- length(genes_A) - a                 # in A, not B
  c <- length(genes_B) - a                 # in B, not A
  d <- n_univ - a - b - c                  # in neither (> 0 with global universe)

  ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")

  data.frame(
    n_G1       = length(genes_A),
    n_G2       = length(genes_B),
    n_overlap  = a,
    n_universe = n_univ,
    odds_ratio = unname(ft$estimate),
    p_value    = ft$p.value,
    stringsAsFactors = FALSE
  )
}

results <- lapply(pairs, function(pair) {
  g1 <- pair[1]; g2 <- pair[2]
  rows <- lapply(all_cts, function(ct) {
    lapply(all_dirs, function(dir) {
      res <- run_fisher(gene_sets[[g1]][[ct]][[dir]],
                        gene_sets[[g2]][[ct]][[dir]],
                        global_universe[[ct]][[dir]])
      if (is.null(res)) return(NULL)
      cbind(data.frame(Group1 = g1, Group2 = g2,
                       Celltype = ct, Direction = dir,
                       stringsAsFactors = FALSE), res)
    })
  })
  dplyr::bind_rows(Filter(Negate(is.null), unlist(rows, recursive = FALSE)))
})

full_table <- dplyr::bind_rows(results)

# BH correction across all tests
full_table$p_adj <- p.adjust(full_table$p_value, method = "BH")

message("  ", nrow(full_table), " tests performed")
message("  ", sum(full_table$p_adj < 0.05, na.rm = TRUE),
        " significant at FDR < 0.05")

# ============================================================
# SUMMARY TABLE: per pair, summarised across celltypes × directions
# ============================================================
summary_table <- full_table |>
  dplyr::group_by(Group1, Group2) |>
  dplyr::summarise(
    n_tests              = dplyr::n(),
    n_sig_fdr05          = sum(p_adj < 0.05, na.rm = TRUE),
    mean_neglog10_padj   = mean(-log10(p_adj + 1e-300), na.rm = TRUE),
    median_neglog10_padj = median(-log10(p_adj + 1e-300), na.rm = TRUE),
    median_odds_ratio    = median(odds_ratio, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# SAVE
# ============================================================
if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)

fwrite(full_table,    file = file.path(SAVE_DIR, "panel_B_fisher_dcopa_full.csv"))
fwrite(summary_table, file = file.path(SAVE_DIR, "panel_B_fisher_dcopa_summary.csv"))
message("Saved: panel_B_fisher_dcopa_full.csv  (", nrow(full_table), " rows)")
message("Saved: panel_B_fisher_dcopa_summary.csv  (", nrow(summary_table), " rows)")

# v2.3: changed so output is only SEA celltypes

library(qs)
library(data.table)
library(tidyverse)

# ============================================================
# PLACEHOLDERS
# ============================================================

MIT_EDGER_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_edgeR_ADvsCon_by_celltype.RDS")
MIT_DESEQ_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "MIT_DE_results/mit_DESeq2_ADvsCon_by_celltype.qs")
SEA_EDGER_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_edgeR_ADvsCon_by_celltype.RDS")
SEA_DESEQ_PATH <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "SEA_DE_results/sea_DESeq2_ADvsCon_by_celltype.qs")

DCOPA_PATH <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v2/dcopa_overlap.csv")
SAVE_DIR   <- file.path(Sys.getenv("DE_DIR", "/home/gugene/code/git/kang-oldham-2026/figures/fig_7/DE_old"), "dcopa_overlap/")

# DE thresholds
FDR_THRESH <- 0.05
LFC_THRESH <- 0      # |logFC| > 1

# ============================================================
# CELLTYPE MAPPING
# ============================================================
# Each element maps one dcopa celltype name (first element) to one or more
# DE result celltype names (remaining elements). When matching, DE genes from
# all mapped DE celltypes are pooled before computing overlap with the dcopa genes.
#
# Set DATASET_FOR_MAPPING to the DE dataset whose celltype names the map applies to.
# Other datasets will be matched by exact name only (no mapping applied).
DATASET_FOR_MAPPING <- "MIT_edgeR"   # one of: MIT_edgeR, MIT_DESeq2, SEA_edgeR, SEA_DESeq2
                                      # (mapping is applied to both MIT datasets together
                                      #  if you set either MIT variant)

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

# Build lookup: dcopa_name -> vector of DE celltype names (for the mapped dataset)
# and reverse: DE celltype name -> dcopa name (for diagnostics)
dcopa_to_de  <- setNames(lapply(map_list, function(x) x[-1]), sapply(map_list, `[`, 1))
de_to_dcopa  <- unlist(lapply(map_list, function(x) setNames(rep(x[1], length(x)-1), x[-1])))

# Lookup: any celltype name (MIT or dcopa) -> canonical SEA name.
# SEA names map to themselves; MIT/dcopa names map via the map_list.
sea_names    <- sapply(map_list, `[`, 1)           # first element = SEA name
to_sea_name  <- c(
  setNames(sea_names, sea_names),                  # SEA -> SEA (identity)
  de_to_dcopa                                       # MIT -> SEA (via map)
  # dcopa names ARE the SEA names (first element), so already covered above
)

# Helper: given a dcopa celltype name and a DE result list, return a pooled
# data frame of all DE results from the mapped DE celltypes (or exact match
# if no mapping exists for this dataset).
# Uses data.table for the pooling/dedup step — much faster than dplyr groupby.
get_de_for_ct <- function(dcopa_ct, de_list, use_mapping) {
  if (use_mapping && dcopa_ct %in% names(dcopa_to_de)) {
    mapped_cts <- dcopa_to_de[[dcopa_ct]]
    parts <- Filter(Negate(is.null), lapply(mapped_cts, function(de_ct) de_list[[de_ct]]))
    if (length(parts) == 0) return(NULL)
    # Pool and keep lowest FDR per gene using data.table (faster than dplyr groupby)
    dt <- data.table::rbindlist(parts)
    dt[dt[, .I[which.min(FDR)], by = gene]$V1]
  } else {
    de_list[[dcopa_ct]]
  }
}

# ============================================================
# LOAD & PREPARE DCOPA GENE LIST
# ============================================================
message("Loading dcopa gene list ...")

dcopa_raw <- fread(DCOPA_PATH, data.table = FALSE)

# Filter to AllADVsCon_DFC and celltypes with "all" suffix only
dcopa <- dcopa_raw |>
  dplyr::filter(
    Comparison == "AllADVsCon_DFC",
    grepl(" all$", Celltype)
  ) |>
  dplyr::mutate(
    # Strip " all" suffix so celltype names match DE result keys
    celltype_clean = sub(" all$", "", Celltype),
    direction      = dplyr::case_when(
      Direction == "Higher in more severe" ~ "higher",
      Direction == "Lower in more severe"  ~ "lower",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(direction))

message("  dcopa: ", length(unique(dcopa$celltype_clean)), " celltypes after filtering")

# Split into per-direction, per-celltype gene vectors
dcopa_higher <- split(dcopa$Gene[dcopa$direction == "higher"],
                      dcopa$celltype_clean[dcopa$direction == "higher"])
dcopa_lower  <- split(dcopa$Gene[dcopa$direction == "lower"],
                      dcopa$celltype_clean[dcopa$direction == "lower"])
dcopa_higher <- lapply(dcopa_higher, unique)
dcopa_lower  <- lapply(dcopa_lower,  unique)

# ============================================================
# LOAD & STANDARDISE DE RESULTS
# ============================================================
message("Loading DE results ...")

# Standardise column names to: gene, logFC, FDR
tidy_edger <- function(res_list) {
  lapply(res_list, function(df) {
    df |>
      dplyr::rename(
        gene  = any_of(c("genes", "gene", "Gene")),
        logFC = any_of(c("logFC")),
        FDR   = any_of(c("FDR", "adj.P.Val"))
      ) |>
      dplyr::select(gene, logFC, FDR)
  })
}

tidy_deseq <- function(res_list) {
  lapply(res_list, function(df) {
    df |>
      dplyr::rename(
        gene  = any_of(c("genes", "gene", "Gene")),
        logFC = any_of(c("logFC", "log2FoldChange")),
        FDR   = any_of(c("FDR", "padj"))
      ) |>
      dplyr::select(gene, logFC, FDR)
  })
}

mit_edger <- tidy_edger(readRDS(MIT_EDGER_PATH))
mit_deseq <- tidy_deseq(qread(MIT_DESEQ_PATH))
sea_edger <- tidy_edger(readRDS(SEA_EDGER_PATH))
sea_deseq <- tidy_deseq(qread(SEA_DESEQ_PATH))

# Determine which datasets get the celltype mapping applied
# (both MIT variants share the same celltype names, so map both if either is specified)
mit_datasets <- c("MIT_edgeR", "MIT_DESeq2")
sea_datasets <- c("SEA_edgeR", "SEA_DESeq2")
mapped_group <- if (DATASET_FOR_MAPPING %in% mit_datasets) mit_datasets else sea_datasets

message("\n--- Celltype name check ---")
message("Mapping applied to: ", paste(mapped_group, collapse = ", "))
message("dcopa celltypes (after stripping ' all'): ",
        paste(sort(unique(dcopa$celltype_clean)), collapse = ", "))
message("MIT edgeR celltypes: ", paste(sort(names(mit_edger)), collapse = ", "))
message("SEA edgeR celltypes: ", paste(sort(names(sea_edger)), collapse = ", "))

dcopa_cts  <- unique(dcopa$celltype_clean)
mapped_de_cts <- unique(unlist(dcopa_to_de))  # all DE names covered by the map

# Check coverage for the mapped datasets
mit_covered <- intersect(names(mit_edger), mapped_de_cts)
sea_covered <- intersect(names(sea_edger), names(dcopa_cts))  # exact match for unmapped

message("\nMIT DE celltypes covered by map: ", paste(sort(mit_covered), collapse = ", "))
mit_unmapped <- setdiff(names(mit_edger), mapped_de_cts)
if (length(mit_unmapped) > 0)
  message("MIT DE celltypes NOT in map (will be matched exactly): ",
          paste(sort(mit_unmapped), collapse = ", "))

dcopa_no_map <- setdiff(dcopa_cts, names(dcopa_to_de))
if (length(dcopa_no_map) > 0)
  message("dcopa celltypes with no map entry (exact match only): ",
          paste(sort(dcopa_no_map), collapse = ", "))
message("---\n")

# ============================================================
# HELPER: extract directional DE gene sets from one DE result
# ============================================================
# Uses base R subsetting — avoids dplyr overhead for a simple filter + split.
get_de_directional <- function(de_df, fdr_thresh, lfc_thresh) {
  keep <- !is.na(de_df$FDR) & !is.na(de_df$logFC) &
          de_df$FDR < fdr_thresh & abs(de_df$logFC) > lfc_thresh
  sig  <- de_df[keep, ]
  list(
    higher = unique(sig$gene[sig$logFC > 0]),
    lower  = unique(sig$gene[sig$logFC < 0])
  )
}

# ============================================================
# MAIN ANALYSIS
# ============================================================
de_sets <- list(
  MIT_edgeR  = mit_edger,
  MIT_DESeq2 = mit_deseq,
  SEA_edgeR  = sea_edger,
  SEA_DESeq2 = sea_deseq
)

# All celltypes present across dcopa and any DE result, expressed as SEA names.
# MIT celltypes are translated via to_sea_name; any name not in the map is kept
# as-is (e.g. a SEA celltype that has no MIT counterpart).
all_cts_raw <- unique(c(
  names(dcopa_higher), names(dcopa_lower),
  names(mit_edger), names(mit_deseq),
  names(sea_edger), names(sea_deseq)
))
all_cts <- sort(unique(ifelse(all_cts_raw %in% names(to_sea_name),
                              to_sea_name[all_cts_raw],
                              all_cts_raw)))

# ── Pre-compute directional gene sets once for every ct x dataset ────────────
# This is the key optimisation: get_de_for_ct (which does pooling/dedup) and
# get_de_directional (which filters) are each called only once per combination
# instead of once per table. Results are stored in a nested list:
#   de_cache[[ds_name]][[ct]] = list(higher = char_vec, lower = char_vec)
message("Pre-computing DE gene sets ...")
de_cache <- lapply(names(de_sets), function(ds_name) {
  use_mapping <- ds_name %in% mapped_group
  setNames(lapply(all_cts, function(ct) {
    de_df <- get_de_for_ct(ct, de_sets[[ds_name]], use_mapping)
    if (is.null(de_df)) return(NULL)
    get_de_directional(de_df, FDR_THRESH, LFC_THRESH)
  }), all_cts)
})
names(de_cache) <- names(de_sets)

# ── TABLE 1: Summary counts and % overlap per celltype ───────────────────────
message("Building summary table ...")

summary_rows <- lapply(all_cts, function(ct) {

  n_dcopa_higher <- length(dcopa_higher[[ct]])
  n_dcopa_lower  <- length(dcopa_lower[[ct]])

  row <- data.frame(celltype = ct,
                    n_dcopa_higher = n_dcopa_higher,
                    n_dcopa_lower  = n_dcopa_lower)

  for (ds_name in names(de_sets)) {
    cached <- de_cache[[ds_name]][[ct]]

    if (is.null(cached)) {
      row[[paste0("n_DE_higher_", ds_name)]] <- NA_integer_
      row[[paste0("n_DE_lower_",  ds_name)]] <- NA_integer_
      row[[paste0("pct_overlap_higher_", ds_name)]] <- NA_real_
      row[[paste0("pct_overlap_lower_",  ds_name)]] <- NA_real_
    } else {
      de_higher <- cached$higher
      de_lower  <- cached$lower
      dcopa_h   <- dcopa_higher[[ct]]
      dcopa_l   <- dcopa_lower[[ct]]

      n_overlap_higher <- length(intersect(de_higher, dcopa_h))
      n_overlap_lower  <- length(intersect(de_lower,  dcopa_l))

      pct_h <- if (length(dcopa_h) > 0) round(100 * n_overlap_higher / length(dcopa_h), 1) else NA_real_
      pct_l <- if (length(dcopa_l) > 0) round(100 * n_overlap_lower  / length(dcopa_l), 1) else NA_real_

      row[[paste0("n_DE_higher_", ds_name)]] <- length(de_higher)
      row[[paste0("n_DE_lower_",  ds_name)]] <- length(de_lower)
      row[[paste0("pct_overlap_higher_", ds_name)]] <- pct_h
      row[[paste0("pct_overlap_lower_",  ds_name)]] <- pct_l
    }
  }
  row
})

summary_table <- dplyr::bind_rows(summary_rows)

# Reorder columns for readability
count_cols   <- as.vector(outer(c("n_DE_higher_", "n_DE_lower_"), names(de_sets), paste0))
dcopa_cols   <- c("n_dcopa_higher", "n_dcopa_lower")
overlap_cols <- as.vector(outer(c("pct_overlap_higher_", "pct_overlap_lower_"), names(de_sets), paste0))

summary_table <- summary_table |>
  dplyr::select(celltype, all_of(count_cols), all_of(dcopa_cols), all_of(overlap_cols))

# ── TABLE 2: Overlapping gene names, long format ─────────────────────────────
message("Building overlap gene table ...")

overlap_rows <- lapply(all_cts, function(ct) {
  lapply(names(de_sets), function(ds_name) {
    cached <- de_cache[[ds_name]][[ct]]
    if (is.null(cached)) return(NULL)

    higher_genes <- intersect(cached$higher, dcopa_higher[[ct]])
    lower_genes  <- intersect(cached$lower,  dcopa_lower[[ct]])

    rows <- list()
    if (length(higher_genes) > 0)
      rows[[1]] <- data.frame(celltype  = ct, dataset = ds_name,
                              direction = "higher_in_AD", gene = higher_genes)
    if (length(lower_genes) > 0)
      rows[[2]] <- data.frame(celltype  = ct, dataset = ds_name,
                              direction = "lower_in_AD",  gene = lower_genes)
    dplyr::bind_rows(rows)
  }) |> dplyr::bind_rows()
}) |> dplyr::bind_rows()

if (nrow(overlap_rows) == 0) {
  warning("No overlapping genes found between DE results and dcopa list. ",
          "Check that celltype names match (see diagnostic output above).")
  overlap_table <- data.frame(celltype  = character(),
                               dataset   = character(),
                               direction = character(),
                               gene      = character())
} else {
  overlap_table <- overlap_rows |>
    dplyr::arrange(celltype, direction, dataset, gene)
}

# ============================================================
# SAVE OUTPUTS
# ============================================================
if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)

# All celltype labels in outputs use the canonical SEA-AD name.
# MIT and dcopa celltypes have been translated via to_sea_name before analysis.
fwrite(summary_table,
       file = file.path(SAVE_DIR, "dcopa_overlap_summary.csv"))
message("Saved: dcopa_overlap_summary.csv")

fwrite(overlap_table,
       file = file.path(SAVE_DIR, "dcopa_overlap_genes.csv"))
message("Saved: dcopa_overlap_genes.csv")

# ============================================================
# PRINT PREVIEW
# ============================================================
message("\n--- Summary table (", nrow(summary_table), " celltypes) ---")
print(as_tibble(summary_table), n = Inf, width = Inf)

message("\n--- Overlap genes (", nrow(overlap_table), " entries across all datasets/celltypes) ---")
print(head(overlap_table, 30))

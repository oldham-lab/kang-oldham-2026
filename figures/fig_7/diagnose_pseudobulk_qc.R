# Diagnostic QC: compare SEA-AD pseudobulk data across regions
#
# Checks per dataset:
#   1. Sample counts and group balance
#   2. Library size distribution (colSums)
#   3. Count integrality (are values whole numbers?)
#   4. Per-gene count range (flag possible normalised data)
#   5. Samples per celltype per group
#
# Run this before interpreting differences in DE gene counts between regions.

library(data.table)
library(tidyverse)

# ============================================================
# PATHS
# ============================================================
datasets <- list(
  SEA_DFC = list(
    expr = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/cell_expression_by_donor_sum_subclass.csv"),
    anno = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/annotations_by_donor_sum_subclass.csv"),
    celltype_col = "Subclass",
    dx_col       = "Dx",
    dx_control   = "Control",
    dx_ad        = "Alzheimers"
  ),
  SEA_MTG = list(
    expr = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/mtg_data_pseudobulk/sea_mtg_pseudobulk.csv.gz"),
    anno = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/mtg_data_pseudobulk/sea_mtg_pseudobulk_annotations.csv"),
    celltype_col = "Subclass",
    dx_col       = "Dx",
    dx_control   = "Control",
    dx_ad        = "Alzheimers"
  )
)

# ============================================================
# HELPERS
# ============================================================
check_integrality <- function(mat, n_genes_sample = 200, n_samples_sample = 20) {
  # Sample a subset to keep this fast on large matrices
  row_idx <- sample(nrow(mat), min(n_genes_sample, nrow(mat)))
  col_idx <- sample(ncol(mat), min(n_samples_sample, ncol(mat)))
  sub     <- mat[row_idx, col_idx]
  frac_integer <- mean(sub == floor(sub), na.rm = TRUE)
  list(
    frac_integer = frac_integer,
    max_value    = max(mat, na.rm = TRUE),
    min_value    = min(mat, na.rm = TRUE),
    any_negative = any(mat < 0, na.rm = TRUE)
  )
}

# ============================================================
# RUN CHECKS
# ============================================================
for (ds_name in names(datasets)) {
  cfg <- datasets[[ds_name]]
  cat("\n", strrep("=", 60), "\n")
  cat("DATASET:", ds_name, "\n")
  cat(strrep("=", 60), "\n\n")

  # --- Load ---
  cat("Loading expression matrix ...\n")
  expr_raw <- fread(cfg$expr, data.table = FALSE)

  # First column is gene names
  genes <- expr_raw[, 1]
  mat   <- as.matrix(expr_raw[, -1])
  rownames(mat) <- genes
  cat("  Dimensions:", nrow(mat), "genes x", ncol(mat), "samples\n")

  cat("Loading annotations ...\n")
  anno <- fread(cfg$anno, data.table = FALSE)
  cat("  Annotation rows:", nrow(anno), "  Columns:", paste(names(anno), collapse = ", "), "\n")

  # --- Check sample alignment ---
  expr_samples <- colnames(mat)
  anno_samples <- anno$label %||% anno[[1]]   # use first col if no 'label' col
  n_match <- sum(expr_samples %in% anno_samples)
  cat("\n--- Sample alignment ---\n")
  cat("  Expr samples:", length(expr_samples),
      " | Anno rows:", nrow(anno),
      " | Matching:", n_match, "\n")
  if (n_match < length(expr_samples))
    cat("  WARNING: not all expression samples are in the annotation file\n")

  # --- Group balance ---
  cat("\n--- Group balance ---\n")
  if (cfg$dx_col %in% names(anno)) {
    print(table(anno[[cfg$dx_col]]))
  } else {
    cat("  WARNING: dx_col '", cfg$dx_col, "' not found in annotations\n", sep = "")
    cat("  Available columns:", paste(names(anno), collapse = ", "), "\n")
  }

  # --- Celltypes ---
  cat("\n--- Celltypes ---\n")
  if (cfg$celltype_col %in% names(anno)) {
    ct_counts <- sort(table(anno[[cfg$celltype_col]]), decreasing = TRUE)
    cat("  ", length(ct_counts), "celltypes\n")
    print(ct_counts)
  } else {
    cat("  WARNING: celltype_col '", cfg$celltype_col, "' not found\n", sep = "")
  }

  # --- Samples per celltype per group ---
  if (cfg$celltype_col %in% names(anno) && cfg$dx_col %in% names(anno)) {
    cat("\n--- Samples per celltype per group ---\n")
    ct_dx <- table(anno[[cfg$celltype_col]], anno[[cfg$dx_col]])
    print(ct_dx)

    low_n <- ct_dx < 3
    if (any(low_n)) {
      cat("  WARNING: celltype×group combinations with <3 samples (will be skipped in DE):\n")
      idx <- which(low_n, arr.ind = TRUE)
      for (i in seq_len(nrow(idx)))
        cat("   ", rownames(ct_dx)[idx[i, 1]], "/", colnames(ct_dx)[idx[i, 2]],
            ":", ct_dx[idx[i, 1], idx[i, 2]], "samples\n")
    }
  }

  # --- Library sizes ---
  cat("\n--- Library sizes (colSums) ---\n")
  lib_sizes <- colSums(mat)
  cat("  Min:", format(min(lib_sizes),  big.mark = ","),
      " Median:", format(median(lib_sizes), big.mark = ","),
      " Max:", format(max(lib_sizes),  big.mark = ","), "\n")
  cat("  Ratio max/min:", round(max(lib_sizes) / max(min(lib_sizes), 1), 1), "\n")
  if (sd(lib_sizes) / mean(lib_sizes) < 0.01)
    cat("  WARNING: library sizes are suspiciously uniform — data may be normalised\n")

  # --- Count integrality ---
  cat("\n--- Count integrality (sampled) ---\n")
  ic <- check_integrality(mat)
  cat("  Fraction of sampled values that are whole numbers:", round(ic$frac_integer, 4), "\n")
  cat("  Value range: [", ic$min_value, ",", format(ic$max_value, big.mark = ","), "]\n")
  if (ic$any_negative)
    cat("  WARNING: negative values found — data may be log-transformed\n")
  if (ic$frac_integer < 0.99)
    cat("  WARNING: many non-integer values — data may be normalised/log-transformed\n")
  if (ic$max_value < 50)
    cat("  WARNING: max value < 50 — data is almost certainly normalised\n")

  # --- Float32 precision risk ---
  cat("\n--- Float32 precision check ---\n")
  # float32 represents integers exactly up to 2^24 = 16,777,216
  float32_limit <- 2^24
  max_val       <- ic$max_value
  cat("  Max count:", format(max_val, big.mark = ","),
      " | float32 exact-integer limit:", format(float32_limit, big.mark = ","), "\n")
  if (max_val > float32_limit)
    cat("  WARNING: max count exceeds float32 precision limit.",
        "Pseudobulk sums accumulated as float32 may have lost precision.\n",
        "  Consider re-running prepare_seaad_mtg_pseudobulk.py with float64 accumulation.\n")
  else
    cat("  OK: max count is within float32 exact-integer range\n")

  # --- Per-gene max count (flag very high counts) ---
  cat("\n--- Top 10 genes by max count across samples ---\n")
  gene_max <- sort(apply(mat, 1, max), decreasing = TRUE)
  print(head(gene_max, 10))

  cat("\n")
}

cat(strrep("=", 60), "\n")
cat("QC complete.\n")

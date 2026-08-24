library(data.table)
library(tibble)

# ============================================================================
# Per-module percentage of genes with a positive control - case difference.
#
# Loops over the same comparisons as the three *AllComparisons*_pval_runcode.R
# drivers (SEAAD, MIT multiome, brainSCOPE), but instead of computing
# euclidean-distance p-values it reports, per module and per cell type:
#
#   100 * (# genes where (control - case) > 0) / (# genes with non-zero diff)
#
# i.e. the fraction of a module's genes that are higher in control than in
# case, computed gene-by-gene on the genome-wide cell-type means. Ties
# (control == case) are excluded from both numerator and denominator.
#
# Output: module x cell-type CSVs, same structure as the old
# *_projdistpvalindiv.csv files, written to the subfolder defined below.
# ============================================================================

filter_under <- 3   # keep modules with > filter_under member genes (as before)

out_subdir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13/projdistpvalindiv/pct_positive_signs")

# ── Module-definition directories (shared across datasets) ──────────────────
mod_dirs <- list(
  bulk_megaset = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"),
  rosmap       = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC"),
  brainseq_scz = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/")
)

sigcount_bonf_file <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv")

# ── Dataset configurations (the three drivers, merged) ──────────────────────
# means_dir is a function of (base_dir, region) because brainSCOPE nests a
# "DFC" level that the AD datasets do not.
datasets <- list(
  list(
    prefix    = "SEAAD_",
    base_dir  = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output"),
    regions   = c("PFC", "MTC"),
    mod_types = c("bulk_megaset", "rosmap"),
    runs      = list(list(case = "allAD", control = "Con", label = "allAD_vs_Con")),
    means_dir = function(base, region) file.path(base, region, "means")
  ),
  list(
    prefix    = "MIT_",
    base_dir  = file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output"),
    regions   = c("PFC", "MTC"),
    mod_types = c("bulk_megaset", "rosmap"),
    runs      = list(list(case = "allAD", control = "Con", label = "allAD_vs_Con")),
    means_dir = function(base, region) file.path(base, region, "means")
  ),
  list(
    prefix    = "brainSCOPE_",
    base_dir  = file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output"),
    regions   = c("CMC", "SZBDMulti-Seq"),
    mod_types = c("bulk_megaset", "brainseq_scz"),
    runs      = list(list(case = "Schizophrenia", control = "control", label = "Schizophrenia_vs_control")),
    means_dir = function(base, region) file.path(base, region, "DFC", "means")
  )
)

# ── Module membership + filtering (same selection as the pval drivers) ──────
# Returns a named list of gene-vectors, keyed by the kept module ids (in
# numeric order), matching the row order of the old *_projdistpvalindiv.csv.
get_module_genes <- function(mod_type) {
  datkme <- fread(data.table = F, file = file.path(mod_dirs[[mod_type]], "kme_tables", "topmodposbc_table.csv"))
  mods   <- split(datkme[, 2], datkme[, 3])               # gene lists keyed by module id
  these  <- as.numeric(names(mods)[lengths(mods) > filter_under])

  if (mod_type == "bulk_megaset") {
    sigcount_bonf <- fread(data.table = F, file = sigcount_bonf_file)
    these <- these[!these %in% which(sigcount_bonf$vals < 2)]
  }

  mods[as.character(these)]
}

# ── Genome-wide cell-type means for one group (NA -> 0, as before) ──────────
load_means <- function(means_dir, group) {
  out <- fread(data.table = F, file = file.path(means_dir, paste0("genomewide_means_", group, ".csv")))
  out <- tibble::column_to_rownames(out, names(out)[1])
  out[is.na(out)] <- 0
  return(out)
}

if (!dir.exists(out_subdir)) dir.create(out_subdir, recursive = TRUE)

# ── Main loop ───────────────────────────────────────────────────────────────
for (ds in datasets) {
  for (region in ds$regions) {
    means_dir <- ds$means_dir(ds$base_dir, region)

    for (run in ds$runs) {
      # Genome-wide means files for this comparison (skip if either is missing)
      case_file    <- file.path(means_dir, paste0("genomewide_means_", run$case,    ".csv"))
      control_file <- file.path(means_dir, paste0("genomewide_means_", run$control, ".csv"))
      missing <- !file.exists(c(case_file, control_file))
      if (any(missing)) {
        cat(sprintf("Skipping %s %s %s - missing file(s):\n  %s\n",
                    ds$prefix, region, run$label,
                    paste(c(case_file, control_file)[missing], collapse = "\n  ")))
        next
      }

      case_m    <- load_means(means_dir, run$case)
      control_m <- load_means(means_dir, run$control)

      # Align genes (rows) and cell types (cols) between the two groups
      genes <- rownames(case_m)
      cts   <- intersect(colnames(case_m), colnames(control_m))
      diff  <- as.matrix(control_m[genes, cts, drop = FALSE]) -
               as.matrix(case_m[genes, cts, drop = FALSE])      # control - case

      for (mod_type in ds$mod_types) {
        cat(sprintf("\n%s\n%s | %s | %s | %s\n%s\n",
                    strrep("=", 60), ds$prefix, run$label, mod_type, region, strrep("=", 60)))

        module_genes <- get_module_genes(mod_type)

        # Per module: % of present, non-tied genes with (control - case) > 0
        pct <- t(sapply(module_genes, function(g) {
          idx <- which(genes %in% g)
          sub <- diff[idx, , drop = FALSE]
          pos     <- colSums(sub > 0)
          nonzero <- colSums(sub != 0)
          out <- 100 * pos / nonzero          # ties excluded from numerator and denominator
          out[!is.finite(out)] <- NA          # modules with no present/non-tied genes -> NA
          out
        }))
        colnames(pct) <- cts
        pct <- as.data.frame(pct)
        rownames(pct) <- seq_len(nrow(pct))   # match old row indexing

        out_file <- file.path(out_subdir,
                              paste0(ds$prefix, region, "_", run$label, "_", mod_type, "_pctposdiff.csv"))
        write.csv(pct, file = out_file)
        cat(sprintf("  Wrote %s  (%d modules x %d cell types)\n", basename(out_file), nrow(pct), ncol(pct)))
      }
    }
  }
}

cat("\nAll runs complete.\n")

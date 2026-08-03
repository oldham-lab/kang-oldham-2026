library(data.table)
library(qs)
library(future.apply)
library(tibble)
library(dplyr)
library(tidyr)
library(CoPA)

plan(multisession, workers = 8)

# ── Run options ────────────────────────────────────────────────────────────────
# All possible combinations (region x label x mod_type):
# list(list(region = "CMC",           label = "SZ_vs_Con", mod_type = "bulk_megaset"),
#      list(region = "CMC",           label = "SZ_vs_Con", mod_type = "brainseq_scz"),
#      list(region = "SZBDMulti-Seq", label = "SZ_vs_Con", mod_type = "bulk_megaset"),
#      list(region = "SZBDMulti-Seq", label = "SZ_vs_Con", mod_type = "brainseq_scz"))

run_only <- NULL
# NULL runs all combinations

save_randinds_for <- NULL
# NULL disables saving for all runs

# ── Configuration ──────────────────────────────────────────────────────────────
regions  <- c("CMC",
              "SZBDMulti-Seq")

base_dir <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output"

mod_configs <- list(
  bulk_megaset = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"),
  brainseq_scz = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/")
)

expr       <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table = F)
bulk_genes <- expr[, 2]

# ── Module filter (computed once per mod_type from its respective module_output_dir) ──
filter_under <- 3
mod_filters <- lapply(mod_configs, function(mod_dir) {
  datkme <- fread(data.table = F, file = file.path(mod_dir, "kme_tables", "topmodposbc_table.csv"))
  if (sum(duplicated(datkme[, 2])) > 0) datkme[, 2] <- make.unique(datkme[, 2])
  mods          <- tapply(datkme[, 2], datkme[, 3], list)
  modulelengths <- unlist(lapply(mods, length))
  as.numeric(names(mods)[which(modulelengths > filter_under)])
})

# ── Helper ─────────────────────────────────────────────────────────────────────
# means_dir is set inside the region loop before this is called
load_means <- function(group) {
  out <- fread(data.table = F, file = file.path(means_dir, paste0("genomewide_means_", group, ".csv")))
  out <- tibble::column_to_rownames(out, names(out)[1])
  out[is.na(out)] <- 0
  return(out)
}

# ── Run definitions ────────────────────────────────────────────────────────────
runs <- list(
  list(case = "Schizophrenia", control = "control", label = "Schizophrenia_vs_control")
)

mod_types <- names(mod_configs)

# ── Source combined function ───────────────────────────────────────────────────

# ── Main loop ──────────────────────────────────────────────────────────────────
for (region in regions) {
  means_dir    <- file.path(base_dir, region, "DFC", "means")
  modmeans_dir <- file.path(base_dir, region, "DFC", "mod_means", "log_native")

  for (run in runs) {
    for (mod_type in mod_types) {

      if (!is.null(run_only) && !any(sapply(run_only, \(x) region == x$region && run$label == x$label && mod_type == x$mod_type))) next

      module_output_dir <- mod_configs[[mod_type]]
      these_mods        <- mod_filters[[mod_type]]

      cat(sprintf("\n%s\nRun: %s | %s | region: %s\n%s\n",
                  strrep("=", 60), run$label, mod_type, region, strrep("=", 60)))

      # Verify both mod_means files exist before proceeding
      case_mod_file    <- file.path(modmeans_dir, paste0("mod_means_", run$case,    "_", mod_type, ".csv"))
      control_mod_file <- file.path(modmeans_dir, paste0("mod_means_", run$control, "_", mod_type, ".csv"))
      missing <- !file.exists(c(case_mod_file, control_mod_file))
      if (any(missing)) {
        cat("  Skipping — missing file(s):\n")
        cat(paste(" ", c(case_mod_file, control_mod_file)[missing], collapse = "\n"), "\n")
        next
      }

      # sn_objs: per-celltype mean expression matrices for case and control
      sn_objs <- list(
        mean = list(
          Subclass = setNames(
            list(load_means(run$case), load_means(run$control)),
            c(run$case, run$control)
          )
        )
      )

      # proj_all: module-level mean projections for case and control, filtered to these_mods
      proj_all <- setNames(
        list(fread(data.table = F, file = case_mod_file),
             fread(data.table = F, file = control_mod_file)),
        c(run$case, run$control)
      )
      proj_all <- lapply(proj_all, \(x) {
        out <- x[these_mods, ]
        out[is.na(out)] <- 0
        return(out)
      })

      # Output saved to: {base_dir}/{region}/DFC/euclidean_distances/{label}_{mod_type}_output_table.csv
      project_rand_and_calculate_euclidean(
        module_output_dir = module_output_dir,
        filter_under      = filter_under,
        do_log            = TRUE,
        bulk_genes        = bulk_genes,
        save_dir1         = file.path(base_dir, region, "DFC", "euclidean_distances"),
        sn_objs           = sn_objs,
        proj_all          = proj_all,
        rand_n            = 10000,
        seed              = 26,
        out_prefix        = paste0(run$label, "_", mod_type, "_"),
        save_randinds     = !is.null(save_randinds_for) && any(sapply(save_randinds_for, \(x) region == x$region && run$label == x$label && mod_type == x$mod_type))
      )
    }
  }
}

cat("\nAll runs complete.\n")

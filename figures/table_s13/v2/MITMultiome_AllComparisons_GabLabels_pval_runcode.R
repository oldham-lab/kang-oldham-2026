library(data.table)
library(qs)
library(future.apply)
library(tibble)
library(dplyr)
library(tidyr)
library(CoPA)

plan(multisession, workers = 8)

# ── Run options ────────────────────────────────────────────────────────────────
#run_only <- list(list(region = "PFC", label = "allAD_vs_Con", mod_type = "bulk_megaset"))
                  # set to a list of combinations to run, e.g.:
                  # list(list(region = "PFC", label = "allAD_vs_Con", mod_type = "bulk_megaset"))
                  # NULL runs all combinations
run_only <- NULL

#save_randinds_for <- list(region = "PFC", label = "allAD_vs_Con", mod_type = "bulk_megaset")  # set to NULL to disable
save_randinds_for <- NULL

# ── Configuration ──────────────────────────────────────────────────────────────
regions  <- c("PFC",
              "MTC")

base_dir <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output"

mod_configs <- list(
  bulk_megaset = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"),
  rosmap       = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC")
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
  out[is.na(out)] <- 0 # for cell types that have zero cells (see APOE44 and Sst Chodl)
  return(out)
}


# ── Run definitions: 4 comparisons x 2 module types = 8 runs ──────────────────
runs <- list(
  list(case = "allAD",   control = "Con",     label = "allAD_vs_Con")#,
  #list(case = "earlyAD", control = "Con",     label = "earlyAD_vs_Con"),
  #list(case = "lateAD",  control = "earlyAD", label = "lateAD_vs_earlyAD"),
  #list(case = "APOE44",  control = "APOE33",  label = "APOE44_vs_APOE33")
)

mod_types <- names(mod_configs)

# ── Source truncated function ──────────────────────────────────────────────────
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13/project_rand_and_calculate_pval_standalone.R"))

# ── Main loop ──────────────────────────────────────────────────────────────────
for (region in regions) {
  means_dir    <- file.path(base_dir, region, "means")
  modmeans_dir <- file.path(base_dir, region, "mod_means", "log_native")

  for (run in runs) {
    for (mod_type in mod_types) {

      if (!is.null(run_only) && !any(sapply(run_only, \(x) region == x$region && run$label == x$label && mod_type == x$mod_type))) next

      module_output_dir <- mod_configs[[mod_type]]
      these_mods        <- mod_filters[[mod_type]]

      if(mod_type == "bulk_megaset"){
        sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
        these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
      }


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

      # Output saved to: {script_dir}/MIT_{region}_{label}_{mod_type}_projdistpvalindiv.csv
      project_rand_and_calculate_pval(
        module_output_dir = module_output_dir,
        filter_under      = filter_under,
        do_log            = TRUE,
        bulk_genes        = bulk_genes,
        save_dir1         = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s13/v2"),
        sn_objs           = sn_objs,
        proj_all          = proj_all,
        rand_n            = 10000,
        seed              = 26,
        out_prefix        = paste0("MIT_", region, "_", run$label, "_", mod_type, "_"),
        save_randinds     = !is.null(save_randinds_for) && region == save_randinds_for$region && run$label == save_randinds_for$label && mod_type == save_randinds_for$mod_type
      )
    }
  }
}

cat("\nAll runs complete.\n")

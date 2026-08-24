library(data.table)
library(qs)
library(future.apply)
options(future.globals.maxSize=1e9)
plan(multisession, workers=10)
library(CoPA)

######### Input data (SEA-AD 2024)
# Load bulk megaset dataset
expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SCZ_SampleNetworks/1_04-04-29/brainseq_samp_filt_SCZ_1_171_outliers_removed_geneSymbolsAdded.csv"), data.table=F)
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/")
sn_summary_object_path <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/brainSCOPE/")
module_output_dir <- save_dir1
##################################
 
COPA(expr = expr, 
     plot = F,
     save_dir1 = save_dir1,
     sn_summary_object_path = sn_summary_object_path)

# COPA_compare() was removed when the CoPA package was folded into the engine.
# It sourced files that had already moved to deprecated/ and shared state via
# source(local = TRUE). project_rand_and_calculate_euclidean() is the
# explicit-parameter successor: it recomputes the module filter itself but takes
# sn_objs and proj_all directly, so build them here the way the old
# calculate_euclidean_distances() did internally.
filter_under <- 3

mod_bc <- fread(data.table = FALSE,
                file = file.path(module_output_dir, "kme_tables", "topmodposbc_table.csv")) |>
  (\(.) tapply(.[, 2], .[, 3], list))()
these_mods <- which(lengths(mod_bc) > filter_under)

# The sn summary objects sit with the single-nucleus inputs, while the projection
# indices are COPA() output under save_dir1. COPA_compare() used one path for
# both and defaulted it to save_dir1, which has no sn_summary_tables/ - that is
# why this call could not have run as written.
sn_objs <- qread(file = file.path(sn_summary_object_path, "sn_summary_tables",
                                  "sn_summary_objects_log.qs"))
snnames <- mapply(\(ct, ct_name) paste0(ct_name, "_", names(ct)),
                  sn_objs$mean, names(sn_objs$mean), SIMPLIFY = FALSE)
index_files <- list.files(file.path(save_dir1, "sn_proj_indices", "log_native"),
                          full.names = TRUE)
index_files <- index_files[grep("indices_over", index_files)]
proj_all <- lapply(unlist(snnames), \(x)
                   fread(data.table = FALSE,
                         file = index_files[grep(x, index_files)])[these_mods, ])

project_rand_and_calculate_euclidean(module_output_dir = module_output_dir,
                                     filter_under      = filter_under,
                                     do_log            = TRUE,
                                     bulk_genes        = expr[, 2],
                                     save_dir1         = save_dir1,
                                     sn_objs           = sn_objs,
                                     proj_all          = proj_all,
                                     rand_n            = 10000)
 
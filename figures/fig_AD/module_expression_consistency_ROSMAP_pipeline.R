library(DESeq2)
library(qs)
library(data.table)
library(AnnotationHub)
library(tidyverse)
library(cowplot)
options(bitmapType = 'cairo')

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/module_expression_consistency_fxns.R"))
 
#########
# ROSMAP (DFC)
#########

# Load expr (ROSMAP bulk AD samples only)
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_AD_only_SampleNetworks/1_10-52-54/ROSMAP_samp_filt_AD_only_1_248_ComBat.csv"), data.table=F)
megaset_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
expr <- expr[match(megaset_expr[,1], expr[,1]),]
expr <- data.frame("ensembl_id"=expr[,1],"Gene"= megaset_expr[,2], expr[,3:ncol(expr)])
expr <- expr[!is.na(expr[,1]),]
bulk_expr <- expr
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
# Load mods (ROSMAP bulk AD samples only)
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Con vs all AD
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024/")
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/rosmapAD_unnormalized/conVsAllAD")
objs4 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/"),
                        return_obj = T) 

## Con vs early AD
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_conVsEarly/")
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/rosmapAD_unnormalized/conVsEarly")
objs5 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/"),
                        return_obj = T)

## early vs late AD
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_earlyVsLate/")
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/rosmapAD_unnormalized/earlyVsLate")
objs6 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/"),
                        return_obj = T)

###############
# ROSMAP on MIT
###############
# Load expr (ROSMAP bulk AD samples only)
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_AD_only_SampleNetworks/1_10-52-54/ROSMAP_samp_filt_AD_only_1_248_ComBat.csv"), data.table=F)
megaset_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
expr <- expr[match(megaset_expr[,1], expr[,1]),]
expr <- data.frame("ensembl_id"=expr[,1],"Gene"= megaset_expr[,2], expr[,3:ncol(expr)])
expr <- expr[!is.na(expr[,1]),]
bulk_expr <- expr
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
# Load mods (ROSMAP bulk AD samples only)
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Con vs all AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_rosmap/conVAll")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome")
objs1 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/"),
                        return_obj = T)

# Con vs early AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_rosmap/conVEarly/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_conVEarly")
objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly/"),
                        return_obj = T)

# Early vs late AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_rosmap/earlyVLate/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_earlyVLate")
objs3 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate/"),
                        return_obj = T)

# # Con vs late AD
# save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_rosmap/conVLate/")
# copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_conVLate")
# objs4 <- find_cons_mods(save_path = save_path,
#                         copa_dir1 = copa_dir1,
#                         sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVLate/"),
#                         return_obj = T)


############
# Overlap (SEAAD2024 PFC, MIT_Multiome PFC)
##############

# Find overlap
paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/rosmapAD_unnormalized/"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_rosmap/"))
pathnames <- c("bulk_megaset", "bulk_MIT")
common_pool <- list(
  c("Endothelial", "SMC", "VLMC", "End", "Per"),
  c("L4 IT", "Exc L4-5 IT-2", "Exc L3-4 IT","Exc L4-5 IT-1"),
  c("L5 ET",  "Exc L5 ET"),
  c("L5 IT", "Exc L4-5 IT-2", "Exc L4-5 IT-1","Exc L3-5 IT", "Exc L5-6 IT"),
  c("L5/6 NP", "Exc L5/6 NP"),
  c("Lamp5", "Inh LAMP5"),
  c("Pvalb", "Inh PVALB"),
  c("Sst", "Inh SST"),
  c("L6 IT", "Exc L5-6 IT"),
  c("L6 IT Car3", "Exc L5/6 IT Car3"),
  c("L6 CT", "Exc L6 CT"),
  c("Pax6","Inh PAX6"),
  c("Astrocyte", "Ast"),
  c("OPC", "OPC"),
  c("Vip",  "Inh VIP"),
  c("L6b", "Exc L6b"),
  c("L2/3 IT", "Exc L2-3 IT")
)

find_overlapping_mods(paths = paths,
                      pathnames = pathnames,
                      common_pool = common_pool,
                      save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_ROSMAP.qs"))
# $conVsAllAD
# [1] 0.2404533

# $conVsEarly
# [1] 2.55855e-08

# $earlyVsLate
# [1] 1



# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

overlap1 <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_ROSMAP.qs"))



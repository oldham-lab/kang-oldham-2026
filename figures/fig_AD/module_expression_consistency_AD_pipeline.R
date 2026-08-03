library(DESeq2)
library(qs)
library(data.table)
library(AnnotationHub)
library(tidyverse)
library(cowplot)
options(bitmapType = 'cairo')

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/module_expression_consistency_fxns.R"))
 

############
# SEAAD2024 (DFC)
##############

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# conVsAll
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized/conVAll")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("AD", "Con"),
               header_names = c("Gabitto et al. 2024"))

# # conVsAll subset
# save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized/conVAll")
# copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
# objs1 <- find_cons_mods(save_path = save_path, 
#                copa_dir1 = copa_dir1,
#                return_obj = T,
#                plot = T,
#                subset = c())

# Load copa_compare data (conVsEarly)
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized/conVEarly/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("Con", "Early"),
                       header_names = c("Gabitto et al. 2024"))

# early vs late
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized/earlyVLate")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/")
objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        return_obj = T,
                        plot = T,
                        comp_names = c("Late", "Early"),
                        header_names = c("Gabitto et al. 2024"))

# # con vs late
# save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized/conVLate")
# copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsLate/")
# objs2 <- find_cons_mods(save_path = save_path,
#                         copa_dir1 = copa_dir1,
#                         return_obj = T)
 
############
# SEAAD2024 (MTG)
##############

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# conVsAll
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MTG_unnormalized/conVAll")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("AD", "Con"),
               header_names = c("Gabitto et al. 2024"))

# Load copa_compare data (conVsEarly)
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MTG_unnormalized/conVEarly/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_conVEarly/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("Con", "Early"),
                       header_names = c("Gabitto et al. 2024"))

# early vs late
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MTG_unnormalized/earlyVLate")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_earlyVLate/")
objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        return_obj = T,
                        plot = T,
                        comp_names = c("Late", "Early"),
                        header_names = c("Gabitto et al. 2024"))

###############
# Morabito 2021
###############

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Load copa_compare data
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/morabito/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/Morabito_ABIanno/")
objs3 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        return_obj = T)



#############################
# MIT AD Multiome Multiregion (MTC)
#############################

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Con vs all AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_MTC/conVAll/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("AD", "Con"),
                       header_names = c("Liu et al. 2025"))

# Con vs early AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_MTC/conVEarly/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_MTC/")
objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_MTC/"),
                        return_obj = T,
                        plot = T,
                        comp_names = c("AD", "Con"),
                       header_names = c("Liu et al. 2025"))

# Early vs late AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_MTC/earlyVLate/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_MTC/")
objs3 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_MTC/"),
                        return_obj = T,
                        plot = T,
                        comp_names = c("AD", "Con"),
                       header_names = c("Liu et al. 2025"))

# # Con vs late AD
# save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_MTC/conVLate/")
# copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVLate_MTC/")
# objs4 <- find_cons_mods(save_path = save_path,
#                         copa_dir1 = copa_dir1,
#                         sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVLate_MTC/"),
#                         return_obj = T,
#                         plot = F,
#                         comp_names = c("AD", "Con"))

#############################
# MIT AD Multiome Multiregion (PFC)
#############################

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Con vs all AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_PFC/conVAll/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("AD", "Con"),
                       header_names = c("Liu et al. 2025"))

# Con vs early AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_PFC/conVEarly/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_PFC/")
objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_PFC/"),
                        return_obj = T,
                        plot = T,
                        comp_names = c("Con", "Early"),
                       header_names = c("Liu et al. 2025"))

# Early vs late AD
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_PFC/earlyVLate/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_PFC/")
objs3 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_PFC/"),
                        return_obj = T,
                        plot = T,
                        comp_names = c("Late", "Early"),
                       header_names = c("Liu et al. 2025"))

############
# Overlap (SEAAD2024 PFC, MIT_Multiome PFC)
##############

# Find overlap
paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_PFC"))
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
                      save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_PFC.qs"))
# $conVAll
# [1] 1.843936e-19

# $conVEarly
# [1] 0.9897561

# $earlyVLate
# [1] 0.5459825


# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

overlap1 <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_PFC.qs"))


# conVsAll
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_PFC/conVAll")
#copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
#copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/")
copa_dir1 <- list(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/"),
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/"))
over <- tapply(overlap1[[1]]$mod, overlap1[[1]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("AD", "Con"),
               subset = over,
               header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))
 
# Load copa_compare data (conVsEarly)
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_PFC/conVEarly/")
copa_dir1 <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/"),
               file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_PFC/"))
over <- tapply(overlap1[[2]]$mod, overlap1[[2]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("Con", "Early"),
                       subset = over,
                       header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))

# early vs late
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_PFC/earlyVLate")
copa_dir1 <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/"),
               file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_PFC/"))
over <- tapply(overlap1[[3]]$mod, overlap1[[3]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        return_obj = T,
                        plot = T,
                        comp_names = c("Late", "Early"),
                        subset = over,
                        header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))

####################################
# Overlap (SEAAD2024 MTG vs MIT MTG)
####################################
# Find overlap
paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MTG_unnormalized"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_MTC"))
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
                      save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_MTG.qs"))
# Fisher P-values:
# $conVAll
# [1] 9.032377e-24
# $conVEarly
# [1] 0.9981529
# $earlyVLate
# [1] 6.966233e-42

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

overlap1 <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_MTG.qs"))

# conVsAll
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_MTG/conVAll")
copa_dir1 <- list(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG/"),
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/"))
over <- tapply(overlap1[[1]]$mod, overlap1[[1]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("AD", "Con"),
               subset = over,
               header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))
 
# Load copa_compare data (conVsEarly)
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_MTG/conVEarly/")
copa_dir1 <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_conVEarly/"),
               file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_MTC/"))
over <- tapply(overlap1[[2]]$mod, overlap1[[2]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       plot = T,
                       comp_names = c("Con", "Early"),
                       subset = over,
                       header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))

# early vs late
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_MTG/earlyVLate")
copa_dir1 <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_earlyVLate/"),
               file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_MTC/"))
over <- tapply(overlap1[[3]]$mod, overlap1[[3]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        return_obj = T,
                        plot = T,
                        comp_names = c("Late", "Early"),
                        subset = over,
                        header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))




#################################################################
# Plot modules that are shared between SEAAD2024 and MIT Multiome
#################################################################
# (Projected onto SEA)

paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_unnormalized"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion"))
pathnames <- c("bulk_megaset", "bulk_MIT")
 
# List of significant modules
pathList <- lapply(paths, \(d){
  compname <- list.files(d)
  modlist <- lapply(compname, \(x){
    path <- file.path(d, x, "modlist_filtered.qs")
    if(file.exists(path)){
      df <- qread(file.path(d, x, "modlist_filtered.qs"))
      outlist <- lapply(df, \(a) a$mod)
      names(outlist) <- c("pos", "neg")
    } else {
      outlist <- list()
    }
    return(outlist)
  })
  names(modlist) <- compname
  return(modlist)
})
names(pathList) <- pathnames

## SEA
# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# conVsAll
subset_vec <- intersect(pathList[[1]]$conVAll |> unlist(), pathList[[2]]$conVAll |> unlist())
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEA_MIT_overlap/conVAll_on_SEA/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               subset = subset_vec)

# conVsEarly
subset_vec <- intersect(pathList[[1]]$conVEarly |> unlist(), pathList[[2]]$conVEarly |> unlist())
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEA_MIT_overlap/conVEarly_on_SEA/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       subset = subset_vec)

# early vs late
subset_vec <- intersect(pathList[[1]]$earlyVLate |> unlist(), pathList[[2]]$earlyVLate |> unlist())
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEA_MIT_overlap/earlyVLate_on_SEA/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/")
objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        return_obj = T,
                        subset = subset_vec)

## MIT
# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Load copa_compare data
subset_vec <- intersect(pathList[[1]]$conVAll |> unlist(), pathList[[2]]$conVAll |> unlist())
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEA_MIT_overlap/conVAll_on_MIT/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/")
objs <- find_cons_mods(save_path = save_path,
                       copa_dir1 = copa_dir1,
                       return_obj = T,
                       subset = subset_vec)

# Con vs early AD
subset_vec <- intersect(pathList[[1]]$conVEarly |> unlist(), pathList[[2]]$conVEarly |> unlist())
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEA_MIT_overlap/conVEarly_on_MIT/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly/")
objs2 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly/"),
                        return_obj = T,
                        subset = subset_vec)

# Early vs late AD
subset_vec <- intersect(pathList[[1]]$earlyVLate |> unlist(), pathList[[2]]$earlyVLate |> unlist())
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEA_MIT_overlap/earlyVLate_on_MIT/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate/")
objs3 <- find_cons_mods(save_path = save_path,
                        copa_dir1 = copa_dir1,
                        sn_summary_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate/"),
                        return_obj = T,
                        subset = subset_vec)

###############
# APOE 3/3 vs 4/4 
################
# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# SEAAD2024, DFC
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_APOE_DFC/44vs33/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_DFC_apoe_44_vs_33/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("44", "33"),
               header_names = c("Gabitto et al. 2024"))

# MIT DFC
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_APOE_DFC/44vs33/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("44", "33"),
               header_names = c("Liu et al. 2025"))

# SEAAD2024, MTG
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_APOE_MTG/44vs33/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_apoe_44_vs_33/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("44", "33"),
               header_names = c("Gabitto et al. 2024"))

# MIT MTG
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_APOE_MTC/44vs33/")
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE_MTC/")
objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("44", "33"),
               header_names = c("Liu et al. 2025"))

############
# Overlap (SEAAD2024 PFC, MIT_Multiome PFC), APOE 4/4 vs 3/3
##############

# Find overlap
paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_APOE_DFC/"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_APOE_DFC/"))
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
                      save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_PFC_APOE.qs"))
# $`44vs33`
# [1] 1

overlap1 <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_PFC_APOE.qs"))

save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_PFC_APOE/44vs33/")
copa_dir1 <- list(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_DFC_apoe_33_vs_44/"),
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE/"))
over <- tapply(overlap1[[1]]$mod, overlap1[[1]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("44", "33"),
               subset = over,
               header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))

############
# Overlap (SEAAD2024 MTG, MIT_Multiome MTC), APOE 4/4 vs 3/3
##############

# Find overlap
paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_APOE_MTG/"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion_APOE_MTC/"))
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
                      save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_MTG_APOE.qs"))
# $`44vs33`
# [1] 1

overlap1 <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_MTG_APOE.qs"))

save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_MTG_APOE/44vs33/")
copa_dir1 <- list(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_apoe_33_vs_44/"),
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE_MTC/"))
over <- tapply(overlap1[[1]]$mod, overlap1[[1]]$type.x, list) |>
  lapply(unique) # names have to be "pos" and "neg"

objs1 <- find_cons_mods(save_path = save_path, 
               copa_dir1 = copa_dir1,
               return_obj = T,
               plot = T,
               comp_names = c("44", "33"),
               subset = over,
               header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))

############
# Overlap APOE 4/4 vs 3/3, region AND dataset (SEA, MIT; DFC ,MTG)
##############

# Find overlap
paths <- c(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_PFC_APOE"),
           file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_MTG_APOE"))
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
                      save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_DFCandMTG_APOE.qs"))
# $`44vs33`
# [1] 5.558864e-05

overlap1 <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/seaad2024_mit_dcopa_overlap_DFCandMTG_APOE.qs"))

# save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEAAD2024_MITMultiome_overlap_APOE_DFCandMTG/44vs33/")
# copa_dir1 <- list(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_apoe_33_vs_44/"),
#                   file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE_MTC/"))
# over <- tapply(overlap1[[1]]$mod, overlap1[[1]]$type.x, list) |>
#   lapply(unique) # names have to be "pos" and "neg"

# objs1 <- find_cons_mods(save_path = save_path, 
#                copa_dir1 = copa_dir1,
#                return_obj = T,
#                plot = T,
#                comp_names = c("44", "33"),
#                subset = over,
#                header_names = c("Gabitto et al. 2024", "Liu et al. 2025"))
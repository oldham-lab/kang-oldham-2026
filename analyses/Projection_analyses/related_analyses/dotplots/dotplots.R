library(tidyverse)
library(data.table)
library(qs)
library(ggpubr)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/calculate_rand_euclidean_distances.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/dotplots/fxns.R"))

# Use new input 

######################################
# Generate summary dotplot for sig mods that overlap between SEAAD2024 and MIT (bulk megaset modules)
######################################
# - Add class info
class_info <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |>
  arrange(Subclass) |>
  mutate(Subclass_fixed = factor(c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
                                   "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip"),
                                 levels = c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
                                             "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
                                             "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC"))) |>
  arrange(Subclass_fixed) |>
  select(Subclass, Class)

# Load module data (topmodposbc - all dCoPA analyses done using topmodposbc definitions)
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
filter_under <- 3
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])


## con vs All (DFC)
# Load output table for shared modules
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon/shared_output_table.csv")) |>
  dplyr::filter(mod %in% these_mods)
# length(unique(paste0(sum_tab$mod, sum_tab$Celltype1)))
# # 190
# length(unique(paste0(sum_tab$mod, sum_tab$Celltype_Liu_2025)))
# # 206
# # There are fewer unique mod-subclass combinations using Gabitto annotations compared to Liu annotations
# Subset to unique mod-subclass combinations (using Gabitto annotations)
# (necessary due to the many-to-many nature of subclass label mapping between Gabitto and Liu)
# (need to add this to pipeline code)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
                           
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("all AD", "con"),
                                 file_suffix = "conVsAllAD",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); control vs all AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_DFC/"))

## early vs con (DFC)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_earlyVsCon/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("early AD", "con"),
                                 file_suffix = "earlyADVscon",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); early AD vs con samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_DFC/"))

## late vs early (DFC)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_lateVsEarly/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("late AD", "early AD"),
                                 file_suffix = "lateVsEarlyAD",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); late vs early AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_DFC/"))

## con vs all (MTG)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("all AD", "con"),
                                 file_suffix = "conVsAllAD",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 MTG (SEAAD2024) and Liu et al. 2025 MTG (MIT_Multiome_Multiregion); con vs all AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_MTG/"))

## early vs con (MTG)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_earlyVsCon/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("early AD", "con"),
                                 file_suffix = "earlyADVscon",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 MTG (SEAAD2024) and Liu et al. 2025 MTG (MIT_Multiome_Multiregion); early AD vs con samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_MTG/"))


## late vs early (MTG)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_lateVsEarly/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("late AD", "early AD"),
                                 file_suffix = "lateVsEarlyAD",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 MTG (SEAAD2024) and Liu et al. 2025 MTG (MIT_Multiome_Multiregion); late vs early AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_MTG/"))

## APOE 3/3 vs 4/4
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_APOE_DFC/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("44", "33"),
                                 file_suffix = "APOE44vs33",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 MTG (SEAAD2024) and Liu et al. 2025 MTG (MIT_Multiome_Multiregion); APOE 4/4 vs 3/3 samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_DFC_APOE/"))

## con vs all (DFC, ROSMAP)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon_ROSMAP/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("all AD", "con"),
                                 file_suffix = "conVsAllAD_DFC",
                                 plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1010 total modules, FDR cutoff",
                                 plot_caption = "Mods (ROSMAP) are significant in both Gabitto et al. 2024 DFC (SEAAD2024) and Liu et al. 2025 DFC (MIT_Multiome_Multiregion); con vs all AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_ROSMAPmods/"))

## early vs con (DFC, ROSMAP)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_earlyVsCon_ROSMAP/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("early AD", "con"),
                                 file_suffix = "earlyADVscon_DFC",
                                 plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1010 total modules, FDR cutoff",
                                 plot_caption = "Mods (ROSMAP) are significant in both Gabitto et al. 2024 DFC (SEAAD2024) and Liu et al. 2025 DFC (MIT_Multiome_Multiregion); early AD vs con samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_ROSMAPmods/"))

## late vs early (DFC, ROSMAP)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_lateVsEarly_ROSMAP/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("late AD", "early AD"),
                                 file_suffix = "lateVsEarlyAD_DFC",
                                 plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1010 total modules, FDR cutoff",
                                 plot_caption = "Mods (ROSMAP) are significant in both Gabitto et al. 2024 DFC (SEAAD2024) and Liu et al. 2025 DFC (MIT_Multiome_Multiregion); late AD vs early AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_ROSMAPmods/"))

## con vs all (MTG, ROSMAP)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon_ROSMAP/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("all AD", "con"),
                                 file_suffix = "conVsAllAD_MTG",
                                 plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1010 total modules, FDR cutoff",
                                 plot_caption = "Mods (ROSMAP) are significant in both Gabitto et al. 2024 MTG (SEAAD2024) and Liu et al. 2025 MTG (MIT_Multiome_Multiregion); con vs all AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_ROSMAPmods/"))

## early vs con (DFC, ROSMAP)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_earlyVsCon_ROSMAP/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("early AD", "con"),
                                 file_suffix = "earlyADVscon_MTG",
                                 plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1010 total modules, FDR cutoff",
                                 plot_caption = "Mods (ROSMAP) are significant in both Gabitto et al. 2024 MTG (SEAAD2024) and Liu et al. 2025 MTG (MIT_Multiome_Multiregion); early AD vs con samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_ROSMAPmods/"))

## late vs early (DFC, ROSMAP)
sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_lateVsEarly_ROSMAP/shared_output_table.csv"))|>
  dplyr::filter(mod %in% these_mods)
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("late AD", "early AD"),
                                 file_suffix = "lateVsEarlyAD_MTG",
                                 plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                                 plot_subtitle = "1010 total modules, FDR cutoff",
                                 plot_caption = "Mods (ROSMAP) are significant in both Gabitto et al. 2024 MTG (SEAAD2024) and Liu et al. 2025 MTG (MIT_Multiome_Multiregion); late AD vs early AD samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_ROSMAPmods/"))


######################################
# Generate summary dotplot for sig mods that overlap between brainSCOPE CMC and SZBD
######################################

# brainSCOPE

class_info <- data.frame("Subclass" = c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
                                        "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
                                        "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro", "Endo", "Immune", "PC", "SMC", "VLMC"),
                               "Class" = c(rep("Glutamatergic", 9), rep("GABAergic", 9), rep("Non-neuronal", 9))) |>
  mutate(Subclass = factor(Subclass, levels = unique(Subclass)))

sum_tab <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/brainSCOPE_CMC_vs_SZBD/shared_output_table.csv"))
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype1)), ] 
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("SCZ", "con"),
                                 file_suffix = "SCZvsCon",
                                 plot_title = "# of significant dCoPA mods shared between brainSCOPE cohorts (CMC vs SZBDMultiseq)",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in two cohorts (CMC, SZBDMultiseq) from Emani et al. 2024 (brainSCOPE); control vs SCZ samples",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/brainSCOPE_CMC_SZBD_shared/"))

#############################
## Create a dotplot for overlap between DFC and MTG, Gabitto and Liu, all AD vs con (4 datasets total)
#############################

# First, find overlap between 4 datasets
map_list <- list(
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

dat1 = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/output_table_Subclass.csv"))
dat2 = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/euclidean_distances/output_table_Subclass.csv"))
dat3 = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG/euclidean_distances/output_table_Subclass.csv"))
dat4 = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/euclidean_distances/output_table_Subclass.csv"))
              
sig_type = "fdr"
dataset_names = c("Gabitto_2024_DFC", "Liu_2025_DFC", "Gabitto_2024_MTC", "Liu_2025_MTC")
save_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFCandMTG_conVAll/")

dat1 <- dat1 |> dplyr::filter(sig_FDR)
dat2 <- dat2 |> dplyr::filter(sig_FDR)
dat3 <- dat3 |> dplyr::filter(sig_FDR)
dat4 <- dat4 |> dplyr::filter(sig_FDR)
dat1 <- dat1 |> dplyr::filter(Consistency %in% c(0, 1))
dat2 <- dat2 |> dplyr::filter(Consistency %in% c(0, 1))
dat3 <- dat3 |> dplyr::filter(Consistency %in% c(0, 1))
dat4 <- dat4 |> dplyr::filter(Consistency %in% c(0, 1))

#dat_out <- dplyr::inner_join(dat1, dat2, by = dplyr::join_by(mod, Direction, Consistency), suffix = dataset_names, relationship = "many-to-many")
dat_out <- Reduce(\(x,y){
  dplyr::inner_join(x, y, by = dplyr::join_by(mod, Direction, Consistency), relationship = "many-to-many")
}, list(dat1, dat2, dat3, dat4))

col_index <- grep("Celltype", colnames(dat_out))
keep_these <- lapply(map_list, \(x){
    which(dat_out[ ,col_index[1]] %in% x & dat_out[ ,col_index[2]] %in% x & dat_out[ ,col_index[3]] %in% x & dat_out[ ,col_index[4]] %in% x)
}) |> unlist() |> unique()
dat_out1 <- dat_out[keep_these, ] |>
    filter(!duplicated(paste(mod, Celltype.x))) |>
    dplyr::arrange(mod) |>
    select(1:9)
fwrite(dat_out1, file = file.path(save_dir, "shared_output_table.csv"))


# - Add class info
class_info <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |>
  arrange(Subclass) |>
  mutate(Subclass_fixed = factor(c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
                                   "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip"),
                                 levels = c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
                                             "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
                                             "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC"))) |>
  arrange(Subclass_fixed) |>
  select(Subclass, Class)

sum_tab <- fread(data.table = F, file = file.path(save_dir, "shared_output_table.csv"))
sum_tab <- sum_tab[!duplicated(paste0(sum_tab$mod, sum_tab$Celltype.x)), ] |>
  rename("Celltype" = "Celltype.x")                    
create_shared_dotplot_summary_v2(sum_tab = sum_tab,
                                 class_info = class_info,
                                 splits = c("all AD", "con"),
                                 file_suffix = "AllADVsCon",
                                 plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025, DFC and MTG, control vs all AD samples ",
                                 plot_subtitle = "1023 total modules, FDR cutoff",
                                 plot_caption = "Mods are significant in both Gabitto et al. 2024 (SEAAD2024) and Liu et al. 2025 (MIT_Multiome_Multiregion); DFC and MTG",
                                 save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_DFCandMTG/"))

#######################################
# For all comparisons, create dotplots and put in single plot
#######################################

# First, load all data:
dir1 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/")
dot_list <- list(
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_DFC.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_earlyADVscon.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_earlyADVscon_DFC.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_lateVsEarlyAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_lateVsEarlyAD_DFC.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_DFC_APOE", "dcopa_scorecard_summary_APOE44vs33.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_MTG.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_earlyADVscon.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_earlyADVscon_MTG.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_lateVsEarlyAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_lateVsEarlyAD_MTG.csv"),
  file.path(dir1, "brainSCOPE_CMC_SZBD_shared", "dcopa_scorecard_summary_SCZvsCon.csv")
) |>
  lapply(\(x) fread(x, data.table = F))
typevec <- c("All AD vs Con (DFC), p = 5e-157",
             "All AD vs Con (ROSMAP, DFC), p = 2e-241",
             "Early AD vs Con (DFC), p = 2e-20",
             "Early AD vs Con (ROSMAP, DFC), p = 4e-24",
             "Late vs Early AD (DFC), p = 1e-11",
             "Late vs Early AD (ROSMAP, DFC), p = 1e-24",
             "APOE 4/4 vs 3/3 (DFC), p = 4e-2",
             "All AD vs Con (MTG), p = 1e-251",
             "All AD vs Con (ROSMAP, MTG), p = 0e-000",
             "Early AD vs Con (MTG), p = 5e-29",
             "Early AD vs Con (ROSMAP, MTG), p = 9e-51",
             "Late vs Early AD (MTG), p = 2e-243",   
             "Late vs Early AD (ROSMAP, MTG), p = 1e-255",
             "SCZ vs Con (DFC), p = 9e-28")

dot_list <- mapply(\(l, name){
  l$comp <- name
  return(l)
}, dot_list, typevec, SIMPLIFY = F)

cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

# p <- dot_list |>
#   do.call(what = "rbind") |>
#   ggplot(aes(x = Celltype, y = type, size = num_sig, fill = Class)) +
#     theme_minimal() + 
#     geom_point(color = "black", pch = 21) +
#     theme(text = element_text(size = 6),
#           axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.2, 0, 0, 0, "cm")),
#           axis.text.y = element_text(size = 6), 
#           strip.text = element_text(size = 7),
#           legend.position = "bottom",
#           legend.spacing.y = unit(0, "mm"),
#           legend.title = element_blank(),
#           legend.margin = margin(-0.5, 0, 0, 0, "cm"),
#           legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()) +
#     labs(x = "", y = "") +#,
#         #title = plot_title,
#         #subtitle = plot_subtitle,
#         #caption =  plot_caption) +
#     scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
#     scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) +
#     facet_wrap(~comp, ncol = 2, scales = "free_y") 

# ggsave(p, file = file.path(Sys.getenv("SCRATCH_DIR", "~/test"), "test1.pdf"), width = 7, height = 6)

p2 <- mapply(\(x, y){
  x |>
    mutate(type = factor(type, levels = rev(unique(type)))) |>
    ggplot(aes(x = Celltype, y = type, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.2, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 6), 
            strip.text = element_text(size = 7),
            legend.position = "bottom",
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 8, hjust = 0.5)) +
      labs(x = "", y = "",
          title = y) +#,
          #subtitle = plot_subtitle,
          #caption =  plot_caption) +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3))  +
      scale_y_discrete(drop = FALSE)
}, dot_list, typevec, SIMPLIFY = F)

pout <- cowplot::plot_grid(plotlist = p2, ncol = 2, byrow = F) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave(pout, file = file.path(dir1, "summary_plots", "combined_summary_dotplots.pdf"), height = 12, width = 10)
ggsave(pout, file = file.path(dir1, "summary_plots", "combined_summary_dotplots.png"), height = 12, width = 10, bg = "white")

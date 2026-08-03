library(tidyverse)
library(data.table)
library(qs)
library(ggpubr)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/calculate_rand_euclidean_distances.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/dotplots/fxns.R"))

###########
# SEAAD2024
###########

sn_anno_subclass <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |> 
  arrange(Subclass) |>
  mutate(Subclass_final = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) 

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")


# conVAll
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 (SEAAD2024) control vs all AD samples",
                  splits = c("all AD", "con"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )
# conVEarly
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_conVsEarly/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 (SEAAD2024) control vs early AD samples",
                  splits = c("con", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = T
                  )

# earlyVLate
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_earlyVsLate/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 (SEAAD2024) early vs late AD samples",
                  splits = c("late AD", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )

###########
# SEAAD2024 (MTC)
###########

sn_anno_subclass <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |> 
  arrange(Subclass) |>
  mutate(Subclass_final = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) 

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")


# conVAll
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MTG/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 MTG (SEAAD2024) control vs all AD samples",
                  splits = c("all AD", "con"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )

# conVEarly
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_conVEarly/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MTG_conVEarly/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 MTG (SEAAD2024) control vs early AD samples",
                  splits = c("con", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = T
                  )

# earlyVLate
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_earlyVLate/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MTG_earlyVLate/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 MTG (SEAAD2024) early vs late AD samples",
                  splits = c("late AD", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )


#################
# SEAAD2024 DFC (ROSMAP modules)
#################
ct_order1 <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")

sn_anno_subclass <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |> 
  arrange(Subclass) |>
  mutate(Subclass_final = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) |>
  arrange(factor(Subclass_final, levels = ct_order1))

ct_order <- sn_anno_subclass$Subclass

# conVAll
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_DFC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_DFC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_DFC_ROSMAPmods/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 (SEAAD2024) control vs all AD samples",
                  splits = c("all AD", "con"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  mod_count = 1010,
                  rev_bool = F
                  )

# conVEarly
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_conVsEarly/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_DFC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_DFC_ROSMAPmods_conVEarly/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 (SEAAD2024) control vs early AD samples",
                  splits = c("con", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  mod_count = 1010,
                  rev_bool = T
                  )

# earlyVsLate
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_DFC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_DFC_ROSMAPmods_earlyVsLate/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 (SEAAD2024) control vs early AD samples",
                  splits = c("late AD", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  mod_count = 1010,
                  rev_bool = F
                  )

#################
# SEAAD2024 (APOE)
#################
ct_order1 <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")

sn_anno_subclass <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |> 
  arrange(Subclass) |>
  mutate(Subclass_final = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) |>
  arrange(factor(Subclass_final, levels = ct_order1))

ct_order <- sn_anno_subclass$Subclass

# DFC
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_DFC_apoe_33_vs_44"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_DFC_APOE/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 DFC (SEAAD2024) APOE 4/4 vs 3/3 samples",
                  splits = c("4/4", "3/3"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  mod_count = 1023,
                  rev_bool = F
                  )

# MTG
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_MTG_apoe_33_vs_44"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MTG_APOE/"),
                  caption1 = "Dataset used: Gabitto et al. 2024 MTG (SEAAD2024) APOE 4/4 vs 3/3 samples",
                  splits = c("4/4", "3/3"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  mod_count = 1023,
                  rev_bool = F
                  )

##########################
# MIT_Multiome_Multiregion (MTC)
##########################

obj <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/sn_summary_tables/sn_summary_objects_log.qs"))

sn_anno_subclass <- fread(data.table = F, file = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv") |> 
  select(RNA.Class, RNA.Subclass) |>
  filter(!duplicated(RNA.Subclass), 
          RNA.Subclass %in% colnames(obj[[1]][[1]][[1]])) |>
  mutate(Class = case_match(RNA.Class, "Exc" ~ "Glutamatergic", "Inh" ~ "GABAergic", .default = "Non-neuronal")) |>
  rename(Subclass = RNA.Subclass) |>
  arrange(factor(Class, levels = c("Glutamatergic", "GABAergic", "Non-neuronal"))) |>
  mutate(Subclass_final = Subclass)

ct_order <- unique(sn_anno_subclass$Subclass)

# conVAll
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITMTC/"),
                  caption1 = "Dataset used: Liu et al. 2025 MTC (MIT_Multiome_Multiregion) control vs all AD samples",
                  splits = c("all AD", "con"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )
# # conVEarly
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_MTC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITMTC_conVEarly/"),
                  caption1 = "Dataset used: Liu et al. 2025 MTC (MIT_Multiome_Multiregion) control vs early AD samples",
                  splits = c("con", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = T
                  )

# # earlyVLate
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_MTC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_MTC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITMTC_earlyVLate/"),
                  caption1 = "Dataset used: Liu et al. 2025 MTC (MIT_Multiome_Multiregion) early vs late AD samples",
                  splits = c("late AD", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )

##########################
# MIT_Multiome_Multiregion (PFC)
##########################

obj <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/sn_summary_tables/sn_summary_objects_log.qs"))

sn_anno_subclass <- fread(data.table = F, file = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv") |> 
  select(RNA.Class, RNA.Subclass) |>
  filter(!duplicated(RNA.Subclass), 
          RNA.Subclass %in% colnames(obj[[1]][[1]][[1]])) |>
  mutate(Class = case_match(RNA.Class, "Exc" ~ "Glutamatergic", "Inh" ~ "GABAergic", .default = "Non-neuronal")) |>
  rename(Subclass = RNA.Subclass) |>
  arrange(factor(Class, levels = c("Glutamatergic", "GABAergic", "Non-neuronal"))) |>
  mutate(Subclass_final = Subclass)

ct_order <- unique(sn_anno_subclass$Subclass)

# conVAll
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITPFC/"),
                  caption1 = "Dataset used: Liu et al. 2025 PFC (MIT_Multiome_Multiregion) control vs all AD samples",
                  splits = c("all AD", "con"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )
# # conVEarly
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_conVEarly_PFC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITPFC_conVEarly/"),
                  caption1 = "Dataset used: Liu et al. 2025 PFC (MIT_Multiome_Multiregion) control vs early AD samples",
                  splits = c("con", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = T
                  )

# # earlyVLate
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVLate_PFC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MITPFC_earlyVLate/"),
                  caption1 = "Dataset used: Liu et al. 2025 PFC (MIT_Multiome_Multiregion) early vs late AD samples",
                  splits = c("late AD", "early AD"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )

#################
# MIT (APOE)
#################
obj <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_PFC/sn_summary_tables/sn_summary_objects_log.qs"))

sn_anno_subclass <- fread(data.table = F, file = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv") |> 
  select(RNA.Class, RNA.Subclass) |>
  filter(!duplicated(RNA.Subclass), 
          RNA.Subclass %in% colnames(obj[[1]][[1]][[1]])) |>
  mutate(Class = case_match(RNA.Class, "Exc" ~ "Glutamatergic", "Inh" ~ "GABAergic", .default = "Non-neuronal")) |>
  rename(Subclass = RNA.Subclass) |>
  arrange(factor(Class, levels = c("Glutamatergic", "GABAergic", "Non-neuronal"))) |>
  mutate(Subclass_final = Subclass)

ct_order <- unique(sn_anno_subclass$Subclass)

# DFC
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MIT_DFC_APOE/"),
                  caption1 = "Dataset used: Liu et al. 2025 PFC (MIT_Multiome_Multiregion) APOE 4/4 vs 3/3",
                  splits = c("4/4", "3/3"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  mod_count = 1023,
                  rev_bool = F
                  )

# MTG
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE_MTC"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/MIT_MTG_APOE/"),
                  caption1 = "Dataset used: Liu et al. 2025 MTG (MIT_Multiome_Multiregion) APOE 4/4 vs 3/3",
                  splits = c("4/4", "3/3"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  mod_count = 1023,
                  rev_bool = F
                  )


############
# brainSCOPE
############

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro", "Endo", "Immune", "PC", "SMC", "VLMC")

sn_anno_subclass <- data.frame("Subclass" = ct_order,
                                "Class" = c(rep("Glutamatergic", 9), rep("GABAergic", 9), rep("Non-neuronal", 9))) |>
  mutate(Subclass_final = Subclass)

# conVSCZ (CMC)
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE_CMC/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/brainSCOPE_CMC/"),
                  caption1 = "Dataset used: Emani et al. 2024 (brainSCOPE CMC); control vs SCZ samples",
                  splits = c("SCZ", "con"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )

# conVSCZ (SZBDMultiseq)
rand_mod_analysis(sn_summary_object_path = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/brainSCOPE_SZBDMultiseq/"),
                  module_output_dir = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024"),
                  save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/brainSCOPE_SZBDMultiseq/"),
                  caption1 = "Dataset used: Emani et al. 2024 (brainSCOPE SZBDMultiseq); control vs SCZ samples",
                  splits = c("SCZ", "con"),
                  calc_rand_mods = F,
                  plot_real_vs_rand = F,
                  sn_anno_subclass = sn_anno_subclass,
                  ct_order = ct_order,
                  rev_bool = F
                  )

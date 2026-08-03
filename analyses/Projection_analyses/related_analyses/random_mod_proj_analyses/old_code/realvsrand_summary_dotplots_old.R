library(tidyverse)
library(data.table)
library(qs)
library(ggpubr)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/calculate_rand_euclidean_distances.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/random_mod_proj_analyses/fxns.R"))

######################################
# Generate summary dotplot for sig mods that overlap between SEAAD2024 and MIT (bulk megaset modules)
######################################
save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared/")
if(!dir.exists(save_dir1)){dir.create(save_dir1)}

# Load sig modules (conVAll)
overlap1 <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/seaad2024_mit_dcopa_overlap_PFC.qs"))
overlap1 <- lapply(overlap1, \(x){
  out <- x |> 
    rename("type" = "type.x", "subclass" = "subclass.x")
}) # use Gabitto subclass names

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
  mutate(Subclass = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) 

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")

create_shared_dotplot_summary(slist = overlap1[[1]],
                              splits = c("all AD", "con"),
                              file_suffix = "conVsAllAD",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                              plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); control vs all AD samples",
                              rev_bool =  F,
                              save_dir1 = save_dir1)

create_shared_dotplot_summary(slist = overlap1[[2]],
                              splits = c("con", "early AD"),
                              file_suffix = "conVsEarly",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                              plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); control vs early AD samples",
                              rev_bool =  T,
                              save_dir1 = save_dir1)

create_shared_dotplot_summary(slist = overlap1[[3]],
                              splits = c("late AD", "early AD"),
                              file_suffix = "earlyVsLate",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                              plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); early vs late AD samples",
                              rev_bool =  F,
                              save_dir1 = save_dir1)

# How many modules are shared in common across shared AD and shared SCZ modules?
overlap_scz <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/brainscope_cmc_szbd_dcopa_overlap.qs"))

length(unique(overlap1[[1]]$mod))
# 89 unique mods for AD conVAll (86 lower in AD, 3 higher in AD)
length(unique(overlap_scz[[1]]$mod))
# 21 unique mods for SCZ (18 lower in SCZ, 3 higher in SCZ)
sum(unique(overlap_scz[[1]]$mod) %in% unique(overlap1[[1]]$mod))
# 13 mods shared between the two (all lower in AD)

all_subclasses <- unique(c(overlap1[[1]]$subclass, overlap_scz[[1]]$subclass))
overlap_join <- full_join(overlap1[[1]][,1:3], overlap_scz[[1]][,1:3], by = join_by(mod))
overlap_join1 <- overlap_join |> filter(!is.na(subclass) & !is.na(subclass.x)) |>
  filter(subclass == subclass.x)
# 7 unique modules that move in same direction and same subclass:
# > overlap_join1
#    mod subclass type subclass.x type.x
# 1   34    L4 IT  neg      L4 IT    neg
# 2  100    L4 IT  neg      L4 IT    neg
# 3  100    L4 IT  neg      L4 IT    neg
# 4  132    L4 IT  neg      L4 IT    neg
# 5  132    L4 IT  neg      L4 IT    neg
# 6  324    L4 IT  neg      L4 IT    neg
# 7  324    L4 IT  neg      L4 IT    neg
# 8  329    L4 IT  neg      L4 IT    neg
# 9  329    L4 IT  neg      L4 IT    neg
# 10 405    L4 IT  neg      L4 IT    neg
# 11 489    L5 IT  neg      L5 IT    neg
# 12 489    L5 IT  neg      L5 IT    neg
# 13 489  L5/6 NP  neg    L5/6 NP    neg

######################################
# Generate summary dotplot for sig mods that overlap between SEAAD2024 and MIT (APOE, DFC)
######################################
save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_DFC_APOE/")
if(!dir.exists(save_dir1)){dir.create(save_dir1)}

# Load sig modules (conVAll)
overlap1 <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/seaad2024_mit_dcopa_overlap_PFC_APOE.qs"))
overlap1 <- lapply(overlap1, \(x){
  out <- x |> 
    rename("type" = "type.x", "subclass" = "subclass.x")
}) # use Gabitto subclass names

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
  mutate(Subclass = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) 

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")

create_shared_dotplot_summary(slist = overlap1[[1]],
                              splits = c("4/4", "3/3"),
                              file_suffix = "APOE",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods shared between Gabitto et al. 2024 and Liu et al. 2025",
                              plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); APOE 4/4 vs 3/3",
                              rev_bool =  F,
                              mod_count = 1023,
                              save_dir1 = save_dir1)


######################################
# Generate summary dotplot for sig mods that overlap between SEAAD2024 and MIT (ROSMAP modules, DFC)
######################################
 
save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/SEAAD2024_MIT_shared_ROSMAPmods/")
if(!dir.exists(save_dir1)){dir.create(save_dir1)}

# Load sig modules (conVAll)
overlap1 <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/seaad2024_mit_dcopa_overlap_ROSMAP.qs"))
overlap1 <- lapply(overlap1, \(x){
  out <- x |> 
    rename("type" = "type.x", "subclass" = "subclass.x")
}) # use Gabitto subclass names

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
  mutate(Subclass = c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip")) 

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC")

create_shared_dotplot_summary(slist = overlap1[[1]],
                              splits = c("all AD", "con"),
                              file_suffix = "conVsAllAD",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                              plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); control vs all AD samples",
                              rev_bool =  F,
                              save_dir1 = save_dir1)

create_shared_dotplot_summary(slist = overlap1[[2]],
                              splits = c("con", "early AD"),
                              file_suffix = "conVsEarly",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                              plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); control vs early AD samples",
                              rev_bool =  T,
                              save_dir1 = save_dir1)

create_shared_dotplot_summary(slist = overlap1[[3]],
                              splits = c("late AD", "early AD"),
                              file_suffix = "earlyVsLate",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods (ROSMAP) shared between Gabitto et al. 2024 and Liu et al. 2025",
                              plot_caption = "Mods are significant in both Gabitto et al. 2024 PFC (SEAAD2024) and Liu et al. 2025 PFC (MIT_Multiome_Multiregion); early vs late AD samples",
                              rev_bool =  F,
                              save_dir1 = save_dir1)

#######################
# Generate summary dotplot for sig mods that overlap between brainSCOPE cohorts (CMC vs SZBDMultiseq)
#######################

save_dir1 = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/brainSCOPE_CMC_SZBD_shared/")
if(!dir.exists(save_dir1)){dir.create(save_dir1)}

# Load sig modules
overlap1 <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/brainscope_cmc_szbd_dcopa_overlap.qs"))
overlap1 <- lapply(overlap1, \(x){
  out <- x |> 
    rename("type" = "type.x", "subclass" = "subclass.x")
}) 

# overlaptest <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/seaad2024_mit_dcopa_overlap_PFC.qs")) # load AD shared mods
# l1 <- unique(overlap1[[1]][,1]) # unique brainscope shared mods
# llist <- lapply(overlaptest, \(x) unique(x[,1])) unique AD shared mods 
# > length(l1)
# [1] 54
# > lapply(llist, length) |> unlist()
#    conVAll  conVEarly earlyVLate 
#         89         17          8 
# > lapply(llist, \(x) sum(x %in% l1)) |> unlist() # of AD shared mods in brainscope shared mods
#    conVAll  conVEarly earlyVLate 
#         24          9          2 
# lapply(llist, \(x) sum(x %in% l1)/length(x)) |> unlist() # % of AD shared mods in brainscope shared mods
#    conVAll  conVEarly earlyVLate 
#  0.2696629  0.5294118  0.2500000 
#  lapply(llist, \(x) sum(l1 %in% x)/length(l1)) |> unlist() # % of brainscope shared mods in AD shared mods
#    conVAll  conVEarly earlyVLate 
# 0.44444444 0.16666667 0.03703704 

# create_shared_dotplot_summary_AD(slist = qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/brainscope_CMC_SZBD_overlap/conVsSCZ/modlist_filtered.qs")),
#                                  comparison = "control vs SCZ samples",
#                                  splits = c("SCZ", "con"),
#                                  file_suffix = "conVsSCZ",
#                                  save_dir1)

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro", "Endo", "Immune", "PC", "SMC", "VLMC")

sn_anno_subclass <- data.frame("Subclass" = ct_order,
                               "Class" = c(rep("Glutamatergic", 9), rep("GABAergic", 9), rep("Non-neuronal", 9)))

create_shared_dotplot_summary(#slist = qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_AD/panels/brainscope_CMC_SZBD_overlap/conVsSCZ/modlist_filtered.qs")),
                              slist = overlap1[[1]],
                              splits = c("SCZ", "con"),
                              file_suffix = "conVsSCZ",
                              sn_anno_subclass = sn_anno_subclass,
                              ct_order = ct_order,
                              plot_title = "# of significant dCoPA mods shared between brainSCOPE cohorts (CMC vs SZBDMultiseq)",
                              plot_caption =  "Mods are significant in two cohorts (CMC, SZBDMultiseq) from Emani et al. 2024 (brainSCOPE); control vs SCZ samples",
                              save_dir1 = save_dir1)

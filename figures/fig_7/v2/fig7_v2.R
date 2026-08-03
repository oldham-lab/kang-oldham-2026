library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v2"
if(!dir.exists(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)))
  dir.create(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version), recursive = T)

##########
# Panel A-B
# Dotplots (Control vs all AD)
##########
# Find overlap between all AllADvsCon comparisons
# v2: find overlap between shared regions (i.e. DFC bulk overlap with DFC ADmods)

# For celltype ordering:
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

# Load dCoPA output files (AllADvsCon DFC, MTG)
file1 <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon"),        
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon_ROSMAP"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon"),        
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon_ROSMAP"))
flist <- lapply(file1, \(x) fread(data.table = F, file = file.path(x, "shared_output_table.csv")) |>
    filter(mod %in% these_mods)
)

# Combined (DFC)
flistjoindfc <- Reduce(\(x, y){
  a1 <- x |> select(c("mod", "Celltype1", "Direction", "Consistency")) 
  a2 <- y |> select(c("mod", "Celltype1", "Direction", "Consistency"))
  out <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype1, Direction, Consistency), relationship = "many-to-many") |>
    filter(mod %in% these_mods)
  return(out)
}, flist[c(1,2)]) |>
  filter(!duplicated(paste0(mod, Celltype1, Direction, Consistency))) |>
    mutate(Direction = factor(Direction, levels = c(1, -1))) |> 
    mutate(Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass))) |>
    select(mod, Celltype1, Direction) |>
    group_by(Celltype1, Direction, .drop = FALSE) |>
    summarise(num_sig = n(), .groups = "drop") |>
    left_join(class_info, by = join_by(Celltype1 == Subclass)) |>
    mutate(Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass))) |>
    arrange(Celltype1) |>
    mutate(num_sig = case_match(num_sig, 
                                0 ~ NA,
                                .default = num_sig
                                ),
           comp = "Combined (DFC)") 

# Combined (ROSMAP modules)
flistjoinmtg <- Reduce(\(x, y){
  a1 <- x |> select(c("mod", "Celltype1", "Direction", "Consistency")) 
  a2 <- y |> select(c("mod", "Celltype1", "Direction", "Consistency"))
  out <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype1, Direction, Consistency), relationship = "many-to-many") |>
    filter(mod %in% these_mods)
  return(out)
}, flist[c(3, 4)]) |>
  filter(!duplicated(paste0(mod, Celltype1, Direction, Consistency))) |>
    mutate(Direction = factor(Direction, levels = c(1, -1))) |> 
    mutate(Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass))) |>
    select(mod, Celltype1, Direction) |>
    group_by(Celltype1, Direction, .drop = FALSE) |>
    summarise(num_sig = n(), .groups = "drop") |>
    left_join(class_info, by = join_by(Celltype1 == Subclass)) |>
    mutate(Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass))) |>
    arrange(Celltype1) |>
    mutate(num_sig = case_match(num_sig, 
                                0 ~ NA,
                                .default = num_sig
                                ),
           comp = "Combined (MTG)") 
flistall <- rbind(flistjoinmtg, flistjoindfc)

## Plot dotplot
dir1 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/")
dot_list <- list(
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_MTG.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_DFC.csv")
) |>
  lapply(\(x){
    fread(x, data.table = F) |> 
      mutate(Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass))) 
  })
names(dot_list) <- c("MTG", "MTG (AD modules)", "DFC", "DFC (AD modules)")

typevec <- c("MTG",
             "MTG (AD modules)",
             "Combined (MTG)",
             "DFC",
             "DFC (AD modules)",
             "Combined (DFC)"
             )

df_hi <- lapply(seq_along(dot_list), \(x){
  dot_list[[x]] |>
    filter(grepl("Higher", type)) |>
    mutate(comp = names(dot_list)[x])
}) |> do.call(what = "rbind") |>
  select(!type)
df_hi <- rbind(df_hi, flistall |> filter(Direction == 1) |> select(!Direction)) |> 
  mutate(comp = factor(comp, levels = rev(typevec)))

df_lo <- lapply(seq_along(dot_list), \(x){
  dot_list[[x]] |>
    filter(grepl("Lower", type)) |>
    mutate(comp = names(dot_list)[x])
}) |> do.call(what = "rbind") |>
  select(!type)
df_lo <- rbind(df_lo, flistall |> filter(Direction == -1) |> select(!Direction)) |> 
  mutate(comp = factor(comp, levels = rev(typevec)))

cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

plo <- df_lo |>
    ggplot(aes(x = Celltype1, y = comp, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 7), 
            legend.direction = "horizontal",    
            legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = unit(4, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(x = "", y = "") +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_A.pdf"), height = 2.6, width = 4.6)
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_A.svg"), height = 2.6, width = 4.6)

phi <- df_hi |>
    ggplot(aes(x = Celltype1, y = comp, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 7), 
            legend.direction = "horizontal",    
            legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = unit(4, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
    #  guides(fill = guide_legend(ncol = 1),
     #        size = guide_legend(ncol = 1)) +
      labs(x = "", y = "") +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(limits = c(1, max(df_lo$num_sig, na.rm = T)), # Manually set the scale to be the same as lo object
                            breaks = scales::pretty_breaks(n = 2)) 
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_B.pdf"), height = 2.6, width = 4.6)
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_B.svg"), height = 2.6, width = 4.6)

#########
# Panel A-B related
# Creating dotplots for shared genes
#########
flist1 <- flist
flist1[[5]] <- Reduce(\(x, y){
  a1 <- x |> select(c("mod", "Celltype1", "Direction", "Consistency")) 
  a2 <- y |> select(c("mod", "Celltype1", "Direction", "Consistency"))
  out <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype1, Direction, Consistency), relationship = "many-to-many") |>
    filter(mod %in% these_mods)
  return(out)
}, flist[c(1,2)]) 

flist1[[6]] <- Reduce(\(x, y){
  a1 <- x |> select(c("mod", "Celltype1", "Direction", "Consistency")) 
  a2 <- y |> select(c("mod", "Celltype1", "Direction", "Consistency"))
  out <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype1, Direction, Consistency), relationship = "many-to-many") |>
    filter(mod %in% these_mods)
  return(out)
}, flist[c(1,2)]) 

names(flist1) <- c("DFC",
                "DFC (AD modules)",
                "MTG",
                "MTG (AD modules)",
                "Combined (DFC)",
                "Combined (MTG)"
                )

foutlo <- mapply(\(x, c){
  y <- -1
  temp <- x |> filter(Direction == y) 
  if(nrow(temp) > 0){
    temp1 <- split(temp[,1], temp[,2])
    tempgenes <- lapply(temp1, \(z){
      mods[z] |> unlist() |> unique()
    })
    return(data.frame("Celltype1" = names(tempgenes),
                      "num_sig" = unlist(lapply(tempgenes, length))) |>
           mutate(Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass))) |>
           complete(Celltype1) |>
           mutate(comp = c))
  } else {
    return(data.frame(Celltype1 =  unique(class_info$Subclass), 
                  num_sig = NA,
                  comp = c))
  }
}, flist1, names(flist1), SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  left_join(class_info, by = join_by(Celltype1 == Subclass)) |>
  mutate(comp = factor(comp, levels = rev(typevec)),
         Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass)))


fouthi <- mapply(\(x, c){
  y <- 1
  temp <- x |> filter(Direction == y) 
  if(nrow(temp) > 0){
    temp1 <- split(temp[,1], temp[,2])
    tempgenes <- lapply(temp1, \(z){
      mods[z] |> unlist() |> unique()
    })
    return(data.frame("Celltype1" = names(tempgenes),
                      "num_sig" = unlist(lapply(tempgenes, length))) |>
        mutate(Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass))) |>
        complete(Celltype1) |>
        mutate(comp = c))
  } else{
    return(data.frame(Celltype1 =  unique(class_info$Subclass), 
                      num_sig = NA,
                      comp = c))
  }
}, flist1, names(flist1), SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  left_join(class_info, by = join_by(Celltype1 == Subclass)) |>
  mutate(comp = factor(comp, levels = rev(typevec)),
         Celltype1 = factor(Celltype1, levels = unique(class_info$Subclass)))

plo <- foutlo |>
    ggplot(aes(x = Celltype1, y = comp, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 7), 
            legend.direction = "horizontal",    
            legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = unit(4, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(x = "", y = "") +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_A_genes.pdf"), height = 2.6, width = 4.6)
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_A_genes.svg"), height = 2.6, width = 4.6)

phi <- fouthi |>
    ggplot(aes(x = Celltype1, y = comp, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 7), 
            legend.direction = "horizontal",    
            legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = unit(4, "mm"),
            legend.title = element_blank(),
            legend.margin = margin(-0.5, 0, 0, 0, "cm"),
            legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
    #  guides(fill = guide_legend(ncol = 1),
     #        size = guide_legend(ncol = 1)) +
      labs(x = "", y = "") +
      scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
      scale_size_continuous(limits = c(1, max(foutlo$num_sig, na.rm = T)), # Manually set the scale to be the same as lo object
                            breaks = scales::pretty_breaks(n = 2)) 
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_B_genes.pdf"), height = 2.6, width = 4.6)
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_B_genes.svg"), height = 2.6, width = 4.6)


# ###########
# # Panel C
# # - Gene overlap barplots (plus signficance w.r.t. hypergeometric p-value)
# ############
# # hypergeom function:
# fisherTest_modoverlap <- function(dat1, # ex. SEAAD2024
#                                   dat2, # ex. MIT
#                                   shared.in.mod, # intersection between SEAAD2024 and MIT
#                                   all){ # Fisher's test function
#   #total.shared = length(intersect(all,dat1))
#   #shared.in.mod = length(intersect(dat2,dat1))
#   shared.out.mod = length(dat1) - shared.in.mod
#   in.mod.not.shared = length(dat2) - shared.in.mod
#   out.mod.not.shared = length(all) - length(dat2) - shared.out.mod
#   fisher.test(matrix(c(shared.in.mod,
#                       in.mod.not.shared,
#                       shared.out.mod,
#                       out.mod.not.shared), ncol=2),alternative="greater")$p.val
# }

# # Load module data (topmodposbc - all dCoPA analyses done using topmodposbc definitions)
# module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
# datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
# if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
# mods <- tapply(datkme[,2], datkme[,3], list)
# modulelengths <- unlist(lapply(mods,length))
# filter_under <- 3
# these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])

# # Load individual output tables
# # con vs All DFC
# olist <- list("g_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv")),
#               "l_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_PFC/euclidean_distances/output_table_Subclass.csv")),
#               "g_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv")),
#               "l_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_MTC/euclidean_distances/output_table_Subclass.csv")),
#               "g_r_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv")),
#               "l_r_dfc" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv")),
#               "g_r_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv")),
#               "l_r_mtg" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv"))
#               ) |> 
#   lapply(\(x){
#     x |> dplyr::filter(sig_FDR, Consistency %in% c(0, 1))
#   })

# # Convert Liu celltypes to Gabitto celltypes
# olist[c(2,4,6,8)] <- olist[c(2,4,6,8)] |>
#   lapply(\(x){
#     x |> 
#       mutate(Celltype = case_match(Celltype,
#         c("SMC", "VLMC", "End", "Per") ~ "Endothelial",
#         c("Exc L4-5 IT-2", "Exc L3-4 IT","Exc L4-5 IT-1") ~ "L4 IT",
#         "Exc L5 ET" ~ "L5 ET",
#         c("Exc L4-5 IT-2", "Exc L4-5 IT-1","Exc L3-5 IT", "Exc L5-6 IT") ~ "L5 IT",
#         "Exc L5/6 NP" ~ "L5/6 NP",
#         "Inh LAMP5" ~ "Lamp5",
#         "Inh PVALB" ~ "Pvalb",
#         "Inh SST" ~ "Sst",
#         "Exc L5-6 IT" ~ "L6 IT",
#         "Exc L5/6 IT Car3" ~ "L6 IT Car3",
#         "Exc L6 CT" ~ "L6 CT",
#         "Inh PAX6" ~ "Pax6",
#         "Ast" ~ "Astrocyte",
#         "Inh VIP" ~ "Vip",
#         "Exc L6b" ~ "L6b",
#         "Exc L2-3 IT" ~ "L2/3 IT",
#         "Exc L6 IT" ~ "L6 IT",
#         "Mic" ~ "Microglia-PVM",
#         "Oli" ~ "Oligodendrocyte",
#         .default = Celltype
#       )) |>
#       filter(!Celltype %in% c("T", "Fib", "Exc EC"))
#   })

# comp_ind <- list(c(1,2), c(3,4), c(5,6), c(7,8))
# comp_names <- c("DFC",
#                 "MTG",
#                 "DFC (ROSMAP)",
#                 "MTG (ROSMAP)")

# outdf <- lapply(seq_along(comp_ind), \(x){ # For each of the 4 comparisons,
#   a1 <- olist[[comp_ind[[x]][1]]]
#   a2 <- olist[[comp_ind[[x]][2]]]
#   ctout <- lapply(class_info$Subclass, \(ct){ # For each celltype,
#     # Find union of mods (per direction)
#     dirout <- lapply(c(1, -1), \(dir){ # For each direction,
#       m1 <- a1 |> filter(Celltype == ct,
#                          Direction == dir)
#       g1 <- unique(unlist(mods[m1$mod]))
#       m2 <- a2 |> filter(Celltype == ct,
#                          Direction == dir)
#       g2 <- unique(unlist(mods[m2$mod]))
  
#       # # Find intersection of mods (per direction)
#       # i1 <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype, Direction, Consistency), relationship = "many-to-many") |>
#       #   filter(Celltype == ct, Direction == dir) 
#       # ig1 <- unique(unlist(mods[m1$mod]))

#       ig <- length(unique(intersect(g1, g2)))
#       fisher_pval <- fisherTest_modoverlap(g1, 
#                                            g2, 
#                                            ig, 
#                                            all = 1:18913)
#       return(data.frame("comp" = comp_names[x],
#                         "ct" = ct,
#                         "direction" = dir,
#                         "length" = ig,
#                         "pval" = fisher_pval))
#     }) |> do.call(what = "rbind")
#     return(dirout)
#   }) |> do.call(what = "rbind")
#   return(ctout)
# }) |> do.call(what = "rbind") |>
#   left_join(class_info, by = join_by(ct == Subclass)) |>
#   mutate(
#     label = case_when(
#       pval > 0.05 ~ "",
#       (0.01 < pval & pval <= 0.05) ~ "*",
#       (0.001 < pval & pval <= 0.01) ~ "**",
#       (0.0001 < pval & pval <= 0.001) ~ "***",
#       pval <= 0.0001 ~ "****",
#       .default = NA
#     ),
#     ct = factor(ct, levels = class_info$Subclass),
#     comp = factor(comp, levels = comp_names)) 


# # Plot (higher in disease)
# p <- outdf |>
#   filter(direction == -1) |>
#   ggplot(aes(x = ct, y = length, fill = Class)) +
#     theme_classic() + 
#     geom_bar(stat = "identity") + 
#     geom_text(aes(label = label), 
#               #y = max(outdf$length[outdf$direction == -1]), 
#               vjust = -0.5) + 
#     facet_wrap(~comp, nrow = 4) +
#     ylim(0, max(outdf$length[outdf$direction == -1]) * 1.2) +
#     labs(x = "", y = "# of overlapping module genes") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#           legend.title = element_blank(),
#           legend.position = "bottom",
#           strip.text = element_blank(),
#       #    strip.background = element_blank()
#       ) +
#     scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) 
# ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v1/panel_C1.pdf"), height = 5, width = 7)

# p <- outdf |>
#   filter(direction == 1) |>
#   ggplot(aes(x = ct, y = length, fill = Class)) +
#     theme_classic() + 
#     geom_bar(stat = "identity") + 
#     geom_text(aes(label = label), 
#               #y = max(outdf$length[outdf$direction == 1]), 
#               vjust = -0.5) + 
#     facet_wrap(~comp, nrow = 4) +
#     ylim(0, max(outdf$length[outdf$direction == 1]) * 1.2) +
#     labs(x = "", y = "# of overlapping module genes") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#           legend.title = element_blank(),
#           legend.position = "bottom",
#           strip.text = element_blank(),
#        #   strip.background = element_blank()
#        ) +
#     scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) 
# ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v1/panel_C2.pdf"), height = 5, width = 7)

###########
# Panel C new
# - Gene overlap barplots (plus signficance w.r.t. hypergeometric p-value)
############
# hypergeom function:
fisherTest_modoverlap <- function(dat1, # ex. SEAAD2024
                                  dat2, # ex. MIT
                                  shared.in.mod, # intersection between SEAAD2024 and MIT
                                  all){ # Fisher's test function
  #total.shared = length(intersect(all,dat1))
  #shared.in.mod = length(intersect(dat2,dat1))
  shared.out.mod = length(dat1) - shared.in.mod
  in.mod.not.shared = length(dat2) - shared.in.mod
  out.mod.not.shared = length(all) - length(dat2) - shared.out.mod
  fisher.test(matrix(c(shared.in.mod,
                      in.mod.not.shared,
                      shared.out.mod,
                      out.mod.not.shared), ncol=2),alternative="greater")$p.val
}



dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/")

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

dat_vec <- c("gabitto_vs_liu_MTG_AllADVsCon",         
             "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP",
             "gabitto_vs_liu_DFC_AllADVsCon",    
             "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP"
             ) 

title_vec <- c("Con vs all AD (MTG)",
               "Con vs all AD (ROSMAP, MTG)",
               "Con vs all AD (DFC)",
               "Con vs all AD (ROSMAP, DFC)"
               )

label_vec <- c("Shared in Gabitto and Liu MTG (con vs all AD)",           
               "Shared in Gabitto and Liu MTG (con vs all AD, ROSMAP mods)",
               "Shared in Gabitto and Liu DFC (con vs all AD)",    
               "Shared in Gabitto and Liu DFC (con vs all AD, ROSMAP mods)"
               ) 

dat_list <- mapply(\(x, y){
  out <- fread(file.path(dir1, x, "shared_output_table.csv"), data.table = F) |>
    filter(mod %in% these_mods)
  return(out)
}, dat_vec, title_vec, SIMPLIFY = F)
names(dat_list) <- label_vec

comp_names <- c("MTG", "DFC")

genedf <- lapply(seq_along(comp_names), \(i){
  d1 <- dat_list[[i*2 - 1]]
  d2 <- dat_list[[i*2]]
  #common_cts <- unique(c(d1[,2], d2[,2]))
  common_cts <- class_info$Subclass

  outdf <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(common_cts)){
      m1 <- d1 |> filter(Celltype1 == common_cts[c],
                        Direction == d)
      g1 <- unique(unlist(mods[m1$mod]))
      m2 <- d2 |> filter(Celltype1 == common_cts[c],
                        Direction == d)
      g2 <- unique(unlist(mods[m2$mod]))

      ngene <- 18913 * 24 * 2
      fisher_pval <- fisherTest_modoverlap(g1, 
                                           g2, 
                                           length(unique(intersect(g1, g2))), # intersection 
                                           all = 1:ngene)

      sublist[[c]] <- data.frame("Comparison" = comp_names[i],
                                 "Direction" = d,
                                 "Celltype" = common_cts[c],
                                 "Length1" = length(g1),
                                 "Length2" = length(g2),
                                 "NMods1" = length(unique(m1$mod)),
                                 "NMods2" = length(unique(m2$mod)),
                                 "LengthIntersectGene" = length(unique(intersect(g1, g2))),
                                 "LengthUnionGene" = length(unique(union(g1, g2))), 
                                 "LengthIntersectMod" = length(unique(intersect(m1$mod, m2$mod))),
                                 "Fisher_p_val" = fisher_pval)
    }
    return(sublist |> do.call(what = "rbind"))
  }) |> do.call(what = "rbind") |>
    mutate(log_pval = -log10(Fisher_p_val),
           log_pval = case_match(log_pval, Inf ~ 324, .default = log_pval),
           Direction = factor(Direction, levels = c("-1", "1")),
           Direction = fct_recode(Direction,
                                  "Lower expression in disease" = "-1",
                                  "Higher expression in disease" = "1"),
           Celltype = factor(Celltype, levels = class_info$Subclass)
           )
  return(outdf)
}) |> do.call(what = "rbind") |>
  left_join(class_info, by = join_by(Celltype == Subclass)) |>
  mutate(
    label = case_when(
      Fisher_p_val > 0.05 ~ "",
      (0.01 < Fisher_p_val & Fisher_p_val <= 0.05) ~ "*",
      (0.001 < Fisher_p_val & Fisher_p_val <= 0.01) ~ "**",
      (0.0001 < Fisher_p_val & Fisher_p_val <= 0.001) ~ "***",
      Fisher_p_val <= 0.0001 ~ "****",
      .default = NA
    ),
    Celltype = factor(Celltype, levels = class_info$Subclass),
    Comparison = factor(Comparison, levels = comp_names))

p <- genedf |>
  filter(Direction == "Lower expression in disease") |>
  ggplot(aes(x = Celltype, y = LengthIntersectGene, fill = Class)) +
    theme_bw() + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = label), 
              #y = max(outdf$length[outdf$direction == -1]), 
              vjust = -0.5) + 
    facet_wrap(~Comparison, nrow = 4) +
    ylim(0, max(genedf$LengthIntersectGene[genedf$Direction == "Lower expression in disease"]) * 1.3) +
    labs(x = "", y = "# of overlapping module genes") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
       #   strip.text = element_blank(),
      #    strip.background = element_blank()
      ) +
    scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) 
ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C1.pdf"), height = 4, width = 7)
ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C1.svg"), height = 4, width = 7)
ggsave(p + theme(legend.position = "none"), 
       file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C1_nolegend.pdf"), height = 3.5, width = 7)
ggsave(p + theme(legend.position = "none"), 
       file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C1_nolegend.svg"), height = 3.5, width = 7)


p <- genedf |>
  filter(Direction == "Higher expression in disease") |>
  ggplot(aes(x = Celltype, y = LengthIntersectGene, fill = Class)) +
    theme_bw() + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = label), 
              #y = max(outdf$length[outdf$direction == -1]), 
              vjust = -0.5) + 
    facet_wrap(~Comparison, nrow = 4) +
    ylim(0, max(genedf$LengthIntersectGene[genedf$Direction == "Higher expression in disease"]) * 1.3) +
    labs(x = "", y = "# of overlapping module genes") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
       #   strip.text = element_blank(),
      #    strip.background = element_blank()
      ) +
    scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) 
ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C2.pdf"), height = 4, width = 7)
ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C2.svg"), height = 4, width = 7)
ggsave(p + theme(legend.position = "none"), 
       file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C2_nolegend.pdf"), height = 3.5, width = 7)
ggsave(p + theme(legend.position = "none"), 
       file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_C2_nolegend.svg"), height = 3.5, width = 7)



############
# Upset plot
############

# DFC modules
d1 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon/shared_output_table.csv"))|>
    filter(mod %in% these_mods)
d2 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon_ROSMAP/shared_output_table.csv"))|>
    filter(mod %in% these_mods)

common_cts <- unique(c(d1[,2], d2[,2]))

# Collect lists of mods per subclass (DFC)
outlist <- lapply(c(-1, 1), \(d){
  sublist <- list()
  for(c in seq_along(common_cts)){
    m1 <- d1 |> filter(Celltype1 == common_cts[c],
                      Direction == d)
    g1 <- unique(unlist(mods[m1$mod]))
    m2 <- d2 |> filter(Celltype1 == common_cts[c],
                      Direction == d)
    g2 <- unique(unlist(mods[m2$mod]))
    
    outvec <- unique(intersect(g1, g2))
    if(length(outvec) > 0){
      sublist[[c]] <- outvec
    } else {
      sublist[[c]] <- "none"
    }
  }
  names(sublist) <- common_cts
  sublist <- sublist[-which(unlist(lapply(sublist, \(x) "none" %in% x)))]
  return(sublist)
}) 
names(outlist) <- c(-1, 1)

# Code for only producing the top barplot of UpSetR, due to the package plotting the barplot overlapping with the upset portion
# Convert outlist to the "expression" input a la UpSetR
# get_disjoint_combinations <- function(sets_list) {
  
#   set_names <- names(sets_list)
#   n         <- length(set_names)
#   all_combos <- list()
  
#   # Generate all non-empty combinations of set names
#   for (i in 1:n) {
#     combos <- combn(set_names, i, simplify = FALSE)
#     all_combos <- c(all_combos, combos)
#   }
  
#   # For each element, find exactly which sets it belongs to
#   all_elements <- unique(unlist(sets_list))
  
#   element_membership <- sapply(all_elements, function(el) {
#     paste(set_names[sapply(sets_list, function(s) el %in% s)], collapse = "&")
#   })
  
#   # Count how many elements belong to each exact combination
#   combo_counts <- table(element_membership)
  
#   # Convert to named numeric vector
#   result <- as.numeric(combo_counts)
#   names(result) <- names(combo_counts)
  
#   return(sort(result, decreasing = TRUE))
# }

# outlistvec <- get_disjoint_combinations(outlist[[1]])

# # Create top half of upset plot in ggplot
# p <- data.frame("which" = names(outlistvec), "count" = outlistvec) |>
#   arrange(desc(count)) |>
#   mutate(which = factor(which, levels = which)) |>
#   ggplot(aes(x = which, y = count)) +
#     theme_minimal() + 
#     theme(text = element_text(size = 12),
#           axis.title.y = element_text(size = 10),
#           axis.text.x = element_blank(),
#           legend.position = "none",
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           plot.margin = margin(1, 0, 0, 0, "cm"),
#           axis.line.y = element_line(colour = "black", linewidth = 0.5),
#           axis.line.x = element_line(colour = "black", linewidth = 0.5),
#           axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
#     geom_bar(stat = 'identity', width = 0.5, fill = "black") +
#     #geom_text(aes(label = count), vjust = -0.5) +
#     labs(x = "", y = "# of overlaps") +
#     scale_x_discrete(expand = c(0, 1)) +
#     scale_y_continuous(limits = c(0, max(outlistvec)), expand = c(0, 1)) 
# ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_D1_top.svg"), bg = "white", width = 4, height = 2)
# ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_D1_top.pdf"), bg = "white", width = 4, height = 2)

# New solution: just save the plot as a pdf and edit the overlap in Inkscape (the "missing" portion of the graph is underneath the upset portion)

# DFC (lower in disease)
p1 <- upset(fromList(outlist[[1]]), 
            sets = names(outlist[[1]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlaps", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.5, 0.5),
            show.numbers = F)

svg(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_D1.svg"), 
    width = 3.5, height = 2)
p1 
#grid.text("All AD vs Con, DFC, lower in severe", x = 0.65, y=0.98, gp=gpar(fontsize=8))
dev.off()

pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_D1.pdf"), 
    width = 3.5, height = 2)
p1 
#grid.text("All AD vs Con, DFC, lower in severe", x = 0.65, y=0.98, gp=gpar(fontsize=8))
dev.off()

# DFC (higher in disease)
# Only astrocytes

# MTG modules
d3 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon/shared_output_table.csv"))|>
    filter(mod %in% these_mods)
d4 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon_ROSMAP/shared_output_table.csv"))|>
    filter(mod %in% these_mods)

common_cts <- unique(c(d3[,2], d4[,2]))

# Collect lists of mods per subclass (DFC)
outlist <- lapply(c(-1, 1), \(d){
  sublist <- list()
  for(c in seq_along(common_cts)){
    m1 <- d3 |> filter(Celltype1 == common_cts[c],
                      Direction == d)
    g1 <- unique(unlist(mods[m1$mod]))
    m2 <- d4 |> filter(Celltype1 == common_cts[c],
                      Direction == d)
    g2 <- unique(unlist(mods[m2$mod]))
    
    outvec <- unique(intersect(g1, g2))
    if(length(outvec) > 0){
      sublist[[c]] <- outvec
    } else {
      sublist[[c]] <- "none"
    }
  }
  names(sublist) <- common_cts
  sublist <- sublist[-which(unlist(lapply(sublist, \(x) "none" %in% x)))]
  return(sublist)
}) 
names(outlist) <- c(-1, 1)

# outlistvec <- get_disjoint_combinations(outlist[[1]])

# p <- data.frame("which" = names(outlistvec), "count" = outlistvec) |>
#   arrange(desc(count)) |>
#   mutate(which = factor(which, levels = which)) |>
#   ggplot(aes(x = which, y = count)) +
#     theme_minimal() + 
#     theme(text = element_text(size = 12),
#           axis.title.y = element_text(size = 10),
#           axis.text.x = element_blank(),
#           legend.position = "none",
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           plot.margin = margin(1, 0, 0, 0, "cm"),
#           axis.line.y = element_line(colour = "black", linewidth = 0.5),
#           axis.line.x = element_line(colour = "black", linewidth = 0.5),
#           axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
#     geom_bar(stat = 'identity', width = 0.5, fill = "black") +
#     #geom_text(aes(label = count), vjust = -0.5) +
#     labs(x = "", y = "# of overlaps") +
#     scale_x_discrete(expand = c(0, 1)) +
#     scale_y_continuous(limits = c(0, max(outlistvec)), expand = c(0, 1)) 
# ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_D2_top.svg"), bg = "white", width = 4, height = 2)
# ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_D2_top.pdf"), bg = "white", width = 4, height = 2)


# MTG (lower in disease)
p1 <- upset(fromList(outlist[[1]]), 
            sets = names(outlist[[1]]), 
            order.by = "freq",
            mainbar.y.label = "# of overlaps", 
            sets.x.label = "# of genes",
            mb.ratio = c(0.4, 0.6),
            show.numbers = F)

svg(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_D2.svg"), 
    width = 4.7, height = 2.7)
p1 
dev.off()

pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_D2.pdf"), 
    width = 4.7, height = 2.7)
p1 
dev.off()

# MTG (higher in disease)
# None

#############
# Map conserved dCoPA modules onto metaclusters
#############

# Clustering of eigenmodules (Jorstad consensusMin)
metacluster_min <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/jorstad_consensusMin_dendro.qs"))

# Eigenmodules (already in the order of metacluster_min)
mod_eig <- fread(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/mod_eig.csv")), data.table = F) |>
  column_to_rownames("V1")

# kME table of eigenmodules
kme <- fread(data.table = F, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v2/Jorstad_DFC_consensusMin_noMerge_Modules/Pearson-no_TO_signum0.703_minSize3_merge_ME_1_1021/kME_table_.csv"))) |>
  mutate(Gene = gsub("Mod", "", Gene) |> as.integer())

# Create mapping table for old to new module indices
modindmap <- data.frame("new" = 1:1023,
                        "old" = these_mods)

# Find list of shared dCoPA mods between bulk mod output and ROSMAP mod output
# AllADVsCon DFC 
dcopa_shared <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon/shared_output_table.csv"))
dcopa_shared_r <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon_ROSMAP/shared_output_table.csv"))
dsr <- dplyr::inner_join(dcopa_shared, dcopa_shared_r, by = dplyr::join_by(mod, Celltype1, Direction, Consistency), relationship = "many-to-many") |>
  select(mod, Celltype1, Direction) |>
  filter(!duplicated(paste0(mod, Celltype1, Direction))) |>
  left_join(modindmap, by = join_by(mod == old)) |> # map old indices (1158) to new indices (1023)
  select(!mod) |> 
  left_join(kme[,c(1,6)], by = join_by(new == Gene)) |> # Add metamodule information (topmodposfdr)
  group_by(Celltype1, TopModPosFDR_0.0782) |>
  summarise(n = n()) |>
  mutate(TopModPosFDR_0.0782 = factor(TopModPosFDR_0.0782, levels = colnames(mod_eig))) |> # Order by eigenmodule dendrogram order
  arrange(TopModPosFDR_0.0782) |> 
  pivot_wider(names_from = TopModPosFDR_0.0782, values_from = n, 
              names_expand = T, values_fill = NA) |> # create wide matrix of celltype by metamodule
  column_to_rownames("Celltype1") 

dsr[is.na(dsr)] <- 0
fwrite(dsr, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/metamod_proj_dfc_alladvscon.csv"), row.names = T)


# p <- Heatmap(dsr,
#              name = "# of modules",
#              cluster_columns = metacluster_min
#              )

# svg(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v2/metamodule_projection_heatmap.svg"))
# draw(p)
# dev.off()

# AllADVsCon MTG
dcopa_shared <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon/shared_output_table.csv"))
dcopa_shared_r <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon_ROSMAP/shared_output_table.csv"))
dsr <- dplyr::inner_join(dcopa_shared, dcopa_shared_r, by = dplyr::join_by(mod, Celltype1, Direction, Consistency), relationship = "many-to-many") |>
  select(mod, Celltype1, Direction) |>
  filter(!duplicated(paste0(mod, Celltype1, Direction))) |>
  left_join(modindmap, by = join_by(mod == old)) |> # map old indices (1158) to new indices (1023)
  select(!mod) |> 
  left_join(kme[,c(1,6)], by = join_by(new == Gene)) |> # Add metamodule information (topmodposfdr)
  group_by(Celltype1, TopModPosFDR_0.0782) |>
  summarise(n = n()) |>
  mutate(TopModPosFDR_0.0782 = factor(TopModPosFDR_0.0782, levels = colnames(mod_eig))) |> # Order by eigenmodule dendrogram order
  arrange(TopModPosFDR_0.0782) |> 
  pivot_wider(names_from = TopModPosFDR_0.0782, values_from = n, 
              names_expand = T, values_fill = NA) |> # create wide matrix of celltype by metamodule
  column_to_rownames("Celltype1") 
dsr[is.na(dsr)] <- 0
fwrite(dsr, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/metamod_proj_mtg_alladvscon.csv"), row.names = T)

##############
# GSEA 
#############

# Different strategy:
# Get all L4 IT genes per comparison, run GSEA

dat_vec <- list(list("gabitto_vs_liu_DFC_AllADVsCon",  
                     "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP"),    
                list("gabitto_vs_liu_DFC_earlyVsCon",       
                     "gabitto_vs_liu_DFC_earlyVsCon_ROSMAP"),  
                list("gabitto_vs_liu_DFC_lateVsEarly",        
                     "gabitto_vs_liu_DFC_lateVsEarly_ROSMAP"),
                list("gabitto_vs_liu_MTG_AllADVsCon",  
                     "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP"),          
                list("gabitto_vs_liu_MTG_earlyVsCon",  
                     "gabitto_vs_liu_MTG_earlyVsCon_ROSMAP"),  
                list("gabitto_vs_liu_MTG_lateVsEarly",        
                     "gabitto_vs_liu_MTG_lateVsEarly_ROSMAP")#,
                #"gabitto_vs_liu_APOE_DFC",
                #"brainSCOPE_CMC_vs_SZBD"
                ) 

dat_vec_names <- c("AllADVsCon_DFC",
                   "earlyVsCon_DFC",
                   "lateVsEarly_DFC",
                   "AllADVsCon_MTG",
                   "earlyVsCon_MTG",
                   "lateVsEarly_MTG")

# Load all relevant dcopa output tables (for all comparisons)
dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared")
dcopa_shared <- lapply(dat_vec, \(x){
  lapply(x, \(y){
    fread(file.path(dir1, y, "shared_output_table.csv"), data.table = F) |>
      filter(mod %in% these_mods)
  })
})

# Function for extracting genes belonging to shared significant dCoPA mods (topmodposbc)
temp_func <- function(d1, d2, 
                      common_cts = NULL,
                      cts){
  if(is.null(common_cts)){
    common_cts <- unique(c(d1[,2], d2[,2]))
  }
  outlist <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(common_cts)){
      m1 <- d1 |> filter(Celltype1 == common_cts[c],
                        Direction == d)
      g1 <- unique(unlist(mods[m1$mod]))
      m2 <- d2 |> filter(Celltype1 == common_cts[c],
                        Direction == d)
      g2 <- unique(unlist(mods[m2$mod]))
      
      outvec <- unique(intersect(g1, g2))
      if(length(outvec) > 0){
        sublist[[c]] <- outvec
      } else {
        sublist[[c]] <- "none"
      }
      sublist <- sublist[common_cts == cts]
    }
    names(sublist) <- cts
    #names(sublist) <- common_cts
    #sublist <- sublist[-which(unlist(lapply(sublist, \(x) "none" %in% x)))]
    return(sublist)
  }) 
  names(outlist) <- c(-1, 1)
  return(outlist)
}

# Using function, create list of genes for all celltypes
dcopa_allct <- lapply(dcopa_shared, \(x){
   temp_func(d1 = x[[1]], d2 = x[[2]], 
             common_cts = unique(class_info$Subclass),
             cts = unique(class_info$Subclass))
}) |> setNames(dat_vec_names)
# Organize list of genes into dataframe and save
dcopa_to_df <- function(dcopa_allct) {                                         
  dir_labels <- c("-1" = "Lower in more severe", "1" = "Higher in more severe")

  lapply(names(dcopa_allct), \(comp) {
    lapply(names(dcopa_allct[[comp]]), \(dir) {
      ct_list <- dcopa_allct[[comp]][[dir]]
      lapply(names(ct_list), \(ct) {
        genes <- ct_list[[ct]]
        if (is.null(genes) || identical(genes, "none")) return(NULL)
        data.frame(
          Comparison = comp,
          Direction  = dir_labels[[dir]],
          Celltype   = ct,
          Gene       = genes,
          stringsAsFactors = FALSE
        )
      }) |> do.call(what = "rbind")
    }) |> do.call(what = "rbind")
  }) |> do.call(what = "rbind")
}

gsea_input <- dcopa_to_df(dcopa_allct)
fwrite(gsea_input, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/Allct_genes_allComparisons.csv"))

# Using function, create list of genes for L4 IT
dcopa_l4it <- lapply(dcopa_shared, \(x){
   temp_func(d1 = x[[1]], d2 = x[[2]], cts = c("L4 IT"))
}) |> setNames(dat_vec_names)

# Organize list of genes into dataframe and save
gsea_input <- mapply(\(x, y){
  t1 <- data.frame("Comparison" = y,
                   "Direction" = "Lower in more severe",
                   "Genes" = if(is.null(x[[1]][[1]])) NA else x[[1]][[1]])
  t2 <- data.frame("Comparison" = y,
                   "Direction" = "Higher in more severe",
                   "Genes" = if(is.null(x[[2]][[1]])) NA else x[[2]][[1]])    
  return(rbind(t1, t2))
}, dcopa_l4it, dat_vec_names, SIMPLIFY = F) |>
  do.call(what = "rbind")

# Write genes for all comparisons
fwrite(gsea_input, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/L4IT_genes_allComparisons.csv"))
# Write genes for L4 IT, AllADVsCon_MTG
fwrite(gsea_input |> 
        filter(Comparison == "AllADVsCon_MTG",
        Direction == "Lower in more severe"),
        file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/L4IT_genes_AllADVsCon_MTG.csv"))

# Run GSEA on all AD vs Con, MTG
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func.R"))

gsea_input1 <- gsea_input |>
  filter(Comparison == "AllADVsCon_MTG",
        Direction == "Lower in more severe") |>  
  pull(Genes) 

b_out <- run_gsea_for_proj(list("L4 IT" = gsea_input1),
                           set_list=NULL,
                           file_desc=NULL,
                           broad = T,
                           save_dir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)
                           )

g_out <- run_gsea_for_proj(list("L4 IT" = gsea_input1),
                           set_list=NULL,
                           file_desc=NULL,
                           broad = F,
                           save_dir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)
                           )

# Organize, filter, and save 
all_out <- rbind(
  b_out |> 
    arrange(L4.IT) |>
    filter(L4.IT < 0.05),
  g_out |> 
    arrange(L4.IT) |>
    filter(L4.IT < 0.05)
)

fwrite(all_out, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v2/gsea_allADvsconMTG.csv"))


#### Find overlap between L4 IT dcopa mods in...
# - lower in severe, early vs con, MTG
# - lower in severe, all AD vs con, MTG
# - lower in severe, late vs early, MTG
# - the same in DFC

str(dcopa_l4it, max.level = 3)
length(intersect(dcopa_l4it$lateVsEarly_MTG$`-1`$`L4 IT`,
                 dcopa_l4it$AllADVsCon_MTG$`-1`$`L4 IT`))
# 239 genes intersect between lower in late AD and lower in all AD, MTG
# There are no overlapping genes for L4 IT for lower in early AD, MTG or DFC

#### Find overlap between SCZ and AD for deep layer neurons

# Load SCZ modules
temp_func2 <- function(d1, d2, cts){
  common_cts <- unique(c(d1[,2], d2[,2]))
  outlist <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(common_cts)){
      m1 <- d1 |> filter(Celltype == common_cts[c],
                        Direction == d)
      g1 <- unique(unlist(mods[m1$mod]))
      m2 <- d2 |> filter(Celltype == common_cts[c],
                        Direction == d)
      g2 <- unique(unlist(mods[m2$mod]))
      
      outvec <- unique(intersect(g1, g2))
      if(length(outvec) > 0){
        sublist[[c]] <- outvec
      } else {
        sublist[[c]] <- "none"
      }
      sublist <- sublist[common_cts == cts]
    }
    names(sublist) <- cts
    #names(sublist) <- common_cts
    #sublist <- sublist[-which(unlist(lapply(sublist, \(x) "none" %in% x)))]
    return(sublist)
  }) 
  names(outlist) <- c(-1, 1)
  return(outlist)
}

d1 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/brainSCOPE_CMC_vs_SZBD/shared_output_table.csv"))|>
    filter(mod %in% these_mods)
d2 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/brainSCOPE_CMC_vs_SZBD_ROSMAP/shared_output_table.csv"))|>
    filter(mod %in% these_mods)

scz_mods <- temp_func2(d1 = d1, d2 = d2, cts = "L4 IT")
# Seems like there are no genes shared between ROSMAP AD and bulk megaset projections

# What if we just compared CMC L4IT genes to L4IT genes from AD
d1m <- d1 |> filter(Direction == -1, Celltype == "L4 IT") |>
  pull(mod)
g1 <- unique(unlist(mods[d1m]))
mtgscz1 <- g1[g1 %in% dcopa_l4it[[4]][[1]][[1]]]
# length(mtgscz1)
# 31 shared with AllADVsCon_MTG
# > mtgscz1
#  [1] "PTCD2"    "NRDC"     "STIM2"    "UBE4A"    "KPNA1"    "ARMC8"   
#  [7] "FBXW2"    "USP9X"    "ECPAS"    "ANAPC1"   "PPP1R21"  "PDCD6IP" 
# [13] "KIAA0232" "PIK3R4"   "ZNF407"   "AKAP11"   "WASF1"    "SGIP1"   
# [19] "PPP6C"    "RAB22A"   "NDFIP1"   "KRAS"     "ACTR2"    "FAM98B"  
# [25] "PITPNB"   "LRRTM3"   "BMPR2"    "MAST4"    "SYNE1"    "NTRK3"   
# [31] "MED12L"  
dfcscz1 <- g1[g1 %in% dcopa_l4it[[1]][[1]][[1]]]
# length(dfcscz1)
# 48 shared with AllADVsCon_DFC
#  [1] "XK"       "UBE2K"    "USP14"    "WASL"     "SHOC2"    "CCDC6"   
#  [7] "CAMSAP2"  "INSIG2"   "MARCHF6"  "LRP12"    "CACNA2D1" "PPP2R5E" 
# [13] "FAM126B"  "PABIR2"   "VCPIP1"   "GPR137C"  "PPTC7"    "ERO1A"   
# [19] "TLK1"     "ZNF670"   "PTCD2"    "NRDC"     "STIM2"    "UBE4A"   
# [25] "KPNA1"    "ARMC8"    "FBXW2"    "USP9X"    "ECPAS"    "ANAPC1"  
# [31] "PPP1R21"  "PDCD6IP"  "KIAA0232" "PIK3R4"   "ZNF407"   "CCNT2"   
# [37] "HACE1"    "ZRANB2"   "UBR3"     "PHOSPHO2" "SSBP2"    "C9orf72" 
# [43] "GNPDA2"   "WDR41"    "TMEM135"  "BRWD1"    "CD47"     "ZFP28"   

# genes shared between CMC_SCZ, alladvscon_DFC, alladvscon_MTG
length(intersect(mtgscz1, dfcscz1))
# 15
# [1] "PTCD2"    "NRDC"     "STIM2"    "UBE4A"    "KPNA1"    "ARMC8"   
#  [7] "FBXW2"    "USP9X"    "ECPAS"    "ANAPC1"   "PPP1R21"  "PDCD6IP" 
# [13] "KIAA0232" "PIK3R4"   "ZNF407"  

########
# GSEA v2
########

# Strategy:
# For all datasets,
# - Create tables for all genes belonging to significant dCoPA modules
# - Run geneset enrichment per celltype and comparison

dat_vec <- list("gabitto_vs_liu_DFC_AllADVsCon",    
                "gabitto_vs_liu_DFC_earlyVsCon",       
                "gabitto_vs_liu_DFC_lateVsEarly",        
                "gabitto_vs_liu_MTG_AllADVsCon",  
                "gabitto_vs_liu_MTG_earlyVsCon",  
                "gabitto_vs_liu_MTG_lateVsEarly",        
                "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP",    
                "gabitto_vs_liu_DFC_earlyVsCon_ROSMAP",  
                "gabitto_vs_liu_DFC_lateVsEarly_ROSMAP",
                "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP",          
                "gabitto_vs_liu_MTG_earlyVsCon_ROSMAP",  
                "gabitto_vs_liu_MTG_lateVsEarly_ROSMAP",
                "brainSCOPE_CMC_vs_SZBD",
                "brainSCOPE_CMC_vs_SZBD_ROSMAP"
                ) 

dat_vec_names <- c("AllADVsCon_DFC",
                   "earlyVsCon_DFC",
                   "lateVsEarly_DFC",
                   "AllADVsCon_MTG",
                   "earlyVsCon_MTG",
                   "lateVsEarly_MTG",
                   "AllADVsCon_DFC_ROSMAP",
                   "earlyVsCon_DFC_ROSMAP",
                   "lateVsEarly_DFC_ROSMAP",
                   "AllADVsCon_MTG_ROSMAP",
                   "earlyVsCon_MTG_ROSMAP",
                   "lateVsEarly_MTG_ROSMAP",
                   "brainSCOPE_CMC_vs_SZBD_DFC",
                   "brainSCOPE_CMC_vs_SZBD_DFC_ROSMAP")

# Load all relevant dcopa output tables (for all comparisons)
dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared")
dcopa_shared <- lapply(dat_vec, \(x){
  fread(file.path(dir1, x, "shared_output_table.csv"), data.table = F)|>
    filter(mod %in% these_mods)
})

# Function for extracting genes belonging to shared significant dCoPA mods (topmodposbc)
extract_genes <- function(d1){
  cts <- unique(d1[,2])

  outlist <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(cts)){
      m1 <- d1 |> filter(Celltype1 == cts[c],
                        Direction == d)
      outvec <- unique(unlist(mods[m1$mod]))
      
      if(length(outvec) > 0){
        sublist[[c]] <- outvec
      } else {
        sublist[[c]] <- "none"
      }
    }
    names(sublist) <- cts
    return(sublist)
  }) 
  names(outlist) <- c(-1, 1)
  return(outlist)
}

# Using function, create list of genes for all celltypes
dcopa_allct <- lapply(dcopa_shared, extract_genes) |> setNames(dat_vec_names)

# Function for creating a list of all intersections from a list of genes by celltypes
# Note that this is inclusive (e.g. genes belonging to group A will represent all genes belonging to group A regardless of whether they are shared or not in group B, C, etc)
all_intersections_fast <- function(lst) {
  # Get all unique elements across all vectors
  all_elements <- unique(unlist(lst))
  n <- length(lst)
  set_names <- if (!is.null(names(lst))) names(lst) else paste0("V", seq_len(n))
  
  # Build binary membership matrix (elements x sets)
  membership <- matrix(0, nrow = length(all_elements), ncol = n,
                       dimnames = list(all_elements, set_names))
  
  for (i in seq_len(n)) {
    membership[all_elements %in% lst[[i]], i] <- 1
  }
  
  # For each unique row pattern, recover the elements and label the combination
  # Convert each row to a string key for fast grouping
  row_keys <- apply(membership, 1, paste, collapse = "")
  
  # Get all 2^n - 1 possible non-zero patterns
  result <- list()
  
  for (i in seq_len(2^n - 1)) {
    # Build the expected bit pattern for this subset
    idx <- which(as.logical(intToBits(i)[seq_len(n)]))
    pattern <- rep(0, n)
    pattern[idx] <- 1
    key <- paste(pattern, collapse = "")
    
    # Find elements that match this exact pattern
    matched_elements <- all_elements[row_keys == key]
    
    if (length(matched_elements) > 0) {
      label <- paste(set_names[idx], collapse = " & ")
      result[[label]] <- matched_elements
    }
  }
  
  result
}

# Create list of overlaps from simple lists
dcopa_allct_int <- lapply(dcopa_allct, \(x){
  lapply(x, all_intersections_fast)
})


# Organize list of genes into dataframe and save
dcopa_to_df <- function(dcopa_allct) {                                         
  dir_labels <- c("-1" = "Lower in more severe", "1" = "Higher in more severe")

  lapply(names(dcopa_allct), \(comp) {
    lapply(names(dcopa_allct[[comp]]), \(dir) {
      ct_list <- dcopa_allct[[comp]][[dir]]
      lapply(names(ct_list), \(ct) {
        genes <- ct_list[[ct]]
        if (is.null(genes) || identical(genes, "none")) return(NULL)
        data.frame(
          Comparison = comp,
          Direction  = dir_labels[[dir]],
          Celltype   = ct,
          Gene       = genes,
          stringsAsFactors = FALSE
        )
      }) |> do.call(what = "rbind")
    }) |> do.call(what = "rbind")
  }) |> do.call(what = "rbind")
}

# Save list off all genes by celltype (no intersections)
gsea_input <- dcopa_to_df(dcopa_allct) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
         "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
         Celltype = paste(Celltype, "all")) |>
  dplyr::select(all_of(c("Comparison", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))

# Save list of all genes by celltype (with intersections)
gsea_input2 <- dcopa_to_df(dcopa_allct_int) |> 
  mutate("Region" = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         "Disease" = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
         "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
         Celltype = ifelse(grepl("&", Celltype), Celltype, paste(Celltype, "unique"))) |>
  dplyr::select(all_of(c("Comparison", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))

# Combine and save all
gsea_input_all <- rbind(gsea_input, gsea_input2) |>
  arrange(Disease, Comparison)
fwrite(gsea_input_all, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/dcopa_genelist.csv"))

# Run gsea
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
ginputlist1 <- split(gsea_input_all$Gene, paste(gsea_input_all$Comparison, gsea_input_all$Direction, gsea_input_all$Celltype))
b_out <- run_gsea_for_proj_optimized(ginputlist1,
                                     set_list=NULL,
                                     file_desc=NULL,
                                     broad = T,
                                     save_dir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)
                                     )

g_out <- run_gsea_for_proj_optimized(ginputlist1,
                           set_list=NULL,
                           file_desc=NULL,
                           broad = F,
                           save_dir = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)
                           )
all_out <- rbind(b_out, g_out)
fwrite(all_out, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/gsea_all.csv"))

# DE results in Gabitto:
# ex. fig 8
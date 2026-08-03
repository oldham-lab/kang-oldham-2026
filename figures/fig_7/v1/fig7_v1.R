library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)

version <- "v1"
if(!dir.exists(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)))
  dir.create(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version), recursive = T)

##########
# Panel A-B
# Dotplots (Control vs all AD)
##########
# Find overlap between all AllADvsCon comparisons
file1 <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared"), full.names = T)
file1 <- file1[c(3:4, 10:11)]
flist <- lapply(file1, \(x) fread(data.table = F, file = file.path(x, "shared_output_table.csv")))
# Combined (Bulk megaset modules)
flistjoin <- Reduce(\(x, y){
  a1 <- x |> select(c("mod", "Celltype_Gabitto_2024", "Direction", "Consistency")) 
  a2 <- y |> select(c("mod", "Celltype_Gabitto_2024", "Direction", "Consistency"))
  out <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype_Gabitto_2024, Direction, Consistency), relationship = "many-to-many") 
  return(out)
}, flist[c(1,3)]) |>
  filter(!duplicated(paste0(mod, Celltype_Gabitto_2024, Direction, Consistency))) |>
    mutate(Direction = factor(Direction, levels = c(1, -1))) |> 
    rename("Celltype" = "Celltype_Gabitto_2024") |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    select(mod, Celltype, Direction) |>
    group_by(Celltype, Direction, .drop = FALSE) |>
    summarise(num_sig = n(), .groups = "drop") |>
    left_join(class_info, by = join_by(Celltype == Subclass)) |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    arrange(Celltype) |>
    mutate(num_sig = case_match(num_sig, 
                                0 ~ NA,
                                .default = num_sig
                                ),
           comp = "Combined") 
# Combined (ROSMAP modules)
flistjoinros <- Reduce(\(x, y){
  a1 <- x |> select(c("mod", "Celltype_Gabitto_2024", "Direction", "Consistency")) 
  a2 <- y |> select(c("mod", "Celltype_Gabitto_2024", "Direction", "Consistency"))
  out <- dplyr::inner_join(a1, a2, by = dplyr::join_by(mod, Celltype_Gabitto_2024, Direction, Consistency), relationship = "many-to-many") 
  return(out)
}, flist[c(2,4)]) |>
  filter(!duplicated(paste0(mod, Celltype_Gabitto_2024, Direction, Consistency))) |>
    mutate(Direction = factor(Direction, levels = c(1, -1))) |> 
    rename("Celltype" = "Celltype_Gabitto_2024") |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    select(mod, Celltype, Direction) |>
    group_by(Celltype, Direction, .drop = FALSE) |>
    summarise(num_sig = n(), .groups = "drop") |>
    left_join(class_info, by = join_by(Celltype == Subclass)) |>
    mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) |>
    arrange(Celltype) |>
    mutate(num_sig = case_match(num_sig, 
                                0 ~ NA,
                                .default = num_sig
                                ),
           comp = "Combined (AD modules)") 
flistall <- rbind(flistjoin, flistjoinros)

## Plot dotplot
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

dir1 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/Projection_analyses/related_analyses/rand_mod_analysis/rand_mod_dcopa_pvals/")
dot_list <- list(
  file.path(dir1, "SEAAD2024_MIT_shared_DFC", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_MTG", "dcopa_scorecard_summary_conVsAllAD.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_DFC.csv"),
  file.path(dir1, "SEAAD2024_MIT_shared_ROSMAPmods", "dcopa_scorecard_summary_conVsAllAD_MTG.csv")
) |>
  lapply(\(x){
    fread(x, data.table = F) |> 
      mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass))) 
  })
names(dot_list) <- c("DFC", "MTG", "DFC (AD modules)", "MTG (AD modules)")
typevec_hi <- c("DFC",
                "MTG",
                "Combined",
                "DFC (AD modules)",
                "MTG (AD modules)",
                "Combined (AD modules)")

typevec_lo <- c("DFC",
                "MTG",
                "Combined",
                "DFC (AD modules)",
                "MTG (AD modules)",
                "Combined (AD modules)")

df_hi <- lapply(seq_along(dot_list), \(x){
  dot_list[[x]] |>
    filter(grepl("Higher", type)) |>
    mutate(comp = names(dot_list)[x])
}) |> do.call(what = "rbind") |>
  select(!type)
df_hi <- rbind(df_hi, flistall |> filter(Direction == 1) |> select(!Direction)) |> 
  mutate(comp = factor(comp, levels = rev(typevec_hi)))

df_lo <- lapply(seq_along(dot_list), \(x){
  dot_list[[x]] |>
    filter(grepl("Lower", type)) |>
    mutate(comp = names(dot_list)[x])
}) |> do.call(what = "rbind") |>
  select(!type)
df_lo <- rbind(df_lo, flistall |> filter(Direction == -1) |> select(!Direction)) |> 
  mutate(comp = factor(comp, levels = rev(typevec_lo)))

cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

phi <- df_hi |>
    ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
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
      scale_size_continuous(breaks = scales::pretty_breaks(n = 3)) 
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_B.pdf"), height = 2.6, width = 5)

plo <- df_lo |>
    ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
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
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_A.pdf"), height = 2.6, width = 5)

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

# Load module data (topmodposbc - all dCoPA analyses done using topmodposbc definitions)
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
filter_under <- 3
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])

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

dat_vec <- c("gabitto_vs_liu_DFC_AllADVsCon",    
             "gabitto_vs_liu_MTG_AllADVsCon",         
             "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP") 

title_vec <- c("Con vs all AD (DFC)",
               "Con vs all AD (MTG)",
               "Con vs all AD (ROSMAP, DFC)",
               "Con vs all AD (ROSMAP, MTG)")

label_vec <- c("Shared in Gabitto and Liu DFC (con vs all AD)",    
               "Shared in Gabitto and Liu MTG (con vs all AD)",           
               "Shared in Gabitto and Liu DFC (con vs all AD, ROSMAP mods)",    
               "Shared in Gabitto and Liu MTG (con vs all AD, ROSMAP mods)") 

dat_list <- mapply(\(x, y){
  out <- fread(file.path(dir1, x, "shared_output_table.csv"), data.table = F)
  return(out)
}, dat_vec, title_vec, SIMPLIFY = F)
names(dat_list) <- label_vec

comp_names <- c("DFC", "MTG")

genedf <- lapply(seq_along(comp_names), \(i){
  d1 <- dat_list[[i]]
  d2 <- dat_list[[i + 2]]
  #common_cts <- unique(c(d1[,2], d2[,2]))
  common_cts <- class_info$Subclass

  outdf <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(common_cts)){
      m1 <- d1 |> filter(Celltype_Gabitto_2024 == common_cts[c],
                        Direction == d)
      g1 <- unique(unlist(mods[m1$mod]))
      m2 <- d2 |> filter(Celltype_Gabitto_2024 == common_cts[c],
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
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v1/panel_C1.pdf"), height = 4, width = 7)

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
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v1/panel_C2.pdf"), height = 4, width = 7)



############
# Upset plot
############

# DFC modules
d1 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon/shared_output_table.csv"))
d2 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon_ROSMAP/shared_output_table.csv"))

common_cts <- unique(c(d1[,2], d2[,2]))

# Collect lists of mods per subclass (DFC)
outlist <- lapply(c(-1, 1), \(d){
  sublist <- list()
  for(c in seq_along(common_cts)){
    m1 <- d1 |> filter(Celltype_Gabitto_2024 == common_cts[c],
                      Direction == d)
    g1 <- unique(unlist(mods[m1$mod]))
    m2 <- d2 |> filter(Celltype_Gabitto_2024 == common_cts[c],
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

# DFC (lower in disease)
p1 <- upset(fromList(outlist[[1]]), 
      sets = names(outlist[[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlaps", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.5, 0.5),
      show.numbers = F)

pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v1/panel_D1.pdf"), 
    width = 3.5, height = 2)
p1 
#grid.text("All AD vs Con, DFC, lower in severe", x = 0.65, y=0.98, gp=gpar(fontsize=8))
dev.off()

# DFC (higher in disease)
# Only astrocytes

# MTG modules
d3 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon/shared_output_table.csv"))
d4 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon_ROSMAP/shared_output_table.csv"))

common_cts <- unique(c(d3[,2], d4[,2]))

# Collect lists of mods per subclass (DFC)
outlist <- lapply(c(-1, 1), \(d){
  sublist <- list()
  for(c in seq_along(common_cts)){
    m1 <- d3 |> filter(Celltype_Gabitto_2024 == common_cts[c],
                      Direction == d)
    g1 <- unique(unlist(mods[m1$mod]))
    m2 <- d4 |> filter(Celltype_Gabitto_2024 == common_cts[c],
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

# MTG (lower in disease)
p1 <- upset(fromList(outlist[[1]]), 
      sets = names(outlist[[1]]), 
      order.by = "freq",
      mainbar.y.label = "# of overlaps", 
      sets.x.label = "# of genes",
      mb.ratio = c(0.4, 0.6),
      show.numbers = F)

pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v1/panel_D3.pdf"), 
    width = 4.7, height = 2.7)
p1 
#grid.text("All AD vs Con, MTG, lower in severe", x = 0.65, y=0.98, gp=gpar(fontsize=8))
dev.off()

# MTG (higher in disease)
# None
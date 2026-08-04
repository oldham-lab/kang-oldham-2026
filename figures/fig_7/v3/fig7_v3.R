library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v3"
if(!dir.exists(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)))
  dir.create(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version), recursive = T)


#########################################
# Upset plot for DE genes (pseudobulk DE)
# Gabitto vs Liu, DFC vs MTG, DESeq2 vs edgeR
# 8 comparisons total
###########################################

# see /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_7/v3/panel_A_upset_DE.R

###########
# Panel B
############
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
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]

dat_vec <- list(
  # con vs All DFC
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv"),
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_PFC/euclidean_distances/output_table_Subclass.csv"),
  # con vs All MTG
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv"),
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_MTC/euclidean_distances/output_table_Subclass.csv"),
  # con vs All DFC ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv"),
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_DFC/euclidean_distances/output_table_Subclass.csv"),
  # all AD vs con MTG ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv"),
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_MTG/euclidean_distances/output_table_Subclass.csv")
)

dat_vec_names <- c("Gabitto_AllADVsCon_DFC",
                   "Liu_AllADVsCon_DFC",
                   "Gabitto_AllADVsCon_MTG",
                   "Liu_AllADVsCon_MTG",
                   "Gabitto_AllADVsCon_DFC_ROSMAP",
                   "Liu_AllADVsCon_DFC_ROSMAP",
                   "Gabitto_AllADVsCon_MTG_ROSMAP",
                   "Liu_AllADVsCon_MTG_ROSMAP")

# Load all relevant dcopa output tables (for all comparisons)
dcopa_shared <- lapply(dat_vec, \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods) |>
    filter(sig_FDR)
})

# Function for extracting genes belonging to shared significant dCoPA mods (topmodposbc)
extract_genes <- function(d1){
  cts <- unique(d1[,2])

  outlist <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(cts)){
      m1 <- d1 |> filter(Celltype == cts[c],
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

# Save list off all genes by celltype 
gsea_input <- dcopa_to_df(dcopa_allct) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
         "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
         "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu")
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))
fwrite(gsea_input, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_B_dcopa_genelist.csv"))

# /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_7/v3/panel_B_upset_dcopa.R

###########
# Panel C-D (old panel A-B)
##########

# Find overlap between all AllADvsCon comparisons
# v3: remove "combined" row

# Load dCoPA output files (AllADvsCon DFC, MTG)
file1 <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon"),        
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon_ROSMAP"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon"),        
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_MTG_AllADVsCon_ROSMAP"))
flist <- lapply(file1, \(x) fread(data.table = F, file = file.path(x, "shared_output_table.csv")) |>
    filter(mod %in% these_mods)
)

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
             "DFC",
             "DFC (AD modules)"
             )

df_hi <- lapply(seq_along(dot_list), \(x){
  dot_list[[x]] |>
    filter(grepl("Higher", type)) |>
    mutate(comp = names(dot_list)[x])
}) |> do.call(what = "rbind") |>
  select(!type) |> 
  mutate(comp = factor(comp, levels = rev(typevec)))

df_lo <- lapply(seq_along(dot_list), \(x){
  dot_list[[x]] |>
    filter(grepl("Lower", type)) |>
    mutate(comp = names(dot_list)[x])
}) |> do.call(what = "rbind") |>
  select(!type) |> 
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
            #legend.position = "bottom",
            legend.position = "none", # separate legend created in inkscape
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
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_C.pdf"), height = 2, width = 5)
ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_C.svg"), height = 2, width = 5)

phi <- df_hi |>
    ggplot(aes(x = Celltype1, y = comp, size = num_sig, fill = Class)) +
      theme_minimal() + 
      geom_point(color = "black", pch = 21) +
      theme(text = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
            axis.text.y = element_text(size = 7), 
            legend.direction = "horizontal",    
            #legend.position = "bottom",
            legend.position = "none", # separate legend created in inkscape
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
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_D.pdf"), height = 2, width = 5)
ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "panel_D.svg"), height = 2, width = 5)


################
# Panel C + D alternate
###############

# Load class info + subclass labels for Liu
obj <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_MTC/sn_summary_tables/sn_summary_objects_log.qs"))
class_info_mit <- fread(data.table = F, file = file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv")) |> 
  select(RNA.Class, RNA.Subclass) |>
  filter(!duplicated(RNA.Subclass), 
          RNA.Subclass %in% colnames(obj[[1]][[1]][[1]])) |>
  mutate(Class = case_match(RNA.Class, "Exc" ~ "Glutamatergic", "Inh" ~ "GABAergic", .default = "Non-neuronal")) |>
  rename(Subclass = RNA.Subclass) |>
  arrange(factor(Class, levels = c("Glutamatergic", "GABAergic", "Non-neuronal"))) |>
  mutate(Subclass_final = Subclass)
rm(obj)

# Names for each dataset
dat_vec_names_short <- c("Gabitto DFC",
                   "Liu DFC",
                   "Gabitto MTG",
                   "Liu MTG",
                   "Gabitto DFC (AD modules)",
                   "Liu DFC (AD modules)",
                   "Gabitto MTG (AD modules)",
                   "Liu MTG (AD modules)")

# Function for processing dcopa output into plotting input
process_dcopa_output_tab <- function(sum_tab, dat_name, class_info_df){
  allsum <- sum_tab |> 
    mutate(Direction = factor(Direction, levels = c(-1, 1)),
           Celltype = factor(Celltype, levels = unique(class_info_df$Subclass))) |>
    select(mod, Celltype, Direction) |>
    group_by(Celltype, Direction, .drop = FALSE) |>
    summarise(num_sig = n(), .groups = "drop") |>
    left_join(class_info_df, by = join_by(Celltype == Subclass)) |>
    arrange(Celltype) |>
    mutate(num_sig = case_match(num_sig, 
                                0 ~ NA,
                                .default = num_sig
                                ),
           comp = dat_name) 
  return(allsum)
}

# Create inputs for class info
plotlist <- list()
class_info_list <- rep(list(class_info), 4)
class_info_list_mit <- rep(list(class_info_mit), 4)

# Process liu dcopa outputs
plotlist[[1]] <- mapply(process_dcopa_output_tab, dcopa_shared[c(2,4,6,8)], dat_vec_names_short[c(2,4,6,8)], class_info_list_mit, SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(Celltype = factor(Celltype, levels = unique(class_info_mit$Subclass)),
         comp = factor(comp, levels = dat_vec_names_short[c(6,2,8,4)]))

# Process gabitto dcopa outputs
plotlist[[2]] <- mapply(process_dcopa_output_tab, dcopa_shared[c(1,3,5,7)], dat_vec_names_short[c(1,3,5,7)], class_info_list, SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass)),
         comp = factor(comp, levels = dat_vec_names_short[c(5,1,7,3)]))

# Set size range for dot legend
limit_vec = c(1, max(c(plotlist[[1]]$num_sig, plotlist[[2]]$num_sig), na.rm = T))

# Plot (4 plots total)
dvec <- c("liu", "gabitto")
for(i in 1:2){
  plo <- plotlist[[i]] |>
      filter(Direction == -1) |>
      ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
        theme_minimal() + 
        geom_point(color = "black", pch = 21) +
        theme(text = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
              axis.text.y = element_text(size = 7), 
              legend.direction = "horizontal",    
              legend.position = "bottom",
              #legend.position = "none", # separate legend created in inkscape
              legend.box = "vertical",
              legend.spacing.y = unit(4, "mm"),
              legend.title = element_blank(),
              legend.margin = margin(-0.5, 0, 0, 0, "cm"),
              legend.text = element_text(margin = margin(0, 0, 0, -0.02, "cm")),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        labs(x = "", y = "") +
        scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
        scale_size_continuous(limits = limit_vec,
                              breaks = scales::pretty_breaks(n = 3)) 
  ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_C_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  #ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_C_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_C_indiv_", dvec[i], "_nolegend.pdf")), height = 2, width = 5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_C_indiv_", dvec[i], "_nolegend.svg")), height = 2, width = 5)

  phi <- plotlist[[i]] |>
      filter(Direction == 1) |>
      ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
        theme_minimal() + 
        geom_point(color = "black", pch = 21) +
        theme(text = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
              axis.text.y = element_text(size = 7), 
              legend.direction = "horizontal",    
              legend.position = "bottom",
              #legend.position = "none", # separate legend created in inkscape
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
        scale_size_continuous(limits = limit_vec, # Manually set the scale to be the same as lo object
                              breaks = scales::pretty_breaks(n = 2)) 
  ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_D_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  #ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_D_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_D_indiv_", dvec[i], "_nolegend.pdf")), height = 2, width = 5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_D_indiv_", dvec[i], "_nolegend.svg")), height = 2, width = 5)
}

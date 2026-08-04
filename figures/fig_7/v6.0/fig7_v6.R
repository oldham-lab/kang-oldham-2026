library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v6.0"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)



##########
# Load initial data
###########
CT_RENAME <- c(
  "Lamp5"     = "LAMP5",
  "Lamp5 Lhx6"= "LAMP5 LHX6",
  "Pax6"      = "PAX6",
  "Pvalb"     = "PVALB",
  "Sncg"      = "SNCG",
  "Sst"       = "SST",
  "Vip"       = "VIP",
  "Sst Chodl" = "SST CHODL"
)

# Load class info + subclass labels for Gabitto
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
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  # con vs All MTG
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"), 
  # con vs All DFC ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"), 
  # all AD vs con MTG ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"), 
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv")
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
    filter(mod %in% these_mods,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency))
})

# Function for extracting genes belonging to shared significant dCoPA mods (topmodposbc)
extract_genes <- function(d1){
  cts <- unique(d1$Celltype)

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

# Organize list of genes into dataframe 
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


################
# Panel A+B
###############

find_output_overlap <- function(dat1, dat2){
  dat_out <- dplyr::inner_join(dat1, dat2, by = dplyr::join_by(mod, Direction, Celltype, Consistency)) |>
    dplyr::arrange(mod) 
  return(dat_out)
}

# Function for finding overlap of overlap of outputs:
find_output_overlap2 <- function(dat1, dat2){
  dat_out <- dplyr::inner_join(dat1, dat2, by = dplyr::join_by(mod, Celltype, Direction, Consistency)) |>
      dplyr::filter(!duplicated(paste0(mod, Celltype, Direction))) |>
      dplyr::arrange(mod) |>
      dplyr::select(c(mod, Celltype, Direction, Consistency))
  return(dat_out)
}

# Calculate overlaps
dcopa_shared_overlaps <- list(
  find_output_overlap(dcopa_shared[[1]], dcopa_shared[[2]]), # Gabitto + Liu DFC
  find_output_overlap(dcopa_shared[[3]], dcopa_shared[[4]]), # Gabitto + Liu MTG
  find_output_overlap(dcopa_shared[[5]], dcopa_shared[[6]]), # Gabitto + Liu DFC (AD modules)
  find_output_overlap(dcopa_shared[[7]], dcopa_shared[[8]]) # Gabitto + Liu DFC (AD modules)
)

dcopa_shared_overlaps[[5]] <- find_output_overlap2(dcopa_shared_overlaps[[1]], dcopa_shared_overlaps[[2]]) # Gabitto/Liu DFC/MTG
dcopa_shared_overlaps[[6]] <- find_output_overlap2(dcopa_shared_overlaps[[3]], dcopa_shared_overlaps[[4]]) # Gabitto/Liu DFC/MTG AD modules 

# Names for each dataset (plus overlaps)
# dat_vec_names_short <- c("Gabitto DFC",
#                          "Liu DFC",
#                          "Gabitto MTG",
#                          "Liu MTG",
#                          "Gabitto DFC (AD modules)",
#                          "Liu DFC (AD modules)",
#                          "Gabitto MTG (AD modules)",
#                          "Liu MTG (AD modules)",
#                          "Gabitto + Liu DFC",
#                          "Gabitto + Liu MTG",
#                          "Gabitto + Liu DFC (AD modules)",
#                          "Gabitto + Liu MTG (AD modules)",
#                          "Gabitto + Liu, DFC + MTG",
#                          "Gabitto + Liu, DFC + MTG (AD modules)")

dat_vec_names_short <- c("CTRL modules | Gabitto SN",
                         "CTRL modules | Liu SN",
                         "CTRL modules | Gabitto SN",
                         "CTRL modules | Liu SN",
                         "AD modules | Gabitto SN",
                         "AD modules | Liu SN",
                         "AD modules | Gabitto SN",
                         "AD modules | Liu SN",
                         "CTRL modules | Gabitto + Liu SN",
                         "CTRL modules | Gabitto + Liu SN",
                         "AD modules | Gabitto + Liu SN",
                         "AD modules | Gabitto + Liu SN",
                         "CTRL modules | Gabitto + Liu SN",
                         "AD modules | Gabitto + Liu SN")                         

# Function for processing dcopa output into plotting input
process_dcopa_output_tab <- function(sum_tab, dat_name){
  class_info_df <- class_info
  allsum <- sum_tab |> 
  mutate(Direction = factor(Direction, levels = c(-1, 1)),
          Celltype = factor(Celltype, levels = unique(class_info_df$Subclass))) |>
  dplyr::select(mod, Celltype, Direction) |>
  group_by(Celltype, Direction, .drop = FALSE) |>
  summarise(num_sig = n(), .groups = "drop") |>
  left_join(class_info_df |> 
  dplyr::select(Subclass, Class), by = join_by(Celltype == Subclass)) |>
  arrange(Celltype) |>
  mutate(num_sig = case_match(num_sig, 
                              0 ~ NA,
                              .default = num_sig
                              ),
          comp = dat_name) 
  return(allsum)
}

# # Process bulk megaset outputs
plotlist <- list()
# MTG
plotlist[[1]] <- mapply(process_dcopa_output_tab, 
                        c(dcopa_shared[c(3:4, 7:8)], dcopa_shared_overlaps[c(2, 4)]), 
                        dat_vec_names_short[c(3:4, 7:8, 10, 12)], 
                        SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass)),
         comp = factor(comp, levels = rev(dat_vec_names_short[c(3,7,4,8,10,12)])))
# DFC
plotlist[[2]] <- mapply(process_dcopa_output_tab, 
                        c(dcopa_shared[c(1:2, 5:6)], dcopa_shared_overlaps[c(1, 3)]), 
                        dat_vec_names_short[c(1:2, 5:6, 9, 11)], 
                        SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass)),
         comp = factor(comp, levels = rev(dat_vec_names_short[c(1,5,2,6,9,11)])))

# Capitalize celltype names
plotlist <- plotlist |> lapply(\(x) x |> mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME)))

# Set size range for dot legend
limit_vec = c(1, max(c(plotlist[[1]]$num_sig, plotlist[[2]]$num_sig), na.rm = T))

# Plot (4 plots total)
dvec <- c("MTG", "DFC") # save suffix
cc_colors <- RColorBrewer::brewer.pal(3, "Set1")
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
  ggsave(plo, file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  #ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_C_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)

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
  ggsave(phi, file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  #ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_D_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)
}

# Save legend as standalone SVG
legend_grob <- cowplot::get_legend(
  plo + theme(
    legend.text       = element_text(size = 5),
    legend.title      = element_text(size = 5),
    legend.key.size   = unit(3, "mm"),
    legend.justification = "center",
    legend.box.just   = "center"
  )
)
p_legend_out <- ggpubr::as_ggplot(legend_grob) +
  theme(plot.margin = margin(0, 0, 0, 0))
ggsave(p_legend_out, file = file.path(save_dir, "panel_A_B_legend.svg"),
        height = 0.8, width = 3, bg = "transparent")


###########
# Panel C-D?
############

# /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_7/v3/panel_B_upset_dcopa.R


###########
# Panel E-F
###########

# Save genelist for N_GL_MTG/DGC and AD_GL_MTG/DFC
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

# Gabitto + Liu MTG (Bulk mods vs AD mods)
for(ct in c("MTG", "DFC")){
  if(ct == "MTG"){
    dcopa_allct <- lapply(dcopa_shared_overlaps[c(2,4)], extract_genes) |> setNames(dat_vec_names_short[c(10, 12)])
  } else if (ct == "DFC"){
    dcopa_allct <- lapply(dcopa_shared_overlaps[c(1,3)], extract_genes) |> setNames(dat_vec_names_short[c(9, 11)])
  }

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
  fwrite(gsea_input, file = file.path(save_dir, paste0("/panel_E_", ct,  "_dcopa_genelist.csv")))
}


#######
# Panel G+
#######

# Load genelists
dfc_list <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/panel_E_DFC_dcopa_genelist.csv"))
mtg_list <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/panel_E_MTG_dcopa_genelist.csv"))

find_comps <- function(dcopa){                                                          
  comps <- unique(dcopa$Comparison)  # exactly 2                  
                                                                  
  genes_a <- dcopa |> filter(Comparison == comps[1]) |>           
  distinct(Celltype, Gene)                                        
  genes_b <- dcopa |> filter(Comparison == comps[2]) |>           
  distinct(Celltype, Gene)                                        
                                                                  
  overlap_table <- inner_join(genes_a, genes_b, by = c("Celltype",
   "Gene")) |>                                                    
    arrange(Celltype, Gene)     
  return(overlap_table)
}

##########
# Start with DFC
g_overlaps <- find_comps(dfc_list)
fwrite(g_overlaps, file = file.path(save_dir, "dfc_overlaps.csv"))

# Run GSEA
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
b_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = T)
g_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = F)
gsea_out <- rbind(b_out, g_out)

# Significant genesets per celltype (p < GSEA_PTHRESH)
GSEA_PTHRESH <- 0.05
gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
  pvals <- gsea_out[[3 + i]]
  gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
}), celltypes)

# ── AD gene database: genes in g_overlaps with documented AD associations ────
# Each entry: gene = list(n_refs = approx. no. relevant papers, ref = bibliography no.)
# Ordered by n_refs descending within each celltype when displayed.
# Bibliography: gsea_summary_bibliography.md
library(gt)

ad_db <- list(
  # Retromer / endosomal sorting
  VPS35    = list(n=8,  ref=1),
  VPS26A   = list(n=6,  ref=1),
  SNX2     = list(n=5,  ref=1),
  APPBP2   = list(n=5,  ref=2),
  APPL1    = list(n=5,  ref=3),
  RANBP9   = list(n=5,  ref=4),
  RAB11FIP2= list(n=3,  ref=1),
  RAB3GAP1 = list(n=3,  ref=1),
  # Autophagy / lysosomal clearance
  PIK3C3   = list(n=10, ref=5),
  ATG2B    = list(n=4,  ref=5),
  SPG11    = list(n=4,  ref=5),
  TRAPPC11 = list(n=4,  ref=5),
  TRAPPC8  = list(n=3,  ref=5),
  TBC1D15  = list(n=3,  ref=5),
  TAX1BP1  = list(n=3,  ref=5),
  # Mitochondrial function / dynamics
  DNM1L    = list(n=8,  ref=6),
  IDH3A    = list(n=4,  ref=7),
  LRPPRC   = list(n=3,  ref=7),
  TOMM70   = list(n=3,  ref=7),
  AFG3L2   = list(n=3,  ref=7),
  # Ubiquitin-proteasome system
  PPP2R2A  = list(n=10, ref=8),
  USP14    = list(n=6,  ref=9),
  UBQLN1   = list(n=6,  ref=10),
  CUL3     = list(n=4,  ref=9),
  CUL1     = list(n=3,  ref=9),
  CUL2     = list(n=3,  ref=9),
  CUL4B    = list(n=3,  ref=9),
  PSMD5    = list(n=3,  ref=9),
  UBE2K    = list(n=3,  ref=9),
  USP15    = list(n=3,  ref=9),
  FBXW11   = list(n=2,  ref=9),
  FBXO11   = list(n=2,  ref=9),
  KLHL7    = list(n=2,  ref=9),
  WWP1     = list(n=2,  ref=9),
  # Tau phosphorylation / modification
  PPM1A    = list(n=5,  ref=8),
  SENP1    = list(n=4,  ref=11),
  MBNL1    = list(n=3,  ref=12),
  MBNL2    = list(n=3,  ref=12),
  # ER stress / ERAD
  EIF2AK2  = list(n=6,  ref=13),
  EDEM3    = list(n=3,  ref=13),
  ERLEC1   = list(n=3,  ref=13),
  ERO1A    = list(n=3,  ref=13),
  # Neuroinflammation
  NFKB1    = list(n=12, ref=14),
  NAPEPLD  = list(n=3,  ref=14),
  CACNA2D1 = list(n=3,  ref=15),
  # Axonal / vesicular transport
  DYNC1LI1 = list(n=6,  ref=16),
  DCTN4    = list(n=4,  ref=16),
  STX7     = list(n=3,  ref=1),
  # Epigenetics
  HDAC2    = list(n=10, ref=17),
  # Apoptosis
  APAF1    = list(n=5,  ref=18),
  # Lipid metabolism
  CERT1    = list(n=3,  ref=19),
  SC5D     = list(n=3,  ref=19),
  DDHD2    = list(n=2,  ref=19),
  SACM1L   = list(n=2,  ref=19),
  # RNA regulation / stress granules
  CAPRIN1  = list(n=3,  ref=20),
  DDX6     = list(n=2,  ref=20),
  YTHDF3   = list(n=2,  ref=20),
  # APP processing / trafficking
  LAMTOR3  = list(n=3,  ref=5),
  SEC23A   = list(n=2,  ref=2),
  SEC24A   = list(n=2,  ref=2),
  WASL     = list(n=2,  ref=3)
)

# ── Build table ─────────────────────────────────────────────────
celltypes  <- unique(g_overlaps$Celltype)

# Significant genesets per celltype (BH-adjusted p < GSEA_FDR_THRESH)
GSEA_FDR_THRESH <- 0.05
gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
  pvals <- gsea_out[[3 + i]]
  fdr   <- p.adjust(pvals, method = "BH")
  gsea_out$SetName[!is.na(fdr) & fdr < GSEA_FDR_THRESH]
}), celltypes)
qsave(gsea_sig, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/gsea_sig.qs"))

CT_RENAME <- c(
  "Lamp5"      = "LAMP5",
  "Lamp5 Lhx6" = "LAMP5 LHX6",
  "Pax6"       = "PAX6",
  "Pvalb"      = "PVALB",
  "Sncg"       = "SNCG",
  "Sst"        = "SST",
  "Vip"        = "VIP",
  "Sst Chodl"  = "SST CHODL"
)

build_row <- function(ct) {
  genes  <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  # Top 3 GO gene sets for this celltype (columns assumed in same order as celltypes)
  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    top_idx  <- order(go_pvals)[1:min(3, sum(!is.na(go_pvals)))]
    top_sets <- gsea_out$SetName[top_idx]
    top_pval <- pvals[top_idx]
    top_str  <- paste0(
      seq_along(top_sets), ". ", top_sets,
      " (p=", formatC(top_pval, format = "e", digits = 2), ")",
      collapse = "<br>"
    )
  }

  # AD genes: match uppercase gene names against ad_db, sort by n_refs descending
  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db)]

  if (length(hits) > 0) {
    n_refs   <- sapply(hits_up, function(g) ad_db[[g]]$n)
    ref_nums <- sapply(hits_up, function(g) ad_db[[g]]$ref)
    ord      <- order(n_refs, decreasing = TRUE)
    hits_s   <- hits[ord]
    refs_s   <- ref_nums[ord]
    ad_str   <- paste0(hits_s, "<sup>", refs_s, "</sup>", collapse = ", ")
  } else {
    ad_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    AD_genes     = ad_str,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(celltypes, build_row)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

# ── Format with gt ─────────────────────────────────────────────────────────────
GT_SCALE <- 0.60  # uniform scale factor for PDF table size (1 = original)

gt_table <- summary_df |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by p-value)",
    AD_genes     = "AD-associated genes"
  ) |>
  fmt_markdown(columns = AD_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    AD_genes     ~ px(262 * GT_SCALE)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = list(cells_body(columns = N_genes),
                     cells_column_labels(columns = N_genes))
  ) |>
  tab_options(
    table.font.names      = c("Arial", "Helvetica", "sans-serif"),
    table.font.size       = 11 * GT_SCALE,
    data_row.padding      = px(6 * GT_SCALE),
    column_labels.padding = px(8 * GT_SCALE)
  )

out_path <- file.path(save_dir, "panel_G_gsea_summary_table.html")
gtsave(gt_table, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table.pdf")
tryCatch({
  gtsave(gt_table, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})

bib_md   <- here::here("Code_for_figures/fig_7/v6/gsea_summary_bibliography.md")
bib_docx <- here::here("Code_for_figures/fig_7/v6/gsea_summary_bibliography_dfc.docx")
tryCatch({
  rmarkdown::pandoc_convert(bib_md, to = "docx", output = bib_docx)
  message("Saved: ", bib_docx)
}, error = function(e) {
  message("Bibliography docx conversion failed (requires pandoc): ", conditionMessage(e))
})


##########dd
# do mtg
g_overlaps <- find_comps(mtg_list)
fwrite(g_overlaps, file = file.path(save_dir, "mtg_overlaps.csv"))

# Run GSEA
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
b_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = T)
g_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = F)
gsea_out <- rbind(b_out, g_out)

# Significant genesets per celltype (p < GSEA_PTHRESH)
GSEA_PTHRESH <- 0.05
gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
  pvals <- gsea_out[[3 + i]]
  gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
}), celltypes)

# ── AD gene database: genes in g_overlaps with documented AD associations ────
# Each entry: gene = list(n_refs = approx. no. relevant papers, ref = bibliography no.)
# Ordered by n_refs descending within each celltype when displayed.
# Bibliography: gsea_summary_bibliography_mtg.md
library(gt)

ad_db <- list(
  # Retromer / endosomal sorting (confirmed in MTG)
  SORL1    = list(n=12, ref=21),
  VPS35    = list(n=8,  ref=1),
  VPS26B   = list(n=8,  ref=29),
  VPS29    = list(n=7,  ref=29),
  APPBP2   = list(n=5,  ref=2),
  APPL1    = list(n=5,  ref=3),
  RANBP9   = list(n=5,  ref=4),
  APBA2    = list(n=6,  ref=28),
  RAB3GAP1 = list(n=3,  ref=1),
  SEC23A   = list(n=2,  ref=2),
  SEC24A   = list(n=2,  ref=2),
  # Autophagy / lysosomal clearance
  PIK3C3   = list(n=10, ref=5),
  SYNJ1    = list(n=9,  ref=25),
  ATG16L1  = list(n=9,  ref=5),
  MAP1LC3B = list(n=8,  ref=5),
  ATP6V1A  = list(n=8,  ref=24),
  RAB7A    = list(n=8,  ref=30),
  ATG2B    = list(n=4,  ref=5),
  TRAPPC11 = list(n=4,  ref=5),
  TRAPPC8  = list(n=3,  ref=5),
  TBC1D15  = list(n=3,  ref=5),
  LAMTOR3  = list(n=3,  ref=5),
  # Mitochondrial function / dynamics
  DNM1L    = list(n=8,  ref=6),
  IDH3A    = list(n=4,  ref=7),
  LRPPRC   = list(n=3,  ref=7),
  TOMM70   = list(n=3,  ref=7),
  AFG3L2   = list(n=3,  ref=7),
  # Tau kinase / phosphatase signalling
  CDK5R1   = list(n=10, ref=22),
  CAMKK2   = list(n=9,  ref=26),
  PPP2CA   = list(n=8,  ref=8),
  PPP2R2B  = list(n=6,  ref=8),
  MAPK1    = list(n=7,  ref=27),
  # Ubiquitin-proteasome system
  UBQLN1   = list(n=6,  ref=10),
  UBQLN2   = list(n=6,  ref=10),
  PSMD7    = list(n=6,  ref=9),
  PSMD1    = list(n=5,  ref=9),
  CUL3     = list(n=4,  ref=9),
  CUL1     = list(n=3,  ref=9),
  CUL2     = list(n=3,  ref=9),
  USP15    = list(n=3,  ref=9),
  FBXW11   = list(n=2,  ref=9),
  FBXO11   = list(n=2,  ref=9),
  # ER stress / ERAD
  EDEM3    = list(n=3,  ref=13),
  ERLEC1   = list(n=3,  ref=13),
  # Neurodegeneration / proteostasis
  TMEM106B = list(n=9,  ref=23),
  DNAJC13  = list(n=6,  ref=36),
  TOLLIP   = list(n=5,  ref=37),
  # Axonal / vesicular transport
  DYNC1LI1 = list(n=6,  ref=16),
  DCTN4    = list(n=4,  ref=16),
  CLTC     = list(n=7,  ref=31),
  NSF      = list(n=7,  ref=34),
  RAB11A   = list(n=8,  ref=35),
  # ER-mitochondria contacts / lipid transfer
  VAPB     = list(n=7,  ref=32),
  DDHD2    = list(n=2,  ref=19),
  # PI3K-Akt-mTOR signalling
  AKT3     = list(n=7,  ref=33),
  # RNA regulation / stress granules
  CAPRIN1  = list(n=3,  ref=20)
)

# ── Build table ─────────────────────────────────────────────────
celltypes  <- unique(g_overlaps$Celltype)

# Significant genesets per celltype (BH-adjusted p < GSEA_FDR_THRESH)
GSEA_FDR_THRESH <- 0.05
gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
  pvals <- gsea_out[[3 + i]]
  fdr   <- p.adjust(pvals, method = "BH")
  gsea_out$SetName[!is.na(fdr) & fdr < GSEA_FDR_THRESH]
}), celltypes)
qsave(gsea_sig, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/gsea_sig_mtg.qs"))

build_row <- function(ct) {
  genes  <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  # Top 3 GO gene sets for this celltype (columns assumed in same order as celltypes)
  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    top_idx  <- order(go_pvals)[1:min(3, sum(!is.na(go_pvals)))]
    top_sets <- gsea_out$SetName[top_idx]
    top_pval <- pvals[top_idx]
    top_str  <- paste0(
      seq_along(top_sets), ". ", top_sets,
      " (p=", formatC(top_pval, format = "e", digits = 2), ")",
      collapse = "<br>"
    )
  }

  # AD genes: match uppercase gene names against ad_db, sort by n_refs descending
  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db)]

  if (length(hits) > 0) {
    n_refs   <- sapply(hits_up, function(g) ad_db[[g]]$n)
    ref_nums <- sapply(hits_up, function(g) ad_db[[g]]$ref)
    ord      <- order(n_refs, decreasing = TRUE)
    hits_s   <- hits[ord]
    refs_s   <- ref_nums[ord]
    ad_str   <- paste0(hits_s, "<sup>", refs_s, "</sup>", collapse = ", ")
  } else {
    ad_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    AD_genes     = ad_str,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(celltypes, build_row)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

# ── Format with gt ─────────────────────────────────────────────────────────────
GT_SCALE <- 0.60  # uniform scale factor for PDF table size (1 = original)

gt_table <- summary_df |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by p-value)",
    AD_genes     = "AD-associated genes"
  ) |>
  fmt_markdown(columns = AD_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    AD_genes     ~ px(262 * GT_SCALE)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = list(cells_body(columns = N_genes),
                     cells_column_labels(columns = N_genes))
  ) |>
  tab_options(
    table.font.names      = c("Arial", "Helvetica", "sans-serif"),
    table.font.size       = 11 * GT_SCALE,
    data_row.padding      = px(6 * GT_SCALE),
    column_labels.padding = px(8 * GT_SCALE)
  )

out_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg.html")
gtsave(gt_table, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg.pdf")
tryCatch({
  gtsave(gt_table, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})

bib_md   <- here::here("Code_for_figures/fig_7/v6/gsea_summary_bibliography_mtg.md")
bib_docx <- here::here("Code_for_figures/fig_7/v6/gsea_summary_bibliography_mtg.docx")
tryCatch({
  rmarkdown::pandoc_convert(bib_md, to = "docx", output = bib_docx)
  message("Saved: ", bib_docx)
}, error = function(e) {
  message("Bibliography docx conversion failed (requires pandoc): ", conditionMessage(e))
})


#  Prompt: Build an AD-associated gene database from a            
#   celltype–gene overlap table                                    
                  
#   I have a CSV file at [PATH] with two columns: Celltype and     
#   Gene.                                                          
#   It contains the overlap genes between two independent analyses
#   for a specific                                                 
#   brain region/comparison.
                                                                 
#   Please do the following:                                       
   
#   1. Read the CSV and extract all unique genes across all        
#   celltypes.      
                                                                 
#   2. For each gene, determine whether it has documented          
#   associations with
#      Alzheimer's disease. Consider any of the following types of 
#   evidence:                                                      
#      - GWAS hits or genetic risk variants
#      - Protein interactions with core AD proteins (APP, BACE1,   
#   PSEN1/2, MAPT, APOE)                                           
#      - Functional studies in AD brain tissue or AD mouse models  
#      - Involvement in AD-relevant pathways (amyloid processing,  
#   tau phosphorylation,                                           
#        autophagy/lysosomal clearance, neuroinflammation, synaptic
#    dysfunction,                                                  
#        mitochondrial dysfunction, UPS impairment, axonal
#   transport, ER stress)                                          
#      - Co-expression or proteomic studies in AD patient brain
                                                                 
#   3. For each AD-associated gene, assign:                        
#      - n: approximate number of relevant publications (use your  
#   knowledge;                                                     
#        higher = stronger/more replicated evidence)
#      - ref: a bibliography number pointing to the single most    
#   informative paper                                              
#        (landmark study, first mechanistic link, or highest-impact
#    GWAS paper)                                                   
                  
#   4. Rank genes by n (descending) and group them into biological 
#   pathway categories.
                                                                 
#   5. Output the result as an R list in this format:              
   
#      ad_db <- list(                                              
#        # [Pathway category]
#        GENE1 = list(n=12, ref=1),                                
#        GENE2 = list(n=8,  ref=2),                                
#        ...
#      )                                                           
                  
#   6. Also produce a companion bibliography Markdown file with one
#    entry per
#      reference number, in this format:                           
                                                                 
#      **[N]** Author et al. (Year). Title. *Journal*.             
#   Volume(Issue):pages. PMID: XXXXXXX                             
#      *(GENE1, GENE2 — one-line summary of the AD-relevant        
#   finding)*                                                      
   
#      Only include genes that were actually confirmed present in  
#   the input CSV.  
#      Do not invent genes or references.                          
                                                                 
#   ---
#   Key things to fill in: [PATH] for the input file. Optionally   
#   add context like the brain region, dataset name, or which      
#   pathway categories matter most for your analysis.

# ● Roughly correlated in rank order, but with systematic          
#   divergences worth knowing about:
                                                                 
#   Where they track well:                                         
#   - High-n genes like SORL1, PIK3C3, CDK5R1 would also have many 
#   PubMed hits for "[gene] AND Alzheimer's" — these are genuinely 
#   well-studied in AD context                                    
#   - The top ~10 genes by n are likely to be in the top ~15 by    
#   PubMed count                                               
                                                                 
#   Where they diverge:
                                                                 
#   Situation: Gene in large AD GWAS papers                        
#   My n: low-moderate                                             
#   PubMed count: high                                             
#   Reason: GWAS papers mention hundreds of genes; raw hits inflate
                                                                 
#     PubMed count                                                 
#   ────────────────────────────────────────
#   Situation: Gene studied deeply in small mechanistic literature 
#   My n: high      
#   PubMed count: moderate                                         
#   Reason: Few focused papers but each is directly relevant
#   ────────────────────────────────────────
#   Situation: Housekeeping gene (e.g. CUL1, CUL2)                 
#   My n: low (2–3)
#   PubMed count: potentially high                                 
#   Reason: Many mentions in broad proteasome/AD reviews
#   ────────────────────────────────────────
#   Situation: Recently implicated genes (TMEM106B, SYNJ1)         
#   My n: moderate
#   PubMed count: moderate-high                                    
#   Reason: Rapidly growing but recent literature

#   The key difference: My n was meant to capture mechanistically  
#   relevant papers — studies where the gene plays a direct role in
#    AD pathology. A raw PubMed search for GENE AND Alzheimer's    
#   picks up any mention, including large -omics datasets where the
#    gene appears in a supplementary table.

#   Practical implication: If you wanted to validate the rankings, 
#   a better PubMed proxy would be something like
#   GENE[Title/Abstract] AND Alzheimer[Title/Abstract]             
#   (title/abstract only, not full text), which filters out
#   supplementary-table mentions. That would likely correlate more
#   closely with my n values.

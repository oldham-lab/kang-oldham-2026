# v7
# - layout finalization
# - major error fix: used proper module list for AD modules (was previously using CTRL modules for both CTRL and AD dCoPA outputs)

library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v1/slide5/"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7_sup/"), version)
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

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "SST CHODL", "SST", "Chandelier", 
              "PVALB", "LAMP5", "LAMP5 LHX6", "PAX6",  "SNCG", 
              "VIP",  "Astrocyte", "Oligodendrocyte", "OPC", "Microglia-PVM", "Endothelial", "VLMC")

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

# Load CTRL module data (topmodposbc - all dCoPA analyses done using topmodposbc definitions)
filter_under <- 3

module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]

# Load AD module data
datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods_case <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods_case,length))
these_mods_case <- as.numeric(names(mods_case)[which(modulelengths>filter_under)])
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_rosmapAD_1127.csv"))
these_mods_case <- these_mods_case[!these_mods_case %in% which(sigcount_bonf$vals < 1)]

dat_vec <- list(
  # con vs All DFC
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/earlyAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/earlyAD_vs_Con_bulk_megaset_output_table.csv"),
  # con vs All MTG
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/earlyAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/earlyAD_vs_Con_bulk_megaset_output_table.csv"), 
  # con vs All DFC ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/earlyAD_vs_Con_rosmap_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/earlyAD_vs_Con_rosmap_output_table.csv"), 
  # all AD vs con MTG ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/earlyAD_vs_Con_rosmap_output_table.csv"), 
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/earlyAD_vs_Con_rosmap_output_table.csv")
)

dat_vec_names <- c("Gabitto_earlyADvsCTRL_DFC",
                   "Liu_earlyADvsCTRL_DFC",
                   "Gabitto_earlyADvsCTRL_MTG",
                   "Liu_earlyADvsCTRL_MTG",
                   "Gabitto_earlyADvsCTRL_DFC_ROSMAP",
                   "Liu_earlyADvsCTRL_DFC_ROSMAP",
                   "Gabitto_earlyADvsCTRL_MTG_ROSMAP",
                   "Liu_earlyADvsCTRL_MTG_ROSMAP")

# Load all relevant dcopa output tables (for all comparisons)
dcopa_sharedctrl <- lapply(dat_vec[1:4], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency))
})

dcopa_sharedcase <- lapply(dat_vec[5:8], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods_case,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency))
})

dcopa_shared <- c(dcopa_sharedctrl, dcopa_sharedcase)

# Function for extracting genes belonging to shared significant dCoPA mods (topmodposbc)
extract_genes <- function(d1, mods){
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
dcopa_allct_ctrl <- lapply(dcopa_sharedctrl, \(x) extract_genes(x, mods)) |> setNames(dat_vec_names[1:4])
dcopa_allct_case <- lapply(dcopa_sharedcase, \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names[5:8])

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
gsea_input_ctrl <- dcopa_to_df(dcopa_allct_ctrl) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
         "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
         "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu")
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))

gsea_input_case <- dcopa_to_df(dcopa_allct_case) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
         "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
         "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu")
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))

gsea_input <- rbind(gsea_input_ctrl, gsea_input_case)

fwrite(gsea_input, file = file.path(save_dir, "/panel_B_dcopa_genelist.csv"))


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
plotlist <- plotlist |> lapply(\(x) x |> 
  mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME),
         Celltype = factor(Celltype, levels = ct_order),
         Class = factor(Class, levels = c("Glutamatergic", "GABAergic", "Non-neuronal"))))

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
  plo 
)
p_legend_out <- ggpubr::as_ggplot(legend_grob) +
  theme(plot.margin = margin(0, 0, 0, 0))
ggsave(p_legend_out, file = file.path(save_dir, "panel_A_B_legend.svg"),
        height = 1.4, width = 3, bg = "transparent")


###########
# Panel C-D
###########

extract_genes <- function(d1, mods){
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
    dcopa_allct1 <- lapply(dcopa_shared_overlaps[c(2)], \(x) extract_genes(x, mods)) |> setNames(dat_vec_names_short[c(10)])
    dcopa_allct2 <- lapply(dcopa_shared_overlaps[c(4)], \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names_short[c(12)])
    dcopa_allct <- c(dcopa_allct1, dcopa_allct2)
  } else if (ct == "DFC"){
    dcopa_allct1 <- lapply(dcopa_shared_overlaps[c(1)], \(x) extract_genes(x, mods)) |> setNames(dat_vec_names_short[c(9)])
    dcopa_allct2 <- lapply(dcopa_shared_overlaps[c(3)], \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names_short[c(11)])
    dcopa_allct <- c(dcopa_allct1, dcopa_allct2)
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
  fwrite(gsea_input, file = file.path(save_dir, paste0("/panel_C_D_", ct,  "_dcopa_genelist.csv")))
}


###########
# Panel E-F
###########
extract_gene_overlaps <- function(d1, d2, dir, mods1, mods2){
  common_cts <- unique(c(d1[,2], d2[,2]))

  outlist <- lapply(dir, \(d){
    sublist <- list()
    for(c in seq_along(common_cts)){
      m1 <- d1 |> filter(Celltype == common_cts[c],
                        Direction == d)
      g1 <- unique(unlist(mods1[m1$mod]))
      m2 <- d2 |> filter(Celltype == common_cts[c],
                        Direction == d)
      g2 <- unique(unlist(mods2[m2$mod]))
      
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
  names(outlist) <- dir
  return(outlist)
}

dir <- -1
for(reg in c("MTG", "DFC")){
  if(reg == "MTG"){
    dcopa_allct <- extract_gene_overlaps(dcopa_shared_overlaps[[c(2)]], dcopa_shared_overlaps[[c(4)]], dir, mods, mods_case)[[1]] 
  } else if (reg == "DFC"){
    dcopa_allct <- extract_gene_overlaps(dcopa_shared_overlaps[[c(1)]], dcopa_shared_overlaps[[c(3)]], dir, mods, mods_case)[[1]] 
  }

  names(dcopa_allct) <- dplyr::recode(names(dcopa_allct), !!!CT_RENAME)

  # svg(file.path(save_dir, paste0("/panel_E_F_", reg, ".svg")), width = 7, height = 4)           
  #   upset(fromList(dcopa_allct), 
  #         sets = names(dcopa_allct), 
  #         order.by = "freq",          
  #         mainbar.y.label = "# of overlaps",
  #         sets.x.label = "# of genes",                    
  #         mb.ratio = c(0.5, 0.5), 
  #         show.numbers = F)                                   
  # dev.off()

  m <- make_comb_mat(dcopa_allct)
  m  <- m[comb_size(m) > 2]     # drop combinations with ≤1 overlapping gene
  cs <- comb_size(m)
  ss <- set_size(m)

  # Identify diagonal cells: single-set combinations, in display order
  comb_ord <- order(cs, decreasing = TRUE)
  set_ord  <- match(names(dcopa_allct), set_name(m))
  n_sets   <- length(set_name(m))
  n_combs  <- length(comb_name(m))
  is_diag  <- matrix(FALSE, nrow = n_sets, ncol = n_combs)
  for (j_disp in seq_len(n_combs)) {
    bits <- as.integer(strsplit(comb_name(m)[comb_ord[j_disp]], "")[[1]])
    if (sum(bits) == 1) {
      i_disp <- which(set_ord == which(bits == 1))
      if (length(i_disp) == 1) is_diag[i_disp, j_disp] <- TRUE
    }
  }

  p <- UpSet(
    m,            
    # sets = names(dcopa_allct) equivalent        
    set_order  = names(dcopa_allct),              
    # order.by = "freq" equivalent                
    comb_order = order(cs, decreasing = TRUE),    
    # show.numbers = FALSE equivalent (no count labels on bars)                                 
    # mb.ratio = c(0.5, 0.5) equivalent: control bar height vs matrix                            
    top_annotation = HeatmapAnnotation(
      "# of overlaps" = anno_barplot(             
        cs,       
        border  = FALSE,                          
        gp      = gpar(fill = "black"),
        height  = unit(3, "cm")   # increase/decrease to adjust ratio               
      ),                                          
      annotation_name_side = "left",              
      annotation_name_rot  = 90,
      annotation_name_offset = unit(10, "mm") 
    ),
    # sets.x.label = "# of genes" equivalent
    right_annotation = rowAnnotation(             
      "# of genes" = anno_barplot(                
        ss,                                       
        border    = FALSE,                        
        gp        = gpar(fill = "black"),
        width     = unit(2, "cm"),
        xlim      = c(max(ss), 0)
      )
    )                                             
  )               

  svglite::svglite(file.path(save_dir,                     
  paste0("/panel_E_F_", reg, ".svg")), width = 4, height = 2.5)                                     
  draw(p)         
  dev.off()
}


#######
# Summary tables (panels E-F)
#######

# Load genelists
dfc_list <- fread(data.table = F, file = file.path(save_dir, "panel_C_D_DFC_dcopa_genelist.csv"))
mtg_list <- fread(data.table = F, file = file.path(save_dir, "/panel_C_D_MTG_dcopa_genelist.csv"))

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

g_overlaps <- find_comps(dfc_list)
fwrite(g_overlaps, file = file.path(save_dir, "dfc_overlaps.csv"))
g_overlaps <- find_comps(mtg_list)
fwrite(g_overlaps, file = file.path(save_dir, "mtg_overlaps.csv"))

####### Start with DFC

g_overlaps <- find_comps(dfc_list)

# Run GSEA
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
b_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = T)
g_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = F)
gsea_out <- rbind(b_out, g_out)

# Significant genesets per celltype (p < GSEA_PTHRESH)
# GSEA_PTHRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
# }), celltypes)

# ── AD gene database: genes in g_overlaps with documented AD associations ────
# Each entry: gene = list(n_refs = approx. no. relevant papers, ref = bibliography no.)
# Ordered alphabetically within each celltype when displayed.
# Bibliography: gsea_summary_bibliography_dfc.md
# NOTE: All genes in dfc_overlaps are ribosomal proteins (RPL*/RPS*, excluded by criteria)
#       or lack any verified AD association (MCMDC2, WDHD1, ZRANB2). No entries qualify.
library(gt)

ad_db <- list()

# ── Build table ─────────────────────────────────────────────────
celltypes  <- unique(g_overlaps$Celltype)

# # Significant genesets per celltype (BH-adjusted p < GSEA_FDR_THRESH)
# GSEA_FDR_THRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   fdr   <- p.adjust(pvals, method = "BH")
#   gsea_out$SetName[!is.na(fdr) & fdr < GSEA_FDR_THRESH]
# }), celltypes)
# qsave(gsea_sig, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/gsea_sig.qs"))

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

  # AD genes: match uppercase gene names against ad_db, sort alphabetically
  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db)]

  if (length(hits) > 0) {
    ref_nums <- sapply(hits_up, function(g) ad_db[[g]]$ref)
    ord      <- order(hits_up)
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

out_path <- file.path(save_dir, "panel_G_gsea_summary_table_dfc.html")
gtsave(gt_table, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_dfc.pdf")
tryCatch({
  gtsave(gt_table, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})

bib_md   <- file.path(save_dir, "gsea_summary_bibliography_dfc.md")
bib_docx <- file.path(save_dir, "gsea_summary_bibliography_dfc.docx")
tryCatch({
  rmarkdown::pandoc_convert(bib_md, to = "docx", output = bib_docx)
  message("Saved: ", bib_docx)
}, error = function(e) {
  message("Bibliography docx conversion failed (requires pandoc): ", conditionMessage(e))
})


##########
# do mtg
g_overlaps <- find_comps(mtg_list)

# Run GSEA
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
b_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = T)
g_out <- run_gsea_for_proj_optimized(tapply(g_overlaps[,2], g_overlaps[,1],list),
                                     broad = F)
gsea_out <- rbind(b_out, g_out)

# Significant genesets per celltype (p < GSEA_PTHRESH)
# GSEA_PTHRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
# }), celltypes)

# ── AD gene database: genes in g_overlaps with documented AD associations ────
# Each entry: gene = list(n_refs = approx. no. relevant papers, ref = bibliography no.)
# Ordered alphabetically within each celltype when displayed.
# Bibliography: gsea_summary_bibliography_mtg.md
# NOTE: All genes in mtg_overlaps are tubulin isoforms (TUBA*/TUBB*, excluded by criteria)
#       or ribosomal proteins (RPL*/RPS*, excluded by criteria). No entries qualify.
library(gt)

ad_db <- list()

# ── Build table ─────────────────────────────────────────────────
celltypes  <- unique(g_overlaps$Celltype)

# Significant genesets per celltype (BH-adjusted p < GSEA_FDR_THRESH)
# GSEA_FDR_THRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   fdr   <- p.adjust(pvals, method = "BH")
#   gsea_out$SetName[!is.na(fdr) & fdr < GSEA_FDR_THRESH]
# }), celltypes)
# qsave(gsea_sig, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v6/gsea_sig_mtg.qs"))

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

  # AD genes: match uppercase gene names against ad_db, sort alphabetically
  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(ad_db)]
  hits_up     <- genes_upper[genes_upper %in% names(ad_db)]

  if (length(hits) > 0) {
    ref_nums <- sapply(hits_up, function(g) ad_db[[g]]$ref)
    ord      <- order(hits_up)
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

bib_md   <- file.path(save_dir, "gsea_summary_bibliography_mtg.md")
bib_docx <- file.path(save_dir, "gsea_summary_bibliography_mtg.docx")
tryCatch({
  rmarkdown::pandoc_convert(bib_md, to = "docx", output = bib_docx)
  message("Saved: ", bib_docx)
}, error = function(e) {
  message("Bibliography docx conversion failed (requires pandoc): ", conditionMessage(e))
})

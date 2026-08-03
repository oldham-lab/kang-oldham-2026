library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v1/"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_SCZ_AD/"), version)
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
#sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_rosmapAD_1127.csv"))
#these_mods_case <- these_mods_case[!these_mods_case %in% which(sigcount_bonf$vals < 1)]

# Load SCZ module data
datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods_scz <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods_scz,length))
these_mods_scz <- as.numeric(names(mods_scz)[which(modulelengths>filter_under)])


dat_vec <- list(
  # con vs All DFC
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv",
  ## con vs SCZ (CTRL modules)
  # CMC
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC/DFC/euclidean_distances/Schizophrenia_vs_control_bulk_megaset_output_table.csv",
  # SZBD
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq/DFC/euclidean_distances/Schizophrenia_vs_control_bulk_megaset_output_table.csv",
  # con vs All DFC ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"),
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv",
  ## con vs SCZ (SCZ modules)
  # CMC
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC/DFC/euclidean_distances/Schizophrenia_vs_control_brainseq_scz_output_table.csv",
  # SZBD
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq/DFC/euclidean_distances/Schizophrenia_vs_control_brainseq_scz_output_table.csv"
)

dat_vec_names <- c("Gabitto_AllADVsCon_DFC",
                   "Liu_AllADVsCon_DFC",
                   "CMC_DFC",
                   "SZBD_DFC",
                   "Gabitto_AllADVsCon_DFC_ROSMAP",
                   "Liu_AllADVsCon_DFC_ROSMAP",
                   "CMC_DFC_SCZmods",
                   "SZBD_DFC_SCZmods")

# Load all relevant dcopa output tables (for all comparisons)
dcopa_sharedctrl <- lapply(dat_vec[1:4], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency))
})

dcopa_sharedctrl[3:4] <- lapply(dcopa_sharedctrl[3:4], \(x){
   x |> 
     dplyr::filter(!Celltype %in% c("Immune", "SMC", "PC")) |>
     dplyr::mutate(Celltype = case_match(
        Celltype,
        "Oligo" ~ "Oligodendrocyte",
        "Endo" ~ "Endothelial",
        "Micro" ~ "Microglia-PVM",
        "Astro" ~ "Astrocyte",
        .default = Celltype
     ))
})

dcopa_sharedcase <- lapply(dat_vec[5:6], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods_case,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency))
})

dcopa_sharedcase[3:4] <- lapply(dat_vec[7:8], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods_scz,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency)) |> 
    dplyr::filter(!Celltype %in% c("Immune", "SMC", "PC")) |>
    dplyr::mutate(Celltype = case_match(
        Celltype,
        "Oligo" ~ "Oligodendrocyte",
        "Endo" ~ "Endothelial",
        "Micro" ~ "Microglia-PVM",
        "Astro" ~ "Astrocyte",
        .default = Celltype
     ))
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
dcopa_allct_case <- lapply(dcopa_sharedcase[1:2], \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names[5:6])
dcopa_allct_case[3:4] <- lapply(dcopa_sharedcase[3:4], \(x) extract_genes(x, mods_scz)) 
names(dcopa_allct_case)[3:4] <- dat_vec_names[7:8]

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
gsea_input_ctrl1 <- dcopa_to_df(dcopa_allct_ctrl[1:2]) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = "AD",
         "Module type" = "Bulk megaset",
         "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu")
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))

gsea_input_ctrl2 <- dcopa_to_df(dcopa_allct_ctrl[3:4]) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = "SCZ",
         "Module type" = "Bulk megaset",
         Dataset = ifelse(grepl("CMC", Comparison), "CMC", "SZBD")
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))


gsea_input_case1 <- dcopa_to_df(dcopa_allct_case[1:2]) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = "AD",
         "Module type" = "ROSMAP AD",
         Dataset = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu")
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))

gsea_input_case2 <- dcopa_to_df(dcopa_allct_case[3:4]) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = "SCZ",
         "Module type" = "SCZmods",
         "Dataset" = ifelse(grepl("CMC", Comparison), "CMC", "SZBD")
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))


gsea_input <- rbind(gsea_input_ctrl1, gsea_input_ctrl2, gsea_input_case1, gsea_input_case2)

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
  find_output_overlap(dcopa_shared[[5]], dcopa_shared[[6]]), # Gabitto + Liu DFC (AD modules)
  find_output_overlap(dcopa_shared[[3]], dcopa_shared[[4]]), # CMC vs SZBD DFC
  find_output_overlap(dcopa_shared[[7]], dcopa_shared[[8]]) # CMC vs SZBD (SCZ modules)
)

dcopa_shared_overlaps[[5]] <- find_output_overlap2(dcopa_shared_overlaps[[1]], dcopa_shared_overlaps[[3]]) # CTRL modules | AD + SCZ


# Names for each dataset (plus overlaps)

# dat_vec_names_short <- c("CTRL modules | Gabitto SN",
#                          "CTRL modules | Liu SN",
#                          "CTRL modules | CMC SN",
#                          "CTRL modules | SZBDMultiseq SN",
#                          "AD modules | Gabitto SN",
#                          "AD modules | Liu SN",
#                          "SCZ modules | CMC SN",
#                          "SCZ modules | SZBDMultiseq SN",
#                          "CTRL modules | Gabitto + Liu SN",
#                          "CTRL modules | CMC + SZBDMultiseq SN",
#                          "AD modules | Gabitto + Liu SN",
#                          "SCZ modules | CMC + SZBDMultiseq SN")   

dat_vec_names_short <- c("CTRL modules | Gabitto + Liu SN",
                         "AD modules | Gabitto + Liu SN",
                         "CTRL modules | CMC + SZBDMultiseq SN",
                         "SCZ modules | CMC + SZBDMultiseq SN",
                         "CTRL modules | AD + SCZ")     

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
plotlist[[1]] <- mapply(process_dcopa_output_tab, 
                        dcopa_shared_overlaps, 
                        dat_vec_names_short, 
                        SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass)),
         comp = factor(comp, levels = rev(dat_vec_names_short)))

# Set size range for dot legend
limit_vec = c(1, max(c(plotlist[[1]]$num_sig), na.rm = T))

# Plot (4 plots total)
dvec <- c("DFC") # save suffix
cc_colors <- RColorBrewer::brewer.pal(3, "Set1")
for(i in 1){
  plo <- plotlist[[i]] |>
      filter(Direction == -1) |>
      ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
        theme_minimal() +
        geom_point(color = "black", pch = 21) +
        theme(text = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
              axis.text.y = element_text(size = 7, colour = c("#B8860B", 
                                                              "black", "black", "black", "black")),
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
  ggsave(plo, file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5.5)
  ggsave(plo, file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], ".svg")), height = 2.5, width = 5.5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)

  phi <- plotlist[[i]] |>
      filter(Direction == 1) |>
      ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
        theme_minimal() +
        geom_point(color = "black", pch = 21) +
        theme(text = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
              axis.text.y = element_text(size = 7, colour = c("#B8860B", "black", "black", "black", "black")),
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
  ggsave(phi, file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5.5)
  ggsave(phi, file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], ".svg")), height = 2.5, width = 5.5)
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


dcopa_allct1 <- lapply(dcopa_shared_overlaps[c(1)], \(x) extract_genes(x, mods)) |> setNames(dat_vec_names_short[c(1)])
dcopa_allct2 <- lapply(dcopa_shared_overlaps[c(3)], \(x) extract_genes(x, mods)) |> setNames(dat_vec_names_short[c(3)])
dcopa_allct <- c(dcopa_allct1, dcopa_allct2)

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
mutate(Region = "DFC",
        Disease = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
        "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
        "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu")
        #Celltype = paste(Celltype, "all")
        ) |>
dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Gene")))
fwrite(gsea_input, file = file.path(save_dir, paste0("/panel_C_D_", ct,  "_dcopa_genelist.csv")))


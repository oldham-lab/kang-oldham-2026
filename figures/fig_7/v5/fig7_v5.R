library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v5"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)


#########################################
# Upset plot for DE genes (pseudobulk DE)
# Gabitto vs Liu, DFC vs MTG, DESeq2 vs edgeR
# 8 comparisons total
###########################################

# see /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_7/v3/panel_A_upset_DE.R

###########
# Panel B
############
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
class_info_mit <- fread(data.table = F, file = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv") |> 
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
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv",
  # con vs All MTG
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv", 
  # con vs All DFC ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"),
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv", 
  # all AD vs con MTG ROSMAP
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"), 
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"
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

# /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_7/v3/panel_B_upset_dcopa.R

################
# Panel C + D alternate
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
      select(c(mod, Celltype, Direction, Consistency))
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
dat_vec_names_short <- c("Gabitto DFC",
                         "Liu DFC",
                         "Gabitto MTG",
                         "Liu MTG",
                         "Gabitto DFC (AD modules)",
                         "Liu DFC (AD modules)",
                         "Gabitto MTG (AD modules)",
                         "Liu MTG (AD modules)",
                         "Gabitto + Liu DFC",
                         "Gabitto + Liu MTG",
                         "Gabitto + Liu DFC (AD modules)",
                         "Gabitto + Liu MTG (AD modules)",
                         "Gabitto + Liu, DFC + MTG",
                         "Gabitto + Liu, DFC + MTG (AD modules)")

# Function for processing dcopa output into plotting input
process_dcopa_output_tab <- function(sum_tab, dat_name){
  class_info_df <- class_info
  allsum <- sum_tab |> 
  mutate(Direction = factor(Direction, levels = c(-1, 1)),
          Celltype = factor(Celltype, levels = unique(class_info_df$Subclass))) |>
  select(mod, Celltype, Direction) |>
  group_by(Celltype, Direction, .drop = FALSE) |>
  summarise(num_sig = n(), .groups = "drop") |>
  left_join(class_info_df |> select(Subclass, Class), by = join_by(Celltype == Subclass)) |>
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
  ggsave(plo, file = file.path(save_dir, paste0("panel_C_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  #ggsave(plo, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_C_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_C_indiv_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_C_indiv_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)

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
  ggsave(phi, file = file.path(save_dir, paste0("panel_D_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  #ggsave(phi, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, paste0("panel_D_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_D_indiv_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_D_indiv_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)
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
ggsave(p_legend_out, file = file.path(save_dir, "panel_C_D_legend.svg"),
        height = 0.8, width = 3, bg = "transparent")
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
# Panel F Upset overlap
#######
d1 <- dcopa_shared_overlaps[[5]]
d2 <- dcopa_shared_overlaps[[6]]
common_cts <- unique(c(d1[,2], d2[,2]))

# Collect lists of genes per subclass (Gabitto+Liu DFC/MTG, Gabitto+Liu DFC/MTG ROSMAP)
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
  }
  names(sublist) <- common_cts
  sublist <- sublist[-which(unlist(lapply(sublist, \(x) "none" %in% x)))]
  return(sublist)
}) 
names(outlist) <- c(-1, 1)

# DFC (lower in disease)
local({
  orig_fn <- getFromNamespace("Make_main_bar", "UpSetR")
  fn_src  <- paste(deparse(body(orig_fn)), collapse = "\n")
  fn_src  <- sub("bottom_margin <- \\(-1\\) \\* 0\\.65", "bottom_margin <- 0", fn_src)
  patched_fn       <- orig_fn
  body(patched_fn) <- parse(text = fn_src)[[1]]
  assignInNamespace("Make_main_bar", patched_fn, "UpSetR")
  on.exit(assignInNamespace("Make_main_bar", orig_fn, "UpSetR"), add = TRUE)

  p1 <- upset(fromList(outlist[[1]]),
              sets = names(outlist[[1]]),
              order.by = "freq",
              mainbar.y.label = "# of overlaps",
              sets.x.label = "# of unique\ndCoPA genes",
              mb.ratio = c(0.5, 0.5),
              show.numbers = F)

  svg(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_F_upset.svg"),
      width = 2.4, height = 2)
  print(p1)
  dev.off()

  pdf(file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "/panel_F_upset.pdf"),
      width = 2.4, height = 2)
  print(p1)
  dev.off()
})

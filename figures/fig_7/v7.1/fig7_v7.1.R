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

version <- "v7.1/"
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

ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "SST CHODL", "SST", "Chandelier", 
              "PVALB", "LAMP5", "LAMP5 LHX6", "PAX6",  "SNCG", 
              "VIP",  "Astrocyte", "Oligodendrocyte", "OPC", "Microglia-PVM", "Endothelial", "VLMC")

# Load class info + subclass labels for Gabitto
# class_info <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
#   select(Subclass, Class) |>
#   filter(!duplicated(Subclass)) |>
#   mutate(Class = case_match(
#     Class, 
#     "Neuronal: GABAergic" ~ "GABAergic",
#     "Neuronal: Glutamatergic" ~ "Glutamatergic",
#     "Non-neuronal and Non-neural" ~ "Non-neuronal"
#   )) |>
#   arrange(Subclass) |>
#   mutate(Subclass_fixed = factor(c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
#                                    "L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","OPC","Oligo",  "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","VLMC","Vip"),
#                                  levels = c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
#                                              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
#                                              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro/PVM", "Endo", "VLMC"))) |>
#   arrange(Subclass_fixed) |>
#   select(Subclass, Class)

class_info <- data.frame(
  "Subclass" = c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT",
                "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6",
                "Pvalb", "Sncg", "Sst", "Vip", "Sst Chodl", "Astrocyte", "Oligodendrocyte",
                "OPC", "Microglia-PVM", "Endothelial", "VLMC"),
  "Class" = c("Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic",
              "Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic",
              "Glutamatergic", "GABAergic",     "GABAergic",     "GABAergic",    
              "GABAergic",     "GABAergic",     "GABAergic",     "GABAergic",    
              "GABAergic",     "GABAergic",     "Non-neuronal",  "Non-neuronal", 
              "Non-neuronal",  "Non-neuronal",  "Non-neuronal",  "Non-neuronal")
)

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
mods_trim <- mods[these_mods]

# Load AD module data
datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods_case <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods_case,length))
these_mods_case <- as.numeric(names(mods_case)[which(modulelengths>filter_under)])
#sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_rosmapAD_1127.csv"))
#these_mods_case <- these_mods_case[!these_mods_case %in% which(sigcount_bonf$vals < 1)]

# Gene -> module-index lookup (each gene belongs to exactly one module).
# CTRL genes come from mods_trim; AD/case genes come from the trimmed case mods.
mods_case_trim <- mods_case[these_mods_case]
make_gene2mod <- function(mod_list) {
  setNames(rep(names(mod_list), lengths(mod_list)), unlist(mod_list, use.names = FALSE))
}
gene2mod_ctrl <- make_gene2mod(mods_trim)
gene2mod_case <- make_gene2mod(mods_case_trim)

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
         "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu"),
         Module = unname(gene2mod_ctrl[Gene])
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Module", "Gene")))

gsea_input_case <- dcopa_to_df(dcopa_allct_case) |>
  mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
         Disease = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
         "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
         "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu"),
         Module = unname(gene2mod_case[Gene])
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Module", "Gene")))

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

# Save overlap tables for shiny app
# Capitalize celltype names
# NOTE: this plotlist recode is applied later (after plotlist is built, ~L349); the
# copy here referenced plotlist before it existed and crashed fresh Rscript runs. Disabled.
# plotlist <- plotlist |> lapply(\(x) x |>
#   mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME),
#          Celltype = factor(Celltype, levels = ct_order),
#          Class = factor(Class, levels = c("Glutamatergic", "GABAergic", "Non-neuronal"))))

dcopa_shared_overlaps_save <- lapply(dcopa_shared_overlaps, \(x){
  x_temp <- x[!duplicated(x$mod), ]
  x_temp <- c(table(x_temp$Direction))
  higher_count <- x_temp[names(x_temp) == 1]
  lower_count <- x_temp[names(x_temp) == -1]
  out <- x |> 
    dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME),
         Celltype = factor(Celltype, levels = ct_order)) |>
    dplyr::mutate(type = case_match(
      Direction,
      -1 ~ paste0("Lower in all AD vs con\n(n = ", lower_count, ")"),
      1 ~ paste0("Higher in all AD vs con\n(n = ", higher_count, ")"),
      .default = NA
  ))
  return(out)
})
fwrite(dcopa_shared_overlaps_save[[1]], file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/ad_dfc.csv"))
fwrite(dcopa_shared_overlaps_save[[2]], file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/ad_mtg.csv"))
fwrite(dcopa_shared_overlaps_save[[3]], file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/ad_dfc_admods.csv"))
fwrite(dcopa_shared_overlaps_save[[4]], file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_svg_map_path/ad_mtg_admods.csv"))


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

# Save tables for dot_data_path for Shiny app
plotlist_save_create <- function(df, comp_var, overlap_df){
  x <- df |> dplyr::filter(comp == comp_var)
  
  overlap_temp <- overlap_df[!duplicated(overlap_df$Direction), ]
  hi_var <- overlap_temp$type[overlap_temp$Direction == 1]
  lo_var <- overlap_temp$type[overlap_temp$Direction == -1]
  if(length(hi_var) == 0) hi_var <- NA
  if(length(lo_var) == 0) lo_var <- NA

  out <- x |> dplyr::mutate(type = case_match(
    Direction,
    "-1" ~ lo_var,
    "1" ~ hi_var,
    .default = NA
  ))
  return(out)
}

fwrite(plotlist_save_create(plotlist[[2]], "CTRL modules | Gabitto + Liu SN", dcopa_shared_overlaps_save[[1]]),
       file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_dot_data_path/ad_dfc.csv")     
) # DFC ctrl mods
fwrite(plotlist_save_create(plotlist[[1]], "CTRL modules | Gabitto + Liu SN", dcopa_shared_overlaps_save[[2]]),
       file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_dot_data_path/ad_mtg.csv")     
) # MTG ctrl mods
fwrite(plotlist_save_create(plotlist[[2]], "AD modules | Gabitto + Liu SN", dcopa_shared_overlaps_save[[3]]),
       file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_dot_data_path/ad_dfc_admods.csv")     
) # DFC AD mods
fwrite(plotlist_save_create(plotlist[[1]],  "AD modules | Gabitto + Liu SN", dcopa_shared_overlaps_save[[4]]),
       file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/dcopa_dot_data_path/ad_mtg_admods.csv")     
) # MTG AD mods

# Set size range for dot legend
limit_vec = c(1, max(c(plotlist[[1]]$num_sig, plotlist[[2]]$num_sig), na.rm = T))

# Plot (4 plots total)
dvec <- c("MTG", "DFC") # save suffix
cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

# Square-bracket layers spanning the two combined "Gabitto + Liu SN" rows (bottom
# two rows, y = 1-2), labelled with the Pearson correlation of dot sizes (num_sig)
# between those two rows across cell types (complete cell types only).
add_r_bracket <- function(df, dir_val) {
  d <- df |>
    dplyr::filter(Direction == dir_val,
                  comp %in% c("CTRL modules | Gabitto + Liu SN",
                              "AD modules | Gabitto + Liu SN")) |>
    dplyr::select(Celltype, comp, num_sig) |>
    tidyr::pivot_wider(names_from = comp, values_from = num_sig)
  r <- cor(d[["CTRL modules | Gabitto + Liu SN"]],
           d[["AD modules | Gabitto + Liu SN"]], use = "complete.obs")
  xb   <- nlevels(df$Celltype) + 1   # bracket vertical bar, just right of last cell type
  tick <- 0.4                        # length of the bracket's horizontal ticks
  list(
    annotate("segment", x = xb, xend = xb,        y = 1, yend = 2),
    annotate("segment", x = xb, xend = xb - tick, y = 1, yend = 1),
    annotate("segment", x = xb, xend = xb - tick, y = 2, yend = 2),
    annotate("text", x = xb + 0.7, y = 1.5,
             label = sprintf("r = %.2f", r), angle = 270, size = 7 / .pt),
    coord_cartesian(clip = "off"),
    theme(plot.margin = margin(2, 16, 2, 2))
  )
}

for(i in 1:2){
  plo <- plotlist[[i]] |>
      filter(Direction == -1) |>
      ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
        theme_minimal() +
        geom_point(color = "black", pch = 21) +
        theme(text = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
              axis.text.y = element_text(size = 7, colour = c("#B8860B", "#B8860B", "black", "black", "black", "black")),
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
                              breaks = scales::pretty_breaks(n = 3)) +
        add_r_bracket(plotlist[[i]], -1)
  ggsave(plo, file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  ggsave(plo, file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)

  phi <- plotlist[[i]] |>
      filter(Direction == 1) |>
      ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
        theme_minimal() +
        geom_point(color = "black", pch = 21) +
        theme(text = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
              axis.text.y = element_text(size = 7, colour = c("#B8860B", "#B8860B", "black", "black", "black", "black")),
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
  ggsave(phi, file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
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

  m  <- make_comb_mat(dcopa_allct)
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
  paste0("/panel_E_F_", reg, ".svg")), width = 6, height = 3.5)                                     
  draw(p)         
  dev.off()
}
# MTG: 368/544 (68%) of genes are unique to a single subclass
# DFC: 165/198 (83%) of genes are unique to a single subclass

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

# Calculate FDR
gsea_pval_to_fdr <- function(gsea_out) {                              
    meta_cols <- c("SetID", "SetName", "SetSize")                       
    pval_cols <- setdiff(colnames(gsea_out), meta_cols)                 
    fdr_out <- gsea_out                                                 
    fdr_out[, pval_cols] <- lapply(gsea_out[, pval_cols, drop = FALSE], 
                                    p.adjust, method = "BH")            
    fdr_out                                                             
  }                                                                     
gsea_out_fdr <- gsea_pval_to_fdr(gsea_out)


# Significant genesets per celltype (p < GSEA_PTHRESH)
# GSEA_PTHRESH <- 0.05
# gsea_sig <- setNames(lapply(seq_along(celltypes), function(i) {
#   pvals <- gsea_out[[3 + i]]
#   gsea_out$SetName[!is.na(pvals) & pvals < GSEA_PTHRESH]
# }), celltypes)

# ── AD gene database: genes in g_overlaps with documented AD associations ────
# Each entry: gene = list(n_refs = approx. no. relevant papers, ref = bibliography no.)
# Ordered alphabetically within each celltype when displayed.
# Bibliography: gsea_summary_bibliography.md
library(gt)

ad_db <- list(
  # PI3K–Akt signalling
  AKT3     = list(n=7,  ref=1),
  # Endosomal APP sorting / Rab5
  APPL1    = list(n=5,  ref=2),
  # APP secretory trafficking (COPII)
  APPBP2   = list(n=5,  ref=3),
  SEC23A   = list(n=2,  ref=3),
  # COPI / APP retrograde Golgi-to-ER trafficking
  COPA     = list(n=3,  ref=4),
  # Axonal transport / dynein-dynactin
  DCTN4    = list(n=4,  ref=5),
  PAFAH1B1 = list(n=6,  ref=5),
  # Mitochondrial fission
  DNM1L    = list(n=8,  ref=6),
  # ER stress / ERAD
  EDEM3    = list(n=3,  ref=7),
  # JNK / tau hyperphosphorylation
  MAPK8    = list(n=4,  ref=8),
  # Mitochondrial complex I / bioenergetics
  NDUFS1   = list(n=5,  ref=9),
  # Pannexin-1 / synaptic plasticity
  PANX1    = list(n=3,  ref=10),
  # Tau phosphatase (PP2A)
  PPP2R2B  = list(n=8,  ref=11),
  # Tau phosphatase (calcineurin)
  PPP3CB   = list(n=5,  ref=12),
  # Proteasome assembly and regulation
  PSMG1    = list(n=3,  ref=13),
  PSMD5    = list(n=3,  ref=13),
  # BACE1 restriction / amyloid processing
  RTN3     = list(n=5,  ref=14),
  # PI(4)P lipid homeostasis / APP trafficking
  SACM1L   = list(n=2,  ref=15),
  # SUMO / tau solubility
  SENP1    = list(n=4,  ref=16),
  # Endophilin A1 / synaptic Aβ injury
  SH3GL2   = list(n=3,  ref=17),
  # Cytoskeletal integrity / calpain cleavage
  SPTAN1   = list(n=3,  ref=18),
  # Store-operated Ca2+ entry / spine stability
  STIM2    = list(n=4,  ref=19),
  # Tau deubiquitination / aggregation
  USP10    = list(n=3,  ref=20)
)

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
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (p=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
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

out_path <- file.path(save_dir, "panel_G_gsea_summary_table.html")
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

# ── FDR-corrected version (DFC) ─────────────────────────────────────────────
build_row_fdr <- function(ct) {
  genes   <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out_fdr[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out_fdr$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out_fdr$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (padj=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
  }

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

summary_df_fdr <- do.call(rbind, lapply(celltypes, build_row_fdr)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

gt_table_fdr <- summary_df_fdr |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by FDR-adjusted p-value)",
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

out_path <- file.path(save_dir, "panel_G_gsea_summary_table_dfc_fdr.html")
gtsave(gt_table_fdr, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_dfc_fdr.pdf")
tryCatch({
  gtsave(gt_table_fdr, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
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

# Calculate FDR
gsea_pval_to_fdr <- function(gsea_out) {                              
    meta_cols <- c("SetID", "SetName", "SetSize")                       
    pval_cols <- setdiff(colnames(gsea_out), meta_cols)                 
    fdr_out <- gsea_out                                                 
    fdr_out[, pval_cols] <- lapply(gsea_out[, pval_cols, drop = FALSE], 
                                    p.adjust, method = "BH")            
    fdr_out                                                             
  }                                                                     
gsea_out_fdr <- gsea_pval_to_fdr(gsea_out)

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
library(gt)

ad_db <- list(
  # PI3K–Akt signalling
  AKT3     = list(n=7,  ref=1),
  PIK3CA   = list(n=5,  ref=1),
  # Endosomal APP sorting / Rab5
  APPL1    = list(n=5,  ref=2),
  # APP cytoplasmic tail / secretory trafficking
  APPBP2   = list(n=5,  ref=3),
  # Mitophagy receptor (PINK1/Parkin-independent)
  BNIP3    = list(n=4,  ref=4),
  # Autophagy / ALS-FTD endosomal pathway
  C9orf72  = list(n=4,  ref=5),
  # Kinesin-1 adaptor / APP axonal transport
  CLSTN1   = list(n=5,  ref=6),
  # Dystrophic neurites / APP vulnerability
  CLSTN3   = list(n=3,  ref=7),
  # Endocytosis / clathrin
  CLTC     = list(n=7,  ref=8),
  # COPI / APP retrograde trafficking
  COPA     = list(n=3,  ref=9),
  # Ubiquitin-proteasome system (SCF-E3 ligase)
  CUL1     = list(n=3,  ref=10),
  KLHL7    = list(n=2,  ref=10),
  # Clathrin uncoating / endocytic recycling
  DNAJC6   = list(n=3,  ref=11),
  # Mitochondrial fission and fusion
  DNM1L    = list(n=8,  ref=12),
  OPA1     = list(n=5,  ref=12),
  # Retrograde axonal transport / dynein
  DYNC1H1  = list(n=5,  ref=13),
  DYNC1I1  = list(n=5,  ref=13),
  PAFAH1B1 = list(n=6,  ref=13),
  # ER stress / ERAD
  EDEM3    = list(n=3,  ref=14),
  # Stress granules / tau condensates
  G3BP2    = list(n=3,  ref=15),
  # Autophagy / LC3-II
  MAP1LC3B = list(n=8,  ref=16),
  # DENN/MADD / neuroprotection
  MADD     = list(n=4,  ref=17),
  # JNK / tau hyperphosphorylation
  MAPK8    = list(n=4,  ref=18),
  # Neprilysin / Abeta degradation
  MME      = list(n=5,  ref=19),
  # mTOR / autophagy / tau clearance
  MTOR     = list(n=5,  ref=20),
  # Neurexin / synaptic organisation
  NRXN1    = list(n=4,  ref=21),
  # Oxidative stress / neuroprotection
  OXR1     = list(n=3,  ref=22),
  # Endolysosomal PI(3,5)P2 synthesis
  PIKFYVE  = list(n=4,  ref=23),
  # Tau phosphatase (PP2A)
  PPP2R2B  = list(n=8,  ref=24),
  # Tau phosphatase (calcineurin)
  PPP3CB   = list(n=5,  ref=25),
  PPP3R1   = list(n=3,  ref=25),
  # BACE1-APP scaffolding
  RANBP9   = list(n=5,  ref=26),
  # Selective vulnerability / NFT-prone neurons
  RORB     = list(n=4,  ref=27),
  # ER morphology / BACE1 restriction
  RTN1     = list(n=3,  ref=28),
  RTN3     = list(n=5,  ref=28),
  # PI(4)P lipid homeostasis / APP trafficking
  SACM1L   = list(n=2,  ref=29),
  # COPII / APP ER-to-Golgi trafficking
  SEC23A   = list(n=2,  ref=30),
  # Endophilin A1 / synaptic Aβ injury
  SH3GL2   = list(n=3,  ref=31),
  # Synaptic vesicle exocytosis / SNARE
  SNAP25   = list(n=5,  ref=32),
  # Cytoskeletal integrity / calpain cleavage
  SPTAN1   = list(n=3,  ref=33),
  # Store-operated Ca2+ entry / spine stability
  STIM2    = list(n=4,  ref=34),
  # Synaptic phosphoinositide turnover
  SYNJ1    = list(n=9,  ref=35),
  # Synaptic vesicle density / AD PET biomarker
  SV2A     = list(n=5,  ref=36),
  # Tau kinase (TAOK / MARK cascade)
  TAOK1    = list(n=5,  ref=37),
  # Neuroinflammation / NF-kB / TLR
  TBK1     = list(n=4,  ref=38),
  # Lysosomal / TDP-43 pathology
  TMEM106B = list(n=9,  ref=39),
  # Selective autophagy / TLR adaptor
  TOLLIP   = list(n=5,  ref=40),
  # Tau deubiquitination / aggregation
  USP10    = list(n=3,  ref=41),
  # Proteasomal tau clearance
  USP14    = list(n=5,  ref=42),
  # ER-mitochondria contact sites (MAM)
  VAPB     = list(n=4,  ref=43),
  # Mitochondria / Abeta interaction
  VDAC1    = list(n=5,  ref=44),
  VDAC3    = list(n=4,  ref=45),
  # Retromer / BACE1 endosomal trafficking
  VPS35    = list(n=5,  ref=46),
  # Ca2+ sensor / AD CSF biomarker
  VSNL1    = list(n=6,  ref=47),
  # 14-3-3 / phospho-tau interaction
  YWHAB    = list(n=5,  ref=48),
  YWHAZ    = list(n=5,  ref=49)
)

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
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (p=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
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

# ── FDR-corrected version (MTG) ─────────────────────────────────────────────
build_row_fdr <- function(ct) {
  genes   <- g_overlaps$Gene[g_overlaps$Celltype == ct]
  n_genes <- length(genes)

  ct_idx  <- which(celltypes == ct)
  pvals   <- gsea_out_fdr[[3 + ct_idx]]
  go_mask <- grepl("^GO", gsea_out_fdr$SetName)
  if (!is.null(pvals)) {
    go_pvals <- ifelse(go_mask, pvals, NA)
    sig_mask <- !is.na(go_pvals) & go_pvals < 0.05
    n_sig    <- sum(sig_mask)
    if (n_sig == 0) {
      top_str <- "—"
    } else {
      top_idx  <- order(go_pvals)[1:min(3, n_sig)]
      top_sets <- gsea_out_fdr$SetName[top_idx]
      top_pval <- pvals[top_idx]
      top_str  <- paste0(
        seq_along(top_sets), ". ", top_sets,
        " (padj=", formatC(top_pval, format = "e", digits = 2), ")",
        collapse = "<br>"
      )
    }
  }

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

summary_df_fdr <- do.call(rbind, lapply(celltypes, build_row_fdr)) |>
  arrange(Celltype) |>
  dplyr::mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

gt_table_fdr <- summary_df_fdr |>
  gt() |>
  cols_label(
    Celltype     = "Cell type",
    N_genes      = "Genes (n)",
    Top_genesets = "Top 3 GO gene sets (by FDR-adjusted p-value)",
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

out_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg_fdr.html")
gtsave(gt_table_fdr, out_path)
message("Saved: ", out_path)

pdf_path <- file.path(save_dir, "panel_G_gsea_summary_table_mtg_fdr.pdf")
tryCatch({
  gtsave(gt_table_fdr, pdf_path)
  message("Saved: ", pdf_path)
}, error = function(e) {
  message("PDF export failed (requires webshot2 + Chrome): ", conditionMessage(e))
})

############
# Additional analyses
############
save_dir_add <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/"), version, "additional")


### Create upset plots for bottom two rows of a-b (4 total)
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

save_suffixes <- c("CTRLmods_DFC", "CTRLmods_MTG", "ADmods_DFC", "ADmods_MTG")
for(i in 1:4){
  if(i %in% 1:2){
    dcopa_allct <- extract_genes(dcopa_shared_overlaps[[c(i)]], mods)[[1]] 
  } else {
    dcopa_allct <- extract_genes(dcopa_shared_overlaps[[c(i)]], mods_case)[[1]] 
  }

  names(dcopa_allct) <- dplyr::recode(names(dcopa_allct), !!!CT_RENAME)

  m  <- make_comb_mat(dcopa_allct)
  #m  <- m[comb_size(m) > 2]     # drop combinations with ≤1 overlapping gene
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

  SET_ORDER <- ct_order[ct_order %in% names(dcopa_allct)]

  p <- UpSet(
    m,            
    # sets = names(dcopa_allct) equivalent        
    set_order  = SET_ORDER,              
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

  if(i %in% c(2,4)){
    plotwidth <- 10
  } else if (i %in% c(1,3)){
    plotwidth <- 6
  }

  svglite::svglite(file.path(save_dir_add,                     
  paste0("/additional_upset_", save_suffixes[i], ".svg")), width = plotwidth, height = 4)                                     
  draw(p)         
  dev.off()
}

######### Create summary tables for bottom two rows of a-b (4 total)
# Load genelists
dfc_list <- fread(data.table = F, file = file.path(save_dir, "panel_C_D_DFC_dcopa_genelist.csv"))
mtg_list <- fread(data.table = F, file = file.path(save_dir, "/panel_C_D_MTG_dcopa_genelist.csv"))

# Save gene tables for each comparison separately
fwrite(dfc_list |> 
  dplyr::filter(Comparison == "CTRL modules | Gabitto + Liu SN"), 
  file = file.path(save_dir_add, "CTRLmods_dfc_overlaps.csv"))
fwrite(dfc_list |> 
  dplyr::filter(Comparison == "AD modules | Gabitto + Liu SN"), 
  file = file.path(save_dir_add, "ADmods_dfc_overlaps.csv"))
fwrite(mtg_list |> 
  dplyr::filter(Comparison == "CTRL modules | Gabitto + Liu SN"), 
  file = file.path(save_dir_add, "CTRLmods_mtg_overlaps.csv"))
fwrite(mtg_list |> 
  dplyr::filter(Comparison == "AD modules | Gabitto + Liu SN"), 
  file = file.path(save_dir_add, "ADmods_mtg_overlaps.csv"))

# Create ad_db objects
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.1/summary_table_fxn.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.1/additional/ad_db_objs.R"))

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


g_overlap_list <- list(
  fread(file.path(save_dir_add, "CTRLmods_dfc_overlaps.csv"), data.table = F),
  fread(file.path(save_dir_add, "CTRLmods_mtg_overlaps.csv"), data.table = F),
  fread(file.path(save_dir_add, "ADmods_dfc_overlaps.csv"), data.table = F),
  fread(file.path(save_dir_add, "ADmods_mtg_overlaps.csv"), data.table = F)
) |> lapply(\(dcopa){                                                          
  dcopa |>           
    distinct(Celltype, Gene) |>                                                                                    
    arrange(Celltype, Gene)     
})

bla <- fread(data.table=F,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.1/mtg_overlaps.csv"))
bla1 <- bla |> filter(Celltype=="L6 IT")
bla2 <- g_overlap_list[[2]] |> filter(Celltype=="L6 IT")

create_summary_table(g_overlaps = g_overlap_list[[1]],
                     ad_db = ad_db_CTRLmods_dfc,
                     column_4_header = "AD-associated genes",
                     out_string = "CTRLmods_dfc",
                     save_dir = save_dir_add)

create_summary_table(g_overlaps = g_overlap_list[[2]],
                     ad_db = ad_db_CTRLmods_mtg,
                     column_4_header = "AD-associated genes",
                     out_string = "CTRLmods_mtg",
                     save_dir = save_dir_add)

create_summary_table(g_overlaps = g_overlap_list[[3]],
                     ad_db = ad_db_ADmods_dfc,
                     column_4_header = "AD-associated genes",
                     out_string = "ADmods_dfc",
                     save_dir = save_dir_add)

create_summary_table(g_overlaps = g_overlap_list[[4]],
                     ad_db = ad_db_ADmods_mtg,
                     column_4_header = "AD-associated genes",
                     out_string = "ADmods_mtg",
                     save_dir = save_dir_add)
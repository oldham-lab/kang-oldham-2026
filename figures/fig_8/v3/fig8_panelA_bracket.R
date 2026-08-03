# v1
# - added representative snapshot for VIP/SST


library(tidyverse)
library(qs)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)
library(showtext)
showtext_auto()

version <- "v3"   # v3 panel-a-only run: outputs (incl. bracketed panel A) go to fig_8/v3
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/"), version)

HIGHLIGHT_GENE <- "RB1CC1"  # gene to highlight in red in the SCZ column; set to "" to disable
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

##########
# Load initial data
###########

# Dataframe for converting SCZ celltypes to Gabitto celltypes
ct_scz_map <- data.frame("scz" = c("Astro", "Endo", "Micro", "Oligo"),
                         "gab" = c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte"))

# Capitalizing certain gabitto celltypes
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

# Final celltype order
ct_order <- c("L2/3 IT", "L4 IT",  "L5 ET", "L5 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "SST CHODL", "SST", "Chandelier", 
              "PVALB", "LAMP5", "LAMP5 LHX6", "PAX6",  "SNCG", 
              "VIP",  "Astrocyte", "Oligodendrocyte", "OPC", "Microglia-PVM", "Endothelial", "VLMC")

# Load class info + subclass labels for Gabitto
class_info <- data.frame(
    Subclass = c(
      "L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
      "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Vip", "Sst Chodl",
      "Astrocyte", "Oligodendrocyte", "OPC", "Microglia-PVM", "Endothelial", "VLMC"
    ),
    Class = c(
      "Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic",
      "Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic",
      "GABAergic", "GABAergic", "GABAergic", "GABAergic", "GABAergic",
      "GABAergic", "GABAergic", "GABAergic", "GABAergic",
      "Non-neuronal", "Non-neuronal", "Non-neuronal", "Non-neuronal", "Non-neuronal", "Non-neuronal"
    )
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

# Load SCZ module data
datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods_case <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods_case,length))
these_mods_case <- as.numeric(names(mods_case)[which(modulelengths>filter_under)])

# Gene -> module-index lookup (each gene belongs to exactly one module).
# CTRL genes come from mods_trim; SCZ/case genes come from the trimmed case mods.
mods_case_trim <- mods_case[these_mods_case]
make_gene2mod <- function(mod_list) {
  setNames(rep(names(mod_list), lengths(mod_list)), unlist(mod_list, use.names = FALSE))
}
gene2mod_ctrl <- make_gene2mod(mods_trim)
gene2mod_case <- make_gene2mod(mods_case_trim)

dat_vec <- list(
  ## con vs SCZ (CTRL modules)
  # CMC
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC/DFC/euclidean_distances/Schizophrenia_vs_control_bulk_megaset_output_table.csv",
  # SZBD
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq/DFC/euclidean_distances/Schizophrenia_vs_control_bulk_megaset_output_table.csv",
  ## con vs SCZ (SCZ modules)
  # CMC
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC/DFC/euclidean_distances/Schizophrenia_vs_control_brainseq_scz_output_table.csv",
  # SZBD
  "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq/DFC/euclidean_distances/Schizophrenia_vs_control_brainseq_scz_output_table.csv"
)

dat_vec_names <- c("CMC",
                   "SZBD",
                   "CMC_SCZmods",
                   "SZBD_SCZmods")

# Load all relevant dcopa output tables (for all comparisons)
dcopa_sharedctrl <- lapply(dat_vec[1:2], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency)) |>
    mutate(Celltype = case_match(                   
      Celltype,                                     
      "Astro"  ~ "Astrocyte",
      "Endo"    ~ "Endothelial",                        
      "Micro"    ~ "Microglia-PVM",                 
      "Oligo" ~ "Oligodendrocyte",                     
      .default = Celltype                      
    )) |>
    filter(!Celltype %in% c("Immune", "PC", "SMC"))
})

dcopa_sharedcase <- lapply(dat_vec[3:4], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods_case,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency)) |>
    mutate(Celltype = case_match(                   
      Celltype,                                     
      "Astro"  ~ "Astrocyte",
      "Endo"    ~ "Endothelial",                        
      "Micro"    ~ "Microglia-PVM",                 
      "Oligo" ~ "Oligodendrocyte",                     
      .default = Celltype                      
    )) |>
    filter(!Celltype %in% c("Immune", "PC", "SMC"))
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
dcopa_allct_ctrl <- lapply(dcopa_sharedctrl, \(x) extract_genes(x, mods)) |> setNames(dat_vec_names[1:2])
dcopa_allct_case <- lapply(dcopa_sharedcase, \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names[3:4])

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
  mutate(Region = "DFC",
         Disease = "SCZ",
         "Module type" = ifelse(grepl("SCZmods", Comparison), "SCZ mods", "Bulk megaset"),
         "Dataset" = ifelse(grepl("CMC", Comparison), "CMC", "SZBD"),
         Module = unname(gene2mod_ctrl[Gene])
         #Celltype = paste(Celltype, "all")
         ) |>
  dplyr::select(all_of(c("Comparison", "Dataset", "Disease", "Region", "Module type", "Direction", "Celltype", "Module", "Gene")))

gsea_input_case <- dcopa_to_df(dcopa_allct_case) |>
  mutate(Region = "DFC",
         Disease = "SCZ",
         "Module type" = ifelse(grepl("SCZmods", Comparison), "SCZ mods", "Bulk megaset"),
         "Dataset" = ifelse(grepl("CMC", Comparison), "CMC", "SZBD"),
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
  find_output_overlap(dcopa_shared[[1]], dcopa_shared[[2]]), # CMC + SZBD DFC
  find_output_overlap(dcopa_shared[[3]], dcopa_shared[[4]]) # CMC + SZBD DFC (SCZ modules)
)

# Save tables for shiny app
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
      -1 ~ paste0("Lower in SCZ vs con\n(n = ", lower_count, ")"),
      1 ~ paste0("Higher in SCZ vs con\n(n = ", higher_count, ")"),
      .default = NA
  ))
  return(out)
})
# (v3 panel-a-only run: shared fig_7 writes disabled to avoid clobbering)
# fwrite(dcopa_shared_overlaps_save[[1]], file = ".../fig_7/dcopa_svg_map_path/scz_dfc.csv")
# fwrite(dcopa_shared_overlaps_save[[2]], file = ".../fig_7/dcopa_svg_map_path/scz_dfc_sczmods.csv")

#dcopa_shared_overlaps[[3]] <- find_output_overlap2(dcopa_shared_overlaps[[1]], dcopa_shared_overlaps[[2]]) # Gabitto/Liu DFC/MTG

# Names for each dataset (plus overlaps)


dat_vec_names_short <- c("CTRL modules | CMC SN",
                         "CTRL modules | SZBDMultiseq SN",
                         "SCZ modules | CMC SN",
                         "SCZ modules | SZBDMultiseq SN",
                         "CTRL modules | CMC + SZBDMultiseq SN",
                         "SCZ modules | CMC + SZBDMultiseq SN")                         

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
# DFC
plotlist[[1]] <- mapply(process_dcopa_output_tab, 
                        c(dcopa_shared, dcopa_shared_overlaps), 
                        dat_vec_names_short, 
                        SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(Celltype = factor(Celltype, levels = unique(class_info$Subclass)),
         comp = factor(comp, levels = rev(dat_vec_names_short[c(1,3,2,4,5,6)])))

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

# (v3 panel-a-only run: shared fig_7 dot-data writes disabled)
# fwrite(plotlist_save_create(plotlist[[1]], "CTRL modules | CMC + SZBDMultiseq SN", dcopa_shared_overlaps_save[[1]]),
#        file = ".../fig_7/dcopa_dot_data_path/scz_dfc.csv") # DFC ctrl mods
# fwrite(plotlist_save_create(plotlist[[1]], "SCZ modules | CMC + SZBDMultiseq SN", dcopa_shared_overlaps_save[[2]]),
#        file = ".../fig_7/dcopa_dot_data_path/scz_dfc_sczmods.csv") # DFC SCZ mods

# Set size range for dot legend
limit_vec = c(1, max(c(plotlist[[1]]$num_sig), na.rm = T))

# Plot (4 plots total)
dvec <- c("DFC") # save suffix
cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

# Square-bracket spanning the two combined "CMC + SZBDMultiseq SN" rows (bottom two
# rows, y = 1-2), labelled with the Pearson correlation of dot sizes (num_sig)
# between those two rows across cell types. Same style as fig_7/v8's add_r_bracket.
add_r_bracket <- function(df, dir_val) {
  d <- df |>
    dplyr::filter(Direction == dir_val,
                  comp %in% c("CTRL modules | CMC + SZBDMultiseq SN",
                              "SCZ modules | CMC + SZBDMultiseq SN")) |>
    dplyr::select(Celltype, comp, num_sig) |>
    tidyr::pivot_wider(names_from = comp, values_from = num_sig)
  r <- cor(d[["CTRL modules | CMC + SZBDMultiseq SN"]],
           d[["SCZ modules | CMC + SZBDMultiseq SN"]], use = "complete.obs")
  xb   <- nlevels(df$Celltype) + 1   # bracket bar, just right of the last cell type
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
for(i in 1){
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
        guides(fill = guide_legend(order = 1),
               size = guide_legend(order = 2)) +
        add_r_bracket(plotlist[[i]], -1)   # fig_7-style r bracket on the combined rows
  ggsave(plo, file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  ggsave(plo, file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(plo + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_A_indiv_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)
  message("v3: saved bracketed panel A; stopping before panel B / GSEA / tables.")
  quit(save = "no", status = 0)

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
        labs(x = "", y = "") +
        scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
        scale_size_continuous(limits = limit_vec, # Manually set the scale to be the same as lo object
                              breaks = scales::pretty_breaks(n = 2)) +
        guides(fill = guide_legend(order = 1),
               size = guide_legend(order = 2))
  ggsave(phi, file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], ".pdf")), height = 2.5, width = 5)
  ggsave(phi, file = file.path(save_dir, paste0("panel_B_indiv_", dvec[i], ".svg")), height = 2.5, width = 5)
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
for(ct in c("DFC")){
  dcopa_allct1 <- lapply(dcopa_shared_overlaps[c(1)], \(x) extract_genes(x, mods)) |> setNames(dat_vec_names_short[c(5)])
  dcopa_allct2 <- lapply(dcopa_shared_overlaps[c(2)], \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names_short[c(6)])
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
           Disease = "SCZ",
           "Module type" = ifelse(grepl("SCZmods", Comparison), "SCZ mods", "Bulk megaset"),
           "Dataset" = ifelse(grepl("CMC", Comparison), "CMC", "SZBD")
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

#for(dir in c(-1, 1)){
 dir <- -1
  for(reg in c("DFC")){
    dcopa_allct <- extract_gene_overlaps(dcopa_shared_overlaps[[c(1)]], dcopa_shared_overlaps[[c(2)]], dir, mods, mods_case)[[1]] 
    
    names(dcopa_allct) <- dplyr::recode(names(dcopa_allct), !!!CT_RENAME)

    m <- make_comb_mat(dcopa_allct)
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

    if(dir==1){
      dirlab <- "higher"
    } else if (dir==-1){
      dirlab <- "lower"
    }

    svglite::svglite(file.path(save_dir,                     
    paste0("/panel_E_F_", reg, 
           "_", dirlab, 
           ".svg")), width = 5.5, height = 3.5)                                     
    draw(p)         
    dev.off()
  }
#}

#######
# Summary tables (panels E-F)
#######

# Load genelists
dfc_list <- fread(data.table = F, file = file.path(save_dir, "panel_C_D_DFC_dcopa_genelist.csv"))

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

####### DFC

g_overlaps <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/dfc_overlaps.csv"), data.table = F)

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

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/scz_db.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/scz_ref_index.R"))
scz_hits    <- fread(file.path(save_dir, "scz_db_summary_table_dfc.csv"))
pubmed_hits <- setNames(as.integer(scz_hits$Pubmed_total_hits), scz_hits$Gene)

# ── SCZ gene database ───────────────────────────────────────────────────────
library(gt)

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

  # SCZ genes: top 10 by PubMed hits, placeholder reference superscripts
  genes_upper <- toupper(genes)
  hits        <- genes[genes_upper %in% names(scz_db_dfc)]
  hits_up     <- genes_upper[genes_upper %in% names(scz_db_dfc)]

  if (length(hits) > 0) {
    hit_counts <- pubmed_hits[hits_up]
    hit_counts[is.na(hit_counts)] <- 0L
    ord    <- order(-hit_counts)[seq_len(min(10L, length(hits)))]
    hits_s <- hits[ord]
    if (nchar(HIGHLIGHT_GENE) > 0)
      hits_s <- ifelse(toupper(hits_s) == toupper(HIGHLIGHT_GENE),
                       paste0('<span style="color:red">', hits_s, "</span>"), hits_s)
    refs_s <- sapply(seq_along(hits_up[ord]), function(i) {
      r <- scz_ref_index[hits_up[ord][i]]
      if (is.na(r)) as.character(i) else as.character(r)
    })
    scz_str <- paste0(hits_s, "<sup>", refs_s, "</sup>", collapse = ", ")
  } else {
    scz_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    SCZ_genes    = scz_str,
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
    SCZ_genes    = "SCZ-associated genes"
  ) |>
  fmt_markdown(columns = SCZ_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    SCZ_genes    ~ px(262 * GT_SCALE)
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
  hits        <- genes[genes_upper %in% names(scz_db_dfc)]
  hits_up     <- genes_upper[genes_upper %in% names(scz_db_dfc)]

  if (length(hits) > 0) {
    hit_counts <- pubmed_hits[hits_up]
    hit_counts[is.na(hit_counts)] <- 0L
    ord    <- order(-hit_counts)[seq_len(min(10L, length(hits)))]
    hits_s <- hits[ord]
    if (nchar(HIGHLIGHT_GENE) > 0)
      hits_s <- ifelse(toupper(hits_s) == toupper(HIGHLIGHT_GENE),
                       paste0('<span style="color:red">', hits_s, "</span>"), hits_s)
    refs_s <- sapply(seq_along(hits_up[ord]), function(i) {
      r <- scz_ref_index[hits_up[ord][i]]
      if (is.na(r)) as.character(i) else as.character(r)
    })
    scz_str <- paste0(hits_s, "<sup>", refs_s, "</sup>", collapse = ", ")
  } else {
    scz_str <- "—"
  }

  data.frame(
    Celltype     = ct,
    N_genes      = n_genes,
    Top_genesets = top_str,
    SCZ_genes    = scz_str,
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
    SCZ_genes    = "SCZ-associated genes"
  ) |>
  fmt_markdown(columns = SCZ_genes) |>
  fmt(columns = Top_genesets, fns = function(x) x) |>
  cols_width(
    Celltype     ~ px(60  * GT_SCALE),
    N_genes      ~ px(45  * GT_SCALE),
    Top_genesets ~ px(450 * GT_SCALE),
    SCZ_genes    ~ px(262 * GT_SCALE)
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


#####
# Genes that show up in VIP and SST:
# c("DLG2", "SCN2A", "TCF4", "PIK3CB", "PIK3R1")
# These genes show up in:
# 67 106 120 245 (1158 index)
# 61  99 113 229 (1016 index)
# Candidates: (gene shows up in top 10 module genes):
# 67 (SST, has good enrichments), 120 (VIP, has two of the genes), TCF4 (SST) (1158 index)

# Do these genes show up?
# c("SETD1A", "GRIN2A", "TRIO", "CACNA1G", "SP4", "RB1CC1", "CUL1", "XPO7", "HERC1", "SRRM2")
# RB1CC1 shows up in dfc_overlaps
bla <- dfc_list |> dplyr::filter(Gene %in% c("SETD1A", "GRIN2A", "TRIO", "CACNA1G", "SP4", "RB1CC1", "CUL1", "XPO7", "HERC1", "SRRM2")) 
# # These genes show up for CTRL mods or SCZ mods:
#                              Comparison Dataset Disease Region  Module type
# 1  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 2  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 3  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 4  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 5  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 6  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 7  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 8  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 9  CTRL modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 10  SCZ modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 11  SCZ modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
# 12  SCZ modules | CMC + SZBDMultiseq SN     CMC     SCZ    DFC Bulk megaset
#               Direction   Celltype   Gene
# 1  Lower in more severe      L6 IT RB1CC1
# 2  Lower in more severe        Sst   XPO7
# 3  Lower in more severe      L5 IT   CUL1
# 4  Lower in more severe      L5 IT RB1CC1
# 5  Lower in more severe    L5/6 NP   CUL1
# 6  Lower in more severe L6 IT Car3   CUL1
# 7  Lower in more severe L6 IT Car3 RB1CC1
# 8  Lower in more severe        L6b   CUL1
# 9  Lower in more severe        L6b RB1CC1
# 10 Lower in more severe      L6 IT   CUL1
# 11 Lower in more severe        L6b   TRIO
# 12 Lower in more severe L6 IT Car3 RB1CC1

# What about in conVsAD data (DFC)?
#                        Comparison Dataset Disease Region  Module type
# 1 CTRL modules | Gabitto + Liu SN Gabitto      AD    MTG Bulk megaset
# 2 CTRL modules | Gabitto + Liu SN Gabitto      AD    MTG Bulk megaset
# 3 CTRL modules | Gabitto + Liu SN Gabitto      AD    MTG Bulk megaset
# 4 CTRL modules | Gabitto + Liu SN Gabitto      AD    MTG Bulk megaset
# 5   AD modules | Gabitto + Liu SN Gabitto      AD    MTG Bulk megaset
# 6   AD modules | Gabitto + Liu SN Gabitto      AD    MTG Bulk megaset
#              Direction Celltype   Gene
# 1 Lower in more severe    L4 IT   TRIO
# 2 Lower in more severe    L4 IT RB1CC1
# 3 Lower in more severe    Lamp5   XPO7
# 4 Lower in more severe  L5/6 NP   CUL1
# 5 Lower in more severe    Lamp5    SP4
# 6 Lower in more severe    Pvalb    SP4


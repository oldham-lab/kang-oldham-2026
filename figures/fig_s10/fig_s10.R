# fig_s10
# ------------------------------------------------------------------------------
# Supplementary counterpart to Fig. 7 panels a-b.
#
# Fig. 7 (v8) shows, for MTG and DFC, the number of significant dCoPA modules with
# *lower* expression in AD (Direction == -1; "panel_A" dot plots). Those down-in-AD
# plots are the "phenotype" the main figure highlights. The mirror-image *higher*-in-AD
# dot plots (Direction == 1; what fig_7/v7.1 called "panel_B") were dropped from the
# main figure. This script regenerates those up-in-AD dot plots as a standalone
# supplementary figure (fig_s10, panels a = MTG, b = DFC + shared legend).
#
# The upstream data loading / dCoPA processing is taken verbatim from
# fig_7/v8/fig7_v7.1.R so the numbers are identical to the main figure; only the
# plotted Direction differs (== 1 here vs == -1 there). Per-figure decision: the
# up-in-AD data is sparse (mostly a few non-neuronal cell types), so the r-bracket
# that v8 added to the down-in-AD panels is intentionally NOT drawn here (it would
# rest on too few complete cell-type pairs to be meaningful) -- this keeps panel B
# faithful to its original fig_7/v7.1 form.
#
# Side-effect writes present in fig7_v7.1.R (shiny dot/svg-map tables, the shared
# panel_B genelist, the GSEA/summary-table panels C-G) are intentionally omitted:
# fig_s10 only re-renders the dot plots and must not overwrite those canonical outputs.

library(tidyverse)
library(qs)
library(data.table)
library(showtext)
showtext_auto()

save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s10")
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


################
# Panel a + b (up-in-AD dCoPA module counts, MTG + DFC)
###############

find_output_overlap <- function(dat1, dat2){
  dat_out <- dplyr::inner_join(dat1, dat2, by = dplyr::join_by(mod, Direction, Celltype, Consistency)) |>
    dplyr::arrange(mod)
  return(dat_out)
}

# Calculate overlaps (Gabitto + Liu, per region / module set)
dcopa_shared_overlaps <- list(
  find_output_overlap(dcopa_shared[[1]], dcopa_shared[[2]]), # Gabitto + Liu DFC
  find_output_overlap(dcopa_shared[[3]], dcopa_shared[[4]]), # Gabitto + Liu MTG
  find_output_overlap(dcopa_shared[[5]], dcopa_shared[[6]]), # Gabitto + Liu DFC (AD modules)
  find_output_overlap(dcopa_shared[[7]], dcopa_shared[[8]])  # Gabitto + Liu MTG (AD modules)
)

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

# Process bulk megaset outputs
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

# Set size range for dot legend (shared across MTG + DFC, both directions -> identical to Fig. 7)
limit_vec = c(1, max(c(plotlist[[1]]$num_sig, plotlist[[2]]$num_sig), na.rm = T))

# Plot (up-in-AD == Direction 1; 2 plots total, one per region)
dvec <- c("MTG", "DFC") # save suffix
cc_colors <- RColorBrewer::brewer.pal(3, "Set1")

for(i in 1:2){
  phi <- plotlist[[i]] |>
      filter(Direction == 1) |>
      ggplot(aes(x = Celltype, y = comp, size = num_sig, fill = Class)) +
        theme_minimal() +
        geom_point(color = "black", pch = 21) +
        theme(text = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, margin = margin(-0.1, 0, 0, 0, "cm")),
              axis.text.y = element_text(size = 7, colour = "black"),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.box = "vertical",
              legend.spacing.y = unit(4, "mm"),
              legend.title = element_blank(),
              legend.margin = margin(-0.5, 0, 0, 0, "cm"),
              legend.text = element_text(size = 7, margin = margin(0, 0, 0, -0.02, "cm")),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        labs(x = "", y = "") +
        scale_fill_manual(values = c("Non-neuronal" = cc_colors[1], "Glutamatergic" = cc_colors[3], "GABAergic" = cc_colors[2])) +
        # Size scale is kept identical to Fig. 7's down-in-AD panels (limits = limit_vec)
        # so a dot of a given size means the same module count in both figures (up-in-AD
        # effects are honestly much smaller). The shared max (158) is driven by the
        # down-in-AD counts, so its default breaks (100) exceed every up-in-AD dot
        # (max 34); set explicit breaks within this figure's range for a usable legend.
        scale_size_continuous(limits = limit_vec,
                              breaks = c(10, 20, 30))
  ggsave(phi, file = file.path(save_dir, paste0("panel_", dvec[i], ".pdf")), height = 2.5, width = 5)
  ggsave(phi, file = file.path(save_dir, paste0("panel_", dvec[i], ".svg")), height = 2.5, width = 5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_", dvec[i], "_nolegend.pdf")), height = 2.3, width = 5.5)
  ggsave(phi + theme(legend.position = "none"), file = file.path(save_dir, paste0("panel_", dvec[i], "_nolegend.svg")), height = 2.3, width = 5.5)
}

# Save legend as standalone SVG (shared by both panels; matches the panel dot scale)
legend_grob <- cowplot::get_legend(
  phi
)
p_legend_out <- ggpubr::as_ggplot(legend_grob) +
  theme(plot.margin = margin(0, 0, 0, 0))
ggsave(p_legend_out, file = file.path(save_dir, "panel_legend.svg"),
        height = 1.4, width = 3, bg = "transparent")

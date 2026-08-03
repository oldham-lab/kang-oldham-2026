# Project all genes 

library(qs)
library(data.table)
library(tidyverse)
library(showtext)
showtext_auto()
options(bitmapType = 'cairo')

version_number  <- "v2"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/"), version_number)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

# Load bulk expression data for genes and save as txt for input into python script
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
fwrite(bulk_expr[2], file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s6/genes.txt"))

# Run python script
# /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_s7/calc_mean_se_per_celltype.py

# Load output files
paths <- c(
  file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v1/expression_DFC_mean_se.csv"),
  file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v1/expression_MTG_mean_se.csv"),
  file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2/SEAAD_A9_RNAseq_final-nuclei.2024-02-13_mean_se.csv"),
  file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_mean_se.csv"),
  file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2/MIT_AD_PFC_mean_se.csv"),
  file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2/MIT_AD_MTC_mean_se.csv"),
  file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v2/morabito_2021_mean_se.csv")
)
path_names <- c(
  "Jorstad et al. 2023 (DFC)",
  "Jorstad et al. 2023 (MTG)",
  "Gabitto et al. 2024 (control DFC)",
  "Gabitto et al. 2024 (control MTG)",
  "Liu et al. 2025 (control DFC)",
  "Liu et al. 2025 (control MTG)",
  "Morabito et al. 2021 (control DFC)"
)
mats <- mapply(\(x, y){
  fread(x, data.table = F) |>
    mutate(dataset = y)
}, paths, path_names, SIMPLIFY = F) |>
  do.call(what = "rbind")

# Change celltype names to match
mats$CellType[mats$CellType == "Astro"] <- "Astrocyte"
mats$CellType[mats$CellType == "Oligo"] <- "Oligodendrocyte"
mats$CellType[mats$CellType == "Micro/PVM"] <- "Microglia-PVM"
mats$CellType[mats$CellType == "Endo"] <- "Endothelial"


allcts <- c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
          "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip",
          "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")

allcts_cap <- c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
        "Chandelier", "LAMP5", "LAMP5 LHX6", "PAX6", "PVALB", "SNCG", "SST", "SST CHODL", "VIP",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")

# Capitalize labels
mats$CellType_cap <- allcts_cap[match(mats$CellType, allcts)]

# Cell-type fill palette, matched to the fig_4 projection barplots (fig_4_v6.R).
# Built in fig_4's cell-type order and named by cell type so each subclass gets the
# identical color in both figures regardless of factor order here.
proj_pal_order <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
                    "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5", "PAX6", "SNCG", "VIP",
                    "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
proj_pal <- setNames(
  colorRampPalette(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                     "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"))(length(proj_pal_order)),
  proj_pal_order
)

# Plot indices
  # plotdf <- data.frame("ct" = rep(allcts, 3),
  #                      "dataset" = c(rep("Jorstad et al. 2023 (PFC)", n), 
  #                                    rep("Gabitto et al. 2024 (control PFC)", n), 
  #                                    rep("Morabito et al. 2021 (control PFC)", n)),
  #                      "ind" = c(unlist(lein[i, ]), 
  #                                unlist(SEAcon[i, ]), 
  #                                unlist(mora[i, ])),
  #                      "ind_se" = c(unlist(lein_se[i, ]), 
  #                                   unlist(SEAcon_se[i, ]), 
  #                                   unlist(mora_se[i, ]))) |>
  #   dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
  #                 ct = factor(ct, levels = allcts))
 
p <- mats |>
  mutate(dataset = factor(dataset, levels = unique(dataset)),
         CellType_cap = factor(CellType_cap, levels = allcts_cap)) |>
  ggplot(aes(x = CellType_cap, y = Mean, fill = CellType_cap)) +
    theme_classic() +
    geom_col(position = position_dodge(), alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean - 2 * SE,
                      ymax = Mean + 2 * SE),
                      width = 0.2,
                      linewidth = 0.3,
                      position = position_dodge(0.5)) +
    scale_fill_manual(values = proj_pal) +
    theme(text = element_text(family = "sans", color = "black", size = 12),
          legend.position = "none", 
          axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(color = "black"),
          strip.background = element_rect(fill = "white")) +
    facet_wrap(~dataset, ncol = 1, nrow = length(path_names)#, 
               #scales = "free_y"
               ) +
    labs(y = "Mean expression (log UMI counts + 1)", x = "") #+
   # scale_y_continuous(breaks = c(0, 0.5, 1))

ggsave(p, file = file.path(save_dir, paste0("panel_1.pdf")), height = 7.5, width = 5.5)
ggsave(p, file = file.path(save_dir, paste0("panel_1.svg")), height = 7.5, width = 5.5)



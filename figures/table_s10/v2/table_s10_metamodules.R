# Table S10 (meta-module tabs)
# Meta-module eigengene values underlying the heatmaps in fig_5 panel C (DFC) and
# fig_s8 panel C (MTG). Each value is a meta-module's eigengene (PC1 of the REIs of
# its constituent modules) in one cell type. Four tabs = 2 datasets (Jorstad,
# Gabitto) x 2 regions (DFC, MTG). Data prep mirrors fig_5/v5/fig_5_v5.R and
# fig_s8/v2/fig_s8_v2.R: load Module_eigengenes.csv, order Gabitto's meta-modules to match
# Jorstad's (jor_order), number meta-modules 1..N as in the figures, label cell
# types via class_info. Output is transposed (rows = meta-module #, cols = cell type)
# and column-ordered to match the per-module REI tabs.

library(data.table)
library(tidyverse)

version <- "v2"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s10/"), version)

# --- Cell type labels for eigengene rows (same construction as the figures) ---
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
  arrange(Subclass_fixed)

allcts <- c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
        "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
allcts_cap <- c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
        "Chandelier", "LAMP5", "LAMP5 LHX6", "PAX6", "PVALB", "SNCG", "SST", "SST CHODL", "VIP",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
ct_labels <- allcts_cap[match(class_info[,1], allcts)]   # cell-type label per eigengene row

# Cell-type column order to match the per-module REI tabs (table_s10.R)
ct_order <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
          "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5",  "PAX6",  "SNCG",  "VIP",
          "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")

# --- Eigengene file paths (consensus-min, no-merge; the networks plotted in the figures) ---
# HGNC v2: the FindModules meta-module networks were regenerated (fig_5/v5 DFC,
# fig_s8/v2 MTG) at shifted signum + node-count suffixes, so resolve the minSize3
# "largest network" (lowest signum) by glob rather than hard-coding, exactly as
# fig_5/v5/fig_5_v5.R and fig_s8/v2/fig_s8_v2.R do.
f5 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v5")
f8 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s8/v2")
find_net <- function(load_dir, prefix){
  hits <- Sys.glob(file.path(load_dir, paste0(prefix, "_consensusMin_noMerge_Modules"),
                             "Pearson-no_TO_signum*_minSize3_merge_ME_1_*"))
  if(length(hits) == 0) stop("No minSize3 ", prefix, " network in ", load_dir)
  signum <- as.numeric(sub(".*signum([0-9.]+)_minSize3.*", "\\1", basename(hits)))
  hits[which.min(signum)]   # largest network = lowest signum
}
paths <- list(
  jorstad_dfc = file.path(find_net(f5, "Jorstad_DFC"), "Module_eigengenes.csv"),
  gabitto_dfc = file.path(find_net(f5, "Gabitto_DFC"), "Module_eigengenes.csv"),
  jorstad_mtg = file.path(find_net(f8, "Jorstad_MTG"), "Module_eigengenes.csv"),
  gabitto_mtg = file.path(find_net(f8, "Gabitto_MTG"), "Module_eigengenes.csv")
)

# rows = cell types, cols = meta-modules. Gabitto reordered to Jorstad's meta-module
# order; meta-modules numbered 1..N. Returns a df: Meta-module # + 24 cell type cols.
build_metamod <- function(path, jor_order = NULL){
  me <- fread(data.table = F, file = path) |> column_to_rownames("Sample")
  if(!is.null(jor_order)) me <- me[, jor_order, drop = FALSE]  # match Jorstad numbering
  ord <- colnames(me)                # meta-module order (= 1..N labels in the figure)
  rownames(me) <- ct_labels          # label rows by cell type
  out <- as.data.frame(t(me))        # rows = meta-modules, cols = cell types
  out <- out[, ct_order, drop = FALSE]
  out <- cbind("Meta-module #" = seq_len(nrow(out)), out)
  rownames(out) <- NULL
  list(df = out, order = ord)
}

# DFC: Jorstad defines the meta-module order; Gabitto follows it
jdfc <- build_metamod(paths$jorstad_dfc)
gdfc <- build_metamod(paths$gabitto_dfc, jor_order = jdfc$order)
# MTG: Jorstad defines the meta-module order; Gabitto follows it
jmtg <- build_metamod(paths$jorstad_mtg)
gmtg <- build_metamod(paths$gabitto_mtg, jor_order = jmtg$order)

fwrite(jdfc$df, file.path(save_dir, "MetaMod_Eigenmod_Jorstad_DFC.csv"))
fwrite(gdfc$df, file.path(save_dir, "MetaMod_Eigenmod_Gabitto_DFC.csv"))
fwrite(jmtg$df, file.path(save_dir, "MetaMod_Eigenmod_Jorstad_MTG.csv"))
fwrite(gmtg$df, file.path(save_dir, "MetaMod_Eigenmod_Gabitto_MTG.csv"))

message("Meta-module tabs: DFC = ", nrow(jdfc$df), " meta-modules; MTG = ", nrow(jmtg$df), " meta-modules.")

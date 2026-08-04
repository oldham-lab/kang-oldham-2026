# panel_S.R
# Standalone fast regenerator for Fig. 4 panel S (projection-correlation heatmap).
# Recomputes only `cormeans` from the projection-index CSVs and renders panel_S.svg
# through the shared make_panel_S() in panel_S_plot.R -- avoiding the full ~10 min
# fig_4_v6.R pipeline (which also recomputes the per-module variance-explained loop
# and the full gene x gene correlation matrix that panel S does not need).

library(qs)
library(data.table)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(ggdendro)
library(showtext)
showtext_auto()
options(bitmapType = 'cairo')

v6_dir   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/v6")
base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")

# --- Module set: topmodposbc modules with >= 2 significant bulk-platform cors ---
module_output_dir <- file.path(base_dir, "SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table = F, file = file.path(module_output_dir, "kme_tables", "topmodposbc_table.csv"))
if (sum(duplicated(datkme[, 2])) > 0) { datkme[, 2] <- make.unique(datkme[, 2]) }
mods <- tapply(datkme[, 2], datkme[, 3], list)
modulelengths <- unlist(lapply(mods, length))
these_mods <- as.numeric(names(mods)[which(modulelengths > filter_under)])
sigcount_bonf <- fread(data.table = F,
  file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods_final <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]

# --- Projection indices (REI), per dataset ------------------------------------
d <- "log_REI"
lein   <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", d, "/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |> dplyr::select(!module)
mitcon <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/mod_means/"), d, "mod_means_Con_bulk_megaset.csv"))
# Gabitto (SEA-AD): HGNC-harmonized Python output (same format as MIT/Liu; no `module` col)
SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/mod_means/", d, "/mod_means_Con_bulk_megaset.csv"))
mora   <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", d, "/indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv")) |> dplyr::select(!module)

# --- Align cell types across datasets (mirrors fig_4_v6.R) --------------------
colnames(lein)[c(1, 3, 15, 16)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
colnames(mora)[c(1, 3, 14, 15)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
allcts <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
        "Chandelier", "Pvalb", "Sst", "Sst Chodl", "Lamp5 Lhx6", "Lamp5", "Pax6", "Sncg", "Vip",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
for (nm in c("lein", "mitcon", "SEAcon", "mora")) {
  df <- get(nm)
  df[setdiff(allcts, colnames(df))] <- NA
  df <- df[, match(allcts, colnames(df))]
  assign(nm, df)
}

# --- Mean pairwise correlation matrix across datasets -------------------------
corlist <- lapply(these_mods_final, \(i){
  tempdf <- data.frame("lein" = unlist(lein[i, ]),
                       "SEA"  = unlist(SEAcon[i, ]),
                       "mit"  = unlist(mitcon[i, ]),
                       "mora" = unlist(mora[i, ]))
  cor(tempdf, use = "pairwise.complete.obs")
})
cormeans <- Reduce("+", lapply(corlist, function(m) { m[is.na(m)] <- 0; m })) /
            Reduce("+", lapply(corlist, function(m) !is.na(m)))
colnames(cormeans) <- rownames(cormeans) <- c("Jorstad 2023", "Gabitto 2024", "Liu 2025", "Morabito 2021")
print(cormeans)

# --- Render panel S -----------------------------------------------------------
source(file.path(v6_dir, "panel_S_plot.R"))
make_panel_S(cormeans, file.path(v6_dir, "panel_S.svg"))
cat("wrote", file.path(v6_dir, "panel_S.svg"), "\n")

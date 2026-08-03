# panel_mp.R
# Standalone fast regenerator for Fig. 4 panels m-p (projection barplots) WITH the
# nested dataset-correlation brackets baked in (panel_X_4_with_brackets.svg). Loads
# only the projection-index means/SEs and renders via the shared builders in
# panel_mp_plot.R -- avoiding the full ~10 min fig_4_v6.R pipeline.

library(qs)
library(data.table)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(showtext)
showtext_auto()
options(bitmapType = 'cairo')

v6_dir   <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/v6")
base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")
mod_id   <- c(112, 669, 681, 1007)   # modules shown in panels m-p (1158 indices)

# --- Projection indices (REI): means + SEs, per dataset ------------------------
d <- "log_REI"
lein   <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", d, "/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |> dplyr::select(!module)
lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", d, "/indices_se_topmodposbc.csv"))
lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]
mitcon <- fread(data.table = F, file = file.path("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/mod_means/", d, "mod_means_Con_bulk_megaset.csv"))
mitse  <- fread(data.table = F, file = file.path("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/se/", d, "se_Con_bulk_megaset.csv"))
# Gabitto (SEA-AD): HGNC-harmonized Python output (same format as MIT/Liu; no `module` col)
SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/mod_means/", d, "/mod_means_Con_bulk_megaset.csv"))
SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/se/", d, "/se_Con_bulk_megaset.csv"))
SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]
mora   <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", d, "/indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv")) |> dplyr::select(!module)
mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", d, "/indices_se_SubclassContopmodposbc.csv"))
mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]

# --- Align cell types across datasets (mirrors fig_4_v6.R, including its SE quirks) -
colnames(lein)[c(1, 3, 15, 16)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
colnames(mora)[c(1, 3, 14, 15)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
colnames(lein_se) <- colnames(lein); colnames(SEAcon_se) <- colnames(SEAcon); colnames(mora_se) <- colnames(mora)
allcts <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
        "Chandelier", "Pvalb", "Sst", "Sst Chodl", "Lamp5 Lhx6", "Lamp5", "Pax6", "Sncg", "Vip",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
allcts_cap <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
        "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5", "PAX6", "SNCG", "VIP",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
fill_match <- function(df) { df[setdiff(allcts, colnames(df))] <- NA; df[, match(allcts, colnames(df))] }
lein <- fill_match(lein); lein_se <- fill_match(lein_se)
mitcon <- fill_match(mitcon); mitse <- fill_match(mitse)
SEAcon <- fill_match(SEAcon); SEAcon_se <- fill_match(SEAcon_se)
mora <- fill_match(mora); mora_se <- fill_match(mora_se)
for (nm in c("lein", "mitcon", "SEAcon", "mora", "lein_se", "mitse", "SEAcon_se", "mora_se")) {
  df <- get(nm); colnames(df) <- allcts_cap; assign(nm, df)
}

# --- Projection barplot palette (matched to fig_3_v4.R, mapped by cell type) ----
proj_pal <- setNames(
  colorRampPalette(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                     "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"))(length(allcts_cap)),
  allcts_cap)

# --- Build each panel + brackets ---------------------------------------------
source(file.path(v6_dir, "panel_mp_plot.R"))
dfs <- list(lein, SEAcon, mitcon, mora)
ses <- list(lein_se, SEAcon_se, mitse, mora_se)
for (j in seq_along(mod_id)) {
  out <- make_proj_panel(mod_id[j], dfs, ses, allcts_cap, proj_pal, "Relative expression index")
  cat(mod_id[j], "pairwise cors:\n"); print(round(out$cormat, 3))
  ggsave(add_corr_brackets(out$plot, out$cormat),
         file = file.path(v6_dir, paste0("panel_", j, "_4_with_brackets.svg")),
         width = 6.6, height = 5)
}
cat("wrote panel_{1..4}_4_with_brackets.svg\n")

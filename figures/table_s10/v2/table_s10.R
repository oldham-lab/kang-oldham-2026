# Table S10
# Relative expression index (REI) values underlying the fig_3 panel F REI barplots,
# for ALL modules (not just the single module shown in the figure).
# Four datasets (Jorstad, Gabitto, Liu, Morabito) + an element-wise average.
#
# Data prep mirrors fig_3/v4/fig_3_v4.R exactly for the REI index (index_dirs "log_REI"),
# so values here match the bar heights in panel F. Modules are the size-filtered set
# (size > 3) with the fig_3 significance (sigmod) filter applied -> 1016 modules,
# renumbered 1..N, matching the figure's module numbering.

library(data.table)
library(tidyverse)

version <- "v2"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s10/"), version)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

# --- Module set + numbering (same as fig_3_v4.R) ---
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])

# Significance (sigmod) filter from fig_3_v4.R: drop modules whose bulk-correlation
# significance count is < 2. Takes the set from 1023 -> 1016 (figure module stays at 354).
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]

# --- Load REI indices for the four datasets (REI = index_dirs "log_REI") ---
base_dir   <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")
index_dir  <- "log_REI"

lein <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dir, "/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |>
  dplyr::select(!module)
mitcon <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/mod_means/"), index_dir, "mod_means_Con_bulk_megaset.csv"))
# Gabitto (SEA-AD) now routes through the HGNC-harmonized Python SEA-AD output
# (SEAAD2024_full_python_output), matching fig_4/v6 + fig_5/v5. Supersedes the pre-fix
# R sn_proj_indices output; the Python format has no `module` column (no select needed).
SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/mod_means/", index_dir, "/mod_means_Con_bulk_megaset.csv"))
mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dir, "/indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv")) |>
  dplyr::select(!module)

# Change celltypes to match (same column fixes as fig_3_v4.R)
colnames(lein)[c(1, 3, 15, 16)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
colnames(mora)[c(1, 3, 14, 15)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")

# Cell type order / labels used in the figure
allcts <- c("L2/3 IT", "L4 IT","L5 IT",  "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
          "Chandelier", "Pvalb", "Sst", "Sst Chodl", "Lamp5 Lhx6", "Lamp5",  "Pax6",  "Sncg",  "Vip",
          "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")

allcts_cap <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
          "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5",  "PAX6",  "SNCG",  "VIP",
          "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")

# Fill out / order cell types, then relabel to the figure's capitalized labels
fill_order <- function(d){
  d[setdiff(allcts, colnames(d))] <- NA
  d <- d[, match(allcts, colnames(d))]
  colnames(d) <- allcts_cap
  d
}
lein   <- fill_order(lein)
mitcon <- fill_order(mitcon)
SEAcon <- fill_order(SEAcon)
mora   <- fill_order(mora)

# --- Subset to the filtered modules and renumber 1..N (figure numbering) ---
build_tab <- function(d){
  out <- d[these_mods, , drop = FALSE]
  out <- cbind("Module #" = seq_along(these_mods), out)
  rownames(out) <- NULL
  out
}
jorstad  <- build_tab(lein)    # Jorstad et al. 2023
gabitto  <- build_tab(SEAcon)  # Gabitto et al. 2024
liu      <- build_tab(mitcon)  # Liu et al. 2025
morabito <- build_tab(mora)    # Morabito et al. 2021

# --- Average at each cell across the four datasets (na.rm: average available datasets) ---
ct_cols <- allcts_cap
stack <- abind::abind(
  as.matrix(jorstad[ , ct_cols]),
  as.matrix(gabitto[ , ct_cols]),
  as.matrix(liu[ , ct_cols]),
  as.matrix(morabito[ , ct_cols]),
  along = 3
)
avg_mat <- apply(stack, c(1, 2), function(x) if(all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
average <- cbind(data.frame("Module #" = seq_along(these_mods), check.names = FALSE),
                 as.data.frame(avg_mat))

# --- Write CSVs (one per tab) ---
fwrite(jorstad,  file.path(save_dir, "REI_Jorstad_2023.csv"))
fwrite(gabitto,  file.path(save_dir, "REI_Gabitto_2024.csv"))
fwrite(liu,      file.path(save_dir, "REI_Liu_2025.csv"))
fwrite(morabito, file.path(save_dir, "REI_Morabito_2021.csv"))
fwrite(average,  file.path(save_dir, "REI_Average.csv"))

message("Wrote ", length(these_mods), " modules x ", length(ct_cols), " cell types for 4 datasets + average.")

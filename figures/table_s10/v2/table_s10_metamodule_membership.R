# Table S10 (meta-module membership tabs)
# Module-membership of meta-modules: the meta-module analogue of Table S9's
# gene-membership of modules. Each meta-module is a FindModules cluster of modules
# (each module a "node") built on module-module REI correlations; data prep mirrors
# fig_5/v5/fig_5_v5.R (DFC) and fig_s8/v2/fig_s8_v2.R (MTG). Four tabs = 2 datasets
# (Jorstad, Gabitto) x 2 regions (DFC, MTG), matching the four Eigenmod tabs.
#
# Each tab is keyed to the figure module numbering (1..1016, the size>3 + significance
# set from table_s10.R / fig_3_v4.R, matching the REI tabs), and gives each module's
# seed and topmodposbc meta-module assignment. Meta-modules are numbered as in the
# Eigenmod tabs (Figure 5 / Figure S8 panel C; Gabitto ordered to match Jorstad).
# Modules with no assignment are left blank, exactly as unassigned genes are in Table S9.

library(data.table)
library(tidyverse)

version <- "v2"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s10/"), version)

# ============================================================
# Figure module numbering (size>3 + significance) -> 1016 modules, numbered 1..N.
# Identical to table_s10.R / fig_3_v4.R so Module # matches the REI tabs.
# ============================================================
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(data.table = F, file = file.path(module_output_dir, "kme_tables", "topmodposbc_table.csv"))
if(sum(duplicated(datkme[, 2])) > 0){ datkme[, 2] <- make.unique(datkme[, 2]) }
mods <- tapply(datkme[, 2], datkme[, 3], list)
modulelengths <- unlist(lapply(mods, length))
size_pass_ids <- as.numeric(names(mods)[which(modulelengths > 3)])   # 1023 module ids (tapply order)
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
final_ids <- size_pass_ids[!size_pass_ids %in% which(sigcount_bonf$vals < 2)]   # 1016 figure modules

# ============================================================
# Mod# (meta-module network node label) -> module id, per region.
#   DFC (fig_5_v5.R): nodes = size>3 + significance set, Mod# = 1..length(final).
#     v5 now applies the significance filter (fig_5 v2/v3 did not), so DFC Mod# matches
#     the figure Module # directly -- same as MTG. (The old 1023-space crosswalk is gone.)
#   MTG (fig_s8/v2/fig_s8_v2.R): nodes = size>3 + significance set, Mod# = 1..length(final)
# Zero-variance nodes were dropped before the network, so those Mod#s never appear in
# the kME table -- the join below simply leaves those modules blank.
# ============================================================
dfc_modk_to_id <- setNames(final_ids, paste0("Mod", seq_along(final_ids)))
mtg_modk_to_id <- setNames(final_ids, paste0("Mod", seq_along(final_ids)))

# ============================================================
# Meta-module color -> number, matching the Eigenmod tabs (build_metamod in
# table_s10_metamodules.R): Jorstad column order defines the numbering; Gabitto is
# reordered to Jorstad so a given color carries the same number in both datasets.
# ============================================================
# HGNC v2: resolve the regenerated minSize3 "largest network" (lowest signum) by glob
# rather than hard-coding the signum/node-count suffix (fig_5/v5 DFC, fig_s8/v2 MTG);
# mirrors fig_5/v5/fig_5_v5.R + table_s10_metamodules.R.
f5 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/v5")
f8 <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s8/v2")
find_net <- function(load_dir, prefix){
  hits <- Sys.glob(file.path(load_dir, paste0(prefix, "_consensusMin_noMerge_Modules"),
                             "Pearson-no_TO_signum*_minSize3_merge_ME_1_*"))
  if(length(hits) == 0) stop("No minSize3 ", prefix, " network in ", load_dir)
  signum <- as.numeric(sub(".*signum([0-9.]+)_minSize3.*", "\\1", basename(hits)))
  hits[which.min(signum)]
}
paths <- list(
  jorstad_dfc = find_net(f5, "Jorstad_DFC"),
  gabitto_dfc = find_net(f5, "Gabitto_DFC"),
  jorstad_mtg = find_net(f8, "Jorstad_MTG"),
  gabitto_mtg = find_net(f8, "Gabitto_MTG")
)

# color -> meta-module # from an eigengene file (optionally reordered to Jorstad)
color_to_metanum <- function(net_dir, jor_order = NULL){
  cols <- colnames(fread(data.table = F, file.path(net_dir, "Module_eigengenes.csv")))[-1]
  if(!is.null(jor_order)) cols <- jor_order
  setNames(seq_along(cols), cols)
}
jor_dfc_order <- colnames(fread(data.table = F, file.path(paths$jorstad_dfc, "Module_eigengenes.csv")))[-1]
jor_mtg_order <- colnames(fread(data.table = F, file.path(paths$jorstad_mtg, "Module_eigengenes.csv")))[-1]
metanum <- list(
  jorstad_dfc = color_to_metanum(paths$jorstad_dfc),
  gabitto_dfc = color_to_metanum(paths$gabitto_dfc, jor_order = jor_dfc_order),
  jorstad_mtg = color_to_metanum(paths$jorstad_mtg),
  gabitto_mtg = color_to_metanum(paths$gabitto_mtg, jor_order = jor_mtg_order)
)

# ============================================================
# Build one membership tab: 1016 figure modules x (seed, topmodposbc) meta-module #.
# ============================================================
build_membership <- function(net_dir, modk_to_id, colormap){
  k <- fread(data.table = F, file.path(net_dir, "kME_table_.csv"))
  bc_col <- grep("^TopModPosBC", colnames(k), value = TRUE)
  stopifnot(length(bc_col) == 1)
  link <- data.frame(
    module_id = modk_to_id[k$Gene],
    seed      = colormap[k$ModSeed],
    topmodposbc = colormap[k[[bc_col]]]
  )
  out <- data.frame("Module #" = seq_along(final_ids), module_id = final_ids,
                    check.names = FALSE) |>
    left_join(link, by = "module_id") |>
    transmute(`Module #`,
              `Meta-module # (seed)` = seed,
              `Meta-module # (topmodposbc)` = topmodposbc)
  out
}

jdfc <- build_membership(paths$jorstad_dfc, dfc_modk_to_id, metanum$jorstad_dfc)
gdfc <- build_membership(paths$gabitto_dfc, dfc_modk_to_id, metanum$gabitto_dfc)
jmtg <- build_membership(paths$jorstad_mtg, mtg_modk_to_id, metanum$jorstad_mtg)
gmtg <- build_membership(paths$gabitto_mtg, mtg_modk_to_id, metanum$gabitto_mtg)

fwrite(jdfc, file.path(save_dir, "MetaMod_Membership_Jorstad_DFC.csv"))
fwrite(gdfc, file.path(save_dir, "MetaMod_Membership_Gabitto_DFC.csv"))
fwrite(jmtg, file.path(save_dir, "MetaMod_Membership_Jorstad_MTG.csv"))
fwrite(gmtg, file.path(save_dir, "MetaMod_Membership_Gabitto_MTG.csv"))

report <- function(name, df) message(
  name, ": ", nrow(df), " modules; assigned (topmodposbc) = ",
  sum(!is.na(df$`Meta-module # (topmodposbc)`)), ", seed = ",
  sum(!is.na(df$`Meta-module # (seed)`)))
report("Jorstad DFC", jdfc); report("Gabitto DFC", gdfc)
report("Jorstad MTG", jmtg); report("Gabitto MTG", gmtg)

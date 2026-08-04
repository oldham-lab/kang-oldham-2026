library(ComplexHeatmap)
library(tidyverse)
library(data.table)
library(dendextend)
library(qs)
library(showtext)
showtext_auto()

version <- "v5"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

###########################
# Panel A: REI matrix schematic
###########################
set.seed(23)

rei1 <- matrix(runif(20, 0, 1), nrow = 5, ncol = 4)
rei2 <- matrix(runif(20, 0, 1), nrow = 5, ncol = 4)

rei1c <- cor(rei1)
rei2c <- cor(rei2)
rei1m <- cor(t(rei1))
rei2m <- cor(t(rei2))

reip1 <- pmin(rei1c, rei2c)
reip2 <- pmin(rei1m, rei2m)

h_col <- circlize::colorRamp2(c(-1, 1), c("white", "grey"))
h_col1 <- circlize::colorRamp2(c(0, 1), c("white", "#420D09"))
h_col2 <- circlize::colorRamp2(c(0, 1), c("white", "#151B54"))

p <- Heatmap(as.matrix(rei1),
             name = "REI",
             col = h_col1,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Modules (REI)",
             column_title = "Celltypes",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 3)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A1.svg"), width = 2.5, height = 3)
draw(p, heatmap_legend_side = "bottom")
dev.off()
pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A1.pdf"), width = 2.5, height = 3)
draw(p, heatmap_legend_side = "bottom")
dev.off()

p <- Heatmap(as.matrix(rei2),
             name = "REI",
             col = h_col2,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Modules (REI)",
             column_title = "Celltypes",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 8)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A2.svg"), width = 2.5, height = 3)
draw(p, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A2.pdf"), width = 2.5, height = 3)
draw(p, heatmap_legend_side = "bottom")
dev.off()

p <- Heatmap(as.matrix(rei1c),
             name = "REI",
             col = h_col1,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Celltypes",
             column_title = "Celltypes",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 10)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A3.svg"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A3.pdf"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

p <- Heatmap(as.matrix(rei2c),
             name = "REI",
             col = h_col2,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Celltypes",
             column_title = "Celltypes",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 10)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A4.svg"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A4.pdf"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

p <- Heatmap(as.matrix(rei1m),
             name = "REI",
             col = h_col1,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Modules",
             column_title = "Modules",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 10)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A5.svg"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A5.pdf"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

p <- Heatmap(as.matrix(rei2m),
             name = "REI",
             col = h_col2,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Modules",
             column_title = "Modules",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 10)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A6.svg"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A6.pdf"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

p <- Heatmap(as.matrix(reip1),
             name = "REI",
             col = h_col,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Celltype",
             column_title = "Celltype",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 10)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A7.svg"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A7.pdf"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()


p <- Heatmap(as.matrix(reip2),
             name = "REI",
             col = h_col,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Module",
             column_title = "Module",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 10)),
             show_heatmap_legend = F,
             width = 5*unit(5, "mm"), 
             height = 5*unit(5, "mm"),
             border = "black",
             column_title_gp = gpar(fontsize = 11),
             row_title_gp = gpar(fontsize = 11))

svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A8.svg"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "/panel_A8.pdf"), width = 2.5, height = 2.5)
draw(p, heatmap_legend_side = "bottom")
dev.off()

###############
# Panel B: Dendrogram - clustering subclasses over all projections
###############
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

class_info[,4] <- allcts_cap[match(class_info[,1], allcts)]


# HGNC fix: Gabitto REI repointed from the stale R sn_proj_indices to the HGNC-fixed
# Python SEA-AD output (region PFC = DFC). Python file has no `module` column;
# select(class_info[,1]) keeps the 24 celltype columns by name, so no other change needed.
rei1 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/mod_means/log_REI/mod_means_Con_bulk_megaset.csv")) |>
  select(class_info[,1])

rei2 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_REI/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |>
  select(class_info[,3])


# Filter modules
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)
# v5: add the significance (sigmod) filter used in fig_3_v4.R / fig_s8, taking the
# module set from 1023 -> 1016 so DFC matches the figure module numbering (REI tabs)
# and the MTG pipeline. (v3 omitted this filter and stayed at 1023.)
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
mod_seed <- mod_seed[these_mods]
mod_bc <- mod_bc[these_mods]
rei1 <- rei1[these_mods, ]
rei2 <- rei2[these_mods, ]
colnames(rei2) <- colnames(rei1)
rownames(rei1) <- 1:nrow(rei1)
rownames(rei2) <- 1:nrow(rei2)

# Remove zero-variance projections.
# HGNC re-run: recompute the zero-variance module positions from the (fixed) data instead of
# hardcoding c(753, 972). The HGNC fix recovered the all-histone module 972 (was all-zero in
# Gabitto pre-fix), so post-fix only Jorstad's protocadherin module (753) is zero-variance.
# Network input = length(these_mods) - length(zerovar_vec).
var1 <- apply(rei1, 1, var)   # Gabitto
var2 <- apply(rei2, 1, var)   # Jorstad
zerovar_vec <- sort(union(which(var1 == 0), which(var2 == 0)))

rei1 <- rei1[-zerovar_vec, ]
rei2 <- rei2[-zerovar_vec, ]

# Calculate mod-mod and celltype-celltype matrices

ct1 <- cor(rei1)
ct2 <- cor(rei2)
m1 <- cor(t(rei1))
m2 <- cor(t(rei2))

# Calculate consensus matrices (mean)

ctc1 <- (ct1 + ct2) / 2
mc1 <- (m1 + m2) / 2
colnames(mc1) <- paste0("Mod", colnames(mc1))
rownames(mc1) <- paste0("Mod", rownames(mc1))

# Calculate consensus matrices (min)

ctm1 <- pmin(ct1, ct2)
mcm1 <- pmin(m1, m2)
colnames(mcm1) <- paste0("Mod", colnames(mcm1))
rownames(mcm1) <- paste0("Mod", rownames(mcm1))

# Create dendrograms from consensus ct matrices
clustermean <- hclust(as.dist(1 - ctc1), method="complete") |>
  as.dendrogram() 
svg(file.path(save_dir, "/panel_B_mean.svg"), width = 6, height = 4)
par(mar = c(8.1, 4.1, 4.1, 2.1))
plot(clustermean, ylab = "1 - cor"#, 
    # main = "Clustering all subclasses over all modules\nusing mean consensus matrix", 
    # cex.main = 1, 
    # font.main = 1
     )
dev.off()

clustermin <- hclust(as.dist(1 - ctm1), method="complete") |>
  as.dendrogram() 
svg(file.path(save_dir, "/panel_B_min.svg"), width = 6, height = 4)
par(mar = c(8.1, 4.1, 4.1, 2.1))
plot(clustermin, ylab = "1 - cor", 
     #main = "Clustering all subclasses over all modules\nusing min consensus matrix", 
     #cex.main = 1, 
     #font.main = 1
     )
dev.off()
 
# ##########
# # Panel C: FM on mod-mod consensus matrix
# ##########

setwd(file.path(Sys.getenv("FINDMODULES_DIR", "/home/gugene/code/git/FindModules"), "FindModules/R/"))
source("FindModules.R")
source("map_identifiers_function.R")
source("FM_helper_fxns.R")
source("FindModules.R")
source("find_seed_genes_greedy_march_megaset.R")
source("similarityType.R")
source("plotting_functions.R")
source("overlapType.R")
source("networkOutputs.R")
source("module_quant_functions.R")
source("iteration_code.R")

setwd(save_dir)

# ### Try no merge

# Gabitto DFC
# HGNC fix: repointed to the HGNC-fixed Python SEA-AD output (PFC = DFC), as for rei1 above.
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/mod_means/log_REI/mod_means_Con_bulk_megaset.csv")) |>
  select(class_info[,1])
expr <- data.frame("Gene" = paste0("Mod", 1:length(these_mods)), proj[these_mods, ]) 
expr <- expr[-zerovar_vec, ] # remove zero-variance REI indices (recomputed from fixed data)

# Gabitto DFC min
FindModules(
  projectname="Gabitto_DFC_consensusMin_noMerge",
  data_cols=expr,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(expr),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = mcm1,
  saveSimMat = FALSE,
  simType = c("Pearson"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3, 5, 8),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.97 ,.96, .95),
  minMEcorvec = c(0.95),
  merge.by = c("ME"),
  merge.param = 1,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

# # Jorstad DFC
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_REI/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |>
  select(class_info[,3])
expr <- data.frame("Gene" = paste0("Mod", 1:length(these_mods)), proj[these_mods, ])
expr <- expr[-zerovar_vec, ] # remove zero-variance REI indices (recomputed from fixed data)

FindModules(
  projectname="Jorstad_DFC_consensusMin_noMerge",
  data_cols=expr,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(expr),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = mcm1,
  saveSimMat = FALSE,
  simType = c("Pearson"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3, 5, 8),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.97 ,.96, .95),
  minMEcorvec = c(0.95),
  merge.by = c("ME"),
  merge.param = 1,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

# With no merging, module consistency is retained between FM runs on different datasets

############
# Cluster modules
############

# v5: load the regenerated DFC networks produced by the FindModules calls above (which now
# write into this v5 dir). FindModules writes one minSize3 network per signum in signumvec;
# the figure plots the "largest network" = lowest signum / most meta-modules. Resolve by glob
# + min-signum rather than hard-coding the auto-selected value. The module-count suffix is also
# wildcarded (data-driven zerovar removal changes it, e.g. 1014 -> 1015).
load_dir <- save_dir
find_net <- function(prefix){
  hits <- Sys.glob(file.path(load_dir, paste0(prefix, "_consensusMin_noMerge_Modules"),
                             "Pearson-no_TO_signum*_minSize3_merge_ME_1_*"))
  if(length(hits) == 0)
    stop("No minSize3 ", prefix, " network in ", load_dir,
         " (run the FindModules calls above first).")
  signum <- as.numeric(sub(".*signum([0-9.]+)_minSize3.*", "\\1", basename(hits)))
  hits[which.min(signum)]  # largest network (lowest signum), as in v2/v3
}
jor_net <- find_net("Jorstad_DFC")
gab_net <- find_net("Gabitto_DFC")

### Plot largest network (consensus min, Jorstad)
mod_eig <- fread(data.table = F,file = file.path(jor_net, "Module_eigengenes.csv")) |>
  column_to_rownames("Sample")
rownames(mod_eig) <- class_info[,4]
jor_order <- colnames(mod_eig)
colnames(mod_eig) <- 1:ncol(mod_eig)

metacluster_min <- hclust(as.dist(1 - cor(mod_eig)), method="complete") |>
  as.dendrogram()
qsave(metacluster_min, file = file.path(save_dir, "jorstad_consensusMin_dendro.qs"))

#mod_eig <- mod_eig[order.dendrogram(clustermin), order.dendrogram(metacluster_min)]
#fwrite(mod_eig, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version, "mod_eig.csv"), row.names = T)

kme <- fread(data.table = F, file = file.path(jor_net, "kME_table_.csv"))
mod_fdr <- tapply(kme[,1], kme[,6], list) |>
  lapply(\(x) as.numeric(gsub("Mod", "", x)))
mod_fdr_gene <- lapply(mod_fdr, \(x) mod_seed[x] |> unlist() |> unique())

# Count mods per meta-module
modcountdf <- data.frame("mod" = names(mod_fdr),
                    "FDR" = lapply(mod_fdr,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "mod_count") |>
  mutate(mod_per = mod_count/length(these_mods) * 100) |>
  arrange(mod_per)
genecountdf <- data.frame("mod" = names(mod_fdr_gene),
                          "FDR" = lapply(mod_fdr_gene,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "gene_count") |>
  mutate(gene_per = gene_count/18913 * 100) |>
  arrange(gene_per)
countdf <- full_join(modcountdf[,-c(2:3)], genecountdf[,-c(2:3)], by = join_by(mod)) |>
  as.data.frame()
countdf <- countdf[match(jor_order, countdf$mod), ]
countdf$mod <- 1:nrow(countdf)

p <- Heatmap(as.matrix(mod_eig),
        name = "Eigenmodule",
        cluster_rows = clustermin,
        cluster_columns = metacluster_min,
        top_annotation = HeatmapAnnotation("% of all mods" = anno_barplot(countdf[,2],
                                                                          axis_param = list(at = c(0, 2))), 
                                           "% of all genes" = anno_barplot(countdf[,3],
                                                                           axis_param = list(at = c(0, 1)))),
        heatmap_legend_param = list(title_gp = gpar(fontface = "plain"),
                                    title_position = "topcenter"),
        column_title_side = "bottom",
        column_title = "Meta-modules",
        row_title = "1 - cor",
        row_names_side = "left",
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 12),
        show_heatmap_legend = FALSE,
        show_column_dend = FALSE   # labeled dendrogram is composited on afterwards
        )

# Draw WITHOUT the plain column dendrogram and measure the column band + barplot-top
# so composite_dendrogram_panelC.py can align the labeled dendrogram onto it.
svg(file.path(save_dir, "/panel_C_consensusMin_Jorstad_noDend.svg"), width = 15, height = 7)
draw(p, padding = unit(c(6, 6, 6, 24), "mm"), heatmap_legend_side = "bottom")
decorate_heatmap_body("Eigenmodule", {
  tl <- deviceLoc(unit(0, "npc"), unit(1, "npc")); tr <- deviceLoc(unit(1, "npc"), unit(1, "npc"))
  x0 <<- convertX(tl$x, "inches", valueOnly = TRUE)
  x1 <<- convertX(tr$x, "inches", valueOnly = TRUE)
  ybody <<- convertY(tl$y, "inches", valueOnly = TRUE)   # heatmap-body top (crop line for panel B)
})
decorate_annotation("% of all mods", {
  tl <- deviceLoc(unit(0, "npc"), unit(1, "npc")); ybar <<- convertY(tl$y, "inches", valueOnly = TRUE)
})
dev.off()
writeLines(c(paste0("W_px=", 15 * 72), paste0("H_px=", 7 * 72),
             paste0("x0=", x0 * 72), paste0("x1=", x1 * 72),
             paste0("ybar_top=", (7 - ybar) * 72),
             paste0("ybody_top=", (7 - ybody) * 72)),
           file.path(save_dir, "jorstad_panel_coords.txt"))

# Save the heatmap colour legend on its own (matches the bottom legend above)
heatmap_legend <- Legend(col_fun = p@matrix_color_mapping@col_fun,
                         title = "Eigenmodules\n(PC1 of REIs\nfor merged\nmodules)",
                         title_gp = gpar(fontface = "plain"),
                         title_position = "topcenter",
                         direction = "vertical")

svg(file.path(save_dir, "/panel_C_consensusMin_Jorstad_legend.svg"), width = 1.5, height = 4)
grid.newpage()
draw(heatmap_legend)
dev.off()

pdf(file.path(save_dir, "/panel_C_consensusMin_Jorstad_legend.pdf"), width = 1.5, height = 4)
grid.newpage()
draw(heatmap_legend)
dev.off()

### Plot largest network (consensus min, Gabitto)
mod_eig <- fread(data.table = F,file = file.path(gab_net, "Module_eigengenes.csv")) |>
  column_to_rownames("Sample") |>
  select(jor_order)
rownames(mod_eig) <- class_info[,4]
colnames(mod_eig) <- 1:ncol(mod_eig)

# metacluster_min <- hclust(as.dist(1 - cor(mod_eig)), method="complete") |>
#   as.dendrogram()

kme <- fread(data.table = F, file = file.path(gab_net, "kME_table_.csv"))
mod_fdr <- tapply(kme[,1], kme[,6], list) |>
  lapply(\(x) as.numeric(gsub("Mod", "", x)))
mod_fdr_gene <- lapply(mod_fdr, \(x) mod_seed[x] |> unlist() |> unique())

# Count mods per meta-module
modcountdf <- data.frame("mod" = names(mod_fdr),
                    "FDR" = lapply(mod_fdr,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "mod_count") |>
  mutate(mod_per = mod_count/length(these_mods) * 100) |>
  arrange(mod_per)
genecountdf <- data.frame("mod" = names(mod_fdr_gene),
                          "FDR" = lapply(mod_fdr_gene,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "gene_count") |>
  mutate(gene_per = gene_count/18913 * 100) |>
  arrange(gene_per) 
countdf <- full_join(modcountdf[,-c(2:3)], genecountdf[,-c(2:3)], by = join_by(mod)) |>
  as.data.frame() 
countdf <- countdf[match(jor_order, countdf$mod), ]
countdf$mod <- 1:nrow(countdf)

p <- Heatmap(as.matrix(mod_eig),
        name = "Eigenmodule",
        cluster_rows = clustermin,
        cluster_columns = metacluster_min, # use the clustering from jorstad
        top_annotation = HeatmapAnnotation("% of all mods" = anno_barplot(countdf[,2],
                                                                          axis_param = list(at = c(0, 2))), 
                                           "% of all genes" = anno_barplot(countdf[,3],
                                                                           axis_param = list(at = c(0, 1)))),
        heatmap_legend_param = list(title_gp = gpar(fontface = "plain"),
                                    title_position = "topcenter"),
        column_title_side = "bottom",
        column_title = "Meta-modules",
        row_title = "1 - cor",
        row_names_side = "left",
        column_dend_height = unit(2, "cm"),
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 12),
        show_heatmap_legend = FALSE#,
       # show_column_dend = FALSE
        )

svg(file.path(save_dir, "/panel_C_consensusMin_Gabitto.svg"), width = 15, height = 7)
draw(p, padding = unit(c(6, 6, 6, 24), "mm"), heatmap_legend_side = "bottom")
decorate_column_dend("Eigenmodule", {
   grid.yaxis(gp = gpar(fontsize = 8)) 
})
dev.off()

pdf(file.path(save_dir, "/panel_C_consensusMin_Gabitto.pdf"), width = 15, height = 7)
draw(p, padding = unit(c(6, 6, 6, 24), "mm"), heatmap_legend_side = "bottom")
decorate_column_dend("Eigenmodule", {
   grid.yaxis(gp = gpar(fontsize = 8))
})
dev.off()


# ########
# Labeling branchpoints of eigenmodule dendrogram (v5 pipeline)
# ########
# v5 port of the v3 labeled-dendrogram pipeline (run_dendrogram.sh ->
# make_labeled_dendrogram_only_v1.2.py), keyed to the v5 Jorstad consensus-min
# network. Produces the two files the renderer needs, then renders:
#   mod_eig.csv                             - meta-module COLOUR-named eigengenes,
#                                             rows/cols ordered by the panel-C dendrograms
#   branchpoint_table_modeig_with_genes.csv - branchpoints with % mods and % genes
# (Formerly generated in v2 from a hand-loaded branchpoint_table_modeig.csv; here the
# branchpoint table is regenerated from the v5 dendrogram via build_branchpoint_table.)

# --- mod_eig.csv: colour-named eigengenes, ordered by the panel-C dendrograms ---
mod_eig_lab <- fread(data.table = F, file = file.path(jor_net, "Module_eigengenes.csv")) |>
  column_to_rownames("Sample")
rownames(mod_eig_lab) <- class_info[, 1]
mod_eig_lab <- mod_eig_lab[order.dendrogram(clustermin), order.dendrogram(metacluster_min)]
fwrite(mod_eig_lab, file = file.path(save_dir, "mod_eig.csv"), row.names = TRUE)

# --- branchpoint_table_modeig.csv: branchpoints of the meta-module dendrogram ---
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s8/build_branchpoint_table.R"))
bp_table <- build_branchpoint_table(
  mod_eig_file = file.path(save_dir, "mod_eig.csv"),
  cut_height   = 0.3,
  output_file  = file.path(save_dir, "branchpoint_table_modeig.csv")
)

# --- Add % genes and % mods per branchpoint (topmodposFDR meta-module membership) ---
bp_table <- bp_table |>
  mutate(Elements = map(Elements, ~ str_split(.x, ",\\s*")[[1]]))
kme <- fread(data.table = F, file = file.path(jor_net, "kME_table_.csv"))
mod_fdr <- tapply(kme[, 1], kme[, 6], list) |>
  lapply(\(x) as.numeric(gsub("Mod", "", x)))
mod_fdr_gene <- lapply(mod_fdr, \(x) mod_seed[x] |> unlist() |> unique())
numfdrgene <- length(unique(unlist(mod_seed)))
nummodgene <- length(mod_seed)
gene_per <- sapply(bp_table$Elements, \(x)
  length(unlist(mod_fdr_gene[names(mod_fdr_gene) %in% x])) / numfdrgene)
bp_table$Pct_of_Total_genes <- signif(gene_per * 100, 3)
mod_per <- sapply(bp_table$Elements, \(x)
  length(unlist(mod_fdr[names(mod_fdr) %in% x])) / nummodgene)
bp_table$Pct_of_Total <- signif(mod_per * 100, 3)
bp_table <- bp_table[, c(1:4, 6, 5)]
fwrite(bp_table, file = file.path(save_dir, "branchpoint_table_modeig_with_genes.csv"))

# --- Render the labeled dendrogram (coords mode) and composite it onto the Jorstad
#     panel_C, replacing that panel's top dendrogram. (matplotlib renderer -> svgutils
#     composite -> rsvg-convert for the PDF.) Only the Jorstad panel is modified. ---
render_py  <- Sys.getenv("PYTHON_BIN", "/home/gugene/miniconda3/bin/python")
render_scr <- file.path(save_dir, "make_labeled_dendrogram_only_v5.py")
composite  <- file.path(save_dir, "composite_dendrogram_panelC.py")
system2(render_py, c(shQuote(render_scr),
  shQuote(file.path(save_dir, "mod_eig.csv")),
  shQuote(file.path(save_dir, "branchpoint_table_modeig_with_genes.csv")),
  "-o", shQuote(file.path(save_dir, "jorstad_dendro_labeled.svg")),
  "--cut-height", "0.3", "--font-size", "14", "--fig-width", "36", "--fig-height", "8",
  "--coords-out", shQuote(file.path(save_dir, "jorstad_dendro_coords.txt"))))
system2(render_py, c(shQuote(composite), shQuote(save_dir)))
system2("rsvg-convert", c("-f", "pdf",
  "-o", shQuote(file.path(save_dir, "panel_C_consensusMin_Jorstad.pdf")),
  shQuote(file.path(save_dir, "panel_C_consensusMin_Jorstad.svg"))))


############
# Figure-assembly panels (see assemble_figure_v5.py)
#   panel_A.svg  - REI schematic (provided)
#   panel_B.svg  - labeled dendrogram + barplots (no heatmap)  [from composite above]
#   panel_C.svg  - Jorstad heatmap (row dendrogram, no top dendrogram/barplots/legend)
#   panel_D.svg  - Gabitto heatmap (same layout)
# C and D are drawn at width 15 (like the noDend panel) so their column band matches
# panel_B's -> the three stay column-aligned when placed at the same width in the pptx.
############
# xaxis: TRUE keeps the "Meta-modules" title + meta-module tick labels (panel D, the
# bottom panel); FALSE drops both (panel C, sandwiched between B and D). Also writes the
# heatmap column band (x0,x1 in px) so the assembly can align each panel's columns
# regardless of differing left/right margins (e.g. the barplot annotation text in B).
panel_heatmap <- function(me, out, coords_out, xaxis = TRUE){
  hm <- Heatmap(as.matrix(me),
    name = "Eigenmodule",
    cluster_rows = clustermin,
    cluster_columns = metacluster_min,
    show_column_dend = FALSE,
    column_title_side = "bottom",
    column_title = if (xaxis) "Meta-modules" else NULL,
    show_column_names = xaxis,
    row_title = "1 - cor",
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 12),
    show_heatmap_legend = FALSE)
  svg(out, width = 15, height = 5)
  draw(hm, padding = unit(c(6, 6, 1, 24), "mm"))       # tight top padding
  decorate_heatmap_body("Eigenmodule", {
    tl <- deviceLoc(unit(0, "npc"), unit(1, "npc")); tr <- deviceLoc(unit(1, "npc"), unit(1, "npc"))
    writeLines(c(paste0("W_px=", 15 * 72),
                 paste0("x0=", convertX(tl$x, "inches", valueOnly = TRUE) * 72),
                 paste0("x1=", convertX(tr$x, "inches", valueOnly = TRUE) * 72)),
               coords_out)
  })
  dev.off()
}

# Jorstad heatmap (reload: mod_eig was reassigned to Gabitto above); no x-axis (panel C)
jor_me <- fread(data.table = F, file = file.path(jor_net, "Module_eigengenes.csv")) |>
  column_to_rownames("Sample")
rownames(jor_me) <- class_info[, 4]; colnames(jor_me) <- 1:ncol(jor_me)
panel_heatmap(jor_me, file.path(save_dir, "panel_C.svg"),
              file.path(save_dir, "panel_C_coords.txt"), xaxis = FALSE)

# Gabitto heatmap (columns ordered to match Jorstad); keep x-axis (panel D, bottom)
gab_me <- fread(data.table = F, file = file.path(gab_net, "Module_eigengenes.csv")) |>
  column_to_rownames("Sample") |>
  select(all_of(jor_order))
rownames(gab_me) <- class_info[, 4]; colnames(gab_me) <- 1:ncol(gab_me)
panel_heatmap(gab_me, file.path(save_dir, "panel_D.svg"),
              file.path(save_dir, "panel_D_coords.txt"), xaxis = TRUE)

message("Figure panels written: panel_B.svg, panel_C.svg, panel_D.svg (panel_A.svg is provided)")



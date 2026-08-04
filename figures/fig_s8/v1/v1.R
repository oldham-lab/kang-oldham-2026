# Same as fig 5 but for MTG

library(ComplexHeatmap)
library(tidyverse)
library(data.table)
library(dendextend)
library(qs)
library(showtext)
showtext_auto()

version <- "v1"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/"), version)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

###############
# Panel A: Dendrogram - clustering subclasses over all projections
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

# Load indices (MTG)
rei1 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")) |>
  select(class_info[,1])

rei2 <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_con.csv")) |>
  select(class_info[,3])


# Filter modules
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
mod_seed <- mod_seed[these_mods]
mod_bc <- mod_bc[these_mods]
rei1 <- rei1[these_mods, ]
rei2 <- rei2[these_mods, ]
colnames(rei2) <- colnames(rei1)
rownames(rei1) <- 1:nrow(rei1)
rownames(rei2) <- 1:nrow(rei2)

# Remove zero variance projections
var1 <- apply(rei1, 1, var)
# which(var1==0)
# 972
var2 <- apply(rei2, 1, var)
# which(var2==0)
# 753

zerovar_vec <- c(753, 972)
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
clustermin <- hclust(as.dist(1 - ctm1), method="complete") |>
  as.dendrogram() 
svg(file.path(save_dir, "/panel_A_min.svg"), width = 6, height = 4)
par(mar = c(8.1, 4.1, 4.1, 2.1))
plot(clustermin, ylab = "1 - cor", 
     #main = "Clustering all subclasses over all modules\nusing min consensus matrix", 
     #cex.main = 1, 
     #font.main = 1
     )
dev.off()

# ##########
# # Panel B: FM on mod-mod consensus matrix
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

# No merge

# Gabitto MTG
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")) |>
  select(class_info[,1])
expr <- data.frame("Gene" = paste0("Mod", 1:length(these_mods)), proj[these_mods, ]) 
expr <- expr[-zerovar_vec, ] # remove zero-variance REI indices

# Gabitto DFC min
FindModules(
  projectname="Gabitto_MTG_consensusMin_noMerge",
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

# Jorstad DFC mean
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_con.csv")) |>
  select(class_info[,3])
expr <- data.frame("Gene" = paste0("Mod", 1:length(these_mods)), proj[these_mods, ]) 
expr <- expr[-zerovar_vec, ] # remove zero-variance REI indices

FindModules(
  projectname="Jorstad_MTG_consensusMin_noMerge",
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

load_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/v1"))

### Plot largest network (consensus min, Jorstad)
mod_eig <- fread(data.table = F,file = file.path(load_dir, "Jorstad_MTG_consensusMin_noMerge_Modules/Pearson-no_TO_signum0.702_minSize3_merge_ME_1_1014/Module_eigengenes.csv")) |>
  column_to_rownames("Sample")
rownames(mod_eig) <- class_info[,1]
jor_order <- colnames(mod_eig)

metacluster_min <- hclust(as.dist(1 - cor(mod_eig)), method="complete") |>
  as.dendrogram() 
qsave(metacluster_min, file = file.path(save_dir, "jorstad_consensusMin_dendro.qs"))

kme <- fread(data.table = F, file = file.path(load_dir, "Jorstad_MTG_consensusMin_noMerge_Modules/Pearson-no_TO_signum0.702_minSize3_merge_ME_1_1014/kME_table_.csv"))
mod_fdr <- tapply(kme[,1], kme[,6], list) |>
  lapply(\(x) as.numeric(gsub("Mod", "", x)))
mod_fdr_gene <- lapply(mod_fdr, \(x) mod_seed[x] |> unlist() |> unique()) 

# Count mods per meta-module
modcountdf <- data.frame("mod" = names(mod_fdr),
                    "FDR" = lapply(mod_fdr,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "mod_count") |>
  mutate(mod_per = mod_count/1023 * 100) |>
  arrange(mod_per) 
genecountdf <- data.frame("mod" = names(mod_fdr_gene),
                          "FDR" = lapply(mod_fdr_gene,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "gene_count") |>
  mutate(gene_per = gene_count/18913 * 100) |>
  arrange(gene_per) 
countdf <- full_join(modcountdf[,-c(2:3)], genecountdf[,-c(2:3)], by = join_by(mod)) |>
  as.data.frame() 
countdf <- countdf[match(colnames(mod_eig), countdf$mod), ]

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
        column_dend_height = unit(2, "cm"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 12),
        show_heatmap_legend = FALSE#,
        #show_column_dend = FALSE
        )

svg(file.path(save_dir, "/panel_C_consensusMin_Jorstad.svg"), width = 15, height = 7)
draw(p, padding = unit(c(6, 6, 6, 24), "mm"), heatmap_legend_side = "bottom")
decorate_column_dend("Eigenmodule", {
   grid.yaxis(gp = gpar(fontsize = 8)) 
})
dev.off()

pdf(file.path(save_dir, "/panel_C_consensusMin_Jorstad.pdf"), width = 15, height = 7)
draw(p, padding = unit(c(6, 6, 6, 24), "mm"), heatmap_legend_side = "bottom")
decorate_column_dend("Eigenmodule", {
   grid.yaxis(gp = gpar(fontsize = 8)) 
})
dev.off()

### Plot largest network (consensus min, Gabitto)
mod_eig <- fread(data.table = F,file = file.path(load_dir, "Gabitto_MTG_consensusMin_noMerge_Modules/Pearson-no_TO_signum0.702_minSize3_merge_ME_1_1014/Module_eigengenes.csv")) |>
  column_to_rownames("Sample") |>
  select(jor_order)
rownames(mod_eig) <- class_info[,1]

# metacluster_min <- hclust(as.dist(1 - cor(mod_eig)), method="complete") |>
#   as.dendrogram() 

kme <- fread(data.table = F, file = file.path(load_dir, "Gabitto_MTG_consensusMin_noMerge_Modules/Pearson-no_TO_signum0.702_minSize3_merge_ME_1_1014/kME_table_.csv"))
mod_fdr <- tapply(kme[,1], kme[,6], list) |>
  lapply(\(x) as.numeric(gsub("Mod", "", x)))
mod_fdr_gene <- lapply(mod_fdr, \(x) mod_seed[x] |> unlist() |> unique()) 

# Count mods per meta-module
modcountdf <- data.frame("mod" = names(mod_fdr),
                    "FDR" = lapply(mod_fdr,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "mod_count") |>
  mutate(mod_per = mod_count/1023 * 100) |>
  arrange(mod_per) 
genecountdf <- data.frame("mod" = names(mod_fdr_gene),
                          "FDR" = lapply(mod_fdr_gene,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "gene_count") |>
  mutate(gene_per = gene_count/18913 * 100) |>
  arrange(gene_per) 
countdf <- full_join(modcountdf[,-c(2:3)], genecountdf[,-c(2:3)], by = join_by(mod)) |>
  as.data.frame() 
countdf <- countdf[match(colnames(mod_eig), countdf$mod), ]

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
        column_names_gp = gpar(fontsize = 10),
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
# # Labeling branchpoints of eigenmodule dendrogram using Claude
# ########

# Write mod_eig.csv to disk
mod_eig <- fread(data.table = F,file = file.path(load_dir, "Jorstad_MTG_consensusMin_noMerge_Modules/Pearson-no_TO_signum0.702_minSize3_merge_ME_1_1014/Module_eigengenes.csv")) |>
  column_to_rownames("Sample")
rownames(mod_eig) <- class_info[,1]
jor_order <- colnames(mod_eig)

metacluster_min <- hclust(as.dist(1 - cor(mod_eig)), method="complete") |>
  as.dendrogram() 

mod_eig <- mod_eig[order.dendrogram(clustermin), order.dendrogram(metacluster_min)]
fwrite(mod_eig, file = file.path(save_dir, "mod_eig.csv"), row.names = T)

# Build branchpoint table
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s7/build_branchpoint_table.R"))
bp_table <- build_branchpoint_table(
  mod_eig_file = file.path(save_dir, "mod_eig.csv"),
  cut_height   = 0.3,
  output_file  = file.path(save_dir, "branchpoint_table_modeig.csv")
)

# Convert the Elements column from a comma-separated string to a list column
bp_table <- bp_table %>%
  mutate(Elements = map(Elements, ~ str_split(.x, ",\\s*")[[1]]))

# Load meta-modules (topmodposfdr) (Jorstad min)
kme <- fread(data.table = F, file = file.path(save_dir, "Jorstad_MTG_consensusMin_noMerge_Modules/Pearson-no_TO_signum0.702_minSize3_merge_ME_1_1014/kME_table_.csv"))
mod_fdr <- tapply(kme[,1], kme[,6], list) |>
  lapply(\(x) as.numeric(gsub("Mod", "", x)))
mod_fdr_gene <- lapply(mod_fdr, \(x) mod_seed[x] |> unlist() |> unique()) 

# Calculate total # of genes and modules represented by modules used in metamodule network
numfdrgene <- length(unique(unlist(mod_seed)))
# 12937
nummodgene <- length(mod_seed)

# Calculate %gene per branchpoint
gene_per <- lapply(bp_table$Elements, \(x){
  gene_temp <- mod_fdr_gene[names(mod_fdr_gene) %in% x] |> unlist() |> length()
  return(gene_temp / numfdrgene)
}) |> unlist()
bp_table$Pct_of_Total_genes <- signif(gene_per * 100, 3)

# Recalculate %mod per branchpoint
mod_per <- lapply(bp_table$Elements, \(x){
  mod_temp <- mod_fdr[names(mod_fdr) %in% x] |> unlist() |> length()
  return(mod_temp / nummodgene)
}) |> unlist()
bp_table$Pct_of_Total <- signif(mod_per * 100, 3)

bp_table <- bp_table[, c(1:4, 6, 5)]
fwrite(bp_table, file = file.path(save_dir, "branchpoint_table_modeig_with_genes.csv"))
# Use this modified bp_table to add %mod and %gene information to the dendrogram


# Run python script
# python make_labeled_dendrogram_only_v1.2.py /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_s7/v1/mod_eig.csv /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_s7//v1/branchpoint_table_modeig_with_genes.csv



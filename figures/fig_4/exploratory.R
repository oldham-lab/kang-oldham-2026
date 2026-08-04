library(tidyverse)
library(data.table)
library(ggrepel)
library(grid)

#######################################
# Load required data
#######################################
# Subclass/Class info
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

# Module data
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])


#######################################
# Run FindModules on module projections 
#######################################
# - Gabitto DFC AllADVsCon, REI, Control
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")) |>
  select(class_info[,1])
expr <- data.frame("Gene" = 1:length(these_mods), proj[these_mods, ]) 

# Celltypes are "samples", modules are "genes"
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

setwd(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/"))

# Bicor
FindModules(
  projectname="modIndices_exploratory_minMEcor0.95_Gabitto_DFC",
  data_cols=expr,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(expr),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(5, 8, 10, 15, 20),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999,.999, .99, .98,.97 ,.96, .95),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.95,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

# Pearson
FindModules(
  projectname="modIndices_exploratory_minMEcor0.95_pearson_Gabitto_DFC",
  data_cols=expr,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(expr),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Pearson"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3, 5, 8, 10, 15, 20),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999,.999, .99, .98,.97 ,.96, .95),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.95,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)
# Pearson modules look much better.

## Summary plots (bicor)
# Plot # of modules associated with each pattern
# Largest network: 0.685_5
kme <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/kME_table_.csv"))
mod_bc <- tapply(kme[,1], kme[,5], list)
mod_fdr <- tapply(kme[,1], kme[,6], list)

outdf <- data.frame("mod" = names(mod_bc),
                    "BC" = lapply(mod_bc, length) |> unlist(),
                    "FDR" = lapply(mod_fdr,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "count") |>
  arrange(count) 
ord <- outdf |> group_by(mod) |> summarise(sum_count = sum(count)) |> arrange(sum_count)
p <- outdf |>
  mutate(mod = factor(mod, levels = ord$mod)) |>
  ggplot(aes(x = mod, y = count, fill = sig_cut)) +
    theme_classic() +
    geom_bar(stat = "identity") +
    labs(x = "Pattern", y = "# of associated mods", title = "# of mods per meta-module pattern") +
    theme(legend.title = element_blank()) +
    coord_flip()
ggsave(p,file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/mod_count_per_pattern.pdf"))

# Plot # of genes per meta-module
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(data.table=F, file = file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)

mod_bc_gene <- lapply(mod_bc, \(x) mods[x] |> unlist() |> unique()) 
mod_fdr_gene <- lapply(mod_fdr, \(x) mods[x] |> unlist() |> unique()) 

outdf <- data.frame("mod" = names(mod_bc_gene),
                    "BC" = lapply(mod_bc_gene, length) |> unlist(),
                    "FDR" = lapply(mod_fdr_gene,length) |> unlist()) |>
  pivot_longer(!mod, names_to = "sig_cut", values_to = "count") |>
  arrange(count) 
ord <- outdf |> group_by(mod) |> summarise(sum_count = sum(count)) |> arrange(sum_count)
p <- outdf |>
  mutate(mod = factor(mod, levels = ord$mod)) |>
  ggplot(aes(x = mod, y = count, fill = sig_cut)) +
    theme_classic() +
    geom_bar(stat = "identity") +
    labs(x = "Pattern", y = "# of associated genes (topmodposbc)", title = "# of genes (topmodposbc) per meta-module pattern") +
    theme(legend.title = element_blank()) +
    coord_flip()
ggsave(p,file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/gene_count_per_pattern.pdf"))

# Pie chart of genes per meta-module
pie_bc <- outdf |> filter(sig_cut == "BC")
pie_bc[36, 1] <- c("Unassigned")
pie_bc[36, 2] <- "BC"
pie_bc[36, 3] <- 18913 - sum(pie_bc[-36, 3])
pie_bc <- pie_bc  |>
  mutate(prop = count/18913,
         lab = paste0(mod," (", sprintf("%.2f", prop * 100), "%)" ))
col_bc <- pie_bc$mod
col_bc[36] <- "white"

pie_fdr <- outdf |> filter(sig_cut == "FDR")
pie_fdr[36, 1] <- c("Unassigned")
pie_fdr[36, 2] <- "FDR"
pie_fdr[36, 3] <- 18913 - sum(pie_fdr[-36, 3])
pie_fdr <- pie_fdr  |>
  mutate(prop = count/18913,
         lab = paste0(mod," (", sprintf("%.2f", prop * 100), "%)" ))
col_fdr <- pie_fdr$mod
col_fdr[36] <- "white"

pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/gene_pie_all.pdf"))
pie(pie_bc$prop , labels = pie_bc$lab, col = col_bc, main = "% of genes per meta-module (BC)", cex = 0.6)
pie(pie_fdr$prop , labels = pie_fdr$lab, col = col_fdr, main = "% of genes per meta-module (FDR)", cex = 0.6)
dev.off()

pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/gene_pie_bc.pdf"))
pie(pie_bc$prop , labels = pie_bc$lab, col = col_bc, main = "% of genes per meta-module (BC)", cex = 0.6)
dev.off()

pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/gene_pie_fdr.pdf"))
pie(pie_fdr$prop , labels = pie_fdr$lab, col = col_fdr, main = "% of genes per meta-module (FDR)", cex = 0.6)
dev.off()

# %VE by each meta-module
eig <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/Module_eigengenes.csv"))
exprt <- t(expr[,-1])
vedf <- lapply(1:length(mod_bc), \(i){
  vevec <- lapply(mod_bc[[i]], \(j){
    summary(lm(exprt[,j] ~ eig[ ,i+1]))$r.squared
  })
  data.frame("metamod" = names(eig)[i + 1], 
             "pcnt_ve" = unlist(vevec))
}) |> do.call(what = "rbind") 
ordve <- vedf |> group_by(metamod) |> summarise(mean_ve = mean(pcnt_ve)) |> arrange(mean_ve) 
vedf$metamod <- factor(vedf$metamod, levels = unique(ordve$metamod))
p <- ggplot(vedf, aes(x = metamod, y = pcnt_ve, fill = metamod)) + 
  theme_classic() +
  #geom_violin() +
  geom_boxplot(width = 0.5, alpha = 0.8, notch = T, outlier.shape = NA) + 
  geom_jitter(alpha = 0.2) + 
  labs(x = "Meta-module", y = "% variance explained", title = "% variance explained by meta-module eigengene") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = unique(ordve$metamod)) +
  coord_flip()
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/pcnt_VE.pdf"))

## Summary plots (Pearson cor mods)
# Cluster subclasses based on metamodules
kme <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Modules/Pearson-no_TO_signum0.85_minSize3_merge_ME_0.95_1157/kME_table_.csv")) 
mod <- tapply(kme[,1], kme[,2], list)
mod_eig <- fread(data.table = F,file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Modules/Pearson-no_TO_signum0.85_minSize3_merge_ME_0.95_1157/Module_eigengenes.csv"))

cormat <- cor(t(mod_eig[, -1]))
rownames(cormat) <- class_info[,1]
colnames(cormat) <- class_info[,1]
library(ComplexHeatmap)
p <- Heatmap(cormat,
             name = "Cor",
             column_title = "Distance between subclasses based on meta-module eigengenes")
pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_Modules/Bicor-no_TO_signum0.685_minSize5_merge_ME_0.95_1157/subclass_heatmap.pdf"))
draw(p)
dev.off()

###################################
# FM networks for different regions
###################################

# Gabitto MTG
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")) |>
  select(class_info[,1])
expr <- data.frame("Gene" = 1:length(these_mods), proj[these_mods, ]) 

FindModules(
  projectname="modIndices_exploratory_minMEcor0.95_pearson_Gabitto_MTG",
  data_cols=expr,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(expr),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Pearson"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3, 5, 8, 10, 15, 20),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999,.999, .99, .98,.97 ,.96, .95),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.95,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

# Jorstad MTG
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_con.csv")) |>
  select(class_info[,3])
expr <- data.frame("Gene" = 1:length(these_mods), proj[these_mods, ]) 

FindModules(
  projectname="modIndices_exploratory_minMEcor0.95_pearson_Jorstad_MTG",
  data_cols=expr,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(expr),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Pearson"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3, 5, 8, 10, 15, 20),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999,.999, .99, .98,.97 ,.96, .95),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.95,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

# Jorstad DFC
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_REI/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |>
  select(class_info[,3])
expr <- data.frame("Gene" = 1:length(these_mods), proj[these_mods, ]) 

FindModules(
  projectname="modIndices_exploratory_minMEcor0.95_pearson_Jorstad_DFC",
  data_cols=expr,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(expr),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Pearson"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3, 5, 8, 10, 15, 20),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999,.999, .99, .98,.97 ,.96, .95),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.95,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

###################################################
# Categorize mods based on coefficient of variation
###################################################

# Possible ways to categorize mods:
# - A spectrum ranging from completely uniform spread to expressed in only a single celltype
# - Distinguishing class-specific modules (excitatory vs inhibitory; neuronal vs non-neuronal)

proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")) |>
  select(!module)

cv_func <- function(x){
  sd(x)/mean(x)
}

out_df <- data.frame("Module" = 1:nrow(proj),
                     "CV" = apply(proj, 1, cv_func)) 
# Higher CV = Expressed in single celltype

# p <- ggplot(out_df, aes(x = CV)) + 
#   geom_density()
# ggsave(p, file = "~/test/test.pdf")
# # Long right tail

out_df |> 
  arrange(desc(CV)) |>
  head()

# Model each module as a function of cell class.
# Load class info:
class_info <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |>
  arrange(Subclass)
class_info <- class_info[match(colnames(proj)[-25], class_info[,1]), ]

# Model each module by class:
class_lm <- c()
for(x in 1:nrow(proj)){
  class_lm[x] <- summary(lm(t(proj[x,]) ~ class_info$Class))$r.squared
}
which(is.na(class_lm)) # 1096: histone genes, all zeroes in Gabitto

# Set colors for each quadrant:
pal3 <- RColorBrewer::brewer.pal(4, "Set1")

out_df2 <- out_df |>
  mutate(class_r2 = class_lm,
         VE = qs::qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/module_VE.qs"))) # % variance explained by PC1 per module
out_df2$quad <- cut(out_df2$CV, breaks = c(0, 1, Inf), labels = c("left", "right"))
out_df2$quad <- interaction(out_df2$quad, cut(out_df2$class_r2, breaks = c(0, .40, 1.00), labels = c("bottom", "top")))

# p3 <- ggplot(out_df2, aes(x = CV, y = VE)) + 
#   geom_point()

# ggsave(p3, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/cv_vs_pcntVE.pdf"))

a <- 0.2
p2 <- ggplot(out_df2, aes(x = CV, y = class_r2 * 100, color = quad)) +
  theme_classic() +
  geom_point(alpha = a, shape = 16) +
  geom_point(data = out_df2[out_df2$Module %in% c(347, 635, 822, 253), ], alpha = 1, color = "black", shape = 16) +
  annotate("label", x = 1, y = 20, label = "Non-specific", color = pal3[1], alpha = a) +
  annotate("label", x = 3.2, y = 20, label = "Celltype-specific", color = pal3[2], alpha = a) +
  annotate("label", x = 1, y = 70, label = "Class-specific", color = pal3[3], alpha = a) +
  geom_text_repel(data = out_df2[out_df2$Module %in% c(347, 635, 822, 253), ], aes(label = Module),
                  box.padding = 0.5, max.overlaps = Inf, color = "Black", force_pull = 0.05) +
  labs(x = "Coefficient of variation (REI)",
       y = "% variance explained by cell class") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("left.bottom" = pal3[1], "right.bottom" = pal3[2], "left.top" = pal3[3], "right.top" = pal3[4]))

ggsave(p2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/cv_vs_ClassR2.pdf"), width = 4, height = 3)

# Ported CoV analysis to figure 4.

#############
# Assess meta-module conservation
############

# Load list of input projection matrices
expr_list <- list(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv"), # Gabitto DFC
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv"), # Gabitto MTG
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_REI/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv"), # Jorstad DFC
                  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_con.csv") # Jorstad MTG
                  ) 
expr_list[1:2] <- lapply(expr_list[1:2], \(x){
    proj <- fread(data.table = F, file = x) |>
      select(class_info[,1])
    expr <- data.frame("Gene" = 1:length(these_mods), proj[these_mods, ]) 
    return(expr)
  })
expr_list[3:4] <- lapply(expr_list[3:4], \(x){
    proj <- fread(data.table = F, file = x) |>
      select(class_info[,3])
    expr <- data.frame("Gene" = 1:length(these_mods), proj[these_mods, ]) 
    return(expr)
  })
dat_names <- c("Gabitto_DFC", "Gabitto_MTG", "Jorstad_DFC", "Jorstad_MTG")

# Paths Gabitto/Jorstad DFC/MTG networks
network_folders <- list.dirs(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory"), recursive = F)[-5] 

minSizevec = c(3, 5, 8, 10, 15, 20)
signumvecr = lapply(c(.9999,.999, .99, .98,.97 ,.96, .95), \(x) rep(x, 6)) |> unlist()

outdf <- lapply(seq_along(network_folders), \(i){ # for each network
  net_names <- list.dirs(network_folders[i], full.names = F, recursive = F)
  nets <- list.dirs(network_folders[i], recursive = F)
  params <- lapply(net_names, \(x){
    split <- strsplit(x, "_")
    out <- data.frame("signum" = gsub("signum", "", split[[1]][3]) |> as.numeric(),
                      "minsize" = gsub("minSize", "", split[[1]][4]) |> as.numeric())
    return(out)
  }) |> do.call(what = "rbind") |> mutate("index" = 1:n()) |> arrange(desc(signum)) |>
    mutate(relsignum = signumvecr) |>
    arrange(index)

  netdf <- lapply(seq_along(nets), \(j){ # for each parameter level
    ## Calculate % VE by meta-module eigengene
    # Load meta-modules (topmodposbc)
    if(file.exists(file.path(nets[j], "kME_table_.csv"))){ # if kme table exists
      kme <- fread(data.table = F, file = file.path(nets[j], "kME_table_.csv"))
      mod_bc <- tapply(kme[,1], kme[,5], list)
      if(length(mod_bc) > 0){ # if there are topmodposbc modules
        # Load meta-module eigengene
        eig <- fread(data.table = F, file = file.path(nets[j], "Module_eigengenes.csv")) |>
          select(!Sample) |>
          select(names(mod_bc))
        # Calculate %VE eexplained
        expr <- expr_list[[i]]
        exprt <- t(expr[,-1])
        vedf <- lapply(1:length(mod_bc), \(a){
          vevec <- lapply(mod_bc[[a]], \(b){
            summary(lm(exprt[,b] ~ eig[ ,a]))$r.squared
          }) |> unlist()
          data.frame("metamod" = names(eig)[a], 
                     "pcnt_ve" = mean(vevec),
                     "metamod_count" = length(mod_bc))
        }) |> do.call(what = "rbind") |>
          mutate("network" = net_names[j],
                 "dataset" = dat_names[i],
                 "relsignum" = params[j, 4],
                 "minsize" = params[j, 2])
        return(vedf) 
      } else {
        return(data.frame())
      }
    } else {
      return(data.frame())
    }
  }) |> do.call(what = "rbind")
  return(netdf)
}) |> do.call(what = "rbind")

p <- outdf |> 
 # filter(minsize==3, relsignum %in% c(0.95, 0.99)) |>
  mutate(relsignum = factor(relsignum, levels = unique(sort(relsignum))),
         dataset = factor(dataset, levels = unique(dataset)),
         minsize = factor(paste0("Minsize ", minsize), levels = paste0("Minsize ", unique(sort(minsize))))) |>
  ggplot(aes(x = dataset, y = pcnt_ve, fill = relsignum)) +
    theme_bw() +
    #geom_jitter() + 
    geom_boxplot(alpha = 0.8) + 
    facet_wrap(~minsize) +
    ylim(0,1) +
    labs(x = "", y = "Mean % variance explained\nby meta-module eigengene")
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/Exploratory_analysis_pcntVE_comparison_across_FM_networks.pdf"), width = 20)

# asdf <- exprt[,mod_bc[[a]]] 
# qwer <- WGCNA::moduleEigengenes(asdf, colors = rep("bla", ncol(asdf)))$eigengenes
#   pc1_bc <-  prcomp(asdf, scale=T)$x[,1]
# WGCNA::moduleEigengenes(asdf, colors = rep("bla", ncol(asdf)))$varExplained

# Use PC1VE from module statistics table instead
outdf <- lapply(seq_along(network_folders), \(i){ # for each network
  net_names <- list.dirs(network_folders[i], full.names = F, recursive = F)
  nets <- list.dirs(network_folders[i], recursive = F)
  params <- lapply(net_names, \(x){
    split <- strsplit(x, "_")
    out <- data.frame("signum" = gsub("signum", "", split[[1]][3]) |> as.numeric(),
                      "minsize" = gsub("minSize", "", split[[1]][4]) |> as.numeric())
    return(out)
  }) |> do.call(what = "rbind") |> mutate("index" = 1:n()) |> arrange(desc(signum)) |>
    mutate(relsignum = signumvecr) |>
    arrange(index)

  netdf <- lapply(seq_along(nets), \(j){ # for each parameter level
    ## Calculate % VE by meta-module eigengene
    # Load meta-modules (topmodposbc)
    if(file.exists(file.path(nets[j], "kME_table_.csv"))){ # if kme table exists
      kme <- fread(data.table = F, file = file.path(nets[j], "kME_table_.csv"))
      mod_bc <- tapply(kme[,1], kme[,5], list)
      if(length(mod_bc) > 0){ # if there are topmodposbc modules
        # Load meta-module eigengene
        modstat <- fread(data.table = F, file = file.path(nets[j], "Module_statistics.csv")) 
        # Calculate %VE eexplained
        vedf <- data.frame("metamod" = names(eig)[a], 
                           "pcnt_ve" = modstat$PC1VE,
                           "metamod_count" = nrow(modstat),
                           "network" = net_names[j],
                           "dataset" = dat_names[i],
                           "relsignum" = params[j, 4],
                           "minsize" = params[j, 2])
        return(vedf) 
      } else {
        return(data.frame())
      }
    } else {
      return(data.frame())
    }
  }) |> do.call(what = "rbind")
  return(netdf)
}) |> do.call(what = "rbind")

p <- outdf |> 
 # filter(minsize==3, relsignum %in% c(0.95, 0.99)) |>
  mutate(relsignum = factor(relsignum, levels = unique(sort(relsignum))),
         dataset = factor(dataset, levels = unique(dataset)),
         minsize = factor(paste0("Minsize ", minsize), levels = paste0("Minsize ", unique(sort(minsize))))) |>
  ggplot(aes(x = dataset, y = pcnt_ve, fill = relsignum)) +
    theme_bw() +
    #geom_jitter() + 
    geom_boxplot(alpha = 0.8) + 
    facet_wrap(~minsize) +
    ylim(0,1) +
    labs(x = "", y = "Mean % variance explained\nby meta-module eigengene")
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/Exploratory_analysis_pcntVE_comparison_across_FM_networks_modStat.pdf"), width = 20)



# # # test
# # source("~/root_dir/home/shared/code/FindModules/FindModules094.R")

# proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")) |>
#   select(class_info[,1])
# expr <- data.frame("Gene" = 1:length(these_mods), proj[these_mods, ]) 

# setwd("~/test/")
# FindModules(
#   projectname="testfm",
#   expr=expr,
#   geneinfo=c(1),
#   sampleindex=c(2:ncol(expr)),
#   samplegroups=NULL,
#   subset=NULL,
#   simMat=NULL,
#   saveSimMat=FALSE,
#   simType="Pearson",
#   beta=1,
#   overlapType="None",
#   TOtype="signed",
#   TOdenom="min",
#   MIestimator="mi.mm",
#   MIdisc="equalfreq",
#   signumType="rel",
#   iterate=TRUE,
#   signumvec=c(.96),
#   minsizevec=c(3),
#   signum=NULL,
#   minSize=NULL,
#   minMEcor=0.95,
#   ZNCcut=2,
#   calcSW=FALSE,
#   loadTree=FALSE,
#   writeKME=TRUE,
#   calcBigModStat=FALSE
# )


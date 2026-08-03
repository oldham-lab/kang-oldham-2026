library(ComplexHeatmap)
library(tidyverse)
library(data.table)
library(dendextend)
library(qs)

version_folder <- "v1.1"

###########################
# Panel A: REI matrix schematic
###########################

rei_proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_REI/indices_over_all_datasets_Cell_Type_1.csv")) |>
  select(!module)

rei_subset <- rei_proj[1:5, 1:5]
h_col <- circlize::colorRamp2(c(0, 1), c("white", "grey"))
p <- Heatmap(as.matrix(rei_subset),
             name = "REI",
             col = h_col,
             cluster_columns = F,
             cluster_rows = F,
             row_title = "Modules (REI)",
             column_title = "Celltypes",
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         title_gp = gpar(fontface = "plain", fontsize = 10)),
             show_heatmap_legend = F,
             width = ncol(rei_subset)*unit(5, "mm"), 
             height = nrow(rei_subset)*unit(5, "mm"))

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, "/panel_A.pdf"), width = 2.5, height = 2.5)
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


rei <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")) |>
  select(class_info[,1])

# Filter modules
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)
rei <- rei[these_mods, ]
#cols <- RColorBrewer::brewer.pal(3, "Set1")


# Cluster subclasses
cluster1 <- hclust(as.dist(1 - cor(rei)), method="complete") |>
  as.dendrogram() #|> 
  # branches_attr_by_labels(colnames(rei)[c(4:12)], cols[3]) |>
  # branches_attr_by_labels(colnames(rei)[c(2, 13:14, 18:23)], cols[2]) |>
  # branches_attr_by_labels(colnames(rei)[c(1, 3, 15:17, 24)], cols[1]) |>  
  # color_labels(labels = colnames(rei)[c(4:12)], col = cols[3]) |>
  # color_labels(labels = colnames(rei)[c(2, 13:14, 18:23)], col = cols[2]) |>
  # color_labels(labels = colnames(rei)[c(1, 3, 15:17, 24)], col = cols[1]) 

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, "/panel_B.pdf"), width = 6, height = 4)
par(mar = c(8.1, 4.1, 4.1, 2.1))
plot(cluster1, ylab = "1 - cor", main = "Clustering all subclasses over all modules", cex.main = 1, font.main = 1)
dev.off()

###############
# Panel C: Metamodule dendrogram + heatmap
###############

mod_eig <- fread(data.table = F,file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Modules/Pearson-no_TO_signum0.85_minSize3_merge_ME_0.95_1157/Module_eigengenes.csv")) |>
  column_to_rownames("Sample")
rownames(mod_eig) <- class_info[,1]

kme <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Modules/Pearson-no_TO_signum0.85_minSize3_merge_ME_0.95_1157/kME_table_.csv"))
mod_fdr <- tapply(kme[,1], kme[,6], list)
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
        cluster_rows = cluster1,
        top_annotation = HeatmapAnnotation("% of all mods" = anno_barplot(countdf[,2],
                                                                          axis_param = list(at = c(0, 5))), 
                                           "% of all genes" = anno_barplot(countdf[,3],
                                                                           axis_param = list(at = c(0, 3)))),
        heatmap_legend_param = list(title_gp = gpar(fontface = "plain"),
                                    title_position = "topcenter"),
        column_title_side = "bottom",
        column_title = "Meta-modules",
        row_title = "1 - cor",
        row_names_side = "left",
        column_dend_height = unit(2, "cm")
        )

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, "/panel_C.pdf"), width = 15, height = 8)
draw(p, padding = unit(c(6, 6, 6, 24), "mm"), heatmap_legend_side = "bottom")
decorate_column_dend("Eigenmodule", {
   grid.yaxis(gp = gpar(fontsize = 8)) 
})
# decorate_row_dend("Eigenmodule", {
#    grid.xaxis(gp = gpar(fontsize = 8)) 
# })
dev.off()

################
# panel D
# (old panel 4Q)
################
mod_id <- c(# 347, # oligo
            # 635, # VIP although VIP is not a top 10 seed gene
            # 822, # Endo
            681, # Microglia
            669, # left top, postsynaptic
            # 253, # mystery module
            1007, # jun/fos, top right
            173 # left bottom, electron transport chain
            )
# Fetch vector of mods filtered by CoPA
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

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
# which(is.na(class_lm)) # 1096: histone genes, all zeroes in Gabitto

# Set colors for each quadrant:
pal3 <- RColorBrewer::brewer.pal(4, "Set1")

out_df2 <- out_df |>
  mutate(class_r2 = class_lm,
         VE = qs::qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/module_VE.qs"))) # % variance explained by PC1 per module
out_df2$quad <- cut(out_df2$CV, breaks = c(0, 1, Inf), labels = c("left", "right"))
out_df2$quad <- interaction(out_df2$quad, cut(out_df2$class_r2, breaks = c(0, .40, 1.00), labels = c("bottom", "top")))
out_df2 <- out_df2[these_mods, ]
a <- 0.9
p2 <- ggplot(out_df2, aes(x = CV, y = class_r2 * 100, color = quad)) +
  theme_classic() +
  geom_point(alpha = 0.2, shape = 16) +
  geom_point(data = out_df2[out_df2$Module %in% mod_id, ], alpha = 1, color = "black", shape = 16) +
  annotate("text", x = 0.5, y = 20, label = "Non-specific", color = pal3[1], alpha = a) +
  annotate("text", x = 3.2, y = 20, label = "Celltype-specific", color = pal3[2], alpha = a) +
  annotate("text", x = 0.5, y = 80, label = "Class\nspecific", color = pal3[3], alpha = a) +
  annotate("text", x = 2, y = 70, label = "Non-neuronal", color = pal3[4], alpha = a) +
  ggrepel::geom_label_repel(data = out_df2[out_df2$Module %in% mod_id, ], aes(label = Module),
                  box.padding = 0.5, max.overlaps = Inf, color = "Black", force_pull = 0.05, alpha = 0.8) +
  labs(x = "Coefficient of variation (REI)",
       y = "% variance explained by cell class") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("left.bottom" = pal3[1], "right.bottom" = pal3[2], "left.top" = pal3[3], "right.top" = pal3[4]))

ggsave(p2, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, "/panel_D.pdf"), width = 4, height = 3)

###########
# Panel E
# - Proportion of mods/genes in each category from panel D
###########
mod_props <- c(table(out_df2$quad)) 
mod_bc_filt <- mod_bc[these_mods]
gene_list_4 <- tapply(out_df2$Module, out_df2$quad, list) |>
  lapply(\(x) unlist(mod_bc_filt[x]) |> unique())
gene_props <- lapply(gene_list_4, length) |> unlist()
gene_props[5] <- 18913 - sum(gene_props[1:4])

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, "/panel_E1.pdf"), width = 7, height = 4, bg = "transparent")
pie(mod_props, 
    labels = c("40%","8%", "50%","2%"),
    col = pal3,
    cex = 1.2)
dev.off()

pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, "/panel_E2.pdf"), width = 7, height = 4, bg = "transparent")
pie(gene_props, 
    labels = c("34%","5%","46%","2%", "Unassigned (13%)"),
    col = c(pal3, "grey"),
    cex = 1.2)
dev.off()

##########
# Panel F
# displaying violin plots of the distributions of REI correlations for each category shown in panel g.  
########## 
proj_filt <- proj[these_mods, ]
proj_cor <- cor(t(proj_filt))
mod_list_4 <- tapply(out_df2$Module, out_df2$quad, list)

corlist4 <- mapply(\(x, y){
  t <- colnames(proj_cor) %in% x
  temp <- proj_cor[t, t]
  outdf <- data.frame("type" = y,
                      "cor" = temp[upper.tri(temp)]) |>
    filter(!is.na(cor))
  return(outdf)
}, mod_list_4, names(mod_list_4), SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate("label" = case_match(type,
    "left.bottom" ~ "Non-specific",
    "right.bottom" ~ "Celltype-specific",
    "left.top" ~ "Class-specific",
    "right.top" ~ "Non-neuronal",
    .default = type))

p <- ggplot(corlist4, aes(x = label, y = cor, fill = type)) + 
  theme_classic() + 
  geom_violin() + 
  labs(x = "", y = "Pairwise correlations\n(REI projections)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("left.bottom" = pal3[1], "right.bottom" = pal3[2], "left.top" = pal3[3], "right.top" = pal3[4])) +
  coord_flip()
ggsave(p, file = file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, "/panel_F.pdf"), width = 4, height = 2) 



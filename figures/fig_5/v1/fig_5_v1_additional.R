library(ComplexHeatmap)
library(tidyverse)
library(data.table)
library(dendextend)
library(qs)
library(showtext)
showtext_auto()
version_folder <- "v2"

# Expanding on Fig. 5 code

######
# Run panel B analysis on 4 datasets (Jorstad or Gabitto, DFC or MTG)
#######
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


rei_list <- list("Gabitto_DFC" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")),
                 "Jorstad DFC" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_proj_indices/log_REI/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")),
                 "Gabitto_MTG" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_Control.csv")),
                 "Jorstad MTG" = fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinMTG/sn_proj_indices/log_REI/indices_over_all_datasets_Subclass_con.csv"))) |>
  lapply(\(x) x |> select(!module))
rei_list[c(2,4)] <- lapply(rei_list[c(2,4)], \(x){
   colnames(x) <- class_info$Subclass[match(colnames(x), class_info$Subclass_fixed)]
   return(x)
})
rei_list <- lapply(rei_list, \(x) x |> select(class_info[,1]))
rei_names <- c("Gabitto_DFC", "Jorstad_DFC", "Gabitto_MTG", "Jorstad_MTG")

cluster_list <- list()
for(i in seq_along(rei_list)){
    rei <- rei_list[[i]]
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
    as.dendrogram() 

    cluster_list[[i]] <- cluster1

    pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, paste0("/panel_B_",rei_names[i], ".pdf")), width = 6, height = 4)
    par(mar = c(8.1, 4.1, 4.1, 2.1))
    plot(cluster1, ylab = "1 - cor", main = "Clustering all subclasses over all modules", cex.main = 1, font.main = 1)
    dev.off()

    svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, paste0("/panel_B_",rei_names[i], ".svg")), width = 6, height = 4)
    par(mar = c(8.1, 4.1, 4.1, 2.1))
    plot(cluster1, ylab = "1 - cor", main = "Clustering all subclasses over all modules", cex.main = 1, font.main = 1)
    dev.off()
}

##########
# Run panel C on all datasets
##########

# Pick largest network from each dataset
mod_eig_list <- list(fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Gabitto_DFC_Modules/Pearson-no_TO_signum0.772_minSize3_merge_ME_0.95_1022/Module_eigengenes.csv")),
                     fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Jorstad_DFC_Modules/Pearson-no_TO_signum0.746_minSize3_merge_ME_0.95_1022/Module_eigengenes.csv")),
                     fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Gabitto_MTG_Modules/Pearson-no_TO_signum0.796_minSize3_merge_ME_0.95_1022/Module_eigengenes.csv")),
                     fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Jorstad_MTG_Modules/Pearson-no_TO_signum0.75_minSize3_merge_ME_0.95_1022/Module_eigengenes.csv"))
                     ) |>
  lapply(\(x){
    out <- x |> column_to_rownames("Sample")
    rownames(out) <- class_info[,1]
    return(out)
  })

kme_list <- list(fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Gabitto_DFC_Modules/Pearson-no_TO_signum0.772_minSize3_merge_ME_0.95_1022/kME_table_.csv")),
                 fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Jorstad_DFC_Modules/Pearson-no_TO_signum0.746_minSize3_merge_ME_0.95_1022/kME_table_.csv")),
                 fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Gabitto_MTG_Modules/Pearson-no_TO_signum0.796_minSize3_merge_ME_0.95_1022/kME_table_.csv")),
                 fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/exploratory/modIndices_exploratory_minMEcor0.95_pearson_Jorstad_MTG_Modules/Pearson-no_TO_signum0.75_minSize3_merge_ME_0.95_1022/kME_table_.csv"))
                 ) 
rei_names <- c("Gabitto_DFC", "Jorstad_DFC", "Gabitto_MTG", "Jorstad_MTG")

for(i in 1:4){
    mod_eig <- mod_eig_list[[i]]
    kme <- kme_list[[i]]
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
            cluster_rows = cluster_list[[i]],
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
            column_dend_height = unit(2, "cm"),
            show_heatmap_legend = F
            )

    pdf(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, paste0("/panel_C_", rei_names[i], ".pdf")), width = 16, height = 8)
    draw(p, padding = unit(c(6, 6, 6, 24), "mm"))#, 
         #heatmap_legend_side = "bottom")
    decorate_column_dend("Eigenmodule", {
    grid.yaxis(gp = gpar(fontsize = 8)) 
    })
    # decorate_row_dend("Eigenmodule", {
    #    grid.xaxis(gp = gpar(fontsize = 8)) 
    # })
    dev.off()

    svg(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_5/"), version_folder, paste0("/panel_C_", rei_names[i], ".svg")), width = 16, height = 8)
    draw(p, padding = unit(c(6, 6, 6, 24), "mm"))#, 
         #heatmap_legend_side = "bottom")
    decorate_column_dend("Eigenmodule", {
    grid.yaxis(gp = gpar(fontsize = 8)) 
    })
    dev.off()
}


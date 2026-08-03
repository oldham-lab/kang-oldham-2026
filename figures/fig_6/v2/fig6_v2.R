library(qs)
library(data.table)
library(scales)
library(ggplot2)
library(cowplot)
library(eulerr)
library(tidyverse)
library(ComplexHeatmap)
library(ggplotify)
library(showtext)
showtext_auto()
options(bitmapType = 'cairo')

# Panel A: Schematic (see v1)
# Panel B: Snapshots
# - Expr line plot, enrichment plot
# - Indices 

# Load bulk-related objects
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table=F, file = file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
mod_seed <- qread(file.path(module_output_dir, "modules/unmerged_modules.qs"))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]
mod_seed <- mod_seed[these_mods]
modulelengths <- unlist(lapply(mods,length))
modulelengths_seed <- unlist(lapply(mod_seed,length))


# Calc % VE (seed)
expr_t <- t(bulk_expr[,3:ncol(bulk_expr)])
colnames(expr_t) <- bulk_expr[,2]
expr_z <- apply(expr_t, 2, function(x) (x-mean(x))/sd(x))
r2_vec_bc_list <- list()
for(j in seq_along(mods)){
  expr_plot <- as.data.frame(expr_z[ ,colnames(expr_z) %in% mods[[j]]])
  pc1_bc <-  prcomp(expr_plot[apply(expr_plot,2,var)>0], scale=T)$x[,1]
  r2_vec_bc <- apply(expr_plot[apply(expr_plot,2,var)>0],2,function(x) summary(lm(x~pc1_bc))$r.squared)
  r2_vec_bc_list[[j]] <- r2_vec_bc
 # cat(j, " ")
}
vevec <- lapply(r2_vec_bc_list, mean) |> unlist()

# Load plot objects
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/")
all_plots_full = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_seed.qs")),
                 qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs")))
plots_bc_full <- list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs")))

# Tweak aesthetics of module graphs
cols2 <- RColorBrewer::brewer.pal(10, "Spectral")
newtheme <- theme_light() + theme(text = element_text(family = "sans"),
                                  axis.text.x = element_text(size = 14), 
                                  axis.text.y = element_text(size = 14), 
                                  axis.title.x = element_text(size = 14), 
                                  axis.title.y = element_text(size = 14), 
                                  plot.subtitle = element_blank(),
                                  legend.key.size = unit(0.2, "cm"), 
                                  legend.title = element_blank())

all_plots <- list()                                  
all_plots[[1]] <- mapply(\(i, j){
  subt <- paste0("# of module genes: ", modulelengths[j], "\n", sprintf("%.0f", vevec[j] * 100), "% variance explained by module eigengene")
  x <- plots_bc_full[[1]][[i]]
  x$layers[[1]] <- NULL

  if(length(unique(plots_bc_full[[1]][[i]]$data$gene)) >= 10){
    legt <- "Top 10 genes:"
  } else {
    legt <- "Genes:"
  }

  x + 
    geom_line(linewidth = 0.2) +
    labs(x = "Bulk sample (n = 100 random samples)") + 
    newtheme + 
    theme(text = element_text(size = 16),
          legend.box.margin = margin(0, 0, 0, -10),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "right") +
    scale_color_manual(values = cols2) +
    labs(title = paste0("Module ", j),
         subtitle = subt,
         colour = legt
         ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))
}, seq_along(plots_bc_full[[1]]), seq_along(mods), SIMPLIFY = F)
 
all_plots[[2]] <- lapply(all_plots_full[[2]], \(x){
  x2 <- x$data |> 
    slice(1:5) |>
    dplyr::mutate(SetName_breaks = factor(SetName_breaks, levels=rev(unique(SetName_breaks)))) |>
    ggplot(aes(x = SetName_breaks, y = pval)) +
        newtheme + 
        geom_bar(stat="identity") 
  x2$layers[[2]] <- x$layers[[2]]
  # Scale y-label size by length of label
  ysizevar <- (max(nchar(as.character(x2$data$SetName))) - 31)/58 * -4 + 10
  ysizevar <- max(ysizevar, 6)
  ysizevar <- min(ysizevar, 10)

  x2 <- x2 + coord_flip() +
    theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5, size = 10, color = "black"),
          axis.text.y = element_text(hjust = 1, vjust = 0.5, size = ysizevar, color = "black"),
          axis.title.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 11, color = "black"),
          plot.title = element_blank()) + #,
          #    plot.margin = margin(0,0,0,-10)) +
    labs(x = "", y = bquote(-log[10]~"(p-val)"), title = "Geneset enrichment") #+
  return(x2)
})

# Load indices (native, REI) and SE
base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")
out_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/")
index_dirs <- c("log_native", "log_REI")
index_list <- list()
index_save_names <- c("native_log", "REI")
index_xaxis_names <- c("Mean expression (log UMI counts)", "Relative expression index")
#for(d in 1:2){
d <- 1 # native_log

# Gabitto + Liu shared dCoPA results (for highlighting grouped barplots)
comparison <- "AllADVsCon_DFC"
dcopa_shared <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared/gabitto_vs_liu_DFC_AllADVsCon/shared_output_table.csv"))

# Gabitto control indices
gabcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv")) |>
  dplyr::select(!module)
gabcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) |>
  select(colnames(gabcon))

# Gabitto AD indices
gabad <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Alzheimers_topmodposbc_mean.csv")) |>
  dplyr::select(!module)
gabad_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassADtopmodposbc.csv")) |>
  select(colnames(gabad))

# Liu control indices
mitcon <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome_AllADVsCon_PFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
  dplyr::select(!module)
mitcon_se <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome_AllADVsCon_PFC/sn_proj_indices/", index_dirs[d], "/indices_se_RNA.SubclassCon.csv")) |>
  select(colnames(mitcon))

# Liu AD indices
mitad <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome_AllADVsCon_PFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Alzheimers.csv")) |>
  dplyr::select(!module)
mitad_se <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome_AllADVsCon_PFC/sn_proj_indices/", index_dirs[d], "/indices_se_RNA.SubclassAD.csv")) |>
  select(colnames(mitcon))

allcts1 <- colnames(gabcon)
allcts2 <- colnames(mitcon)

# Plot projection indices
projplotlist1 <- lapply(unique(dcopa_shared$mod), \(i){ # Gabitto
  # Fetch projection index data for the module
  plotdf1 <- data.frame("ct" = rep(allcts1, 2),
                       "dataset" = c(rep("Gabitto et al. 2024 (DFC)", length(allcts1))),
                       "ind" = c(unlist(gabcon[i, ]), 
                                 unlist(gabad[i, ])),
                       "ind_se" = c(unlist(gabcon_se[i, ]), 
                                    unlist(gabad_se[i, ])),
                       "dx" = c(rep("Con", length(allcts1)),
                                rep("AD", length(allcts1)))) |>
    dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                  ct = factor(ct, levels = allcts1))

  # Fetch dCoPA significance data for given module
  highlight_mat1 <- dcopa_shared |> 
    filter(mod == i, sig_FDR_Gabitto_2024 == T) |>
    select(Celltype_Gabitto_2024, P_value_Gabitto_2024) |>
    filter(!duplicated(Celltype_Gabitto_2024))
  
  # Plot grouped barplot of projection indices with asterisks for dCoPA significance
  plotdf1 <- plotdf1 |>
    left_join(highlight_mat, by = join_by(ct == Celltype_Gabitto_2024)) |>
    mutate(
      label = case_when(
        P_value_Gabitto_2024 > 0.05 ~ "",
        (0.01 < P_value_Gabitto_2024 & P_value_Gabitto_2024 <= 0.05) ~ "*",
        (0.001 < P_value_Gabitto_2024 & P_value_Gabitto_2024 <= 0.01) ~ "**",
        (0.0001 < P_value_Gabitto_2024 & P_value_Gabitto_2024 <= 0.001) ~ "***",
        P_value_Gabitto_2024 <= 0.0001 ~ "****",
        .default = NA
      )) |>
    select(!P_value_Gabitto_2024)

  plotdf2 <- data.frame(#"ct" = rep(commonct, 4),
                        "ct" = rep(allcts2, 2),
                        "dataset" = c(rep("Liu et al. 2025 (DFC)", 2 * length(allcts2))),
                        "ind" = c(unlist(mitcon[i, ]),
                                  unlist(mitad[i, ])),
                        "ind_se" = c(unlist(mitcon_se[i, ]),
                                     unlist(mitad_se[i, ])),
                        "dx" = c(rep("Con", length(allcts2)),
                                 rep("AD", length(allcts2)))) |>
    dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                  ct = factor(ct, levels = allcts2))
  
  # Fetch dCoPA significance data for given module
  highlight_mat2 <- dcopa_shared |> 
    filter(mod == i, sig_FDR_Liu_2025 == T) |>
    select(Celltype_Liu_2025, P_value_Liu_2025) |>
    filter(!duplicated(Celltype_Liu_2025))

  plotdf2 <- plotdf2 |>
      left_join(highlight_mat2, by = join_by(ct == Celltype_Liu_2025)) |>
      mutate(
        label = case_when(
          P_value_Liu_2025 > 0.05 ~ "",
          (0.01 < P_value_Liu_2025 & P_value_Liu_2025 <= 0.05) ~ "*",
          (0.001 < P_value_Liu_2025 & P_value_Liu_2025 <= 0.01) ~ "**",
          (0.0001 < P_value_Liu_2025 & P_value_Liu_2025 <= 0.001) ~ "***",
          P_value_Liu_2025 <= 0.0001 ~ "****",
          .default = NA
        )) |>
    select(!P_value_Liu_2025)

  plotdf <- rbind(plotdf1, plotdf2)  

  # Remove duplicate labels (assign label to whichever dx has highest value)
  plotdf_label <- plotdf |> 
    filter(!is.na(label)) |>
    group_by(ct, dataset) |>
    slice_max(ind, n = 1)

  dodge_width <- 0.3
  #plot_max_y <- plotdf$ind[which.max(plotdf$ind)] + plotdf$ind_se[which.max(plotdf$ind)] * 2
  p <- ggplot(plotdf, aes(x = ct, y = ind, fill = dx)) +
    theme_classic() +
    geom_col(position = position_dodge(width = dodge_width), alpha = 0.5) +
    geom_text(data = plotdf_label, 
              aes(x = ct, y = ind, label = label), 
              vjust = -0.5#,
              #y = plot_max_y + 0.1
              ) + 
    geom_errorbar(aes(ymin = ind - 2 * ind_se,
                      ymax = ind + 2 * ind_se),
                      width = 0.2,
                      linewidth = 0.3,
                      position = position_dodge(width = dodge_width)) +
    scale_fill_manual(values = c("#90D5FF", "#FFA500")) +
    theme(text = element_text(family = "sans", color = "black", size = 12),
          legend.position = "bottom", 
          legend.title = element_blank(),
          legend.margin = margin(-1, 0, 0, 0, "cm"),
          axis.title.y = element_text(size = 10, margin = margin(0, 5, 0, 0)),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(color = "black"),
          strip.background = element_rect(fill = "white"),
          plot.margin = margin(5, 0, 0, 0)) +
    labs(y = index_xaxis_names[d], x = "") +
    scale_y_continuous(breaks = c(0, 0.5, 1)#, 
                       #limits = c(0, plot_max_y + 0.2)
                       ) +
    facet_wrap(~dataset, nrow = 2, ncol = 1, scales = "free")

    return(p)
})

## Plot boxplots comparing individual module gene expression changes
# Gabitto mean expr by celltype tables (for boxplots)
gab_sum <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_summary_tables/sn_summary_objects_log.qs"))$mean$Subclass
liu_sum <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_PFC/sn_summary_tables/sn_summary_objects_log.qs"))$mean$Subclass

boxplotlist_gab <- lapply(unique(dcopa_shared$mod), \(i){
  highlight_mat <- dcopa_shared |> 
    filter(mod == i, sig_FDR_Liu_2025 == T) |>
    select(Celltype_Gabitto_2024, P_value_Gabitto_2024)
  summary_mat <- gab_sum

  p <- mapply(\(x, y){
    df <- x[rownames(x) %in% mods[[which(these_mods == i)]], colnames(x) %in% highlight_mat[,1], drop = F] 
    colnames(df) <- gsub(" ", "_", colnames(df))
    out <- data.frame("type" = y, "gene" = rownames(df), df) |>
      pivot_longer(!type:gene, names_to = "ct", values_to = "vals") |>
      mutate(ct = gsub("_", " ", ct))
    return(out)
  }, summary_mat[1:2], names(summary_mat)[1:2], SIMPLIFY = F) |>
        do.call(what = "rbind") |>
        ggplot(aes(x = type, y = vals)) +
          theme_classic() + 
          geom_point(color = "lightgrey", alpha = 0.7, size = 0.4) + 
          geom_line(aes(group = gene), color = "lightgrey", alpha = 0.7, linewidth = 0.2) +
          geom_boxplot(aes(color = type), linewidth = 0.4, notch = F, outlier.size = 0.4, alpha = 0.7, fill = NA) + 
          facet_wrap(~ct, nrow = 1) + 
          labs(x = "", y = "Mean expression\n(log UMI counts)") +
          scale_color_manual(values = c("#90D5FF", "#FFA500")) +
          theme(legend.position = "none", 
                axis.text.x = element_text(size = 12, angle = 30, hjust = 1, vjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
                #strip.text = element_text(size = 12, margin = margin(-5, 0, -5, 0))
                ) 
  return(p)
})

boxplotlist_liu <- lapply(unique(dcopa_shared$mod), \(i){
  highlight_mat <- dcopa_shared |> 
    filter(mod == i, sig_FDR_Liu_2025 == T) |>
    select(Celltype_Liu_2025, P_value_Liu_2025)
  summary_mat <- liu_sum

  p <- mapply(\(x, y){
    df <- x[rownames(x) %in% mods[[which(these_mods == i)]], colnames(x) %in% highlight_mat[,1], drop = F] 
    colnames(df) <- gsub(" ", "_", colnames(df))
    out <- data.frame("type" = y, "gene" = rownames(df), df) |>
      pivot_longer(!type:gene, names_to = "ct", values_to = "vals") |>
      mutate(ct = gsub("_", " ", ct))
    return(out)
  }, summary_mat[1:2], names(summary_mat)[1:2], SIMPLIFY = F) |>
        do.call(what = "rbind") |>
        ggplot(aes(x = type, y = vals)) +
          theme_classic() + 
          geom_point(color = "lightgrey", alpha = 0.7, size = 0.4) + 
          geom_line(aes(group = gene), color = "lightgrey", alpha = 0.7, linewidth = 0.2) +
          geom_boxplot(aes(color = type), linewidth = 0.4, notch = F, outlier.size = 0.4, alpha = 0.7, fill = NA) + 
          facet_wrap(~ct, nrow = 1) + 
          labs(x = "", y = "Mean expression\n(log UMI counts)") +
          scale_color_manual(values = c("#90D5FF", "#FFA500")) +
          theme(legend.position = "none", 
                axis.text.x = element_text(size = 12, angle = 30, hjust = 1, vjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
                #strip.text = element_text(size = 12, margin = margin(-5, 0, -5, 0))
                ) 
  return(p)
})

# Construct and save module snapshot

out_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/dCoPA_snapshots/"), comparison)
if(!dir.exists(out_dir))
  dir.create(out_dir, recursive = T)

for(i in seq_along(unique(dcopa_shared$mod))){
  boxplot_combined <- plot_grid(boxplotlist_gab[[i]], boxplotlist_liu[[i]], nrow = 1)
  pall <- suppressMessages(plot_grid(all_plots[[1]][[which(these_mods == unique(dcopa_shared$mod)[i])]],
                                     all_plots[[2]][[which(these_mods == unique(dcopa_shared$mod)[i])]],
                                     projplotlist1[[i]],
                                     boxplot_combined,
                                     nrow = 4, 
                                     rel_heights = c(0.5, 0.3, 1, 0.5)))
  ggsave(pall, 
         file = file.path(out_dir, paste0(unique(dcopa_shared$mod)[i], ".svg")), 
         device = svglite::svglite, 
         width = 6, 
         height = 12, 
         bg = "white")
}

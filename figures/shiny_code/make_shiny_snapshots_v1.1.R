# Lein DFC
# MIT AD Multiome/Multiregion
# brainSCOPE
# starting server: sudo shiny-server

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


# Load module graphs
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
all_plots_full = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_seed.qs")),
                 qread(file.path(save_dir1,"sn_proj_objects","bulkcor_object.qs")),
                 qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs")))
plots_bc_full <- list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs")))

# for indexing and file names
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
modulelengths_seed <- unlist(lapply(mod_seed, length))
mod_id <- 1:1158
mod_index <- lapply(mod_id, \(x) which(these_mods %in% x)) |> unlist()

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
  #subt <- as.list(as.character(all_plots[[1]][[i]]$labels$subtitle)) |> do.call(what = paste0)
  subt <- paste0(modulelengths_seed[these_mods[j]], " seed genes | ", modulelengths[these_mods[j]], " topmodposbc genes", "\n", sprintf("%.0f", vevec[these_mods[j]] * 100), "% variance explained by module eigengene")
  x <- all_plots_full[[1]][[i]]
  x$layers[[1]] <- NULL
  if(length(unique(all_plots_full[[1]][[i]]$data$gene)) >= 10){
    legt <- "Top 10 genes:"
  } else {
    legt <- "Genes:"
  }
  x + 
    geom_line(linewidth = 0.2) +
    labs(x = "Human dorsolateral prefrontal cortex samples\n(n = 100 random samples)") + 
    newtheme + 
    theme(text = element_text(size = 16),
          legend.box.margin = margin(0, 0, 0, -10),
          #plot.title = element_text(size = 14, hjust = 0.5),# margin = margin(0, 0, -0.1, 0)),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "right") +
    scale_color_manual(values = cols2) +
    labs(#title = paste0("Module ", these_mods[i]),
         title = paste0("Module ", i), # Numbering post-filter
         subtitle = subt,
         colour = legt
         #subtitle = bquote(subt)
         ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))#, 
                                #nrow = 5, 
                                #theme = theme(legend.byrow = TRUE)))
}, seq_along(all_plots_full[[1]]), mod_index, SIMPLIFY = F)

plots_bc <- list()
plots_bc[[1]] <- mapply(\(i, j){
    #subt <- as.list(as.character(all_plots[[1]][[i]]$labels$subtitle)) |> do.call(what = paste0)
  subt <- paste0("# of module genes: ", modulelengths[these_mods[j]], "\n", sprintf("%.0f", vevec[these_mods[j]] * 100), "% variance explained by module eigengene")
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
          #plot.title = element_text(size = 14, hjust = 0.5),# margin = margin(0, 0, -0.1, 0)),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "right") +
    scale_color_manual(values = cols2) +
    labs(#title = paste0("Module ", these_mods[i]),
         title = paste0("Module ", i),
         subtitle = subt,
         colour = legt
         #subtitle = bquote(subt)
         ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))#, 
                                #nrow = 5, 
                                #theme = theme(legend.byrow = TRUE)))
}, seq_along(plots_bc_full[[1]]), mod_index, SIMPLIFY = F)

all_plots[[2]] <- lapply(all_plots_full[[2]], \(x){
  x + newtheme +
    theme(legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.margin = margin(-10, 0, 0, 0),
          plot.margin = margin(1, 0.3, 0.6, 0.3, "cm")) +
    guides(fill = guide_legend(override.aes = list(shape = 22)))
})

all_plots[[3]] <- lapply(all_plots_full[[3]], \(x){
  x2 <- x$data |> 
    #dplyr::filter(pval > 0) |>
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


# mean, median
# REI, log native
# seed, topmodposbc

# Load indices (native, REI) and SE
base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")
index_dirs <- c("log_native", "log_REI")
index_save_names <- c("native_log", "REI")
index_xaxis_names <- c("Mean expression (log UMI counts)", "Relative expression index")
for(def in c(#"seed"#, 
             "topmodposbc"
             )){
  cat(def, "\n")
  for(sum_type in c("mean", 
                    "median"
                   )){
    for(d in 1:2){ # native log vs REI
      lein <- fread(data.table = F, file = file.path(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], paste0("/indices_over_all_datasets_Cell_Type_1_", def, "_", sum_type,".csv"))) |>
        dplyr::select(!module)
      SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], paste0("/indices_over_all_datasets_Subclass_Control_", def, "_", sum_type,".csv"))) |>
        dplyr::select(!module)
      mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], paste0("/indices_over_all_datasets_Subclass_Control_", def, "_", sum_type,".csv"))) |>
        dplyr::select(!module)

      lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_se_", def, ".csv")) 
      SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon", def,".csv")) 
      mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon", def, ".csv")) 

      lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]  
      SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]
      mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]


      # Change celltypes to match

      colnames(lein)[15] <- "Micro"
      colnames(mora)[c(1, 3, 14, 15)] <- c("Astro", "Endo", "Micro", "Oligo")
      colnames(SEAcon)[c(1, 3, 15, 16)] <- c("Astro", "Endo", "Micro", "Oligo")

      colnames(lein_se) <- colnames(lein)
      colnames(SEAcon_se) <- colnames(SEAcon)
      colnames(mora_se) <- colnames(mora)

      allcts <- c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
                "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip",
                "Astro", "Endo", "Micro", "Oligo", "OPC", "VLMC")
      
      # Fill out celltypes
      lein[setdiff(allcts, colnames(lein))] <- NA
      lein <- lein[, match(allcts, colnames(lein))]
      lein_se[setdiff(allcts, colnames(lein_se))] <- NA
      lein_se <- lein_se[, match(allcts, colnames(lein_se))]

      SEAcon[setdiff(allcts, colnames(SEAcon))] <- NA
      SEAcon <- SEAcon[, match(allcts, colnames(SEAcon))]
      SEAcon_se[setdiff(allcts, colnames(SEAcon_se))] <- NA
      SEAcon_se <- SEAcon_se[, match(allcts, colnames(SEAcon_se))]

      mora[setdiff(allcts, colnames(mora))] <- NA
      mora <- mora[, match(allcts, colnames(mora))]
      mora_se[setdiff(allcts, colnames(mora_se))] <- NA
      mora_se <- mora_se[, match(allcts, colnames(mora_se))]

      # Plot projection indices
      projplotlist <- list()
      for(i in 1:nrow(lein)){
        n <- length(allcts)
        plotdf <- data.frame("ct" = rep(allcts, 3),
                            "dataset" = c(rep("Jorstad et al. 2023 (normal human dorsolateral prefrontal cortex)", n), 
                                          rep("Gabitto et al. 2024 (normal human dorsolateral prefrontal cortex)", n), 
                                          rep("Morabito et al. 2021 (normal human dorsolateral prefrontal cortex)", n)),
                            "ind" = c(unlist(lein[i, ]), 
                                      unlist(SEAcon[i, ]), 
                                      unlist(mora[i, ])),
                            "ind_se" = c(unlist(lein_se[i, ]), 
                                          unlist(SEAcon_se[i, ]), 
                                          unlist(mora_se[i, ]))) |>
          dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                        ct = factor(ct, levels = allcts))

        projplotlist[[i]] <- ggplot(plotdf, aes(x = ct, y = ind, fill = ct)) +
          theme_classic() +
          geom_col(position = position_dodge(), alpha = 0.5) +
          geom_errorbar(aes(ymin = ind - 2 * ind_se,
                            ymax = ind + 2 * ind_se),
                            width = 0.2,
                            linewidth = 0.3,
                            position = position_dodge(0.5)) +
          scale_fill_manual(values = hue_pal()(n)) +
          theme(text = element_text(family = "sans", color = "black", size = 14),
                legend.position="none", 
                axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
                axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                strip.text = element_text(color = "black"),
                strip.background = element_rect(fill = "white")) +
          labs(y = index_xaxis_names[d], x = "")
          
          if(d == 1){
            projplotlist[[i]] <- projplotlist[[i]] +       
              facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free") 
          } else {
            projplotlist[[i]] <- projplotlist[[i]] +       
              facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_x") 
          }
      }

      # Plot individual cor heatmaps for each module
      corlistind <- list()
      for(i in 1:nrow(lein)){
        tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), 
                            "Gabitto 2024" = unlist(SEAcon[i, ]), 
                            "Morabito 2021" = unlist(mora[i, ]))
        tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
        colnames(tempcor) <- c("Jorstad 2023", "Gabitto 2024", "Morabito 2021")
        rownames(tempcor) <- c("Jorstad 2023", "Gabitto 2024", "Morabito 2021")
        testvec <- apply(tempcor, 2, \(x) sum(is.na(x)))
        if(any(testvec > 0)){
          rem <- which(testvec > 0)
          tempcor <- tempcor[-rem, -rem]
        } 
        col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        p <-  Heatmap(tempcor, 
                      column_title = "Correlations of module projections",
                      column_title_gp = gpar(fontsize = 14),
                      heatmap_legend_param = list(
                        direction = "horizontal",
                        title_gp = gpar(fontface = "plain"),
                        title_position = "topcenter"
                      ),
                      name = "Correlation", 
                      column_names_rot = 90,  
                      col = col_fun,    
                      row_names_gp = gpar(fontsize = 12),
                      row_names_side = "left",
                      row_dend_side = "right",
                      width = unit(3, "cm"), 
                      height = unit(3, "cm"),
                      show_column_names = F)
        corlistind[[i]] <- as.ggplot(p)
        
      #  cat(i, " ")
      }

      # Arrange snapshots and save as SVGs
      out_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/shiny_code/shiny_snapshots/"), index_save_names[d], def, sum_type)
      if(!dir.exists(out_dir)){
        dir.create(out_dir, recursive = T)
      }
      for(j in seq_along(all_plots[[1]])){
        rightside <- plot_grid(projplotlist[[these_mods[j]]], corlistind[[these_mods[j]]], nrow = 2, rel_heights = c(3, 1))
    
        if(def == "seed"){
          leftsideall <- suppressMessages(plot_grid(all_plots[[1]][[j]],
                                                    all_plots[[2]][[j]],
                                                    all_plots[[3]][[j]],
                                                    nrow = 3, rel_heights = c(1,1,1)))
        } else {
          leftsideall <- suppressMessages(plot_grid(plots_bc[[1]][[j]],
                                                    all_plots[[2]][[j]],
                                                    all_plots[[3]][[j]],
                                                    nrow = 3, rel_heights = c(1,1,1)))
        }

        p <- suppressMessages(plot_grid(leftsideall, rightside, ncol=2))
        ggsave(p, file = file.path(out_dir, paste0(j, ".svg")), device = svglite::svglite, width = 13, height = 9, bg = "white")
      }
      cat(def, sum_type, index_dirs[d], "done\n")
      timestamp()
    }
  }
}
 
# ## Summarize the results:
# # Find cors and plot as boxplots
# corlist <- list()
# for(i in 1:nrow(lein)){
#   tempdf <- data.frame("lein" = unlist(lein[i, ]), 
#                        "mit" = unlist(mit[i, ]), 
#                       # "scope" = unlist(scope[i, ]), 
#                        "SEA" = unlist(SEAcon[i, ]), 
#                        "mora" = unlist(mora[i, ]))
#   tempcor <- cor(tempdf)
#   corlist[[i]] <- data.frame("mod" = i, "cors" = tempcor[upper.tri(tempcor)])
# }
# cordf <- do.call(rbind, corlist)
# cordforder <- cordf |> group_by(mod) |> summarize(mean = mean(cors)) |> arrange(mean)
# cordf$mod <- factor(cordf$mod, levels = cordforder$mod)

# p <- ggplot(cordf, aes(x = mod, y = cors)) + 
#   theme_light() +
#   geom_boxplot(linewidth = 0.1, outlier.shape = NA, width = 0.2) +
#   labs(x = "Module", y = "Pairwise correlation", title = "Pairwise correlations between 5 dataset projections") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/5datasetproj_pairwiseCors.svg"), width = 12, height = 5)

# # Plot as means (points) plus error bars instead
# cordfmean <- cordf |> group_by(mod) |> summarize(mean = mean(cors),
#                                                  se = sd(cors)/sqrt(n())) |> arrange(mean)
# cordfmean$mod <- factor(cordfmean$mod, levels = cordfmean$mod)
# p <- ggplot(cordfmean, aes(x = mod, y = mean)) + 
#   theme_light() +
#   geom_pointrange(aes(ymin = mean - 2*se, ymax = mean + 2*se), linewidth = 0.2) +
#   labs(x = "Module (n = 1023)", y = "Mean pairwise correlation", title = "Pairwise correlations between 4 dataset projections") +
#   theme(text = element_text(size = 20),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/4datasetproj_pairwiseCors_meanSE.svg"), width = 12, height = 5)


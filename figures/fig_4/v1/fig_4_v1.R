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
# Create snapshots for 4 modules



# Load bulk-related objects
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024")
filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
modulelengths <- unlist(lapply(mods,length))
modulelengths_seed <- unlist(lapply(mod_seed,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
#mod_id <- 20
mod_id <- c(347, # oligo
            635, # VIP although VIP is not a top 10 seed gene
            822, # Endo
            253 # mystery module
            )
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

# Load plot objects
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
all_plots_full = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_seed.qs"))[mod_index],
                 qread(file.path(save_dir1,"sn_proj_objects","bulkcor_object.qs"))[mod_index],
                 qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs"))[mod_index])
plots_bc_full <- list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs"))[mod_index])


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
  subt <- paste0("# of module genes: ", modulelengths[j], "\n", sprintf("%.0f", vevec[j] * 100), "% variance explained by module eigengene")
  x <- all_plots_full[[1]][[i]]
  x$layers[[1]] <- NULL

  if(length(unique(all_plots_full[[1]][[i]]$data$gene)) >= 10){
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
          plot.title = element_blank(),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "right") +
    scale_color_manual(values = cols2) +
    labs(subtitle = subt,
         colour = legt
         #subtitle = bquote(subt)
         ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))#, 
                                #nrow = 5, 
                                #theme = theme(legend.byrow = TRUE)))
}, seq_along(all_plots_full[[1]]), mod_id, SIMPLIFY = F)
plots_bc <- list()
plots_bc[[1]] <- mapply(\(i, j){
    #subt <- as.list(as.character(all_plots[[1]][[i]]$labels$subtitle)) |> do.call(what = paste0)
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
          #plot.title = element_text(size = 14, hjust = 0.5),# margin = margin(0, 0, -0.1, 0)),
          plot.title = element_blank(),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "right") +
    scale_color_manual(values = cols2) +
    labs(subtitle = subt,
         colour = legt
         #subtitle = bquote(subt)
         ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))#, 
                                #nrow = 5, 
                                #theme = theme(legend.byrow = TRUE)))
}, seq_along(plots_bc_full[[1]]), mod_id, SIMPLIFY = F)
all_plots[[2]] <- lapply(all_plots_full[[2]], \(x){
  x$data <- x$data |> mutate(split = case_match(split, 
    "Brainseq (n=182)" ~ "Brainseq",
    "BrainGVEX (n=255)" ~ "BrainGVEX",
    "CMC (n=274)" ~ "CMC",
    "CMC_HBCC (n=145)" ~ "CMC_HBCC",
    "GTEx (n=190)" ~ "GTEx",
    "NABEC (n=63)" ~ "NABEC",
    "ROSMAP (n=347)" ~ "ROSMAP",
    .default = split
  ))
  return(x + newtheme +
    theme(text = element_text(color = "black"),
          axis.text.x = element_text(color =  "black", angle = 30, hjust = 1, vjust = 1),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.margin = margin(-10, 0, 0, 0),
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 14),
          plot.title = element_blank(),
          plot.margin = margin(1,1,1,1,"cm")) +
    guides(fill = guide_legend(override.aes = list(shape = 22))))
})

all_plots[[3]] <- lapply(all_plots_full[[3]], \(x){
  x2 <- x$data |> 
    slice(1:5) |>
    #dplyr::filter(pval > 0) |>
    dplyr::mutate(SetName=factor(SetName, levels=rev(unique(SetName)))) |>
    ggplot(aes(x = SetName, y = pval)) +
        newtheme + 
        geom_bar(stat="identity") #+
        # theme(axis.text.x = element_text(hjust=1,vjust = 0.5, size = 10),
        #       axis.text.y = element_text(hjust=1,vjust = 0.5, size = 10),
        #       axis.title.y = element_text(size = 11),
        #       axis.title.x = element_text(size = 11),
        #       plot.title = element_blank()) + #,
        #   #    plot.margin = margin(0,0,0,-10)) +
        # labs(x = "", y = bquote(-log[10]~"(p-val)"), title = "Geneset enrichment") #+
        #scale_x_discrete(limits = rev(levels(SetName))) +
        #geom_hline(yintercept = gsea_cutFDR, color = "red")
  x2$layers[[2]] <- x$layers[[2]]
  x2 <- x2 + coord_flip() +
    theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5, size = 10, color = "black"),
          axis.text.y = element_text(hjust = 1,vjust = 0.5, size = 10, color = "black"),
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
    d=2

  lein <- fread(data.table = F, file = paste0(base_dir, "/LeinMTG/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_con.csv")) |>
    dplyr::select(!module)
  lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinMTG/sn_proj_indices/", index_dirs[d], "/indices_se_Cell_Type.csv")) 
  lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]


  mit <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
    dplyr::select(!module)
  mit_se <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_se_RNA.SubclassCon.csv")) 
  mit_se <- mit_se[, match(colnames(mit), colnames(mit_se))]

  SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
    dplyr::select(!module)
  SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon.csv")) 
  SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]

  mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
    dplyr::select(!module)
  mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon.csv")) 
  mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]

  # Change celltypes to match

  colnames(lein)[15] <- "Micro"
  colnames(mit) <- gsub("Exc ", "", colnames(mit))
  colnames(mit) <- gsub("Inh ", "", colnames(mit))
  colnames(mit) <- gsub("-", "/", colnames(mit))
  colnames(mit)[c(2, 3, 4, 6, 11, 12, 16, 21, 25, 27)] <- c("Micro", "Vip", "Astro", "Oligo", "Pvalb", "Endo", "Sst", "Lamp5", "Pax6", "L6 IT Car3")
  colnames(mora)[c(1, 3, 14, 15)] <- c("Astro", "Endo", "Micro", "Oligo")
  colnames(SEAcon)[c(1, 3, 15, 16)] <- c("Astro", "Endo", "Micro", "Oligo")

  colnames(lein_se) <- colnames(lein)
  colnames(mit_se) <- colnames(mit)
  colnames(SEAcon_se) <- colnames(SEAcon)
  colnames(mora_se) <- colnames(mora)

  # Find intersection of celltypes and create euler diagram
  commonct <- Reduce(intersect, list(colnames(lein), colnames(mit), 
                   # colnames(scope), 
                    colnames(SEAcon), colnames(mora)))
  allcts <- unique(c(colnames(lein), colnames(mit), 
                    #colnames(scope), 
                    colnames(SEAcon), colnames(mora)))

  # Hardcode order of cts
  #  "excitatory neurons, inhibitory neurons, astrocytes, oligodendrocytes, OPCs, 
  #  microglia, vascular cells (including fibroblasts), T cells."
  # (remove EC - entorhinal cortex-specific celltype from Liu 2025)
  allcts <- c("L2/3 IT", "L3/4 IT", "L3/5 IT", "L4 IT", "L4/5 IT/1", "L4/5 IT/2", "L5 ET", "L5 IT",  "L5/6 IT", "L5/6 NP", 
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", 
              "Sst", "Vip", "Sst Chodl", "Astro", "Oligo", "OPC", "Micro", "Endo", "Fib", "Per", "SMC", "VLMC", "T")

 
  # Fill out celltypes
  lein[setdiff(allcts, colnames(lein))] <- NA
  lein <- lein[, match(allcts, colnames(lein))]
  lein_se[setdiff(allcts, colnames(lein_se))] <- NA
  lein_se <- lein_se[, match(allcts, colnames(lein_se))]

  mit[setdiff(allcts, colnames(mit))] <- NA
  mit <- mit[, match(allcts, colnames(mit))]
  mit_se[setdiff(allcts, colnames(mit_se))] <- NA
  mit_se <- mit_se[, match(allcts, colnames(mit_se))]

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
  for(i in mod_id){
 # i = mod_id
    #n <- length(commonct)
    n <- length(allcts)
    plotdf <- data.frame(#"ct" = rep(commonct, 4),
                         "ct" = rep(allcts, 4),
                         "dataset" = c(rep("Jorstad et al. 2023 (MTG)", n), rep("Liu et al. 2025 (control MTC)", n), 
                                       #rep("brainSCOPE", n), 
                                       rep("Gabitto et al. 2024 (control MTG)", n), rep("Morabito et al. 2021 (control PFC)", n)),
                         "ind" = c(unlist(lein[i, ]), unlist(mit[i, ]), 
                                   #unlist(scope[i, ]), 
                                   unlist(SEAcon[i, ]), unlist(mora[i, ])),
                         "ind_se" = c(unlist(lein_se[i, ]), unlist(mit_se[i, ]), 
                                      #unlist(scope_se[i, ]), 
                                      unlist(SEAcon_se[i, ]), unlist(mora_se[i, ]))) |>
      dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                    ct = factor(ct, levels = allcts))

    projplotlist[[which(mod_id == i)]] <- ggplot(plotdf, aes(x = ct, y = ind, fill = ct)) +
      theme_classic() +
      geom_col(position = position_dodge(), alpha = 0.5) +
      geom_errorbar(aes(ymin = ind - 2 * ind_se,
                        ymax = ind + 2 * ind_se),
                        width = 0.2,
                        linewidth = 0.3,
                        position = position_dodge(0.5)) +
      scale_fill_manual(values = hue_pal()(n)) +
      theme(text = element_text(family = "sans", color = "black", size = 12),
            legend.position="none", 
            axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
            axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
            strip.text = element_text(color = "black"),
            strip.background = element_rect(fill = "white")) +
      facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_y") +
      labs(y = index_xaxis_names[d], x = "") +
      scale_y_continuous(breaks = c(0, 0.5, 1))
      
      # if(d == 1){
      #   projplotlist[[i]] <- projplotlist[[i]] +       
      #     facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free") 
      # } else {
      #   projplotlist[[i]] <- projplotlist[[i]] +       
      #     facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_x") 
      # }
      cat(i, " ")
  }

  # Plot individual cor heatmaps for each module
  corlistind <- list()
  for(i in mod_id){
    tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), "Liu 2025" = unlist(mit[i, ]), 
                         #"brainSCOPE_con" = unlist(scope[i, ]), 
                         "Gabitto 2024" = unlist(SEAcon[i, ]), "Morabito 2021" = unlist(mora[i, ]))
    tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
    colnames(tempcor) <- c("Jorstad 2023",  "Liu 2025", "Gabitto 2024", "Morabito 2021")
    rownames(tempcor) <- c("Jorstad 2023",  "Liu 2025", "Gabitto 2024", "Morabito 2021")
    testvec <- apply(tempcor, 2, \(x) sum(is.na(x)))
    if(any(testvec > 0)){
      rem <- which(testvec > 0)
      tempcor <- tempcor[-rem, -rem]
    } 
    col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    p <-  Heatmap(tempcor, 
                  heatmap_legend_param = list(
                    direction = "horizontal",
                    title_gp = gpar(fontface = "plain"),
                    title_position = "topcenter"
                  ),
                  name = "Correlation", 
                  column_names_rot = 90,  
                  col = col_fun,    
                  row_names_gp = gpar(fontsize = 22),
                  row_names_side = "left",
                  row_dend_side = "right",
                  #column_names_gp = gpar(fontsize = 24),
                  width = unit(4, "cm"), 
                  height = unit(4, "cm"),
                  show_column_names = F)
    #corlistind[[which(mod_id == i)]] <- as.ggplot(p) 
    #corlistind[[which(mod_id == i)]] <- p
    
    pdf(file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/panel_"),i, "_5.pdf"), height = 4, width = 6)
    draw(p, merge_legend = TRUE, heatmap_legend_side = "bottom")
    dev.off()
    cat(i, " ")
  }

  ## Plot individual panels
outindices <- 1:length(all_plots[[1]])
for(j in outindices){
  # Expression plot:
  ggsave(plots_bc[[1]][[j]], 
        file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/panel_"),mod_id[j], "_1.pdf"),
        width = 6, height = 3)


  # Bulk cor plot:
  ggsave(all_plots[[2]][[j]],
          file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/panel_"), mod_id[j], "_2.pdf"),
          width = 6, height = 4)

  # GSEA plot:
  max_setname <- lapply(all_plots[[3]][[1]]$data$SetName |> as.character(), nchar) |> unlist() |> max()
  ggsave(all_plots[[3]][[j]], file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/panel_"), mod_id[j], "_3.pdf"),
          height = 1.5, width = (3 + max_setname/30 * 2.5))

  # Projection indices:
  ggsave(projplotlist[[j]], file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/panel_"), mod_id[j],"_4.pdf"), height = 5, width = 5.5)

  # Index correlation heatmap
  #ggsave(corlistind[[j]], file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/panel_"),mod_id[j], "_5.pdf"), height = 4, width = 6)
  
}

  
# # Arrange snapshots and save as SVGs
# out_plot <- list()
# outindices <- 1:length(all_plots[[1]])

# for(j in outindices){
#   #rightside <- plot_grid(projplotlist[[these_mods[j]]],
#   #                       align = "hv")

#   p1 <- plot_grid(all_plots[[1]][[j]],
#                   all_plots[[2]][[j]])
#   p <- suppressMessages(plot_grid(p1,
#                                    all_plots[[3]][[j]],
#                                    projplotlist[[j]], 
#                                    corlistind[[j]],
#                                    nrow = 5, rel_heights = c(1,0.6,1,1,0.7)))

#   ggsave(p, file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/"), j, ".pdf"), width = 7, height = 15, bg = "white")
#   cat(j, " ")
# }









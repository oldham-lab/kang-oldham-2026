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
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024")
all_plots = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_seed.qs")),
                 qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs")))

# for indexing and file names
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024")
filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])


# Tweak aesthetics of module graphs
cols2 <- RColorBrewer::brewer.pal(10, "Spectral")
newtheme <- theme_light() + theme(text = element_text(family = "sans"),
                                  axis.text.x = element_text(size = 12), 
                                  axis.text.y = element_text(size = 12), 
                                  axis.title.x = element_text(size = 12), 
                                  axis.title.y = element_text(size = 12), 
                                  plot.subtitle = element_blank(),
                                  legend.key.size = unit(0.2, "cm"), 
                                  legend.title = element_blank())
all_plots[[1]] <- lapply(seq_along(all_plots[[1]]), \(i){
  subt <- as.list(as.character(all_plots[[1]][[i]]$labels$subtitle)) |> do.call(what = paste0)
  x <- all_plots[[1]][[i]]
  x$layers[[1]] <- NULL

  x + 
    geom_line(linewidth = 0.2) +
    labs(x = "Sample") + 
    newtheme + 
    theme(legend.box.margin = margin(0, 0, 0, -10),
          plot.title = element_text(size = 20, hjust = 0.5, margin = margin(0, 0, -0.1, 0)),
          legend.text = element_text(size = 18)) +
    scale_color_manual(values = cols2) +
    labs(title = paste0("Module ", these_mods[i]),
         subtitle = bquote(subt))
})
all_plots[[2]] <- lapply(all_plots[[2]], \(x){
  x2 <- x$data |> 
    #dplyr::filter(pval > 0) |>
    dplyr::mutate(SetName=factor(SetName, levels=rev(unique(SetName)))) |>
    ggplot(aes(x = SetName, y = pval)) +
        newtheme + 
        geom_bar(stat="identity") +
        theme(axis.text.x = element_text(hjust=1,vjust = 0.5, size = 10),
              axis.text.y = element_text(hjust=1,vjust = 0.5, size = 10),
              plot.title = element_blank()) + #,
          #    plot.margin = margin(0,0,0,-10)) +
        labs(x = "", y = bquote(-log[10](p-val)), title = "Geneset enrichment") #+
        #scale_x_discrete(limits = rev(levels(SetName))) +
        #geom_hline(yintercept = gsea_cutFDR, color = "red")
  x2$layers[[2]] <- x$layers[[2]]
  x2 <- x2 + coord_flip()
  return(x2)
})

# Load indices (native) for SEAAD2024 and MIT_Multiome
base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")

SEA <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Alzheimers.csv")) |>
  dplyr::select(!module)
SEA_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassAD.csv")) 
SEA_se <- SEA_se[, match(colnames(SEA), colnames(SEA_se))]

SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
  dplyr::select(!module)
SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon.csv")) 
SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]

mit <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Alzheimers.csv")) |>
  dplyr::select(!module)
mit_se <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_se_RNA.SubclassAD.csv")) 
mit_se <- mit_se[, match(colnames(mit), colnames(mit_se))]

mitcon <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
  dplyr::select(!module)
mitcon_se <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_se_RNA.SubclassCon.csv")) 
mitcon_se <- mitcon_se[, match(colnames(mitcon), colnames(mitcon_se))]

# Plot projection indices
projplotlist <- list()
for(i in 1:nrow(lein)){
  #n <- length(commonct)
  n <- length(allcts)
  plotdf <- data.frame(#"ct" = rep(commonct, 4),
                        "ct" = rep(allcts, 4),
                        "dataset" = c(rep("Lein DLPFC", n), rep("MIT MTC", n), 
                                      #rep("brainSCOPE", n), 
                                      rep("SEAAD2024", n), rep("Morabito 2021", n)),
                        "ind" = c(unlist(lein[i, ]), unlist(mit[i, ]), 
                                  #unlist(scope[i, ]), 
                                  unlist(SEAcon[i, ]), unlist(mora[i, ])),
                        "ind_se" = c(unlist(lein_se[i, ]), unlist(mit_se[i, ]), 
                                    #unlist(scope_se[i, ]), 
                                    unlist(SEAcon_se[i, ]), unlist(mora_se[i, ])))  
  plotdf$dataset <- factor(plotdf$dataset, levels = unique(plotdf$dataset))

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
    #facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_y") +
    labs(y = index_xaxis_names[d], x = "")
    
    if(d == 1){
      projplotlist[[i]] <- projplotlist[[i]] +       
        facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_y") 
    } else {
      projplotlist[[i]] <- projplotlist[[i]] +       
        facet_wrap(~dataset, ncol = 1, nrow = 5) 
    }
    cat(i, " ")
}

# Plot individual cor heatmaps for each module
corlistind <- list()
for(i in 1:nrow(lein)){
  tempdf <- data.frame("Lein" = unlist(lein[i, ]), "MIT_Multiome_con" = unlist(mit[i, ]), 
                        #"brainSCOPE_con" = unlist(scope[i, ]), 
                        "SEAAD2024_con" = unlist(SEAcon[i, ]), "Morabito_con" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
  testvec <- apply(tempcor, 2, \(x) sum(is.na(x)))
  if(any(testvec > 0)){
    rem <- which(testvec > 0)
    tempcor <- tempcor[-rem, -rem]
  } 
  col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  p <-  Heatmap(tempcor, 
                heatmap_legend_param = list(
                  direction = "horizontal",
                  title_gp = gpar(fontface = "plain")
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
  corlistind[[i]] <- as.ggplot(p)
  
  cat(i, " ")
}

# Arrange snapshots and save as SVGs
out_plot <- list()
outindices <- 1:length(all_plots[[1]])
if(!dir.exists(paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/module_reproducibility_"), index_save_names[d], "/"))){
  dir.create(paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/module_reproducibility_"), index_save_names[d], "/"), recursive = T)
}
for(j in outindices){
  #rightside <- plot_grid(projplotlist[[these_mods[j]]],
  #                       align = "hv")
  rightside <- plot_grid(projplotlist[[these_mods[j]]], corlistind[[these_mods[j]]], nrow=2,rel_heights=c(3,1))

  leftsideall <- suppressMessages(plot_grid(all_plots[[1]][[j]],
                                            all_plots[[2]][[j]],
                                            all_plots[[3]][[j]],
                                            nrow = 3, rel_heights = c(1,1,1)))

  p <- suppressMessages(plot_grid(leftsideall, rightside, ncol=2))
  ggsave(p, file = paste0(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/module_reproducibility_"), index_save_names[d], "/", paste0(j, ".svg")), device = svglite::svglite, width = 13, height = 9, bg = "white")
  cat(j, " ")
}




# > colnames(lein)
#  [1] "Astro"      "Chandelier" "Endo"       "L2/3 IT"    "L4 IT"     
#  [6] "L5 ET"      "L5 IT"      "L5/6 NP"    "L6 CT"      "L6 IT"     
# [11] "L6 IT Car3" "L6b"        "Lamp5"      "Lamp5 Lhx6" "Micro/PVM" 
# [16] "Oligo"      "OPC"        "Pax6"       "Pvalb"      "Sncg"      
# [21] "Sst"        "Sst Chodl"  "Vip"        "VLMC"       "module"   

# > colnames(mit)
#  [1] "Exc L2-3 IT"      "Mic"              "Inh VIP"          "Ast"             
#  [5] "Exc L3-5 IT"      "Oli"              "Exc L4-5 IT-2"    "Exc L6 IT"       
#  [9] "Exc L5/6 NP"      "Exc EC"           "Inh PVALB"        "End"             
# [13] "Exc L3-4 IT"      "OPC"              "Exc L5-6 IT"      "Inh SST"         
# [17] "Exc L4-5 IT-1"    "Fib"              "Exc L6b"          "Per"             
# [21] "Inh LAMP5"        "SMC"              "Exc L5 ET"        "Exc L6 CT"       
# [25] "Inh PAX6"         "T"                "Exc L5/6 IT Car3" "module" 

# > colnames(scope)
#  [1] "Astro"      "Chandelier" "Endo"       "Immune"     "L2/3 IT"   
#  [6] "L4 IT"      "L5 ET"      "L5 IT"      "L5/6 NP"    "L6 CT"     
# [11] "L6 IT"      "L6 IT Car3" "L6b"        "Lamp5"      "Lamp5 Lhx6"
# [16] "Micro"      "Oligo"      "OPC"        "Pax6"       "PC"        
# [21] "Pvalb"      "SMC"        "Sncg"       "Sst"        "Sst Chodl" 
# [26] "Vip"        "VLMC"       "module" 





## Summarize the results:
# Find cors and plot as boxplots
corlist <- list()
for(i in 1:nrow(lein)){
  tempdf <- data.frame("lein" = unlist(lein[i, ]), 
                       "mit" = unlist(mit[i, ]), 
                      # "scope" = unlist(scope[i, ]), 
                       "SEA" = unlist(SEAcon[i, ]), 
                       "mora" = unlist(mora[i, ]))
  tempcor <- cor(tempdf)
  corlist[[i]] <- data.frame("mod" = i, "cors" = tempcor[upper.tri(tempcor)])
}
cordf <- do.call(rbind, corlist)
cordforder <- cordf |> group_by(mod) |> summarize(mean = mean(cors)) |> arrange(mean)
cordf$mod <- factor(cordf$mod, levels = cordforder$mod)

p <- ggplot(cordf, aes(x = mod, y = cors)) + 
  theme_light() +
  geom_boxplot(linewidth = 0.1, outlier.shape = NA, width = 0.2) +
  labs(x = "Module", y = "Pairwise correlation", title = "Pairwise correlations between 5 dataset projections") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/5datasetproj_pairwiseCors.svg"), width = 12, height = 5)

# Plot as means (points) plus error bars instead
cordfmean <- cordf |> group_by(mod) |> summarize(mean = mean(cors),
                                                 se = sd(cors)/sqrt(n())) |> arrange(mean)
cordfmean$mod <- factor(cordfmean$mod, levels = cordfmean$mod)
p <- ggplot(cordfmean, aes(x = mod, y = mean)) + 
  theme_light() +
  geom_pointrange(aes(ymin = mean - 2*se, ymax = mean + 2*se), linewidth = 0.2) +
  labs(x = "Module (n = 1023)", y = "Mean pairwise correlation", title = "Pairwise correlations between 4 dataset projections") +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/4datasetproj_pairwiseCors_meanSE.svg"), width = 12, height = 5)


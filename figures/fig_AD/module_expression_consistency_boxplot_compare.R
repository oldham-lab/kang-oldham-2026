library(ggpubr)
library(DESeq2)
library(qs)
library(data.table)
library(AnnotationHub)
library(tidyverse)
library(cowplot)
options(bitmapType = 'cairo')

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_expression_consistency_fxns.R"))

# Load modules (bulk megaset, BC)
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Load lists of modules of interest
sealist <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/SEA_adVsCon_modList_filtered.qs"))
mitlist <- qread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/MIT_Multiome_MultiRegion/MIT_adVsCon_modlist_filtered.qs"))
titlevec <- c("SEA", "MIT")

# Load sn dataset means
copa_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
sea_means1 <- qread(paste0(copa_dir1, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
sea_means1[[2]] <- sea_means1[[2]][match(rownames(sea_means1[[1]]), rownames(sea_means1[[2]])), ]
sea_means1 <- lapply(sea_means1, as.data.frame)
copa_dir2 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/MIT_ADMultiome/")
sea_means2 <- qread(paste0(copa_dir2, "/sn_summary_tables/sn_summary_objects_log.qs"))$mean[[1]]
sea_means2[[2]] <- sea_means2[[2]][match(rownames(sea_means2[[1]]), rownames(sea_means2[[2]])), ]
sea_means2 <- lapply(sea_means2, as.data.frame)

# Initialize objects
highlight_mat <- data.frame("ct1" = colnames(sea_means1[[1]]),
                            "ct2" = c("Ast", NA, "End", "Exc L2-3 IT",
                                      NA, "Exc L5 ET", "Exc L5-6 IT", "Exc L5/6 NP",
                                      "Exc L6 CT", "Exc L5-6 IT", "Exc L5/6 IT Car3", "Exc L6b",
                                      "Inh LAMP5", NA, "Mic", "Oli",
                                      "OPC", "Inh PAX6", "Inh PVALB", NA,
                                      "Inh SST", NA, "Inh VIP", NA))

# Plotting

plotfunc <- function(indsmat,
                     file_name = "sea_vs_mit_seaPos.png",
                     save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels"),
                     hcol){
  inds <- indsmat$mod
  plots <- list()
  for(i in 1:length(inds)){
    if(hcol == 1){
      highlight_ct1 <- indsmat$subclass[i]
      highlight_ct2 <- highlight_mat[which(highlight_mat[, hcol] == highlight_ct1), 2]
    } else {
      highlight_ct2 <- indsmat$subclass[i]
      highlight_ct1 <- highlight_mat[which(highlight_mat[, hcol] == highlight_ct2), 1]
    }
    if(length(highlight_ct1) == 0 | length(highlight_ct2) == 0)
      next

    plots1 <- mapply(\(x, y){
      df <- x[rownames(x) %in% mod_bc[[inds[i]]], colnames(x) == highlight_ct1] 
      out <- data.frame("type" = y, "vals" = df, "gene" = 1:length(df))
      return(out)
    }, sea_means1[1:2], names(sea_means1)[1:2], SIMPLIFY = F) |>
      do.call(what = "rbind") |>
      ggpaired(x = "type", 
              y = "vals", 
              color = "type", 
              line.color = "gray", 
              line.size = 0.4,
              notch = T) +
        labs(x = "", y = paste0("Mean expression\nin ", highlight_ct1), title = paste0(titlevec[1], "_", indsmat[i, 1])) +
        theme(legend.position = "none", 
          axis.text.x = element_text(size = 4),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=5),
          axis.text.y = element_text(size=6),
          plot.title = element_text(hjust = 0.5))  

    plots2 <- mapply(\(x, y){
      df <- x[rownames(x) %in% mod_bc[[inds[i]]], colnames(x) == highlight_ct2] 
      out <- data.frame("type" = y, "vals" = df, "gene" = 1:length(df))
      return(out)
    }, sea_means2[1:2], names(sea_means2)[1:2], SIMPLIFY = F) |>
      do.call(what = "rbind") |>
      ggpaired(x = "type", 
              y = "vals", 
              color = "type", 
              line.color = "gray", 
              line.size = 0.4,
              notch = T) +
              labs(x = "", y = paste0("Mean expression\nin ", highlight_ct2), title = paste0(titlevec[2], "_", indsmat[i, 1])) +
              theme(legend.position = "none", 
                    axis.text.x = element_text(size = 4),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(size=5),
                    axis.text.y = element_text(size=6),
                    plot.title = element_text(hjust = 0.5))  
      plots[[i]] <- cowplot::plot_grid(plots1, plots2, nrow = 1, ncol = 2)
  }                    

  outall2 <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "h", axis = "bt")
  ggsave(outall2, file = file.path(save_path, file_name),  height = length(inds) * 1, width = 3, bg = "white", limitsize = F)
}

plotfunc(indsmat = sealist[[1]],
         file_name = "sea_vs_mit_seaPos.png",
         hcol = 1)
plotfunc(indsmat = sealist[[2]],
         file_name = "sea_vs_mit_seaNeg.png",
         hcol = 1)
plotfunc(indsmat = mitlist[[1]],
         file_name = "sea_vs_mit_mitPos.png",
         hcol = 2)
plotfunc(indsmat = mitlist[[2]],
         file_name = "sea_vs_mit_mitNeg.png",
         hcol = 2)




# Plot diff for all means (astrocyte)
diff_df <- data.frame("seadiff" = sea_means1[[1]]$Astrocyte - sea_means1[[2]]$Astrocyte,
                      "mitdiff" = sea_means2[[1]]$Ast - sea_means2[[2]]$Ast)
# cor(diff_df) 
# -0.18
library(ggplot2)
p <- ggplot(diff_df, aes(x = seadiff, y = mitdiff)) + 
  geom_point() +
  geom_abline()
ggsave(p, file = "~/test/test.png")
# Not much of a pattern between SEA and MIT
# What if we look at direction only and disregard magnitude?
diff_df2 <- data.frame("seadiff" = lapply(mod_bc, \(x){
    bla <- sea_means1[[1]]$Astrocyte[rownames(sea_means1[[1]]) %in% x] - sea_means1[[2]]$Astrocyte[rownames(sea_means1[[2]]) %in% x]
    return((sum(bla) > 0)/length(bla))  
  }) |> unlist(), 
  "mitdiff" = lapply(mod_bc, \(x){
    bla <- sea_means2[[1]]$Ast[rownames(sea_means2[[1]]) %in% x] - sea_means2[[2]]$Ast[rownames(sea_means2[[2]]) %in% x]
    return((sum(bla) > 0)/length(bla))
  }) |> unlist())
# cor(diff_df2, use = 'pairwise.complete.obs')
# -0.0368
p <- ggplot(diff_df2, aes(x = seadiff, y = mitdiff)) + 
  geom_point() +
  geom_abline()
ggsave(p, file = "~/test/test2.png")
sum(diff_df2[,1]==0 & diff_df2[,2] == 0, na.rm = T)
# 691
# Most modules (n=691) are lower in AD astrocytes in both SEA and MIT
# How do these overlap with mods of interest?
lapply(sealist, \(x) sum(x$mod[x$subclass == "Astrocyte"] %in% which(diff_df2[,1]==0 & diff_df2[,2] == 0)))
# 0 (out of 4 AD astrocyte), 0 (out of 0 con astrocyte)
lapply(mitlist, \(x) sum(x$mod[x$subclass == "Ast"] %in% which(diff_df2[,1]==0 & diff_df2[,2] == 0)))
# 0 (out of 0 AD ast), 0 (out of 89 con ast)
# i.e. there are no astrocyte mods that share directionality in both SEA and MIT


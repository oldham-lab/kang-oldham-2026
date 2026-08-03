# v2.1
# - Adding all comparisons in a loop (in addition to All AD vs. Con DFC)

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

# 0. Create various vectors and lists for all comparisons
save_suffix_vec <- c("AllADVsCon_DFC",
                     "EarlyADVsCon_DFC",
                     "LateVsEarlyAD_DFC",
                     "AllADVsCon_MTG",
                     "EarlyADVsCon_MTG",
                     "LateVsEarlyAD_MTG",
                     "AllADVsCon_DFC_ROSMAP",
                     "EarlyADVsCon_DFC_ROSMAP",
                     "LateVsEarlyAD_DFC_ROSMAP",
                     "AllADVsCon_MTG_ROSMAP",
                     "EarlyADVsCon_MTG_ROSMAP",
                     "LateVsEarlyAD_MTG_ROSMAP",
                     "APOE_DFC"#,
                    # "CMC_vs_SZBD"
                     )

dat_vec <- c("gabitto_vs_liu_DFC_AllADVsCon",    
             "gabitto_vs_liu_DFC_earlyVsCon",        
             "gabitto_vs_liu_DFC_lateVsEarly",        
             "gabitto_vs_liu_MTG_AllADVsCon",         
             "gabitto_vs_liu_MTG_earlyVsCon",   
             "gabitto_vs_liu_MTG_lateVsEarly",        
             "gabitto_vs_liu_DFC_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_DFC_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_DFC_lateVsEarly_ROSMAP",
             "gabitto_vs_liu_MTG_AllADVsCon_ROSMAP",    
             "gabitto_vs_liu_MTG_earlyVsCon_ROSMAP", 
             "gabitto_vs_liu_MTG_lateVsEarly_ROSMAP",
             "gabitto_vs_liu_APOE_DFC"#,
             #"brainSCOPE_CMC_vs_SZBD"
             ) 

title_vec <- c("Con vs all AD (DFC)",
             "Early AD vs Con (DFC)",
             "Late vs Early AD (DFC)",
             "Con vs all AD (MTG)",
             "Early AD vs Con (MTG)",
             "Late vs Early AD (MTG)",
             "Con vs all AD (ROSMAP, DFC)",
             "Early AD vs Con (ROSMAP, DFC)",
             "Late vs Early AD (ROSMAP, DFC)",
             "Con vs all AD (ROSMAP, MTG)",
             "Early AD vs Con (ROSMAP, MTG)",
             "Late vs Early AD (ROSMAP, MTG)",
             "APOE 4/4 vs 3/3 (DFC)"#,
             #"SCZ vs con (DFC)"
             )

copa_dir_strings <- list(
  c("finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC",
  "finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_PFC"),
  c("finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsCon_DFC",
  "finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVsCon_PFC"),
  c("finalNonNorm_minsize10_unmerged/SEAAD2024_lateVsEarly_DFC",
  "finalNonNorm_minsize10_unmerged/MIT_ADMultiome_lateVsEarly_PFC"),
  c("finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_MTG",
  "finalNonNorm_minsize10_unmerged/MIT_ADMultiome_AllADVsCon_MTC"),
  c("finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsCon_MTG",
  "finalNonNorm_minsize10_unmerged/MIT_ADMultiome_earlyVsCon_MTC"),
  c("finalNonNorm_minsize10_unmerged/SEAAD2024_lateVsEarly_MTG",
  "finalNonNorm_minsize10_unmerged/MIT_ADMultiome_lateVsEarly_MTC"),
  c("rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC",
  "rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_DFC"),
  c("rosmap_AD/rosmap_AD_SEAAD2024_earlyVsCon_DFC",
  "rosmap_AD/rosmap_AD_MITMultiome_earlyVsCon_DFC"),
  c("rosmap_AD/rosmap_AD_SEAAD2024_lateVsEarly_DFC",
  "rosmap_AD/rosmap_AD_MITMultiome_lateVsEarly_DFC"),
  c("rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_MTG",
  "rosmap_AD/rosmap_AD_MITMultiome_AllADVsCon_MTG"),
  c("rosmap_AD/rosmap_AD_SEAAD2024_earlyVsCon_MTG",
  "rosmap_AD/rosmap_AD_MITMultiome_earlyVsCon_MTG"),
  c("rosmap_AD/rosmap_AD_SEAAD2024_lateVsEarly_MTG",
  "rosmap_AD/rosmap_AD_MITMultiome_lateVsEarly_MTG"),
  c("finalNonNorm_minsize10_unmerged/SEAAD2024_DFC_apoe_44_vs_33",
  "finalNonNorm_minsize10_unmerged/MIT_ADMultiome_APOE")#,
  # c("brainSCOPE_CMC",
  # "brainSCOPE_SZBDMultiseq")
)

dx_list <- list(
  c("AD", "Con"),
  c("Early", "Con"),
  c("Late", "Early"),
  c("AD", "Con"),
  c("Early", "Con"),
  c("Late", "Early"),
  c("AD", "Con"),
  c("Early", "Con"),
  c("Late", "Early"),
  c("AD", "Con"),
  c("Early", "Con"),
  c("Late", "Early"),
  c("44", "33")
)

index_suffixes <- list(
  c("indices_over_all_datasets_Subclass_Alzheimers_topmodposbc_mean.csv",
    "indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv",
    "indices_se_SubclassADtopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Alzheimers.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassAD.csv",
    "indices_se_RNA.SubclassCon.csv"),
  c("indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Con.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassEarly.csv",
    "indices_se_RNA.SubclassCon.csv"
   ),
  c("indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_SubclassLatetopmodposbc.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_RNA.SubclassLate.csv",
    "indices_se_RNA.SubclassEarly.csv"
    ),
  c("indices_over_all_datasets_Subclass_Alzheimers.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_SubclassADtopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Alzheimers.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassAD.csv",
    "indices_se_RNA.SubclassCon.csv"
    ),
  c("indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Con.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassEarly.csv",
    "indices_se_RNA.SubclassCon.csv"    
    ),
  c("indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_SubclassLatetopmodposbc.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",    
    "indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_RNA.SubclassLate.csv",
    "indices_se_RNA.SubclassEarly.csv"        
    ),
  c("indices_over_all_datasets_Subclass_Alzheimers.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_SubclassADtopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Alzheimers.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassAD.csv",
    "indices_se_RNA.SubclassCon.csv"    
    ),
  c("indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Con.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassEarly.csv",
    "indices_se_RNA.SubclassCon.csv"    
    ),
  c("indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_SubclassLatetopmodposbc.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",  
    "indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_RNA.SubclassLate.csv",
    "indices_se_RNA.SubclassEarly.csv"        
    ),
  c("indices_over_all_datasets_Subclass_Alzheimers.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_SubclassADtopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Alzheimers.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassAD.csv",
    "indices_se_RNA.SubclassCon.csv"        
    ),
  c("indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Con.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",
    "indices_se_SubclassContopmodposbc.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_over_all_datasets_Subclass_Control.csv",
    "indices_se_RNA.SubclassEarly.csv",
    "indices_se_RNA.SubclassCon.csv"    
    ),
  c("indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_SubclassLatetopmodposbc.csv",
    "indices_se_SubclassEarlytopmodposbc.csv",      
    "indices_over_all_datasets_Subclass_Late.csv",
    "indices_over_all_datasets_Subclass_Early.csv",
    "indices_se_RNA.SubclassLate.csv",
    "indices_se_RNA.SubclassEarly.csv"       
    ),
  c("indices_over_all_datasets_Subclass_44.csv",
    "indices_over_all_datasets_Subclass_33.csv",
    "indices_se_Subclass44topmodposbc.csv",
    "indices_se_Subclass33topmodposbc.csv",
    "indices_over_all_datasets_Subclass_44.csv",
    "indices_over_all_datasets_Subclass_33.csv",
    "indices_se_RNA.Subclass44.csv",
    "indices_se_RNA.Subclass33.csv"
   )#,
  # c("indices_over_all_datasets_Subclass_Schizophrenia.csv",
  #   "indices_over_all_datasets_Subclass_control.csv",
    
    
  #  "indices_over_all_datasets_Subclass_Schizophrenia.csv",
  #  "indices_over_all_datasets_Subclass_control.csv",

   
  #  )
)

dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/dcopa_shared")
# 1. Gabitto + Liu shared dCoPA results (for highlighting grouped barplots)
dcopa_shared <- mapply(\(x, y){
  fread(file.path(dir1, x, "shared_output_table.csv"), data.table = F)
}, dat_vec, title_vec, SIMPLIFY = F)

proj_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output")
# 2. Gabitto mean expr by celltype tables (for boxplots)
mean_sum_list <- lapply(copa_dir_strings, \(x){
  lapply(x, \(y){
    qread(file.path(proj_dir, y, "sn_summary_tables/sn_summary_objects_log.qs"))$mean$Subclass
  })
}) |>
  setNames(save_suffix_vec)

# 3. Load projection indices (mean, SE)
index_list <- mapply(\(x, y){
  lapply(seq_along(x), \(i){
    if(i == 1){
      ind <- 1:4
    } else {
      ind <- 5:8
    }
    z <- x[i]
    out_list <- list()
    out_list[[1]] <- fread(data.table = F, file = file.path(proj_dir, z, "/sn_proj_indices/", index_dirs[d], y[ind[1]])) 
    if("module" %in% colnames(out_list[[1]])){
      out_list[[1]] <- out_list[[1]] |> dplyr::select(!module)
    }
    out_list[[2]] <- fread(data.table = F, file = file.path(proj_dir, z, "/sn_proj_indices/", index_dirs[d], y[ind[2]])) |>
      select(colnames(out_list[[1]]))
    out_list[[3]] <- fread(data.table = F, file = file.path(proj_dir, z, "/sn_proj_indices/", index_dirs[d], y[ind[3]])) 
    if("module" %in% colnames(out_list[[3]])){
      out_list[[3]] <- out_list[[3]] |> dplyr::select(!module)
    }
    out_list[[4]] <- fread(data.table = F, file = file.path(proj_dir, z, "/sn_proj_indices/", index_dirs[d], y[ind[4]])) |>
      select(colnames(out_list[[3]]))
    return(out_list)
  })
}, copa_dir_strings, index_suffixes, SIMPLIFY = F)


##### Run loop for all comparisons
for(i in 1:length(dcopa_shared)){

  # Gabitto AD indices
  gabad <- index_list[[i]][[1]][[1]]
  gabad_se <- index_list[[i]][[1]][[3]]

  # Gabitto control indices
  gabcon <- index_list[[i]][[1]][[2]]
  gabcon_se <- index_list[[i]][[1]][[4]]

   # Liu AD indices
  mitad <- index_list[[i]][[2]][[1]]
  mitad_se <- index_list[[i]][[2]][[3]]

   # Liu control indices
  mitcon <- index_list[[i]][[2]][[2]]
  mitcon_se <- index_list[[i]][[2]][[4]]

  allcts1 <- intersect(colnames(gabcon), colnames(gabad))
  allcts2 <- intersect(colnames(mitcon), colnames(mitad))

  # Align columns
  mitad     <- mitad     |> select(allcts2)
  mitad_se  <- mitad_se  |> select(allcts2)
  mitcon    <- mitcon    |> select(allcts2)
  mitcon_se <- mitcon_se |> select(allcts2)

  # Plot projection indices
  projplotlist1 <- lapply(unique(dcopa_shared[[i]]$mod), \(j){ # Gabitto

    # Fetch projection index data for the module
    plotdf1 <- data.frame("ct" = rep(allcts1, 2),
                        "dataset" = c(rep(paste0("Gabitto ", title_vec[i]), 2 * length(allcts1))),
                        "ind" = c(unlist(gabcon[j, ]), 
                                  unlist(gabad[j, ])),
                        "ind_se" = c(unlist(gabcon_se[j, ]), 
                                      unlist(gabad_se[j, ])),
                        "dx" = c(rep(dx_list[[i]][2], length(allcts1)),
                                  rep(dx_list[[i]][1], length(allcts1)))) |>
      dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                    ct = factor(ct, levels = allcts1))

    # Fetch dCoPA significance data for given module
    highlight_mat1 <- dcopa_shared[[i]] |> 
      filter(mod == j, sig_FDR_Gabitto_2024 == T) |>
      select(Celltype_Gabitto_2024, P_value_Gabitto_2024) |>
      filter(!duplicated(Celltype_Gabitto_2024))
    
    # Plot grouped barplot of projection indices with asterisks for dCoPA significance
    plotdf1 <- plotdf1 |>
      left_join(highlight_mat1, by = join_by(ct == Celltype_Gabitto_2024)) |>
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

    plotdf2 <- data.frame("ct" = rep(allcts2, 2),
                          "dataset" = c(rep(paste0("Liu ", title_vec[i]), 2 * length(allcts2))),
                          "ind" = c(unlist(mitcon[j, ]),
                                    unlist(mitad[j, ])),
                          "ind_se" = c(unlist(mitcon_se[j, ]),
                                      unlist(mitad_se[j, ])),
                          "dx" = c(rep(dx_list[[i]][2], length(allcts2)),
                                  rep(dx_list[[i]][1], length(allcts2)))) |>
      dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                    ct = factor(ct, levels = allcts2))
    
    # Fetch dCoPA significance data for given module
    highlight_mat2 <- dcopa_shared[[i]] |> 
      filter(mod == j, sig_FDR_Liu_2025 == T) |>
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

    dodge_width <- 0.6
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
  boxplotlist_gab <- lapply(unique(dcopa_shared[[i]]$mod), \(j){
    highlight_mat <- dcopa_shared[[i]] |> 
      filter(mod == j, sig_FDR_Liu_2025 == T) |>
      select(Celltype_Gabitto_2024, P_value_Gabitto_2024)
    summary_mat <- mean_sum_list[[i]][[1]]

    p <- mapply(\(x, y){
      df <- x[rownames(x) %in% mods[[which(these_mods == j)]], colnames(x) %in% highlight_mat[,1], drop = F] 
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

  boxplotlist_liu <- lapply(unique(dcopa_shared[[i]]$mod), \(j){
    highlight_mat <- dcopa_shared[[i]] |> 
      filter(mod == j, sig_FDR_Liu_2025 == T) |>
      select(Celltype_Liu_2025, P_value_Liu_2025)
    summary_mat <- mean_sum_list[[i]][[2]]

    p <- mapply(\(x, y){
      df <- x[rownames(x) %in% mods[[which(these_mods == j)]], colnames(x) %in% highlight_mat[,1], drop = F] 
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

  out_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/dCoPA_snapshots/"), save_suffix_vec[i])
  if(!dir.exists(out_dir))
    dir.create(out_dir, recursive = T)

  for(j in seq_along(unique(dcopa_shared[[i]]$mod))){
    boxplot_combined <- plot_grid(boxplotlist_gab[[j]], boxplotlist_liu[[j]], nrow = 1)
    pall <- suppressMessages(plot_grid(all_plots[[1]][[which(these_mods == unique(dcopa_shared[[i]]$mod)[j])]],
                                      all_plots[[2]][[which(these_mods == unique(dcopa_shared[[i]]$mod)[j])]],
                                      projplotlist1[[j]],
                                      boxplot_combined,
                                      nrow = 4, 
                                      rel_heights = c(0.5, 0.3, 1, 0.5)))
    ggsave(pall, 
          file = file.path(out_dir, paste0(unique(dcopa_shared[[i]]$mod)[j], ".svg")), 
          device = svglite::svglite, 
          width = 6, 
          height = 12, 
          bg = "white")
    cat(j, " ")
  }
  cat(i, "\n")
}

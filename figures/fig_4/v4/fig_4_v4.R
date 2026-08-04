# v4:
# adding Liu projections back, changed first column module, changed bottom

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

version <- "v4"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/"), version)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

# Load bulk-related objects
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table=F, file = file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
mod_seed <- qread(file.path(module_output_dir, "modules/unmerged_modules.qs"))
modulelengths <- unlist(lapply(mods,length))
modulelengths_seed <- unlist(lapply(mod_seed,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])


# # mod_ids are out of 1158
# mod_id <- c(# 347, # oligo
#             # 635, # VIP although VIP is not a top 10 seed gene
#             # 822, # Endo
#             681, # Microglia
#             669, # left top, postsynaptic
#             # 253, # mystery module
#             1007, # jun/fos, top right
#             #173 # left bottom, electron transport chain
#             59
#             )

# Other candidates for first row: (1023 index)
# 76 (spliceosome)


mod_id <- c(82, 669, 681, 1007)  #1158 indices

# Figure out the final mod labels 
# - with filtering by topmodposbc AND >=2 significant bulk plat cors
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods_final <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
mod_labels_final <- which(these_mods_final %in% mod_id)
these_mods_final_df <- data.frame("old" = these_mods_final, "new" = 1:length(these_mods_final))

# candidates for first column (ids subject to change):
# - mod 17: dnm1l, tomm70, vps35, napg
# mod 21: atp5f1a, pgam1, pgk1, pfkm (housekeeping genes)
# mod 27: ubqln1, cul3, cul1, fbxo11 (cellcycle, apoptosis)
# mod 43: psmd14, psmc2, cops5
# 59: mitochondrial micu1, ube2q1, ndufa9

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
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/")
all_plots_full = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_seed.qs"))[mod_index],
                 qread(file.path(save_dir1,"sn_proj_objects","bulkcor_object.qs"))[mod_index],
                 qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs"))[mod_index])
plots_bc_full <- list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs"))[mod_index])

all_plots_full[[3]][[2]]$data$SetName_breaks <- recode(all_plots_full[[3]][[2]]$data$SetName_breaks, "GOMF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY_INVOLVED_IN_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL" = "GOMF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY_INVOLVED_\nIN_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL")

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
          plot.title = element_blank(),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "right") +
    scale_color_manual(values = cols2) +
    labs(subtitle = subt,
         colour = legt
         ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))
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
    dplyr::mutate(SetName_breaks = factor(SetName_breaks, levels=rev(unique(SetName_breaks)))) |>
    ggplot(aes(x = SetName_breaks, y = pval)) +
        newtheme + 
        geom_bar(stat="identity") 
  x2$layers[[2]] <- x$layers[[2]]
  # Calculate string lengths of SetName_breaks
  breaks_lengths <- regexpr("\n", x2$data$SetName_breaks) - 1
  # Scale y-label size by length of label
  maxnamesize <- min(max(nchar(as.character(x2$data$SetName))), max(breaks_lengths)) # take the minimum of the maximum of SetName and SetName_breaks
  #ysizevar <- (max(nchar(as.character(x2$data$SetName))) - 31)/58 * -4 + 10
  ysizevar <- (maxnamesize - 31)/58 * -4 + 10
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
d=2

lein <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Cell_Type_1.csv")) |>
  dplyr::select(!module)
lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_se_topmodposbc.csv")) 
lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]

mitcon <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/mod_means/"), index_dirs[d], "mod_means_Con_bulk_megaset.csv"))
mitse <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/se/"), index_dirs[d], "se_Con_bulk_megaset.csv"))

SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
  dplyr::select(!module)
SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) 
SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]

mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
  dplyr::select(!module)
mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) 
mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]

# Change celltypes to match

colnames(lein)[c(1, 3, 15, 16)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
colnames(mora)[c(1, 3, 14, 15)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")

colnames(lein_se) <- colnames(lein)
colnames(SEAcon_se) <- colnames(SEAcon)
colnames(mora_se) <- colnames(mora)

allcts <- c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
        "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")

allcts_cap <- c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
        "Chandelier", "LAMP5", "LAMP5 LHX6", "PAX6", "PVALB", "SNCG", "SST", "SST CHODL", "VIP",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")


# Fill out celltypes
lein[setdiff(allcts, colnames(lein))] <- NA
lein <- lein[, match(allcts, colnames(lein))]
lein_se[setdiff(allcts, colnames(lein_se))] <- NA
lein_se <- lein_se[, match(allcts, colnames(lein_se))]

mitcon[setdiff(allcts, colnames(mitcon))] <- NA
mitcon <- mitcon[, match(allcts, colnames(mitcon))]
mitse[setdiff(allcts, colnames(mitse))] <- NA
mitse <- mitse[, match(allcts, colnames(mitse))]

SEAcon[setdiff(allcts, colnames(SEAcon))] <- NA
SEAcon <- SEAcon[, match(allcts, colnames(SEAcon))]
SEAcon_se[setdiff(allcts, colnames(SEAcon_se))] <- NA
SEAcon_se <- SEAcon_se[, match(allcts, colnames(SEAcon_se))]

mora[setdiff(allcts, colnames(mora))] <- NA
mora <- mora[, match(allcts, colnames(mora))]
mora_se[setdiff(allcts, colnames(mora_se))] <- NA
mora_se <- mora_se[, match(allcts, colnames(mora_se))]

# Switch to capitalized labels
colnames(lein) <- allcts_cap
colnames(mitcon) <- allcts_cap
colnames(SEAcon) <- allcts_cap
colnames(mora) <- allcts_cap
colnames(lein_se) <- allcts_cap
colnames(mitse) <- allcts_cap
colnames(SEAcon_se) <- allcts_cap
colnames(mora_se) <- allcts_cap

# Plot projection indices
projplotlist <- lapply(mod_id, \(i){
  n <- length(allcts)
  plotdf <- data.frame("ct" = rep(allcts_cap, 4),
                       "dataset" = c(rep("Jorstad et al. 2023 (PFC)", n), 
                                     rep("Gabitto et al. 2024 (control PFC)", n), 
                                     rep("Liu et al. 2025 (normal human dorsolateral prefrontal cortex)", n), 
                                     rep("Morabito et al. 2021 (control PFC)", n)),
                       "ind" = c(unlist(lein[i, ]), 
                                 unlist(SEAcon[i, ]), 
                                 unlist(mitcon[i, ]), 
                                 unlist(mora[i, ])),
                       "ind_se" = c(unlist(lein_se[i, ]), 
                                    unlist(SEAcon_se[i, ]), 
                                    unlist(mitse[i, ]), 
                                    unlist(mora_se[i, ]))) |>
    dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                  ct = factor(ct, levels = allcts_cap))

  # Calculate and print pairwise cor mat
  cat(i, "pairwise cors: \n")
  tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), 
                       "Gabitto 2024" = unlist(SEAcon[i, ]), 
                       "Liu 2025" = unlist(mitcon[i, ]),
                       "Morabito 2021" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
  colnames(tempcor) <- c("Jorstad 2023","Gabitto 2024", "Liu 2025", "Morabito 2021")
  rownames(tempcor) <- c("Jorstad 2023", "Gabitto 2024", "Liu 2025", "Morabito 2021")
  print(tempcor)
 
  return(ggplot(plotdf, aes(x = ct, y = ind, fill = ct)) +
    theme_classic() +
    geom_col(position = position_dodge(), alpha = 0.5) +
    geom_errorbar(aes(ymin = ind - 2 * ind_se,
                      ymax = ind + 2 * ind_se),
                      width = 0.2,
                      linewidth = 0.3,
                      position = position_dodge(0.5)) +
    scale_fill_manual(values = scales::hue_pal()(n)) +
    theme(text = element_text(family = "sans", color = "black", size = 12),
          legend.position="none", 
          axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(color = "black"),
          strip.background = element_rect(fill = "white")) +
    facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_y") +
    labs(y = index_xaxis_names[d], x = "") +
    scale_y_continuous(breaks = c(0, 0.5, 1)))
})

## Plot individual panels
for(j in seq_along(all_plots[[1]])){
  # Expression plot:
  ggsave(plots_bc[[1]][[j]], 
        file = file.path(save_dir, paste0("panel_", j, "_1.pdf")),
        width = 6, height = 3)
  ggsave(plots_bc[[1]][[j]], 
        file = file.path(save_dir, paste0("panel_", j, "_1.svg")),
        width = 6, height = 3)

  # Bulk cor plot:
  ggsave(all_plots[[2]][[j]],
          file = file.path(save_dir, paste0("panel_", j, "_2.pdf")),
          width = 6, height = 4)
    ggsave(all_plots[[2]][[j]],
          file = file.path(save_dir, paste0("panel_", j, "_2.svg")),
          width = 6, height = 4)

  # GSEA plot:
  #max_setname <- lapply(all_plots[[3]][[1]]$data$SetName |> as.character(), nchar) |> unlist() |> max()
  max_setname <- 35
  ggsave(all_plots[[3]][[j]], file = file.path(save_dir, paste0("panel_", j, "_3.pdf")),
          height = 1.5, width = (3 + max_setname/30 * 2.5))
  ggsave(all_plots[[3]][[j]], file = file.path(save_dir, paste0("panel_", j, "_3.svg")),
          height = 1.5, width = (3 + max_setname/30 * 2.5))

  # Projection indices:
  ggsave(projplotlist[[j]], file = file.path(save_dir, paste0("panel_", j, "_4.pdf")), height = 5, width = 5.5)
  ggsave(projplotlist[[j]], file = file.path(save_dir, paste0("panel_", j, "_4.svg")), height = 5, width = 5.5)


}

########## 
# Panel Q:
# Stacked plots (x axis is ordered modules):
# 1. Module size (seed + topmodposbc)
# 2. Median seed gene cor
# 3. PC1 var explained
# 4. Mean REI correlation
##########

input_df <- data.frame("mod index" = 1:length(these_mods_final),
                       "Seed genes" = modulelengths_seed[these_mods_final], 
                       "Topmodposbc genes" = modulelengths[these_mods_final], 
                       "VE" = vevec[these_mods_final])

# Gather mean REI
cordf <- lapply(these_mods_final, \(i){
  tempdf <- data.frame("lein" = unlist(lein[i, ]), 
                       "SEA" = unlist(SEAcon[i, ]), 
                       "mora" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = "pairwise.complete.obs")
  return(data.frame("mod" = i, "cors" = tempcor[upper.tri(tempcor)]))
})
input_df$REI_mean_cor <- cordf |> do.call(what = "rbind") |> group_by(mod) |> summarise(mean = mean(cors)) |> pull(mean)

# Gather median pairwise cor of seed genes

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
expr_t <- t(expr[,3:ncol(expr)])
colnames(expr_t) <- expr[,2]
simMat <- cor(expr_t)

pair_cor <- lapply(mod_seed[these_mods_final], \(x){
  temp <- simMat[rownames(simMat) %in% x, colnames(simMat) %in% x]
  return(median(temp[upper.tri(temp)]))
})

input_df$pair_cor <- unlist(pair_cor)
fwrite(input_df, file = file.path(save_dir, "panel_Q_statistics.csv"))

# faceted ggplot
input_df_long <- pivot_longer(input_df, cols = -mod.index, names_to = "variable",
   values_to = "value") |>                              
    mutate(facet = case_when(
      variable %in% c("length_seed", "length_bc") ~ "lengths",             
      TRUE ~ variable                                              
    ))     


plotlist <- 
  p <- ggplot(input_df_long, aes(x = mod.index, y = value, color = variable)) +  
      theme_classic() +                    
      geom_line() +
      facet_wrap(~ facet, ncol = 1, scale = "free_y") +
      scale_color_manual(values = c("length_seed" = "steelblue", "length_bc" =   
      "tomato", "VE" = "black", "REI_mean_cor" = "black", "pair_cor" = "black"),
                        breaks = c("length_seed", "length_bc")  ) +
      theme( 
        legend.title = element_blank(),         
        legend.position      = c(0.05, 0.95),  # top-right of plot area
        legend.justification = c(0, 1),         # anchor point is top-right of legend box                                          
        legend.background    = element_rect(fill = alpha("white", 0.7))
      ) +
      labs(x = "Module indices (n = 1016)")                                                     
                  
ggsave(p, file = file.path(save_dir, "panel_Q.pdf"), height = 9)
ggsave(p, file = file.path(save_dir, "panel_Q.svg"), height = 9)


# patchwork
library(patchwork)                                             

# Length
plist <- list()

plist[[1]] <- input_df[,c(1:3)] |>
  pivot_longer(!mod.index, names_to = "lengthtype", values_to = "value") |>
  ggplot(aes(x = mod.index, y = value, color = lengthtype)) +  
    theme_classic() +                    
    geom_line(linewidth = 0.1) +
    scale_color_manual(values = c("Seed.genes" = "steelblue", "Topmodposbc.genes" = "tomato"),
                       labels = c("Seed.genes" = "Seed genes", "Topmodposbc.genes" = "Topodposbc genes")) +
    theme( 
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.title = element_blank(),         
      legend.position      = c(0.05, 0.95),  # top-right of plot area
      legend.justification = c(0, 1),         # anchor point is top-right of legend box                                          
      legend.background    = element_rect(fill = alpha("white", 0.7))
    ) +
    labs(x = "Module indices (n = 1016)", 
         y = "# of module genes",
         title = "Size of modules") +
    guides(color = guide_legend(override.aes = list(linewidth = 2))) 

ggsave(plist[[1]], file = file.path(Sys.getenv("SCRATCH_DIR", "~/test"), "test.svg"))

# %VE
plist[[2]] <- input_df[,c(1, 4)] |>
  ggplot(aes(x = mod.index, y = VE)) +  
    theme_classic() +                    
    geom_line(linewidth = 0.1) +
    theme(      
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    labs(x = "Module indices (n = 1016)", 
         y = "% Variance\nexplained",
         title = "% Variance explained by PC1 of topmodposbc genes")   

# REI cor
plist[[3]] <- input_df[,c(1, 5)] |>
  ggplot(aes(x = mod.index, y = REI_mean_cor)) +  
    theme_classic() +                    
    geom_line(linewidth = 0.1) +
    theme( 
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    labs(x = "Module indices (n = 1016)", 
         y = "Correlation",
         title = "Mean pairwise correlation of REI indices")   

# 
plist[[4]] <- input_df[,c(1, 6)] |>
  ggplot(aes(x = mod.index, y = pair_cor)) +  
    theme_classic() +                    
    geom_line(linewidth = 0.1) +
    theme(      
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    labs(x = "Module indices (n = 1016)", 
         y = "Correlation",
         title = "Median pairwise correlation of seed genes")  

pall <- wrap_plots(plist, ncol = 1) +                      
  plot_layout(axes = "collect_x") &                     
  labs(x = "Module indices (n = 1016)") 
ggsave(pall, file = file.path(save_dir, "panel_Q.pdf"), height = 7)
ggsave(pall, file = file.path(save_dir, "panel_Q.svg"), height = 7)

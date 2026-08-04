# v5:
# changed first column, changed bottom

library(qs)
library(data.table)
library(scales)
library(ggplot2)
library(cowplot)
library(eulerr)
library(tidyverse)
library(ggdendro)
library(ggplotify)
library(ggh4x)
library(showtext)
showtext_auto()
options(bitmapType = 'cairo')
# Create snapshots for 4 modules

version <- "v6"
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


mod_id <- c(112, 669, 681, 1007)  #1158 indices

# # 112
#               Jorstad 2023 Gabitto 2024  Liu 2025 Morabito 2021                                                                                                            
# Jorstad 2023     1.0000000    0.8597202 0.7662443     0.7934046                                                                                                            
# Gabitto 2024     0.8597202    1.0000000 0.8335956     0.8961764
# Liu 2025         0.7662443    0.8335956 1.0000000     0.5909439                                                                                                            
# Morabito 2021    0.7934046    0.8961764 0.5909439     1.0000000 

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
    labs(x = "", y = bquote(-log[10]~"(p-val)"), title = "Geneset enrichment") + 
    force_panelsizes(rows = unit(1, "in"), cols = unit(1, "in"))  
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

lein <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |>
  dplyr::select(!module)
lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_se_topmodposbc.csv"))
lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]

mitcon <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/mod_means/"), index_dirs[d], "mod_means_Con_bulk_megaset.csv"))
mitse <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/se/"), index_dirs[d], "se_Con_bulk_megaset.csv"))

# Gabitto (SEA-AD) now routes through the HGNC-harmonized Python SEA-AD output
# (SEAAD2024_full_python_output), matching the means/se format of the MIT/Liu path
# above (no `module` column). Supersedes the pre-fix R sn_proj_indices output.
SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/mod_means/", index_dirs[d], "/mod_means_Con_bulk_megaset.csv"))
SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/se/", index_dirs[d], "/se_Con_bulk_megaset.csv"))
SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]

mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv")) |>
  dplyr::select(!module)
mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) 
mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]

# Change celltypes to match

colnames(lein)[c(1, 3, 15, 16)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
colnames(mora)[c(1, 3, 14, 15)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")

colnames(lein_se) <- colnames(lein)
colnames(SEAcon_se) <- colnames(SEAcon)
colnames(mora_se) <- colnames(mora)

allcts <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
        "Chandelier", "Pvalb", "Sst", "Sst Chodl", "Lamp5 Lhx6", "Lamp5", "Pax6", "Sncg", "Vip",
        "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")

allcts_cap <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
        "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5", "PAX6", "SNCG", "VIP",
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

# Projection barplot palette (matched to fig_3_v4.R; mapped by cell type, in fig_3's order)
proj_pal_order <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
                    "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5", "PAX6", "SNCG", "VIP",
                    "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
proj_pal <- setNames(
  colorRampPalette(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                     "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"))(length(proj_pal_order)),
  proj_pal_order
)

# Plot projection indices (shared builder; also used by the standalone panel_mp.R).
# make_proj_panel() returns the barplot + its 4x4 pairwise-cor matrix; add_corr_brackets()
# overlays the nested dataset-correlation brackets (see panel_mp_plot.R).
source(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/v6"), "panel_mp_plot.R"))
proj_dfs <- list(lein, SEAcon, mitcon, mora)
proj_ses <- list(lein_se, SEAcon_se, mitse, mora_se)
proj_panels  <- lapply(mod_id, \(i) make_proj_panel(i, proj_dfs, proj_ses, allcts_cap, proj_pal, index_xaxis_names[d]))
projplotlist <- lapply(proj_panels, `[[`, "plot")
proj_cormats <- lapply(proj_panels, `[[`, "cormat")
for (j in seq_along(mod_id)) { cat(mod_id[j], "pairwise cors:\n"); print(proj_cormats[[j]]) }

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
  max_setname <- 45
  ggsave(all_plots[[3]][[j]], file = file.path(save_dir, paste0("panel_", j, "_3.pdf")),
#          height = 1.5, width = (3 + max_setname/30 * 2.5))
          height = 2, width = 8)

  ggsave(all_plots[[3]][[j]], file = file.path(save_dir, paste0("panel_", j, "_3.svg")),
#          height = 1.5, width = (3 + max_setname/30 * 2.5))
          height = 2, width = 8)


  # Projection indices:
  ggsave(projplotlist[[j]], file = file.path(save_dir, paste0("panel_", j, "_4.pdf")), height = 5, width = 5.5)
  ggsave(projplotlist[[j]], file = file.path(save_dir, paste0("panel_", j, "_4.svg")), height = 5, width = 5.5)
  # ... and with nested dataset-correlation brackets (used in the assembled figure):
  ggsave(add_corr_brackets(projplotlist[[j]], proj_cormats[[j]]),
         file = file.path(save_dir, paste0("panel_", j, "_4_with_brackets.svg")), height = 5, width = 6.6)
}

########## 
# Panel Q:
# Stacked plots (x axis is ordered modules):
# 1. Module size (seed + topmodposbc)
# 2. Median seed gene cor
# 3. PC1 var explained
# 4. Mean REI correlation
##########

input_df <- data.frame("mod_index" = 1:length(these_mods_final),
                       "Seed genes" = modulelengths_seed[these_mods_final], 
                       "Topmodposbc genes" = modulelengths[these_mods_final])

input_df_long <- input_df |> pivot_longer(!mod_index, names_to = "mod_type", values_to = "length")

p <- ggplot(input_df_long, aes(x = length, color = mod_type)) +
  theme_classic() + 
  scale_color_manual(values = c("#E63946", "#457B9D"), 
                     labels = c("Seed genes", "Expanded module genes\n(topmodposbc)")) +
  geom_density(key_glyph = "path") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = c(0.95, 0.95),                                           
        legend.justification = c(1, 1),                                            
        legend.key = element_blank(),
        legend.title = element_blank()) +                                            
  guides(color = guide_legend(override.aes = list(fill = "NA",linetype = 1, linewidth = 1))) +                                                     
  labs(x = "Module size",
       y = "Density" 
       #title = "Distribution of module lengths"
       )

ggsave(p, file = file.path(save_dir, "panel_Q.svg"), height = 2, width = 3.5)

######
# Panel R
##########
 
# Gather median pairwise cor of seed genes

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
expr_t <- t(expr[,3:ncol(expr)])
colnames(expr_t) <- expr[,2]
simMat <- cor(expr_t)

# ## line plot of median cors
# pair_cor <- lapply(mod_seed[these_mods_final], \(x){
#   temp <- simMat[rownames(simMat) %in% x, colnames(simMat) %in% x]
#   return(median(temp[upper.tri(temp)]))
# })

# plotdf <- data.frame("mod.index" = 1:length(pair_cor),
#                      "pair_cor" = unlist(pair_cor))

# p <- plotdf |> ggplot(aes(x = mod.index, y = pair_cor)) +  
#     theme_classic() +                    
#     geom_line(linewidth = 0.1) +
#     theme(      
#       text = element_text(size = 12),
#       plot.title = element_text(hjust = 0.5, size = 12)
#     ) +
#     labs(x = "Module indices (n = 1016)", 
#          y = "Correlation",
#          title = "Median pairwise correlation of seed genes")  
# ggsave(p, file = file.path(save_dir, "panel_R.svg"), height = 2, width = 7)

## mean cors plus error bars
pair_cor <- lapply(mod_seed[these_mods_final], \(x){
  temp <- simMat[rownames(simMat) %in% x, colnames(simMat) %in% x]
  temp <- temp[upper.tri(temp)]
  outdf <- data.frame("mean" = mean(temp),
                      "se" = sd(temp)/sqrt(length(temp)))
  return(outdf)
})

plotdf <- pair_cor |>
  do.call(what = "rbind") |>
  dplyr::mutate(mod.index = 1:n())

p <- ggplot(plotdf, aes(x = mod.index, y = mean)) + 
  theme_classic() +
  geom_pointrange(aes(ymin = mean - 2*se, ymax = mean + 2*se), 
                  linewidth = 0.2,
                  size = 0.05,
                  alpha = 0.3,
                  shape = 16) +
  geom_smooth(color = "darkred",
              alpha = 0.3) +
    labs(x = "Module index (n = 1016)",
         y = "Correlation"#,
         #title = "Mean pairwise correlation of seed genes"
         ) +
  scale_x_continuous(breaks = c(1, 1016)) +
  ylim(0, 1) +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(0, 8, 0, 0)),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.justification = c("right", "bottom"),
        legend.box.background = element_rect(color = "black", linewidth = 1)) 
ggsave(p, file = file.path(save_dir, "panel_R.svg"), height = 2, width = 7)


##########
# Panel S
#########

# Gather mean REI
corlist <- lapply(these_mods_final, \(i){
  tempdf <- data.frame("lein" = unlist(lein[i, ]), 
                       "SEA" = unlist(SEAcon[i, ]), 
                       "mit" = unlist(mitcon[i, ]),
                       "mora" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = "pairwise.complete.obs")
  return(tempcor)
})
cormeans <- Reduce("+", lapply(corlist,
  function(m) { m[is.na(m)] <- 0; m })) /           
              Reduce("+", lapply(corlist,           
  function(m) !is.na(m)))
colnames(cormeans) <- c("Jorstad 2023", "Gabitto 2024", "Liu 2025", "Morabito 2021")
rownames(cormeans) <- c("Jorstad 2023", "Gabitto 2024", "Liu 2025", "Morabito 2021")

# > cormeans
#               Jorstad 2023 Gabitto 2024  Liu 2025 Morabito 2021
# Jorstad 2023     1.0000000    0.8970314 0.8807404     0.8134280
# Gabitto 2024     0.8970314    1.0000000 0.9144771     0.8600619
# Liu 2025         0.8807404    0.9144771 1.0000000     0.8153554
# Morabito 2021    0.8134280    0.8600619 0.8153554     1.0000000

# Panel S heatmap is rendered by the shared make_panel_S() (see panel_S_plot.R),
# which is also used by the standalone panel_S.R for fast iteration.
source(file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/v6"), "panel_S_plot.R"))
make_panel_S(cormeans, file.path(save_dir, "panel_S.svg"))


######
# find module to feature
######
proj <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/mod_means/log_REI/mod_means_Con_bulk_megaset.csv"))
proj <- proj[these_mods_final, ]

cv_func <- function(x){
  sd(x)/mean(x)
}

out_df <- data.frame("Module" = 1:nrow(proj),
                     "CV" = apply(proj, 1, cv_func)) 
cordf <- lapply(these_mods_final, \(i){
  tempdf <- data.frame("lein" = unlist(lein[i, ]), 
                       "SEA" = unlist(SEAcon[i, ]), 
                       "mora" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = "pairwise.complete.obs")
  return(data.frame("mod" = i, "cors" = tempcor[upper.tri(tempcor)]))
})
out_df$REI_mean_cor <- cordf |> do.call(what = "rbind") |> group_by(mod) |> summarise(mean = mean(cors)) |> pull(mean)

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

out_df$lm <- class_lm
out_df$ind_1158 <- these_mods_final
out_df$ind_1023 <- c(1:1023)[match(out_df$ind_1158, these_mods)]

cand <- out_df[out_df$CV < 0.15, ] |>
   arrange(desc(REI_mean_cor)) 

# 1023 indices)
# 284
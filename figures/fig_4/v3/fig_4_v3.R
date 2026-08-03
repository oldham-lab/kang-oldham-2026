# v2: 
# - Remove Liu projections
# - Remove individual cor heatmaps

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

version <- "v3"
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

mod_id <- c(65, 669, 681, 1007) #1158 indices
# Other candidates for first row: (1023 index)
# 76 (spliceosome)

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

# all_plots[[3]] <- lapply(all_plots_full[[3]], \(x){
#   x2 <- x$data |> 
#     slice(1:5) |>
#     dplyr::mutate(SetName=factor(SetName, levels=rev(unique(SetName)))) |>
#     ggplot(aes(x = SetName, y = pval)) +
#         newtheme + 
#         geom_bar(stat="identity")
#   x2$layers[[2]] <- x$layers[[2]]
#   x2 <- x2 + coord_flip() +
#     theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5, size = 10, color = "black"),
#           axis.text.y = element_text(hjust = 1,vjust = 0.5, size = 10, color = "black"),
#           axis.title.y = element_text(size = 11, color = "black"),
#           axis.title.x = element_text(size = 11, color = "black"),
#           plot.title = element_blank()) +
#     labs(x = "", y = bquote(-log[10]~"(p-val)"), title = "Geneset enrichment") #+
#   return(x2)
# })
 
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

# mit <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
#   dplyr::select(!module)
# mit_se <- fread(data.table = F, file = paste0(base_dir, "/MIT_ADMultiome/sn_proj_indices/", index_dirs[d], "/indices_se_RNA.SubclassCon.csv")) 
# mit_se <- mit_se[, match(colnames(mit), colnames(mit_se))]

SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
  dplyr::select(!module)
SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) 
SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]

mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control.csv")) |>
  dplyr::select(!module)
mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) 
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
projplotlist <- lapply(mod_id, \(i){
  n <- length(allcts)
  plotdf <- data.frame("ct" = rep(allcts, 3),
                       "dataset" = c(rep("Jorstad et al. 2023 (PFC)", n), 
                                     rep("Gabitto et al. 2024 (control PFC)", n), 
                                     rep("Morabito et al. 2021 (control PFC)", n)),
                       "ind" = c(unlist(lein[i, ]), 
                                 unlist(SEAcon[i, ]), 
                                 unlist(mora[i, ])),
                       "ind_se" = c(unlist(lein_se[i, ]), 
                                    unlist(SEAcon_se[i, ]), 
                                    unlist(mora_se[i, ]))) |>
    dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                  ct = factor(ct, levels = allcts))

  # Calculate and print pairwise cor mat
  cat(i, "pairwise cors: \n")
  tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), 
                        "Gabitto 2024" = unlist(SEAcon[i, ]), "Morabito 2021" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
  colnames(tempcor) <- c("Jorstad 2023","Gabitto 2024", "Morabito 2021")
  rownames(tempcor) <- c("Jorstad 2023", "Gabitto 2024", "Morabito 2021")
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
  max_setname <- lapply(all_plots[[3]][[1]]$data$SetName |> as.character(), nchar) |> unlist() |> max()
  ggsave(all_plots[[3]][[j]], file = file.path(save_dir, paste0("panel_", j, "_3.pdf")),
          height = 1.5, width = (3 + max_setname/30 * 2.5))
  ggsave(all_plots[[3]][[j]], file = file.path(save_dir, paste0("panel_", j, "_3.svg")),
          height = 1.5, width = (3 + max_setname/30 * 2.5))

  # Projection indices:
  ggsave(projplotlist[[j]], file = file.path(save_dir, paste0("panel_", j, "_4.pdf")), height = 4, width = 5.5)
  ggsave(projplotlist[[j]], file = file.path(save_dir, paste0("panel_", j, "_4.svg")), height = 4, width = 5.5)


}

########## 
# Panel Q:
# Module size histogram
##########

modlengthdf <- data.frame("length" = modulelengths[these_mods_final], "VE" = vevec[these_mods_final])

p <- ggplot(modlengthdf, aes(x = length)) +
  theme_classic() +
  geom_histogram() +
  labs(x = "Module size (topmodposbc)", y = "Count") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(p, file = file.path(save_dir, "panel_Q_size_hist.svg"), bg = "white", width = 4, height = 2.5)

##########
# Panel R:
# Module size vs VE
#########

p <- ggplot(modlengthdf, aes(x = length, y = VE)) +
  theme_classic() +
  geom_point(size = 0.4) +
  labs(x = "Module size (topmodposbc)", y = "% Variance explained by PC1") #+
  #scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(p, file = file.path(save_dir, "panel_R_size_vs_VE.svg"), bg = "white", width = 4, height = 2.5)


##########
# Panel S: 
# Summary figure for pairwise correlations of indices (n = 3 datasets) for all modules
# v3: remove coloring
##########
# Gather pairwise cors
cordf <- lapply(1:nrow(lein), \(i){
  tempdf <- data.frame("lein" = unlist(lein[i, ]), 
                      # "mit" = unlist(mit[i, ]), 
                      # "scope" = unlist(scope[i, ]), 
                       "SEA" = unlist(SEAcon[i, ]), 
                       "mora" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = "pairwise.complete.obs")
  return(data.frame("mod" = i, "cors" = tempcor[upper.tri(tempcor)]))
})
cordforder <- cordf |> do.call(what = "rbind") |> group_by(mod) |> summarise(mean = mean(cors)) |> arrange(mean) # Arrange by mean pairwise cor

# Plot as means (points) plus error bars 
cordfmean <- cordf |> do.call(what = "rbind") |> 
  left_join(these_mods_final_df, by = join_by(mod == old)) |>
  filter(mod %in% these_mods_final) |>
  group_by(new) |> 
  summarise(mean = mean(cors),
            se = sd(cors)/sqrt(n())) |> arrange(mean) |>
  #inner_join(out_df2 |> select(Module, quad), by = join_by(mod == Module)) |>
  mutate(new = factor(new, levels = new))  #|>
  #filter(!is.na(quad)) # 1 mod
p <- ggplot(cordfmean, aes(x = new, y = mean#, color = quad
      )) + 
  theme_classic() +
  geom_pointrange(aes(ymin = mean - 2*se, ymax = mean + 2*se), 
                  linewidth = 0.2,
                  size = 0.05,
                  alpha = 0.3) +
  ggrepel::geom_label_repel(data = cordfmean[cordfmean$new %in% mod_labels_final, ], aes(label = new),
                  box.padding = 0.5, max.overlaps = Inf, color = "Black", force_pull = 0.05, size = 6) +
  labs(x = paste0("Module (n = (", length(these_mods_final), ")"), 
       y = "Mean pairwise correlation"#, 
       #title = "Pearson correlation"
       ) +
  ylim(0, 1) +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(0, 8, 0, 0)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.justification = c("right", "bottom"),
        legend.box.background = element_rect(color = "black", linewidth = 1)) +
  #scale_color_manual(labels = c("left.bottom" = "Non-specific", "right.bottom" = "Celltype-specific", "left.top" = "Class\nspecific", "right.top" = "Non-neuronal"),
   #                  values = c("left.bottom" = pal3[1], "right.bottom" = pal3[2], "left.top" = pal3[3], "right.top" = pal3[4])) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave(p, file = file.path(save_dir,"panel_S.svg"), width = 6, height = 4, bg = "white")
ggsave(p, file = file.path(save_dir, "panel_S.pdf"), width = 6, height = 4, bg = "white")

# # Save table of CV and class r2
# out_df3 <- out_df2 |>
#   mutate(label = case_match(quad,
#     "left.top" ~ "Class-specific",
#     "left.bottom" ~ "Non-specific",
#     "right.bottom" ~ "Celltype-specific",
#     "right.top" ~ "Non-neuronal",
#     .default = NA))
# colnames(out_df3)[4] <- "pcnt_VE_by_mod_eig"
# fwrite(out_df3, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_4/CV_table.csv"))









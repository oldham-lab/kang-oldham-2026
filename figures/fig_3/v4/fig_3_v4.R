# v4

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

version_number  <- "v4"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_3/"), version_number)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

## Overlay nested dataset-correlation brackets + r-value labels in the right
## margin of a projection-index panel.
##   p      : the projection ggplot (saved at 6.5 x 6 in -> 468 x 432 pt)
##   cormat : 4x4 correlation matrix in dataset order
##            (Jorstad, Gabitto, Liu, Morabito) -- e.g. tempcor below
## One bracket is drawn per dataset pair (C(4,2)=6), nested. The bracket
## geometry was recovered from the original Inkscape annotation in group-local
## coords; matrix(0.75,0,0,0.75,609.7666,-162.0675) maps it into panel E's
## 468x432 pt plot space. base_width_in is the panel's own saved width before
## the annotation margin (6.5in for panel E, 6.8in for panel F); a panel saved
## wider has its plot edge further right by exactly (base_width_in-6.5)*72 pt,
## so the brackets are shifted by the same amount (the facets are identically
## positioned vertically, so y is unchanged). The panel keeps its base width;
## the extra (width_in - base_width_in) is blank margin holding the annotation.
## Returns a cowplot ggdraw object; save it with width = width_in, height = height_in.
add_corr_brackets <- function(p, cormat, base_width_in = 6.5, width_in = 7.5, height_in = 6) {
  W <- width_in * 72; H <- height_in * 72            # canvas in pt
  right_margin_pt <- W - base_width_in * 72 + 5.5    # keep panel, add blank right
  x_shift <- base_width_in * 72 - 468                # plot-edge offset vs panel E
  gap <- 2; font <- 10                               # label gap (pt) and size
  bxs <- 4   # the rightmost two brackets (+ their labels) shift this many pt right,
             # opening space to the 3rd-from-right label without moving any text
             # relative to its own bracket (all labels stay equally spaced).
  bl <- list(  # mx,my = top-left; dx,dy = serif/spine lengths; dx2 = lower serif;
               # bx = per-bracket horizontal shift in pt (bracket + its own label move together)
    list(mx=-189.22643, my=299.20186, dx=4.99878, dy= 94.89874, dx2=5.41535, i=1, j=2, bx=0),
    list(mx=-189.55725, my=401.63598, dx=5.21998, dy=108.57642, dx2=5.65498, i=2, j=3, bx=0),
    list(mx=-190.02781, my=517.12821, dx=5.23063, dy=105.59946, dx2=5.66652, i=3, j=4, bx=0),
    list(mx=-167.87949, my=299.14475, dx=4.94123, dy=210.68423, dx2=5.35301, i=1, j=3, bx=0),
    list(mx=-147.07679, my=402.34114, dx=4.92456, dy=220.06389, dx2=5.33494, i=2, j=4, bx=bxs),  # 2nd from right
    list(mx=-126.99299, my=299.18036, dx=4.76309, dy=323.60958, dx2=5.16003, i=1, j=4, bx=bxs))  # rightmost
  sx <- function(lx) 0.75 * lx + 609.7666 + x_shift  # local -> pt
  sy <- function(ly) 0.75 * ly - 162.0675
  xn <- function(xpt) xpt / W                        # pt -> npc (ggdraw coords)
  yn <- function(ypt) 1 - ypt / H
  g <- ggdraw(p + theme(plot.margin = margin(5.5, right_margin_pt, 5.5, 5.5, "pt")))
  for (b in bl) {
    px <- sx(c(b$mx, b$mx + b$dx, b$mx + b$dx, b$mx + b$dx - b$dx2)) + b$bx
    py <- sy(c(b$my, b$my, b$my + b$dy, b$my + b$dy))
    g <- g + draw_line(xn(px), yn(py), color = "black", size = 0.4)
    # truncate toward zero so a near-perfect r (e.g. 0.996) shows "0.99", not "1.00"
    r_disp <- trunc(cormat[b$i, b$j] * 100 + sign(cormat[b$i, b$j]) * 1e-9) / 100
    g <- g + draw_label(sprintf("r = %.2f", r_disp),
                        x = xn(sx(b$mx + b$dx) + gap + font / 2 + b$bx),
                        y = yn(sy(b$my + b$dy / 2)),
                        angle = 270, size = font, fontfamily = "sans")
  }
  g
}

# Load bulk-related objects
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]
modulelengths <- modulelengths[these_mods]

mod_index <- 354 # 1023

# sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
# these_mods_final <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
# 1016 index is also 354

# Calc % VE
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
  x <- all_plots_full[[1]][[i]]
  x$layers[[1]] <- NULL

  x + 
    geom_line(linewidth = 0.2) +
    labs(x = "Sample") + 
    newtheme + 
    theme(text = element_text(size = 14),
          legend.box.margin = margin(0, 0, 0, -10),
          plot.title = element_text(size = 14, hjust = 0.5),# margin = margin(0, 0, -0.1, 0)),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_color_manual(values = cols2) +
    labs(title = paste0("Module ", mod_index),
         subtitle = j,
         colour = "Top 10 genes:"
         #subtitle = bquote(subt)
         )
}, seq_along(all_plots_full[[1]]), paste0("# of module genes: ", modulelengths[mod_index], "\n", sprintf("%.0f", vevec[mod_index] * 100), "% variance explained by module eigengene"), SIMPLIFY = F)

plots_bc <- list()
plots_bc[[1]] <- mapply(\(i, j){
  #subt <- as.list(as.character(all_plots[[1]][[i]]$labels$subtitle)) |> do.call(what = paste0)
  x <- plots_bc_full[[1]][[i]]
  x$layers[[1]] <- NULL

  x + 
    geom_line(linewidth = 0.2) +
    labs(x = "Sample") + 
    newtheme + 
    theme(legend.box.margin = margin(0, 0, 0, -10),
          plot.title = element_text(size = 14, hjust = 0.5),# margin = margin(0, 0, -0.1, 0)),
          plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_color_manual(values = cols2) +
    labs(title = paste0("Module ", mod_index),
         subtitle = j,
         colour = "Top 10 genes:"
         #subtitle = bquote(subt)
         )
}, seq_along(plots_bc_full[[1]]), paste0("# of module genes: ", modulelengths[mod_index], "\n", sprintf("%.0f", vevec[mod_index] * 100), "% variance explained by module eigengene"), SIMPLIFY = F)



all_plots[[2]] <- lapply(all_plots_full[[2]], \(x){
   x$data <- x$data |> 
    mutate(split = case_match(split,
      "BrainGVEX (n=255)" ~ "BrainGVEX (n=275)",
      "Brainseq (n=182)" ~ "Brainseq (n=190)",
      "CMC (n=274)" ~ "CMC (n=285)",
      "CMC_HBCC (n=145)" ~ "CMC_HBCC (n=162)",
      "NABEC (n=63)" ~ "NABEC (n=69)",
      .default = split))
  
  sig_df <- x$data |>
    group_by(split) |>
    summarise(
        p = wilcox.test(cors[type == levels(factor(type))[1]],
                        cors[type == levels(factor(type))[2]])$p.value,
        .groups = "drop"
    ) |>
    mutate(
        label = case_when(
            p < 0.0001 ~ "****",
            p < 0.001 ~ "***",
            p < 0.01  ~ "**",
            p < 0.05  ~ "*",
            TRUE      ~ "ns"
        ),
        y = 1.07   # x-axis position after coord_flip; adjust to sit above violins
    )
  dodge <- position_dodge(width = 0.6)
  x$data |>
    ggplot(aes(x = split, y = cors, fill = type)) + theme_classic() +
    geom_violin(position = dodge) +
    geom_boxplot(notch = T, width = 0.1, position = dodge, show.legend = F, outlier.shape = NA) +
    geom_text(data = sig_df, aes(x = split, y = y, label = label),
              inherit.aes = FALSE, size = 5) +
    ggtitle("Distributions of bulk pairwise correlations") +
    ylab("Pearson correlation") +
    xlab("") +
    scale_x_discrete(labels = label_wrap(10)) +
    theme(plot.title = element_text(size=13, hjust=0.5),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.direction = "vertical",
          legend.margin = margin(-4, 0, 0, 0),
          legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.4, "cm"),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=11),
          axis.text.y = element_text(size=11),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c("salmon", "grey")) +
    ylim(NA, 1.1)# +
    #coord_flip()
})

all_plots[[3]] <- lapply(all_plots_full[[3]], \(x){
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

  x2 <- x2 + 
    coord_flip() +
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
index_dirs <- c("log_native", "log_REI")
index_list <- list()
index_save_names <- c("native_log", "REI")
index_xaxis_names <- c("Mean expression (log UMI counts + 1)", "Relative expression index (REI)\n(Normalized to mean genome-wide\nexpression in each subclass)")
for(d in 1:2){
  lein <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Cell_Type_1_topmodposbc_mean.csv")) |>
    dplyr::select(!module)
  lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_se_topmodposbc.csv")) 
  lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]

  mitcon <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/mod_means/"), index_dirs[d], "mod_means_Con_bulk_megaset.csv"))
  mitse <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/se/"), index_dirs[d], "se_Con_bulk_megaset.csv"))

  SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv")) |>
    dplyr::select(!module)
  SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_AllADVsCon_DFC/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) 
  SEAcon_se <- SEAcon_se[, match(colnames(SEAcon), colnames(SEAcon_se))]

  mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control_topmodposbc_mean.csv")) |>
    dplyr::select(!module)
  mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassContopmodposbc.csv")) 
  mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]

  # Change celltypes to match
  colnames(lein)[c(1, 3, 15, 16)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
  colnames(mora)[c(1, 3, 14, 15)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")

  colnames(lein_se) <- colnames(lein)
  # colnames(SEAcon_se) <- colnames(SEAcon) # already matching
  colnames(mora_se) <- colnames(mora)

  allcts <- c("L2/3 IT", "L4 IT","L5 IT",  "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3", 
            "Chandelier", "Pvalb", "Sst", "Sst Chodl", "Lamp5 Lhx6", "Lamp5",  "Pax6",  "Sncg",  "Vip",
            "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")


  allcts_cap <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3",
            "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5",  "PAX6",  "SNCG",  "VIP",
            "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")
  proj_pal <- setNames(
    colorRampPalette(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                       "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"))(length(allcts_cap)),
    allcts_cap
  )

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
  projplotlist <- list()
  i = mod_index
  n <- length(allcts)
  plotdf <- data.frame("ct" = rep(allcts_cap, 4),
                        "dataset" = c(rep("Jorstad et al. 2023 (normal human dorsolateral prefrontal cortex)", n), 
                                      rep("Gabitto et al. 2024 (normal human dorsolateral prefrontal cortex)", n), 
                                      rep("Liu et al. 2025 (normal human dorsolateral prefrontal cortex)", n), 
                                      rep("Morabito et al. 2021 (normal human dorsolateral prefrontal cortex)", n)),
                        "ind" = c(unlist(lein[these_mods[i], ]), 
                                  unlist(SEAcon[these_mods[i], ]), 
                                  unlist(mitcon[these_mods[i], ]), 
                                  unlist(mora[these_mods[i], ])),
                        "ind_se" = c(unlist(lein_se[these_mods[i], ]), 
                                     unlist(SEAcon_se[these_mods[i], ]), 
                                     unlist(mitse[these_mods[i], ]), 
                                     unlist(mora_se[these_mods[i], ]))) |>
    dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                  ct = factor(ct, levels = allcts_cap))

  projplotlist[[i]] <- ggplot(plotdf, aes(x = ct, y = ind, fill = ct)) +
    theme_classic() +
    geom_col(position = position_dodge(), alpha = 0.8) +
    geom_errorbar(aes(ymin = ind - 2 * ind_se,
                      ymax = ind + 2 * ind_se),
                      width = 0.2,
                      linewidth = 0.3,
                      position = position_dodge(0.5)) +
    scale_fill_manual(values = proj_pal, na.value = "grey70") +
    theme(text = element_text(family = "sans", color = "black", size = 12),
          legend.position="none", 
          axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(color = "black"),
          strip.background = element_rect(fill = "white")) +
    labs(y = index_xaxis_names[d], x = "")
    
    if(d == 1){
      projplotlist[[i]] <- projplotlist[[i]] +       
        facet_wrap(~dataset, ncol = 1, nrow = 4, scales = "free_y") 
    } else {
      projplotlist[[i]] <- projplotlist[[i]] +       
        facet_wrap(~dataset, ncol = 1, nrow = 4) 
    }
    cat(i, " ")

  ## Calculate pairwise cors
  tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), 
                       "Gabitto 2024" = unlist(SEAcon[i, ]), 
                       "Liu 2025" = unlist(mitcon[i, ]),
                       "Morabito 2021" = unlist(mora[i, ]))
  tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
  colnames(tempcor) <- c("Jorstad 2023", "Gabitto 2024", "Liu 2025", "Morabito 2021")
  rownames(tempcor) <- c("Jorstad 2023", "Gabitto 2024", "Liu 2025", "Morabito 2021")
  print(tempcor)

  ## Plot individual panels
  outindices <- 1:length(all_plots[[1]])
  for(j in outindices){
    # Expression plot:
    ggsave(all_plots[[1]][[j]] +
           labs(x = "Bulk sample (n = 100 random samples)") +
           theme(legend.position = "bottom") +
           guides(color = guide_legend(title.position = "top", 
                                       title.hjust = 0.5,
                                       override.aes = list(linewidth = 2))), 
           file = file.path(save_dir, "panel_B_seed.svg"),
           width = 5.5, height = 4)

    ggsave(plots_bc[[1]][[j]] +
           labs(x = "Bulk sample (n = 100 random samples)") +
           theme(legend.position = "bottom") +
           guides(color = guide_legend(title.position = "top", 
                                       title.hjust = 0.5,
                                       override.aes = list(linewidth = 2))), 
           file = file.path(save_dir, "panel_B_bc.svg"),
           width = 5.5, height = 4)

    # Bulk cor plot:
    ggsave(all_plots[[2]][[j]] +
             theme(plot.title = element_blank()#,
                   #plot.margin = margin(1,2,1,1,"cm")
                   ) +
             guides(fill = guide_legend(override.aes = list(shape = 22))), 
           file = file.path(save_dir, "panel_C.svg"),
           width = 5.5, height = 3)

    # GSEA plot:
    max_setname <- lapply(all_plots[[3]][[1]]$data$SetName |> as.character(), nchar) |> unlist() |> max()
    ggsave(all_plots[[3]][[j]],# +
            # theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1)), 
           file = file.path(save_dir, "panel_D.svg"),
           height = 1.5, width = (2 + max_setname/30 * 2))

    # Projection indices:
    if(d==1){
      ggsave(projplotlist[[mod_index]],
             file = file.path(save_dir, paste0("panel_E_", index_save_names[d], ".svg")),
             height = 6,
             width = 6.5)
      ## Same panel with nested dataset-correlation brackets + r labels baked in.
      ## Correlations are taken live from tempcor (computed just above).
      ggsave(add_corr_brackets(projplotlist[[mod_index]], tempcor),
             file = file.path(save_dir, paste0("panel_E_", index_save_names[d], "_with_brackets.svg")),
             height = 6,
             width = 7.5)
    } else if(d==2){
      ggsave(projplotlist[[mod_index]],
             file = file.path(save_dir, paste0("panel_F_", index_save_names[d], ".svg")),
             height = 6,
             width = 6.8)
      ## Same panel with nested dataset-correlation brackets + r labels baked in.
      ggsave(add_corr_brackets(projplotlist[[mod_index]], tempcor, base_width_in = 6.8, width_in = 7.8),
             file = file.path(save_dir, paste0("panel_F_", index_save_names[d], "_with_brackets.svg")),
             height = 6,
             width = 7.8)
    }
  }
}

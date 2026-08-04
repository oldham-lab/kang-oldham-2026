# v1.3
# - prepping code to create snapshots directly from data
#   instead of loading from premade .qs objects
# - to be incorporated into shiny app (app.R)

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

## Load objects
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")

# Load individual bulk cor data
bulk_plat_cors <- fread(file.path(module_output_dir,"bulk_cors_by_dataset","bulk_plat_cors.csv"), data.table=F)
bulk_rand_cors <- qread(file.path(module_output_dir,"bulk_cors_by_dataset","bulk_rand_cors.qs"))
bulk_rand_cors$type <- gsub("^Random genes, n = (.+)$", "Random genes\n(n = \\1)", bulk_rand_cors$type)

# Load GSEA-related data
gsea_cutFDR <- 4.695451
gsea_plot_dfs <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/shiny_code/gsea_plot_dfs.qs"))

# save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")

# Load module information
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods_1023 <- as.numeric(names(mods)[which(modulelengths>filter_under)])
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
ind_1016 <- which(!these_mods_1023 %in% which(sigcount_bonf$vals < 2))
these_mods <- these_mods_1023[ind_1016]
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
modulelengths_seed <- unlist(lapply(mod_seed, length))
mod_id <- 1:1158

# # Calc % VE (seed)
# expr_t <- t(bulk_expr[,3:ncol(bulk_expr)])
# colnames(expr_t) <- bulk_expr[,2]
# expr_z <- apply(expr_t, 2, function(x) (x-mean(x))/sd(x))
# r2_vec_bc_list <- list()
# for(j in seq_along(mods)){
#   expr_plot <- as.data.frame(expr_z[ ,colnames(expr_z) %in% mods[[j]]])
#   pc1_bc <-  prcomp(expr_plot[apply(expr_plot,2,var)>0], scale=T)$x[,1]
#   r2_vec_bc <- apply(expr_plot[apply(expr_plot,2,var)>0],2,function(x) summary(lm(x~pc1_bc))$r.squared)
#   r2_vec_bc_list[[j]] <- r2_vec_bc
#  # cat(j, " ")
# }
# vevec <- lapply(r2_vec_bc_list, mean) |> unlist()
# qsave(vevec, file = "/home/gugene/ShinyApps/CoPA/www/mod_eig_var_explained.qs")
vevec <- qread(file.path(Sys.getenv("SHINYAPP_DIR", "/home/gugene/ShinyApps/CoPA"), "www/CoPA_files/mod_eig_var_explained.qs"))

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

# all_plots[[2]] <- lapply(all_plots_full[[2]], \(x){
#   x + newtheme +
#     theme(legend.text = element_text(size = 8),
#           plot.title = element_text(hjust = 0.5, size = 14),
#           axis.text.x = element_text(size = 12),
#           legend.position = "bottom",
#           legend.direction = "vertical",
#           legend.margin = margin(-10, 0, 0, 0),
#           plot.margin = margin(1, 0.3, 0.6, 0.3, "cm")) +
#     guides(fill = guide_legend(override.aes = list(shape = 22)))
# })


# Produce expression line plot input dataframes
expr_t <- t(expr[,3:ncol(expr)])
colnames(expr_t) <- expr[,2]
expr_z <- apply(expr_t, 2, function(x) (x-mean(x))/sd(x))

expr_plots <- list()
expr_plots_seed <- list()
for(i in seq_along(these_mods)){
  j <- these_mods[i]
  #cat(j, " ")
  # collect top 10 genes by kme (topmodposbc)
  mod_kme <- datkme[, c(1:3, j+4)]
  mod_kme <- mod_kme[order(mod_kme[, 4], decreasing=T), ]
  mod_kme <- mod_kme[which(mod_kme$topmodposbc == j), ]
  
  # collect top 10 seed genes ordered by kme
  seed_kme <- datkme[,c(1:2,j+4)]
  seed_kme <- seed_kme[which(seed_kme[,2] %in% mod_seed[[j]]),]
  seed_kme <- seed_kme[order(seed_kme[,3], decreasing=T),]
  
  # Select random samples to graph
  samp_ind1 <- sample(1:nrow(expr_z), 100)
  expr_plot <- as.data.frame(expr_z[samp_ind1, colnames(expr_z) %in% mods[[which(names(mods)==j)]]])
  expr_plot_seed <- as.data.frame(expr_z[samp_ind1, colnames(expr_z) %in% mod_seed[[j]]])
  while(all(apply(expr_plot, 2, var) == 0) | all(apply(expr_plot_seed, 2, var) == 0)){
    samp_ind1 <- sample(1:nrow(expr_z), 100)
    expr_plot <- as.data.frame(expr_z[samp_ind1, colnames(expr_z) %in% mods[[which(names(mods)==j)]]])
    expr_plot_seed <- as.data.frame(expr_z[samp_ind1,colnames(expr_z) %in% mod_seed[[j]]])
  }
  
  # Select gene expr vectors and calculate mean r2 for topmodposbc graph
  pc1_bc <-  prcomp(expr_plot[apply(expr_plot,2,var)>0], scale=T)$x[,1]
  r2_vec_bc <- apply(expr_plot[apply(expr_plot,2,var)>0],2,function(x) summary(lm(x~pc1_bc))$r.squared)
  mean_r2_bc <- paste0(signif(mean(r2_vec_bc),2) *100, "%")
  expr_plot <- expr_plot[,colnames(expr_plot) %in% mod_kme[1:min(10, nrow(seed_kme)),2]]
  expr_plot$sample <- c(1:nrow(expr_plot))
  expr_plot <- pivot_longer(as.data.frame(expr_plot),cols=!sample, names_to = "gene", values_to = "value")
  expr_plot$gene <- factor(expr_plot$gene, levels=mod_kme[1:length(unique(expr_plot$gene)),2])
  
  # Select gene expr vectors and calculate mean r2 for seed gene graph
  pc1 <-  prcomp(expr_plot_seed[apply(expr_plot_seed,2,var)>0], scale=T)$x[,1]
  r2_vec <- apply(expr_plot_seed[apply(expr_plot_seed,2,var)>0],2,function(x) summary(lm(x~pc1))$r.squared)
  mean_r2_seed <- paste0(signif(mean(r2_vec),2) *100, "%")
  expr_plot_seed <- expr_plot_seed[,colnames(expr_plot_seed) %in% seed_kme[1:min(10, nrow(seed_kme)),2]]
  expr_plot_seed$sample <- c(1:nrow(expr_plot_seed))
  expr_plot_seed <- pivot_longer(as.data.frame(expr_plot_seed),cols=!sample, names_to = "gene", values_to = "value")
  expr_plot_seed$gene <- factor(expr_plot_seed$gene, levels=seed_kme[1:length(unique(expr_plot_seed$gene)),2])
  
  expr_plots[[i]] <- expr_plot
  expr_plots_seed[[i]] <- expr_plot_seed
  if (j %% 100 == 0) print(j)
}

saveoutstring <- file.path(file.path(Sys.getenv("SHINYAPP_DIR", "/home/gugene/ShinyApps/CoPA"), "www/CoPA_files/expr_line/"), "topmodposbc.qs")
qsave(expr_plots, file = saveoutstring)

saveoutstring <- file.path(file.path(Sys.getenv("SHINYAPP_DIR", "/home/gugene/ShinyApps/CoPA"), "www/CoPA_files/expr_line/"), "seed.qs")
qsave(expr_plots_seed, file = saveoutstring)


# Precalculate bulk cor dfs
bulk_plat_split <- split(bulk_plat_cors, bulk_plat_cors$mod)
bulk_rand_split  <- split(bulk_rand_cors,  bulk_rand_cors$mod)
mod_len          <- setNames(modulelengths, names(mods))
bulk_cor_list <- pbmcapply::pbmclapply(
  setNames(these_mods, these_mods),
  function(i) {
    bulk_cor <- bulk_plat_split[[as.character(i)]]
    bulk_cor$type <- paste0("Module genes\n(n = ", mod_len[[as.character(i)]],")")
    all_megacors <- bulk_rand_split[[as.character(i)]] %>%
      dplyr::group_by(split, type) %>%
      dplyr::slice_sample(n = 1000, replace = FALSE) %>%
      dplyr::ungroup()
    bulk_cor <- rbind(bulk_cor, all_megacors)
    bt1 <- unique(bulk_cor$type)
    bulk_cor$type <- factor(bulk_cor$type, levels = c(bt1[grep("Module", bt1)],
                                                      bt1[grep("Random", bt1)]))
    bulk_cor
  },
  mc.cores = parallel::detectCores() - 1
)
bulk_cor_list <- lapply(bulk_cor_list, \(x){
  x |> mutate(split = case_match(split,
    "BrainGVEX (n=255)" ~ "BrainGVEX (n=275)",
    "Brainseq (n=182)" ~ "Brainseq (n=190)",
    "CMC (n=274)" ~ "CMC (n=285)",
    "CMC_HBCC (n=145)" ~ "CMC_HBCC (n=162)",
    "NABEC (n=63)" ~ "NABEC (n=69)",
    .default = split))
})

# Pre-compute per-split Wilcoxon significance labels so the Shiny app
# can use geom_text instead of running stat_compare_means on every render
bulk_cor_list <- lapply(bulk_cor_list, function(df) {
  sig_df <- df |>
    dplyr::group_by(split) |>
    dplyr::summarise(
      p_value = wilcox.test(
        cors[grepl("Module", type)],
        cors[grepl("Random", type)],
        alternative = "greater"
      )$p.value,
      .groups = "drop"
    ) |>
    dplyr::mutate(sig_label = dplyr::case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )) |>
    dplyr::select(split, sig_label)
  dplyr::left_join(df, sig_df, by = "split")
})

qsave(bulk_cor_list, file = file.path("/home", "gugene", "ShinyApps", "CoPA", "www", "CoPA_files", "bulk_cor", "bulk_cor_list.qs"))
# NOTE: the app no longer loads this .qs (or the gene-projection / core_gbmap
# CSVs) into RAM. After these source files are in place, run build_fst_store.R
# to package the on-disk fst store (www/fst/) that app.R actually queries.


# Produce projection and correlation matrices for input
base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")
index_dirs <- c("log_native", "log_REI")
index_save_names <- c("native_log", "REI")
index_xaxis_names <- c("Mean expression (log UMI counts)", "Relative expression index")
for(def in c("seed", 
             "topmodposbc"
             )){
  modtoken <- if (def == "seed") "bulk_megaset_seed" else "bulk_megaset"
  for(d in 1:2){ # native log vs REI
    dstring <- index_save_names[d]
 
    ## Load projection data
    lein <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Cell_Type_1_", def, "_mean.csv")) |>
      dplyr::select(!module)
    lein_se <- fread(data.table = F, file = paste0(base_dir, "/LeinDFC/sn_proj_indices/", index_dirs[d], "/indices_se_", def, ".csv"))
    lein_se <- lein_se[, match(colnames(lein), colnames(lein_se))]

    mitcon <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/mod_means/"), index_dirs[d], paste0("mod_means_Con_", modtoken, ".csv")))
    mitse <- fread(data.table = F, file = file.path(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/se/"), index_dirs[d], paste0("se_Con_", modtoken, ".csv")))

    SEAcon <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/mod_means/", index_dirs[d], "/mod_means_Con_", modtoken, ".csv"))
    SEAcon_se <- fread(data.table = F, file = paste0(base_dir, "/SEAAD2024_full_python_output/PFC/se/", index_dirs[d], "/se_Con_", modtoken, ".csv"))

    mora <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_over_all_datasets_Subclass_Control_", def, "_mean.csv")) |>
      dplyr::select(!module)
    mora_se <- fread(data.table = F, file = paste0(base_dir, "/Morabito_ABIanno/sn_proj_indices/", index_dirs[d], "/indices_se_SubclassCon", def, ".csv"))
    mora_se <- mora_se[, match(colnames(mora), colnames(mora_se))]

    # Change celltypes to match

    colnames(lein)[c(1, 3, 15, 16)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")
    colnames(mora)[c(1, 3, 14, 15)] <- c("Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte")

    colnames(lein_se) <- colnames(lein)
    colnames(SEAcon_se) <- colnames(SEAcon)
    colnames(mora_se) <- colnames(mora)

    allcts <- c("L2/3 IT", "L4 IT","L5 IT",  "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3", 
            "Chandelier", "Pvalb", "Sst", "Sst Chodl", "Lamp5 Lhx6", "Lamp5",  "Pax6",  "Sncg",  "Vip",
            "Astrocyte", "Endothelial", "Microglia-PVM", "Oligodendrocyte", "OPC", "VLMC")


    allcts_cap <- c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 CT", "L6b", "L6 IT Car3", 
            "Chandelier", "PVALB", "SST", "SST CHODL", "LAMP5 LHX6", "LAMP5",  "PAX6",  "SNCG",  "VIP",
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
    projplotlist <- list()
    for(j in seq_along(these_mods)){
      i <- these_mods[j]
      n <- length(allcts)
      plotdf <- data.frame("ct" = rep(allcts_cap, 4),
                    "dataset" = c(rep("Jorstad et al. 2023 (normal human dorsolateral prefrontal cortex)", n), 
                                  rep("Gabitto et al. 2024 (normal human dorsolateral prefrontal cortex)", n), 
                                  rep("Liu et al. 2025 (normal human dorsolateral prefrontal cortex)", n), 
                                  rep("Morabito et al. 2021 (normal human dorsolateral prefrontal cortex)", n)),
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
      projplotlist[[j]] <- plotdf
      if (j %% 100 == 0) print(j)
    }
    saveoutstring_proj <- file.path(file.path(Sys.getenv("SHINYAPP_DIR", "/home/gugene/ShinyApps/CoPA"), "www/CoPA_files/projs/"), paste0(def, "_", dstring, ".qs"))
    qsave(projplotlist, file = saveoutstring_proj)

    # Plot individual cor heatmaps for each module
    corlistind <- list()
    for(j in seq_along(these_mods)){
      i <- these_mods[j]
      tempdf <- data.frame("Jorstad 2023" = unlist(lein[i, ]), 
                    "Gabitto 2024" = unlist(SEAcon[i, ]), 
                    "Liu 2025" = unlist(mitcon[i, ]),
                    "Morabito 2021" = unlist(mora[i, ]))
      tempcor <- cor(tempdf, use = 'pairwise.complete.obs')
      colnames(tempcor) <- c("Jorstad 2023","Gabitto 2024", "Liu 2025", "Morabito 2021")
      rownames(tempcor) <- c("Jorstad 2023", "Gabitto 2024", "Liu 2025", "Morabito 2021")
      testvec <- apply(tempcor, 2, \(x) sum(is.na(x)))
      if(any(testvec > 0)){
        rem <- which(testvec > 0)
        tempcor <- tempcor[-rem, -rem]
      } 
      corlistind[[j]] <- tempcor
      if (j %% 100 == 0) print(j)
    }

    saveoutstring <- file.path(file.path(Sys.getenv("SHINYAPP_DIR", "/home/gugene/ShinyApps/CoPA"), "www/CoPA_files/proj_cors/"), paste0(def, "_", dstring, ".qs"))
    qsave(corlistind, file = saveoutstring)
    
  }
}

 



# ## Produce full snapshots
# out_dir <- "/home/gugene/test/shiny_test"

# base_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/")
# index_dirs <- c("log_native", "log_REI")
# index_save_names <- c("native_log", "REI")
# index_xaxis_names <- c("Mean expression (log UMI counts)", "Relative expression index")
# for(def in c("seed", 
#              "topmodposbc"
#              )){
#   cat(def, "\n")
#   for(d in 1:2){ # native log vs REI
#     if(d == 1){
#       dstring <- "native_log"
#     } else {
#       dstring <- "REI"
#     }

#     # Load expression line plot
#     if(def == "seed"){
#       expr_line_df <- qread(file.path("/home/gugene/ShinyApps/CoPA/www/CoPA_files/expr_line/", "seed.qs"))
#     } else {
#       expr_line_df <- qread(file.path("/home/gugene/ShinyApps/CoPA/www/CoPA_files/expr_line/", "topmodposbc.qs"))
#     }

#     # Load appropriate bulk correlation data
#     bulk_cor_list <- qread(file = file.path("/home", "gugene", "ShinyApps", "CoPA", "www", "CoPA_files", "bulk_cor", "bulk_cor_list.qs"))

#     # Load appropriate projection data
#     proj_path_proj <- file.path("/home/gugene/ShinyApps/CoPA/www/CoPA_files/projs", paste0(def, "_", dstring, ".qs"))
#     projplotdf <- qread(proj_path_proj)
    
#     # Load appropriate projection correlation data
#     proj_path <- file.path("/home/gugene/ShinyApps/CoPA/www/CoPA_files/proj_cors", paste0(def, "_", dstring, ".qs"))
#     corlistdf <- qread(proj_path)

#     # Loop through modules

#     # expr_plots <- list()
#     # bulk_cor_plots <- list()
#     # gsea_plots <- list()
#     # projplotlist <- list()
#     # corlistind <- list()
#     for(i in seq_along(projplotdf)){

#       # Expression line plots
#       subt <- paste0(modulelengths_seed[these_mods[i]], " seed genes | ", modulelengths[these_mods[i]], " topmodposbc genes", "\n", sprintf("%.0f", vevec[these_mods[i]] * 100), "% variance explained by module eigengene")
#       if(length(expr_line_df[[i]]$gene) >= 10){
#         legt <- "Top 10 genes:"
#       } else {
#         legt <- "Genes:"
#       }
#       expr_plots <- ggplot(expr_line_df[[i]], aes(x=sample,y=value, color=gene )) + 
#         geom_line(linewidth = 0.2) +
#         labs(x = "Human prefrontal cortex samples\n(n = 100 random samples)",
#              y = "Expression z-score") + 
#         newtheme + 
#         theme(text = element_text(size = 16),
#               legend.box.margin = margin(0, 0, 0, -10),
#               plot.title = element_text(hjust = 0.5),
#               plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 10, 0)),
#               legend.title = element_text(size = 14),
#               legend.text = element_text(size = 14),
#               legend.position = "right") +
#         scale_color_manual(values = cols2) +
#         labs(#title = paste0("Module ", these_mods[i]),
#             title = paste0("Module ", i), # Numbering post-filter
#             subtitle = subt,
#             colour = legt
#             #subtitle = bquote(subt)
#             ) +
#         guides(color = guide_legend(override.aes = list(linewidth = 2)))

#       # bulk cors
#       bulk_cor_plots <- ggplot(bulk_cor_list[[i]], aes(x = split, y=cors, fill = type)) + 
#         theme_classic() +
#         geom_violin(position = position_dodge(width = 0.6)) +
#         geom_boxplot(notch=T,width=0.1, position = position_dodge(width = 0.6), show.legend = F, outlier.shape = NA) +
#         ggtitle("Distributions of bulk pairwise correlations") +
#         #geom_hline(color = "black", yintercept=0) +
#         ylab("Pearson correlation") +
#         xlab("") +
#         scale_x_discrete(labels = label_wrap(10)) +
#         # theme(plot.title = element_text(size=13,hjust=0.5),
#         #       plot.subtitle = element_text(size=13,hjust=0.5),
#         #       legend.position = "bottom",
#         #       legend.title=element_blank(),
#         #       axis.text.x = element_text(size=10),
#         #       axis.title.y = element_text(size=11),
#         #       axis.text.y = element_text(size=11),
#         #       axis.ticks.x=element_blank()) +
#         scale_fill_manual(values=c("salmon", "grey")) +
#         ylim(NA,1.1) +
#         ggpubr::stat_compare_means(aes(group = type), label = "p.signif",method.args = list(alternative = "less"), label.y=0.95) +
#         newtheme +
#         theme(legend.text = element_text(size = 8),
#               plot.title = element_text(hjust = 0.5, size = 14),
#               axis.text.x = element_text(size = 12),
#               legend.position = "bottom",
#               legend.direction = "vertical",
#               legend.title=element_blank(),
#               legend.margin = margin(-10, 0, 0, 0),
#               plot.margin = margin(1, 0.3, 0.6, 0.3, "cm")) +
#         guides(fill = guide_legend(override.aes = list(shape = 22)))
#       if(length(unique(bulk_plat_cors$split))==1){
#         bulk_cor_plots <- bulk_cor_plots +
#           theme(axis.text.x = element_blank(),
#                 axis.ticks.x=element_blank())
#       } else {
#         bulk_cor_plots <- bulk_cor_plots +
#           theme(axis.text.x = element_text(size=10))
#       }

#       # GSEA
#       ysizevar <- (max(nchar(as.character(gsea_plot_dfs[[i]]$SetName))) - 31)/58 * -4 + 10
#       ysizevar <- max(ysizevar, 6)
#       ysizevar <- min(ysizevar, 10)
      
#       gsea_plots <- gsea_plot_dfs[[i]] |>
#         dplyr::mutate(SetName_breaks = factor(SetName_breaks, levels=rev(unique(SetName_breaks)))) |>
#         ggplot(aes(x = SetName_breaks, y = pval)) +
#           newtheme + 
#           geom_bar(stat="identity") +
#           labs(x = "", y = bquote(-log[10]~"(p-val)"), title = "Geneset enrichment") +
#           theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5, size = 10, color = "black"),
#               axis.text.y = element_text(hjust = 1, vjust = 0.5, size = ysizevar, color = "black"),
#               axis.title.y = element_text(size = 11, color = "black"),
#               axis.title.x = element_text(size = 11, color = "black"),
#               plot.title = element_blank()) + 
#           coord_flip() +
#           #scale_x_discrete(limits = rev(levels(gsea_plot$SetName))) +
#           geom_hline(yintercept = gsea_cutFDR, color = "red")
      
#       # Projections
#       n <- length(unique(projplotdf[[i]]$ct))

#       projplotlist <- ggplot(projplotdf[[i]], aes(x = ct, y = ind, fill = ct)) +
#         theme_classic() +
#         geom_col(position = position_dodge(), alpha = 0.5) +
#         geom_errorbar(aes(ymin = ind - 2 * ind_se,
#                           ymax = ind + 2 * ind_se),
#                           width = 0.2,
#                           linewidth = 0.3,
#                           position = position_dodge(0.5)) +
#         scale_fill_manual(values = hue_pal()(n)) +
#         theme(text = element_text(family = "sans", color = "black", size = 14),
#               legend.position="none", 
#               axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
#               axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
#               strip.text = element_text(color = "black"),
#               strip.background = element_rect(fill = "white")) +
#         labs(y = index_xaxis_names[d], x = "") +       
#         facet_wrap(~dataset, ncol = 1, nrow = 5, scales = "free_y") 

#       # Projection cor heatmap
#       col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#       p <-  Heatmap(corlistdf[[i]], 
#                     column_title = "Correlations of module projections",
#                     column_title_gp = gpar(fontsize = 14),
#                     heatmap_legend_param = list(
#                       direction = "horizontal",
#                       title_gp = gpar(fontface = "plain"),
#                       title_position = "topcenter"
#                     ),
#                     name = "Correlation", 
#                     column_names_rot = 90,  
#                     col = col_fun,    
#                     row_names_gp = gpar(fontsize = 12),
#                     row_names_side = "left",
#                     row_dend_side = "right",
#                     width = unit(3, "cm"), 
#                     height = unit(3, "cm"),
#                     show_column_names = F)
#       corlistind <- as.ggplot(p)

#       # Arrange snapshots and save as SVGs
#       rightside <- plot_grid(projplotlist, corlistind, nrow = 2, 
#                              rel_heights = c(3, 1))
   
#       leftsideall <- suppressMessages(plot_grid(expr_plots,
#                                                 bulk_cor_plots,
#                                                 gsea_plots,
#                                                 nrow = 3, rel_heights = c(1,1,1)))

#       p <- suppressMessages(plot_grid(leftsideall, rightside, ncol=2))
#       ggsave(p, file = file.path(out_dir, paste0(i, ".svg")), device = svglite::svglite, width = 13, height = 9, bg = "white")

#     }


#     cat(def, index_dirs[d], "done\n")
#     timestamp()
#   }
# }

 

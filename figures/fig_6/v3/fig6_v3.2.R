# v3.2
# Change MIT paths to new outputs based on gabitto metacell mapping

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

# version <- "v3.2"
# save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/"), version)
# if(!dir.exists(save_dir))
#   dir.create(save_dir, recursive = F)

# Load bulk-related objects
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
filter_under <- 3

# Load bulk ctrl modules
datkme <- fread(data.table=F, file = file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]

sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods_final <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]

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
index_dirs <- c("log_native", "log_REI")
index_list <- list()
index_save_names <- c("native_log", "REI")
index_xaxis_names <- c("Mean expression (log UMI counts + 1)", "Relative expression index")
#for(d in 1:2){
d <- 1 # native_log

#### 0. Create various vectors and lists for all comparisons

# Helper functions for loading
make_name <- function(path, base, stem_prefix) {
  rel   <- sub(paste0(base, "/"), "", path)
  parts <- strsplit(rel, "/")[[1]]
  region <- parts[1]                                    # PFC or MTC
  stem   <- tools::file_path_sans_ext(basename(path))   # e.g. genomewide_means_allAD
  stem   <- sub(paste0("^", stem_prefix, "_"), "", stem) # strip leading prefix
  paste(region, stem, sep = "_")
}

make_se_name <- function(path, base) {
  rel    <- sub(paste0(base, "/"), "", path)
  parts  <- strsplit(rel, "/")[[1]]
  region    <- parts[1]                                   # PFC or MTC
  subfolder <- parts[3]                                   # log_native / log_normByMean / log_REI
  stem      <- tools::file_path_sans_ext(basename(path)) # se_{group}_{modtype}
  stem      <- sub("^se_", "", stem)                      # {group}_{modtype}
  paste(region, stem, subfolder, sep = "_")
}

make_dcopa_name <- function(path, base) {
  rel    <- sub(paste0(base, "/"), "", path)
  parts  <- strsplit(rel, "/")[[1]]
  region <- parts[1]                                     # PFC or MTC
  stem   <- tools::file_path_sans_ext(basename(path))    # e.g. allAD_vs_Con_bulk_megaset_output_table
  stem   <- sub("_output_table$", "", stem)              # allAD_vs_Con_bulk_megaset
  paste(region, stem, sep = "_")
}

### Load all SEA data
sea_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output")
sea_files <- list.files(sea_dir,
                        pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
# Exclude seed-module outputs (e.g. *_bulk_megaset_seed.csv): shiny-app only,
# would otherwise be swept up here and double-match the grepl(keys) selection.
sea_files <- sea_files[!grepl("_seed", sea_files)]
# SEA genomewide means
sea_means_files <- sea_files[grepl("/means/", sea_files)]
sea_means_list <- setNames(
  lapply(sea_means_files, \(x) fread(x, data.table = F)),
  sapply(sea_means_files, make_name, base = sea_dir, stem_prefix = "genomewide_means")
)

# SEA log native projections
sea_mod_means_files <- sea_files[grepl("/mod_means/log_native/", sea_files)]
sea_mod_means_list <- setNames(
  lapply(sea_mod_means_files, \(x) fread(x, data.table = F)),
  sapply(sea_mod_means_files, make_name, base = sea_dir, stem_prefix = "mod_means")
)

# SEA log native projections (SE)
sea_se_files <- sea_files[grepl("/se/log_native/", sea_files)]
sea_se_list <- setNames(
  lapply(sea_se_files, \(x) fread(x, data.table = F)),
  sapply(sea_se_files, make_se_name, base = sea_dir)
)

# SEA dcopa output
sea_dcopa_files <- sea_files[grepl("/euclidean_distances/", sea_files)]
sea_dcopa_list <- setNames(
  lapply(sea_dcopa_files, \(x) fread(x, data.table = F)),
  sapply(sea_dcopa_files, make_dcopa_name, base = sea_dir)
)

### Load all MIT data
mit_dir <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output"

mit_files <- list.files(mit_dir,
                        pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
# Exclude seed-module outputs (e.g. *_bulk_megaset_seed.csv): they share this
# tree but are only for the shiny app, and would otherwise be swept up here and
# double-match the grepl(keys) selection below.
mit_files <- mit_files[!grepl("_seed", mit_files)]

# MIT genomewide means
mit_means_files <- mit_files[grepl("/means/", mit_files)]
mit_means_list <- setNames(
  lapply(mit_means_files, \(x) fread(x, data.table = F)),
  sapply(mit_means_files, make_name, base = mit_dir, stem_prefix = "genomewide_means")
)

# MIT log native projections
mit_mod_means_files <- mit_files[grepl("/mod_means/log_native/", mit_files)]
mit_mod_means_list <- setNames(
  lapply(mit_mod_means_files, \(x) fread(x, data.table = F)),
  sapply(mit_mod_means_files, make_name, base = mit_dir, stem_prefix = "mod_means")
)

# MIT log native projections (SE)
mit_se_files <- mit_files[grepl("/se/log_native/", mit_files)]
mit_se_list <- setNames(
  lapply(mit_se_files, \(x) fread(x, data.table = F)),
  sapply(mit_se_files, make_se_name, base = mit_dir)
)

# MIT dcopa output
mit_dcopa_files <- mit_files[grepl("/euclidean_distances/", mit_files)]
mit_dcopa_list <- setNames(
  lapply(mit_dcopa_files, \(x) fread(x, data.table = F)),
  sapply(mit_dcopa_files, make_dcopa_name, base = mit_dir)
)

keys <- list(
  c("PFC", "allAD", "Con", "bulk_megaset"),
  c("PFC", "earlyAD", "Con", "bulk_megaset"),
  c("PFC", "lateAD", "earlyAD", "bulk_megaset"),
  c("PFC", "APOE44", "APOE33", "bulk_megaset"),
  c("MTC", "allAD", "Con", "bulk_megaset"),
  c("MTC", "earlyAD", "Con", "bulk_megaset"),
  c("MTC", "lateAD", "earlyAD", "bulk_megaset"),
  c("MTC", "APOE44", "APOE33", "bulk_megaset"),
  c("PFC", "allAD", "Con", "rosmap"),
  c("PFC", "earlyAD", "Con", "rosmap"),
  c("PFC", "lateAD", "earlyAD", "rosmap"),
  c("PFC", "APOE44", "APOE33", "rosmap"),
  c("MTC", "allAD", "Con", "rosmap"),
  c("MTC", "earlyAD", "Con", "rosmap"),
  c("MTC", "lateAD", "earlyAD", "rosmap"),
  c("MTC", "APOE44", "APOE33", "rosmap")
)

title_vec <- c("CTRL vs AD (DFC)",
               "Early AD vs CTRL (DFC)",
               "Late vs Early AD (DFC)",
               "APOE 4/4 vs 3/3 (DFC)",
               "CTRL vs AD (MTG)",
               "Early AD vs CTRL (MTG)",
               "Late vs Early AD (MTG)",
               "APOE 4/4 vs 3/3 (MTG)",
               "CTRL vs AD (ROSMAP, DFC)",
               "Early AD vs CTRL (ROSMAP, DFC)",
               "Late vs Early AD (ROSMAP, DFC)",
               "APOE 4/4 vs 3/3 (ROSMAP, DFC)",
               "CTRL vs AD (ROSMAP, MTG)",
               "Early AD vs CTRL (ROSMAP, MTG)",
               "Late vs Early AD (ROSMAP, MTG)",
            #  "SCZ vs CTRL (DFC)",
               "APOE 4/4 vs 3/3 (ROSMAP, MTG)")

title_vec_dataset <- list(
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
#  c("CMC", "SZBDMulti-Seq"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025"),
  c("Gabitto 2024", "Liu 2025")
)

dx_list <- list(
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
  c("APOE 4/4", "APOE 3/3"),
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
  c("APOE 4/4", "APOE 3/3"),
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
  c("APOE 4/4", "APOE 3/3"),
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
 # c("SCZ", "Con"),
  c("APOE 4/4", "APOE 3/3")
)

# Save file suffix
save_suffix_vec <- c("AllADVsCon_DFC",
                     "EarlyADVsCon_DFC",
                     "LateVsEarlyAD_DFC",
                     "APOE44_vs_33_DFC",
                     "AllADVsCon_MTG",
                     "EarlyADVsCon_MTG",
                     "LateVsEarlyAD_MTG",
                     "APOE44_vs_33_MTG",
                     "AllADVsCon_DFC_ROSMAP",
                     "EarlyADVsCon_DFC_ROSMAP",
                     "LateVsEarlyAD_DFC_ROSMAP",
                     "APOE44_vs_33_DFC_ROSMAP",
                     "AllADVsCon_MTG_ROSMAP",
                     "EarlyADVsCon_MTG_ROSMAP",
                     "LateVsEarlyAD_MTG_ROSMAP",
                     #"CMC_vs_SZBD",
                     "APOE44_vs_33_MTG_ROSMAP"
                     )


# failed at i = 12

##### Run loop for all comparisons
#for(i in 1:length(keys)){
for(i in 1:8){ # Need to create separate objects for ROSMAP (new plot objects, new input set of mods)

  # Load genomewide means by ct
  mean_sum_gab <- list(
    sea_means_list[[which(grepl(keys[[i]][1], names(sea_means_list)) & grepl(keys[[i]][2], names(sea_means_list)))]],
    sea_means_list[[which(grepl(keys[[i]][1], names(sea_means_list)) & grepl(keys[[i]][3], names(sea_means_list)))]]
  ) |> 
    lapply(\(x) x |> tibble::column_to_rownames("gene_ids") |> as.data.frame()) |>
    setNames(dx_list[[i]])

  mean_sum_mit <- list(
    mit_means_list[[which(grepl(keys[[i]][1], names(mit_means_list)) & grepl(keys[[i]][2], names(mit_means_list)))]],
    mit_means_list[[which(grepl(keys[[i]][1], names(mit_means_list)) & grepl(keys[[i]][3], names(mit_means_list)))]]
  ) |> 
    lapply(\(x) x |> tibble::column_to_rownames("gene_sym") |> as.data.frame()) |>
    setNames(dx_list[[i]])


  # Load dcopa data
  dcopa_gab <- sea_dcopa_list[[which(grepl(keys[[i]][1], names(sea_dcopa_list)) & grepl(keys[[i]][2], names(sea_dcopa_list)) & grepl(keys[[i]][3], names(sea_dcopa_list)) & grepl(keys[[i]][4], names(sea_dcopa_list)))]]
  dcopa_mit <- mit_dcopa_list[[which(grepl(keys[[i]][1], names(mit_dcopa_list)) & grepl(keys[[i]][2], names(mit_dcopa_list)) & grepl(keys[[i]][3], names(mit_dcopa_list)) & grepl(keys[[i]][4], names(mit_dcopa_list)))]]
  dcopa_shared <- dplyr::inner_join(dcopa_gab, dcopa_mit, by = dplyr::join_by(mod, Direction, Consistency), suffix = c("1", "2"), relationship = "many-to-many") |>
      dplyr::filter(Celltype1 == Celltype2,
                    sig_FDR1,
                    sig_FDR2,
                    Consistency %in% c(0, 1)) |>
      dplyr::arrange(mod)

  ## Load indices
  # Gabitto AD indices
  gabad <- sea_mod_means_list[[which(grepl(keys[[i]][1], names(sea_mod_means_list)) & grepl(keys[[i]][2], names(sea_mod_means_list)) & grepl(keys[[i]][4], names(sea_mod_means_list)))]]
  gabad_se <- sea_se_list[[which(grepl(keys[[i]][1], names(sea_se_list)) & grepl(keys[[i]][2], names(sea_se_list)) & grepl(keys[[i]][4], names(sea_se_list)))]]
 
  # Gabitto control indices
  gabcon <- sea_mod_means_list[[which(grepl(keys[[i]][1], names(sea_mod_means_list)) & grepl(keys[[i]][3], names(sea_mod_means_list)) & grepl(keys[[i]][4], names(sea_mod_means_list)))]]
  gabcon_se <- sea_se_list[[which(grepl(keys[[i]][1], names(sea_se_list)) & grepl(keys[[i]][3], names(sea_se_list)) & grepl(keys[[i]][4], names(sea_se_list)))]]
 
  # Liu AD indices
  mitad <- mit_mod_means_list[[which(grepl(keys[[i]][1], names(mit_mod_means_list)) & grepl(keys[[i]][2], names(mit_mod_means_list)) & grepl(keys[[i]][4], names(mit_mod_means_list)))]]
  mitad_se <- mit_se_list[[which(grepl(keys[[i]][1], names(mit_se_list)) & grepl(keys[[i]][2], names(mit_se_list)) & grepl(keys[[i]][4], names(mit_se_list)))]]
 
  # Liu control indices
  mitcon <- mit_mod_means_list[[which(grepl(keys[[i]][1], names(mit_mod_means_list)) & grepl(keys[[i]][3], names(mit_mod_means_list)) & grepl(keys[[i]][4], names(mit_mod_means_list)))]]
  mitcon_se <- mit_se_list[[which(grepl(keys[[i]][1], names(mit_se_list)) & grepl(keys[[i]][3], names(mit_se_list)) & grepl(keys[[i]][4], names(mit_se_list)))]]
 
  # Find intersection of all cts over all datasets
  allcts <- intersect(intersect(colnames(gabcon), colnames(gabad)), intersect(colnames(mitcon), colnames(mitad)))

  # Align columns
  gabad     <- gabad     |> select(all_of(allcts))
  gabad_se  <- gabad_se  |> select(all_of(allcts))
  gabcon    <- gabcon    |> select(all_of(allcts))
  gabcon_se <- gabcon_se |> select(all_of(allcts))
  mitad     <- mitad     |> select(all_of(allcts))
  mitad_se  <- mitad_se  |> select(all_of(allcts))
  mitcon    <- mitcon    |> select(all_of(allcts))
  mitcon_se <- mitcon_se |> select(all_of(allcts))

  # Capitalize celltypes
  CT_RENAME <- c(
    "Lamp5"     = "LAMP5",
    "Lamp5 Lhx6"= "LAMP5 LHX6",
    "Pax6"      = "PAX6",
    "Pvalb"     = "PVALB",
    "Sncg"      = "SNCG",
    "Sst"       = "SST",
    "Vip"       = "VIP",
    "Sst Chodl" = "SST CHODL"
  )

  allcts_cap <- allcts |> dplyr::recode(!!!CT_RENAME)

  # Capitalize index titles
  colnames(gabad)     <- allcts_cap
  colnames(gabad_se)  <- allcts_cap
  colnames(gabcon)    <- allcts_cap
  colnames(gabcon_se) <- allcts_cap
  colnames(mitad)     <- allcts_cap
  colnames(mitad_se)  <- allcts_cap
  colnames(mitcon)    <- allcts_cap
  colnames(mitcon_se) <- allcts_cap

  # Capitalize dcopa_shared celltype columns
  dcopa_shared <- dcopa_shared |> 
    mutate(Celltype1 = dplyr::recode(Celltype1, !!!CT_RENAME),
           Celltype2 = dplyr::recode(Celltype2, !!!CT_RENAME))
  dcopa_gab <- dcopa_gab |> 
    mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))
  dcopa_mit <- dcopa_mit |> 
    mutate(Celltype = dplyr::recode(Celltype, !!!CT_RENAME))

  # Capitalize mean_sum_* column names
  mean_sum_gab <- lapply(mean_sum_gab, \(x){
    colnames(x) <- dplyr::recode(colnames(x), !!!CT_RENAME)
    return(x)
  })

  mean_sum_mit <- lapply(mean_sum_mit, \(x){
    colnames(x) <- dplyr::recode(colnames(x), !!!CT_RENAME)
    return(x)
  })

  # Plot projection indices
  projplotlist1 <- lapply(unique(dcopa_shared$mod), \(j){ # Gabitto

    # Fetch projection index data for the module
    plotdf1 <- data.frame("ct" = rep(allcts_cap, 2),
                          "dataset" = c(rep(paste0(title_vec_dataset[[i]][1], ", ", title_vec[i]), 2 * length(allcts))),
                          "ind" = c(unlist(gabcon[j, ]), 
                                    unlist(gabad[j, ])),
                          "ind_se" = c(unlist(gabcon_se[j, ]), 
                                       unlist(gabad_se[j, ])),
                          "dx" = c(rep(dx_list[[i]][2], length(allcts)),
                                   rep(dx_list[[i]][1], length(allcts)))) |>
      dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                    ct = factor(ct, levels = allcts_cap))

    # Fetch dCoPA significance data for given module
    highlight_mat1 <- dcopa_shared |> 
      filter(mod == j, sig_FDR1 == T) |>
      select(Celltype1, P_value1) |>
      filter(!duplicated(Celltype1))
    
    # Plot grouped barplot of projection indices with asterisks for dCoPA significance
    plotdf1 <- plotdf1 |>
      left_join(highlight_mat1, by = join_by(ct == Celltype1)) |>
      mutate(
        label = case_when(
          P_value1 > 0.05 ~ "",
          (0.01 < P_value1 & P_value1 <= 0.05) ~ "*",
          (0.001 < P_value1 & P_value1 <= 0.01) ~ "**",
          (0.0001 < P_value1 & P_value1 <= 0.001) ~ "***",
          P_value1 <= 0.0001 ~ "****",
          .default = NA
        )) |>
      select(!P_value1)

    plotdf2 <- data.frame("ct" = rep(allcts_cap, 2),
                          "dataset" = c(rep(paste0(title_vec_dataset[[i]][2], ", ", title_vec[i]), 2 * length(allcts))),
                          "ind" = c(unlist(mitcon[j, ]),
                                    unlist(mitad[j, ])),
                          "ind_se" = c(unlist(mitcon_se[j, ]),
                                      unlist(mitad_se[j, ])),
                          "dx" = c(rep(dx_list[[i]][2], length(allcts)),
                                   rep(dx_list[[i]][1], length(allcts)))) |>
      dplyr::mutate(dataset = factor(dataset, levels = unique(dataset)),
                    ct = factor(ct, levels = allcts_cap))
    
    # Fetch dCoPA significance data for given module
    highlight_mat2 <- dcopa_shared |> 
      filter(mod == j, sig_FDR2 == T) |>
      select(Celltype2, P_value2) |>
      filter(!duplicated(Celltype2))

    plotdf2 <- plotdf2 |>
        left_join(highlight_mat2, by = join_by(ct == Celltype2)) |>
        mutate(
          label = case_when(
            P_value2 > 0.05 ~ "",
            (0.01 < P_value2 & P_value2 <= 0.05) ~ "*",
            (0.001 < P_value2 & P_value2 <= 0.01) ~ "**",
            (0.0001 < P_value2 & P_value2 <= 0.001) ~ "***",
            P_value2 <= 0.0001 ~ "****",
            .default = NA
          )) |>
      select(!P_value2)

    plotdf <- rbind(plotdf1, plotdf2) |>
      mutate(dx = factor(dx, levels = rev(dx_list[[i]])))

    # Remove duplicate labels (assign label to whichever dx has highest value)
    plotdf_label <- plotdf |> 
      filter(!is.na(label)) |>
      group_by(ct, dataset) |>
      slice_max(ind, n = 1)

    dodge_width <- 0.6
    plot_max_y <- plotdf$ind[which.max(plotdf$ind)] + plotdf$ind_se[which.max(plotdf$ind)] * 2
    p <- ggplot(plotdf, aes(x = ct, y = ind, fill = dx)) +
      theme_classic() +
      geom_col(position = position_dodge(width = dodge_width), alpha = 0.5) +
      # geom_text(data = plotdf_label, 
      #           aes(x = ct, y = ind, label = label), 
      #           vjust = -0.5#,
      #           #y = plot_max_y + 0.1
      #           ) + 
      geom_errorbar(aes(ymin = ind - 2 * ind_se,
                        ymax = ind + 2 * ind_se),
                        width = 0.2,
                        linewidth = 0.3,
                        position = position_dodge(width = dodge_width)) +
      scale_fill_manual(values = c("#90D5FF", "#FFA500")) +
      theme(text = element_text(family = "sans", color = "black", size = 12),
            legend.position = "bottom", 
            legend.title = element_blank(),
            legend.margin = margin(-0.8, 0, 0, 0, "cm"),
            axis.title.y = element_text(size = 10, margin = margin(0, 5, 0, 0)),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
            strip.text = element_text(color = "black"),
            strip.background = element_rect(fill = "white"),
            plot.margin = margin(5, 0, 0, 0)) +
      labs(y = index_xaxis_names[d], x = "") +
      scale_y_continuous(#breaks = c(0, 0.5, 1)#, 
                         limits = c(0, plot_max_y + 0.27)
                        ) +
      facet_wrap(~dataset, nrow = 2, ncol = 1, scales = "free_y") + 
      geom_text(      
        data = plotdf_label,                                 
        aes(x = ct, y = ind + 2 * ind_se + 0.15, label =     
      "↓"),                                                  
        position = position_dodge(width = dodge_width),      
        size = 3.5,                                          
        vjust = 0                                            
      ) #+
      # geom_segment(                                      
      #   data = plotdf_label,
      #   aes(x = ct, xend = ct,                               
      #       y = ind + 2 * ind_se + 0.06,   # start above error bar                                              
      #       yend = ind + 2 * ind_se + 0.02), # end closer to bar                                                    
      #   position = position_dodge(width = dodge_width),
      #   arrow = arrow(length = unit(0.15, "cm"), type = "closed"),                                             
      #   linewidth = 0.4,                                     
      #   color = "black"                                      
      # )           
      return(p)
  })

  ## Plot boxplots comparing individual module gene expression changes
  boxplotlist_gab <- lapply(unique(dcopa_shared$mod), \(j){
    cat(j, " ")
    highlight_mat <- dcopa_shared |> 
      filter(mod == j, sig_FDR2 == T) |>
      select(Celltype1, P_value1)
    summary_mat <- mean_sum_gab

    if(length(unique(highlight_mat[,1]) > 2)){
      stripsize <- 8
    } else {
      stripsize <- 10
    }
   
    out_df <- mapply(\(x, y){
      df <- x[rownames(x) %in% mods[[which(these_mods == j)]], colnames(x) %in% highlight_mat[,1], drop = F] 
      out <- data.frame("type" = y, "gene" = rownames(df), df, check.names = F) |>
        pivot_longer(!type:gene, names_to = "ct", values_to = "vals")
      return(out)
    }, summary_mat[1:2], names(summary_mat)[1:2], SIMPLIFY = F) |>
      do.call(what = "rbind") |>
      mutate(type = factor(type, levels = rev(dx_list[[i]]))) |>
      left_join(highlight_mat, by = join_by(ct == Celltype1)) |>
      mutate(
        label = case_when(
          P_value1 > 0.05 ~ "",
          (0.01 < P_value1 & P_value1 <= 0.05) ~ "*",
          (0.001 < P_value1 & P_value1 <= 0.01) ~ "**",
          (0.0001 < P_value1 & P_value1 <= 0.001) ~ "***",
          P_value1 <= 0.0001 ~ "****",
          .default = NA
        )) |>
      select(!P_value1)

    outdf_label <- out_df |> 
      filter(!is.na(label)) |>
      group_by(ct) |>
      slice_max(vals, n = 1)

    outdf_label <- out_df |>                               
      group_by(ct) |>                                      
      summarise(                                           
        vals = max(vals, na.rm = TRUE),
        label = unique(label[!is.na(label) & label != ""])[1]                                                
      ) |>
      filter(!is.na(label))  
    
    p <- ggplot(out_df, aes(x = type, y = vals)) +
            theme_classic() + 
            geom_point(color = "lightgrey", alpha = 0.7, size = 0.4) + 
            geom_line(aes(group = gene), color = "lightgrey", alpha = 0.7, linewidth = 0.2) +
            geom_boxplot(aes(color = type), linewidth = 0.4, notch = F, outlier.size = 0.4, alpha = 0.7, fill = NA) + 
            geom_text(      
              data = outdf_label,
              aes(y = vals + 0.1, label = label),                                            
              x = 1.5,
              size = 3,                                            
              vjust = 0                                            
            ) +
            facet_wrap(~ct, nrow = 1) + 
            labs(x = "", y = "Mean expression\n(log UMI counts + 1)", title = "Gabitto et al. 2024") +
            scale_color_manual(values = c("#90D5FF", "#FFA500")) +
            scale_y_continuous(limits = c(NA, max(out_df$vals, na.rm = TRUE) + 0.4)) +
            theme(legend.position = "none", 
                  plot.title = element_text(size = 12, hjust = 0.5),
                  axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 1),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.text = element_text(size = stripsize)
                  ) 
    return(p)
  })

  boxplotlist_liu <- lapply(unique(dcopa_shared$mod), \(j){
    highlight_mat <- dcopa_shared |> 
      filter(mod == j, sig_FDR2 == T) |>
      select(Celltype2, P_value2)
    summary_mat <- mean_sum_mit
    
    if(length(unique(highlight_mat[,1]) > 2)){
      stripsize <- 8
    } else {
      stripsize <- 10
    }

    out_df <- mapply(\(x, y){
      df <- x[rownames(x) %in% mods[[which(these_mods == j)]], colnames(x) %in% highlight_mat[,1], drop = F] 
      out <- data.frame("type" = y, "gene" = rownames(df), df, check.names = F) |>
        pivot_longer(!type:gene, names_to = "ct", values_to = "vals")
      return(out)
    }, summary_mat[1:2], names(summary_mat)[1:2], SIMPLIFY = F) |>
      do.call(what = "rbind") |>
      mutate(type = factor(type, levels = rev(dx_list[[i]]))) |>
      left_join(highlight_mat, by = join_by(ct == Celltype2)) |>
      mutate(
        label = case_when(
          P_value2 > 0.05 ~ "",
          (0.01 < P_value2 & P_value2 <= 0.05) ~ "*",
          (0.001 < P_value2 & P_value2 <= 0.01) ~ "**",
          (0.0001 < P_value2 & P_value2 <= 0.001) ~ "***",
          P_value2 <= 0.0001 ~ "****",
          .default = NA
        )) |>
      select(!P_value2)

    outdf_label <- out_df |> 
      filter(!is.na(label)) |>
      group_by(ct) |>
      slice_max(vals, n = 1)

    outdf_label <- out_df |>                               
      group_by(ct) |>                                      
      summarise(                                           
        vals = max(vals, na.rm = TRUE),
        label = unique(label[!is.na(label) & label != ""])[1]                                                
      ) |>
      filter(!is.na(label))  
          
    p <- ggplot(out_df, aes(x = type, y = vals)) +
      theme_classic() + 
      geom_point(color = "lightgrey", alpha = 0.7, size = 0.4) + 
      geom_line(aes(group = gene), color = "lightgrey", alpha = 0.7, linewidth = 0.2) +
      geom_boxplot(aes(color = type), linewidth = 0.4, notch = F, outlier.size = 0.4, alpha = 0.7, fill = NA) + 
      geom_text(      
        data = outdf_label,
        aes(y = vals + 0.1, label = label),                                            
        x = 1.5,
        size = 3,                                            
        vjust = 0                                            
      ) +
      facet_wrap(~ct, nrow = 1) + 
      labs(x = "", y = "Mean expression\n(log UMI counts + 1)", title = "Liu et al. 2025") +
      scale_color_manual(values = c("#90D5FF", "#FFA500")) +
      scale_y_continuous(limits = c(NA, max(out_df$vals, na.rm = TRUE) + 0.4)) +
      theme(legend.position = "none", 
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = stripsize)
            ) 

    return(p)
  })

  # Construct and save module snapshot

  out_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/dCoPA_snapshots/"), save_suffix_vec[i])
  if(!dir.exists(out_dir))
    dir.create(out_dir, recursive = T)

  for(j in seq_along(unique(dcopa_shared$mod))){
    boxplot_combined <- plot_grid(boxplotlist_gab[[j]], boxplotlist_liu[[j]], nrow = 1)
    pall <- suppressMessages(plot_grid(all_plots[[1]][[which(these_mods == unique(dcopa_shared$mod)[j])]],
                                      all_plots[[2]][[which(these_mods == unique(dcopa_shared$mod)[j])]],
                                      projplotlist1[[j]],
                                      NULL,
                                      boxplot_combined,
                                      nrow = 5, 
                                      rel_heights = c(0.5, 0.3, 0.8, 0.05, 0.5)))
    ggsave(pall, 
          file = file.path(out_dir, paste0(unique(dcopa_shared$mod)[j], ".svg")), 
          device = svglite::svglite, 
          width = 6, 
          height = 12, 
          bg = "white")
    cat(j, " ")
  }
  cat(i, "\n")
}


######
# Mods to feature in Fig 6
#####

# Genes implicated in L4IT (see fig 7, L4IT_MTG_AD_gene_analysis.docx, ranked in order of relevance)
genetest <- c("MAPT", "ICA1", "DLG2", "ANKH", "HTT", "TANC2", "MAPK8", "MSH2", "TOLLIP", "NUS1", "ATP6V0A1", "PLBD2")
# Find mods that contain these genes
modhits <- lapply(mods, \(x) any(x %in% genetest)) |> unlist() |> which()
sum(these_mods[modhits] %in% dcopa_shared[[1]][,1])
# 7
# Find which mod has the most hits
modhitcount <- lapply(mods, \(x) sum(x %in% genetest)) |> unlist()
modhitcount[which(these_mods %in% dcopa_shared[[1]][,1])]
# Mods 54, 99, 337 (1023 indices) have more than two hits; others have 1 or 0
# 1158 indices: 59, 106, 358
# 1016 indices: 54, 99, 337
# Which mod has MAPT?
lapply(mods, \(x) any(x %in% "MAPT")) |> unlist() |> which()
# 82 (1023 indices)

# Use:
# 82 (mod that contains MAPT, top hit in L4IT_MTG_AD_gene_analysis.docx)
#  - 12 topmodposbc genes total, MAPT is top 12, maybe include all 12 genes
# 99 (3rd and 4th ranked hits)  
# 54 (7th and 8th hits both in top 10 genes)

# Run loop for specific modules:
plot_these_mods <- these_mods[c(54, 99)] # need to redo the top panel for 82 if i want to use it


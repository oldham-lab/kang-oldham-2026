# v4
# generalize for ad/scz modules

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
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/gsea_func_optimized.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/fig6_fxns.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/full_pipeline.R"))

# expr line plots + GSEA depend only on module gene sets + bulk expression, so they are
# INVARIANT to the HGNC projection fix. Cache them per (module set, bulk) to skip the
# expensive recompute on re-runs. Delete gsea_expr_cache/ if module defs or bulk change.
GSEA_CACHE_DIR <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/v4/gsea_expr_cache")
dir.create(GSEA_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
cached_all_plots <- function(mods, bulk_expr, these_mods, datkme) {
  key <- paste(length(mods), sum(lengths(mods)), nrow(bulk_expr), ncol(bulk_expr),
               sum(these_mods), sep = "_")
  cp  <- file.path(GSEA_CACHE_DIR, paste0("allplots_", key, ".qs"))
  if (file.exists(cp)) { message("[cache] load ", basename(cp)); return(qread(cp)) }
  message("[cache] compute all_plots (", length(mods), " modules)")
  ap <- list()
  ap[[1]] <- make_expr_line_plots(expr = bulk_expr, mods = mods, these_mods = these_mods, datkme = datkme)
  gsea_broad <- run_gsea_for_proj_optimized(mods, broad = TRUE)
  gsea       <- run_gsea_for_proj_optimized(mods, broad = FALSE)
  ap[[2]] <- make_gsea_plots(gsea = gsea, gsea_broad = gsea_broad, these_mods = these_mods)
  qsave(ap, cp); ap
}

# version <- "v4"
# save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/"), version)
# if(!dir.exists(save_dir))
#   dir.create(save_dir, recursive = F)

################### AD (CTRL modules)
# Load bulk expression matrix
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)

# Load modules
filter_under <- 3
datkme <- fread(data.table=F, file = file.path(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"),"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]

sea_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output")
mit_dir <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output")

snapshot_out_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/dCoPA_snapshots/")

keys <- list(
  c("PFC", "allAD", "Con", "bulk_megaset"),
  c("PFC", "earlyAD", "Con", "bulk_megaset"),
  c("PFC", "lateAD", "earlyAD", "bulk_megaset"),
  c("PFC", "APOE44", "APOE33", "bulk_megaset"),
  c("MTC", "allAD", "Con", "bulk_megaset"),
  c("MTC", "earlyAD", "Con", "bulk_megaset"),
  c("MTC", "lateAD", "earlyAD", "bulk_megaset"),
  c("MTC", "APOE44", "APOE33", "bulk_megaset")
)

title_vec <- c("CTRL vs AD (DFC)",
               "Early AD vs CTRL (DFC)",
               "Late vs Early AD (DFC)",
               "APOE 4/4 vs 3/3 (DFC)",
               "CTRL vs AD (MTG)",
               "Early AD vs CTRL (MTG)",
               "Late vs Early AD (MTG)",
               "APOE 4/4 vs 3/3 (MTG)")

title_vec_dataset <- c("Gabitto et al. 2024", "Liu et al. 2025")

dx_list <- list(
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
  c("APOE 4/4", "APOE 3/3"),
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
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
                     "APOE44_vs_33_MTG"
                     )

################################

# Create expr + GSEA objects (cached; HGNC-invariant, keyed by module set + bulk)
all_plots <- cached_all_plots(mods = mods, bulk_expr = bulk_expr,
                              these_mods = these_mods, datkme = datkme)
 
fig6_full_pipeline(all_plots = all_plots,
                   mods = mods,
                   these_mods = these_mods,
                   sea_dir = sea_dir,
                   mit_dir = mit_dir,
                   snapshot_out_dir = snapshot_out_dir,
                   keys = keys,
                   title_vec = title_vec,
                   title_vec_dataset = title_vec_dataset,
                   dx_list = dx_list,
                   save_suffix_vec = save_suffix_vec)

# OPTIMIZATION: cache block-1 expr+GSEA (SEAAD modules) for reuse in the SCZ CTRL-mods block below.
all_plots_seaad <- all_plots

################### AD (ROSMAP)
# Load bulk expression matrix
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_AD_only_SampleNetworks/1_10-52-54/ROSMAP_samp_filt_AD_only_1_248_ComBat.csv"), data.table=F)
megaset_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulk_expr <- bulk_expr[match(megaset_expr[,1],bulk_expr[,1]),]
bulk_expr <- data.frame("ensembl_id" = bulk_expr[,1],"Gene"= megaset_expr[,2], bulk_expr[,3:ncol(bulk_expr)])
bulk_expr <- bulk_expr[!is.na(bulk_expr[,1]),]

# Load modules
filter_under <- 3
datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]

sea_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output")
mit_dir <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output")

snapshot_out_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/dCoPA_snapshots/")

keys <- list(
  c("PFC", "allAD", "Con", "rosmap"),
  c("PFC", "earlyAD", "Con", "rosmap"),
  c("PFC", "lateAD", "earlyAD", "rosmap"),
  c("PFC", "APOE44", "APOE33", "rosmap"),
  c("MTC", "allAD", "Con", "rosmap"),
  c("MTC", "earlyAD", "Con", "rosmap"),
  c("MTC", "lateAD", "earlyAD", "rosmap"),
  c("MTC", "APOE44", "APOE33", "rosmap")
)

title_vec <- c("CTRL vs AD (ROSMAP, DFC)",
               "Early AD vs CTRL (ROSMAP, DFC)",
               "Late vs Early AD (ROSMAP, DFC)",
               "APOE 4/4 vs 3/3 (ROSMAP, DFC)",
               "CTRL vs AD (ROSMAP, MTG)",
               "Early AD vs CTRL (ROSMAP, MTG)",
               "Late vs Early AD (ROSMAP, MTG)",
               "APOE 4/4 vs 3/3 (ROSMAP, MTG)")

title_vec_dataset <- c("Gabitto et al. 2024", "Liu et al. 2025")

dx_list <- list(
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
  c("APOE 4/4", "APOE 3/3"),
  c("AD", "CTRL"),
  c("Early AD", "CTRL"),
  c("Late AD", "Early AD"),
  c("APOE 4/4", "APOE 3/3")
)

# Save file suffix
save_suffix_vec <- c("AllADVsCon_DFC_ROSMAP",
                     "EarlyADVsCon_DFC_ROSMAP",
                     "LateVsEarlyAD_DFC_ROSMAP",
                     "APOE44_vs_33_DFC_ROSMAP",
                     "AllADVsCon_MTG_ROSMAP",
                     "EarlyADVsCon_MTG_ROSMAP",
                     "LateVsEarlyAD_MTG_ROSMAP",
                     "APOE44_vs_33_MTG_ROSMAP"
                     )

################################

# Create expr + GSEA objects (cached; HGNC-invariant, keyed by module set + bulk)
all_plots <- cached_all_plots(mods = mods, bulk_expr = bulk_expr,
                              these_mods = these_mods, datkme = datkme)

fig6_full_pipeline(all_plots = all_plots,
                   mods = mods,
                   these_mods = these_mods,
                   sea_dir = sea_dir,
                   mit_dir = mit_dir,
                   snapshot_out_dir = snapshot_out_dir,
                   keys = keys,
                   title_vec = title_vec,
                   title_vec_dataset = title_vec_dataset,
                   dx_list = dx_list,
                   save_suffix_vec = save_suffix_vec)

################## SCZ (CTRL mods)
# Load bulk expression matrix
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)

# Load modules
filter_under <- 3
datkme <- fread(data.table=F, file = file.path(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"), "kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]

sea_dir <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC")
mit_dir <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq")

snapshot_out_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/dCoPA_snapshots/")

keys <- list(
  c("DFC", "Schizophrenia", "control", "bulk_megaset")
)

title_vec <- c("CTRL vs SCZ (DFC)")

title_vec_dataset <- c("CMC", "SZBDMulti-Seq")

dx_list <- list(
  c("SCZ", "CTRL")
)

# Save file suffix
save_suffix_vec <- c("SCZvsCon_DFC_bulkmegaset")
################################

# OPTIMIZATION: reuse block-1 expr+GSEA plots. The SCZ CTRL-mods block uses the SAME
# SEAAD modules + combined_FCX bulk as the AD bulk_megaset block, so make_expr_line_plots
# + GSEA are identical -> skip recomputing (avoids a full run_gsea_for_proj_optimized pass).
all_plots <- all_plots_seaad

fig6_full_pipeline(all_plots = all_plots,
                   mods = mods,
                   these_mods = these_mods,
                   sea_dir = sea_dir,
                   mit_dir = mit_dir,
                   snapshot_out_dir = snapshot_out_dir,
                   keys = keys,
                   title_vec = title_vec,
                   title_vec_dataset = title_vec_dataset,
                   dx_list = dx_list,
                   save_suffix_vec = save_suffix_vec)

################## SCZ (SCZ mods)
# Load bulk expression matrix
bulk_expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SCZ_SampleNetworks/1_04-04-29/brainseq_samp_filt_SCZ_1_171_outliers_removed_geneSymbolsAdded.csv"), data.table=F)

# Load modules
filter_under <- 3
datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]

sea_dir <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC")
mit_dir <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq")

snapshot_out_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/dCoPA_snapshots/")

keys <- list(
  c("DFC", "Schizophrenia", "control", "brainseq_scz")
)

title_vec <- c("CTRL vs SCZ (Brainseq DFC)")

title_vec_dataset <- c("CMC", "SZBDMulti-Seq")

dx_list <- list(
  c("SCZ", "CTRL")
)

# Save file suffix
save_suffix_vec <- c("SCZvsCon_DFC_bulkmegaset_Brainseq")
################################

# Create expr + GSEA objects (cached; HGNC-invariant, keyed by module set + bulk)
all_plots <- cached_all_plots(mods = mods, bulk_expr = bulk_expr,
                              these_mods = these_mods, datkme = datkme)

fig6_full_pipeline(all_plots = all_plots,
                   mods = mods,
                   these_mods = these_mods,
                   sea_dir = sea_dir,
                   mit_dir = mit_dir,
                   snapshot_out_dir = snapshot_out_dir,
                   keys = keys,
                   title_vec = title_vec,
                   title_vec_dataset = title_vec_dataset,
                   dx_list = dx_list,
                   save_suffix_vec = save_suffix_vec)


######
# Mods to feature in Fig 6
#####

# Load summary table genes (claude code literature search pipeline) for ad mtg, dfc
ad_db_mtg <- fread(data.table = F, file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.2/ad_db_summary_table_mtg.csv"))
ad_db_dfc <- fread(data.table = F, file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v7.2/ad_db_summary_table_dfc.csv"))

# Load module data
filter_under <- 3
datkme <- fread(data.table=F, file = file.path(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"),"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods_final <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]

# Find which mods contain summary table genes
hits <- lapply(mods, function(x) sum(x %in% ad_db_mtg[,1])) |> unlist()
these_hits <- which(hits > 0)

# Load lists of dCoPA mods
dmod_gab <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv")) |>
    filter(mod %in% these_mods,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    dplyr::select(c(mod, Celltype, Direction, Consistency))
dmod_liu <- fread(data.table=F, file=file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv")) |>
    filter(mod %in% these_mods,
           sig_FDR,
           Consistency %in% c(0, 1)) |>
    dplyr::select(c(mod, Celltype, Direction, Consistency))

# Create vector of most relevant AD genes (via chatgpt - /home/gugene/code/git/Consensus-analysis/Code_for_figures/fig_7/v7.2/agent_artifacts/alzheimers_genes_ranked.txt)
cg <- c("APOE","PSEN1", "APP", "PSEN2", "TREM2", "SORL1", "ABCA7", "BIN1", "CLU", "PICALM", "CR1", "CD33", "MS4A4A", "MS4A6A", "EPHA1", "INPP5D", "MEF2C", "FERMT2", "HLA-DRB1") 

# Check which mods have most highly relevant AD genes
hits2 <- lapply(mods, function(x) sum(x %in% cg)) |> unlist()
these_hits2 <- which(hits2 > 0)

# Check which mods have dCoPA genes AND most highly relevant AD genes
# > intersect(these_hits, these_hits2)
# [1]  83 153 # Two mods show up (index type: 1158)

# Analyze the two mods that contain both dCoPA genes and highly relevant AD genes
mods[[83]]
cg[cg %in% mods[[83]]]
# "MEF2C" (MEF2C not in summary table; NRXN1 is also in mod and in summary table)
# Shows up in top 10 genes

mods[[153]]
cg[cg %in% mods[[153]]]
# "SORL1"
# Doesn't show up in top 10 genes

# Find modules that contain dCoPA genes (DFC) and check which ones also have most highly relevant AD genes
hits_dfc <- lapply(mods, function(x) sum(x %in% ad_db_dfc[,1])) |> unlist()
these_hits_dfc <- which(hits > 0)
intersect(these_hits_dfc, these_hits2)
# 83 153 # Same modules as MTG

# Conclusion for which mods to feature in fig 6:
# use 83 for NRXN1 (MEF2C also shows up); use 17 for VPS35 (L6 IT) [1158 mod indices]

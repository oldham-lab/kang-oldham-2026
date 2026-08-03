# fig6_specific_modules.R
# Generate Fig 6 module-snapshot SVGs for a USER-SPECIFIED set of modules and a
# USER-SPECIFIED key (region / groupA / groupB / bulk dataset), instead of every
# dCoPA-shared module as in fig6_v4.R.
#
# Defaults: context "AD", modules 83 and 147, key c("MTC", "allAD", "Con", "bulk_megaset"),
#           highlight gene TMEM106B (red in the line-plot legend).
#
# Usage (either edit the parameter block below and source(), or pass args):
#   Rscript fig6_specific_modules.R
#   Rscript fig6_specific_modules.R 59,106
#   Rscript fig6_specific_modules.R 59,106 MTC,allAD,Con,bulk_megaset
#   Rscript fig6_specific_modules.R 59,106 DFC,Schizophrenia,control,bulk_megaset SCZ_ctrlmods
#   Rscript fig6_specific_modules.R 59,106 MTC,allAD,Con,bulk_megaset AD NRXN1,MEF2C
#
# arg 4 (optional): comma-separated gene names to colour red in the line-plot legend.
#
# context (arg 3) selects the data sources / dataset titles / default key:
#   "AD"           - SEAAD control modules, projected onto Gabitto + Liu  (default)
#   "SCZ_ctrlmods" - SEAAD control modules, projected onto CMC + SZBDMulti-Seq
#   "SCZ_sczmods"  - Brainseq SCZ modules,  projected onto CMC + SZBDMulti-Seq
#
# Output: <this folder>/<mod>_<key>.svg  (gitignored)

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

# Patch make_expr_line_plots in memory: fig6_fxns.R titles the snapshot by the
# loop index j (paste0("Module ", j)) rather than the true module ID, which
# mislabels sparse module IDs (e.g. module 992 shows as "Module 891"). Fix the
# title to the real ID here only — without editing the shared fig6_fxns.R.
local({
  .b   <- deparse(body(make_expr_line_plots))
  .pat <- 'paste0\\("Module ", *j\\)'
  if (!any(grepl(.pat, .b)))
    stop("title-label patch failed: 'paste0(\"Module \", j)' not found in make_expr_line_plots")
  .b <- gsub(.pat, 'paste0("Module ", these_mods[[j]])', .b)
  body(make_expr_line_plots) <<- parse(text = paste(.b, collapse = "\n"))[[1]]
})

# -------------------- user parameters --------------------
# Modules, key, and context. CLI args (if given) override these.
args    <- commandArgs(trailingOnly = TRUE)
context <- if(length(args) >= 3 && nzchar(args[3])) args[3] else "SCZ_ctrlmods"
# ---------------------------------------------------------

save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v4")

# Context switch: data sources, dataset titles, and default key per context.
if(context == "AD"){
  bulk_expr_path    <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv")
  datkme_path       <- file.path(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"), "kme_tables", "topmodposbc_table.csv")
  sea_dir           <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output")
  mit_dir           <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output"
  title_vec_dataset <- c("Gabitto et al. 2024", "Liu et al. 2025")
  default_key       <- c("MTC", "allAD", "Con", "bulk_megaset")
} else if(context == "SCZ_ctrlmods"){
  bulk_expr_path    <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv")
  datkme_path       <- file.path(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC"), "kme_tables", "topmodposbc_table.csv")
  sea_dir           <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC"
  mit_dir           <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq"
  title_vec_dataset <- c("CMC", "SZBDMulti-Seq")
  default_key       <- c("DFC", "Schizophrenia", "control", "bulk_megaset")
} else if(context == "SCZ_sczmods"){
  bulk_expr_path    <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SCZ_SampleNetworks/1_04-04-29/brainseq_samp_filt_SCZ_1_171_outliers_removed_geneSymbolsAdded.csv")
  datkme_path       <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/kme_tables/topmodposbc_table.csv")
  sea_dir           <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/CMC"
  mit_dir           <- "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/brainscope_means_SE_output/SZBDMulti-Seq"
  title_vec_dataset <- c("CMC", "SZBDMulti-Seq")
  default_key       <- c("DFC", "Schizophrenia", "control", "brainseq_scz")
} else {
  stop("Unknown context '", context, "'. Use one of: AD, SCZ_ctrlmods, SCZ_sczmods")
}

target_mods     <- if(length(args) >= 1 && nzchar(args[1])) as.numeric(strsplit(args[1], ",")[[1]]) else c(992)
key             <- if(length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]]            else default_key
# Genes to colour red in the line-plot legend (NULL = none). Edit here or pass as arg 4.
# fig8 panel f: no red highlight.
highlight_genes <- if(length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]]            else NULL

# Friendly labels for known tokens (cosmetic: legends, titles, output file name)
region_lookup <- c(PFC = "DFC", MTC = "MTG")
dx_lookup     <- c(allAD = "AD", Con = "CTRL", earlyAD = "Early AD", lateAD = "Late AD",
                   APOE44 = "APOE 4/4", APOE33 = "APOE 3/3",
                   Schizophrenia = "SCZ", control = "CTRL")
map_lab <- function(tok, lookup) if(tok %in% names(lookup)) lookup[[tok]] else tok

dx                <- c(map_lab(key[2], dx_lookup), map_lab(key[3], dx_lookup))  # c(condition, control)
keys              <- list(key)
title_vec         <- paste0(dx[2], " vs ", dx[1], " (", map_lab(key[1], region_lookup), ")")
dx_list           <- list(dx)
save_suffix_vec   <- paste(key, collapse = "_")

cat("Context:", context, "\n")
cat("Target modules:", paste(target_mods, collapse = ", "), "\n")
cat("Key:", paste(key, collapse = ", "), "\n")
cat("Highlight genes:", if(is.null(highlight_genes)) "(none)" else paste(highlight_genes, collapse = ", "), "\n")
cat("Output:", file.path(save_dir, paste0("<mod>_", save_suffix_vec, ".svg")), "\n\n")

################################

# Load bulk expression matrix
bulk_expr <- fread(bulk_expr_path, data.table = F)

# Load modules
filter_under <- 3
datkme <- fread(data.table = F, file = datkme_path)
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
mods <- mods[these_mods]

# Warn about any requested module not present in the module set
missing_in_set <- setdiff(target_mods, these_mods)
if(length(missing_in_set) > 0)
  warning("Requested module(s) absent from module set (size <= ", filter_under, " or not found): ",
          paste(missing_in_set, collapse = ", "))

# Create expr + GSEA plot objects.
# Computed over ALL modules so the GSEA FDR cutoff line is identical to the full figure.
all_plots <- list()
all_plots[[1]] <- make_expr_line_plots(expr = bulk_expr,
                                       mods = mods,
                                       these_mods = these_mods,
                                       datkme = datkme,
                                       highlight_genes = highlight_genes)

gsea_broad <- run_gsea_for_proj_optimized(mods, broad = T)
gsea       <- run_gsea_for_proj_optimized(mods, broad = F)

all_plots[[2]] <- make_gsea_plots(gsea = gsea,
                                  gsea_broad = gsea_broad,
                                  these_mods = these_mods)

# Generate snapshots for the target modules only, written flat into v4
fig6_full_pipeline(all_plots = all_plots,
                   mods = mods,
                   these_mods = these_mods,
                   sea_dir = sea_dir,
                   mit_dir = mit_dir,
                   snapshot_out_dir = save_dir,
                   keys = keys,
                   title_vec = title_vec,
                   title_vec_dataset = title_vec_dataset,
                   dx_list = dx_list,
                   save_suffix_vec = save_suffix_vec,
                   target_mods = target_mods,
                   flat_output = TRUE)

cat("\nDone.\n")

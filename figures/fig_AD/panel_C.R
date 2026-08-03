library(data.table)
library(tidyverse)
library(qs)
library(ComplexHeatmap)
options(bitmapType = 'cairo')
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "generic_enrichment_function.r"))


# Load subclass models
sub_lm <- qread(paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_topmodposFDR_Subclass.qs")))

# Load modules
modfdr <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv"))
mods <- tapply(modfdr$Gene, modfdr$topmodposfdr, list)

# Load vector of all genes used to produce modules
allgenes <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)[,2]

# Load significant modules
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_negative.qs")))
euc_fdr <- unlist(euc_dist[[3]][[3]])
euc_fdr_ct <- unlist(lapply(names(euc_fdr)[!duplicated(euc_fdr)], function(x) unlist(strsplit(x, "[.]"))[1]))

# Run enrichment
gshg_results <- GSHG_custom(mods[unique(euc_fdr)],
                            sub_lm[[1]], # donor 1
                            allgenes)

gshg_plot <- gshg_results[,-1]
rownames(gshg_plot) <- gsub("`", "", gshg_results[,1])
gshg_plot <- apply(gshg_plot, 2, function(x) -log10(x))

png(file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panel_C.png"), width=3000,height=1000)
Heatmap(as.matrix(gshg_plot), column_split = euc_fdr_ct,column_title_gp=gpar(fontsize=30), row_names_gp=gpar(fontsize=30), column_title_rot=45, show_column_names=F)
dev.off()

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))
library(cowplot)
library(tidyverse)
library(qs)
library(ggrepel)
options(bitmapType = 'cairo') 
library(doParallel)
library(rhdf5)
registerDoParallel(cores=8)

####################
# Load required data 
####################
cell_annoall <- fread(data.table=F,file=file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_V1.csv"))
cell_annoall <- cell_annoall[,c(3,1,2,4:8)] # just for consistency
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_V1.RDS")))
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
cell_exprall <- cell_exprall[rownames(cell_exprall) %in% genemap[,2],]
donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")

#san_mean <- apply(cell_exprall,1,mean)
#which.min(san_mean[san_mean>=1])
# TMCO3
# 2786

# calculate means by donor
san_mean_list <- list()
for(i in seq_along(donorvec)){
  df_donor <- cell_exprall[,colnames(cell_exprall) %in% cell_annoall$Cell_ID[cell_annoall$Donor==donorvec[i]]]
  san_mean_list[[i]] <- apply(df_donor, 1, mean)
}
qsave(san_mean_list,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_v1_10x/gene_count_means_byDonor.qs"))


###############################################
# create pseudobulk from individual lein donors (DFC)
###############################################

create_donorspecific_pseudobulk(cell_exprall = cell_exprall,
                                cell_annoall = cell_annoall,
                                donorvec = c("H18.30.002", "H19.30.001", "H19.30.002"), # vector of all unique donors
                                save_dirall = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/V1_indiv_donor/"),
                                pcnt.varvec = c(0), # vector of values for pcnt.var
                                no.samples = ncol(megaexpr)-2,
                                cell.nameindex = 3, # column index of cell_annoall that indicates cell name (matches colnames of cell_exprall)
                                filter = megaexpr[,1] # vector of gene symbols to filter pseudobulk matrix 
                                )

#################################################
# For each donor:
# Load Lein pseudobulk 
# Normalize by sample and scale by gene
# Model using subclasses and supertypes
# DLPFC
#################################################

model_donorspecific_pseudobulk(dirlist = c(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/V1_indiv_donor/")),
                               cell_annoall=cell_annoall,
                               cellannocolvec = c(1,4), # subclass, supertype,
                               cellIDcol = 3, # new_labs, not Cell_ID
                               celltypenames = c("subclass", "supertype"))






source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))
library(cowplot)
library(tidyverse)
library(qs)
library(ggrepel)
options(bitmapType = 'cairo') 
library(doParallel)
registerDoParallel(cores=8)

####################
# Load required data (scvi)
####################
save_dirall <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI/")
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
cell_exprall <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/scVI/denoised.csv"),data.table=F) %>% t
#cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
colnames(cell_exprall) <- cell_annoall$Cell_ID
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
orig_genes <- rownames(cell_exprall)
donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")

###############################################
# create pseudobulk from individual lein donors (DFC)
###############################################

create_donorspecific_pseudobulk(cell_exprall = cell_exprall,
                                cell_annoall = cell_annoall,
                                donorvec = c("H18.30.002", "H19.30.001", "H19.30.002"), # vector of all unique donors
                                save_dirall = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI/"),
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

model_donorspecific_pseudobulk(dirlist = c(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI/")),
                               cell_annoall=cell_annoall,
                               cellannocolvec = c(1,4), # subclass, supertype
                               cellIDcol = 3, # new_labs, not Cell_ID
                               celltypenames = c("subclass", "supertype"))

 
bla <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_ALRA/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_05-16-18_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS"))



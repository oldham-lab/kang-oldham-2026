# Create pb for LeinDFC_CB_scVI from all donors

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))
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
# Load required data (scvi)
####################
save_dirall <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_all_donors_scVI_cellbender/")
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC_cellbender.csv"), data.table = FALSE)
cell_exprall <- fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/scVI/denoised_cellbender.csv"),data.table=F) %>% t
# Add cell labels to cell_exprall
colnames(cell_exprall) <- cell_annoall$new_labs
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
orig_genes <- rownames(cell_exprall)

###############################################
# create pseudobulk from individual lein donors (DFC)
###############################################

no.samples = ncol(megaexpr)-2
cell.nameindex = 14 # column index of cell_annoall that indicates cell name (matches colnames of cell_exprall)
filter = megaexpr[,1] # vector of gene symbols to filter pseudobulk matrix 
orig_genes <- rownames(cell_exprall)
pcnt.var <- 0

if(!dir.exists(save_dirall)){dir.create(save_dirall,recursive=T)}
setwd(save_dirall)
makeSyntheticDatasets(
    cell_exprall,
    sampleindex=c(1:ncol(cell_exprall)),
    cell.info=cell_annoall,
    cell.name=cell.nameindex, # new labs
    cell.type=NULL,
    cell.frac=NULL,
    pcnt.cells=10,
    pcnt.var=pcnt.var,
    no.samples=no.samples, # size of megaset
    no.datasets=1
)
    
paths <- list.files(paste0(save_dirall,"/SyntheticDatasets/"),full.names=T)
paths <- paths[grep(paste0(pcnt.var,"pcntVar"),paths)]
paths <- paths[grep(paste0(no.samples,"samples"),paths)]
    
expr_lp_path <- paths[-grep("legend",paths)]
expr_lp_path <- expr_lp_path[grep(".csv",expr_lp_path)]
sifpath <- paths[grep("legend",paths)]
    
process_lein_pb(
    expr_lp_path, 
    orig_genes,
    filter = megaexpr[,1],
    filtername = "BULKGENESUBSET",
    save_simMat=F
)


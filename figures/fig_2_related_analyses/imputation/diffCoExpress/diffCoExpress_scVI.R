# Run differential co-expression FM between original Lein DFC and cellbender-scVI-corrected Lein DFC

library(tidyverse)
library(qs)
library(data.table)
library(reticulate)
numpy <- import("numpy")

######## SN data
# Load bulk megaset, lein pb, and lein cellbender-scVI-imputed pb
mega <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
leinpb <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_06-13-16_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS"))[[1]]
scvi <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI_cellbender/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_04-15-41_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS"))[[1]]

# Align genes 
common_genes <- intersect(intersect(mega[,2], leinpb[,2]), scvi[,2])
mega <- mega[mega[,2] %in% common_genes,]
mega <- mega[match(common_genes, mega[,2]),]
leinpb <- leinpb[leinpb[,2] %in% common_genes,]
leinpb <- leinpb[match(common_genes, leinpb[,2]),]
scvi <- scvi[scvi[,2] %in% common_genes,]
scvi <- scvi[match(common_genes, scvi[,2]),]

# Calculate adjacency matrix for differential co-expression
calcadj <- function(expr){
  expr <- expr[,3:ncol(expr)]
  simMat <- numpy$corrcoef(expr)
  colnames(simMat) <- rownames(expr)
  rownames(simMat) <- rownames(expr)
  simMat <- (simMat+1)/2
  diag(simMat) <- 0
  return(simMat)
}
simMatmega <- calcadj(mega)
qsave(simMatmega, file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/megamat_adjMat.qs"))
simMatleinpb <- calcadj(leinpb)
qsave(simMatleinpb, file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/leinpb_donor1_adjMat.qs"))
simMatscvi <- calcadj(scvi)
qsave(simMatscvi, file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/scvi_donor1_adjMat.qs"))


# Run FM
setwd(file.path(Sys.getenv("FINDMODULES_DIR", "/home/gugene/code/git/FindModules"), "FindModules/R/"))
source("FindModules.R")
source("map_identifiers_function.R")
source("FM_helper_fxns.R")
source("FindModules.R")
source("find_seed_genes_greedy_march_megaset.R")
source("similarityType.R")
source("plotting_functions.R")
source("overlapType.R")
source("networkOutputs.R")
source("module_quant_functions.R")
source("iteration_code.R")
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"))

# Run original bulk minus lein pseudobulk (find which modules are present in bulk but not lein pseudobulk)
simMat <- simMatmega - simMatleinpb
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

FindModules(
  projectname="donor1_megaMinusLeinPB",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)

# Run lein pseudobulk minus bulk megaset (for comparisons sake) 
simMat <- simMatleinpb - simMatmega 
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"))
FindModules(
  projectname="donor1_LeinPBMinusMega",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)

# Run real bulk minus cellbender+scVI (see which modules are recovered by denoise+imputation)
simMat <- simMatmega - simMatscvi
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"))
FindModules(
  projectname="donor1_megaMinus_cbscVI",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)


# Run cellbender+scVI minus real bulk (for comparisons sake)
simMat <- simMatscvi - simMatmega
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"))
FindModules(
  projectname="donor1_cbscVI_minus_mega",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)

# Run geneset enrichment
source(file.path(Sys.getenv("GSEA_GENERIC_DIR", "/home/gugene/code/git/GSEA_generic"), "GSEAfxsV3_nonpar_temp.r"))
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/"))
MyGSHGloop(kmecut1="seed")
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/"))
MyGSHGloop(kmecut1="seed")

# Collect p-values for AOMN
source(file.path(Sys.getenv("GSEA_GENERIC_DIR", "/home/gugene/code/git/GSEA_generic"), "GSEAfxsV3.r"))
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/"))
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7025"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7026"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7027"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7028"), 
                    whichKmeCut="seed",
                    write=TRUE)                    
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/"))
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7025"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7026"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7027"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7028"), 
                    whichKmeCut="seed",
                    write=TRUE)  

# Run FM with minMEcor = 0.9
simMat <- simMatmega - simMatleinpb
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

FindModules(
  projectname="donor1_megaMinusLeinPB_minMEcor0.9",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.9),
  merge.by = c("ME"),
  merge.param = 0.9,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)

simMat <- simMatmega - simMatscvi
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"))
FindModules(
  projectname="donor1_megaMinus_cbscVI_minMEcor0.9",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.9),
  merge.by = c("ME"),
  merge.param = 0.9,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)

# Run geneset enrichment
dir1 <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_minMEcor0.9_Modules/")
dir2 <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_minMEcor0.9_Modules/")
source(file.path(Sys.getenv("GSEA_GENERIC_DIR", "/home/gugene/code/git/GSEA_generic"), "GSEAfxsV3_nonpar_temp.r"))
setwd(dir1)
MyGSHGloop(kmecut1="seed")
setwd(dir2)
MyGSHGloop(kmecut1="seed")

# Collect p-values for AOMN
source(file.path(Sys.getenv("GSEA_GENERIC_DIR", "/home/gugene/code/git/GSEA_generic"), "GSEAfxsV3.r"))
setwd(dir1)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7025"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7026"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7027"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7028"), 
                    whichKmeCut="seed",
                    write=TRUE)                    
setwd(dir2)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7025"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7026"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7027"), 
                    whichKmeCut="seed",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7028"), 
                    whichKmeCut="seed",
                    write=TRUE)  


# Run PB minus PB_CB_scVI (and reverse)
simMat <- simMatscvi - simMatleinpb
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"))
FindModules(
  projectname="donor1_cbscVI_minus_leinpb",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)


# Run cellbender+scVI minus real bulk (for comparisons sake)
simMat <- simMatleinpb - simMatscvi
rownames(simMat) <- common_genes
colnames(simMat) <- common_genes
orig_mat1 <- mega[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"))
FindModules(
  projectname="donor1_leinpb_minus_cbscVI",
  data_cols=orig_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(orig_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = simMat,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(3,5,7),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.999, .99, .98,.97 ,.96, .95,.94, .93, .92,.91 ,.90),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default"
)

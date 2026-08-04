# Run regular FM on Bulk megaset, Lein PB, Lein PB + CB + scVI

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
# calcadj <- function(expr){
#   expr <- expr[,3:ncol(expr)]
#   simMat <- numpy$corrcoef(expr)
#   colnames(simMat) <- rownames(expr)
#   rownames(simMat) <- rownames(expr)
#   simMat <- (simMat+1)/2
#   diag(simMat) <- 0
#   return(simMat)
# }
# simMatmega <- calcadj(mega)
# qsave(simMatmega, file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/megamat_adjMat.qs"))
# simMatleinpb <- calcadj(leinpb)
# qsave(simMatleinpb, file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/leinpb_donor1_adjMat.qs"))
# simMatscvi <- calcadj(scvi)
# qsave(simMatscvi, file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/scvi_donor1_adjMat.qs"))

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

# Run FM
simMatmega <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/megamat_adjMat.qs"))
simMat <- simMatmega*2 - 1
diag(simMat) <- 0
colnames(simMat) <- mega[,2]
rownames(simMat) <- mega[,2]
orig_mat1 <- mega[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/"))
FindModules(
  projectname="mega",
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
  minSizevec = c(5,8,10,15),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.99, .98,.97 ,.96, .95, .93, .9),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

simMatleinpb <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/leinpb_donor1_adjMat.qs"))
simMat <- simMatleinpb*2 - 1
diag(simMat) <- 0
colnames(simMat) <- leinpb[,2]
rownames(simMat) <- leinpb[,2]
orig_mat1 <- leinpb[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/"))
FindModules(
  projectname="leinpb_donor1",
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
  minSizevec = c(5,8,10,15),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.99, .98,.97 ,.96, .95, .93, .9),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

simMatscvi <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/scvi_donor1_adjMat.qs"))
simMat <- simMatscvi*2 - 1
diag(simMat) <- 0
colnames(simMat) <- scvi[,2]
rownames(simMat) <- scvi[,2]
orig_mat1 <- scvi[,-1]

setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/"))
FindModules(
  projectname="pb_cb_scvi_donor1",
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
  minSizevec = c(5,8,10,15),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.99, .98,.97 ,.96, .95, .93, .9),
  minMEcorvec = c(0.85),
  merge.by = c("ME"),
  merge.param = 0.85,
  export.merge.comp = TRUE,
  loadTree = FALSE,
  writeKME = TRUE,
  writeModSnap = TRUE,
  modSnapExprVal = c("meanexpr"),
  floor = "default",
  prompt_zero_values=F
)

# Run geneset enrichment
source(file.path(Sys.getenv("GSEA_GENERIC_DIR", "/home/gugene/code/git/GSEA_generic"), "GSEAfxsV3_nonpar_temp.r"))
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/mega_Modules/"))
MyGSHGloop(kmecut1="topmodposbc")
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/mega_Modules/"))
MyGSHGloop(kmecut1="topmodposfdr")
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/leinpb_donor1_Modules/"))
MyGSHGloop(kmecut1="topmodposbc")
#setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/leinpb_donor1_Modules/"))
#MyGSHGloop(kmecut1="topmodposfdr")
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/pb_cb_scvi_donor1_Modules"))
MyGSHGloop(kmecut1="topmodposbc")
#setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/pb_cb_scvi_donor1_Modules"))
#MyGSHGloop(kmecut1="topmodposfdr")

# Collect p-values for AOMN
source(file.path(Sys.getenv("GSEA_GENERIC_DIR", "/home/gugene/code/git/GSEA_generic"), "GSEAfxsV3.r"))
dirlist <- c(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/mega_Modules/"),
             #file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/mega_Modules/"),
             file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/leinpb_donor1_Modules/"),
             #file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/leinpb_donor1_Modules/"),
             file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/pb_cb_scvi_donor1_Modules"))
             #file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/pb_cb_scvi_donor1_Modules"))
for(dir in dirlist){
  setwd(dir)
  getMySetEnrichments(simType="any",
                      whichSet=c("MOSET7025"), 
                      whichKmeCut="topmodposbc",
                      write=TRUE)
  getMySetEnrichments(simType="any",
                      whichSet=c("MOSET7026"), 
                      whichKmeCut="topmodposbc",
                      write=TRUE)
  getMySetEnrichments(simType="any",
                      whichSet=c("MOSET7027"), 
                      whichKmeCut="topmodposbc",
                      write=TRUE)
  getMySetEnrichments(simType="any",
                      whichSet=c("MOSET7028"), 
                      whichKmeCut="topmodposbc",
                      write=TRUE)    
} 

# Compare enrichment p-values between datasets
fmdirs <- list.dirs(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM"), recursive=F, full.names=T)
fmtypes <- c("PB", "Bulk", "PB_CB_scVI")
mosetvec <- c("MOSET7025","MOSET7026","MOSET7027","MOSET7028")
mosetvecnames <- c("MOSET7025\n(Astro)","MOSET7026\n(Oligo)","MOSET7027\n(Micro)","MOSET7028\n(Neuron)")
outlist <- list()
for(i in seq_along(fmdirs)){
  outvec <- list()
  for(j in seq_along(mosetvec)){
    enrich <- list.files(fmdirs[i], full.names=T)[grep(mosetvec[j], list.files(fmdirs[i]))] %>% fread(., data.table=F)
    outvec[[j]] <- data.frame("type"=fmtypes[i], "pval"=enrich$Pvalue[1], "moset"=mosetvec[j], "mosetname"=mosetvecnames[j])
  }
  outlist[[i]] <- do.call(rbind, outvec)
}
outdf <- do.call(rbind, outlist) %>%
  mutate(logpval = -log10(pval),
         type=factor(type, levels=c("PB_CB_scVI", "PB", "Bulk")))

p <- ggplot(outdf, aes(y=type, x=logpval)) + 
  theme_bw() +
  geom_bar(stat="identity") + 
  facet_wrap(~mosetname, ncol=2,nrow=2) + 
  theme(text=element_text(size=30)) +
  labs(x=bquote(-log[10]~"P-value"), y="")
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_FM/aomn_enrichment_comparison_bonf.png"), height=6, width=8)



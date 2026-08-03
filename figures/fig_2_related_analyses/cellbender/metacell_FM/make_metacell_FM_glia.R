setwd("/home/gugene/code/git/FindModules/FindModules/R/")
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

library(tidyverse)
library(qs)
library(data.table)

####################
# Load required data (DFC)
####################
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))

# Subset to glial cells only
cell_annoall <- cell_annoall %>% dplyr::filter(Class=="Non-neuronal")
cell_exprall <- cell_exprall[,colnames(cell_exprall) %in% cell_annoall$Cell_ID]

# Create metacells from subtypes
subtype_list <- tapply(cell_annoall$Cell_ID, cell_annoall$Cluster,list)

metacell_list <- lapply(subtype_list, function(x){
    if(sum(colnames(cell_exprall) %in% x)==1){
      return(cell_exprall[, colnames(cell_exprall) %in% x])
    } else {
      return(rowMeans(cell_exprall[, colnames(cell_exprall) %in% x]))
    }  
})
metacell_mat <- do.call(cbind, metacell_list)
qsave(metacell_mat,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/gluta_metacell_expr_mat_nonNeuronal.qs"))


# Run FM
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/"))
metacell_mat1 <- qread(file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/gluta_metacell_expr_mat_nonNeuronal.qs"))
metacell_mat1 <- data.frame("Gene"=rownames(metacell_mat1), metacell_mat1)
# Filter genes to a reasonable subset (to remove singlets)
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
metacell_mat1 <- metacell_mat1[metacell_mat1$Gene %in% megaexpr$Gene,]

FindModules(
  projectname="metacell_nonNeuronal",
  data_cols=metacell_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(metacell_mat1),
  sampleGroups = NULL,
  subset = c(1:nrow(metacell_mat1)),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  #minSizevec = c(8, 10, 12, 15),
  minSizevec = c(5, 8, 10, 12, 15),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999, .999, .99, .98,.97 ,.96, .95),
  #signumvec = c(.95),
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
#.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4/")
WD <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/metacell_nonNeuronal_Modules")
source("/home/gugene/code/git/GSEA_generic/GSEAfxsV3_nonpar_temp.r")
setwd(WD)
#MyGSHGloop(kmecut1="seed")
# setwd(WD)
MyGSHGloop(kmecut1="topmodposbc")
setwd(WD)
MyGSHGloop(kmecut1="topmodposfdr")
setwd(WD)

source("/home/gugene/code/git/GSEA_generic/GSEAfxsV3.r")
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET7", "MOSET8", "MOSET6808", "MOSET6837"), 
                    whichKmeCut="topmodposbc",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET9"), 
                    whichKmeCut="topmodposbc",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET6808"), 
                    whichKmeCut="topmodposbc",
                    write=TRUE)
getMySetEnrichments(simType="any",
                    whichSet=c("MOSET6837"), 
                    whichKmeCut="topmodposbc",
                    write=TRUE)

# Load seed genes of largest network
seedkme <- fread(data.table=F, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/metacell_gluta_Modules/Bicor-no_TO_signum0.64_minSize8_merge_ME_0.85_17675/kME_table_.csv"))
seedkme <- tapply(seedkme[,1], seedkme[,2],list)
# Load module statistics of largest network
modstats <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/metacell_gluta_Modules/Bicor-no_TO_signum0.64_minSize8_merge_ME_0.85_17675/Module_statistics.csv"))
exclude_mods <- modstats$Module[modstats$ZMeanExpr<0]
exclude_genes <- unlist(seedkme[names(seedkme) %in% exclude_mods])

# Run FM again
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/"))
metacell_mat1 <- qread(file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/gluta_metacell_expr_mat.qs"))
metacell_mat1 <- data.frame("Gene"=rownames(metacell_mat1), metacell_mat1)
# Filter genes to a reasonable subset (to remove singlets)
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
metacell_mat1 <- metacell_mat1[metacell_mat1$Gene %in% megaexpr$Gene,]

FindModules(
  projectname="metacell_gluta_filterLowMods",
  data_cols=metacell_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(metacell_mat1),
  sampleGroups = NULL,
  subset = c(which(!metacell_mat1[,1] %in% exclude_genes)),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  #minSizevec = c(8, 10, 12, 15),
  minSizevec = c(5, 8, 10, 12, 15),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999, .999, .99, .98,.97 ,.96, .95),
  #signumvec = c(.95),
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


####################
# Load required data (DFC Cellbender)
####################
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
cell_exprall <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Raw_data_bag_4_Human_Cross_Areal_raw_10x/data/DLPFC/aligned_data/lein_2023_dfc_cellbender_expr.qs"))
# Create new cell IDs for cellbender data so we can match to existing cell_anno data
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Human_Cross_Areal_raw_10x_manifest.tsv"), data.table=F) %>%
  dplyr::filter(Anatomical_site=="DLPFC") %>%
  dplyr::select(File_name,Sample_ID)
sif <- sif[1:12,]
sif$File_name <- gsub("_S01_L003_I1_001.fastq.gz","", sif$File_name)
colnames(cell_exprall) <- lapply(colnames(cell_exprall) ,function(x){
  fn <- strsplit(x,"_")[[1]]
  fn <- paste0(fn[1],"_", fn[2])
  y <- strsplit(x, "-")[[1]][[2]]
  z <- strsplit(y, "_")[[1]][[2]]
  return(paste0(z, "-", sif$Sample_ID[sif$File_name==fn]))
}) %>% unlist
cell_annoall$new_labs <- cell_annoall$Cell_ID %>%
  lapply(function(x){
    y <- strsplit(x, "-")[[1]][[1]]
    z <- strsplit(x, "-")[[1]][[2]]
    z <- stringr::str_extract(z, "L8TX.*")
    return(paste0(y, "-", z))
  }) %>% unlist
cell_exprall <- cell_exprall[,colnames(cell_exprall) %in% cell_annoall$new_labs]
cell_annoall <- cell_annoall[match(colnames(cell_exprall), cell_annoall$new_labs),]
fwrite(cell_annoall,file=file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC_cellbender.csv"))

# Subset to excitatory neurons
cell_annoall <- cell_annoall %>% dplyr::filter(Class=="Glutamatergic")
cell_exprall <- cell_exprall[,colnames(cell_exprall) %in% cell_annoall$new_labs]
cell_exprall <- as.matrix(cell_exprall)

# Create metacells from subtypes
subtype_list <- tapply(cell_annoall$new_labs, cell_annoall$Cluster,list)

metacell_list <- lapply(subtype_list, function(x){
    if(sum(colnames(cell_exprall) %in% x)==1){
      return(cell_exprall[, colnames(cell_exprall) %in% x])
    } else if(sum(colnames(cell_exprall) %in% x)>1){
      return(rowMeans(cell_exprall[, colnames(cell_exprall) %in% x]))
    } else {
      return(data.frame())
    }
})
metacell_mat <- do.call(cbind, metacell_list)
qsave(metacell_mat,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/gluta_metacell_expr_mat_cellbender.qs"))


# Run FM
setwd(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/"))
metacell_mat1 <- qread(file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/gluta_metacell_expr_mat_cellbender.qs"))
metacell_mat1 <- data.frame("Gene"=rownames(metacell_mat1), metacell_mat1)
# Filter to genes in uncorrected Lein DFC
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
metacell_mat1 <- metacell_mat1[metacell_mat1[,1] %in% genemap[,2],]
# Filter genes to a reasonable subset (to remove singlets)
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
metacell_mat1 <- metacell_mat1[metacell_mat1$Gene %in% megaexpr$Gene,]

FindModules(
  projectname="metacell_gluta_cellbender",
  data_cols=metacell_mat1,
  genes=NULL,
  metadata_cols=1,
  sampleIndex=2:ncol(metacell_mat1),
  sampleGroups = NULL,
  subset = c("None"),
  simMat = NULL,
  saveSimMat = FALSE,
  simType = c("Bicor"),
  overlapTO = FALSE,
  TOtype = c("default"),
  TOdenom = c("default"),
  beta = 1,
  iterate = TRUE,
  minSizevec = c(5, 8, 10, 12, 15),
  greedyMarch = FALSE,
  signumType = c("rel"),
  signumvec = c(.9999, .999, .99, .98,.97 ,.96, .95),
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
#.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4/")
WD <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_metacell_FM/metacell_gluta_cellbender_Modules")
source("/home/gugene/code/git/GSEA_generic/GSEAfxsV3_nonpar_temp.r")
setwd(WD)
#MyGSHGloop(kmecut1="seed")
# setwd(WD)
MyGSHGloop(kmecut1="topmodposbc")
setwd(WD)
MyGSHGloop(kmecut1="topmodposfdr")
setwd(WD)

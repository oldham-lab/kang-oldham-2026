# This approach doesn't really work since sctransform produces negative values after the scaling

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))
library(cowplot)
library(tidyverse)
library(qs)
library(sctransform)
library(Seurat)
options(bitmapType = 'cairo') 

####################
# Load required data
####################
save_dirall <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_sct/")
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)

# Run SCTransform on Lein DFC data
#cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
#normalized_data <- cell_exprall %>% CreateSeuratObject() %>% SCTransform() 
#qsave(normalized_data[["SCT"]]@scale.data, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_sct/lein_dfc_sct.qs"))

cell_exprall <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_sct/lein_dfc_sct.qs"))
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")

# san_mean <- apply(cell_exprall,1,mean)
# which.min(san_mean[san_mean>=1])
# STAMBPL1 
# 3262 

###############################################
# create pseudobulk from individual lein donors
###############################################

for(i in 1:3){ # for each donor
  for(j in c(0)){
    save_dir <- paste0(save_dirall,"/donor",i)
    if(!dir.exists(save_dir)){dir.create(save_dir, recursive=T)}
    cell_expr<- cell_exprall[,cell_annoall$Donor==donorvec[i]]
    cell_anno <- cell_annoall[cell_annoall$Donor==donorvec[i],]
    
    pcnt.cells=10
    pcnt.var=j
    no.samples=ncol(megaexpr)-2
    genes=megaexpr[,1]
    setwd(save_dir)
    makeSyntheticDatasets(
      cell_expr,
      sampleindex=1:ncol(cell_expr),
      cell.info=cell_anno,
      cell.name=3,
      cell.type=NULL,
      cell.frac=NULL,
      pcnt.cells=pcnt.cells,
      pcnt.var=pcnt.var,
      no.samples=no.samples, # size of megaset
      no.datasets=1
    )
    
    paths <- list.files(paste0(save_dir,"/SyntheticDatasets/"),full.names=T)
    paths <- paths[grep(paste0(pcnt.var,"pcntVar"),paths)]
    paths <- paths[grep(paste0(no.samples,"samples"),paths)]
    
    expr_lp_path <- paths[-grep("legend",paths)]
    expr_lp_path <- expr_lp_path[grep(".csv",expr_lp_path)]
    sifpath <- paths[grep("legend",paths)]
    
    # match expr feature names to lein 10x
    expr_lp <- fread(expr_lp_path,data.table=F)
    expr_lp <- expr_lp[match(genemap[,2], rownames(cell_exprall)),]

    process_lein_pb(
      expr_lp, 
      filter = genes,
      filtername = "BULKGENESUBSET",
      save_simMat=F
    )
    cat(j, " ")
  }
  cat(i, "i done\n")
}

#################################################
# For each donor:
# Load Lein pseudobulk 
# Normalize by sample and scale by gene
# Model using subclasses and supertypes
#################################################

plotdflist <- list()
cellannocolvec <- c(1,4) # subclass, supertype
for(cellannocol in seq_along(cellannocolvec)){ 
for(i in 1:3){
  exdir <- list.files(paste0(save_dirall,"/donor",i,"/SyntheticDatasets/"),full.names=T)
  exdirexpr <- exdir[grep("EXPRLIST",exdir)]
  exdirsif <- exdir[grep("samples_legend",exdir)]
  exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[9]))
  exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  
  exdirlist <- mapply(function(x,y,z){
    obj <- readRDS(x) 
    obj <- norm_pb_samples(obj[[1]], index=3)
    exprt <- t(obj[,3:ncol(obj)])
    escale <- apply(exprt, 2, scale)
    rownames(escale) <- rownames(exprt)
    escaleperm <- apply(escale, 2,function(x) sample(x, length(x))) # permuted matrix
    genemeansfilt <- apply(exprt,2,function(x) log10(mean(x)))
    
    sif <- fread(data.table=F,file=y)
    sampct <- t(calc_CT_counts(sif, cell_annoall, ct_col=cellannocolvec[cellannocol], cell_id_col = 3))
    sampctmean <- apply(sampct, 2, mean)
    sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
    
    # Run modeling
    mod_list <- lapply(as.list(as.data.frame(escale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
    mod_r2a <- unlist(lapply(mod_list, function(x) summary(x)$r.squared))
    plotdf <- data.frame("adj_r2" = mod_r2a, "mean" = genemeansfilt,"type"="Real",
                         "pcnt.var"=z)
    return(plotdf)
  }, exdirexpr,exdirsif,exdirnames, SIMPLIFY=F)
  plotdflist[[i]] <- do.call(rbind,exdirlist) %>% remove_rownames
} # for i
if(cellannocol==1){
  saveRDS(plotdflist,file=paste0(save_dirall,"/fig1bdflist_subclass.RDS"))
} else {
  saveRDS(plotdflist,file=paste0(save_dirall,"/fig1bdflist_supertype.RDS"))
}

} # for cellannocol





# Before creating pseudobulk, convert UMI counts to CPM
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))
library(cowplot)
library(tidyverse)
library(qs)
options(bitmapType = 'cairo') 

library(doParallel)
registerDoParallel(cores=8)

####################
# Load required data
####################
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")

# Convert cell expression matrix to cpm
cell_exprall <- apply(cell_exprall,2,function(x){
  x/sum(x) * 1e6
})

###############################################
# create pseudobulk from individual lein donors using 50:50 neuron:glia ratio
###############################################

save_dirall <- paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_cpm/"))

for(i in 1:3){ # for each donor
#for(i in seq_along(donorvec)){ # for each donor
  for(j in c(0)){
    save_dir <- paste0(save_dirall,"/donor",i)
    if(!dir.exists(save_dir)){dir.create(save_dir, recursive=T)}
    cell_expr <- cell_exprall[,cell_annoall$Donor==unique(cell_annoall$Donor)[i]]
    cell_anno <- cell_annoall[cell_annoall$Donor==unique(cell_annoall$Donor)[i],]
    
    pcnt.cells=10
    pcnt.var=j
    no.samples=ncol(megaexpr)-2
    genes=megaexpr[,1]
    setwd(save_dir)
    makeSyntheticDatasets(
      cell_expr,
      sampleindex=c(1:ncol(cell_expr)),
      cell.info=cell_anno,
      cell.name=3,
      cell.type=NULL, # Class
      cell.frac=NULL, # ct_per
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
    
    process_lein_pb(
      expr_lp_path, 
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

dirlist <- c(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_cpm/"))
 
for(save_dirall in dirlist){
  plotdflist <- list()
  cellannocolvec <- c(1,4) # subclass, supertype
  for(cellannocol in seq_along(cellannocolvec)){ 
    for(i in 1:3){
      exdir <- list.files(paste0(save_dirall,"/donor",i,"/SyntheticDatasets/"),full.names=T)
      exdirexpr <- exdir[grep("EXPRLIST",exdir)]
      exdirexpr <- exdirexpr[grep("0pcntVar", exdirexpr)]
      exdirsif <- exdir[grep("samples_legend",exdir)]
      exdirsif <- exdirsif[grep("0pcntVar", exdirsif)]
      #exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[9]))
      #exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
      exdirnames <- 0

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
        #mod_list <- lapply(as.list(as.data.frame(escale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
        mod_r2a <- foreach(l=1:ncol(escale)) %dopar% {
          summary(lm(escale[,l] ~ ., data=as.data.frame(sampct)))$adj.r.squared
        } %>% unlist
        
        #mod_r2a <- unlist(lapply(mod_list, function(x) summary(x)$r.squared))
      
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
  cat(save_dirall, " done\n\n\n")
} # for save_dirall






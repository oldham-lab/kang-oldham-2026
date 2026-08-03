source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))
library(cowplot)
library(tidyverse)
library(qs)
options(bitmapType = 'cairo') 
library(doParallel)
library(Seurat)
library(SeuratObject)
registerDoParallel(cores=8)

####################
# Load required data (DFC, cellbender corrected)
####################
save_dirall <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_cellbender/")
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
cell_exprall <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Raw_data_bag_4_Human_Cross_Areal_raw_10x/data/DLPFC/aligned_data/lein_2023_dfc_cellbender_expr.qs"))
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
#san_mean <- apply(cell_exprall,1,mean)
#which.min(san_mean[san_mean>=1])
# SHLD2
# 2203
cell_exprall <- cell_exprall[rownames(cell_exprall) %in% genemap[,2],]
donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")
orig_genes <- rownames(cell_exprall)

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

# Create hdf5 object for imputation algorithms (scVI in particular)
#library(SeuratDisk)
#dfc <- CreateSeuratObject(counts = cell_exprall, project = "leindfc", min.cells = 3, min.features = 200)
#SaveH5Seurat(dfc, filename = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Raw_data_bag_4_Human_Cross_Areal_raw_10x/data/DLPFC/aligned_data/lein_2023_dfc_cellbender_expr_matchedLabels.h5Seurat"))
#Convert(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Raw_data_bag_4_Human_Cross_Areal_raw_10x/data/DLPFC/aligned_data/lein_2023_dfc_cellbender_expr_matchedLabels.h5Seurat"), dest = "h5ad")

#fwrite(as.matrix(cell_exprall),row.names=T,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Raw_data_bag_4_Human_Cross_Areal_raw_10x/data/DLPFC/aligned_data/lein_2023_dfc_cellbender_expr_matchedLabels.csv"))

# calculate means by donor
san_mean_list <- list()
for(i in seq_along(donorvec)){
  df_donor <- cell_exprall[,colnames(cell_exprall) %in% cell_annoall$new_labs[cell_annoall$Donor==donorvec[i]]]
  san_mean_list[[i]] <- apply(df_donor, 1, mean)
}
qsave(san_mean_list,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/cellbender/gene_count_means_byDonor.qs"))



###############################################
# create pseudobulk from individual lein donors (DFC)
###############################################

for(i in 1:3){ # for each donor
  for(j in c(0)){
    save_dir <- paste0(save_dirall,"/donor",i)
    if(!dir.exists(save_dir)){dir.create(save_dir, recursive=T)}
    cell_expr <- cell_exprall[,cell_annoall$Donor==donorvec[i]]
    cell_anno <- cell_annoall[cell_annoall$Donor==donorvec[i],]

    pcnt.var=j
    no.samples=ncol(megaexpr)-2
    
    setwd(save_dir)
    makeSyntheticDatasets(
      cell_expr,
      sampleindex=c(1:ncol(cell_expr)),
      cell.info=cell_anno,
      cell.name=14, # new labs
      cell.type=NULL,
      cell.frac=NULL,
      pcnt.cells=10,
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
      orig_genes,
      filter = megaexpr[,1],
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
# DLPFC
#################################################

dirlist <- c(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_cellbender/"))
 
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
        sampct <- t(calc_CT_counts(sif, cell_annoall, ct_col=cellannocolvec[cellannocol], cell_id_col = 14))
        sampctmean <- apply(sampct, 2, mean)
        sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
        
        # Run modeling
        mod_list <- foreach(l=1:ncol(escale)) %dopar% {
          lm(escale[,l] ~ ., data=as.data.frame(sampct))
        } 
 
        mod_r2a <- unlist(lapply(mod_list, function(x) summary(x)$adj.r.squared))
        mod_rmse <- unlist(lapply(mod_list, function(x) sqrt(mean(x$residuals^2))))
        rm(mod_list)
        gc()
                
        plotdf <- data.frame("adj_r2" = mod_r2a,"rmse"=mod_rmse, "mean" = genemeansfilt,"type"="Real",
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






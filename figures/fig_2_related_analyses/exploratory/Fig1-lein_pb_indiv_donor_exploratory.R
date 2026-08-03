# making pb with same settings but zero var

save_dirall <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_testing/")
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_BrainGVEX_SampleNetworks/1_11-32-03/combined_FCX_minus_Brainseq_BrainGVEX_1_1323_ComBat.csv"), data.table=F)
origdemat2 <- fread(data.table=F, file="~/!softlinks/ANALYSES_for_bulk_megaset_and_snrnaseq/6-SN-pseudobulk/lein_orig_de_genes.csv")
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)

# Split cell expr by donor and create pseudobulk

for(i in 1:length(unique(cell_annoall$Donor))){
  save_dir <- paste0(save_dirall,"/donor",i)
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  cell_expr<- cell_exprall[,cell_annoall$Donor==unique(cell_annoall$Donor)[i]]
  cell_anno <- cell_annoall[cell_annoall$Donor==unique(cell_annoall$Donor)[i],]
  # pb_wrapper(savedir1,
  #            exprtrim,
  #            siftrim,
  #            pcnt.cells=10,
  #            pcnt.var=23,
  #            no.samples=75,
  #            genes=megaexpr[,1],
  #            demat = origdemat2,
  #            search_by="no.samples")
  
  pcnt.cells=10
  pcnt.var=23
  no.samples=1788
  genes=megaexpr[,1]
  demat=origdemat2
  search_by="no.samples"
  setwd(save_dir)
  makeSyntheticDatasets(
    cell_expr,
    sampleindex=c(1:ncol(cell_expr)),
    cell.info=cell_anno,
    cell.name=3,
    cell.type=NULL,
    cell.frac=NULL,
    pcnt.cells=pcnt.cells,
    pcnt.var=0,
    no.samples=no.samples, # size of megaset
    no.datasets=1
  )
  
  paths <- list.files(paste0(save_dir,"/SyntheticDatasets/"),full.names=T)
  if(search_by=="pcnt.var"){
    paths <- paths[grep(paste0(pcnt.var,"pcntVar"),paths)]
  } else {
    paths <- paths[grep(paste0(no.samples,"samples"),paths)]
  }
  expr_lp_path <- paths[-grep("legend",paths)]
  expr_lp_path <- expr_lp_path[grep(".csv",expr_lp_path)]
  sifpath <- paths[grep("legend",paths)]
  
  process_lein_pb(
    expr_lp_path, 
    filter = genes,
    filtername = "BULKGENESUBSET"
  )
}




i=1
exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
exdirexpr <- exdir[grep("EXPRLIST",exdir)]
lein75 <- readRDS(exdirexpr)
lein75 <- norm_pb_samples(lein75[[1]], index=3)
exprt <- t(lein75[,3:ncol(lein75)]) 
escale <- apply(exprt, 2, scale)
rownames(escale) <- rownames(exprt)
escaleperm <- apply(escale, 2,function(x) sample(x, length(x))) # permuted matrix
genemeansfilt <- apply(exprt,2,function(x) log10(mean(x)))

# Load cell abundance for lein pseudobulk (Supertype)
exdirsif <- exdir[grep("samples_legend",exdir)]
lein75sif <- fread(data.table=F,file=exdirsif)
sampct <- t(calc_CT_counts(lein75sif, cell_anno, ct_col=4, cell_id_col = 3))
sampctmean <- apply(sampct, 2, mean)
sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subclass

leintest <- readRDS(exdirexpr)[[1]]
leintestscale <- apply(t(leintest[,3:ncol(leintest)]),2,scale)
rownames(leintestscale) <- rownames(exprt)

# Run modeling
leinmod_list <- lapply(as.list(as.data.frame(escale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
leinmod_r2a <- unlist(lapply(leinmod_list, function(x) summary(x)$r.squared))

test1 <- lapply(as.list(as.data.frame(leintestscale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
test2 <- unlist(lapply(test1, function(x) summary(x)$r.squared))


i=1
exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_testing/donor"),i,"/SyntheticDatasets/"),full.names=T)
exdirexpr <- exdir[grep("EXPRLIST",exdir)]
lein75 <- readRDS(exdirexpr)
lein75 <- norm_pb_samples(lein75[[1]], index=3)
exprt <- t(lein75[,3:ncol(lein75)]) 
escale <- apply(exprt, 2, scale)
rownames(escale) <- rownames(exprt)
escaleperm <- apply(escale, 2,function(x) sample(x, length(x))) # permuted matrix
genemeansfilt <- apply(exprt,2,function(x) log10(mean(x)))

# Load cell abundance for lein pseudobulk (Supertype)
exdirsif <- exdir[grep("samples_legend",exdir)]
lein75sif <- fread(data.table=F,file=exdirsif)
sampct <- t(calc_CT_counts(lein75sif, cell_anno, ct_col=4, cell_id_col = 3))
sampctmean <- apply(sampct, 2, mean)
sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subclass

leintest <- readRDS(exdirexpr)[[1]]
leintestscale <- apply(t(leintest[,3:ncol(leintest)]),2,scale)
rownames(leintestscale) <- rownames(exprt)

# Run modeling
leinmodbla_list <- lapply(as.list(as.data.frame(escale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
leinmodbla_r2a <- unlist(lapply(leinmodbla_list, function(x) summary(x)$r.squared))

test1bla <- lapply(as.list(as.data.frame(leintestscale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
test2bla <- unlist(lapply(test1bla, function(x) summary(x)$r.squared))

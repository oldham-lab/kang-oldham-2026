# Concatenate and process all datasets together
library(AnnotationHub)


dat_names <- c("Brainseq",
               "GTEx",
               "NABEC",
               "ROSMAP",
               "BrainGVEX",
               "CMC",
               "CMC_HBCC")

exprs <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SampleNetworks/1_11-53-31/brainseq_samp_filt_1_190_outliers_removed.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/GTEx_samp_filt_SampleNetworks/1_09-10-16/GTEx_samp_filt_1_190_outliers_removed.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/NABEC_samp_filt_SampleNetworks/1_04-53-24/NABEC_samp_filt_1_69_outliers_removed.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_SampleNetworks/1_08-13-25/ROSMAP_samp_filt_1_347_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/BrainGVEX_samp_filt_SampleNetworks/1_08-25-36/BrainGVEX_samp_filt_1_275_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_samp_filt_SampleNetworks/1_11-19-13/CMC_samp_filt_1_285_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_HBCC_samp_filt_SampleNetworks/1_10-40-08/CMC_HBCC_samp_filt_1_162_ComBat.csv"))
expr_list <- list()
for(i in 1:length(exprs)){
 expr_list[[i]] <- fread(exprs[i], data.table=F)
}

# Find common features
for(i in 1:length(expr_list)){
  if(i==1){
    common_features <- expr_list[[i]][,1]
  } else {
    common_features <- intersect(common_features, expr_list[[i]][,1])
  }
}
common_features <- sort(common_features)

# Filter expr matrices to common features
for(i in 1:length(expr_list)){
  expr_list[[i]] <- expr_list[[i]][expr_list[[i]][,1] %in% common_features,]
}

# Order expr matrices by feature
for(i in 1:length(expr_list)){
  expr_list[[i]] <- expr_list[[i]][match(common_features,expr_list[[1]][,1]),]
}

# Combine expr matrices and create sample info mat
for(i in 1:length(expr_list)){
  if(i==1){
    mega_expr <- expr_list[[i]]
    mega_sif <- data.frame("Dataset" = dat_names[i],
                           "samp_id" = colnames(expr_list[[i]])[2:ncol(expr_list[[i]])])
  } else {
    mega_expr <- cbind(mega_expr, expr_list[[i]][,2:ncol(expr_list[[i]])])
    mega_sif <- rbind(mega_sif, data.frame("Dataset" = dat_names[i],
                                           "samp_id" = colnames(expr_list[[i]])[2:ncol(expr_list[[i]])]))
  }
}

colnames(mega_expr)[1] <- "gene_id"

# Add gene symbol to expr matrix
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, c("GENEID", "GENENAME"), "TXNAME")

genes <- ahEdbgene[ahEdbgene[,2] %in% mega_expr[,1],c(2,3)]
genes <- genes[!duplicated(genes[,1]),]
# Many ENSG ids have no associated gene symbol; remove these
genes <- genes[genes[,2] != "",]
genes <- genes[!duplicated(genes[,2]),]
genes <- genes[order(genes[,1]),]

mega_expr <- mega_expr[mega_expr[,1] %in% genes[,1],]
mega_expr <- mega_expr[match(genes[,1], mega_expr[,1]),]
mega_expr <- cbind(genes, mega_expr[,2:ncol(mega_expr)])
colnames(mega_expr)[1:2] <- c("ensembl_id", "Gene")

fwrite(mega_expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined.csv"))
fwrite(mega_sif, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined.csv"))


source(file.path(Sys.getenv("SAMPLENETWORK_DIR", "/home/gugene/code/labcode_old/SampleNetwork"), "SampleNetwork_1.08.r"))

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/"))
mega_expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined.csv"), data.table=F)
mega_sif$grouplabels1 <- 1

SampleNetwork(
  datExprT=mega_expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=2,
  indices1=list(seq(3,ncol(mega_expr))),
  subgroup1=1, # color samples by Dataset
  sampleinfo1=mega_sif,
  samplelabels1=2,
  grouplabels1=which(colnames(mega_sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(1), # Dataset
  trait1=NULL,
  asfactors1=c(1),
  projectname1="combined_FCX_final",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

# library(edgeR)
# 
# # apply normalization
# expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
# tmmf <- calcNormFactors(expr[,3:ncol(expr)], method = "TMM")
# 
# expr2 <- expr
# for(i in 3:ncol(expr)){
#   expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
# }
# colnames(expr2)[1] <- "gene_id"
# 
# fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat_ctf_asinh.csv"))
# 
# 
# # scale first then apply normalization
# expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
# si <- 3:ncol(expr)
# expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
# tmmf <- calcNormFactors(expr[,si], method = "TMM")
# expr2 <- expr
# for(i in 3:ncol(expr)){
#   expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
# }
# colnames(expr2)[1] <- "gene_id"
# fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat_scale_ctf_asinh.csv"))
# 
# expr3 <- expr
# for(i in 3:ncol(expr)){
#   expr3[,i] <- expr3[,i]/tmmf[i-2]
# }
# colnames(expr3)[1] <- "gene_id"
# fwrite(expr3, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat_scale_ctf.csv"))

# Filter samples again
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
mega_sif <- fread(data.table=F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined.csv"))
mega_sif$grouplabels1 <- 1

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/"))
SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=2,
  indices1=list(seq(3,ncol(expr))),
  subgroup1=1, # color samples by Dataset
  sampleinfo1=mega_sif,
  samplelabels1=2,
  grouplabels1=which(colnames(mega_sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(1), # Dataset
  trait1=NULL,
  asfactors1=c(1),
  projectname1="combined_FCX_final",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# No samples were filtered
# calculate and save simMat
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-48-30/combined_FCX_final_1_1518_outliers_removed.csv"), data.table=F)
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined.csv"), data.table=F)
sif <- sif[sif$samp_id %in% colnames(expr),]
sif[,1] <- as.factor(sif[,1])

simMat <- bicor(t(expr[,3:ncol(expr)]), use="p")
diag(simMat) <- 0
colnames(simMat) <- expr[,1]
rownames(simMat) <- colnames(simMat)
save(simMat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_combined_final.Rdata"))

# 1/22/2026: I'm not sure where 1_10-48-30 went, so I redid the refiltering step.
# To my suprise, 127 samples (of 1518) were filtered.
# New data at 1_10-13-40


##################################
######## try while removing Brainseq/NABEC
###################################

# Concatenate and process all datasets together
library(AnnotationHub)


dat_names <- c("GTEx",
               "ROSMAP",
               "BrainGVEX",
               "CMC",
               "CMC_HBCC")

exprs <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/GTEx_samp_filt_SampleNetworks/1_09-10-16/GTEx_samp_filt_1_190_outliers_removed.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_SampleNetworks/1_08-59-36/ROSMAP_samp_filt_1_617_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/BrainGVEX_samp_filt_SampleNetworks/1_08-25-36/BrainGVEX_samp_filt_1_275_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_samp_filt_SampleNetworks/1_11-19-13/CMC_samp_filt_1_285_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_HBCC_samp_filt_SampleNetworks/1_10-40-08/CMC_HBCC_samp_filt_1_162_ComBat.csv"))
expr_list <- list()
for(i in 1:length(exprs)){
  expr_list[[i]] <- fread(exprs[i], data.table=F)
}

# Find common features
for(i in 1:length(expr_list)){
  if(i==1){
    common_features <- expr_list[[i]][,1]
  } else {
    common_features <- intersect(common_features, expr_list[[i]][,1])
  }
}
common_features <- sort(common_features)

# Filter expr matrices to common features
for(i in 1:length(expr_list)){
  expr_list[[i]] <- expr_list[[i]][expr_list[[i]][,1] %in% common_features,]
}

# Order expr matrices by feature
for(i in 1:length(expr_list)){
  expr_list[[i]] <- expr_list[[i]][match(common_features,expr_list[[1]][,1]),]
}

# Combine expr matrices and create sample info mat
for(i in 1:length(expr_list)){
  if(i==1){
    mega_expr <- expr_list[[i]]
    mega_sif <- data.frame("Dataset" = dat_names[i],
                           "samp_id" = colnames(expr_list[[i]])[2:ncol(expr_list[[i]])])
  } else {
    mega_expr <- cbind(mega_expr, expr_list[[i]][,2:ncol(expr_list[[i]])])
    mega_sif <- rbind(mega_sif, data.frame("Dataset" = dat_names[i],
                                           "samp_id" = colnames(expr_list[[i]])[2:ncol(expr_list[[i]])]))
  }
}

colnames(mega_expr)[1] <- "gene_id"

# Add gene symbol to expr matrix
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, c("GENEID", "GENENAME"), "TXNAME")

genes <- ahEdbgene[ahEdbgene[,2] %in% mega_expr[,1],c(2,3)]
genes <- genes[!duplicated(genes[,1]),]
# Many ENSG ids have no associated gene symbol; remove these
genes <- genes[genes[,2] != "",]
genes <- genes[!duplicated(genes[,2]),]
genes <- genes[order(genes[,1]),]

mega_expr <- mega_expr[mega_expr[,1] %in% genes[,1],]
mega_expr <- mega_expr[match(genes[,1], mega_expr[,1]),]
mega_expr <- cbind(genes, mega_expr[,2:ncol(mega_expr)])
colnames(mega_expr)[1:2] <- c("ensembl_id", "Gene")

fwrite(mega_expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC.csv"))
fwrite(mega_sif, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined_minus_Brainseq_NABEC.csv"))

# Scale and save
si <- 3:ncol(mega_expr)
mega_expr[,si] <- mega_expr[,si]+abs(min(mega_expr[,si],na.rm=T))
fwrite(mega_expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_scaled.csv"))

source(file.path(Sys.getenv("SAMPLENETWORK_DIR", "/home/gugene/code/labcode_old/SampleNetwork"), "SampleNetwork_1.08.r"))

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/"))
#mega_expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC.csv"), data.table=F)
mega_expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_scaled.csv"), data.table=F)
mega_sif <- fread(data.table=F, file= file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined_minus_Brainseq_NABEC.csv"))
mega_sif$grouplabels1 <- 1

SampleNetwork(
  datExprT=mega_expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=2,
  indices1=list(seq(3,ncol(mega_expr))),
  subgroup1=1, # color samples by Dataset
  sampleinfo1=mega_sif,
  samplelabels1=2,
  grouplabels1=which(colnames(mega_sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(1), # Dataset
  trait1=NULL,
  asfactors1=c(1),
  projectname1="combined_FCX_minus_Brainseq_NABEC_scaled",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

library(edgeR)

# apply normalization
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat.csv"), data.table=F)
tmmf <- calcNormFactors(expr[,3:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 3:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
}
colnames(expr2)[1] <- "gene_id"

fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_ctf_asinh.csv"))

# scale first then apply normalization
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat.csv"), data.table=F)
si <- 3:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
tmmf <- calcNormFactors(expr[,si], method = "TMM")
expr2 <- expr
for(i in 3:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
}
colnames(expr2)[1] <- "gene_id"
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_scale_ctf_asinh.csv"))

expr3 <- expr
for(i in 3:ncol(expr)){
  expr3[,i] <- expr3[,i]/tmmf[i-2]
}
colnames(expr3)[1] <- "gene_id"
fwrite(expr3, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_scale_ctf.csv"))

# scale -> filt -> combat -> asinh
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat.csv"), data.table=F)
expr3 <- expr
for(i in 3:ncol(expr)){
  expr3[,i] <- asinh(expr3[,i])
}
colnames(expr3)[1] <- "gene_id"
fwrite(expr3, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat_asinh.csv"))


# scale -> filt -> combat -> ctf -> asinh
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat.csv"), data.table=F)
si <- 3:ncol(expr)
tmmf <- calcNormFactors(expr[,si], method = "TMM")
expr2 <- expr
for(i in 3:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
}
colnames(expr2)[1] <- "gene_id"
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat_ctf_asinh.csv"))

# scale -> filt -> combat -> scale -> ctf -> asinh
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat.csv"), data.table=F)
si <- 3:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
tmmf <- calcNormFactors(expr[,si], method = "TMM")
expr2 <- expr
for(i in 3:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
}
colnames(expr2)[1] <- "gene_id"
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat_scale_ctf_asinh.csv"))

# filt -> combat -> ctf -> scale just for kicks
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat.csv"), data.table=F)
tmmf <- calcNormFactors(expr[,3:ncol(expr)], method = "TMM")
expr2 <- expr
for(i in 3:ncol(expr)){
  expr2[,i] <- expr2[,i]/tmmf[i-2]
}
colnames(expr2)[1] <- "gene_id"
si <- 3:ncol(expr2)
expr2[,si] <- expr2[,si]+abs(min(expr2[,si],na.rm=T))

fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_ctf_scale.csv"))


# Filter samples
mega_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_BrainGVEX_SampleNetworks/1_11-32-03/combined_FCX_minus_Brainseq_BrainGVEX_1_1323_ComBat_ctf_asinh.csv"), data.table=F)

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/"))

SampleNetwork(
  datExprT=mega_expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=2,
  indices1=list(seq(3,ncol(mega_expr))),
  subgroup1=1, # color samples by Dataset
  sampleinfo1=mega_sif,
  samplelabels1=2,
  grouplabels1=which(colnames(mega_sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(1), # Dataset
  trait1=NULL,
  asfactors1=c(1),
  projectname1="combined_FCX_minus_Brainseq_BrainGVEX",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

# calculate and save simMat

library(WGCNA)
#expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_BrainGVEX_SampleNetworks/1_09-32-32/combined_FCX_minus_Brainseq_BrainGVEX_1_1279_outliers_removed.csv"), data.table=F)
#sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined_minus_Brainseq_BrainGVEX.csv"), data.table=F)
# expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_11-42-16/combined_FCX_final_1_1788_ComBat.csv"), data.table=F)
# sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined.csv"), data.table=F)
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_BrainGVEX_SampleNetworks/1_11-32-03/combined_FCX_minus_Brainseq_BrainGVEX_1_1323_ComBat.csv"), data.table=F)
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined_minus_Brainseq_BrainGVEX.csv"), data.table=F)



sif <- sif[sif$samp_id %in% colnames(expr),]
sif[,1] <- as.factor(sif[,1])

simMat <- bicor(t(expr[,3:ncol(expr)]), use="p")
diag(simMat) <- 0
colnames(simMat) <- expr[,1]
rownames(simMat) <- colnames(simMat)
save(simMat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_combined_minus_Brainseq_BrainGVEX_no_norm.Rdata"))



##############################
######## try while removing Brainseq/NABEC and starting with normalized data
#### individual dataset: filt/combat by library batch/scale/ctf
#### then cat, filt, combat by dataset, asinh
############################

# Concatenate and process all datasets together
library(AnnotationHub)


dat_names <- c("GTEx",
               "ROSMAP",
               "BrainGVEX",
               "CMC",
               "CMC_HBCC")

exprs <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_GTEx_scale_ctf.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_ROSMAP_scale_ctf.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_BrainGVEX_scale_ctf.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_CMC_scale_ctf.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_CMC_HBCC_scale_ctf.csv"))
expr_list <- list()
for(i in 1:length(exprs)){
  expr_list[[i]] <- fread(exprs[i], data.table=F)
}

# Find common features
for(i in 1:length(expr_list)){
  if(i==1){
    common_features <- expr_list[[i]][,1]
  } else {
    common_features <- intersect(common_features, expr_list[[i]][,1])
  }
}
common_features <- sort(common_features)

# Filter expr matrices to common features
for(i in 1:length(expr_list)){
  expr_list[[i]] <- expr_list[[i]][expr_list[[i]][,1] %in% common_features,]
}

# Order expr matrices by feature
for(i in 1:length(expr_list)){
  expr_list[[i]] <- expr_list[[i]][match(common_features,expr_list[[1]][,1]),]
}

# Combine expr matrices and create sample info mat
for(i in 1:length(expr_list)){
  if(i==1){
    mega_expr <- expr_list[[i]]
    mega_sif <- data.frame("Dataset" = dat_names[i],
                           "samp_id" = colnames(expr_list[[i]])[2:ncol(expr_list[[i]])])
  } else {
    mega_expr <- cbind(mega_expr, expr_list[[i]][,2:ncol(expr_list[[i]])])
    mega_sif <- rbind(mega_sif, data.frame("Dataset" = dat_names[i],
                                           "samp_id" = colnames(expr_list[[i]])[2:ncol(expr_list[[i]])]))
  }
}

colnames(mega_expr)[1] <- "gene_id"

# Add gene symbol to expr matrix
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, c("GENEID", "GENENAME"), "TXNAME")

genes <- ahEdbgene[ahEdbgene[,2] %in% mega_expr[,1],c(2,3)]
genes <- genes[!duplicated(genes[,1]),]
# Many ENSG ids have no associated gene symbol; remove these
genes <- genes[genes[,2] != "",]
genes <- genes[!duplicated(genes[,2]),]
genes <- genes[order(genes[,1]),]

mega_expr <- mega_expr[mega_expr[,1] %in% genes[,1],]
mega_expr <- mega_expr[match(genes[,1], mega_expr[,1]),]
mega_expr <- cbind(genes, mega_expr[,2:ncol(mega_expr)])
colnames(mega_expr)[1:2] <- c("ensembl_id", "Gene")

fwrite(mega_expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf.csv"))
fwrite(mega_sif, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined_minus_Brainseq_NABEC_firstScaledCtf.csv"))

# Save scaled and ctf 
expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf.csv"), data.table=F)
expr4 <- expr
si <- 3:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_scaled.csv"))
tmmf <- calcNormFactors(expr[,si], method = "TMM")
tmmf2 <- calcNormFactors(expr4[,si], method = "TMM")
expr2 <- expr
expr3 <- expr
for(i in 3:ncol(expr2)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
  expr3[,i] <- expr3[,i]/tmmf[i-2]
  expr4[,i] <- expr4[,i]/tmmf2[i-2]
}
colnames(expr2)[1] <- "gene_id"
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_scaled_ctf_asinh.csv"))
colnames(expr3)[1] <- "gene_id"
fwrite(expr3, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_scaled_ctf.csv"))
colnames(expr4)[1] <- "gene_id"
fwrite(expr4, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_ctf.csv"))


# Perform Combat by dataset
source(file.path(Sys.getenv("SAMPLENETWORK_DIR", "/home/gugene/code/labcode_old/SampleNetwork"), "SampleNetwork_1.08.r"))

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/"))
mega_expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf.csv"), data.table=F)
mega_sif$grouplabels1 <- 1

SampleNetwork(
  datExprT=mega_expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=2,
  indices1=list(seq(3,ncol(mega_expr))),
  subgroup1=1, # color samples by Dataset
  sampleinfo1=mega_sif,
  samplelabels1=2,
  grouplabels1=which(colnames(mega_sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(1), # Dataset
  trait1=NULL,
  asfactors1=c(1),
  projectname1="combined_FCX_minus_Brainseq_NABEC_firstScaledCtf",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

library(edgeR)

# apply asinh 
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat.csv"), data.table=F)
#tmmf <- calcNormFactors(expr[,3:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 3:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i])
}
colnames(expr2)[1] <- "gene_id"

fwrite(expr2, file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_asinh.csv"))

# scale first then apply normalization
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat.csv"), data.table=F)
si <- 3:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale.csv"))

tmmf <- calcNormFactors(expr[,si], method = "TMM")
expr2 <- expr
for(i in 3:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-2])
}
colnames(expr2)[1] <- "gene_id"
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale_ctf_asinh.csv"))

expr3 <- expr
for(i in 3:ncol(expr)){
  expr3[,i] <- expr3[,i]/tmmf[i-2]
}
colnames(expr3)[1] <- "gene_id"
fwrite(expr3, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale_ctf.csv"))

expr4 <- expr
for(i in 3:ncol(expr)){
  expr4[,i] <- asinh(expr4[,i])
}
colnames(expr4)[1] <- "gene_id"
fwrite(expr4, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale_asinh.csv"))


# Filter samples
mega_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_BrainGVEX_SampleNetworks/1_11-32-03/combined_FCX_minus_Brainseq_BrainGVEX_1_1323_ComBat_ctf_asinh.csv"), data.table=F)

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/"))

SampleNetwork(
  datExprT=mega_expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=2,
  indices1=list(seq(3,ncol(mega_expr))),
  subgroup1=1, # color samples by Dataset
  sampleinfo1=mega_sif,
  samplelabels1=2,
  grouplabels1=which(colnames(mega_sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(1), # Dataset
  trait1=NULL,
  asfactors1=c(1),
  projectname1="combined_FCX_minus_Brainseq_BrainGVEX",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

# calculate and save simMat

library(WGCNA)
#expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_BrainGVEX_SampleNetworks/1_09-32-32/combined_FCX_minus_Brainseq_BrainGVEX_1_1279_outliers_removed.csv"), data.table=F)
#sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined_minus_Brainseq_BrainGVEX.csv"), data.table=F)
# expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_11-42-16/combined_FCX_final_1_1788_ComBat.csv"), data.table=F)
# sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined.csv"), data.table=F)
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_BrainGVEX_SampleNetworks/1_11-32-03/combined_FCX_minus_Brainseq_BrainGVEX_1_1323_ComBat.csv"), data.table=F)
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/sif_combined_minus_Brainseq_BrainGVEX.csv"), data.table=F)



sif <- sif[sif$samp_id %in% colnames(expr),]
sif[,1] <- as.factor(sif[,1])

simMat <- bicor(t(expr[,3:ncol(expr)]), use="p")
diag(simMat) <- 0
colnames(simMat) <- expr[,1]
rownames(simMat) <- colnames(simMat)
save(simMat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_combined_minus_Brainseq_BrainGVEX_no_norm.Rdata"))



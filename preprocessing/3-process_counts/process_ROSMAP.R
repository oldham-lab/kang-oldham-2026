library(data.table)
library(tximport)
library(AnnotationHub)


#dirs <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/kal_quant/"), full.names=T)
# 
# genes <- c()
# tpm_mat <- list()
# count_mat <- list()
# for(i in 1:length(dirs)){
#   temp <- fread(paste0(dirs[i], "/abundance.tsv"))
#   if(i==1){
#     genes <- temp[,1]
#   }
#   
#   tpm_mat[[i]] <- temp[,5]
#   count_mat[[i]] <- temp[,4]
#   cat(i, "done\n")
# }
# 
# tpm_mat <- as.data.frame(do.call(cbind, tpm_mat))
# tpm_mat <- cbind(genes, tpm_mat)
# count_mat <- as.data.frame(do.call(cbind, count_mat))
# count_mat <- cbind(genes, count_mat)
# 
# # gather sample names
# samp_names <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/kal_quant/"), full.names=F)
# samp_names <- gsub(".trimmed_transcripts_quant_kal", "", samp_names)
# 
# colnames(tpm_mat)[2:ncol(tpm_mat)] <- samp_names
# colnames(count_mat)[2:ncol(count_mat)] <- samp_names
# 
# temp <- sif2[1:197,]
# 
# dir.create(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/kal_quant/count_mats/"))
# fwrite(tpm_mat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/kal_quant/count_mats/rosmap_tpm.csv"))
# fwrite(count_mat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/kal_quant/count_mats/rosmap_raw_counts.csv"))

# retrieve latest Ensembl annotation
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
#gns <- genes(ahEdb)
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, "GENEID", "TXNAME")

dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/fastq/trimmed/kal_quant/")
dirs <- list.files(dir, full.names=T)
dirs <- dirs[-grep("count_mats", dirs)]
files <- file.path(dirs, "abundance.h5")

txi.kal <- tximport(files, type = "kallisto", tx2gene = ahEdbgene, ignoreTxVersion = T)

prot_gene_ID <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/3-process_counts/protein_coding_gene_IDs.csv"), header=F, data.table=F)

# Save raw counts
expr <- txi.kal[[2]]
expr <- expr[rownames(expr) %in% prot_gene_ID[,1],]
samp_names <- list.files(dir)
samp_names <- samp_names[grep(".trimmed_transcripts",samp_names)]
samp_names <- gsub(".trimmed_transcripts_quant_kal", "", samp_names)
colnames(expr) <- samp_names
expr <- cbind(rownames(expr), data.frame(expr))
colnames(expr)[1] <- "gene_id"
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_ROSMAP_raw_counts.csv"))

## SampleNetwork


library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))
source(file.path(Sys.getenv("SAMPLENETWORK_DIR", "/home/gugene/code/labcode_old/SampleNetwork"), "SampleNetwork_1.08.r"))


expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_ROSMAP_raw_counts.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_ROSMAP.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
sif <- sif[sif$excludeReason=="",]
# filter expr to samples in sif
samp_names <- gsub("X", "", colnames(expr))
expr <- expr[,c(1, which(samp_names %in% sif$specimenID))]
samp_names <- samp_names[samp_names %in% sif$specimenID]
# order sif to match expr
sif <- sif[sif$specimenID %in% samp_names,]
sif <- sif[match(samp_names, sif$specimenID),]
sif[,2] <- colnames(expr)[2:ncol(expr)]
# remove sample 319 (unique library batch)
expr <- expr[-which(colnames(expr)==sif[319,2])]
sif <- sif[-319,]

sif[,24] <- as.factor(sif[,24])

# Remove AD samples
sif_clin <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_clinical.csv"))
sif <- left_join(sif,sif_clin, by=join_by("individualID"))
sif <- sif[sif$cogdx %in% c(1,2),]
expr <- expr[,c(1,which(colnames(expr) %in% sif$specimenID))]

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=24, # color samples by libraryBatch
  sampleinfo1=sif,
  samplelabels1=2,
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(24), #libraryBatch
  trait1=NULL,
  asfactors1=c(24), #libraryBatch
  projectname1="ROSMAP_samp_filt",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# For some reason SampleNetwork fails if I try to do sample filtering and Combat.
# Do sample filtering only first, then run SN again for ComBat only

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_SampleNetworks/1_07-59-12/ROSMAP_samp_filt_1_347_outliers_removed.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- sif[sif[,2] %in% colnames(expr),]
# run above SN code again



# apply CTF, asinh normalization to columns
library(edgeR)

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_SampleNetworks/1_08-13-25/ROSMAP_samp_filt_1_347_ComBat.csv"), data.table=F)
si <- 2:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
tmmf <- calcNormFactors(expr[,2:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 2:ncol(expr)){
  #expr2[,i] <- asinh(expr2[,i]/tmmf[i-1])
  expr2[,i] <- expr2[,i]/tmmf[i-1]
}
colnames(expr2)[1] <- "gene_id"

#fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_ROSMAP_ctf_asinh.csv"))
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_ROSMAP_scale_ctf.csv"))

##########################
# Process AD samples
##########################

library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_ROSMAP_raw_counts.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_ROSMAP.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
sif <- sif[sif$excludeReason=="",]
# filter expr to samples in sif
samp_names <- gsub("X", "", colnames(expr))
expr <- expr[,c(1, which(samp_names %in% sif$specimenID))]
samp_names <- samp_names[samp_names %in% sif$specimenID]
# order sif to match expr
sif <- sif[sif$specimenID %in% samp_names,]
sif <- sif[match(samp_names, sif$specimenID),]
sif[,2] <- colnames(expr)[2:ncol(expr)]
# remove sample 319 (unique library batch)
expr <- expr[-which(colnames(expr)==sif[319,2])]
sif <- sif[-319,]
sif[,24] <- as.factor(sif[,24])

# Keep only AD samples
sif_clin <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_clinical.csv"))
sif <- left_join(sif,sif_clin, by=join_by("individualID"))
sif <- sif[sif$cogdx %in% c(4,5),]
expr <- expr[,c(1,which(colnames(expr) %in% sif$specimenID))]

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=24, # color samples by libraryBatch
  sampleinfo1=sif,
  samplelabels1=2,
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(24), #libraryBatch
  trait1=NULL,
  asfactors1=c(24), #libraryBatch
  projectname1="ROSMAP_samp_filt_AD_only",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# For some reason SampleNetwork fails if I try to do sample filtering and Combat.
# Do sample filtering only first, then run SN again for ComBat only

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_AD_only_SampleNetworks/1_10-49-39/ROSMAP_samp_filt_AD_only_1_248_outliers_removed.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- sif[sif[,2] %in% colnames(expr),]
# run above SN code again

#####################
# Process all samples
#####################

library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_ROSMAP_raw_counts.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_ROSMAP.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
sif <- sif[sif$excludeReason=="",]
# filter expr to samples in sif
samp_names <- gsub("X", "", colnames(expr))
expr <- expr[,c(1, which(samp_names %in% sif$specimenID))]
samp_names <- samp_names[samp_names %in% sif$specimenID]
# order sif to match expr
sif <- sif[sif$specimenID %in% samp_names,]
sif <- sif[match(samp_names, sif$specimenID),]
sif[,2] <- colnames(expr)[2:ncol(expr)]
# remove sample 319 (unique library batch)
expr <- expr[-which(colnames(expr)==sif[319,2])]
sif <- sif[-319,]
sif[,24] <- as.factor(sif[,24])

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=24, # color samples by libraryBatch
  sampleinfo1=sif,
  samplelabels1=2,
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(24), #libraryBatch
  trait1=NULL,
  asfactors1=c(24), #libraryBatch
  projectname1="ROSMAP_samp_filt_all_samples",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# For some reason SampleNetwork fails if I try to do sample filtering and Combat.
# Do sample filtering only first, then run SN again for ComBat only

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_all_samples_SampleNetworks/1_10-45-09/ROSMAP_samp_filt_all_samples_1_617_outliers_removed.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- sif[sif[,2] %in% colnames(expr),]
# run above SN code again


##########################
# Process CON samples only
##########################

library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_ROSMAP_raw_counts.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_ROSMAP.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
sif <- sif[sif$excludeReason=="",]
# filter expr to samples in sif
samp_names <- gsub("X", "", colnames(expr))
expr <- expr[,c(1, which(samp_names %in% sif$specimenID))]
samp_names <- samp_names[samp_names %in% sif$specimenID]
# order sif to match expr
sif <- sif[sif$specimenID %in% samp_names,]
sif <- sif[match(samp_names, sif$specimenID),]
sif[,2] <- colnames(expr)[2:ncol(expr)]
# remove sample 319 (unique library batch)
expr <- expr[-which(colnames(expr)==sif[319,2])]
sif <- sif[-319,]
sif[,24] <- as.factor(sif[,24])

# Keep only con samples
sif_clin <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_clinical.csv"))
sif <- left_join(sif,sif_clin, by=join_by("individualID"))
sif <- sif[!sif$cogdx %in% c(4,5),]
expr <- expr[,c(1,which(colnames(expr) %in% sif$specimenID))]

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=24, # color samples by libraryBatch
  sampleinfo1=sif,
  samplelabels1=2,
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(24), #libraryBatch
  trait1=NULL,
  asfactors1=c(24), #libraryBatch
  projectname1="ROSMAP_samp_filt_con_only",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# For some reason SampleNetwork fails if I try to do sample filtering and Combat.
# Do sample filtering only first, then run SN again for ComBat only

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_con_only_SampleNetworks/1_01-03-44/ROSMAP_samp_filt_con_only_1_369_outliers_removed.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- sif[sif[,2] %in% colnames(expr),]
# run above SN code again

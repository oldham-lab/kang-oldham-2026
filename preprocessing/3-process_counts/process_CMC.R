library(data.table)
library(tximport)
library(AnnotationHub)


# retrieve latest Ensembl annotation
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
#gns <- genes(ahEdb)
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, "GENEID", "TXNAME")

dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/RNAseq/fastq/kal_quant/")
dirs <- list.files(dir, full.names=T)
files <- file.path(dirs, "abundance.h5")

txi.kal <- tximport(files, type = "kallisto", tx2gene = ahEdbgene, ignoreTxVersion = T)

prot_gene_ID <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/3-process_counts/protein_coding_gene_IDs.csv"), header=F, data.table=F)

# Save raw counts
expr <- txi.kal[[2]]
expr <- expr[rownames(expr) %in% prot_gene_ID[,1],]
samp_names <- list.files(dir)
samp_names <- gsub("_trimmed.gz_transcripts_quant_kal", "", samp_names)
colnames(expr) <- samp_names
expr <- cbind(rownames(expr), data.frame(expr))
colnames(expr)[1] <- "gene_id"
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_CMC_raw_counts.csv"))

## SampleNetwork


library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))
source("/home/gugene/code/SampleNetwork/SampleNetwork_1.08.r")


expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_CMC_raw_counts.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_CMC.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable

temp <- colnames(expr)[2:ncol(expr)]
temp <- gsub(".fastq_trimmed_transcripts_quant_kal", "", temp)
sif <- sif[match(temp, sif[,2]),]
sif$Library_Batch <- as.factor(sif$Library_Batch)
colnames(expr)[2:ncol(expr)] <- sif$SampleID

# Split between CMC and CMC_HBCC

sif1 <- sif[sif$Study=="CMC",]
sif2 <- sif[sif$Study=="CMC_HBCC",]
expr1 <- expr[,c(1, which(temp %in% sif1$SampleID)+1)]
expr2 <- expr[,c(1, which(temp %in% sif2$SampleID)+1)]
fwrite(expr1, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_CMC_only_raw_counts.csv"))
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_CMC_HBCC_raw_counts.csv"))

# Z.K < -4

# CMC

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr1,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr1))),
  subgroup1=which(colnames(sif1)=="Library_Batch"), # color samples by libraryBatch
  sampleinfo1=sif1,
  samplelabels1=which(colnames(sif1)=="SampleID"),
  grouplabels1=which(colnames(sif1)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(which(colnames(sif1)=="Library_Batch")), #libraryBatch
  trait1=NULL,
  asfactors1=c(which(colnames(sif1)=="Library_Batch")), #libraryBatch
  projectname1="CMC_samp_filt",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# For some reason SampleNetwork fails if I try to do sample filtering and Combat.
# Do sample filtering only first, then run SN again for ComBat only

expr1 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_samp_filt_SampleNetworks/1_11-13-38/CMC_samp_filt_1_285_outliers_removed.csv"), data.table=F)
colnames(expr1)[1] <- "Gene"
sif1 <- sif1[sif1$SampleID %in% colnames(expr1),]
# run above SN code again


# CMC HBCC

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr2,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr2))),
  subgroup1=which(colnames(sif2)=="Library_Batch"), # color samples by libraryBatch
  sampleinfo1=sif2,
  samplelabels1=which(colnames(sif2)=="SampleID"),
  grouplabels1=which(colnames(sif2)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(which(colnames(sif2)=="Library_Batch")), #libraryBatch
  trait1=NULL,
  asfactors1=c(which(colnames(sif2)=="Library_Batch")), #libraryBatch
  projectname1="CMC_HBCC_samp_filt",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# For some reason SampleNetwork fails if I try to do sample filtering and Combat.
# Do sample filtering only first, then run SN again for ComBat only

expr2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_HBCC_samp_filt_SampleNetworks/1_10-36-40/CMC_HBCC_samp_filt_1_162_outliers_removed.csv"), data.table=F)
colnames(expr2)[1] <- "Gene"
sif2 <- sif2[sif2$SampleID %in% colnames(expr2),]
# run above SN code again


# apply CTF, asinh normalization to columns
library(edgeR)

# CMC
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_samp_filt_SampleNetworks/1_11-19-13/CMC_samp_filt_1_285_ComBat.csv"), data.table=F)
si <- 2:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
tmmf <- calcNormFactors(expr[,2:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 2:ncol(expr)){
  #expr2[,i] <- asinh(expr2[,i]/tmmf[i-1])
  expr2[,i] <- expr2[,i]/tmmf[i-1]
}
colnames(expr2)[1] <- "gene_id"

#fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_CMC_ctf_asinh.csv"))
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_CMC_scale_ctf.csv"))


# CMC HBCC
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_HBCC_samp_filt_SampleNetworks/1_10-40-08/CMC_HBCC_samp_filt_1_162_ComBat.csv"), data.table=F)
si <- 2:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
tmmf <- calcNormFactors(expr[,2:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 2:ncol(expr)){
  #expr2[,i] <- asinh(expr2[,i]/tmmf[i-1])
  expr2[,i] <- expr2[,i]/tmmf[i-1]
}
colnames(expr2)[1] <- "gene_id"

#fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_CMC_HBCC_ctf_asinh.csv"))
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_CMC_HBCC_scale_ctf.csv"))

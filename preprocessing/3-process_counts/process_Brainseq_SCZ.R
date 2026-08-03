# Format quantified reads into counts matrix (Brainseq)

### SampleNetwork
library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_Brainseq_raw_counts.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_Brainseq_SCZ.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
# filter expr to samples in sif
expr <- expr[,c(1, which(colnames(expr) %in% sif[,1]))]
# order sif to match expr
sif <- sif[match(colnames(expr)[2:ncol(expr)], sif[,1]),]

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=1, # color samples by donor_name (individual)
  sampleinfo1=sif,
  samplelabels1=1,
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(5,6,7), # Age, Sex, Race
  trait1=NULL,
  asfactors1=c(6,7),
  projectname1="brainseq_samp_filt_SCZ",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

# Add gene symbols
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SCZ_SampleNetworks/1_04-04-29/brainseq_samp_filt_SCZ_1_171_outliers_removed.csv"), data.table=F)

library(AnnotationHub)
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, c("GENEID", "GENENAME"), "TXNAME")

genes <- ahEdbgene[ahEdbgene[,2] %in% expr[,1], c(2,3)]
genes <- genes[!duplicated(genes[,1]),]
# Many ENSG ids have no associated gene symbol; remove these
genes <- genes[genes[,2] != "",]
genes <- genes[!duplicated(genes[,2]),]
genes <- genes[order(genes[,1]),]

expr <- expr[expr[,1] %in% genes[,1],]
expr <- expr[match(genes[,1], expr[,1]),]
expr <- cbind(genes, expr[,2:ncol(expr)])
colnames(expr)[1:2] <- c("ensembl_id", "Gene")
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SCZ_SampleNetworks/1_04-04-29/brainseq_samp_filt_SCZ_1_171_outliers_removed_geneSymbolsAdded.csv"))


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

dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/GTEx/fastq/kal_quant/")
dirs <- list.files(dir, full.names=T)
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
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_GTEx_raw_counts.csv"))


## SampleNetwork

library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_GTEx_raw_counts.csv"), data.table=F,)
colnames(expr)[1] <- "Gene"

sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_GTEx.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
sif$specimen_id <- gsub("-", ".", sif$specimen_id)

# order sif to match expr
sif <- sif[sif$specimen_id %in% colnames(expr),]
sif <- sif[match(colnames(expr)[2:ncol(expr)], sif$specimen_id),]

# # Z.K < -4

# setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))
# 
# SampleNetwork(
#   datExprT=expr,
#   method1="correlation",
#   impute1=FALSE,
#   subset1=NULL,
#   skip1=1,
#   indices1=list(seq(2,ncol(expr))),
#   subgroup1=8, # color samples by collection site
#   sampleinfo1=sif,
#   samplelabels1=18, # specimen_id
#   grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
#   fitmodels1=TRUE,
#   whichmodel1="univariate",
#   whichfit1="pc1",
#   btrait1=c(8), # collection site
#   trait1=NULL,
#   asfactors1=c(8), # collection site
#   projectname1="GTEx_samp_filt",
#   cexlabels=0.7,
#   normalize1=F,
#   replacenegs1=FALSE,
#   exportfigures1=TRUE,
#   verbose=TRUE
# )

## correct expr matrix according to library batch

# there are few samples from some of the collection sites so remove those
sif <- sif[sif$bss_collection_site %in% c("B1, A1", "C1, A1"),]
expr <- expr[, c(1, which(colnames(expr) %in% sif$specimen_id))]

sif[,45] <- as.factor(sif[,45])

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=45, # color samples by library batch (date_nucleic_acid_isolation)
  sampleinfo1=sif,
  samplelabels1=18, # specimen_id
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(45), # date_nucleic_acid_isolation
  trait1=NULL,
  asfactors1=c(45), # date_nucleic_acid_isolation
  projectname1="GTEx_samp_filt_combat_library_date_batch",
  cexlabels1=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# # run samplenetwork to remove outliers then run for combat
# expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/GTEx_samp_filt_combat_library_date_batch_SampleNetworks/1_03-33-03/GTEx_samp_filt_combat_library_date_batch_1_183_outliers_removed.csv"), data.table=F)
# colnames(expr)[1] <- "Gene"
# sif <- sif[sif$specimen_id %in% colnames(expr),]

# # apply CTF, asinh normalization to columns
# library(edgeR)

# expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/GTEx_samp_filt_combat_library_date_batch_SampleNetworks/1_03-44-38/GTEx_samp_filt_combat_library_date_batch_1_183_ComBat.csv"), data.table=F)
# si <- 2:ncol(expr)
# expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
# tmmf <- calcNormFactors(expr[,2:ncol(expr)], method = "TMM")

# expr2 <- expr
# for(i in 2:ncol(expr)){
#   #expr2[,i] <- asinh(expr2[,i]/tmmf[i-1])
#   expr2[,i] <- expr2[,i]/tmmf[i-1]
# }
# colnames(expr2)[1] <- "gene_id"

# #fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_GTEx_ctf_asinh.csv"))
# fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_GTEx_scale_ctf.csv"))


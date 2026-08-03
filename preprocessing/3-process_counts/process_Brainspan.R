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

dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainspan/fastq/trimmed/kal_quant/")
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
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_Brainspan_raw_counts.csv"))


# SampleNetwork
library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_Brainspan_raw_counts.csv"), data.table=F,)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_Brainspan.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
# filter expr to samples in sif
expr <- expr[,c(1, which(colnames(expr) %in% sif$Run))]
# order sif to match expr
sif <- sif[match(colnames(expr)[2:ncol(expr)], sif$Run),]

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=1, # color samples by sample name
  sampleinfo1=sif,
  samplelabels1=1,
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(32), # sex
  trait1=NULL,
  asfactors1=c(32), # sex
  projectname1="Brainspan_samp_filt",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

# apply CTF, asinh normalization to columns
library(edgeR)

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/Brainspan_samp_filt_SampleNetworks/1_11-32-00/Brainspan_samp_filt_1_40_outliers_removed.csv"), data.table=F)
tmmf <- calcNormFactors(expr[,2:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 2:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-1])
}
colnames(expr2)[1] <- "gene_id"

fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_Brainspan_ctf_asinh.csv"))

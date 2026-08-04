library(data.table)
library(tximport)
library(AnnotationHub)
library(gridExtra)


# retrieve latest Ensembl annotation
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
#gns <- genes(ahEdb)
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, "GENEID", "TXNAME")

dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/RNAseq/BrainGVEX/trimmed/kal_quant/")
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
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_BrainGVEX_raw_counts.csv"))

## SampleNetwork


library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))
source(file.path(Sys.getenv("SAMPLENETWORK_DIR", "/home/gugene/code/labcode_old/SampleNetwork"), "SampleNetwork_1.08.r"))


expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_BrainGVEX_raw_counts.csv"), data.table=F,)
colnames(expr)[1] <- "Gene"
sif <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_BrainGVEX.csv"))
#sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_BrainGVEX.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
# match colnames of expr and sif samp names
sif$RNAseq_filename <- gsub("-", ".",sif$RNAseq_filename)
temp <- lapply(colnames(expr)[2:ncol(expr)], function(x) strsplit(x, "_")[[1]][1:7])
temp <- unlist(lapply(temp, function(x) paste(x, collapse="_")))
sif <- sif[match(temp, sif$RNAseq_filename),]
sif$RNAseq_filename <- colnames(expr)[2:ncol(expr)]

sif$LibraryBatch <- as.factor(sif$LibraryBatch)

# > which(colnames(sif)=="LibraryBatch")
# [1] 9

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=which(colnames(sif)=="LibraryBatch"), # color samples by libraryBatch
  sampleinfo1=sif,
  samplelabels1=which(colnames(sif)=="RNAseq_filename"),
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(which(colnames(sif)=="LibraryBatch")), #libraryBatch
  trait1=NULL,
  asfactors1=c(which(colnames(sif)=="LibraryBatch")), #libraryBatch
  projectname1="BrainGVEX_samp_filt",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
# For some reason SampleNetwork fails if I try to do sample filtering and Combat.
# Do sample filtering only first, then run SN again for ComBat only

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/BrainGVEX_samp_filt_SampleNetworks/1_08-20-48/BrainGVEX_samp_filt_1_275_outliers_removed.csv"), data.table=F)
colnames(expr)[1] <- "Gene"
sif <- sif[sif$RNAseq_filename %in% colnames(expr),]
# run above SN code again



# apply CTF, asinh normalization to columns
library(edgeR)

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/BrainGVEX_samp_filt_SampleNetworks/1_08-25-36/BrainGVEX_samp_filt_1_275_ComBat.csv"), data.table=F)
#y <- DGEList(counts=expr[,2:ncol(expr)], genes=expr[,1])
si <- 2:ncol(expr)
expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
tmmf <- calcNormFactors(expr[,2:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 2:ncol(expr)){
  #expr2[,i] <- asinh(expr2[,i]/tmmf[i-1])
  expr2[,i] <- expr2[,i]/tmmf[i-1]
}
colnames(expr2)[1] <- "gene_id"

#fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_BrainGVEX_ctf_asinh.csv"))
fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_BrainGVEX_scale_ctf.csv"))







############
# test covariates



expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_BrainGVEX_raw_counts.csv"), data.table=F,)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_BrainGVEX.csv"), data.table=F)
sif$grouplabels1 <- 1 # dummy column for "grouplabels1" variable
# match colnames of expr and sif samp names
sif$RNAseq_filename <- gsub("-", ".",sif$RNAseq_filename)
temp <- lapply(colnames(expr)[2:ncol(expr)], function(x) strsplit(x, "_")[[1]][1:7])
temp <- unlist(lapply(temp, function(x) paste(x, collapse="_")))
sif <- sif[match(temp, sif$RNAseq_filename),]
sif$RNAseq_filename <- colnames(expr)[2:ncol(expr)]

sif$LibraryBatch <- as.factor(sif$LibraryBatch)

# > which(colnames(sif)=="LibraryBatch")
# [1] 9

# Z.K < -4

setwd(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/"))

these_indices <- c(1:72)

SampleNetwork(
  datExprT=expr,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(seq(2,ncol(expr))),
  subgroup1=which(colnames(sif)=="LibraryBatch"), # color samples by libraryBatch
  sampleinfo1=sif,
  samplelabels1=which(colnames(sif)=="RNAseq_filename"),
  grouplabels1=which(colnames(sif)=="grouplabels1"), # Group variable
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=c(these_indices), 
  trait1=NULL,
  asfactors1=c(these_indices), 
  projectname1="BrainGVEX_samp_filt_testing",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)


# calculate PCs and model each gene by PCs
expr_t <- t(expr[,2:ncol(expr)])
colnames(expr_t) <- expr[,1]
expr_t <- apply(expr_t,2,scale)
exclude_these <- which(is.na(apply(expr_t,2,mean)))
expr_t <- expr_t[,-exclude_these]
PCs <- prcomp(expr_t, scale=F)
pc_10 <- PCs$x[,1:10]

model_list <- list()
for(i in 1:10){
  vec <- c()
  for(j in 1:ncol(expr_t)){
    vec <- c(vec, summary(lm(expr_t[,j] ~ pc_10[,i]))$r.squared)
  }
  model_list[[i]] <- vec
}

# Graph modeling results

plot_df <- list()
for(i in 1:10){
  plot_df[[i]] <- data.frame("pc" = rep(i, length(model_list[[i]])),
                             "r2" = model_list[[i]])
}
plot_df <- do.call(rbind, plot_df)

plot_df[,1] <- factor(plot_df[,1], levels=1:10)

p1 <- ggplot(plot_df, aes(x = pc, y = r2)) + 
  geom_boxplot(notch=T, outlier.size=0.5) +
  ggtitle("Variance explained by PCs 1-10 in BrainGVEX data",
          sub = "19548 genes, 276 samples") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5)) +
  ylab(bquote(r^2)) +
  xlab("PC")

# plot cors of covariates with pcs
sif_pc_cors <- list()
for(i in 1:ncol(sif)){
  temp <- sif[,i]
  if(class(temp) %in% c("character", "factor")){
    temp <- as.numeric(as.factor(temp))
  }
  sif_pc_cors[[i]] <- cor(temp, pc_10)
}

plot_df2 <- list()
for(i in 1:length(sif_pc_cors)){
  plot_df2[[i]] <- data.frame("covariate" = colnames(sif)[i],
                         "cor" = sif_pc_cors[[i]][1,])
}
plot_df2 <- do.call(rbind, plot_df2)

p2 <- ggplot(plot_df2, aes(x = covariate, y = cor)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.4) +
  xlab("") +
  ylab("Pearson correlation") +
  ggtitle("Correlation of PCs 1-10 with study covariates",
          sub = "BrainGVEX data") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# transpose
sif_pc_cors <- list()
for(i in 1:ncol(pc_10)){
  temp_vec <- c()
  for(j in 1:ncol(sif)){
    temp <- sif[,j]
    if(class(temp) %in% c("character", "factor")){
      temp <- as.numeric(as.factor(temp))
    }
    temp_vec[j] <- cor(pc_10[,i], temp)
  }
  sif_pc_cors[[i]] <- temp_vec
}

plot_df2 <- list()
for(i in 1:length(sif_pc_cors)){
  plot_df2[[i]] <- data.frame("pc" = i,
                              "cor" = sif_pc_cors[[i]])
}
plot_df2 <- do.call(rbind, plot_df2)

plot_df2[,1] <- factor(plot_df2[,1], levels=1:10)

p3 <- ggplot(plot_df2, aes(x = pc, y = cor)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.4) +
  ggtitle("Correlation of study covariates with PCs 1-10",
          sub="BrainGVEX data") +
  xlab("PC") +
  ylab("Pearson correlation") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5)) 

lay = rbind(c(1,1),
            c(2,3))
p_all <- grid.arrange(grobs=list(p1,p2,p3), 
                      #nrow=3,
                      layout_matrix=lay
                      )
ggsave(p_all, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/4-analyses/BrainGVEX_PC_covariate_analysis.png"),
       width=15, height=8)

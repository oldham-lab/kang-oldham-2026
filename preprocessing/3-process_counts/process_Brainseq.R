# Format quantified reads into counts matrix (Brainseq)
library(tidyverse)
library(data.table)
library(tximport)


dirs <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainseq/Salmon_aligned/decoy/"), full.names=T)

samp_list <- list()
genes <- c()
for(i in 1:length(dirs)){
# not entirely sure what the difference is between quant.sf and quant.genes.sf besides order
  test <- fread(list.files(dirs[i], full.names=T)[[6]], data.table=F)
  
  if(i==1){
    genes <- test[,1]
  }
  ind <- unlist(lapply(test[,1], function(x) grepl("protein_coding", x)))
  test <- test[ind,]
  
  if(i==1){
    genes <- test[,1]
  }
  
  samp_list[[i]] <- test[,4]
  cat(i, "done\n")
}

expr <- as.data.frame(do.call(cbind, samp_list))
expr <- cbind(genes, expr)

# add sample names
samp <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainseq/Salmon_aligned/decoy/"), full.names=F)
samp <- gsub(".trimmed_transcripts_quant", "", samp)

colnames(expr)[1] <- "Genes"
colnames(expr)[2:ncol(expr)] <- samp

# write expr matrix (transcript-level)
fwrite(expr, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainseq/counts/expr_TPM.csv"))

# Alternative method: use tximport to concatenate data at gene-level
# retrieve latest Ensembl annotation
library(AnnotationHub)
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
#gns <- genes(ahEdb)
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, "GENEID", "TXNAME")

dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainseq/Salmon_aligned/decoy/")
files <- file.path(list.files(dir, full.names=T), "quant.sf")

txi.salmon <- tximport(files, type = "salmon", tx2gene = ahEdbgene, ignoreAfterBar = T, ignoreTxVersion = T)

# load table of protein coding gene IDS
prot_gene_ID <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/3-process_counts/protein_coding_gene_IDs.csv"), header=F, data.table=F)

# Save raw counts
expr <- txi.salmon[[2]]
expr <- expr[rownames(expr) %in% prot_gene_ID[,1],]
samp_names <- list.files(dir)
samp_names <- gsub(".trimmed_transcripts_quant", "", samp_names)
colnames(expr) <- samp_names
fwrite(data.frame(expr), file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_Brainseq_raw_counts.csv"), row.names=T)

# # testing
# # how many transcripts and genes are there in Salmon output?
# test <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainseq/Salmon_aligned/decoy/R2809.trimmed_transcripts_quant/quant.sf"), data.table=F) # first sample
# test_t <- unlist(lapply(test[,1], function(x) strsplit(x, "[|]")[[1]][1])) # transcript ids with version #
# test_g <- unlist(lapply(test[,1], function(x) strsplit(x, "[|]")[[1]][2])) # gene ids with version #
# test_tt <- unlist(lapply(test_t, function(x) strsplit(x, "[.]")[[1]][1])) # transcript ids without version #
# test_gt <- unlist(lapply(test_g, function(x) strsplit(x, "[.]")[[1]][1])) # gene ids without version #
# # > length(unique(test_t))
# # [1] 243577
# # > length(unique(test_g))
# # [1] 61114
# # > length(unique(test_tt))
# # [1] 243577
# # > length(unique(test_gt))
# # [1] 61114
# # > sum(rownames(txi.salmon[[1]]) %in% test_gt)
# # [1] 60199
# 
# # how many protein coding genes are in Salmon output?
# test_p <- unlist(lapply(test[,1], function(x) strsplit(x, "[|]")[[1]][8]))
# test_pb <- test_p == "protein_coding"
# # > sum(test_pb)
# # [1] 86927
# test_gtp <- test_gt[test_pb]
# # > sum(rownames(txi.salmon[[1]]) %in% test_gtp)
# # [1] 19621
# 
# # > names(txi.salmon)
# # [1] "abundance"           "counts"              "length"             
# # [4] "countsFromAbundance"
# 
# # how are abundance/counts calculated?
# 
# test_gene <- unlist(lapply(test_gt, function(x) grepl("ENSG00000000003", x)))
# test_gene_mat <- test[test_gene,]
# 
# # > test_gene_mat[,2:5]
# #        Length EffectiveLength      TPM NumReads
# # 239603   3768        3618.478 1.719853  253.405
# # 239604   3796        3567.000 0.000000    0.000
# # 239605    900         671.000 0.000000    0.000
# # 239606   1025         717.289 0.705142   20.595
# # 239607    820         591.000 0.000000    0.000
# 
# # > txi.salmon[[1]][1,1] # abundance = sum of TPM of all transcripts
# # ENSG00000000003 
# # 2.424995 
# # > txi.salmon[[2]][1,1] # counts from abundance = sum of NumReads of all transcripts
# # ENSG00000000003 
# # 274 
# # > txi.salmon[[3]][1,1] # length
# # ENSG00000000003 
# # 2774.868 


### SampleNetwork
library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_Brainseq_raw_counts.csv"), data.table=F,)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_Brainseq.csv"), data.table=F)
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
  projectname1="brainseq_samp_filt",
  cexlabels=0.7,
  normalize1=F,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

# apply CTF, asinh normalization to columns
library(edgeR)

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SampleNetworks/1_11-53-31/brainseq_samp_filt_1_190_outliers_removed.csv"), data.table=F)
tmmf <- calcNormFactors(expr[,2:ncol(expr)], method = "TMM")

expr2 <- expr
for(i in 2:ncol(expr)){
  expr2[,i] <- asinh(expr2[,i]/tmmf[i-1])
}
colnames(expr2)[1] <- "gene_id"

fwrite(expr2, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/expr_Brainseq_ctf_asinh.csv"))



### Genomewide cor dist is very right-shifted
### Examine other covariates to investigate

### SampleNetwork
library(tidyverse)
library(data.table)

source(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "code/SampleNetworks/SampleNetwork_1.08.r"))

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_Brainseq_raw_counts.csv"), data.table=F,)
colnames(expr)[1] <- "Gene"
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_Brainseq.csv"), data.table=F)
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
  btrait1=c(1:21), # Age, Sex, Race
  trait1=NULL,
  asfactors1=c(1,2,3,6,7,8,14,15,16,17,18,19,20),
  projectname1="brainseq_samp_filt_testing",
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
  ggtitle("Variance explained by PCs 1-10 in Brainseq data",
          sub = "19356 genes, 204 samples") +
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
          sub = "Brainseq data") +
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
          sub="Brainseq data") +
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
ggsave(p_all, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/4-analyses/Brainseq_PC_covariate_analysis.png"),
       width=15, height=8)

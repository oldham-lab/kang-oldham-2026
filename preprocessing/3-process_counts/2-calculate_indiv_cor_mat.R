# Producing and saving individual simMats for each bulk study

library(WGCNA)
library(data.table)
library(tidyverse)

path <- c(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SampleNetworks/1_11-53-31/brainseq_samp_filt_1_190_outliers_removed.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/GTEx_samp_filt_SampleNetworks/1_09-10-16/GTEx_samp_filt_1_190_outliers_removed.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/NABEC_samp_filt_SampleNetworks/1_04-53-24/NABEC_samp_filt_1_69_outliers_removed.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_SampleNetworks/1_08-13-25/ROSMAP_samp_filt_1_347_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/BrainGVEX_samp_filt_SampleNetworks/1_08-25-36/BrainGVEX_samp_filt_1_275_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_samp_filt_SampleNetworks/1_11-19-13/CMC_samp_filt_1_285_ComBat.csv"),
           file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_HBCC_samp_filt_SampleNetworks/1_10-40-08/CMC_HBCC_samp_filt_1_162_ComBat.csv"))

#path <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/"), full.names=T)
# Remove Brainspan
#path <- path[-grep("Brainspan", path)]
#path <- path[c(1,2,4,5,6,9,10)]
# 
dat_names <- c("Brainseq (n=190)",
               "GTEx (n=190)",
               "NABEC (n=69)",
               "ROSMAP (n=347)",
               "BrainGVEX (n=275)",
               "CMC (n=285)",
               "CMC HBCC (n=162)"
               )
# path <- path[c(2,3,4,6,7)]
# dat_names <- c("GTEx (n=190)",
#                "NABEC (n=63)",
#                "ROSMAP (n=597)",
#                "CMC (n=274)",
#                "CMC HBCC (n=145)"
#                )

# remove Brainseq/BrainGVEX
#path <- path[-c(1,2)]
#dat_names <- dat_names[-c(1,2)]

# load expr matrices
mats <- list()
for(i in 1:length(path)){
  mats[[i]] <- fread(path[i], data.table=F)
}

# align features
features <- mats[[1]][,1]
for(i in 2:length(mats)){
  features <- intersect(features, mats[[i]][,1])
}
for(i in 1:length(mats)){
  mats[[i]] <- mats[[i]][mats[[i]][,1] %in% features,]
  mats[[i]] <- mats[[i]][match(features, mats[[i]][,1]),]
}

# calculate cor mats
cor_mat <- list()
t1 <- timestamp()
for(i in 1:length(mats)){
  cor_mat[[i]] <- bicor(t(mats[[i]][2:ncol(mats[[i]])]), use = "p")
  diag(cor_mat[[i]]) <- 0
  rownames(cor_mat[[i]]) <- mats[[i]][,1]
  colnames(cor_mat[[i]]) <- mats[[i]][,1]
  timestamp()
  cat(i, " ")
}
t2 <- timestamp()

# Convert to gene symbols
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
cor_mat <- lapply(cor_mat, function(x){
  x[rownames(x) %in% expr[,1], colnames(x) %in% expr[,1]]
  rownames(x) <- expr[,2]
  colnames(x) <- expr[,2]
  return(x)
})

# save cor mats
#names(cor_mat) <- dat_names
#save(cor_mat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_bicor_individual_5_datasets.Rdata"))

# save individual mats
#cor_mat_f <- cor_mat[-c(8,9)]
#save(cor_mat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_bicor_individual_minus_Brainseq_BrainGVEX.Rdata"))
#save(cor_mat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_bicor_individual_non_norm_minus_Brainseq_BrainGVEX.Rdata"))
saveRDS(cor_mat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_bicor_individual_final_non_norm.RDS"))


# save transposed expr matrices
# mats_t <- lapply(mats, function(x) t(x[,2:ncol(x)]))
# for(i in 1:length(mats_t)){
#   colnames(mats_t[[i]]) <- mats[[1]][,1]
# }
# save(mats_t, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_transposed_list.Rdata"))

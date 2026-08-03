# Calculate genomewide correlation distributions for each dataset

library(WGCNA)
library(data.table)
library(tidyverse)

path <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/processed_mats/"), full.names=T)
# Remove Brainspan
#path <- path[-grep("Brainspan", path)]
path <- path[c(1,2,4,5,6,9,10)]
path <- c(path, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_12-23-52/combined_FCX_final_1_1706_outliers_removed.csv"))
path <- c(path, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/BrainGVEX_samp_filt_SampleNetworks/1_08-25-36/BrainGVEX_samp_filt_1_275_ComBat.csv"))
path <- c(path, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SampleNetworks/1_11-53-31/brainseq_samp_filt_1_190_outliers_removed.csv"))
dat_names <- c("BrainGVEX (n=255)",
               "Brainseq (n=182)", 
               "CMC (n=274)",
               "CMC HBCC (n=145)",
               "GTEx (n=190)", 
               "NABEC (n=63)", 
               "ROSMAP (n=597)",
               "combined (n=1708)",
               "BrainGVEX, no norm (n=255)",
               "Brainseq, no norm (n=182)")

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
for(i in 1:length(mats)){
  cor_mat[[i]] <- bicor(t(mats[[i]][2:ncol(mats[[i]])]), use = "p")
  diag(cor_mat[[i]]) <- 0
  rownames(cor_mat[[i]]) <- mats[[i]][,1]
  colnames(cor_mat[[i]]) <- mats[[i]][,1]
}

# save cor mats
#names(cor_mat) <- dat_names
#save(cor_mat, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_bicor_individual_5_datasets.Rdata"))

# save individual mats
#cor_mat_f <- cor_mat[-c(8,9)]
#save(cor_mat_f, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/simMat_bicor_individual_final.Rdata"))


# save transposed expr matrices
# mats_t <- lapply(mats, function(x) t(x[,2:ncol(x)]))
# for(i in 1:length(mats_t)){
#   colnames(mats_t[[i]]) <- mats[[1]][,1]
# }
# save(mats_t, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_transposed_list.Rdata"))

# plot distributions of pairwise correlations
set.seed(121)
plot_df <- list()
for(i in 1:length(cor_mat)){
  plot_df[[i]] <- data.frame("Dataset" = rep(dat_names[i], 10000),
                             sample(cor_mat[[i]][upper.tri(cor_mat[[i]])], 10000))
}

plot_df <- data.frame(do.call(rbind, plot_df))
colnames(plot_df)[2] <- "cor"

p <- ggplot(plot_df, aes(x = cor, color = Dataset)) +
  geom_density() +
  xlab("Bicor") +
  ggtitle("Genomewide pairwise correlations", 
          sub = paste0(length(features), " features")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/4-analyses/1-genomewide_bicor_dist.png"))

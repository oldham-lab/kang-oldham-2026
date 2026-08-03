# Plot cor dists of all 7 datasets used in bulk megaset analysis

library(ggridges)
library(edgeR)
library(data.table)
library(tidyverse)
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))

############
# Functions
############

get_10k_cors <- function(expr, 
                         ind=2,
                         no_cors=10000,
                         seed=3){ 
  set.seed(seed)
  # transpose and remove zero var genes
  expr_t <- t(expr[,ind:ncol(expr)])
  #expr_var <- apply(expr_t,2,var)
  #expr_t <- expr_t[,expr_var>0]
  #expr <- expr[expr_var>0,]
  # calc similarities and sample 10k of them
  simMat <- cor(expr_t)
  raw_cors <- simMat[upper.tri(simMat)]
  raw_cors <- sample(raw_cors,no_cors)
  return(raw_cors)
}

get_10k_corsv2 <- function(expr_t, 
                         ind=2,
                         no_cors=10000,
                         seed=3){ 
  set.seed(seed)
#  expr_t <- t(expr[,ind:ncol(expr)])
  #expr_var <- apply(expr_t,2,var)
  #expr_t <- expr_t[,expr_var>0]
  #expr <- expr[expr_var>0,]
  # calc similarities and sample 10k of them
  simMat <- cor(expr_t)
  raw_cors <- simMat[upper.tri(simMat)]
  raw_cors <- sample(raw_cors,no_cors)
  return(raw_cors)
}

get_all_cors <- function(expr,
                         ind=2){
  # transpose and remove zero var genes
  expr_t <- t(expr[,ind:ncol(expr)])
  expr_var <- apply(expr_t,2,var)
  expr_t <- expr_t[,expr_var>0]
  expr <- expr[expr_var>0,]
  # calc similarities and sample 10k of them
  simMat <- cor(expr_t)
  raw_cors <- simMat[upper.tri(simMat)]
  return(raw_cors)
}

scale_pos <- function(expr,
                      si=NULL #sample index
                      ){
  if(is.null(si)){si <- 2:ncol(expr)}
  expr[,si] <- expr[,si]+abs(min(expr[,si],na.rm=T))
  return(expr)
}

ctf <- function(expr,
                si = 2:ncol(expr), # sample index
                gi = 1 # gene index
                ){
  y <- calcNormFactors(expr[,si], method = "TMM")
  expr[,si] <- expr[,si]/y
  return(expr)
} 

lognorm <- function(expr,
                    si=2:ncol(expr)){
  expr[,si] <- log10(expr[,si]+1)
  return(expr)
}

asinhnorm <- function(expr,
                      si=2:ncol(expr)){
  expr[,si] <- asinh(expr[,si])
  return(expr)
}

# format cors into df
make_plot_df <- function(cor_list, #list
                         names=NULL, #vec
                         dataset
){
  if(is.null(names)){
    names <- names(cor_list)
  }
  dfs <- list()
  for(i in 1:length(cor_list)){
    dfs[[i]] <- data.frame("type"=names[i],
                           "vals"=cor_list[[i]],
                           "dataset"=dataset)
  }
  dfs <- do.call(rbind, dfs)
  dfs$type <- factor(dfs$type, levels=rev(names))
  return(dfs)
}

# calculate all different combinations of normalization methods
calc_exprs <- function(expr){
  expr_list <- list()
  expr_list[[1]] <- expr 
  expr_list[[2]] <- expr %>% ctf
  expr_list[[3]] <- expr %>% asinhnorm
  expr_list[[4]] <- expr %>% ctf %>% asinhnorm
  expr_list[[5]] <- expr %>% scale_pos %>% ctf
  expr_list[[6]] <- expr %>% scale_pos %>% asinhnorm
  expr_list[[7]] <- expr %>% scale_pos %>% ctf %>% asinhnorm
  expr_list[[8]] <- expr %>% scale_pos %>% lognorm
  expr_list[[9]] <- expr %>% scale_pos %>% ctf %>% lognorm
  return(expr_list)
}

# wrapper function
# start with expr, end with dataframe of cors
get_cors_full <- function(expr_list,
                          save_name,
                          n,
                          subset=T
                          ){
  full_list <- unlist(lapply(expr_list, calc_exprs), recursive=F)
  comb_vec <- c("CTF","asinh","CTF + asinh","Scale + CTF","Scale + asinh","Scale + CTF + asinh","Scale + log","Scale + CTF + log")
  names(full_list) <- c("Raw counts",paste0("Raw counts + ", comb_vec), 
                        "Filt + ComBat", paste0("Filt + ComBat + ", comb_vec))
  
  if(subset){
    full_cors <- lapply(full_list, get_10k_cors)
    saveRDS(full_cors,file=paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "0-preprocess_cor_dists/norm_compare_cor_sample_10k_"),save_name,".RDS"))
  } else {
    full_cors <- lapply(full_list, get_all_cors)
    saveRDS(full_cors,file=paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "0-preprocess_cor_dists/norm_compare_cor_all_"),save_name,".RDS"))
  }
  
  plot_df <- make_plot_df(full_cors,
                          dataset=paste0(save_name," (n=",n,")"))
  return(plot_df)
}

##########
# Code
##########

expr_list <- list("BrainGVEX" = list("Raw counts"= fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_BrainGVEX_raw_counts.csv"),data.table=F),
                                     "Filt + ComBat" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/BrainGVEX_samp_filt_SampleNetworks/1_08-25-36/BrainGVEX_samp_filt_1_275_ComBat.csv"),data.table=F)),
                  "Brainseq" = list("Raw counts" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_Brainseq_raw_counts.csv"), data.table=F),
                                    "Filt + ComBat" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SampleNetworks/1_11-53-31/brainseq_samp_filt_1_190_outliers_removed.csv"), data.table=F)),
                  "CMC" = list("Raw counts" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_CMC_only_raw_counts.csv"), data.table=F),
                               "Filt + ComBat" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_samp_filt_SampleNetworks/1_11-19-13/CMC_samp_filt_1_285_ComBat.csv"), data.table=F)),
                  "CMC_HBCC" = list("Raw counts" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_CMC_HBCC_raw_counts.csv"), data.table=F),
                                    "Filt + ComBat" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/CMC_HBCC_samp_filt_SampleNetworks/1_10-40-08/CMC_HBCC_samp_filt_1_162_ComBat.csv"), data.table=F)),
                  "GTEx" = list("Raw counts" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_GTEx_raw_counts.csv"), data.table=F),
                                "Filt + ComBat" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/GTEx_samp_filt_combat_library_date_batch_SampleNetworks/1_03-44-38/GTEx_samp_filt_combat_library_date_batch_1_183_ComBat.csv"), data.table=F)),
                  "NABEC" = list("Raw counts" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_NABEC_raw_counts.csv"), data.table=F),
                                 "Filt + ComBat" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/NABEC_samp_filt_SampleNetworks/1_04-53-24/NABEC_samp_filt_1_69_outliers_removed.csv"), data.table=F)),
                  "ROSMAP" = list("Raw counts" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_ROSMAP_raw_counts.csv"), data.table=F),
                                  "Filt + ComBat" = fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_SampleNetworks/1_08-59-36/ROSMAP_samp_filt_1_617_ComBat.csv"), data.table=F))
)

samp_sizes <- unlist(lapply(expr_list, function(x) ncol(x[[1]])-1))
raw_common_genes <- Reduce(intersect, lapply(expr_list, function(x) x[[1]][,1]))
filt_common_genes <- Reduce(intersect, lapply(expr_list, function(x) x[[2]][,1]))
for(i in 1:length(expr_list)){
  expr_list[[i]][[1]] <- expr_list[[i]][[1]][expr_list[[i]][[1]][,1] %in% raw_common_genes,]
  expr_list[[i]][[2]] <- expr_list[[i]][[2]][expr_list[[i]][[2]][,1] %in% filt_common_genes,]
}

plot_list <- mapply(get_cors_full, expr_list, names(expr_list), samp_sizes, SIMPLIFY=FALSE)
saveRDS(plot_list, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "0-preprocess_cor_dists/7_dataset_10k_cor_sample.RDS"))

plot_list <- readRDS(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "0-preprocess_cor_dists/7_dataset_10k_cor_sample.RDS"))
plot_df <- do.call(rbind, plot_list)
plot_df$dataset <- gsub("[(]", " (", plot_df$dataset)

# Calculate cors of rows (cors between each dataset for each preprocessing pipeline)
plot_df2 <- list()
for(i in 1:length(unique(plot_df$type))){
  bla <- plot_df %>% dplyr::filter(type==unique(plot_df$type)[i]) 
  bla <- tapply(bla$vals, bla$dataset, list)
  bla <- as.data.frame(do.call(cbind,bla))
  bla <- cor(bla, use='pairwise.complete.obs')
  bla <- bla[upper.tri(bla)]
  bla <- data.frame("type"=unique(plot_df$type)[i],
                    "vals"=bla,
                    "dataset"="Row cors")
  plot_df2[[i]] <- bla
}
plot_df2 <- as.data.frame(do.call(rbind,plot_df2))
plot_df2 %>% group_by(type) %>% summarise(mean=mean(vals)) %>% arrange(desc(mean))
# type                                 mean
# <fct>                               <dbl>
# 1 Filt + ComBat + Scale + log         0.669
# 2 Filt + ComBat + Scale + asinh       0.669
# 3 Filt + ComBat + asinh               0.667
# 4 Filt + ComBat                       0.656
# 5 Filt + ComBat + CTF + asinh         0.651
# 6 Raw counts + Scale + log            0.634
# 7 Raw counts                          0.631
# 8 Filt + ComBat + CTF                 0.630
# 9 Raw counts + Scale + asinh          0.630
# 10 Raw counts + asinh                  0.630
# 11 Filt + ComBat + Scale + CTF + log   0.627
# 12 Filt + ComBat + Scale + CTF + asinh 0.626
# 13 Raw counts + Scale + CTF + log      0.618
# 14 Raw counts + Scale + CTF + asinh    0.614
# 15 Raw counts + CTF + asinh            0.614
# 16 Filt + ComBat + Scale + CTF         0.612
# 17 Raw counts + Scale + CTF            0.603
# 18 Raw counts + CTF                    0.603

plotdf <- rbind(plot_df,plot_df2)
plotdf$dataset <- factor(plotdf$dataset, levels=c(unique(plot_df$dataset),"Row cors"))
p<-ggplot(plotdf, aes(x=vals,y=type)) +
  geom_density_ridges(alpha=0.3) +
  geom_vline(xintercept=0,color="red",alpha=0.6) +
  th +
  xlab("") +
  ylab("") +
  ggtitle("Pairwise cor distributions for bulk megaset datasets",
          sub="19573 genes, 10000 cors sampled") +
  facet_wrap(~dataset, ncol = 8)
ggsave(p,width=14,height=7, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/4-analyses/Full_comparison_of_preprocessing_normalization_strategies.png"))


p2<-ggplot(plot_df2, aes(x=vals)) +
  geom_density(alpha=0.3) +
  geom_vline(xintercept=0,color="red",alpha=0.6) +
  th +
  xlab("") +
  ylab("") +
  ggtitle("Pairwise cor distributions for bulk megaset datasets",
          sub="19573 genes, 10000 cors sampled") +
  facet_wrap(~type)

bla <- plot_df2 %>% group_by(type,dataset) %>% summarise("mean"=mean(vals)) %>% arrange(desc(mean)) 
# 1 Filt + ComBat + Scale + CTF         Cors of row 0.00483 
# 2 Filt + ComBat + Scale + CTF + log   Cors of row 0.00440 
# 3 Filt + ComBat + Scale + CTF + asinh Cors of row 0.00437 
# 4 Filt + ComBat                       Cors of row 0.00395 
# 5 Filt + ComBat + asinh               Cors of row 0.00392 
# 6 Filt + ComBat + Scale + asinh       Cors of row 0.00379 
# 7 Filt + ComBat + Scale + log         Cors of row 0.00372 
# 8 Filt + ComBat + CTF + asinh         Cors of row 0.00369 
# 9 Filt + ComBat + CTF                 Cors of row 0.00366 
# 10 Raw counts + Scale + asinh          Cors of row 0.00175 
# 11 Raw counts + asinh                  Cors of row 0.00175 
# 12 Raw counts + Scale + log            Cors of row 0.00161 
# 13 Raw counts + Scale + CTF + asinh    Cors of row 0.00113 
# 14 Raw counts + CTF + asinh            Cors of row 0.00113 
# 15 Raw counts + Scale + CTF + log      Cors of row 0.000994
# 16 Raw counts                          Cors of row 0.000878
# 17 Raw counts + Scale + CTF            Cors of row 0.000690
# 18 Raw counts + CTF                    Cors of row 0.000690


# Calculate cors between pairwise cors
test <- plot_list[[1]]
test$index <- rep(1:10000, 18)
test <- pivot_wider(test, names_from="type", values_from="vals")
testcors <- cor(test[,3:ncol(test)],use='pairwise.complete.obs')
summary(testcors[upper.tri(testcors)])

##########
# Plot megaset cors
##########

expr_list <- list("No normalization" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined.csv"), data.table=F),
                  "Filt + ComBat" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_11-42-16/combined_FCX_final_1_1788_ComBat.csv"),data.table=F),
                  "Filt + ComBat + CTF + asinh" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_11-42-16/combined_FCX_final_1_1788_ComBat_ctf_asinh.csv"),data.table=F),
                  "Filt + ComBat + Scale + CTF" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_11-42-16/combined_FCX_final_1_1788_ComBat_scale_ctf.csv"),data.table=F),
                  "Filt + ComBat + Scale + CTF + asinh" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_11-42-16/combined_FCX_final_1_1788_ComBat_scale_ctf_asinh.csv"),data.table=F),
                  "5 dataset No normalization" = fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC.csv"), data.table=F),
                  "5 dataset Filt + ComBat" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24//combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat.csv"),data.table=F),
                  "5 dataset Filt + ComBat + CTF + asinh" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_ctf_asinh.csv"),data.table=F),
                  "5 dataset Filt + ComBat + Scale + CTF + asinh" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_scale_ctf.csv"),data.table=F),
                  "5 dataset Filt + ComBat + Scale + CTF + asinh" = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_scale_ctf_asinh.csv"),data.table=F)
                  )



expr_cors <- lapply(expr_list, get_10k_cors, ind=3)
names(expr_cors) <- c("(7 dataset) No normalization","(7 dataset) Filt + ComBat","(7 dataset) Filt + ComBat + CTF + asinh","(7 dataset) Filt + ComBat + Scale + CTF", "(7 dataset) Filt + ComBat + Scale + CTF + asinh",
                      "(5 dataset) No normalization","(5 dataset) Filt + ComBat","(5 dataset) Filt + ComBat + CTF + asinh","(5 dataset) Filt + ComBat + Scale + CTF", "(5 dataset) Filt + ComBat + Scale + CTF + asinh")
temp <- fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_asinh.csv"), data.table=F)
expr_cors[[11]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[11] <- "(5 dataset scale + CTF) Filt + ComBat + asinh"
temp <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_scaled.csv"), data.table=F)
expr_cors[[12]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[12] <- "(5 dataset) scale"
temp <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat.csv"), data.table=F)
expr_cors[[13]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[13] <- "(5 dataset) scale + Filt + ComBat"
temp <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat_asinh.csv"), data.table=F)
expr_cors[[14]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[14] <- "(5 dataset) scale + Filt + ComBat + asinh"
temp <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat_ctf_asinh.csv"), data.table=F)
expr_cors[[15]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[15] <- "(5 dataset) scale + Filt + ComBat + ctf + asinh"
temp <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_scaled_SampleNetworks/1_07-00-42/combined_FCX_minus_Brainseq_NABEC_scaled_1_1529_ComBat_scale_ctf_asinh.csv"), data.table=F)
expr_cors[[16]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[16] <- "(5 dataset) scale + Filt + ComBat + scale + ctf + asinh"
temp <- fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat.csv"), data.table=F)
expr_cors[[17]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[17] <- "(5 dataset scale + CTF) Filt + ComBat"
temp <- fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale_ctf.csv"), data.table=F)
expr_cors[[18]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[18] <- "(5 dataset scale + CTF) Filt + ComBat + scale + ctf"
temp <- fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale_ctf_asinh.csv"), data.table=F)
expr_cors[[19]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[19] <- "(5 dataset scale + CTF) Filt + ComBat + scale + ctf + asinh"
temp <- fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale_asinh.csv"), data.table=F)
expr_cors[[20]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[20] <- "(5 dataset scale + CTF) Filt + ComBat + scale + asinh"
temp <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf.csv"), data.table=F)
expr_cors[[21]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[21] <- "(5 dataset scale + CTF) None"
temp <- fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_SampleNetworks/1_11-06-47/combined_FCX_minus_Brainseq_NABEC_firstScaledCtf_1_1522_ComBat_scale.csv"), data.table=F)
expr_cors[[22]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[22] <- "(5 dataset scale + CTF) Filt + ComBat + scaled"
temp <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_scaled.csv"), data.table=F)
expr_cors[[23]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[23] <- "(5 dataset scale + CTF) scaled"
temp <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_scaled_ctf.csv"), data.table=F)
expr_cors[[24]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[24] <- "(5 dataset scale + CTF) scaled + CTF"
temp <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_scaled_ctf_asinh.csv"), data.table=F)
expr_cors[[25]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[25] <- "(5 dataset scale + CTF) scaled + CTF + asinh"
temp <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC_firstScaledCtf_ctf.csv"), data.table=F)
expr_cors[[26]] <- get_10k_cors(temp, ind=3)
names(expr_cors)[26] <- "(5 dataset scale + CTF) CTF"

temp = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_ctf_asinh.csv"),data.table=F)
expr_cors[[8]] <- get_10k_cors(temp,ind=3)
names(expr_cors)[[8]] <- "(5 dataset) Filt + ComBat + CTF + asinh"

temp = fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_minus_Brainseq_NABEC_SampleNetworks/1_01-50-24/combined_FCX_minus_Brainseq_NABEC_1_1529_ComBat_ctf_scale.csv"),data.table=F)
expr_cors[[27]] <- get_10k_cors(temp,ind=3)
names(expr_cors)[[27]] <- "(5 dataset) Filt + ComBat + CTF + scale"


saveRDS(expr_cors, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "0-preprocess_cor_dists/megaset_10k_cors_sample.RDS"))



expr_cors <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "0-preprocess_cor_dists/megaset_10k_cors_sample.RDS"))
plot_df <- make_plot_df(expr_cors, dataset="all")

p<-ggplot(plot_df, aes(x=vals,y=type)) +
  geom_density_ridges(alpha=0.3) +
  geom_vline(xintercept=0,color="red",alpha=0.6) +
  th +
  xlab("") +
  ylab("") +
  ggtitle("Pairwise cor distributions for bulk megaset dataset",
          sub="10000 cors sampled") 
ggsave(p, width=8, height=20, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/4-analyses/Full_comparison_of_preprocessing_normalization_strategies_combined_megaset.png"))


# Some testing

# take 100 genes
test <- expr_list[[6]][sample(3:nrow(expr_list[[6]]),1000),3:ncol(expr_list[[6]])]
cortest <- cor(t(test))
plot(density(cortest[upper.tri(cortest)]))

test2 <- test+abs(min(test,na.rm=T))
cortest2 <- cor(t(test2))
plot(density(cortest2[upper.tri(cortest2)]))

test3 <- test+max(test)
cortest3 <- cor(t(test3))
plot(density(cortest3[upper.tri(cortest3)]))


expr_corsDF <- do.call(cbind, expr_cors)
expr_corsDF <- expr_corsDF[,grep("(5 dataset)", colnames(expr_corsDF))]
expr_cdc <- cor(expr_corsDF)
plot(density(expr_cdc[upper.tri(expr_cdc)]))
expr_cdc[,1]

# (5 dataset) No normalization 
# 1.00000000 
# (5 dataset) Filt + ComBat 
# 0.56826827 
# (5 dataset) Filt + ComBat + CTF + asinh 
# 0.55744956 
# (5 dataset) Filt + ComBat + Scale + CTF 
# 0.04512828 
# (5 dataset) Filt + ComBat + Scale + CTF + asinh 
# 0.04509318 
# (5 dataset scale + CTF) Filt + ComBat + asinh 
# 0.02329546 
# (5 dataset) scale 
# 1.00000000 
# (5 dataset) scale + Filt + ComBat 
# 0.56826827 
# (5 dataset) scale + Filt + ComBat + asinh 
# 0.56745089 
# (5 dataset) scale + Filt + ComBat + ctf + asinh 
# 0.05641044 
# (5 dataset) scale + Filt + ComBat + scale + ctf + asinh 
# 0.04509340 
# (5 dataset scale + CTF) Filt + ComBat 
# 0.02318566 
# (5 dataset scale + CTF) Filt + ComBat + scale + ctf 
# 0.01736384 
# (5 dataset scale + CTF) Filt + ComBat + scale + ctf + asinh 
# 0.01740283 
# (5 dataset scale + CTF) Filt + ComBat + scale + asinh 
# 0.02322829 
# (5 dataset scale + CTF) None 
# 0.01349468 
# (5 dataset scale + CTF) Filt + ComBat + scaled 
# 0.02318566 
# (5 dataset scale + CTF) scaled 
# 0.01349468 
# (5 dataset scale + CTF) scaled + CTF 
# 0.01931604 
# (5 dataset scale + CTF) scaled + CTF + asinh 
# 0.01923321 
# (5 dataset scale + CTF) CTF 
# 0.01931604 
# (5 dataset) Filt + ComBat + CTF + scale 
# 0.56963302 

# What about cors with individual datasets
# load braingvex
exprBG <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_BrainGVEX_raw_counts.csv"),data.table=F)
exprcomb <-  fread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/expr_combined_minus_Brainseq_NABEC.csv"), data.table=F)
exprBG <- exprBG[exprBG[,1] %in% exprcomb[,1],]
exprBGc <- get_10k_cors(exprBG,ind=3)

cor(expr_corsDF, exprBGc)

# seems like the combination of scale + ctf greatly distorts the cor distribution
# (even though it produces the most normal looking distributions)

                                     


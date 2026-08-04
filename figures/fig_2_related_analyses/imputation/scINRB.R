# scINRB: single-cell gene expression imputation with network regularization and bulk RNA-seq data Open Access 
# https://academic.oup.com/bib/article/25/3/bbae148/7642689
# doi.org/10.1093/bib/bbae148
# https://github.com/JGuan-lab/scINRB

# Run on Lein DFC with bulk megaset as bulk input

library(data.table)
library(tidyverse)
library(qs)
library(MASS)
source(file.path(Sys.getenv("SCINRB_DIR", "/home/gugene/code/git_other/scINRB"), "code/scINRB.R"))
source(file.path(Sys.getenv("SCINRB_DIR", "/home/gugene/code/git_other/scINRB"), "code/functions.R"))

# Load bulk data
bulkmat <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"))
# Load cellbender-corrected lein DFC
lein <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Raw_data_bag_4_Human_Cross_Areal_raw_10x/data/DLPFC/aligned_data/lein_2023_dfc_cellbender_expr_matchedLabels.csv")) %>%
  column_to_rownames(var="V1")
# Calculate bulk means and match bulk genes with sn genes
bulkmeans <- rowMeans(bulkmat[,-c(1,2)]) %>% data.frame
rownames(bulkmeans) <- bulkmat$Gene
common_genes <- intersect(bulkmat$Gene, rownames(lein))
lein <- lein[rownames(lein) %in% common_genes,]
lein <- lein[match(common_genes, rownames(lein)),]
bulkmeans <- bulkmeans[rownames(bulkmeans) %in% common_genes,,drop=F]
bulkmeans <- bulkmeans[match(common_genes, rownames(bulkmeans)),,drop=F]

# Prepare data
result <- preprocess(as.matrix(lein),as.matrix(bulkmeans))
data_sc <- result[[1]]
data_bulk <- result[[2]]

# Run with default parameters
parameter <- c(0.001,0.001,1) 
r <- 200
t1 <- timestamp()
result <- scINRB(data_sc, data_bulk, parameter, r)
t2 <- timestamp()
# timed out after 3-4 days
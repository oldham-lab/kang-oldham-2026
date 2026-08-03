# https://github.com/satijalab/seurat-wrappers/blob/master/docs/alra.md

library(Seurat)
library(SeuratWrappers)
library(qs)
library(tidyverse)
library(data.table)

dat <- readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/DFC.rds")) %>% 
  NormalizeData()

dat <- RunALRA(dat)
qsave(dat, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/ALRA/lein_DFC_seurat_obj_with_ALRA_run_normalized.qs"))
fwrite(dat@assays$alra@data)

# Output for ALRA run (normalized):
# Rank k = 64
# Identifying non-zero values
# Computing Randomized SVD
# Find the 0.001000 quantile of each gene
# Thresholding by the most negative value of each gene
# Scaling all except for 7943 columns
# 1.02% of the values became negative in the scaling process and were set to zero
# The matrix went from 18.12% nonzero to 49.48% nonzero
# Warning: Layer counts isn't present in the assay object; returning NULL
# Setting default assay as alra


# Run Alra for cellbender data
datmat <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/lein_2023_raw/Raw_data_bag_4_Human_Cross_Areal_raw_10x/data/DLPFC/aligned_data/lein_2023_dfc_cellbender_expr_matchedLabels.csv"))
rownames(datmat) <- datmat[,1]
datmat <- datmat[,-1]
dat <- CreateSeuratObject(counts = datmat, project = "leindfc", min.cells = 3, min.features = 200) %>% 
  NormalizeData()
dat <- RunALRA(dat)
qsave(dat, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/ALRA/lein_DFC_seurat_obj_with_ALRA_run_cellbender.qs"))
# Rank k = 51
# Identifying non-zero values
# Computing Randomized SVD
# Find the 0.001000 quantile of each gene
# Thresholding by the most negative value of each gene
# Scaling all except for 0 columns
# 0.73% of the values became negative in the scaling process and were set to zero
# The matrix went from 23.85% nonzero to 61.77% nonzero
# Warning: Layer counts isn't present in the assay object; returning NULL
# Setting default assay as alra
# overview of imputation methods: 
# https://doi.org/10.1186/s13059-020-02132-x
# https://doi.org/10.1186/s12859-023-05417-7

library(SAVER)
library(qs)

cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
#min.cells = 3, min.features = 200
mincells <- apply(cell_exprall,1,function(x) sum(x>0)) # 8122 features expressed in fewer than 3 cells
cell_exprall <- cell_exprall[mincells>=3,]
#minfeatures <- apply(cell_exprall,2,function(x) sum(x>0)) # no cells expressed in fewer than 200 features

#cortex.saver <- saver(cell_exprall, ncores = 8)
#qsave(cortex.saver, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/SAVER/SAVER_output.qs"))
# this crashes

#saver1 <- saver(cell_exprall, pred.genes = 1:500, pred.genes.only = TRUE, 
#                do.fast = FALSE)
# still takes way too long

t1 <- timestamp()
saver1 <- saver(cell_exprall, pred.genes = 1:100, pred.genes.only = TRUE, 
                do.fast = FALSE, ncores=6)
t2 <- timestamp()
# 8 hrs for 100 genes, prohibitively slow
qsave(saver1, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/SAVER/output_genes_1_to_100.qs"))

t1 <- timestamp()
saver1 <- saver(cell_exprall, pred.genes = 101:200, pred.genes.only = TRUE, 
                do.fast = TRUE, ncores=6)
t2 <- timestamp()
qsave(saver1, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/SAVER/output_genes_101_to_200.qs"))
# 5 hours, still too slow

inputgenevec1 = c("AIF1","ALDH1L1","MOG","SLC17A7","GAD1")
genes.ind <- which(rownames(cell_exprall) %in% inputgenevec1)

# Generate predictions for those genes and return entire dataset
cortex.saver.genes <- saver(cell_exprall, pred.genes = genes.ind, 
                            estimates.only = TRUE, ncores=5)
qsave(cortex.saver.genes, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/SAVER/output_aif1_aldh1l1_mog_slc17a7_gad1.qs"))

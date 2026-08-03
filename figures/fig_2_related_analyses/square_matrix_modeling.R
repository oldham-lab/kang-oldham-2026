# Recapitulate Rebecca's square matrix modeling in Lein DFC
# see: /mnt/bdata/rebecca_home/SCSN_meta_analysis/frontal_cortex/analyses/04_square_matrix/square_matrix_fxn.R
# also see R01 #3

library(tidyverse)
library(data.table)
library(qs)
options(bitmapType = 'cairo')

# Load Lein DFC
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
n_resamples <- 10
sample_size <- 2000

## Sample n genes and nuclei from expression data:
mat_list <- lapply(1:n_resamples, FUN=function(n){
    set.seed(n); row_index <- sample(1:nrow(cell_exprall), size=sample_size)
    set.seed(n); col_index <- sample(1:ncol(cell_exprall), size=sample_size)
    expr1 <- cell_exprall[row_index, col_index]
    ## Avoid zero variance genes:
    zero_var <- rowSums(expr1) == 0
    iter <- 1e6
    while(any(zero_var)){
        row_index <- row_index[!zero_var]
        set.seed(iter)
        row_index <- c(row_index, sample(setdiff(1:nrow(cell_exprall), row_index), size=sum(zero_var)))
        expr1 <- cell_exprall[row_index, col_index]
        zero_var <- rowSums(expr1) == 0
        iter <- iter + 1
    }
    return(expr1)
})
qsave(mat_list, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/square_mats/matlist.qs"))


# Load iterations
mat_list <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/square_mats/matlist.qs"))
# For each iteration, calculate pc1 in each dimension and model
setwd(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/square_mats/"))
#source("/mnt/bdata/rebecca_home/SCSN_meta_analysis/frontal_cortex/analyses/04_square_matrix/analyses/01_variance_explained/variance_explained_fxn.R")
library(dplyr) 
library(gmodels) ## fast.prcomp()
#library(flexiblas)
library(future.apply)

#flexiblas_switch(invisible(flexiblas_load_backend("NETLIB")))

options(future.globals.maxSize=Inf)
plan(multicore, workers=5)

PCs_variance_explained <- function(mat_list,
                                   data_type,
                                   expr_type,
                                   sample_size,
                                   n_resamples){
  
    var_expl <- future_lapply(mat_list, FUN=function(expr){
      col_pcs <- gmodels::fast.prcomp(expr, scale.=T, retx=F)
      col_var_expl <- col_pcs$sdev^2/(sum(col_pcs$sdev^2))
      row_pcs <- gmodels::fast.prcomp(t(expr), scale.=T, retx=F)
      row_var_expl <- row_pcs$sdev^2/(sum(row_pcs$sdev^2))
      return(data.frame(iteration=i,
                        PC=1:sample_size, 
                        Row_VE=row_var_expl, 
                        Col_VE=col_var_expl))
    })
  
  file_path <- paste0("square_matrix_PCA_variance_explained_", data_type, 
                      "_", expr_type, "_", sample_size, "_samples_", n_resamples, "_resamples_", 
                      length(mat_list), "_datasets.qs")
  qsave(var_expl, file=file_path)
  
}

data_type <- "author_data"
expr_type <- "umi_counts"

PCs_variance_explained(mat_list,
                       data_type,
                       expr_type,
                       sample_size,
                       n_resamples)

# Plot row vs col variance
varobj <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/square_mats/square_matrix_PCA_variance_explained_author_data_umi_counts_2000_samples_10_resamples_10_datasets.qs"))

plotobj <- lapply(varobj, function(x){
    out <- data.frame(type=c("Gene\n(n=10)", "Cell\n(n=10)"),
                      ve=c(x[1,3], x[1,4]))
}) %>% do.call(rbind, .)

p <- ggplot(plotobj, aes(x=type, y=ve)) +
  theme_bw() +
  geom_jitter() +
  geom_boxplot(width=0.2,fill="white") +
  labs(x="", y="% Variance explained") +
  theme(text=element_text(size=30))
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/square_mat_pc1.png"), width=4,height=8)

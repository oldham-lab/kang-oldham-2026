# Evaluate whether cells for a given celltype are more similar to each other than expected by chance in Lein DFC

library(qs)
library(data.table)
library(tidyverse)
library(flexiblas)
library(reticulate)
library(qvalue)
library(Seurat)
numpy <- import("numpy")
options(bitmapType = 'cairo')


cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
#genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)

dfc_seurat <- CreateSeuratObject(counts = cell_exprall, project = "pbmc3k", min.cells = 3, min.features = 200) %>%
    NormalizeData %>%
    FindVariableFeatures(nfeatures=3000, selection.method="vst")
var_features <- VariableFeatures(dfc_seurat)
match_features <- gsub("_", "-", rownames(cell_exprall))
cell_exprall_hivar <- cell_exprall[match_features %in% var_features,]

# subset <- cell_exprall[, colnames(cell_exprall) %in% cell_annoall$Cell_ID[cell_annoall[,1]=="Sst Chodl"]] %>% t
# subset_var <- apply(subset,2,var)
# subset <- subset[,subset_var>0]
# system.time(bla1 <- cor(subset))
# #   user  system elapsed 
# # 52.885   0.441  53.333
# system.time(bla2 <- numpy$corrcoef(t(subset)))
# #   user  system elapsed 
# # 16.884   1.476  11.165

cts <- c(table(cell_annoall[,1])) %>% sort

############ Top 3000 hvgs
corpctg <- list()
plotlist <- list()
for(i in seq_along(cts)){
    cells <- cell_annoall$Cell_ID[cell_annoall[,1] == "VLMC"]
    subset <- cell_exprall_hivar[, colnames(cell_exprall_hivar) %in% cells]
    subset_var <- apply(subset,1,var)
    subset <- subset[subset_var>0,]

    real_cors <- numpy$corrcoef(subset, rowvar=F)
    diag(real_cors) <- 0
    
    set.seed(i)
    subset_perm <- apply(subset,2,function(x){
        sample(x, length(x), replace=F)
    })
    rand_cors <- numpy$corrcoef(subset_perm, rowvar=F)
    diag(rand_cors) <- 0

    # Calculate cutoffs
    r2 <- rand_cors[upper.tri(rand_cors)]
    norm_cut <- quantile(r2, 0.95)
    r1 <- real_cors[upper.tri(real_cors)]
    corpctg[[i]] <- sum(r1>norm_cut)/length(r1)

    plotlist[[i]] <- data.frame("cors"=c(sample(r1, min(length(r1),10000), replace=F),sample(r2,min(length(r2),10000), replace=F)),
                                "type"=c(rep("real",min(length(r1),10000)), rep("permuted",min(length(r2),10000))),
                                "region"=names(cts)[i],
                                "cut"=norm_cut)

    cat(i, " ")
}
unlist(corpctg)
#  [1] 0.4293040 0.4357143 0.4369963 0.4201465 0.4236264 0.4313187 0.4417582
#  [8] 0.4316850 0.4040293 0.4362637 0.4478022 0.4368132 0.4219780 0.4265568
# [15] 0.4443223 0.4082418 0.4161172 0.3983516 0.4225275 0.4166667 0.4309524
# [22] 0.4313187 0.4377289 0.4293040

p <- ggplot(do.call(rbind, plotlist), aes(x=cors, color=type)) +
  theme_bw() +
  geom_density() +
  geom_vline(aes(xintercept=cut)) +
  facet_wrap(~region) +
  labs(title="Cell-cell correlations for Lein DFC supertypes", subtitle="Line represents nominal significance (95th percentile of permuted cors)", 
       x="Cell-cell correlation\n(3000 HVGs)") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        text=element_text(size=30),
        plot.subtitle=element_text(size=18))
ggsave(p,file=(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/mislabeling_analyses/cell_cell_cors_top3000hvg.png")),width=14)

########### all genes
corpctg <- list()
plotlist <- list()
for(i in seq_along(cts)){
    cells <- cell_annoall$Cell_ID[cell_annoall[,1] == "VLMC"]
    subset <- cell_exprall[, colnames(cell_exprall) %in% cells]
    subset_var <- apply(subset,1,var)
    subset <- subset[subset_var>0,]

    real_cors <- numpy$corrcoef(subset, rowvar=F)
    diag(real_cors) <- 0
    
    set.seed(i)
    subset_perm <- apply(subset,2,function(x){
        sample(x, length(x), replace=F)
    })
    rand_cors <- numpy$corrcoef(subset_perm, rowvar=F)
    diag(rand_cors) <- 0

    # Calculate cutoffs
    r2 <- rand_cors[upper.tri(rand_cors)]
    norm_cut <- quantile(r2, 0.95)
    r1 <- real_cors[upper.tri(real_cors)]
    corpctg[[i]] <- sum(r1>norm_cut)/length(r1)

    plotlist[[i]] <- data.frame("cors"=c(sample(r1, min(length(r1),10000), replace=F),sample(r2,min(length(r2),10000), replace=F)),
                                "type"=c(rep("real",min(length(r1),10000)), rep("permuted",min(length(r2),10000))),
                                "region"=names(cts)[i],
                                "cut"=norm_cut)

    cat(i, " ")
}

p <- ggplot(do.call(rbind, plotlist), aes(x=cors, color=type)) +
  theme_bw() +
  geom_density() +
  geom_vline(aes(xintercept=cut)) +
  facet_wrap(~region) +
  labs(title="Cell-cell correlations for Lein DFC supertypes", subtitle="Line represents nominal significance (95th percentile of permuted cors)", 
       x="Cell-cell correlation\n(All genes)") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        text=element_text(size=30),
        plot.subtitle=element_text(size=18))
ggsave(p,file=(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/mislabeling_analyses/cell_cell_cors_allgenes.png")),width=14)

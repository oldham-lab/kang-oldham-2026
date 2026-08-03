# Of the modules that are expressed differently in control vs disease, what module genes are driving this?
# DEseq2 documentation: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library(DESeq2)
library(qs)
library(data.table)
library(AnnotationHub)
library(tidyverse)
options(bitmapType = 'cairo')

##########################
# Load objects of interest
##########################

# Load bulk megaset modules and expr
bulk_expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
bulkt <- t(bulk_expr[,-c(1,2)]) |> as.data.frame() |> setNames(bulk_expr[,2])
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/modules/unmerged_modules.qs"))
mod_bc <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv")) |>
  (\(x) tapply(x[, 2], x[, 3], list))()
mod_bc_lengths <- lapply(mod_bc, length) |> unlist()
these_mods <- which(mod_bc_lengths > 3)

# Load control and AD samples, not ComBat processed (DEseq2 does not accept negative values)
expr_SN <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/ROSMAP_samp_filt_all_samples_SampleNetworks/1_10-57-08/ROSMAP_samp_filt_all_samples_1_617_ComBat.csv"))
expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/expr_ROSMAP_raw_counts.csv"), data.table=F)
sif <- fread(data.table=F,file=file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_ROSMAP.csv"))
sif <- sif[sif$excludeReason=="",]
# filter expr to samples in sif
samp_names <- gsub("X", "", colnames(expr))
expr <- expr[,c(1, which(samp_names %in% sif$specimenID))]
samp_names <- samp_names[samp_names %in% sif$specimenID]
# order sif to match expr
sif <- sif[sif$specimenID %in% samp_names,]
sif <- sif[match(samp_names, sif$specimenID),]
sif[,2] <- colnames(expr)[2:ncol(expr)]
# remove sample 319 (unique library batch)
expr <- expr[-which(colnames(expr)==sif[319,2])]
sif <- sif[-319,]
sif[,24] <- as.factor(sif[,24])
# Add gene symbol to expr matrix
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 104))
ahEdb <- ahDb[[1]]
k <- keys(ahEdb, keytype = "TXNAME")
ahEdbgene <- select(ahEdb, k, c("GENEID", "GENENAME"), "TXNAME")
genes <- ahEdbgene[ahEdbgene[,2] %in% expr[,1],c(2,3)] |> 
  dplyr::filter(!duplicated(GENEID), GENENAME != "", !duplicated(GENENAME))
genes <- genes[order(genes[,1]),]
expr <- expr[expr[,1] %in% genes[,1],]
expr <- expr[match(genes[,1], expr[,1]),]
expr <- cbind(genes, expr[,2:ncol(expr)])
colnames(expr)[1:2] <- c("ensembl_id", "Gene")
sif_clin <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_clinical.csv"))
sif2 <- dplyr::left_join(sif, sif_clin, by = dplyr::join_by("individualID")) |> 
  dplyr::mutate(alzdx = NA)
sif2$alzdx[sif2$cogdx == 4] <- "AD"
sif2$alzdx[sif2$cogdx == 1] <- "Con"
sif2 <- sif2 |> dplyr::filter(!is.na(alzdx))

# Prepare for calculating DE genes using DEseq2
expr2 <- expr[,-c(1:2)]
rownames(expr2) <- expr$Gene
expr2 <- expr2[,colnames(expr2) %in% sif2$specimenID]
expr2 <- expr2[,colnames(expr2) %in% colnames(expr_SN)] # filter to samples kept after SampleNetwork filtering
sif2 <- sif2[sif2$specimenID %in% colnames(expr2),]
rownames(sif2) <- sif2$specimenID
sif2 <- sif2 |> dplyr::select(alzdx, libraryBatch)
sif2$alzdx <- as.factor(sif2$alzdx)
sif2$libraryBatch <- as.factor(sif2$libraryBatch)

# Not super kosher but DESeqDatSetFromMatrix does not accept non-integers 
# would have to create a txi object (tximport) to use kallisto data as is
expr3 <- apply(expr2,2,floor)

# Split and transpose rosmap expr for correlations
expr3t <- t(expr3)
expr3tad <- expr3 |> as.data.frame() |> dplyr::select(rownames(sif2)[sif2$alzdx == "AD"]) |> t()
expr3tcon <- expr3 |> as.data.frame() |> dplyr::select(rownames(sif2)[sif2$alzdx == "Con"]) |> t()

# Load copa_compare data
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_negative.qs")))

sea_means <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/sn_summary_tables/allmlist_log.qs"))
sub_diff <- sea_means[[1]] - sea_means[[2]]
seade <- fread(data.table=F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/DE/DE_gene_list_dx_blockSubclass_allgenes_p0.05.csv"))
seapval <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/p_values_all_Subclass.csv"))
seadist <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_distances_all_Subclass.csv"))
rownames(seadist) <- these_mods
rownames(seapval) <- these_mods
cell_anno_pb <- fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/annotations_by_donor_sum_subclass.csv"))
desubclasses <- unique(cell_anno_pb$Subclass)


## Calculate DE

# dds <- DESeqDataSetFromMatrix(countData = expr3,
#                               colData = sif2,
#                               design= ~ alzdx + libraryBatch)
# dds <- DESeq(dds)
# res <- results(dds, name="alzdx_Con_vs_AD")
# qsave(res, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "rosmap_bulk/DEseq2_results_cogdx1vs4.qs"))
# or to shrink log fold changes association with condition (for visualization):
# res2 <- lfcShrink(dds, coef="alzdx_Con_vs_AD", type="apeglm") # for visualization
res <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "rosmap_bulk/DEseq2_results_cogdx1vs4.qs"))
resdf <- as.data.frame(res@listData) |>
  `rownames<-`(res@rownames)
#resdf |> dplyr::arrange(padj) |> head()
siggenes <- rownames(resdf)[resdf$padj < 0.05]
#> length(siggenes)
#[1] 5396
# 5396 genes are differentially expressed between control and AD (cogdx=1 vs cogdx=4)


# How many DE genes are represented in module genes (bc)?
# sum(siggenes %in% unlist(mod_bc)) / length(siggenes)
# [1] 0.7813195

# How many mods contain DE genes?
# has <- lapply(mod_bc, \(x) sum(x %in% siggenes)) |> unlist()
# sum(has > 0)
# 943 mods contain at least 1 DE gene
# sum(has > 2)
# 543

# How many mods are differentially expressed?
# Create module "samples" by averaging module genes

expr3_mod <- lapply(mod_bc, \(x){
  if(length(x) > 1){
    return(colMeans(expr3[rownames(expr3) %in% x, ]))
  } else {
    return(expr3[rownames(expr3) %in% x, ])
  }
})
expr3_mod <- do.call(rbind, expr3_mod)
expr3_mod <- apply(expr3_mod, 2, floor)

# dds_mod <- DESeqDataSetFromMatrix(countData = expr3_mod,
#                                   colData = sif2,
#                                   design= ~ alzdx + libraryBatch)
# dds_mod <- DESeq(dds_mod)
# res_mod <- results(dds_mod, name="alzdx_Con_vs_AD")
# qsave(res_mod, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "rosmap_bulk/DEseq2_results_cogdx1vs4_modLevel.qs"))
res_mod <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "rosmap_bulk/DEseq2_results_cogdx1vs4_modLevel.qs"))
resdf_mod <- as.data.frame(res_mod@listData) |>
  `rownames<-`(res_mod@rownames)
resdf_mod |> dplyr::arrange(padj) |> head()
# Log2 fold change ranges between -0.4 and 0.3

# Which mods are "significant" in bulk ROSMAP?
sigmods <- rownames(resdf_mod)[resdf_mod$padj < 0.05] |> as.numeric()
# 406 mods marked as significant
sigmodsneg <- rownames(resdf_mod)[resdf_mod$padj < 0.05 & resdf_mod$log2FoldChange > 0] # higher in con
sigmodspos <- rownames(resdf_mod)[resdf_mod$padj < 0.05 & resdf_mod$log2FoldChange < 0] # higher in ad

# Load significant modules (projections)
euc_all_fdr <- unique(unlist(euc_dist[[1]][[3]]))
# length(euc_all_fdr)
# 287 mods marked as significant by copa_compare (subclass, fdr)
# sum(euc_all_fdr %in% sigmods)
# 120 of those mods are marked as significant by DEseq2
euc_pos_fdr <- unique(unlist(euc_dist[[2]][[3]]))
euc_neg_fdr <- unique(unlist(euc_dist[[3]][[3]]))
# > sum(sigmodsneg %in% euc_pos_fdr)
# [1] 18
# > sum(sigmodsneg %in% euc_neg_fdr)
# [1] 63
# > sum(sigmodspos %in% euc_pos_fdr)
# [1] 31
# > sum(sigmodspos %in% euc_neg_fdr)
# [1] 20
# Generally seems to correspond

######################
# Exploratory analyses
######################

# What % of each module's genes is significantly DE in ROSMAP?
# temp <- lapply(mod_bc, \(mod){
#   df <- resdf[rownames(resdf) %in% mod, ]
#   return(sum(df$pvalue < 0.05, na.rm = T)/nrow(df))
# }) |> unlist()
# summary(temp)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.1667  0.3333  0.3554  0.5000  1.0000 
# sum(temp == 1)
# 34

# temp3 <- lapply(mod_bc[euc_all_fdr], \(mod){
#   df <- resdf[rownames(resdf) %in% mod, ]
#   return(sum(df$pvalue < 0.05, na.rm = T)/nrow(df))
# }) |> unlist()
# summary(temp3)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.2500  0.3750  0.4033  0.5156  1.0000 
# sum(temp3 == 1)
# 3

# temp2 <- lapply(mod_bc[sigmods], \(mod){
#   df <- resdf[rownames(resdf) %in% mod, ]
#   return(sum(df$pvalue < 0.05, na.rm = T)/nrow(df))
# }) |> unlist() 
# summary(temp2)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.4286  0.5522  0.5797  0.7143  1.0000 
# sum(temp2 == 1)
# 32
# Mean % of significantly DE genes per module is 35% over all modules, 40% in copa_compare modules, 58% in ROSMAP DE mods
# comp1 <- c(1:1158)[temp == 1]
# comp2 <- c(1:1158)[sigmods][temp2 == 1]
# length(intersect(comp1, comp2))
# 32
# 32 significant modules have 100% significantly DE genes
# 2 modules are not significant but have 100% significantly DE genes

# For modules marked as significant by DEseq2, how many individual module genes match directionality of overall module?
# consis <- lapply(mod_bc, \(x){
#   higher_in_con <- sum(resdf[rownames(resdf) %in% x, 2] > 0)
#   return(higher_in_con/length(x))
# }) |> unlist()
#> summary(consis)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.2500  0.5455  0.5312  0.8333  1.0000       3
# About 50/50 when looking at all modules
#sum(consis == 1 | consis == 0, na.rm = T)
# 251 mods are completely coherent
#sum(consis == 1, na.rm = T)
# 144 uniformly higher in control
#sum(consis == 0, na.rm = T)
# 107 uniformly higher in AD
# summary(consis[sigmods])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.1057  0.7778  0.5627  0.9444  1.0000
# More genes higher in control, but with some low outliers, when looking at all DE mods
# sum(consis[sigmods] == 1 | consis[sigmods] == 0, na.rm = T)
# 156 mods are completely coherent and significant
# 251-156=95 mods are completely coherent but not significant
# summary(consis[sigmodspos])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.07692 0.09587 0.15035 0.50000 
# As expected most genes higher in AD when looking at mods significantly higher in AD
# Interesting that there's 5 AD mods with ~40-50% mod genes higher in control
#summary(consis[sigmodsneg])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3333  0.8462  0.9231  0.9024  1.0000  1.0000
# > head(sort(consis[sigmodsneg]), 20)                                                    
#       834       457       767       311       318       662       845       931 
# 0.3333333 0.5000000 0.5000000 0.6000000 0.6666667 0.6666667 0.6666667 0.6666667 
#       317       612       929       371       259       856       225       291 
# 0.6875000 0.6875000 0.6875000 0.6923077 0.7037037 0.7200000 0.7272727 0.7333333         
#       434      1136       177       436 
# 0.7391304 0.7391304 0.7500000 0.7500000 
# As expected most genes higher in con when looking at mods significantly higher in con
# Interesting that there's 4 con mods with <66% genes higher in con (one at 33%)

# Analyze the outliers
# Mod 834: higher in con overall but only has 33% genes higher in con
# resdf[rownames(resdf) %in% mod_bc[[834]], ]
#             baseMean log2FoldChange      lfcSE        stat       pvalue         padj
# TSPOAP1 3233.1905234    0.138969225 0.02833691  4.90417771 9.381952e-07 5.047792e-05
# PRKACG     0.2513604   -0.154997133 0.31654165 -0.48965794 6.243760e-01           NA
# GABBR1  1067.8024473   -0.002191446 0.04336142 -0.05053907 9.596928e-01 9.797024e-01
# One highly expressed gene higher in con, one very low expression gene, one moderately expressed gene very slightly higher in AD

# resdf[rownames(resdf) %in% mod_bc[[457]], ]
#         baseMean log2FoldChange      lfcSE       stat       pvalue        padj
# SPCS2   1393.323   -0.007473515 0.01484526 -0.5034277 0.6146635823 0.764280029
# HNRNPA1 6610.648    0.063332953 0.01846439  3.4300044 0.0006035714 0.006025055
# HNRNPD  2514.646   -0.009292088 0.01585556 -0.5860461 0.5578445347 0.720041515
# HNRNPA3 4550.071    0.009410615 0.01413869  0.6655930 0.5056712672 0.680109419
# Largely driven by one high-expression gene (HNRNPA1)

# resdf[rownames(resdf) %in% mod_bc[[767]], ]
# Largely driven by one high-expression gene (S100B)

# resdf[rownames(resdf) %in% mod_bc[[311]], ]
# Largely driven by one high-expression gene (TSPYL1)

#resdf[rownames(resdf) %in% mod_bc[[806]], ]
# Largely driven by one gene but all genes are in same direction, BUT only one gene is significant

# How many modules have genes where all genes are higher in control AND all genes are significant (nominal)?
consis_strict_con <- lapply(mod_bc, \(x){
  higher_in_con <- sum(resdf[rownames(resdf) %in% x, 2] > 0)
  sig <- sum(resdf$padj[rownames(resdf) %in% x] < 0.05, na.rm=T)
  if(higher_in_con == length(x) & sig == length(x)){
    return(T)
  } else {
    return(F)
  }
}) |> unlist()
# sum(consis_strict_con)
#[1] Of 144 mods with all genes higher in con, 10 have all genes significant as well
# which(consis_strict_con)
#  228  280  286  411  437  625  725  763  764  866  960 1008 1057 1067 1134 
#  228  280  286  411  437  625  725  763  764  866  960 1008 1057 1067 1134 
# mod_bc_lengths[which(consis_strict_con)]
#  228  280  286  411  437  625  725  763  764  866  960 1008 1057 1067 1134 
#    1    1    1    1    1    2    1    1    4    1    1    2    4    1    3 
# Most of these 15 are singlets
# Let's look at the exceptions:
# resdf[rownames(resdf) %in% mod_bc[[764]], ]
#         baseMean log2FoldChange      lfcSE     stat       pvalue        padj
# RHOBTB2 2920.710     0.11852303 0.03263731 3.631519 0.0002817575 0.003441209
# PRRT3    860.143     0.08259637 0.02860422 2.887559 0.0038824417 0.022638815
# HR      1028.107     0.11017785 0.03825086 2.880402 0.0039716864 0.023022844
# LYNX1   5603.292     0.08905526 0.03450436 2.580985 0.0098518741 0.044339562
# resdf[rownames(resdf) %in% mod_bc[[1057]], ]
#        baseMean log2FoldChange      lfcSE     stat       pvalue         padj
# GLP2R  50.95098      0.2096318 0.08395521 2.496948 1.252673e-02 5.273834e-02
# ASB2   90.39403      0.2033964 0.06589831 3.086519 2.025150e-03 1.412140e-02
# PCDH8 550.05510      0.2415359 0.08730875 2.766456 5.666915e-03 2.980609e-02
# FFAR4 101.10004      0.2507085 0.04369355 5.737884 9.586681e-09 1.734039e-06

# How many modules have genes where all genes are higher in AD AND all genes are significant?
consis_strict_AD <- lapply(mod_bc, \(x){
  higher_in_AD <- sum(resdf[rownames(resdf) %in% x, 2] < 0)
  sig <- sum(resdf$padj[rownames(resdf) %in% x] < 0.05, na.rm=T)
  if(higher_in_AD == length(x) & sig == length(x)){
    return(T)
  } else {
    return(F)
  }
}) |> unlist()
#sum(consis_strict_AD)
# Of 107 mods with all genes higher in AD, 11 mods are significant as well
# which(consis_strict_AD)
#  168  381  506  563  565  583  593  671  737  741  819  951  989 1024 1029 1033 1058 1111 
#  168  381  506  563  565  583  593  671  737  741  819  951  989 1024 1029 1033 1058 1111 
# mod_bc_lengths[which(consis_strict_AD)]
#  168  381  506  563  565  583  593  671  737  741  819  951  989 1024 1029 1033 1058 1111 
#    1    6    2    3    3    3    9    1    6    9    3    2    2    1    1    2    1    1 
# resdf[rownames(resdf) %in% mod_bc[[381]], ]
#         baseMean log2FoldChange      lfcSE      stat       pvalue         padj
# TF     5451.5265     -0.2324616 0.05120270 -4.540027 5.624705e-06 1.951304e-04
# PLP1  21785.6841     -0.2231637 0.05761462 -3.873386 1.073338e-04 1.710531e-03
# ENPP2  1043.5214     -0.2052983 0.06448420 -3.183699 1.454062e-03 1.124939e-02
# GPR37   663.8323     -0.2783600 0.06534336 -4.259959 2.044642e-05 5.087137e-04
# UGT8    789.2100     -0.2984713 0.05957382 -5.010108 5.439951e-07 3.258207e-05
# THBS2   438.6921     -0.1882696 0.05332945 -3.530312 4.150696e-04 4.548435e-03
# resdf[rownames(resdf) %in% mod_bc[[593]], ]
#           baseMean log2FoldChange      lfcSE      stat       pvalue
# PLEKHB1 10766.4919     -0.2051835 0.04124518 -4.974725 6.534022e-07
# HMG20B    503.6095     -0.1973205 0.04336119 -4.550625 5.348684e-06
# PADI2    2437.3719     -0.3061306 0.04793693 -6.386111 1.701576e-10                                                                                                              
# HEY2      141.9927     -0.1504472 0.04748500 -3.168310 1.533277e-03
# DOCK1    1292.4756     -0.2430116 0.04764842 -5.100097 3.394801e-07
# CSRP1   11348.8635     -0.3099799 0.04790450 -6.470789 9.749276e-11
# ZCCHC24  3399.7607     -0.1746190 0.04946360 -3.530253 4.151629e-04
# APLN      496.9169     -0.5753698 0.07018821 -8.197527 2.453830e-16
# TAX1BP3   605.4237     -0.1285894 0.04238058 -3.034159 2.412072e-03
# resdf[rownames(resdf) %in% mod_bc[[737]], ]
#          baseMean log2FoldChange      lfcSE      stat       pvalue         padj
# WNK1   4790.03578    -0.11911144 0.03487802 -3.415086 6.376197e-04 6.245509e-03
# TJP1   1179.77294    -0.16471428 0.04501906 -3.658767 2.534314e-04 3.214634e-03
# KIF13A  899.53842    -0.12548544 0.02988509 -4.198931 2.681785e-05 6.148051e-04
# NKD1    357.24374    -0.10368875 0.04656895 -2.226564 2.597643e-02 8.963403e-02
# ZEB2   2929.92245    -0.07161567 0.03232733 -2.215329 2.673750e-02 9.163089e-02
# CNTF     11.82111    -0.37717725 0.07209387 -5.231752 1.679110e-07 1.355881e-05
# resdf[rownames(resdf) %in% mod_bc[[741]], ]
#          baseMean log2FoldChange      lfcSE      stat       pvalue         padj
# RAI14    184.1310    -0.15308106 0.04331324 -3.534279 4.088899e-04 4.502348e-03
# TGFB3    526.9555    -0.15173036 0.05298090 -2.863869 4.185013e-03 2.390228e-02
# KIF1C   5429.3169    -0.17515300 0.05093215 -3.438948 5.839796e-04 5.881416e-03
# ANP32B  1388.0703    -0.10341413 0.03504793 -2.950649 3.171072e-03 1.954644e-02
# ANKRD40 2704.4795    -0.21588973 0.04115570 -5.245682 1.557053e-07 1.280181e-05
# MCM7     855.2434    -0.18541746 0.03997582 -4.638240 3.513880e-06 1.378420e-04
# KIF5B   6017.5733    -0.20871676 0.03741375 -5.578612 2.424460e-08 3.045391e-06
# SMARCC1  874.8956    -0.07588522 0.03500848 -2.167624 3.018730e-02 9.976193e-02
# FOXO4   1603.2719    -0.23964125 0.04609182 -5.199214 2.001329e-07 1.551953e-05

# Let's take a step back and categorize modules as follows:
# Two axes can be defined:
# - "Coherence": do all genes move in same direction or not
# - Significance: how many genes are significantly different
# Are non-significant genes worth considering regardless of their direction?
# Most DE modules are coherent, i.e. mod genes move in the same direction. 

# Start with significance first. For DE modules, how many genes are significant (FDR)?
# signif_per <- lapply(mod_bc, \(x){
#   count <- sum(resdf$padj[rownames(resdf) %in% x] < 0.05, na.rm=T)
#   return(count/length(x))
# }) |> unlist()
# > summary(signif_per)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.07143 0.20000 0.24269 0.34582 1.00000
# summary(signif_per[sigmods])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.2913  0.4074  0.4472  0.5714  1.0000 
# Seems like for most DE mods, only 1/3 to 1/2 of the genes are significant
# signif_count <- lapply(mod_bc, \(x){
#   count <- sum(resdf$padj[rownames(resdf) %in% x] < 0.05, na.rm=T)
#   return(count)
# }) |> unlist()
# > summary(signif_count[sigmods])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.000   2.000   5.000   6.123   9.000  29.000 

# Test coherence of only the sig genes (FDR) for each DE mod
sig_coh_con <- lapply(mod_bc, \(x){
  sigs <- which(resdf$padj[rownames(resdf) %in% x] < 0.05)
  mod <- resdf$log2FoldChange[rownames(resdf) %in% x] 
  count_con <- sum(mod[sigs] > 0)
  return(count_con/length(sigs))
}) |> unlist()
# summary(sig_coh_con[sigmods])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.0000  1.0000  0.5731  1.0000  1.0000       6 
# NA values represent modules where genes are not present in ROSMAP
#sum(sig_coh_con[sigmods] == 1, na.rm=T)
# 218
#sum(sig_coh_con[sigmods] == 0, na.rm=T)
# 158
# length(sig_coh_con[sigmods][sig_coh_con[sigmods] < 1 & sig_coh_con[sigmods] > 0])
# For 30 of the 406 sig mods, the sig genes are not completely coherent (for the rest the sig genes are completely coherent)
# summary(sig_coh_con[sigmods][sig_coh_con[sigmods] < 1 & sig_coh_con[sigmods] > 0])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.05882 0.10833 0.26667 0.46872 0.85337 0.95454       6 
# sig_coh_con[sigmods][sig_coh_con[sigmods] < 1 & sig_coh_con[sigmods] > 0]
#         26       <NA>        152        165        190        303        367 
# 0.95454545         NA 0.77777778 0.05882353 0.84615385 0.06666667 0.91666667 
#        383       <NA>        429        434        439        468        499 
# 0.92307692         NA 0.83333333 0.80000000 0.90000000 0.20000000 0.10000000 
#        560        568        585       <NA>        729       <NA>        752 
# 0.16666667 0.87500000 0.06666667         NA 0.09090909         NA 0.11111111 
#        804       <NA>        871        930       <NA>        992       1007 
# 0.08333333         NA 0.80000000 0.12500000         NA 0.33333333 0.12500000 
#       1011       1082 
# 0.16666667 0.92857143 
# Coherence of these 30 modules is high. Mostly one or two genes go against the grain.
# summary(unlist(lapply(mod_bc[sigmods][sig_coh_con[sigmods] < 1 & sig_coh_con[sigmods] > 0], length)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0.00   13.25   20.00   21.03   29.50   59.00 
# These mods are mostly pretty large

# resdf[rownames(resdf) %in% mod_bc[[26]] & resdf$padj < 0.05, ]
# Mostly coherent (22 genes higher in control) with 1 exception (PPP1R1A)
# expr3t |> as.data.frame() |> dplyr::select(mod_bc[[26]]) |> cor() |> 
#   (\(x) apply((x + 1) / 2, 2, sum))() |> sort()
# expr3tad |> as.data.frame() |> dplyr::select(mod_bc[[26]]) |> cor() |> 
#   (\(x) apply((x + 1) / 2, 2, sum))() |> sort()
# expr3tcon |> as.data.frame() |> dplyr::select(mod_bc[[26]]) |> cor() |> 
#   (\(x) apply((x + 1) / 2, 2, sum))() |> sort()
# Seems that PPP1R1A has the lowest connectivity across all samples and in AD samples but is second lowest ahead of WFDC2 in control samples

# > resdf[rownames(resdf) %in% mod_bc[[152]] & resdf$padj < 0.05, ]
#            baseMean log2FoldChange      lfcSE      stat       pvalue
# TMEM38A   898.21750     0.13573550 0.03495827  3.882786 1.032664e-04
# OLFM1   11420.39077     0.11278009 0.03559870  3.168096 1.534409e-03
# CYP1A1     11.49024    -0.31508807 0.09057433 -3.478779 5.037039e-04
# ENSA     3741.31272     0.06721965 0.02580942  2.604462 9.201875e-03
# FNDC5     488.58311     0.17558351 0.03741749  4.692552 2.698181e-06
# GAP43    7436.87227     0.16215508 0.04453814  3.640814 2.717777e-04
# RIIAD1     43.93498     0.21465188 0.05172253  4.150065 3.323805e-05
# NKX1-2     14.90145    -0.19663988 0.07510874 -2.618069 8.842885e-03
# LYPD8      23.29614     0.19732247 0.07472270  2.640730 8.272767e-03

# In short, most of the 406 DE modules are driven by a majority (in some cases all) of the genes. 
# What do "DE modules" really represent to begin with?
# Not sure many of the modules offer any real information over just a traditional DE analysis...
# - ...with the exception of the completely coherent modules

## What does coherence look like for copa_compare mods?
# lapply(euc_dist$pos[[3]]$Astrocyte, \(mod){
#   temp <- sub_diff$Astrocyte[rownames(sub_diff) %in% mod_bc[[mod]]]
#   return(sum(temp > 0) / length(temp))
# }) |> unlist() 
#       329       348       356       402       421       474       507       514 
# 0.7500000 1.0000000 0.9250000 0.6333333 0.6153846 1.0000000 0.8333333 0.9047619 
#       540       566       593       606       646       693       786       792 
# 0.7500000 0.7619048 0.8888889 0.9090909 0.9000000 1.0000000 0.5384615 0.7777778 
#       804       849       874       877       894      1031      1035      1037 
# 0.7777778 0.9090909 0.8947368 0.6153846 0.8750000 0.5909091 0.8333333 0.8636364 
#      1125 
# 0.7142857 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.5385  0.7500  0.8333  0.8105  0.9048  1.0000 
# Of the modules that show higher expression in AD astrocytes, at least 50% of genes are coherent (median 83%)

# lapply(euc_dist$neg[[3]]$Astrocyte, \(mod){
#   temp <- sub_diff$Astrocyte[rownames(sub_diff) %in% mod_bc[[mod]]]
#   return(sum(temp < 0) / length(temp))
# }) |> unlist() 
#        20        25        83        85       104       112       119       131 
# 0.8400000 1.0000000 0.8235294 1.0000000 0.7894737 0.8684211 0.8214286 0.9166667 
#       149       156       164       192       198       260       295       301 
# 0.9375000 0.7037037 0.8400000 0.9615385 0.7777778 0.8888889 1.0000000 0.9411765 
#       376       382       386       686       709       716       724       767 
# 0.9411765 0.7857143 1.0000000 0.8571429 0.8888889 1.0000000 0.7058824 0.6666667 
#       812       948       972      1049      1100      1146 
# 0.7142857 0.8333333 0.8000000 0.7142857 0.5937500 0.9090909
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.5938  0.7867  0.8486  0.8507  0.9403  1.0000 
# Same situation with the modules with higher expression in Con: majority are coherent
# Module 1100 stands out as an outlier



################################
# Coherence in SEAAD2024 modules
################################

# How many modules per subclass are completely coherent (in SEAAD data)?
# Positive modules (higher in AD):
copacohpos <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
    return(sum(temp > 0) / length(temp))
  }) |> unlist() 
  return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
}, euc_dist$pos[[3]], names(euc_dist$pos[[3]]), SIMPLIFY = F) #|> do.call(what = "rbind")# |>
#   ggplot(aes(x = subclass, y = pcnt)) + 
#     theme_bw() + geom_violin() + coord_flip() + theme(text = element_text(size = 12))
#lapply(copacohpos, \(x) sum(x$pcnt == 1)) |> unlist()
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#               4               1               8               1               0 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               1               0               0               0               1 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#               0               1               0               0               0 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#               2               0               0               0               0 
#             Sst       Sst Chodl             Vip            VLMC 
#               1               0               1               0 
# Anywhere from 0-4 modules are completely coherent and higher in con for each subclass.
# How many overlap with mods significant in ROSMAP?
#test1 <- lapply(copacohpos, \(x) rownames(x)[x$pcnt==1]) |> unlist() |> as.numeric()
#sum(test1 %in% which(consis_strict_con))
# 0
# How many of these modules are specific to a subclass
copacohpos_complete <- lapply(copacohpos, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
tapply(copacohpos_complete[ ,2], copacohpos_complete[ ,1], list)
# $`1055`                                                                                                                                                         
# [1] "Endothelial"                                                                                                                               
# $`13`                                                                                                                                                           
# [1] "Endothelial"                                                                                                                                                                                                                                                                                                             
# $`149`                                                                                                                                                          
# [1] "Endothelial"     "Oligodendrocyte"                                                                                                                                                                                        
# $`181`                                                                          
# [1] "L2/3 IT"                                                                                                                                                                                                                                  
# $`182`                                                                                                                                                          
# [1] "Endothelial"                                                                                                       
# $`348`                                                                          
# [1] "Astrocyte"                                                                                                        
# $`474`                               
# [1] "Astrocyte"
# $`514`
# [1] "Chandelier" "L5 ET"      "L6 IT"      "L6b"        "Sst"       
# [6] "Vip"       
# $`540`
# [1] "Oligodendrocyte"
# $`693`
# [1] "Astrocyte"
# $`716`
# [1] "Endothelial"
# $`882`
# [1] "Endothelial"
# $`947`
# [1] "Endothelial"
# $`948`
# [1] "Endothelial"
# $`970`
# [1] "Astrocyte"
# The majority of modules are only completely coherent in one subclass.
# 4 modules are completely coherent in two subclasses
# 1 module (514) is completely coherent in 6 subclasses - this is a mitochondrial module

#  modules (higher in con):
copacohneg <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], colnames(sub_diff) == subclass_name]
    return(sum(temp < 0) / length(temp))
  }) |> unlist() 
  return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
}, euc_dist$neg[[3]], names(euc_dist$neg[[3]]), SIMPLIFY = F) 
lapply(copacohneg, \(x) sum(x$pcnt == 1)) |> unlist()
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#               5               0               3               5              13 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               0              10              17               2               9 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#               3               2              24               2               0 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#               1               0               2               6               1 
#             Sst       Sst Chodl             Vip            VLMC 
#              20               4               5               0 
# More coherent modules that are higher in control than there are in AD.
# - Around 4-6 on average with max of 23 (Lamp5) are completely coherent and higher in AD.
# How many overlap with mods significant in ROSMAP?
#test1 <- lapply(copacohneg, \(x) rownames(x)[x$pcnt==1]) |> unlist() |> as.numeric()
#sum(test1 %in% which(consis_strict_AD))
# 0 
# How many of these modules are specific to a single subclass?
copacohneg_complete <- lapply(copacohneg, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
copacohneg_count <- tapply(copacohneg_complete[ ,2], copacohneg_complete[ ,1], list) |>
  lapply(length) |> unlist()
sum(copacohneg_count > 1)
# 26
sum(copacohneg_count == 1)
# 31
# 31 of the 57 completely coherent modules higher in AD are specific to a single subclass
summary(copacohneg_count)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.000   1.000   1.000   2.351   2.000  16.000 
# mod 149 is an outlier that is coherent in 16 subclasses (all neuronal)
# - this module is a ribosomal module that is higher in AD in all 
#test <- tapply(copacohneg_complete[ ,2], copacohneg_complete[ ,1], list) 
#names(euc_dist$neg[[3]])[!names(euc_dist$neg[[3]]) %in% test$`149`]
# "Astrocyte"       "Chandelier"      "Endothelial"     "L5 ET"          
# [5] "Microglia-PVM"   "Oligodendrocyte" "OPC"             "VLMC"  
# mod 149 is not present in glial celltypes and L5ET
# Other mods present in many subclasses:
# mod 182 (9 subclasses, neuronal): electron transport related
# mod 716 (9 subclasses, neuronal + astro): ribosomal related
# mod 295 (7 subclasses, neuronal + astro): ribosomal related
# mod 25 (7 subclasses, neuronal + astro): ETC related
# mod 416 (6 subclasses, neuronal): not sure
# mod 85 (6 subclasses, neuronal + astro): ribosomal
# mod 131 (5 subclasses, excitatory + VIP): no idea
# mod 441 (5 subclasses, neuronal): tubulins
# mod 451 (5 subclasses, excitatory): no idea (protein/AA modification, some mito-related?)

# For majority of subclasses, coherence in AD modules ranges from 50% to 100%, median roughly 80-90% ish
summary(mod_bc_lengths[euc_pos_fdr])
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  4.00   12.00   18.00   20.13   26.00   53.00 
# All the AD modules fairly large modules as well

# Can we attach significance to this coherence analysis?
# Load DE genes (AD vs control) for SEAAD2024 for each individual subclass
delistsea <- list.files(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/DE_by_subclass_fast/data/"), full.names = T) |>
  (\(.) grep(".RDS", ., value=T))() |> 
  (\(.) grep("gene_list", ., value=T))() |> 
  lapply(readRDS) 
# lapply(delistsea, \(x) sum(x[[1]]$FDR < 0.05)) |> unlist() |> setNames(desubclasses)
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#               8               2               0               1               0 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               0               2               3               1               0 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#               0               0               0               1               0 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#               0               1               0               2               1 
#             Sst       Sst Chodl             Vip            VLMC 
#               0               1               0               0 
# there are very few differentially expressed genes for any subclass (FDR cutoff)
#astrode <- delistsea[[1]][[1]]$genes[delistsea[[1]][[1]]$FDR < 0.05] # 8 genes total for astrocyte
#sum(astrode %in% mod_bc[euc_dist$all[[3]][[1]]] |> unlist())
# 0 of the 8 DE Astro genes show up in any of the modules. 

# Let's look at astrocyte modules marked by copa_compare:
# test <- lapply(mod_bc[euc_dist$all[[3]][[1]]], \(x){
#   delistsea[[1]][[1]][delistsea[[1]][[1]]$genes %in% x, ]
# })
#sum(test$`1100`$logFC > 0)
# 15
#mean(test$`1100`$logFC[test$`1100`$logFC > 0])
# [1] 0.2852956
#mean(test$`1100`$logCPM[test$`1100`$logFC > 0])
# [1] 2.093802
#sum(test$`1100`$logFC < 0)
# 14
#mean(test$`1100`$logFC[test$`1100`$logFC < 0])
# [1] -0.3376981
#mean(test$`1100`$logCPM[test$`1100`$logFC < 0])
# [1] 1.595362
#sub_diff$Astrocyte[rownames(sub_diff) %in% mod_bc[[1100]]] |> setNames(rownames(sub_diff)[rownames(sub_diff) %in% mod_bc[[1100]]]) |> sort()
#        HSPA1B        HSPA1A         HSPB1          BAG3      SERPINH1 
# -1.070043e-01 -1.000480e-01 -9.605719e-02 -4.193802e-02 -3.597231e-02 
#        TENT5A       C7orf61           CKM           LBP         ESPNL 
# -4.907107e-03 -4.652720e-03 -2.030609e-03 -1.706353e-03 -4.073580e-04 
#         CCL26         TPRX1      CATSPER4        NPC1L1        MOGAT2 
# -2.775541e-04 -2.564176e-04 -2.425647e-04 -1.610975e-04 -7.690338e-05 
#           LTA           LHB       LGALS9C         APOC3        SPACA4 
# -6.467522e-05 -4.510484e-05 -3.592024e-05 -3.196762e-05  0.000000e+00 
#        GAGE2A         AGTR2         IL17C          CRNN          VHLL 
#  0.000000e+00  0.000000e+00  2.303125e-05  6.718434e-05  1.055974e-04 
#         LCN10      PLA2G12B          TLR9        NANOS2          OIT3 
#  1.645277e-04  2.112246e-04  2.241259e-04  3.275207e-04  3.441542e-04 
#         FOXJ1         MKNK2 
#  2.618733e-03  1.909237e-02 
#test$`1100` |> arrange(logFC)
# Load p-values for astrocyte modules
#astrosig <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/p_values_all_Subclass.csv"))
#astrosig$Astrocyte[these_mods %in% euc_dist$all[[3]][[1]]]
# The copa_compare pvalues for all mods are very similar.

# None of the copa_compare module genes are significant at FDR cutoff.
# Majority of modules seem mostly coherent, with a few notable exceptions.
# For example, module 1100 is all over the place.
# - Mod 1100 is higher in control (according to copa_compare)
# - The # of mod genes higher in control vs ad is about half half
# - the genes higher in con are higher expression but lower fold change than genes higher in AD

## How many copa_compare module genes are significant DE genes in ROSMAP data?
# lapply(mod_bc[euc_all_fdr], \(mod){
#   df <- resdf[rownames(resdf) %in% mod, ]
#   return(sum(df$padj < 0.05, na.rm = T)/nrow(df))
# }) |> unlist() |> summary()
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.1250  0.2414  0.2791  0.3750  1.0000 
# The vast majority of copa_compare modules do not have more than 40% of genes that are significantly DE in ROSMAP

## How coherent are significant copa_compare module genes in ROSMAP data?
copa_all_coh <- lapply(mod_bc[euc_pos_fdr], \(x){
  sigs <- which(resdf$padj[rownames(resdf) %in% x] < 0.05)
  mod <- resdf$log2FoldChange[rownames(resdf) %in% x] 
  count_con <- sum(mod[sigs] > 0)
  return(count_con/length(sigs))
}) |> unlist() 
#summary(copa_all_coh)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  # 0.000   0.000   0.000   0.354   1.000   1.000      13 
#sum(copa_all_coh == 0, na.rm = T)
# 56
#sum(copa_all_coh > 0 & copa_all_coh < 1, na.rm = T)
# 21 
#sum(copa_all_coh == 1, na.rm = T)
# 27
#sum(is.na(copa_all_coh))
# 13, # NA indicates that none of the genes are significant (or none of the module genes are present in ROSMAP)
copa_all_neg_coh <- lapply(mod_bc[euc_neg_fdr], \(x){
  sigs <- which(resdf$padj[rownames(resdf) %in% x] < 0.05)
  mod <- resdf$log2FoldChange[rownames(resdf) %in% x] 
  count_con <- sum(mod[sigs] > 0)
  return(count_con/length(sigs))
}) |> unlist() 
#sum(copa_all_neg_coh == 0, na.rm = T)
# 47
#sum(copa_all_neg_coh > 0 & copa_all_neg_coh < 1, na.rm = T)
# 43
#sum(copa_all_neg_coh == 1, na.rm = T)
# 96
#sum(is.na(copa_all_neg_coh), na.rm = T)
# 14
#bla <- resdf[rownames(resdf) %in% mod_bc[[48]],]
# Of the copa_compare module genes that are significantly DE, the majority (>50%) are completely coherent in the right direction, but ~20% are 
#  completely coherent in the wrong direction and ~20% are not coherent (according to ROSMAP bulk data). 
# But keep in mind that most of the module genes are not significant. 

## How does coherence look if we consider all genes and not just the significant ones? (considering that there are not many significant genes to begin with)
copa_coh <- lapply(mod_bc[euc_all_fdr], \(x){
  mod <- resdf$log2FoldChange[rownames(resdf) %in% x] 
  count_con <- sum(mod > 0)
  return(count_con/length(mod))
}) |> unlist() 
# > summary(copa_coh)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.1300  0.3333  0.4407  0.7948  1.0000       1 
#sum(copa_coh == 0, na.rm = T)
# 18 higher in con
#sum(copa_coh == 1, na.rm = T)
# 28 higher in AD
#sum(copa_coh > 0 & copa_coh < 1, na.rm = T)
# 240
# p <- data.frame("pcnt" = copa_coh[!is.na(copa_coh)]) |>
#   ggplot(aes(x = pcnt)) + geom_density()
# ggsave(p, file = "~/test/test.png")
# The distribution of percentages is notably bimodal.
#sum(copa_coh > 0.15 & copa_coh < .85, na.rm = T)
# 147

# Do completely coherent mods in SEAAD overlap with those in ROSMAP?
seapos <- lapply(copacohpos, \(x) rownames(x)[x$pcnt == 1]) |> unlist()
seapos <- seapos[!duplicated(seapos)]
#length(seapos)
# 20 completely coherent mods higher in AD total (but some are significant for multiple subclasses)
seaneg <- lapply(copacohneg, \(x) rownames(x)[x$pcnt == 1]) |> unlist()
seaneg <- seaneg[!duplicated(seaneg)]
#length(seaneg)
# 66 completely coherent mods higher in con total (but some are significant for multiple subclasses)

rospos <- names(copa_coh)[copa_coh == 0]
rosneg <- names(copa_coh)[copa_coh == 1]
#sum(seapos %in% rospos)
# 3
#sum(seapos %in% rosneg)
# 2
#sum(seaneg %in% rospos)
# 5
#sum(seaneg %in% rosneg)
# 6
# Overall a bit of overlap, but mostly no overlap
# meaning that for most copa_compare modules, the mod genes are not 100% coherent across all cells. they are only completely coherent for certain subclasses

#####################################
# Differential co-expression analysis
#####################################

# Load ROSMAP Con-AD simMat
simMatAD <- qread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data_other/ROSMAP_bulk/adjmatAD.qs")) * 2 - 1
simMatCon <- qread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data_other/ROSMAP_bulk/adjmatCon.qs")) *2 - 1
#summary(simMatAD[upper.tri(simMatAD)])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.72667  0.07124  0.24100  0.28396  0.48335  1.00000 
#summary(simMatCon[upper.tri(simMatCon)])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.71039  0.05681  0.21407  0.26492  0.45907  1.00000
# The correlations in AD samples are slightly higher than in con samples. 

# Compare correlations:
# Method 1: wilcox test between cors in normal vs cors in AD
wilcp <- lapply(mod_bc, \(mod){
  if(length(mod) > 1){
    these <- colnames(simMatAD) %in% mod
    temp1 <- simMatAD[these, these]
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- simMatCon[these, these]
    temp2 <- temp2[upper.tri(temp2)]
    out <- wilcox.test(temp1, temp2, alternative = "greater")
    out2 <- wilcox.test(temp1, temp2, alternative = "less")
    return(data.frame("pval_AD" = out$p.value, 
                      "pval_con" = out2$p.value))
  } else {
    return(data.frame("pval_AD" = 1, "pval_con" = 1))
  }
}) |> do.call(what = "rbind")
#wilcp |> arrange(pval_AD) |> head()
#sum(wilcp$pval_AD < 0.05/1158)
# 34 modules are significant (significantly differentially co-expressed); cors are higher in AD
#sum(wilcp$pval_con < 0.05/1158)
# 16 modules are significant; cors are higher in con
wilcp_sig <- wilcp |> dplyr::filter(wilcp$pval_AD < 0.05/1158 | wilcp$pval_con < 0.05/1158)

# Is there overlap between the differentially co-expressed modules and the AD ones?
a1 <- names(consis_strict_AD)[consis_strict_AD] 
c1 <- names(consis_strict_con)[consis_strict_con] 
#sum(rownames(wilcp_sig) %in% a1)
# none
#sum(rownames(wilcp_sig) %in% c1)
# none

# How does directionality work with correlations?
direccor <- lapply(mod_bc, \(mod){
  if(length(mod) > 1){
    these <- colnames(simMatAD) %in% mod
    temp1 <- simMatAD[these, these]
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- simMatCon[these, these]
    temp2 <- temp2[upper.tri(temp2)]
    out <- data.frame("pval_AD" = temp1, 
                      "pval_con" = temp2) |>
      dplyr::mutate(diff = sign(pval_AD - pval_con))
    return(sum(out$diff == 1) / nrow(out))
  } else {
    return(NA)
  }
}) |> unlist()
allADcor <- which(direccor == 1)
#length(allADcor)
# 49 mods have pairwise cors that are all higher in AD
#sum(allADcor %in% as.numeric(rownames(wilcp_sig)))
# 0, none of the 25 are significantly different
#allconcor <- which(direccor == 0)
#length(allconcor)
# 25 mods have pairwise cors that are all higher in AD
#sum(allconcor %in% as.numeric(rownames(wilcp_sig)))
# 0, none of the 49 are significantly different
#summary(direccor[as.numeric(rownames(wilcp_sig))])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01515 0.42253 0.70544 0.60691 0.78788 0.90476 
# Most of the significant modules have some pairwise cors that move up and some that move down
# - suggesting that there are a subset of pairs that are responsible for driving significance
# What do the completely concordant differentially co-expressed modules look like?
#allADcor[allADcor %in% as.numeric(a1)]
# 819  951  989 1033 
# there are 4 module that has higher pairwise cors in AD and is coherent for AD
#allADcor[allADcor %in% as.numeric(c1)]
# 1008 1134
#allconcor[allconcor %in% as.numeric(c1)]
# 625
#allconcor[allconcor %in% as.numeric(a1)]
# 506
# These are all very small mods (<=3 mod genes)
# seems that complete coherence for pairwise cors is not a viable direction.

#summary(direccor[as.numeric(rownames(wilcp_sig[wilcp_sig$sign == 1, ]))])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.6287  0.7044  0.7812  0.7730  0.8208  0.9048 
#summary(direccor[as.numeric(rownames(wilcp_sig[wilcp_sig$sign == -1, ]))])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01515 0.09518 0.24544 0.23907 0.38742 0.44156 
# For significant modules, the directionality is fairly strong for each direction (~75% of pairwise cors are coherent)

#summary(mod_bc_lengths[as.numeric(rownames(wilcp_sig))])
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 12.00   29.00   40.00   41.47   51.00   99.00 
# The significant modules are also fairly large.
#names(copa_coh)[copa_coh == 1]

#########
# Summary
#########
# Intuitively understanding AD-associated modules:
# - Modules where all genes have higher/lower expression in disease
#   - completely coherent
# - Modules where pairwise correlations are significantly lower in disease
#   - probably not completely coherent
# - Modules that are not implicated at bulk level but are at single-cell level
#   - can do an expression-based coherence analysis but correlation-based analysis probably not reliable due to nature of SN data

# Goal: choose candidate AD-associated modules based on several criteria:
# - Coherence of differential expression in ROSMAP bulk
# - Coherence of differential expression in SN data (COPA_compare)
# - Differential co-expression of module genes in bulk vs AD
#   - using ROSMAP bulk, which module pairwise cors are lower as a whole in AD? 

# Of the 1158 bulk megaset greedy march modules, 
#  - 33 mods consist entirely of significantly DE genes in the same direction (bulk ROSMAP),
#    - 10 modules are completely coherent and significantly higher in con
#    - 11 modules are completely coherent and significantly higher in AD
#  - 50 modules are significantly differentially co-expressed (bonf, bulk ROSMAP),
#    - 34 with higher pairwise cors in AD,
#    - 16 with higher pairwise cors in con.  
#    - These do not overlap at all with the 21 modules that are completely coherent w.r.t. DE genes
#  - 287 mods marked as significant by COPA_compare 
#    - Coherence is somewhat less preserved when looking at COPA_compare modules in SEAAD2024 data, but still avg ~80%
#    - 0-4 modules per subclass are completely coherent (higher in AD)
#    - 4-6 modules per subclass are completely coherent (higher in CON)
#    - can't attach significance to this as most genes are not significantly DE per subclass

#####################
# Outline of figures:
#####################

# - DE genes: % of significant mod genes vs logFC?
plotdf <- data.frame(
  "mean_logFC" = lapply(mod_bc, \(mod){
    mean_logFC <- mean(resdf$log2FoldChange[rownames(resdf) %in% mod])
    return(mean_logFC)
    }) |> unlist(),
  "pcnt_sig" = lapply(mod_bc, \(mod){
    sigcount <- sum(resdf$pvalue[rownames(resdf) %in% mod] < 0.05)
    return(sigcount / length(mod) * 100)
    }) |> unlist(),
  "mod_size" = mod_bc_lengths,
  "colorvec" = "black"
)
plotdf$colorvec[consis_strict_con] <- "blue"
plotdf$colorvec[consis_strict_AD] <- "red"

p <- ggplot(plotdf, aes(x = mean_logFC, y = pcnt_sig, color = colorvec, alpha = colorvec, size = mod_size)) +
  theme_classic() +
  geom_point() +
  theme(text = element_text(size = 20),
        legend.position = "right",
        legend.margin = margin(-5, 0, 0, 0),
        legend.title = element_text(size = 14, angle = -90, hjust = 0.5),
        legend.text = element_text(size = 10, margin = margin(0, 0, 0, -0.1))) +
  labs(size = bquote(italic("Module size")), 
       x = bquote("Mean log"[2]~"FC of mod genes"),
       y = "% of mod genes that are DE") +
  scale_size(range = c(0.1, 5), breaks = c(1, 25, 50, 75, 99)) +
  scale_color_manual(values = c("black" = "black", "red" = scales::hue_pal()(3)[1], "blue" = scales::hue_pal()(3)[3])) +
  scale_alpha_manual(values = c(0.1, 1, 1)) +
  guides(color = "none", alpha = "none", size = guide_legend(title.position = "right"))
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_A.png"), bg = "white", width = 6, height = 5)

# What is relationship between logFC and coherence
plotdf2 <- data.frame(
  "mean_logFC" = lapply(mod_bc, \(mod){
    mean_logFC <- mean(resdf$log2FoldChange[rownames(resdf) %in% mod])
    return(mean_logFC)
    }) |> unlist(),
  "pcnt_coh" = lapply(mod_bc, \(mod){
    sigcount <- sum(resdf$log2FoldChange[rownames(resdf) %in% mod] > 0)
    return(sigcount / length(mod) * 100)
    }) |> unlist(),
  "mod_size" = mod_bc_lengths,
  "colorvec" = "black"
)
plotdf2$colorvec[consis_strict_con] <- "blue"
plotdf2$colorvec[consis_strict_AD] <- "red"

p2cor <- cor(plotdf2[, 1], plotdf2[, 2], use = "pairwise.complete.obs") |> signif(2)
p2 <- ggplot(plotdf2, aes(x = mean_logFC, y = pcnt_coh, color = colorvec, alpha = colorvec, size = mod_size)) +
  theme_classic() +
  geom_point() +
  theme(text = element_text(size = 20),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, -10),
        legend.title = element_text(size = 14, angle = -90, hjust = 0.5),
        legend.text = element_text(size = 10, margin = margin(0, 0, 0, -0.1)), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16)) +
  labs(size = bquote(italic("Module size")), 
       x = bquote("Mean log"[2]~"FC of mod genes"),
       y = bquote("% of mod genes where log"[2]~"FC > 0")) +
  scale_size(range = c(0.1, 5), breaks = c(1, 25, 50, 75, 99)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  annotate("text", label = bquote(italic(r)~"="~.(p2cor)), x = -0.25, y = 95, size = 6) +
  scale_color_manual(values = c("black" = "black", "red" = scales::hue_pal()(3)[1], "blue" = scales::hue_pal()(3)[3])) +
  scale_alpha_manual(values = c(0.05, 1, 1)) +
  guides(color = "none", alpha = "none", size = guide_legend(title.position = "right"))
ggsave(p2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_test.png"), bg = "white", width = 5, height = 5)
 
# Mean correlations in con vs AD
meanmodcors1 <- lapply(mod_bc, \(mod){
  these <- colnames(simMatAD) %in% mod
  m1 <- simMatAD[these, these]
  m2 <- simMatCon[these, these]
  return(data.frame("AD" = mean(m1[upper.tri(m1)]), "Con" = mean(m2[upper.tri(m2)])))
}) |> do.call(what = "rbind") |>
  mutate("mod_size" = mod_bc_lengths,
         "colorvec" = "black") 
meanmodcors1$colorvec[wilcp$pval_con < 0.05 / 1158] <- "blue"
meanmodcors1$colorvec[wilcp$pval_AD < 0.05 / 1158] <- "red"

p3 <- ggplot(meanmodcors1, aes(x = Con, y = AD, size = mod_size, color = colorvec, alpha = colorvec)) + 
  theme_classic() + 
  geom_point() + 
  theme(text = element_text(size = 20),
      legend.position = "right",
      legend.margin = margin(0, 0, 0, -10),
      legend.title = element_text(size = 14, angle = -90, hjust = 0.5),
      legend.text = element_text(size = 10, margin = margin(0, 0, 0, -0.1)), 
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16)) +
  labs(size = bquote(italic("Module size")), 
    x = bquote("Pairwise cors of mod genes (Con)"),
    y = bquote("Pairwise cors of mod genes (AD)")) +    
  scale_size(range = c(0.1, 5), breaks = c(1, 25, 50, 75, 99)) +
  scale_color_manual(values = c("black" = "black", "red" = scales::hue_pal()(3)[1], "blue" = scales::hue_pal()(3)[3])) +
  scale_alpha_manual(values = c(0.05, 1, 1)) +
  guides(color = "none", alpha = "none", size = guide_legend(title.position = "right"))

ggsave(p3, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_cor.png"), bg = "white", width = 5, height = 5)
# There is a clear outlier here. What is it? 
meanmodcors1 <- meanmodcors1 |> mutate("diff" = Con - AD)
# meanmodcors1 |> arrange(desc(diff)) |> head()
# Module 1149 is the outlier
# > mod_bc[[1149]]
#  [1] "SELE"      "CXCL2"     "CSF3"      "CCL7"      "CCL8"      "IL6"      
#  [7] "VCAM1"     "CXCL3"     "CXCL1"     "CXCL8"     "C2CD4B"    "C20orf141"
# Many chemokines. interesting!
# This module has significantly higher correlations in control than AD samples. 
# What if we plotted this?
sif2_order <- sif2 |> arrange(alzdx)
expr2_order <- expr2[ ,match(rownames(sif2_order), colnames(expr2))] |>
  #apply(2, \(col) log2(col + 1)) |> 
  t() |> as.data.frame() |> dplyr::select(any_of(mod_bc[[1149]])) |>
  rownames_to_column(var = "sample") |>
  pivot_longer(!sample, names_to = "gene", values_to = "Expression") |>
  mutate("dx" = ifelse(sample %in% rownames(sif2_order)[sif2_order$alzdx == "AD"], "AD", "Con"))

p <- ggplot(expr2_order, aes(x = sample, y = Expression, group = gene, color = gene)) +
  theme_classic() + 
  geom_line() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~dx, ncol = 2, nrow = 1) +
  ggtitle("Module 1149 expression in ROSMAP Con vs AD")
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_cor_outlier.png"), width = 8, height = 3)
outliers <- expr2_order[expr2_order[,3] > 1000, ]
# There are basically 5 samples where expression of several genes (CSF3 especially, also SELE and CXCL1) spikes dramatically
sif2[rownames(sif2) %in% outliers[[1]], ]
sif[sif[ ,2] %in% outliers[[1]], ]

### Plot the results of the seattle analysis
thesead <- unique(copacohpos_complete[,1]) |> as.numeric()
# Which of the coherent mods align with copa mods unique for 1 subclass
eucposcommon <- unlist(euc_dist$pos[[3]])
eucposcommon <- unique(eucposcommon[duplicated(eucposcommon)])
eucposu <- lapply(euc_dist$pos[[3]], \(x){
  return(x[!x %in% eucposcommon])
})
sum(thesead %in% unlist(eucposu))
# 12 of the 20 coherent modules are significantly unique for a single celltype
thesead2 <- thesead[thesead %in% unlist(eucposu)]

copacohpos_complete$copa_unique <- F
for(i in 1:nrow(copacohpos_complete)){
  temp <- names(euc_dist$pos[[3]])[unlist(lapply(euc_dist$pos[[3]], \(x) copacohpos_complete$mod[i] %in% x))]
  if(length(temp) == 1){
    copacohpos_complete$copa_unique[i] <- temp
  } else if(length(temp) == 2){
    copacohpos_complete$copa_unique[i] <- paste(temp, collapse = " ")
  } else {
    copacohpos_complete$copa_unique[i] <- "multiple"
  }
}
copacohpos_complete |> arrange(as.numeric(mod))



# Plot differences in mean expr between ad and control for coherent modules
adconmeans <- lapply(thesead2, \(mod){
  mapply(\(sea, name){
    sea[rownames(sea) %in% mod_bc[[mod]], ] |>
      pivot_longer(everything(), names_to = "subclass", values_to = "mean") |>
      mutate("type" = name)
  }, sea_means[1:2], c("AD", "Con"), SIMPLIFY = F) |>
    do.call(what = "rbind") |>
    mutate(mod = mod)
}) |> do.call(what = "rbind")

p <- ggplot(adconmeans, aes(x = subclass, y = mean, fill = type)) + 
  theme_classic() + 
  geom_boxplot() +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  facet_wrap(~mod, ncol = 1)
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_sea.png"), height = 30, width = 10)

# Plot percentage differences
sub_diff_per <- ((sea_means[[1]] - sea_means[[2]]) / sea_means[[1]]) 
sub_diff_mod <- lapply(thesead2, \(mod){
  moddf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "pcnt_diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(pcnt_diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_mod)){
  if(sub_diff_mod$subclass[i] %in% copacohpos_complete$subclass[copacohpos_complete$mod == sub_diff_mod$modno[i]]){
    sub_diff_mod$highlight[i] <- T
  }
}  

p2 <- ggplot(sub_diff_mod, aes(x = subclass, y = meanpd, fill = highlight)) +
  theme_classic() + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  facet_wrap(~modno, ncol = 1)
ggsave(p2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_sea_pcntDiff.png"), height = 30, width = 10)

# Plot absolute difference
sub_diff_mod1 <- lapply(thesead2, \(mod){
  moddf <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_mod1)){
  if(sub_diff_mod1$subclass[i] %in% copacohpos_complete$subclass[copacohpos_complete$mod == sub_diff_mod1$modno[i]]){
    sub_diff_mod1$highlight[i] <- T
  }
}  

p3 <- ggplot(sub_diff_mod1, aes(x = subclass, y = meanpd, fill = highlight)) +
  theme_classic() + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  facet_wrap(~modno, ncol = 1)
ggsave(p3, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_sea_absDiff.png"), height = 30, width = 10)

# Plot neg modules
theseadneg <- unique(copacohneg_complete[,1]) |> as.numeric()
# Which of the coherent mods align with copa mods unique for 1 subclass
eucnegcommon <- unlist(euc_dist$neg[[3]])
eucnegcommon <- unique(eucnegcommon[duplicated(eucnegcommon)])
eucnegu <- lapply(euc_dist$neg[[3]], \(x){
  return(x[!x %in% eucnegcommon])
})
sum(theseadneg %in% unlist(eucnegu))
# 11 of the 20 coherent modules are significantly unique for a single celltype
theseadneg2 <- theseadneg[theseadneg %in% unlist(eucnegu)]

# Plot percentage differences (neg)
sub_diff_modneg <- lapply(theseadneg2, \(mod){
  moddf <- sub_diff_per[rownames(sub_diff_per) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "pcnt_diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(pcnt_diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_modneg)){
  if(sub_diff_modneg$subclass[i] %in% copacohneg_complete$subclass[copacohneg_complete$mod == sub_diff_modneg$modno[i]]){
    sub_diff_modneg$highlight[i] <- T
  }
}  

p2 <- ggplot(sub_diff_modneg, aes(x = subclass, y = meanpd, fill = highlight)) +
  theme_classic() + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  facet_wrap(~modno, ncol = 1)
ggsave(p2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_sea_pcntDiff_neg.png"), height = 30, width = 10)

# plot abs diff (neg)

sub_diff_modneg1 <- lapply(theseadneg2, \(mod){
  moddf <- sub_diff[rownames(sub_diff) %in% mod_bc[[mod]], ] |>
    pivot_longer(everything(), names_to = "subclass", values_to = "diff") |>
    mutate("modno" = mod) 
  return(moddf)
}) |> do.call(what = "rbind") |>
  group_by(subclass, modno) |>
  summarise("meanpd" = mean(diff)) |>
  mutate("highlight" = F)
for(i in 1:nrow(sub_diff_modneg1)){
  if(sub_diff_modneg1$subclass[i] %in% copacohneg_complete$subclass[copacohneg_complete$mod == sub_diff_modneg1$modno[i]]){
    sub_diff_modneg1$highlight[i] <- T
  }
}  

p4 <- ggplot(sub_diff_modneg1, aes(x = subclass, y = meanpd, fill = highlight)) +
  theme_classic() + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  facet_wrap(~modno, ncol = 1)
ggsave(p4, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_sea_absDiff_neg.png"), height = 30, width = 10)

sub_diff_modneg1[sub_diff_modneg1$highlight, ]

# Plot coherence vs copa significance? 

seacohall <- lapply(colnames(sub_diff), \(subclass_name){
  lapply(mod_bc, \(mod){
    temp <- sub_diff[rownames(sub_diff) %in% mod, colnames(sub_diff) == subclass_name]
  return(sum(temp > 0) / length(temp))
  }) |> unlist()
}) |> do.call(what = "cbind")

subdiffmodall <- rbind(sub_diff_mod1, 
                       sub_diff_modneg1)

# pval vs euclidean distance (highlighted)
plotdf <- lapply(which(colnames(seadist) %in% subdiffmodall[[1]]), \(col){
  sc <- colnames(seadist)[col]
  outdf <- cbind(seadist[, col], seapval[, col]) |> as.data.frame() |>
    mutate(highlight = "black", ct = sc)
  outdf$highlight[rownames(seadist) %in% subdiffmodall$modno[subdiffmodall$subclass == sc & subdiffmodall$meanpd > 0]] <- "red"
  outdf$highlight[rownames(seadist) %in% subdiffmodall$modno[subdiffmodall$subclass == sc & subdiffmodall$meanpd < 0]] <- "blue"
  return(outdf)
}) |> do.call(what = "rbind")
p <- ggplot(plotdf, aes(x = V1, y = -log10(V2), color = highlight, alpha = highlight)) +
  theme_classic() + 
  geom_point() + 
  labs(x = "Module euclidean distance (AD - Con)", 
       y = bquote(-log[10]~"p-value")) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "black" = "black"),
                     labels = c("red" = "higher in AD", "blue" = "higher in Con", "black" = "Other")) +
  scale_alpha_manual(values = c("red" = 2, "blue" = 2, "black" = 0.05),
                     labels = c("red" = "higher in AD", "blue" = "higher in Con", "black" = "Other")) +
  theme(text = element_text(size = 20), 
        legend.title = element_blank()) +
  ylim(NA, 3.5) + 
  facet_wrap(~ct, nrow = 3, ncol = 4, scales = "free_x")
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_sea_indivCTfacets.png"), height = 10, width = 10)

# pval vs euclidean distance (highlighted)
# For sake of graphing, set all pvalues that are zero to 0.005

plotdf <- lapply(which(colnames(seadist) %in% subdiffmodall[[1]]), \(col){
  sc <- colnames(seadist)[col]
  outdf <- cbind(seacohall[these_mods, col], seapval[, col]) |> as.data.frame() |>
    mutate(highlight = "black", ct = sc)
  outdf$highlight[rownames(seadist) %in% subdiffmodall$modno[subdiffmodall$subclass == sc & subdiffmodall$meanpd > 0]] <- "red"
  outdf$highlight[rownames(seadist) %in% subdiffmodall$modno[subdiffmodall$subclass == sc & subdiffmodall$meanpd < 0]] <- "blue"
  return(outdf)
}) |> do.call(what = "rbind") 
plotdf[plotdf[ ,2] == 0, 2] <- 0.005
p <- ggplot(plotdf, aes(x = V1, y = -log10(V2), color = highlight, alpha = highlight)) +
  theme_classic() + 
  geom_point() + 
  labs(x = "% of mod genes higher in AD", 
       y = bquote(-log[10]~"p-value (COPA)")) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "black" = "black"),
                     labels = c("red" = "All genes higher in AD", "blue" = "All genes higher in Con", "black" = "Other"),
                     breaks = c("red", "blue", "black")) +
  scale_alpha_manual(values = c("red" = 2, "blue" = 2, "black" = 0.05),
                     labels = c("red" = "All genes higher in AD", "blue" = "All genes higher in Con", "black" = "Other"),
                     breaks = c("red", "blue", "black")) +
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 12), 
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical") +
  ylim(NA, 3.5) + 
  facet_wrap(~ct, nrow = 3, ncol = 4) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.1)) +
  guides(alpha = "none")
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels/panel_sea_indivCTfacets_coh.png"), height = 6, width = 9)

### Plot modules
# source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/module_plotting_fxn.R"))
# plot_mods(ind = subdiffmodall$modno,
#           save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels"),
#           file_name = "mod_list.png")
subdiffallper <- rbind(sub_diff_mod, sub_diff_modneg)
#inds = unique(sub_diff_mod$modno)
inds = unique(subdiffallper$modno)
save_path = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_AD/panels")
file_name = "mod_list_all.png"

   save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/")
    module_output_dir <- save_dir1

  all_plots = list(qread(file.path(save_dir1,"sn_proj_objects","expr_object_bc.qs")), 
                qread(file.path(save_dir1,"sn_proj_objects","gsea_object.qs")))

  # for indexing and file names
  datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
  if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
  mods <- tapply(datkme[,2], datkme[,3], list)
  modulelengths <- unlist(lapply(mods,length))
  filter_under <- 3
  these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])

  cols <- RColorBrewer::brewer.pal(6, "Paired")
  cols2 <- RColorBrewer::brewer.pal(10, "Spectral")

  out_plot <- lapply(seq_along(inds), \(i){
    j <- which(these_mods == inds[i])
    newtheme <- theme_light() + theme(axis.text.x = element_text(size = 6), 
                                   axis.text.y = element_text(size = 6), 
                                   axis.title.x = element_text(size = 6), 
                                   axis.title.y = element_text(size = 6), 
                                   legend.text = element_text(size = 6),
                                   plot.title = element_blank(),
                                   plot.subtitle = element_blank(),
                                   legend.key.size = unit(0.2, "cm"), 
                                   legend.title = element_blank())
    plots <- lapply(all_plots, \(x) x[[j]] + newtheme)

    plots[[1]]$layers[[1]] <- NULL
    plots[[1]] <- plots[[1]] + 
        geom_line(linewidth = 0.2) +
        labs(x = "Sample") + 
        theme(axis.title.x = element_text(size = 6), 
              axis.title.y = element_text(size = 4),
              legend.text = element_text(size = 7), 
              legend.box.margin = margin(0, 0, 0, -10)) +
        scale_color_manual(values = cols2)
        

    plots[[2]] <- all_plots[[2]][[j]]$data |> 
        dplyr::filter(pval>0) |>
        dplyr::mutate(SetName=factor(SetName, levels=rev(unique(SetName)))) |>
        ggplot(aes(x = SetName, y = pval)) +
            theme_light() +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(size=6,hjust=1,vjust = 0.5),
                axis.title.x = element_text(size=6),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size=3),
                axis.text.y = element_text(size=4),
                plot.margin = margin(0,0,0,-20)) +
            labs(x="", y=bquote(-log[10](p-val))) #+
            #scale_x_discrete(limits = rev(levels(SetName))) +
            #geom_hline(yintercept = gsea_cutFDR, color = "red")
    plots[[2]]$layers[[2]] <- all_plots[[2]][[j]]$layers[[2]]
    plots[[2]] <- plots[[2]] + coord_flip()
    #p <- plot_grid(plotlist = plots, ncol = 4, nrow = 1, align = "h", axis = "bt", rel_widths=c(0.5, 0.85, 0.75, 1, 0.8))
    #return(p)
    plots[[3]] <- subdiffallper |> 
      dplyr::filter(modno == inds[i]) |>
      ggplot(aes(x = subclass, y = meanpd, fill = highlight)) +
        theme_classic() + 
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=6),
              axis.text.y = element_text(size=6),
              #legend.text = element_text(size = 6),
              #legend.key.size = unit(0.2, "cm"),
              legend.position = "none")  +
        labs(y = "% change in expr")

    return(plots)
  })

    out_plot2 <- lapply(c(1,2,3), \(x){
        plist <- lapply(out_plot, \(y) y[[x]])
        if(x != 2){
          return(plot_grid(plotlist = plist, nrow = length(plist), align = "v", axis = "rl"))
        } else {
          return(plot_grid(plotlist = plist, nrow = length(plist)))
        }
    })

    outall2 <- plot_grid(plotlist = out_plot2, ncol = 3, align = "h", axis = "bt", rel_widths = c(0.6,0.4,0.8))
    ggsave(outall2, file = file.path(save_path, file_name),  height = length(inds) * 1, width = 6, bg = "white", limitsize = F)


# Let's add early/late ad modules to this calculus
euc_dist_el <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/euclidean_distances/euclidean_sigmods_all.qs")),
                    "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/euclidean_distances/euclidean_sigmods_positive.qs")),
                    "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/euclidean_distances/euclidean_sigmods_negative.qs")))
sea_means_el <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_earlyVsLate/sn_summary_tables/allmlist_log.qs"))
sub_diff_el <- sea_means_el[[1]] - sea_means_el[[2]]

ccpos <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff_el[rownames(sub_diff_el) %in% mod_bc[[mod]], colnames(sub_diff_el) == subclass_name]
    return(sum(temp > 0) / length(temp))
  }) |> unlist() 
  if(length(outvec) > 0){
    return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
  } else {
    return(data.frame("subclass" = subclass_name, "pcnt" = NA))
  }
}, euc_dist_el$pos[[3]], names(euc_dist_el$pos[[3]]), SIMPLIFY = F)
# How many of these overlap with the conVsAD modules?
t1 <- lapply(ccpos, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # earlyvslate
t2 <- lapply(copacohpos, \(x) rownames(x[x$pcnt == 1, ])) |> unlist() # conVsAD
unique(t1)[which(unique(t1) %in% unique(t2))]
# 149 716
unique(t2)[which(unique(t2) %in% unique(t1))]
# 149 716

# Which are specific to a subclass?
ccpos_com <- lapply(ccpos, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
#tapply(ccpos_com[ ,2], ccpos_com[ ,1], list)
ccpos_u <- tapply(ccpos_com[ ,2], ccpos_com[ ,1], list)
ccpos_u <- ccpos_u[(lapply(ccpos_u, length) |> unlist()) == 1] |> unlist()
    #      25         295         301         441         716         744 
    # "L6 CT"     "L6 CT"     "L6 CT"      "Sncg"     "L6 CT"       "Vip" 
    #      87         909 
    #  "VLMC" "Astrocyte" 

# Negative (higher in con)
ccneg <- mapply(\(subclass, subclass_name){
  outvec <- lapply(subclass, \(mod){
    temp <- sub_diff_el[rownames(sub_diff_el) %in% mod_bc[[mod]], colnames(sub_diff_el) == subclass_name]
    return(sum(temp < 0) / length(temp))
  }) |> unlist() 
  return(data.frame("subclass" = subclass_name, "pcnt" = outvec))
}, euc_dist_el$neg[[3]], names(euc_dist_el$neg[[3]]), SIMPLIFY = F) 
lapply(ccneg, \(x) sum(x$pcnt == 1)) |> unlist()
#       Astrocyte      Chandelier     Endothelial         L2/3 IT           L4 IT 
#               2               7               1               1              37 
#           L5 ET           L5 IT         L5/6 NP           L6 CT           L6 IT 
#               0               5               1               2               3 
#      L6 IT Car3             L6b           Lamp5      Lamp5 Lhx6   Microglia-PVM 
#               7               2              18               2               0 
# Oligodendrocyte             OPC            Pax6           Pvalb            Sncg 
#              15               3               1               1               1 
#             Sst       Sst Chodl             Vip            VLMC 
#               0               0               1               0 
# How many of these modules are specific to a single subclass?
ccneg_com <- lapply(ccneg, \(x){
  x[x$pcnt == 1, ] |>
    rownames_to_column(var = "mod")
}) |> do.call(what = "rbind")
ccneg_u <- tapply(ccneg_com[ ,2], ccneg_com[ ,1], list)
ccneg_u <- ccneg_u[(lapply(ccneg_u, length) |> unlist()) == 1] |> unlist()
#                10               101              1035              1049                                                                                         
#           "L4 IT"           "L4 IT"           "L6 CT"       "Astrocyte" 
#              1086                11                13               131 
#           "L4 IT"           "L4 IT"           "L4 IT"           "L4 IT" 
#                14               140               141               147 
#           "L4 IT"           "Lamp5"           "L4 IT"       "Astrocyte" 
#               157                17               173               184 
#           "L4 IT"           "Lamp5" "Oligodendrocyte"           "Lamp5" 
#               192                21               227                25 
# "Oligodendrocyte"           "L4 IT"           "L4 IT" "Oligodendrocyte" 
#               260               272                28               316 
# "Oligodendrocyte"           "L4 IT"           "L4 IT"           "L4 IT" 
#               329               334                35               442 
#      "L6 IT Car3"           "L4 IT"           "L4 IT"         "L2/3 IT" 
#               452               453                48               489 
#           "Lamp5"            "Sncg"           "L4 IT"           "Lamp5" 
#               493               517               518               526 
#      "L6 IT Car3"           "Lamp5"      "L6 IT Car3" "Oligodendrocyte" 
#                58               614               623               639 
#           "Lamp5"           "L6 IT" "Oligodendrocyte"      "Chandelier" 
#                66               678               716               721 
#           "L4 IT"           "L4 IT" "Oligodendrocyte"      "L6 IT Car3" 
#               754                82                88               893 
# "Oligodendrocyte"           "L4 IT"      "Chandelier"     "Endothelial" 
#                90               910               968               984 
#           "L4 IT"           "L6 CT"           "L4 IT"         "L5/6 NP" 
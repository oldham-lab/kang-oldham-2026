library(tidyverse)
library(data.table)
library(qs)
library(doParallel)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
registerDoParallel(cores=6)

## Load datasets of interest
# Load bulk megaset, lein pb, and lein cellbender-scVI-imputed pb
mega <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
leinpb <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_06-13-16_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS"))[[1]]
scvi <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI_cellbender/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_04-15-41_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS"))[[1]]

# Align genes 
common_genes <- intersect(intersect(mega[,2], leinpb[,2]), scvi[,2])
mega <- mega[mega[,2] %in% common_genes,]
mega <- mega[match(common_genes, mega[,2]),]
leinpb <- leinpb[leinpb[,2] %in% common_genes,]
leinpb <- leinpb[match(common_genes, leinpb[,2]),]
scvi <- scvi[scvi[,2] %in% common_genes,]
scvi <- scvi[match(common_genes, scvi[,2]),]

# Transpose expression matrices
dat_names <- c("mega", "pb", "scvi")
graph_names <- c("Bulk", "PB", "PB_CB_scVI")
graph_names2 <- list(mega="Bulk", pb="PB", scvi="PB_CB_scVI")
matlist <- list()
matlist[[1]] <- t(mega[,-c(1:2)])
matlist[[2]] <- t(leinpb[,-c(1:2)])
matlist[[3]] <- t(scvi[,-c(1:2)])
names(matlist) <- dat_names
matlist <- lapply(matlist, function(x) apply(x,2,scale))

# Load adjacency matrices and convert to correlation matrices
simlist <- list()
simlist[[1]] <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/megamat_adjMat.qs"))
simlist[[2]] <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/leinpb_donor1_adjMat.qs"))
simlist[[3]] <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/scvi_donor1_adjMat.qs"))
simlist <- lapply(simlist, function(x) 2*x-1)
names(simlist) <- dat_names

# Load means for original SN data
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means_byDonor.qs"))
san1 <- san_mean[[1]][names(san_mean[[1]]) %in% common_genes]
san1 <- san1[match(scvi[,2], names(san1))]
# Calculate pcntiles for san1
calc_pcntile <- function(vec){
    (1-(rank(-vec,na.last="keep")/max(rank(-vec,na.last="keep"))))*100
}
sandf <- data.frame("Gene"=names(san1),
                    "Pcntile"=calc_pcntile(san1),
                    "Bracket"=NA)
sandf$Bracket[sandf$Pcntile<20] <- "0-20"
sandf$Bracket[sandf$Pcntile>=20 & sandf$Pcntile < 40] <- "20-40"
sandf$Bracket[sandf$Pcntile>=40 & sandf$Pcntile < 60] <- "40-60"
sandf$Bracket[sandf$Pcntile>=60 & sandf$Pcntile < 80] <- "60-80"
sandf$Bracket[sandf$Pcntile>=80 & sandf$Pcntile < 100] <- "80-99"

# Order knn genes for each gene in each dataset
# knnlist <- lapply(simlist, function(a) apply(a,2,function(x) common_genes[order(x,decreasing=T)]))

# For each dataset, model each gene using knn genes

#k_vec <- 1:20
# r2a_list <- list()
# rmse_list <- list()
# for(k in seq_along(k_vec)){
#     temp1 <- list()
#     temp2 <- list()
#     for(i in 1:3){
#         bla1 <- list()
#         bla2 <- list()
#         for(j in c(1:3)){
#             mod_list1 <- foreach(l=1:ncol(matlist[[i]])) %dopar% {
#                 lm(matlist[[j]][,l] ~ ., data=as.data.frame(matlist[[j]][,common_genes %in% knnlist[[i]][1:k_vec[k],l]]))
#             }
#             bla1[[j]] <- data.frame("modeling"=dat_names[j],"model_with"=dat_names[i],"r2"=unlist(lapply(mod_list1, function(x) summary(x)$adj.r.squared)))
#             bla2[[j]] <- data.frame("modeling"=dat_names[j],"model_with"=dat_names[i],"rmse"=unlist(lapply(mod_list1, function(x) sqrt(mean(x$residuals^2)))))
#             rm(mod_list1)
#             gc()
#         }
#         temp1[[i]] <- bla1 %>% do.call(rbind,.)
#         temp2[[i]] <- bla2 %>% do.call(rbind,.)
#     }
#     r2a_list[[k]] <- temp1 %>% do.call(rbind,.)
#     rmse_list[[k]] <- temp2 %>% do.call(rbind,.)
#     cat(k, " ")
# }
# names(r2a_list) <- paste0("k=", 1:20)
# names(rmse_list) <- paste0("k=", 1:20)
# qsave(r2a_list,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/r2a_list.qs"))
# qsave(rmse_list,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/rmse_list.qs"))
#r2a_list <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/r2a_list.qs"))
#rmse_list <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/rmse_list.qs"))

# Plot violin plots over all values of k (non-reciprocal only, i.e. bulk predicting bulk, pb predicting pb)
#brewer_colors <- brewer.pal(6,"Paired")[c(2,4,6)]
#comps <- list(c("Bulk", "PB"), c("Bulk", "PB_CB_scVI"), c("PB", "PB_CB_scVI"))

# p <- r2a_list %>%
#     mapply(function(dat,name){
#         out <- dat %>% mutate("k"=gsub("k=", "", name))
#     },.,names(.), SIMPLIFY=F) %>% do.call(rbind,.) %>%
#     dplyr::filter(modeling==model_with) %>% 
#     mutate(k=as.numeric(k)) %>%
#     mutate(modeling=factor(modeling, labels=graph_names2)) %>%
#     ggplot(aes(x=k, y=r2, group=k)) + 
#         theme_bw() +
#         geom_violin(aes(fill=k)) +
#         geom_boxplot(fill="white", width=0.1, notch=T, outlier.shape=NA) + 
#         theme(text=element_text(size=30),
#             plot.margin=margin(1,1,1,1,"cm"),
#             legend.position="none",
#             axis.text.x=element_text(size=20)) +
#         labs(x="# of nearest neighbors", y=bquote("Adjusted"~R^2)) +
#         #scale_fill_brewer() +
#         #stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=4)  +
#         facet_wrap(~modeling, ncol=3)
# ggsave(p,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn_"), k_vec[1],"To",k_vec[length(k_vec)],"_R2.png"), width=14, height=5)

# p <- rmse_list %>%
#     mapply(function(dat,name){
#         out <- dat %>% mutate("k"=gsub("k=", "", name))
#     },.,names(.), SIMPLIFY=F) %>% do.call(rbind,.) %>%
#     dplyr::filter(modeling==model_with) %>% 
#     mutate(k=as.numeric(k)) %>%
#     mutate(modeling=factor(modeling, labels=graph_names2)) %>%
#     ggplot(aes(x=k, y=rmse, group=k)) + 
#         theme_bw() +
#         geom_violin(aes(fill=k)) +
#         geom_boxplot(fill="white", width=0.1, notch=T, outlier.shape=NA) + 
#         theme(text=element_text(size=30),
#             plot.margin=margin(1,1,1,1,"cm"),
#             legend.position="none",
#             axis.text.x=element_text(size=20)) +
#         labs(x="# of nearest neighbors", y="RMSE") +
#         #scale_fill_brewer() +
#         #stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=4)  +
#         facet_wrap(~modeling, ncol=3)
# ggsave(p,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn_"), k_vec[1],"To",k_vec[length(k_vec)],"_RMSE.png"), width=14, height=5)

# Plot modeling values vs san means

# p <- r2a_list %>%
#     mapply(function(dat,name){
#         out <- dat %>% mutate("k"=name)
#     },.,names(.), SIMPLIFY=F) %>% do.call(rbind,.) %>%
#     dplyr::filter(modeling==model_with) %>% 
#     mutate(modeling=factor(modeling, labels=graph_names2)) %>%
#     mutate(sanmean = rep(san1, 60)) %>%
#     dplyr::filter(k %in% c("k=5","k=10","k=15","k=20")) %>%
#     mutate(k=factor(k, levels=unique(k))) %>%
#     ggplot(aes(x=r2, y=sanmean, color=k)) + 
#         theme_bw() +
#         geom_point(size=0.1, alpha=0.2) +
#         #geom_smooth() +
#         theme(text=element_text(size=30),
#             plot.margin=margin(1,1,1,1,"cm"),
#             legend.position="bottom",
#             axis.text.x=element_text(size=20)) +
#         labs(y="Mean UMI counts (Lein DFC)", x=bquote("Adjusted"~R^2)) +
#         scale_color_brewer() +
#         #stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=4)  +
#         facet_wrap(vars(modeling), ncol=3, nrow=4) +
#         guides(color=guide_legend(override.aes=list(size=2, alpha=1)))
# ggsave(p,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn_"), k_vec[1],"To",k_vec[length(k_vec)],"_vs_mean_R2.png"), width=16, height=6)

# Plot modeling values vs san means (filtered for k=10)

# p <- r2a_list %>%
#     mapply(function(dat,name){
#         out <- dat %>% mutate("k"=name)
#     },.,names(.), SIMPLIFY=F) %>% do.call(rbind,.) %>%
#     dplyr::filter(modeling==model_with) %>% 
#     mutate(modeling=factor(modeling, labels=graph_names2)) %>%
#     dplyr::filter(k =="k=10") %>%
#     mutate(sanmean = rep(san1, 3)) %>%
#     ggplot(aes(x=r2, y=sanmean)) + 
#         theme_bw() +
#         geom_line(linewidth=0.05) +
#         theme(text=element_text(size=30),
#             plot.margin=margin(1,1,1,1,"cm"),
#             legend.position="bottom",
#             axis.text.x=element_text(size=20),
#             axis.title.y=element_text(margin=margin(0,10,0,0))) +
#         labs(y="Mean UMI counts (Lein DFC)", x=bquote("Adjusted"~R^2)) +
#         #scale_color_brewer() +
#         #stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=4)  +
#         facet_wrap(vars(modeling), nrow=3) #+
#         #guides(color=guide_legend(override.aes=list(linewidth=2)))
# ggsave(p,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn10_vs_mean_R2.png")), width=8, height=10)


# Plot violin plots over all values of k and all combinations of datasets (reciprocal included)

# p <- r2a_list %>%
#     mapply(function(dat,name){
#         out <- dat %>% mutate("k"=gsub("k=", "", name))
#     },.,names(.), SIMPLIFY=F) %>% do.call(rbind,.) %>%
#     mutate(modeling=factor(modeling, labels=graph_names2)) %>%
#     mutate(model_with=factor(model_with, labels=graph_names2)) %>%
#     mutate(model_lab=paste0(modeling, " modeled by ", model_with) %>% as.factor)  %>% 
#     mutate(k=as.numeric(k)) %>%
#     ggplot(aes(x=k, y=r2, group=k)) + 
#         theme_bw() +
#         geom_violin(aes(fill=k)) +
#         geom_boxplot(fill="white", width=0.1, notch=T, outlier.shape=NA) + 
#         theme(text=element_text(size=30),
#             plot.margin=margin(1,1,1,1,"cm"),
#             legend.position="none",
#             axis.text.x=element_text(size=20)) +
#         labs(x="# of nearest neighbors", y=bquote("Adjusted"~R^2)) +
#         facet_wrap(~model_lab, ncol=3, nrow=3)
# ggsave(p,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn_"), k_vec[1],"To",k_vec[length(k_vec)],"_R2_all_combinations.png"), width=20, height=8)


# Repeat modeling but use average of kme genes rather than doing multiple linear regression

# k_vec <- 1:20
# r2a_1_list <- list()
# rmse_1_list <- list()
# for(k in seq_along(k_vec)){
#     temp1 <- list()
#     temp2 <- list()
#     for(i in 1:3){
#         bla1 <- list()
#         bla2 <- list()
#         for(j in c(1:3)){
#             mod_list1 <- foreach(l=1:ncol(matlist[[i]])) %dopar% {
#                 lm(matlist[[j]][,l] ~ rowMeans(as.data.frame(matlist[[j]][,common_genes %in% knnlist[[i]][1:k_vec[k],l]])))
#             }
#             bla1[[j]] <- data.frame("modeling"=dat_names[j],"model_with"=dat_names[i],"r2"=unlist(lapply(mod_list1, function(x) summary(x)$adj.r.squared)))
#             bla2[[j]] <- data.frame("modeling"=dat_names[j],"model_with"=dat_names[i],"rmse"=unlist(lapply(mod_list1, function(x) sqrt(mean(x$residuals^2)))))
#             rm(mod_list1)
#             gc()
#         }
#         temp1[[i]] <- bla1 %>% do.call(rbind,.)
#         temp2[[i]] <- bla2 %>% do.call(rbind,.)
#     }
#     r2a_1_list[[k]] <- temp1 %>% do.call(rbind,.)
#     rmse_1_list[[k]] <- temp2 %>% do.call(rbind,.)
#     cat(k, " ")
# }
# names(r2a_1_list) <- paste0("k=", 1:20)
# names(rmse_1_list) <- paste0("k=", 1:20)
# qsave(r2a_1_list,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/r2a_list_1predictor.qs"))
# qsave(rmse_1_list,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/rmse_list_1predictor.qs"))
r2a_1_list <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/r2a_list_1predictor.qs"))
rmse_1_list <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/rmse_list_1predictor.qs"))

# Plot violin plots over all values of k and all combinations of datasets (reciprocal included)

p <- r2a_1_list %>%
    mapply(function(dat,name){
        out <- dat %>% mutate("k"=gsub("k=", "", name))
    },.,names(.), SIMPLIFY=F) %>% do.call(rbind,.) %>%
    mutate(modeling=factor(modeling, labels=graph_names2)) %>%
    mutate(model_with=factor(model_with, labels=graph_names2)) %>%
    mutate(model_lab=paste0(modeling, " modeled by ", model_with) %>% as.factor)  %>% 
    mutate(k=as.numeric(k)) %>%
    ggplot(aes(x=k, y=r2, group=k)) + 
        theme_bw() +
        geom_violin(aes(fill=k)) +
        geom_boxplot(fill="white", width=0.1, notch=T, outlier.shape=NA) + 
        theme(text=element_text(size=30),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="none",
            axis.text.x=element_text(size=20)) +
        labs(x="# of nearest neighbors", y=bquote("Adjusted"~R^2)) +
        facet_wrap(~model_lab, ncol=3, nrow=3)
ggsave(p,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn_"), k_vec[1],"To",k_vec[length(k_vec)],"_R2_all_combinations_1_predictor.png"), width=20, height=8)

# Plot modeling values vs san means (filtered for k=10)

p <- r2a_1_list %>%
    mapply(function(dat,name){
        out <- dat %>% mutate("k"=name)
    },.,names(.), SIMPLIFY=F) %>% do.call(rbind,.) %>%
    dplyr::filter(k =="k=10") %>%
    mutate(modeling=factor(modeling, labels=graph_names2)) %>%
    mutate(model_with=factor(model_with, labels=graph_names2)) %>%
    mutate(model_lab=paste0(modeling, " modeled by ", model_with) %>% as.factor) %>% 
    mutate(bracket = rep(sandf$Bracket, 9)) %>%
    ggplot(aes(y=r2, x=bracket)) + 
        theme_bw() +
        geom_violin(aes(fill=bracket)) +
        geom_boxplot(fill="white", width=0.2, notch=T, outlier.shape=NA) + 
        theme(text=element_text(size=30),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="none",
            axis.text.x=element_text(size=20),
            axis.title.y=element_text(margin=margin(0,10,0,0))) +
        labs(x="Mean UMI percentile", y=bquote("Adjusted"~R^2)) +
        scale_fill_brewer() +
        #stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=4)  +
        facet_wrap(~model_lab,nrow=3, ncol=3) #+
        #guides(color=guide_legend(override.aes=list(linewidth=2)))
ggsave(p,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn10_vs_meanBracket_R2_1predictor.png")), width=20, height=10)

## Which markers stand out as outliers?
# Filter to k=10
k10 <- r2a_1_list$`k=10` %>%
  mutate(Gene = rep(common_genes,9),
         category = paste0(k10[,1], "_", k10[,2])) %>%
  filter(category %in% c("mega_mega", "mega_scvi", "scvi_mega", "scvi_scvi"))

#topgenes <- k10 %>% group_by(category) %>% arrange(desc(r2), .by_group=T) %>% group_map(~head(.x,10L))

k10w <- pivot_wider(k10[,3:5], names_from=category, values_from=r2)
fwrite(k10w,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn_modeling_megaset_vs_scvi.csv")))

combns <- combn(2:5,2)
combnnames <- combn(colnames(k10w)[-1], 2)
plist <- list()
for(i in 1:ncol(combns)){
    tempdf <- k10w[,c(1, combns[1:2,i])]
    colnames(tempdf)[2:3] <- c("xdata", "ydata")
    plist[[i]] <- ggplot(tempdf, aes(x=xdata, y=ydata)) + 
      theme_bw() + 
      geom_point(alpha=0.5,shape=16, fill="black", size=0.1) + 
      #stat_density_2d(geom = "polygon", contour = TRUE,
      #                aes(fill = after_stat(level)), colour = "black",bins = 5) +
      #scale_fill_distiller(palette = "Blues", direction = 1) +
      theme(text=element_text(size=30),
            plot.margin=margin(1,1,1,1,"cm"),
            legend.position="none") +
      labs(x=combnnames[1,i],y=combnnames[2,i]) +
      xlim(0,1) +
      ylim(0,1)
}
p <- plot_grid(plotlist=plist, ncol=3,nrow=2)
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/knn_modeling_megaset_vs_scvi_scatter.png"), width=20,height=10)

# Genes well modeled by scvi neighbors in scvi but poorly in bulk
ktest <- k10w[,c(1,4,5)] %>%
  mutate(diff=scvi_scvi-mega_scvi) %>%
  arrange(desc(diff)) %>% 
  left_join(sandf, by=join_by(Gene)) %>%
  as.data.frame
head(ktest,40)
# Mostly low percentile expression genes

# Genes well modeled by bulk neighbors in bulk but poorly in scvi
ktest2 <- k10w[,c(1,2,3)] %>%
  mutate(diff=mega_mega-scvi_mega) %>%
  arrange(desc(diff)) %>% 
  left_join(sandf, by=join_by(Gene)) %>%
  as.data.frame
head(ktest2,40)
# many highly expressed genes
# NRXN1
#  Gene  mega_mega scvi_mega mega_scvi scvi_scvi
#   <chr>     <dbl>     <dbl>     <dbl>     <dbl>
# 1 NRXN1     0.940   0.00237     0.410     0.659

kallp <- lapply(r2a_1_list,function(x){
  x <- x %>% mutate(Gene = rep(common_genes,9),
             category = paste0(modeling, "_", model_with))
})
k10w <- kallp$`k=10` 
k10w <- pivot_wider(k10w[,3:5], names_from=category, values_from=r2)

# Which genes are downgraded or upgraded by imputation
sum(k10w$scvi_mega < k10w$pb_mega & k10w$scvi_mega < k10w$mega_mega)
# 3164 - downgraded by imputation
sum(k10w$scvi_mega > k10w$pb_mega & k10w$scvi_mega < k10w$mega_mega)
# 11037 - upgraded by imputation but still worse than self modeling
sum(k10w$scvi_mega > k10w$pb_mega & k10w$scvi_mega > k10w$mega_mega)
# 2984 - upgraded by imputation compared to both datasets

sum(k10w$scvi_mega == k10w$pb_mega)
# 0
sum(k10w$scvi_mega < k10w$pb_mega)
# 3172 mega neighbors perform better pre-correction
sum(k10w$scvi_mega > k10w$pb_mega)
# 14021 mega neighbors perform better post-correction

sum(k10w$pb_pb < k10w$scvi_pb)
# 10753 pb neighbors perform worse after imputation
sum(k10w$pb_pb > k10w$scvi_pb)
# 6400 pb neighbors perform better after imputation
sum(k10w$pb_pb == k10w$scvi_pb)
# 0

sum(k10w$pb_pb < k10w$scvi_scvi)
# 17014
sum(k10w$pb_pb == k10w$scvi_scvi)
# 0
sum(k10w$pb_pb > k10w$scvi_scvi)
# 179
# Vast majority of genes are better explained by neighbors after imputation

sum(k10w$mega_mega < k10w$scvi_scvi)
# 12362
sum(k10w$mega_mega > k10w$scvi_scvi)
# 4831

sum(k10w$mega_mega < k10w$mega_scvi)
# 11
sum(k10w$mega_mega > k10w$mega_scvi)
# 17178

sum(k10w$scvi_mega < k10w$scvi_scvi)
# 17176
sum(k10w$scvi_mega > k10w$scvi_scvi)
# 13

sum(k10w$mega_scvi < k10w$scvi_scvi)
# 16642
sum(k10w$mega_scvi > k10w$scvi_scvi)
# 551
sum(k10w$pb_scvi < k10w$scvi_scvi)
# 17104
sum(k10w$pb_scvi > k10w$scvi_scvi)
# 89

# The corrected neighbors are internally consistent
# but what about the difference between mega neighbors predicting pb vs corrected
ktest <- data.frame("Gene"=k10w$Gene, "pb_vs_corrected"=k10w$pb_mega-k10w$scvi_mega) %>% arrange(desc(pb_vs_corrected)) %>%
  left_join(sandf,by=join_by(Gene))

p <- ggplot(ktest, aes(x=pb_vs_corrected, y=Pcntile)) +
  theme_bw() + 
  geom_point(alpha=0.5, size=0.2) +
  theme(text=element_text(size=30)) +
  geom_smooth()
ggsave(p,file="~/test/test.png")


# Load zero percentages
zero_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/pctage_of_zeros_per_gene.qs"))
zero1 <- zero_mean[names(zero_mean) %in% common_genes]
zero1 <- zero1[match(common_genes, names(zero1))]
# Calculate pcntiles for san1
calc_pcntile <- function(vec){
    (1-(rank(-vec,na.last="keep")/max(rank(-vec,na.last="keep"))))*100
}
zerodf <- data.frame("Gene"=names(zero1),
                    "Pcntile"=calc_pcntile(zero1),
                    "Bracket"=NA)
zerodf$Bracket[zerodf$Pcntile<20] <- "0-20"
zerodf$Bracket[zerodf$Pcntile>=20 & zerodf$Pcntile < 40] <- "20-40"
zerodf$Bracket[zerodf$Pcntile>=40 & zerodf$Pcntile < 60] <- "40-60"
zerodf$Bracket[zerodf$Pcntile>=60 & zerodf$Pcntile < 80] <- "60-80"
zerodf$Bracket[zerodf$Pcntile>=80 & zerodf$Pcntile < 100] <- "80-99"

# Plot by zero percentage
ktest <- data.frame("Gene"=k10w$Gene, "pb_vs_corrected"=k10w$pb_mega-k10w$scvi_mega) %>% arrange(desc(pb_vs_corrected)) %>%
  left_join(zerodf,by=join_by(Gene))

p <- ggplot(ktest, aes(x=pb_vs_corrected, y=Pcntile)) +
  theme_bw() + 
  geom_point(alpha=0.5, size=0.2) +
  theme(text=element_text(size=30)) +
  geom_smooth()
ggsave(p,file="~/test/test.png")
# pretty much the same pattern as with mean expression

# Calculate mean expr percentiles from lein pb
leinpbmean <- rowMeans(leinpb[,-c(1,2)])
# Calculate pcntiles for san1
leinpbmeandf <- data.frame("Gene"=common_genes,
                    "Pcntile"=calc_pcntile(leinpbmean),
                    "Bracket"=NA)
leinpbmeandf$Bracket[leinpbmeandf$Pcntile<20] <- "0-20"
leinpbmeandf$Bracket[leinpbmeandf$Pcntile>=20 & leinpbmeandf$Pcntile < 40] <- "20-40"
leinpbmeandf$Bracket[leinpbmeandf$Pcntile>=40 & leinpbmeandf$Pcntile < 60] <- "40-60"
leinpbmeandf$Bracket[leinpbmeandf$Pcntile>=60 & leinpbmeandf$Pcntile < 80] <- "60-80"
leinpbmeandf$Bracket[leinpbmeandf$Pcntile>=80 & leinpbmeandf$Pcntile < 100] <- "80-99"

ktest <- data.frame("Gene"=k10w$Gene, "pb_vs_corrected"=k10w$pb_mega-k10w$scvi_mega) %>% arrange(desc(pb_vs_corrected)) %>%
  left_join(leinpbmeandf, by=join_by(Gene))

p <- ggplot(ktest, aes(x=pb_vs_corrected, y=Pcntile)) +
  theme_bw() + 
  geom_point(alpha=0.5, size=0.2) +
  theme(text=element_text(size=30)) +
  geom_smooth()
ggsave(p,file="~/test/test.png")

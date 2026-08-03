# A code document for follow-up analyses after producing figure 1 using CompareMarkers

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
library(cowplot)
library(tidyverse)
library(qs)
library(ggrepel)
library(ComplexHeatmap)
options(bitmapType = 'cairo')

library(doParallel)
registerDoParallel(cores=6)

####################
# Load required data
####################
save_dirall <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/")
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
origdemat2 <- fread(data.table=F, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein_orig_de_genes.csv"))
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)

###########
# Testing
###########

# Load model results
modelres <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_topmodposFDR_Subclass.qs"))

# Load t-values for all 3 donors
# ctstring <- "Subclass"
# ze=1
# modtlist <- list()
# modrlist <- list()
# for(i in 1:3){
#   exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
#   exdirexpr <- exdir[grep("EXPRLIST",exdir)]
#   exdirsif <- exdir[grep("samples_legend",exdir)]
#   exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
#   exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))

#   x <- exdirexpr[5]
#   y <- exdirsif[5]
#   z <- exdirnames[5]
#   obj <- readRDS(x) 
#   obj <- norm_pb_samples(obj[[1]], index=3)
#   exprt <- t(obj[,3:ncol(obj)])
#   escale <- apply(exprt, 2, scale)
#   rownames(escale) <- rownames(exprt)
  
#   sif <- fread(data.table=F,file=y)
#   sampct <- t(calc_CT_counts(sif, cell_anno, ct_col=ze, cell_id_col = 3))
#   sampctmean <- apply(sampct, 2, mean)
#   sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
  
#   # Run modeling
    
#   mod_list <- foreach(l=1:ncol(escale)) %dopar% {
#     lm(escale[,l] ~ ., data=as.data.frame(sampct))
#   }
#   modtlist[[i]] <- foreach(l=seq_along(mod_list)) %dopar% {
#     summary(mod_list[[l]])$coefficients[,3:4]
#   }
#   names(modtlist[[i]]) <- obj$Gene
  
#   modrlist[[i]] <- foreach(l=seq_along(mod_list)) %dopar% {
#     summary(mod_list[[l]])$adj.r.squared
#   }
#   names(modrlist[[i]]) <- obj$Gene

#   rm(mod_list)
#   gc()
# }
# qsave(modtlist, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/subclass_tvalues.qs"))
# qsave(modrlist, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/subclass_r2values.qs"))

modtlist <- qread(file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/subclass_tvalues.qs"))
modrlist <- qread(file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/subclass_r2values.qs"))

#####
# Calculate specificity and sensitivity for each gene over all subclasses
####

donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")

# collect lists of cell ids per subclass for each donor
sublist_all <- lapply(donorvec, function(donor){
  cell_anno1 <- cell_anno %>% dplyr::filter(Donor==donor)
  return(tapply(cell_anno1[,3], cell_anno1[,1], list)) 
})

# fidlist <- mapply(function(sublist, donor){
#   # Calculate sensitivity
#   sensi <- lapply(sublist, function(sub){
#     cell_exprall[,colnames(cell_exprall) %in% sub] %>%
#       apply(.,1, function(r) sum(r>0)/length(r))
#   }) %>% do.call(cbind,.) %>% t

#   # Calculate specificity
#   speci <- lapply(sublist, function(x) rowSums(cell_exprall[,colnames(cell_exprall) %in% x])) %>%
#     do.call(cbind,.) %>% apply(., 1, function(x) x/sum(x)) 
  
#   # Calculate fidelity
#   fid <- mapply(function(se, sp) se*sp, as.list(as.data.frame(sensi)), as.list(as.data.frame(speci)), SIMPLIFY=F) %>% do.call(cbind,.)

#   return(list("sensitivity"=sensi, "specificity"=speci, "fidelity"=fid))
# }, sublist_all, donorvec, SIMPLIFY=F)
# names(fidlist) <- c("Donor1", "Donor2", "Donor3")
# fidlist <- lapply(fidlist, function(x){
#   lapply(x, function(y){
#     out <- y %>% as.data.frame
#     rownames(out) <- names(sublist_all[[1]])
#     return(out)
#   })
# })
# qsave(fidlist,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/specificity_sensitivity_fidelity.qs"))

fidlist <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/specificity_sensitivity_fidelity.qs"))
# list of 3: donor1-3
# - list of 3: sensitivity, specificity, fidelity

lapply(fidlist[[1]], function(x){
  x %>% dplyr::select("AIF1") %>% arrange(desc(AIF1))
})
# interesting how specificity is decent across the board but sensitivity of badly modeled markers like AIF1 is low

# load kevin's aomn hi-fi genes
fidaomn <- fread(data.table=F,file="~/code/git/CompareMarkers/data/bulkFidelity_Hs.FCX_AOMN.csv")

# What is the r2 vs sensitivity for DE genes?

leindegenes <- tapply(origdemat2$genes, origdemat2$celltype, list)

modr2proc <- lapply(modrlist, function(x){
  out <- t(as.data.frame(x)) %>% as.data.frame %>% rownames_to_column(var="Gene")
  colnames(out)[2] <- "r2"
  return(out)
})

modtproc <- lapply(modtlist[[1]], function(x){
  x[,1]
}) %>% do.call(cbind,.) %>% t %>% as.data.frame
colnames(modtproc)[1] <- "Sst Chodl"
colnames(modtproc) <- gsub("`", "", colnames(modtproc))
modtproc <- modtproc[,match(rownames(fidlist[[1]][[1]]), colnames(modtproc))]
modtproc <- modtproc %>% rownames_to_column(var="Gene")

fidsenslist <- list()
for(subclass in seq_along(leindegenes)){
  modttrim <- modtproc[,c(1,subclass+1)]
  colnames(modttrim)[2] <- "tval"
  fidsenslist[[subclass]] <- fidlist[[1]][[1]][subclass,] %>%
    dplyr::select(leindegenes[[subclass]]) %>% t %>% as.data.frame %>% 
    mutate("ct"=names(leindegenes)[subclass]) %>% 
    rownames_to_column(var="Gene") %>%
    left_join(modr2proc[[1]]) %>% 
    left_join(modttrim)
  colnames(fidsenslist[[subclass]])[2] <- "sens"
}
plotdf <- do.call(rbind, fidsenslist)
plotdford <- plotdf %>% group_by(ct) %>% summarise(meansens=mean(sens)) %>% arrange(meansens)
plotdf$ct <- factor(plotdf$ct, levels=plotdford$ct)

fidspeclist <- list()
for(subclass in seq_along(leindegenes)){
  modttrim <- modtproc[,c(1,subclass+1)]
  colnames(modttrim)[2] <- "tval"
  fidspeclist[[subclass]] <- fidlist[[1]][[2]][subclass,] %>%
    dplyr::select(leindegenes[[subclass]]) %>% t %>% as.data.frame %>% 
    mutate("ct"=names(leindegenes)[subclass]) %>% 
    rownames_to_column(var="Gene") %>%
    left_join(modr2proc[[1]]) %>% 
    left_join(modttrim)
  colnames(fidspeclist[[subclass]])[2] <- "spec"
}
plotdf2 <- do.call(rbind, fidspeclist)
#plotdford2 <- plotdf2 %>% group_by(ct) %>% summarise(meanspec=mean(spec)) %>% arrange(meanspec)
plotdf2$ct <- factor(plotdf2$ct, levels=plotdford$ct)

fidfidlist <- list()
for(subclass in seq_along(leindegenes)){
  modttrim <- modtproc[,c(1,subclass+1)]
  colnames(modttrim)[2] <- "tval"
  fidfidlist[[subclass]] <- fidlist[[1]][[3]][subclass,] %>%
    dplyr::select(leindegenes[[subclass]]) %>% t %>% as.data.frame %>% 
    mutate("ct"=names(leindegenes)[subclass]) %>% 
    rownames_to_column(var="Gene") %>%
    left_join(modr2proc[[1]]) %>% 
    left_join(modttrim)
  colnames(fidfidlist[[subclass]])[2] <- "fid"
}
plotdf3 <- do.call(rbind, fidfidlist)
#plotdford3 <- plotdf3 %>% group_by(ct) %>% summarise(meanfid=mean(fid)) %>% arrange(meanfid)
plotdf3$ct <- factor(plotdf3$ct, levels=plotdford$ct)

# dotplot of sens/spec/fid vs r2

p <- ggplot(plotdf, aes(x=sens, y=r2)) + 
  theme_bw() + 
  geom_point() + 
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Sensitivity\n(single cell)", y=bquote("Adjusted R"^2)) +
  facet_wrap(~ct, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/sens_vs_r2_leinDEgenes.png"), width=13, height=10)

p <- ggplot(plotdf2, aes(x=spec, y=r2)) + 
  theme_bw() + 
  geom_point() + 
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Specificity\n(single cell)", y=bquote("Adjusted R"^2)) +
  facet_wrap(~ct, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/spec_vs_r2_leinDEgenes.png"), width=13, height=10)

p <- ggplot(plotdf3, aes(x=fid, y=r2)) + 
  theme_bw() + 
  geom_point() + 
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Accuracy\n(Sensitivity x Specificity)", y=bquote("Adjusted R"^2)) +
  facet_wrap(~ct, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/fid_vs_r2_leinDEgenes.png"), width=13, height=10)

p <- ggplot(plotdf, aes(x=sens, y=tval)) + 
  theme_bw() + 
  geom_point() + 
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Sensitivity\n(single cell)", y=bquote("T-value")) +
  facet_wrap(~ct, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/sens_vs_tval_leinDEgenes.png"), width=13, height=10)

p <- ggplot(plotdf2, aes(x=spec, y=tval)) + 
  theme_bw() + 
  geom_point() + 
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Specificity\n(single cell)", y=bquote("T-value")) +
  facet_wrap(~ct, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/spec_vs_tval_leinDEgenes.png"), width=13, height=10)

p <- ggplot(plotdf3, aes(x=fid, y=tval)) + 
  theme_bw() + 
  geom_point() + 
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="Accuracy\n(Sensitivity x Specificity)", y=bquote("T-value")) +
  facet_wrap(~ct, nrow=4, ncol=6)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/fid_vs_tval_leinDEgenes.png"), width=13, height=10)


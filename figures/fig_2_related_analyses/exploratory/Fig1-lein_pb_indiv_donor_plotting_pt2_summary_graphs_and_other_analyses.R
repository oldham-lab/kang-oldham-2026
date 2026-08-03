source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))
library(cowplot)
library(tidyverse)
library(qs)
library(ggpubr)
options(bitmapType = 'cairo')

####################
# Load required data
####################
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)

#############################
# Load modeling data and plot
#############################
plotdflist_sub <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/fig1bdflist_subclass.RDS"))
plotdflist_super <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/fig1bdflist_supertype.RDS"))

# Calculate SANITY threshold for each of the three donors
genemeansfiltvec <- list()
for(i in 1:3){
  exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
  exdirexpr <- exdir[grep("EXPRLIST",exdir)]
  exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
  exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  
  genemeansfiltvec[[i]] <- do.call(rbind,mapply(function(x,y){
    obj <- readRDS(x) 
    obj <- norm_pb_samples(obj[[1]], index=3)
    exprt <- t(obj[,3:ncol(obj)])
    colnames(exprt) <- obj[,2]
    genemeansfilt <- apply(exprt,2,function(x) log10(mean(x)))
    return(data.frame("pcnt.var"=y,
                      "mean"=genemeansfilt[names(genemeansfilt)=="STAMBPL1"]))
  }, exdirexpr, exdirnames, SIMPLIFY=F))
  rownames(genemeansfiltvec[[i]]) <- NULL
  genemeansfiltvec[[i]] <- genemeansfiltvec[[i]] %>% dplyr::filter(pcnt.var==0)
}

#######
# How well are celltype-specific genes (DE genes, genes with sig t-values) modeled?
#######

# Load DE genes from Lein et al
origdemat <- fread(data.table=F, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein_orig_de_genes.csv"))
origdelist <- tapply(origdemat[,1], origdemat[,2], list)
allorigdegenes <- unique(unlist(origdelist))

# Collect gene symbols from each donor (pcnt.var=0)
genesymlist <- list()
for(i in 1:3){
  exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
  exdirexpr <- exdir[grep("EXPRLIST",exdir)]
  obj <- readRDS(exdirexpr[[5]]) 
  genesymlist[[i]] <- obj[[1]][,2]
}

# Add gene symbols to donor meanR2 mat, then format, then indicate which genes are DE genes
zerodonormat <- mapply(function(x,y,z){
  x %>% dplyr::filter(pcnt.var==0) %>%
    mutate(Gene=y, Donor=z)
}, plotdflist_sub, genesymlist,1:3, SIMPLIFY=F) %>% do.call(rbind, .) %>%
  mutate(yes_de=ifelse(Gene %in% allorigdegenes, "DE\ngenes", "non-DE\ngenes"), Donor=paste0("Donor ", Donor)) %>%
  mutate(yes_de=factor(yes_de,levels=c("non-DE\ngenes","DE\ngenes")))
# top genes are NPY, PLP1, ADARB2, SLC6A1

# plot distribution of modeling r2 of DE genes vs non DE genes
p <- ggplot(zerodonormat, aes(x=yes_de, y=adj_r2)) +
  theme_bw() + 
  geom_violin(aes(fill=yes_de),alpha=0.9) +
  geom_boxplot(notch=T,width=0.1) +
  theme(text=element_text(size=30),
        legend.position="none") +
  labs(x="", y=bquote("Adj R"^2)) +
  ylim(0,1) +
  facet_wrap(~Donor, nrow=3) +
  stat_compare_means(label.y=0.85)
  #stat_compare_means(label = "p.signif",method.args = list(alternative = "less"))

ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_DEvsnonDE.png"), width=4, height=9)


p <- ggplot(zerodonormat, aes(x=yes_de, y=mean)) +
  theme_bw() + 
  geom_violin(aes(fill=yes_de),alpha=0.9) +
  geom_boxplot(notch=T,width=0.1) +
  theme(text=element_text(size=30),
        legend.position="none") +
  labs(x="", y=bquote("Mean expression ("~log[10]~")")) +
  ylim(0,1) +
  facet_wrap(~Donor, nrow=3) +
  stat_compare_means(label.y=0.85)
  #stat_compare_means(label = "p.signif",method.args = list(alternative = "less"))

ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_DEvsnonDE_meanexpr.png"), width=4, height=9)

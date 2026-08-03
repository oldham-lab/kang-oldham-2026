source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
library(cowplot)
library(tidyverse)
library(qs)
options(bitmapType = 'cairo')

####################
# Load required data
####################
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)

#############################
# Load modeling data and plot
#############################

dirlist <- c(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_5050/"),
             file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_6040/"),
             file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_7030/"),
             file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/"))

meansvec_sub <- list()
meansvec_super <- list()
for(homedir in seq_along(dirlist)){ 
  plotdflist_sub <- readRDS(paste0(dirlist[homedir],"/fig1bdflist_subclass.RDS"))
  plotdflist_super <- readRDS(paste0(dirlist[homedir],"/fig1bdflist_supertype.RDS"))

  # # Calculate SANITY threshold for each of the three donors
  # genemeansfiltvec <- list()
  # for(i in 1:3){
  #   exdir <- list.files(paste0(homedir,"/donor",i,"/SyntheticDatasets/"),full.names=T)
  #   exdirexpr <- exdir[grep("EXPRLIST",exdir)]
  #   exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
  #   exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
    
  #   genemeansfiltvec[[i]] <- do.call(rbind,mapply(function(x,y){
  #     obj <- readRDS(x) 
  #     obj <- norm_pb_samples(obj[[1]], index=3)
  #     exprt <- t(obj[,3:ncol(obj)])
  #     colnames(exprt) <- obj[,2]
  #     genemeansfilt <- apply(exprt,2,function(x) log10(mean(x)))
  #     return(data.frame("pcnt.var"=y,
  #                       "mean"=genemeansfilt[names(genemeansfilt)=="STAMBPL1"]))
  #   }, exdirexpr, exdirnames, SIMPLIFY=F))
  #   rownames(genemeansfiltvec[[i]]) <- NULL
  #   genemeansfiltvec[[i]] <- genemeansfiltvec[[i]] %>% dplyr::filter(pcnt.var==0)
  # }

  # Load gene symbols for each donor
  genesymlist <- list()
  for(i in 1:3){
    exdir <- list.files(paste0(dirlist[homedir],"/donor",i,"/SyntheticDatasets/"),full.names=T)
    exdirexpr <- exdir[grep("EXPRLIST",exdir)]
    obj <- readRDS(exdirexpr[grep("0pcntVar", exdirexpr)]) 
    genesymlist[[i]] <- obj[[1]][,2]
  }

  # Calculate mean r2 for subclasses
  meansvec_sub[[homedir]] <- mapply(function(x, gene){
    df <- x %>% dplyr::filter(pcnt.var==0) %>%
      mutate("Gene"=gene)
    return(df)
  }, plotdflist_sub, genesymlist, SIMPLIFY=F)

  # Calculate mean r2 for supertypes
  meansvec_super[[homedir]] <- lapply(plotdflist_super, function(x){
    df <- x %>% dplyr::filter(pcnt.var==0)
  })
}


# Plot r2 for DE genes for different ratios
origdemat2 <- fread(data.table=F, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein_orig_de_genes.csv"))
degenes <- tapply(origdemat2$genes, origdemat2$celltype,list)
ratiovec <- c("50:50", "60:40", "70:30", "80:20")
r2df <- mapply(function(x,ratio){
  mapply(function(y,yname){
    data.frame("adj.r2"=x[[1]]$adj_r2[x[[1]]$Gene %in% y],
               "ct"=yname)
  },degenes,names(degenes),SIMPLIFY=F) %>%
    do.call(rbind,.) %>% mutate("ratio"=ratio)
}, meansvec_sub, ratiovec, SIMPLIFY=F) %>% do.call(rbind,.)

p <- ggplot(r2df, aes(x=ratio, y=adj.r2)) +
  theme_bw() + 
  geom_boxplot(notch=T) + 
  theme(text=element_text(size=30),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="Ratio of neurons to non-neurons", 
       y=bquote("Adjusted R"^2)) +
  facet_wrap(~ct, nrow=4, ncol=6)

ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/ratio_analysis_neurons_to_nonneurons.png"), width=13, height=10)

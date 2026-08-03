source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))
library(cowplot)
library(tidyverse)
library(qs)
options(bitmapType = 'cairo')

####################
# Load required data
####################
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_zeroing/")

#############################
# Load modeling data and plot
#############################
plotdflist_sub <- readRDS(paste0(homedir,"/fig1bdflist_subclass.RDS"))
plotdflist_super <- readRDS(paste0(homedir,"/fig1bdflist_supertype.RDS"))

# Calculate SANITY threshold for each of the three donors
sanmeanvec <- c(1.760, 1.709, 1.691)
genemeansfiltvec <- list()
for(i in 1:3){
  exdir <- list.files(paste0(homedir,"/donor",i,"/SyntheticDatasets/"),full.names=T)
  exdirexpr <- exdir[grep("EXPRLIST",exdir)]
  #exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
  #exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  exdirnames <- 0

  genemeansfiltvec[[i]] <- do.call(rbind,mapply(function(x,y){
    obj <- readRDS(x) 
    obj <- norm_pb_samples(obj[[1]], index=3)
    exprt <- t(obj[,3:ncol(obj)])
    colnames(exprt) <- obj[,2]
    genemeansfilt <- apply(exprt,2,function(x) log10(mean(x)))
    return(data.frame("pcnt.var"=y,
                      "mean"=sanmeanvec[i]))
  }, exdirexpr, exdirnames, SIMPLIFY=F))
  rownames(genemeansfiltvec[[i]]) <- NULL
  genemeansfiltvec[[i]] <- genemeansfiltvec[[i]] %>% dplyr::filter(pcnt.var==0)
}

# Calculate mean r2 for subclasses
meansvec_sub <- lapply(plotdflist_sub, function(x){
  df <- x %>% group_by(pcnt.var) %>%
    summarise(mean=mean(adj_r2)) %>%
    dplyr::filter(pcnt.var==0) %>% 
    mutate(mean_label=paste0("Mean adj. r2 = ",signif(mean,2)),
           pcnt.var.label=paste0("Pcnt. var = ", pcnt.var)) %>%
    mutate(pcnt.var=as.factor(pcnt.var.label))
})

# Calculate mean r2 for supertypes
meansvec_super <- lapply(plotdflist_super, function(x){
  df <- x %>% group_by(pcnt.var) %>%
    summarise(mean=mean(adj_r2)) %>%
    dplyr::filter(pcnt.var==0) %>%
    mutate(mean_label=paste0("Mean adj. r2 = ",signif(mean,2)),
           pcnt.var.label=paste0("Pcnt. var = ", pcnt.var)) %>%
    mutate(pcnt.var=as.factor(pcnt.var.label))
})

# Load DE genes from Lein et al
origdemat <- fread(data.table=F, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein_orig_de_genes.csv"))
origdelist <- tapply(origdemat[,1], origdemat[,2], list)
allorigdegenes <- unique(unlist(origdelist))

# Collect gene symbols from each donor (pcnt.var=0)
genesymlist <- list()
for(i in 1:3){
  exdir <- list.files(paste0(homedir,"/donor",i,"/SyntheticDatasets/"),full.names=T)
  exdirexpr <- exdir[grep("EXPRLIST",exdir)]
  obj <- readRDS(exdirexpr) 
  genesymlist[[i]] <- obj[[1]][,2]
}

# Add gene symbols to plotdflist objects
plotdflist_sub <- mapply(function(x,y,z){
  x %>% dplyr::filter(pcnt.var==0) %>%
    mutate(Gene=y)
}, plotdflist_sub, genesymlist,1:3, SIMPLIFY=F)

plotdflist_super <- mapply(function(x,y,z){
  x %>% dplyr::filter(pcnt.var==0) %>%
    mutate(Gene=y)
}, plotdflist_super, genesymlist,1:3, SIMPLIFY=F)

# Plot first two columns of figure 1
plot_cols_facetwrap <- function(a,b,c){
  meansvec <- do.call(rbind, b) %>%
    mutate(donor=c("Donor 1", "Donor 2", "Donor 3"))
  meansvec$mean_label <- gsub("adj. r2 ", "", meansvec$mean_label)
  out <- mapply(function(i, plotdflist){
    plotdflist[plotdflist$pcnt.var==0,] %>%
      mutate(donor=paste0("Donor ", i))
  }, c(1:3), a, SIMPLIFY=F) %>% do.call(rbind,.) %>%
    mutate(de_gene=Gene %in% allorigdegenes)
  
  out <- out %>%
    ggplot(aes(y=adj_r2, x=mean)) +
      theme_bw() +
      geom_point(alpha=0.05,size=0.5) +
      stat_density_2d(geom = "polygon", contour = TRUE,
                    aes(fill = after_stat(level)), colour = "black",
                    bins = 5) +
      scale_fill_distiller(palette = "Blues", direction = 1) +
      xlab(bquote("Mean expression ("~log[10]~")")) +
      ylab(bquote("Adjusted R"^2)) +
      geom_hline(data=meansvec, aes(yintercept=mean),linetype="dashed",color="red", linewidth=2, alpha=0.5) +
      ylim(0,1) +
      ggtitle(c) +
      geom_vline(data=genemeansfiltvec[[i]],aes(xintercept=mean),linetype="dashed",color="blue",linewidth=2, alpha=0.5) +
      geom_text(data=meansvec,aes(label=mean_label, y=mean+0.1), x = -1, color="red", size=8, alpha=0.8) +
      theme(text=element_text(size=30),
            legend.position="none",
            axis.title.x=element_text(size=30,margin=margin(10,0,0,0)),
            axis.title.y=element_text(margin=margin(0,10,0,0)),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=24),
            plot.margin=margin(1,1,1,1,"cm")) +
      facet_wrap(~donor, nrow=3)
  return(out)
}

plot_cols_facetwrap_deGeneColor <- function(a,b,c){
  meansvec <- do.call(rbind, b) %>%
    mutate(donor=c("Donor 1", "Donor 2", "Donor 3"))
  meansvec$mean_label <- gsub("adj. r2 ", "", meansvec$mean_label)
  out <- mapply(function(i, plotdflist){
    plotdflist[plotdflist$pcnt.var==0,] %>%
      mutate(donor=paste0("Donor ", i))
  }, c(1:3), a, SIMPLIFY=F) %>% do.call(rbind,.) %>%
    mutate(de_gene=Gene %in% allorigdegenes)
  
  out <- out %>%
    ggplot(aes(y=adj_r2, x=mean, color=de_gene)) +
      theme_bw() +
      geom_point(alpha=0.05,size=0.5) +
      geom_point() + 
      xlab(bquote("Mean expression ("~log[10]~")")) +
      ylab(bquote("Adjusted R"^2)) +
      geom_hline(data=meansvec, aes(yintercept=mean),linetype="dashed",color="red", linewidth=2, alpha=0.5) +
      ylim(0,1) +
      ggtitle(c) +
      geom_vline(data=genemeansfiltvec[[i]],aes(xintercept=mean),linetype="dashed",color="blue",linewidth=2, alpha=0.5) +
      geom_text(data=meansvec,aes(label=mean_label, y=mean+0.1), x = -1, color="red", size=8, alpha=0.8) +
      theme(text=element_text(size=30),
            legend.position="none",
            axis.title.x=element_text(size=30,margin=margin(10,0,0,0)),
            axis.title.y=element_text(margin=margin(0,10,0,0)),
            plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=24),
            plot.margin=margin(1,1,1,1,"cm")) +
      facet_wrap(~donor, nrow=3)
  return(out)
}

p1 <- plot_cols_facetwrap(plotdflist_sub, meansvec_sub,"Subclass (n=24)")
p2 <- plot_cols_facetwrap(plotdflist_super, meansvec_super, "Supertype (n=136)")

fig1p1 <- plot_grid(p1, p2,
                  ncol=2,
                  labels=c("B","C"),
                  label_size=42)
png(file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/zeroing/fig1_part1.png"),width=1000,height=1000)
plot(fig1p1)
dev.off()

p1 <- plot_cols_facetwrap_deGeneColor(plotdflist_sub, meansvec_sub,"Subclass (n=24)")
p2 <- plot_cols_facetwrap_deGeneColor(plotdflist_super, meansvec_super, "Supertype (n=136)")

fig1p1 <- plot_grid(p1, p2,
                  ncol=2,
                  labels=c("B","C"),
                  label_size=42)
png(file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/zeroing/fig1_part1_colorByDE.png"),width=1000,height=1000)
plot(fig1p1)
dev.off()

# Plot markers (one that is well modeled, two that are not)
fid <- fread(data.table=F,file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "bulk_fidelity/bulkFidelity_Hs.FCX_AOMN.csv"))
ftop20 <- lapply(2:5, function(i){
  df <- fid[order(fid[,i], decreasing=T),]
  return(df[1:20,1])
})

# Gather all gene models from pcnt.var=0
exdirlist <- list()
exprZero <- list()
for(i in 1:3){
  exdir <- list.files(paste0(homedir,"/donor",i,"/SyntheticDatasets/"),full.names=T)
  exdirexpr <- exdir[grep("EXPRLIST",exdir)]
  exdirexpr <- exdirexpr[grep("0pcntVar", exdirexpr)]
  exdirsif <- exdir[grep("samples_legend",exdir)]
  exdirsif <- exdirsif[grep("0pcntVar", exdirsif)]
  #exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
  #exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  exdirnames <- 0

  obj <- readRDS(exdirexpr) 
  obj <- norm_pb_samples(obj[[1]], index=3)
  exprt <- t(obj[,3:ncol(obj)])
  escale <- apply(exprt, 2, scale)
  rownames(escale) <- rownames(exprt)
  colnames(escale) <- genemap[colnames(exprt),2]
  exprZero[[i]] <- escale
  #genemeansfilt <- apply(exprt,2,function(x) log10(mean(x)))
  sif <- fread(data.table=F,file=exdirsif)
  sampct <- t(calc_CT_counts(sif, cell_anno, ct_col=1, cell_id_col = 3))
  sampctmean <- apply(sampct, 2, mean)
  sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
    
  # Run modeling
  exdirlist[[i]] <- lapply(as.list(as.data.frame(escale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
}

# Create graphs for rightmost columns of figure 1
graph_marker_model <- function(gene){
  plotdf <- do.call(rbind,mapply(function(x,y,z){
    index <- which(colnames(y)==gene)
    df <- data.frame("pred"=predict(x[[index]]),"actual"=y[,index], "gene"=gene, "donor"=paste0("Donor ",z))
  }, exdirlist, exprZero,seq_along(exdirlist), SIMPLIFY=F))
  return(plotdf)
}

plot_gene_plots_facetwrap <- function(x){
  df <- graph_marker_model(x) 
  dfcor <- df %>% group_by(donor) %>% summarise(cors=cor(pred,actual)[,1]) %>%
    mutate(corlabel=paste0("~R^{2} == ", signif(cors^2,2)))

  p <- ggplot(df, aes(x=pred, y=actual)) + 
    theme_bw() +
    geom_density_2d_filled() +
    geom_abline() +
    theme(text=element_text(size=30),
          #axis.title.x=element_text(margin=margin(10,0,0,0)),
          #axis.title.y=element_text(margin=margin(0,10,0,0)),
          plot.subtitle=element_text(size=20, color="red"),
          strip.text=element_text(size=24),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position="none") +
          #plot.margin=margin(1,1,1,1,"cm")) +
    labs(title=x, x="", y="") +
    geom_text(data=dfcor,aes(label=corlabel), y=Inf, x = -Inf,fontface="bold", color="white", size=12, 
      alpha=0.8, vjust=2,hjust=-0.1,parse=T) +
    facet_wrap(~donor, nrow=3)
  return(p)
}

inputgenevec1 <- c("MOG", "ALDH1L1", "TYROBP", "VIP", "VAT1L")
pgenelist1 <- lapply(inputgenevec1, plot_gene_plots_facetwrap)

title1 <- ggdraw() + draw_label("\tPredicted expression",size=42)
title2 <- ggdraw() + draw_label("Actual expression",size=42, angle=90)

pglist1 <- plot_grid(plotlist=pgenelist1, ncol=length(pgenelist1))
fig1p2 <- plot_grid(pglist1, 
                  ncol=1,
                  label_size=42)
fig1p2 <- plot_grid(fig1p2,title1,NULL, nrow=3, rel_heights=c(1,0.005,0.05))
fig1p2 <- plot_grid(title2, fig1p2, ncol=2, rel_widths=c(0.03,1))
png(file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/zeroing/fig1_part2.png"),width=1800,height=1200)
plot(fig1p2)
dev.off()





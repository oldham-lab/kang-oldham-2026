source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))
library(cowplot)
library(tidyverse)
library(qs)
library(ggrepel)
options(bitmapType = 'cairo') 
library(doParallel)
registerDoParallel(cores=8)

####################
# Load required data (SAVER)
####################
save_dirall <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/")
cell_annoall <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
genes_of_interest <- c("AIF1","ALDH1L1","MOG","SLC17A7","GAD1")
cell_exprall <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/impute/SAVER/output_aif1_aldh1l1_mog_slc17a7_gad1.qs"))
cell_exprall <- cell_exprall[rownames(cell_exprall) %in% genes_of_interest,]
megaexpr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
orig_genes <- rownames(cell_exprall)
donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")

# genes expressed in less than 3 cells
#bla <- apply(cell_exprall,1,function(x) sum(x>0))
#low_genes <- names(bla)[bla<3]
#qsave(low_genes, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/genes_expressed_in_lessThan3_cells.qs"))

#san_mean <- apply(cell_exprall,1,mean)
# which.min(san_mean[san_mean>=1])
# STAMBPL1 
# 3262 
#qsave(san_mean,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means.qs"))
# other genes around STAMBPL1 (in order):
# [1] "RNF213-AS1" "PABPC1"     "ZNF654"     "C12orf40"   "LINC01331" 
# [6] "STK4"       "CCND2"      "STAMBPL1"   "TNK2"       "FAM78B"    
#[11] "SLC1A1"
# calculate means by donor
san_mean_list <- list()
for(i in seq_along(donorvec)){
  df_donor <- cell_exprall[,colnames(cell_exprall) %in% cell_annoall$Cell_ID[cell_annoall$Donor==donorvec[i]]]
  san_mean_list[[i]] <- apply(df_donor, 1, mean)
}
qsave(san_mean_list,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means_byDonor_SAVER.qs"))

###############################################
# create pseudobulk from individual lein donors (DFC)
###############################################

create_donorspecific_pseudobulk(cell_exprall = cell_exprall,
                                cell_annoall = cell_annoall,
                                donorvec = c("H18.30.002", "H19.30.001", "H19.30.002"), # vector of all unique donors
                                save_dirall = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_SAVER/"),
                                pcnt.varvec = c(0), # vector of values for pcnt.var
                                no.samples = ncol(megaexpr)-2,
                                cell.nameindex = 3, # column index of cell_annoall that indicates cell name (matches colnames of cell_exprall)
                                filter = megaexpr[,1] # vector of gene symbols to filter pseudobulk matrix 
                                )

#################################################
# For each donor:
# Load Lein pseudobulk 
# Normalize by sample and scale by gene
# Model using subclasses and supertypes
# DLPFC
#################################################

model_donorspecific_pseudobulk(dirlist = c(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_SAVER/")),
                               cell_annoall=cell_annoall,
                               cellannocolvec = c(1,4), # subclass, supertype
                               cellIDcol = 3, # new_labs, not Cell_ID
                               celltypenames = c("subclass", "supertype"))

 
######
# Load modeling result and compare to pre-scVI 
######
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))

homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/")
plotdflist_sub <- qread(paste0(homedir,"/fig1bdflist_subclass.qs"))
plotdflist_super <- qread(paste0(homedir,"/fig1bdflist_supertype.qs"))
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")

# Collect gene symbols from each donor (pcnt.var=0)
 genesymlist <- list()
 for(i in 1:3){
      exdir <- list.files(paste0(homedir,"/donor",i,"/SyntheticDatasets/"),full.names=T)
      exdirexpr <- exdir[grep("EXPRLIST",exdir)]
      exdirexpr <- exdirexpr[grep("0pcntVar", exdirexpr)]
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
 
homedir_saver <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_SAVER/")
saver_sub <- qread(paste0(homedir_saver,"/fig1bdflist_subclass.qs"))
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means_byDonor.qs"))

plotdf_sub1 <- mapply(function(x,y,z,a){
  out <- left_join(x,y[,c(1,6)],by=join_by("Gene"=="Gene")) %>%
     mutate("Donor"=z)
  a <- data.frame("san_mean"=a) %>% rownames_to_column(var="Gene")
  out <- left_join(out, a, by=join_by("Gene"=="Gene"))
  return(out)
 }, saver_sub, plotdflist_sub, donornames,san_mean, SIMPLIFY=F) %>% do.call(rbind,.) %>%
   mutate(adj_r2_diff = adj_r2.x-adj_r2.y, san_mean=log2(san_mean))
  
plotdf_sub <- mapply(function(x,y,z,a){
  y <- y[y$Gene %in% x$Gene,]
  x <- x[x$Gene %in% y$Gene,]
  x$mean <- y$mean[match(x$Gene,y$Gene)]
  out <- rbind(x %>% mutate("type"="SAVER"), y %>% mutate("type"="Original")) %>% 
    mutate(Donor=z, type=factor(type,levels=c("Original","SAVER")))
  a <- data.frame("san_mean"=a) %>% rownames_to_column(var="Gene")
  out <- left_join(out, a, by=join_by("Gene"=="Gene"))
  return(out)
}, saver_sub, plotdflist_sub, donornames,san_mean, SIMPLIFY=F) %>% do.call(rbind,.)
  
# Plot adj r2 vs mean expr
p <- ggplot(plotdf_sub1, aes(x=adj_r2_diff, y=san_mean)) + 
  theme_bw() +
  geom_point(aes(color=type)) +
  geom_smooth() +
  facet_wrap(~Donor) +
  labs(x=bquote(Delta~" Adjusted R"^2~" (scVI - Original)"), y=bquote("Mean expression ("~log[2]~")")) +
  theme(text=element_text(size=30),
        legend.position="none") +
  scale_x_continuous(labels = function(x) format(x, nsmall = 2)) 
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/SAVER/adjr2_vs_mean.png"), width=14)

# Plot after filtering to genes of interest
p <- plotdf_sub %>% 
  ggplot(aes(x=type, y=adj_r2)) +
    theme_bw() +
    geom_violin(aes(fill=type), alpha=0.2) + 
    geom_line(aes(group=Gene), alpha=0.2) +
    geom_label_repel(aes(label=Gene), max.overlaps=15) +
    geom_point() +
    facet_wrap(~Donor) +
    labs(x="", y=bquote("Adjusted R"^2)) +
    theme(text=element_text(size=30),
          legend.position="none")
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/SAVER/genes_of_interest.png"),width=10)

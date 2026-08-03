# A code document for follow-up analyses after producing figure 1 using CompareMarkers

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
source(file.path(Sys.getenv("PSEUDOBULK_DIR", "/home/gugene/code/git/Pseudobulk-from-SC-SN-data"), "makeSyntheticDatasets_0.51.r"))
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2/fxns.R"))
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

######################################################
# Calculate greedy mods from pseudobulk for each donor
######################################################

# source("/home/gugene/code/git/COPA/kme_topmodposbc.R")
# source("/home/gugene/code/git/COPA/greedy_march_COPA_20240802.R")
# library(CoPA)
# for(i in 1:3){
#   exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
#   exdirexpr <- exdir[grep("EXPRLIST",exdir)]
#   exdirsif <- exdir[grep("samples_legend",exdir)]
#   exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
#   exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  
#   x <- exdirexpr[5]
#   obj <- readRDS(x) 
#   expr <- norm_pb_samples(obj[[1]], index=3)
#   save_dir1 <- paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/donor"),i,"_greedy/")
#   find_unmerged_greedy_mods(minSizevec=10)
# }

# # Calculate topmodposbc genes for each donor network
# for(i in 1:3){
#   exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
#   exdirexpr <- exdir[grep("EXPRLIST",exdir)]
#   exdirsif <- exdir[grep("samples_legend",exdir)]
#   exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
#   exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  
#   x <- exdirexpr[5]
#   obj <- readRDS(x) 
#   expr <- norm_pb_samples(obj[[1]], index=3)
#   save_dir1 <- paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/donor"),i,"_greedy/")
#   produce_kme_topmodposbc_tables()
# }

############################################################################################
# Create tables of genes positively and significantly associated with each subclass/supertype
############################################################################################

for(ze in c(1, 4)){
  if(ze==1){
    ctstring <- "Subclass"
  } else if(ze==4){
    ctstring <- "Supertype"
  } else {
    ctstring < ""
  }
  modlist <- list()
  modlistU <- list()
  for(i in 1:3){
    exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
    exdirexpr <- exdir[grep("EXPRLIST",exdir)]
    exdirsif <- exdir[grep("samples_legend",exdir)]
    exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
    exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  
    x <- exdirexpr[5]
    y <- exdirsif[5]
    z <- exdirnames[5]
    obj <- readRDS(x) 
    obj <- norm_pb_samples(obj[[1]], index=3)
    exprt <- t(obj[,3:ncol(obj)])
    escale <- apply(exprt, 2, scale)
    rownames(escale) <- rownames(exprt)
    
    sif <- fread(data.table=F,file=y)
    sampct <- t(calc_CT_counts(sif, cell_anno, ct_col=ze, cell_id_col = 3))
    sampctmean <- apply(sampct, 2, mean)
    sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
    
    # Run modeling
     
    mod_list <- foreach(l=1:ncol(escale)) %dopar% {
      lm(escale[,l] ~ ., data=as.data.frame(sampct))
    }
    modt <- foreach(l=seq_along(mod_list)) %dopar% {
      summary(mod_list[[l]])$coefficients[,3:4]
    }
    names(modt) <- obj$Gene
    rm(mod_list)
    gc()
 
    # Collect which genes have positive t-values
    modtlistbool <- lapply(modt, function(x) x[,1]>0)
    # Collect p-values in list
    modplist <- lapply(modt, function(x) x[,2])
    # Order p-values
    modpord <- do.call(cbind, modplist) %>% t %>% as.data.frame %>% as.list %>% 
      lapply(., function(x) data.frame("gene"=names(modplist), "pval"=x) %>% arrange(pval))
    # Calculate FDR cutoff
    modq <- qvalue(unlist(modplist))
    modqcut <- max(modq$pvalues[modq$qvalues < 0.05])
    # Collect which genes have p-values below FDR cutoff
    modplistbool <- lapply(modt, function(x) x[,2] < modqcut)
    # Collect which genes have positive t-values and p-values below FDR cutoff
    modgeneboth <- mapply(function(t, p){
        t & p
    }, modtlistbool, modplistbool, SIMPLIFY=F) %>% as.data.frame %>% t
    # Collect whether gene is max in a given celltype
    modpdf <- modplist %>% as.data.frame %>%
      apply(., 2, function(x) x == min(x)) %>% t
    # Collect "topmodposfdr" genes
    modgenebothU <- mapply(function(a,b){
        a & b
    }, as.list(as.data.frame(modgeneboth)), as.list(as.data.frame(modpdf)), SIMPLIFY=F ) %>% as.data.frame 
    colnames(modgenebothU) <- colnames(modgeneboth)
    # Create list of genes for each subclass
    modlist[[i]] <- apply(modgeneboth,2, function(x) obj[x,2])
    modlist[[i]] <- mapply(function(x,y) y$gene[y$gene %in% x], modlist[[i]], modpord, SIMPLIFY=F)
    modlistU[[i]] <- apply(modgenebothU,2, function(x) obj[x,2])
    modlistU[[i]] <- mapply(function(x,y) y$gene[y$gene %in% x], modlistU[[i]], modpord, SIMPLIFY=F)
  } # for i
  qsave(modlist, file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_FDR_"),ctstring,".qs"))
  qsave(modlistU, file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_topmodposFDR_"),ctstring,".qs"))
}
 
# Load gene symbol vector from pseudobulk object
genevecobj <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor2/SyntheticDatasets//SyntheticDataset1_10pcntCells_0pcntVar_1518samples_06-06-50_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS"))[[1]][,2]

# For each donor, plot enrichment vs # of DE genes from PB modeling
  # Gene set enrichment function
GSHG_custom <- function(allModules, # modules to be analyzed
                 mySets, # a list of all genesets
                 allgenes # all genes in dataset being analyzed
                 ){
  fisherTest <- function(set,mod,all){
    total.shared <- length(intersect(all,set))
    shared.in.mod <- length(intersect(mod,set))
    shared.out.mod <- total.shared-shared.in.mod
    in.mod.not.shared <- length(mod)-shared.in.mod
    out.mod.not.shared <- length(all)-length(mod)-shared.out.mod
    fisher.test(matrix(c(shared.in.mod,in.mod.not.shared,shared.out.mod,out.mod.not.shared),ncol=2),
                alternative="greater")$p.val
  }  
  
  allnetGenes <- unique(unlist(allModules))
  checkOverlap=function(x,y){
    length(intersect(x,y))
  }
  
  intsct <- sapply(mySets,checkOverlap,allnetGenes)
  mySets <- mySets[intsct > 0]
 # cat("Kept",length(mySets), "sets")
  
  mySetNames <- names(mySets)
  
  GSHGresults=matrix(nrow=length(mySets),
                     ncol=length(allModules),
                     data=-8888)
  
  for(i in c(1:length(allModules))){
    allModGenes <- unlist(allModules[[i]])
    GSHGresults[,i] <- unlist(lapply(mySets,
                                     function(aSet)
                                       fisherTest(aSet,allModGenes,all=allgenes)))
  #  print(paste("Finished module",i))
  }

  for(pla in seq_along(mySets)){
    fisherTest(mySets[[pla]], allModGenes, all=allgenes)
  }

  colnames(GSHGresults) <- names(allModules)
  datout <- data.frame("sets"=mySetNames,
                          GSHGresults)
  return(datout)
}


modlist_paths <- c(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_FDR_Subclass.qs"),
                   file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_FDR_Supertype.qs"),
                   file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_topmodposFDR_Subclass.qs"),
                   file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_topmodposFDR_Supertype.qs"))
modlist_desc <- c("FDR_Subclass", "FDR_Supertype", "topmodposFDR_Subclass", "topmodposFDR_Supertype")
for(z in c(1,3)){
  modlist <- qread(modlist_paths[z])
  moddesc <- modlist_desc[z]
  if(z %in% c(1,3)){
    ct_freq <- data.frame(table(cell_anno$Cell_Type))
  } else {
    ct_freq <- data.frame(table(cell_anno$Cluster))
  }
for(i in 1:3){
  mods <- fread(data.table=F,file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/donor"),i,"_greedy/kme_tables/topmodposbc_table.csv"))
  mods <- tapply(mods$Gene, mods$topmodposfdr, list)
  # 874 modules

  lein_gshg <- GSHG_custom(mods,
                           modlist[[i]],
                           allgenes=genevecobj)

  gshgt <- t(lein_gshg[,-1])
  sigmodlistind <- as.list(as.data.frame(apply(gshgt,2, order)))
  wminp <- apply(gshgt,2,which.min)
      
  find_unique_inds <- function(wminp){
    minp <- mapply(function(x,y) x[y], as.list(as.data.frame(gshgt)), wminp)
    unique_inds <- c()
    for(i in 1:length(wminp)){
      if(length(sigmodlistind[[i]])<=1){
        unique_inds[i] <- wminp[i]
      } else {
        these <- which(wminp == wminp[i])
      smallest <- which.min(minp[these])
        if(i==these[smallest]){
          unique_inds[i] <- wminp[i]
        } else {
          unique_inds[i] <- sigmodlistind[[i]][which(sigmodlistind[[i]]==wminp[i])+1]
        }
        } 
      }
      return(unique_inds)
    }
      
    new_inds <- find_unique_inds(wminp)
    old_inds <- wminp
    while(sum(old_inds==new_inds)<length(old_inds)){
      old_inds <- new_inds
      new_inds <- find_unique_inds(new_inds)
    }
    moduplength <- unlist(lapply(mods,length))[new_inds]
    pvls <- mapply(function(x,y) x[y], as.list(as.data.frame(gshgt)), new_inds)    
    mod_counts <- unlist(lapply(modlist[[i]], length))
    de_df <- data.frame("CT" = lein_gshg[,1],
                        "values" = pvls,
                        "which" = new_inds,
                        "mod_length" = moduplength,
                        "de_count" = mod_counts[names(mod_counts) %in% lein_gshg[,1]],
                        "dup" = duplicated(new_inds)) %>%
      dplyr::filter(CT!="(Intercept)") %>%
      mutate(CT = gsub("`","",CT)) %>%
      left_join(ct_freq, by=c("CT"="Var1")) %>%
      mutate(Celltype_freq=Freq/sum(Freq))

  p <- de_df %>%
    ggplot(aes(x=de_count, y=-log10(values))) +
      geom_point(aes(size=Celltype_freq)) + 
      geom_label_repel(aes(label=CT), size=5) +
      theme_classic() +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            plot.title=element_text(hjust=0.5,size=16),
            legend.title=element_text(size=12),
            plot.subtitle=element_text(hjust=0.5,size=14),
            legend.text=element_text(size=12)) +
      geom_smooth(alpha=0.2) +
      labs(title="Geneset enrichment of pseudobulk modules",
           subtitle=paste0("Donor ", i, ", r=",signif(cor(de_df$de_count, -log10(de_df$values)),2)),
            x="Number of DE genes (from pseudobulk modeling)",
            y="Enrichment in pseudobulk greedy modules\n-log10(p-value)") 

  ggsave(p, file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup_donor"),i,"_",moddesc, ".png"),width=12,height=8)

  p2 <- de_df %>%
    dplyr::filter(CT != "L2/3 IT") %>%
    ggplot(aes(x=de_count, y=Celltype_freq)) +
      geom_point() + 
      geom_label_repel(aes(label=CT), size=5) +
      theme_classic() +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            plot.title=element_text(hjust=0.5,size=16),
            legend.title=element_text(size=12),
            plot.subtitle=element_text(hjust=0.5,size=14),
            legend.text=element_text(size=12)) +
      geom_smooth(alpha=0.2) +
      labs(title="Geneset enrichment of pseudobulk modules",
           subtitle=paste0("Donor ", i, ", r=",signif(cor(de_df$de_count, -log10(de_df$values)),2)),
            x="Number of DE genes (from pseudobulk modeling)",
            y="Frequency of cell type in original data") 

  ggsave(p2, file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/numberofDEgenes_vs_enrichment/fig1_followup_donor"),i,"_",moddesc, "_countVsfreq.png"),width=12,height=8)

    p3 <- de_df %>%
    dplyr::filter(CT != "L2/3 IT") %>%
    ggplot(aes(x=-log10(values), y=Celltype_freq)) +
      geom_point() + 
      geom_label_repel(aes(label=CT), size=5) +
      theme_classic() +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            plot.title=element_text(hjust=0.5,size=16),
            legend.title=element_text(size=12),
            plot.subtitle=element_text(hjust=0.5,size=14),
            legend.text=element_text(size=12)) +
      geom_smooth(alpha=0.2) +
      labs(title="Geneset enrichment of pseudobulk modules",
           subtitle=paste0("Donor ", i, ", r=",signif(cor(de_df$de_count, -log10(de_df$values)),2)),
            x="Enrichment p-value (-log10)",
            y="Frequency of cell type in original data") 

  ggsave(p3, file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup_donor"),i,"_",moddesc, "_pvalVsfreq.png"),width=12,height=8)
}
}

###########
# Testing
###########

# Load model results
modelres <- qread(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/models_topmodposFDR_Subclass.qs"))

#### Why L5ET such an outlier?
# Get fdr genes for L5ET
othert <- list()
for(i in seq_along(modelres[[1]])){
  othertvals <- c()
  for(j in seq_along(modelres[[1]][[i]])){
    othertvals <- c(othertvals, modt[[which(names(modt)==modelres[[1]][[i]][[j]])]][-i,1])
  }
  othert[[i]] <- data.frame("ct"=names(modelres[[1]])[i],"other_tvals"=othertvals)
}
othert <- do.call(rbind, othert) %>% as.data.frame
othert$ct <- gsub("`", "", othert$ct)
othert <- othert %>% dplyr::filter(ct!="(Intercept)")
othert2 <- left_join(othert,de_df,by=join_by("ct"=="CT"))
sortbymean <- othert2 %>% group_by(ct) %>% summarise("mean"=mean(other_tvals)) %>% arrange(mean)
othert2$ct <- factor(othert2$ct, levels=sortbymean$ct)
p <- ggplot(othert2,aes(x=other_tvals, y=ct, color=de_count)) +
  geom_boxplot() +
  xlab("Modeling t-values") +
  ylab("") +
  theme(axis.text.y=element_text(size=24),
        plot.title=element_text(size=24,hjust=0.5),
        axis.title.x=element_text(size=24),
        axis.text.x=element_text(size=16)) + 
  scale_color_gradient2() +
  ggtitle("Modeling t-values of a given cell-type's marker genes\n(topmodposFDR) by other CT abundance vectors")
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/tvals_of_other_cts_for_a_given_cts_marker_genes.png"),width=12,height=9)

# How well are cts named after markers modeled by those markers?
graph_marker_enrichment_in_ct <- function(markervec, 
                                          #cthighlightvec, 
                                          save_file){
  mvlist <- list()
  for(i in seq_along(markervec)){
    templist <- list()
    for(j in 1:3){
      modt <- modtlist[[j]]
      this_df <- modt[[which(names(modt)==markervec[i])]]
      templist[[j]] <- data.frame("ct"= rownames(this_df),
                              "tval"= this_df[,1],
                              "marker"=markervec[i],
                              "donor"=j)
    }
    mvlist[[i]] <- do.call(rbind, templist) %>% group_by(ct) %>% summarise("mean_tval"=mean(tval), "se_tval"=(sd(tval)/sqrt(3)), "marker"=unique(marker)) %>% as.data.frame
    mvlist[[i]]$ct <- gsub("`","",mvlist[[i]]$ct)
    #mvlist[[i]]$highlight <- mvlist[[i]]$ct == cthighlightvec[i]
  }
  mvlist <- do.call(rbind, mvlist) %>% as.data.frame %>%
    mutate(marker=factor(marker, levels=markervec))

  p <- ggplot(mvlist, aes(x=ct,y=mean_tval))+
    theme_bw() +
    geom_bar(stat="identity") +
    geom_errorbar(aes(x=ct,ymin=mean_tval-2*se_tval, ymax=mean_tval+2*se_tval), width=0.3, alpha=0.5) +
    ylab("Modeling t-value") +
    xlab("") +
    geom_vline(xintercept=0,color="black") +
    #ggtitle("Modeling t-values for\ncell-type marker genes of interest") +
    scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = "none" ) +
    theme(axis.title.y=element_text(size=24),
          plot.title=element_text(size=24,hjust=0.5),
          axis.title.x=element_text(size=24),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=14,hjust=1,vjust=1,angle=45),
          strip.text=element_text(size=24),
          panel.spacing = unit(1, "lines")) +
    facet_wrap(~marker, nrow=length(markervec), scale="free_y") 

  ggsave(p,file=save_file,height=12,width=10)
}

# Load t-values for all 3 donors
  ctstring <- "Subclass"
  ze=1
  modtlist <- list()
  for(i in 1:3){
    exdir <- list.files(paste0(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/donor"),i,"/SyntheticDatasets/"),full.names=T)
    exdirexpr <- exdir[grep("EXPRLIST",exdir)]
    exdirsif <- exdir[grep("samples_legend",exdir)]
    exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
    exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
  
    x <- exdirexpr[5]
    y <- exdirsif[5]
    z <- exdirnames[5]
    obj <- readRDS(x) 
    obj <- norm_pb_samples(obj[[1]], index=3)
    exprt <- t(obj[,3:ncol(obj)])
    escale <- apply(exprt, 2, scale)
    rownames(escale) <- rownames(exprt)
    
    sif <- fread(data.table=F,file=y)
    sampct <- t(calc_CT_counts(sif, cell_anno, ct_col=ze, cell_id_col = 3))
    sampctmean <- apply(sampct, 2, mean)
    sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
    
    # Run modeling
     
    mod_list <- foreach(l=1:ncol(escale)) %dopar% {
      lm(escale[,l] ~ ., data=as.data.frame(sampct))
    }
    modtlist[[i]] <- foreach(l=seq_along(mod_list)) %dopar% {
      summary(mod_list[[l]])$coefficients[,3:4]
    }
    names(modtlist[[i]]) <- obj$Gene
    rm(mod_list)
    gc()
  }


markervec <- c("SST", "PVALB", "VIP", "LAMP5", "LHX6", "PAX6") # CAR3 is not even in the list, investigate this
#cthighlightvec <- c("Sst", "Pvalb", "Vip", "Lamp5", "Lamp5 Lhx6", "Pax6")
graph_marker_enrichment_in_ct(markervec, #CAR3 is not even in the list, investigate this
                              #cthighlightvec,
                              save_file = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/tvals_of_marker_genes_of_interest.png"))
# None of these appear to be very specific

# Sanity check using glial markers
# How well are cts named after markers modeled by those markers?

graph_marker_enrichment_in_ct(markervec = c("AIF1","TYROBP", "MOG", "NTRK2", "PON2", "GFAP"), # CAR3 is not even in the list, investigate this
                              #cthighlightvec = c("Micro/PVM","Micro/PVM", "Oligo", "Astro", "Astro", "Astro"),
                              save_file = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/tvals_of_marker_genes_of_interest_pt2.png"))
# very odd results

graph_marker_enrichment_in_ct(markervec = c("MOG", "AIF1", "GFAP", "GAD1", "SLC17A7"), # CAR3 is not even in the list, investigate this
                              #cthighlightvec = c("Oligo","Micro/PVM", "Astro", "Pvalb", "L2/3 IT"),
                              save_file = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/tvals_of_marker_genes_of_interest_pt3.png"))

# project L5ET genes onto lein dfc
allmlist <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_summary_tables/allmlist_log.qs"))

markervec <-  head(modelres[[1]][[22]],10) # donor1, L5ET
cthighlightvec <- rep("L5 ET", 10)

  mvlist <- list()
  for(i in seq_along(markervec)){
    this_df <- allmlist[[1]][which(rownames(allmlist[[1]])==markervec[i]),]
    mvlist[[i]] <- data.frame("ct"= names(this_df),
                              "tval"= this_df,
                              "marker"=markervec[i])
    mvlist[[i]]$ct <- gsub("`","",mvlist[[i]]$ct)
    mvlist[[i]]$highlight <- mvlist[[i]]$ct == cthighlightvec[i]
  }
  mvlist <- do.call(rbind, mvlist) %>% as.data.frame %>% mutate(marker=factor(marker,levels=markervec))

  p <- ggplot(mvlist, aes(x=tval,y=ct, fill=highlight))+
    theme_bw() +
    geom_bar(stat="identity") +
    xlab("Expression (log10)") +
    ylab("") +
    geom_vline(xintercept=0,color="black") +
    ggtitle("Expression (log10) for\ncell-type marker genes of interest") +
    scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = "none" ) +
    theme(axis.title.y=element_text(size=24),
          plot.title=element_text(size=24,hjust=0.5),
          axis.title.x=element_text(size=24),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=14,hjust=1,vjust=1,angle=45),
          strip.text=element_text(size=24),
          panel.spacing = unit(1, "lines")) +
    coord_flip() +
    facet_wrap(~marker, nrow=length(markervec), scale="free_y") 

  ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/expression_L5ETmarkerGenes.png"),height=16,width=10)

#####################
# project mean of all ct markers onto lein dfc
######################
markervec <- modelres[[1]] # donor1, L5ET
cthighlightvec <- gsub("`","",names(modelres[[1]]))
cthighlightvec[1] <- "Sst Chodl"

  mvlist <- list()
  for(i in seq_along(markervec)){
    this_df <- allmlist[[1]][which(rownames(allmlist[[1]]) %in% markervec[[i]]),]
    if(!is.null(nrow(this_df))){
      this_df <- colMeans(this_df)
    }
    mvlist[[i]] <- data.frame("ct"= names(this_df),
                              "tval"= this_df,
                              "ct_genes"=cthighlightvec[i])
    mvlist[[i]]$highlight <- mvlist[[i]]$ct == cthighlightvec[i]
  }
  mvlist <- do.call(rbind, mvlist) %>% as.data.frame
  mvlist$ct_genes[mvlist$ct_genes=="(Intercept)"] <- "Sst Chodl"
  p <- ggplot(mvlist, aes(x=tval,y=ct, fill=highlight))+
    theme_gray() +
    geom_bar(stat="identity") +
    xlab("Expression (log10)") +
    ylab("") +
    geom_vline(xintercept=0,color="black") +
    ggtitle("Expression (log10) for\ncell-type marker genes of interest") +
    scale_fill_manual( values = c( "TRUE"="tomato", "FALSE"="gray" ), guide = "none" ) +
    theme(axis.title.y=element_text(size=24),
          plot.title=element_text(size=24,hjust=0.5),
          axis.title.x=element_text(size=24),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=14,hjust=1,vjust=1,angle=45),
          strip.text=element_text(size=24),
          panel.spacing = unit(1, "lines")) +
    coord_flip() +
    facet_wrap(~ct_genes, nrow=length(markervec)/2, ncol=2, scale="free") 

  ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/expression_markerGeneAverage_subclass.png"),height=32,width=14)

# Plot adjacency matrix
cormat <- cor(allmlist[[1]])
adjmat <- (cormat+1)/2

#########
# Calculate mean expression matrix over subclasses
######
donorvec <- c("H18.30.002", "H19.30.001", "H19.30.002")
markervec = c("MOG", "AIF1","GFAP", "GAD1", "SLC17A7")

# collect lists of cell ids per subclass for each donor
sublist_all <- lapply(donorvec, function(donor){
  cell_anno1 <- cell_anno %>% dplyr::filter(Donor==donor)
  return(tapply(cell_anno1[,3], cell_anno1[,1], list)) 
})

#### calculate mean expression over all genes for eac subclass
plotdf <- mapply(function(sublist, donor){
  summat <- lapply(sublist, function(x) rowMeans(cell_exprall[,colnames(cell_exprall) %in% x])) %>%
    do.call(cbind,.)
  out <- t(summat[rownames(summat) %in% markervec,]) %>% as.data.frame %>%
    rownames_to_column(var="ct") %>%
    pivot_longer(!ct, names_to="Gene") %>% 
    mutate(Gene=factor(Gene,levels=markervec), Donor=donor) %>% as.data.frame
  return(out)
}, sublist_all, donorvec, SIMPLIFY=F) %>% do.call(rbind,.) %>% group_by(ct,Gene) %>%
  summarise(mean_value=mean(value), se_value=sd(value)/sqrt(3))

p <- ggplot(plotdf,aes(x=ct,y=mean_value)) + 
  theme_bw() + 
  geom_bar(stat="identity") + 
  geom_errorbar(aes(x=ct,ymin=mean_value-2*se_value, ymax=mean_value+2*se_value), width=0.3, alpha=0.5) +
  theme(axis.title.y=element_text(size=24),
        plot.title=element_text(size=24,hjust=0.5),
        axis.title.x=element_text(size=24),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=14,hjust=1,vjust=1,angle=45),
        strip.text=element_text(size=24),
        panel.spacing = unit(1, "lines")) +
  labs(x="", y="Mean UMI") +
  facet_wrap(~Gene, nrow=6, scale="free")

ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/expression_of_marker_genes_of_interest.png"), height=12, width=10)

# scale each ct by mean expression in each ct
plotdf <- mapply(function(sublist, donor){
  summat <- lapply(sublist, function(x) rowMeans(cell_exprall[,colnames(cell_exprall) %in% x])) %>%
    do.call(cbind,.)
  summat2 <- apply(summat, 2, function(x) x/mean(x))
  summat3 <- t(summat2[rownames(summat2) %in% markervec,]) %>% as.data.frame %>%
    rownames_to_column(var="ct") %>%
    pivot_longer(!ct, names_to="Gene") %>% 
    mutate(Gene=factor(Gene,levels=markervec), Donor=donor) %>% as.data.frame
}, sublist_all, donorvec, SIMPLIFY=F) %>% do.call(rbind,.) %>% 
  group_by(ct,Gene) %>%
  summarise(mean_value=mean(value), se_value=sd(value)/sqrt(3))

p2 <- ggplot(plotdf,aes(x=ct,y=mean_value)) + 
  theme_bw() + 
  geom_bar(stat="identity") + 
  geom_errorbar(aes(x=ct,ymin=mean_value-2*se_value, ymax=mean_value+2*se_value), width=0.3, alpha=0.5) +
  theme(axis.title.y=element_text(size=24),
        plot.title=element_text(size=24,hjust=0.5),
        axis.title.x=element_text(size=24),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=14,hjust=1,vjust=1,angle=45),
        strip.text=element_text(size=24),
        panel.spacing = unit(1, "lines")) +
  labs(x="", y="Scaled mean UMI") +
  facet_wrap(~Gene, nrow=6, scale="free")

ggsave(p2, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/expression_of_marker_genes_of_interest_scaled.png"), height=12, width=10)

#### Calculate proportion of gene expression coming from each subclass

plotdf <- mapply(function(sublist, donor){
  summat <- lapply(sublist, function(x) rowSums(cell_exprall[,colnames(cell_exprall) %in% x])) %>%
    do.call(cbind,.)
  summat2 <- apply(summat, 1, function(x) x/sum(x))
  summat3 <- summat2[,colnames(summat2) %in% markervec] %>% as.data.frame %>%
    rownames_to_column(var="ct") %>%
    pivot_longer(!ct, names_to="Gene") %>% 
    mutate(Gene=factor(Gene,levels=markervec), Donor=donor) %>% as.data.frame
}, sublist_all, donorvec, SIMPLIFY=F) %>% do.call(rbind,.) %>%
  group_by(ct,Gene) %>%
  summarise(mean_value=mean(value), se_value=sd(value)/sqrt(3))

p <- ggplot(plotdf2,aes(x=ct,y=mean_value)) + 
  theme_bw() + 
  geom_bar(stat="identity") + 
  geom_errorbar(aes(x=ct,ymin=mean_value-2*se_value, ymax=mean_value+2*se_value), width=0.3, alpha=0.5) +
  theme(axis.title.y=element_text(size=24),
        plot.title=element_text(size=24,hjust=0.5),
        axis.title.x=element_text(size=24),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=14,hjust=1,vjust=1,angle=45),
        strip.text=element_text(size=24),
        panel.spacing = unit(1, "lines")) +
  labs(x="", y="Proportion of counts") +
  facet_wrap(~Gene, nrow=6, scale="free")

ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/fig1_followup/proportion_of_marker_genes_of_interest.png"), height=12, width=10)



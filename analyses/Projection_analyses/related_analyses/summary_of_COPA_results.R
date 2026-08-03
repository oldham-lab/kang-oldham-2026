# summarise COPA results
# - count number of sig mods per dataset

library(eulerr)

library(CoPA)
options(bitmapType = 'cairo')

# Paths to bulk megaset projections
paths_all <- list.dirs(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/"), full.names=T,recursive=F)[c(4,6,7)] 
paths_all <- c(paths_all, list.dirs(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/"), full.names=T,recursive=F)[c(4,5)])
paths_all <- c(paths_all, list.dirs(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_diffCoExpress/"), full.names=T,recursive=F))
dataset_names_bulk <- c("bulk_Mathys_2019", "bulk_Morabito", "bulkSEAAD2024", "rosmap_AD_SEAAD2024", "rosmap_con_SEAAD2024", "rosmap_diffCoExpress_conMinusAD", "rosmap_diffCoExpress_ADMinusCon")

# Create a matrix of counts for number of significantly differentially projected modules
produce_count_df <- function(COPA_output_paths,dataset_names){
  zlist <- list("all"=lapply(COPA_output_paths, function(x) qread(paste0(x, "/euclidean_distances/euclidean_sigmods_all.qs"))),
              "pos"=lapply(COPA_output_paths, function(x) qread(paste0(x, "/euclidean_distances/euclidean_sigmods_positive.qs"))),
              "neg"=lapply(COPA_output_paths, function(x) qread(paste0(x, "/euclidean_distances/euclidean_sigmods_negative.qs"))))
  zcounts <- lapply(zlist, function(direc){
    p <- lapply(direc, function(dataset){
      lapply(dataset, function(cutoff){
          unlist(lapply(cutoff, function(celltype) length(celltype))) %>% data.frame("counts"=.) %>% rownames_to_column("celltype")
      }) %>% map_df(~as.data.frame(.x), .id="cutoff")
    })
    names(p) <- dataset_names
    p <- p %>% map_df(~as.data.frame(.x), .id="dataset")
  }) %>% map_df(~as.data.frame(.x), .id="direction") 
  # Split cutoff column into celltype_scope and cutoff
  zcounts$celltype_scope <- ifelse(grepl("Subclass", zcounts$cutoff), "Subclass", "Supertype")
  zcounts$cutoff_type <- ifelse(grepl("FDR", zcounts$cutoff), "FDR", "nominal")
  zcounts$cutoff_type[grep("bonf",zcounts$cutoff)] <- "bonf"
  zcounts <- zcounts[,-3]
  return(zcounts)
}

bulkcounts <- produce_count_df(paths_all, dataset_names_bulk)

# Plot total counts per dataset
p <- bulkcounts %>% 
  mutate(direction=factor(direction, labels=list(all="All modules", pos="Higher in AD", neg="Higher in Con"))) %>%
  ggplot(aes(x=dataset, y=counts)) + 
    theme_bw() + 
    #geom_violin(aes(fill=celltype_scope)) +
    geom_boxplot(aes(fill=celltype_scope),notch=F) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=14),
          axis.text.y = element_text(size=14),
          axis.title.y=element_text(size=20, margin=margin(0,10,0,0)),
          strip.text=element_text(size=20),
          plot.title = element_text(size=24, hjust=0.5),
          plot.margin=margin(1,1,1,1,"cm")) + 
    facet_wrap(vars(direction,cutoff_type), scales="free_y") + 
    labs(title="Total number of significant modules per dataset", y="Number of significant modules", x="")
ggsave(p,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/copa_compare_summary.png"),width=16,height=12)

# bulk mods: 1158
# rosmapAD: 1127
# rosmapCON: 1186





## Plot venn diagrams of shared significant modules (pairwise)
COPA_output_paths <- paths_all[1:3]
dataset_names <- dataset_names_bulk[1:3]
zlist <- list("all"=lapply(COPA_output_paths, function(x) qread(paste0(x, "/euclidean_distances/euclidean_sigmods_all.qs"))),
              "pos"=lapply(COPA_output_paths, function(x) qread(paste0(x, "/euclidean_distances/euclidean_sigmods_positive.qs"))),
              "neg"=lapply(COPA_output_paths, function(x) qread(paste0(x, "/euclidean_distances/euclidean_sigmods_negative.qs"))))

# Fisher's exact test function
fisherTest <- function(group1,group2,all){
  y1y2 <- length(intersect(group1,group2))
  y1n2 <- length(group1)-y1y2
  n1y2 <- length(group2)-y1y2
  n1n2 <- length(all)-length(union(group1,group2))
  fisher.test(matrix(c(y1y2,n1y2,y1n2,n1n2),ncol=2),alternative="greater")$p.val
}
comparisons <- combn(length(zlist[[1]]),2)

# Calculate pairwise Fisher's exact test p-values
for(z in 1:3){ #all, pos, neg
  if(z==1){
    typetag <- "all"
  } else if (z==2){
    typetag <- "pos"
  } else {
    typetag <- "neg"
  }
  comp_outs <- zlist[[z]]
  fisherP_pairwise <- apply(comparisons,2,function(i){
    x <- comp_outs[[i[1]]]
    y <- comp_outs[[i[2]]]
    # match cts
    common <- unique(c(names(x), names(y)))
    x <- x[names(x) %in% common]
    y <- y[names(y) %in% common]
    return(lapply(seq_along(x), function(j){
        lapply(seq_along(x[[j]]), function(k){
          xsub <- x[[j]][[k]]
          ysub <- y[[j]][[k]]
          fisherTest(xsub, ysub, 1:1158)
       })         
    }))
  })

  # Count the # of overlapping modules (pairwise comparison)
  pairwise_overlap <- apply(comparisons,2,function(i){
    x <- comp_outs[[i[1]]]
    y <- comp_outs[[i[2]]]
    # match cts
    common <- unique(c(names(x), names(y)))
    x <- x[names(x) %in% common]
    y <- y[names(y) %in% common]
    return(lapply(seq_along(x), function(j){
      lapply(seq_along(x[[j]]), function(k) length(intersect(x[[j]][[k]], y[[j]][[k]])))
    }))
  })
  
  #if(z %in% c(2,3)){
    # plot euler diagrams (two datasets - Morabito vs SEA) (Morabito more deeply sequenced than Mathys 2019)
    for(g in c(2,3,5,6)){
      b <- comp_outs[[2]][[g]]
      c <- comp_outs[[3]][[g]]
      # match cts
      common <- unique(c(names(b), names(c)))
      b <- b[names(b) %in% common]
      c <- c[names(c) %in% common]
      
      eulermats <- list()
      for(j in 1:length(b)){
        commonmods <- unique(c(b[[j]],c[[j]]))
        eulermat <- data.frame("Morabito"=commonmods %in% b[[j]], 
                               "SEAAD2024"=commonmods %in% c[[j]])  
        eulermats[[j]] <- plot(euler(eulermat),
          main=paste0(names(b)[j], ",\np=",signif(fisherP_pairwise[[3]][[g]][[j]],2)),quantities=T,
            labels = list(font = 24))
      }

      keep_these <- which(unlist(fisherP_pairwise[[3]][[g]]) < 0.05)
      
      width1 <- 1200
      height1 <- 1000
      
      outdir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary_euler/")

      if(g==3){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEA_subclassFDR_"),typetag,".png"), width=width1, height=height1)
      } else if(g==6){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEAsupertypeFDR_"),typetag,".png"), width=width1, height=height1)
      } else if(g==2){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEA_subclassBon_"),typetag,".png"), width=width1, height=height1)
      } else if(g==5){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEAsupertypeBon_"),typetag,".png"), width=width1, height=height1)
      }
      gridExtra::grid.arrange(grobs=eulermats)
      dev.off()

      if(g==3){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEA_subclassFDR_"),typetag,"_sig.png"), width=width1, height=height1)
      } else if(g==6){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEAsupertypeFDR_"),typetag,"_sig.png"), width=width1, height=height1)
      } else if(g==2){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEA_subclassBon_"),typetag,"_sig.png"), width=width1, height=height1)
      } else if(g==5){
        png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_MorabitoVsSEAsupertypeBon_"),typetag,"_sig.png"), width=width1, height=height1)
      }
      gridExtra::grid.arrange(grobs=eulermats[keep_these], nrow=3)
      dev.off()
    }
  #}
  
  # plot euler diagrams (all three)
  for(g in c(3,6)){
    a <- comp_outs[[1]][[g]]
    b <- comp_outs[[2]][[g]]
    c <- comp_outs[[3]][[g]]
    # match cts
    common <- unique(c(names(a), names(b), names(c)))
    a <- a[names(a) %in% common]
    b <- b[names(b) %in% common]
    c <- c[names(c) %in% common]
    eulermats <- list()
    for(j in 1:length(a)){
      commonmods <- unique(c(a[[j]],b[[j]],c[[j]]))
      eulermat <- data.frame("Mathys"=commonmods %in% a[[j]], 
                             "Morabito"=commonmods %in% b[[j]], 
                             "SEAAD2024"=commonmods %in% c[[j]])  
      eulermats[[j]] <- plot(euler(eulermat), main=names(a)[j],quantities=T)
    }
    
    if(g==3){
      png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_3dataset_subclassFDR_"),typetag,".png"), width=width1, height=height1)
    } else if(g==6){
      png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_3dataset_supertypeFDR_"),typetag,".png"), width=width1, height=height1)
    } else if(g==2){
      png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_3dataset_subclassBon_"),typetag,".png"), width=width1, height=height1)
    } else if(g==5){
      png(paste0(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/copa_compare_summary/euler/copa_compare_summary_euler_3dataset_supertypeBon_"),typetag,".png"), width=width1, height=height1)
    }
    gridExtra::grid.arrange(grobs=eulermats)
    dev.off()
  }

}

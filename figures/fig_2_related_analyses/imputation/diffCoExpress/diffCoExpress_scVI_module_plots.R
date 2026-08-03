# Plot diffCoExpress module seed genes across both datasets 

library(tidyverse)
library(qs)
library(data.table)
library(cowplot)
library(RColorBrewer)
library(ggpubr)

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
megat <- t(mega[,-c(1:2)])
leinpbt <- t(leinpb[,-c(1:2)])
scvit <- t(scvi[,-c(1:2)])

# Convert megaset to positive values (same as FM)
meganonzero <- mega[,3:ncol(mega)] + abs(min(mega[,3:ncol(mega)], na.rm = TRUE)) + 1

# Load gene means and percent zeroes for Lein DFC (SN data)
zero_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/pctage_of_zeros_per_gene.qs"))
zero_mean <- zero_mean[names(zero_mean) %in% common_genes]
zero_mean <- zero_mean[match(common_genes, names(zero_mean))]
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means.qs"))
san_mean <- san_mean[names(san_mean) %in% common_genes]
san_mean <- san_mean[match(common_genes, names(san_mean))]

# Load adjacency matrices and convert to correlation matrices
simMega <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/megamat_adjMat.qs"))
simMega <- simMega*2-1
simPB <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/leinpb_donor1_adjMat.qs"))
simPB <- simPB*2-1
simSCVI <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/scvi_donor1_adjMat.qs"))
simSCVI <- simSCVI*2-1

# Function for plotting module snapshots over samples from all three datasets
plot_snapshots_allSamples <- function(modules_dir){
  # Load modules
  dirs <- list.dirs(modules_dir,full.names=T)[-1]
  for(i in dirs){
    cat("starting ", i, "\n")
    if(length(grep("kME_table", list.files(i)))==0){
      next
    } else {
      # - Megaset minus Lein pb
      mods <- fread(data.table=F, file=list.files(i,full.names=T)[grep("kME_table",list.files(i))])
      mods <- tapply(mods[,1], mods[,2], list)
      # take top 15 genes for display purposes
      mods <- lapply(mods, function(x){
          if(length(x)>15){
            return(x[1:15])
          } else {
            return(x)
          }
      })

      # log counts
      mergedats <- cbind(meganonzero, leinpb[,-c(1:2)], scvi[,-c(1:2)])
      colnames(mergedats) <- 1:ncol(mergedats)
      mergedatsmods <- lapply(mods, function(x){
          mergedats[common_genes %in% x,] %>% 
            mutate(Gene=x) %>%
            pivot_longer(!Gene, names_to="Sample", values_to="counts") %>%
            mutate("log_counts"=log2(counts+1),
                  "Sample"=factor(Sample, levels=unique(Sample)))
      })

      # z score
      megaz <- apply(mega[,-c(1,2)],2,scale)
      leinpbz <- apply(leinpb[,-c(1,2)],2,scale)
      scviz <- apply(scvi[,-c(1,2)],2,scale)
      mergedatsz <- cbind(megaz, leinpbz,scviz) %>% as.data.frame
      colnames(mergedatsz) <- 1:ncol(mergedatsz)
      mergedatsmodsz <- lapply(mods, function(x){
          mergedatsz[common_genes %in% x,] %>% 
            mutate(Gene=x) %>%
            pivot_longer(!Gene, names_to="Sample", values_to="z_counts") %>%
            mutate("Sample"=factor(Sample,levels=unique(Sample)))
      })

      # Determine dimensions of plot
      # find factors for n where n = the total number of modules
      factor_a_number <- function(x) {
        x <- as.integer(x)
        div <- seq_len(abs(x))
        factors <- div[x %% div == 0L]
        return(factors)
      }

      if(length(mods) < 40){
        factors <- factor_a_number(length(mods))
        if(length(factors)<=2){
          factors <- factor_a_number(length(mods)+1)
        }
        dim1 <- factors[which.min(abs(factors-sqrt(length(mods))))]
        dim2 <- ceiling(length(mods)/dim1)
      } else {
        num_figures <- ceiling(length(mods)/40)
        lengthvec <- rep(ceiling(length(mods)/num_figures), num_figures-1)
        lengthvec <- c(lengthvec, length(mods) - sum(lengthvec))
        factorlist <- lapply(lengthvec, function(x){
            out <- factor_a_number(x)
            if(length(out)<=2){
              out <- factor_a_number(x+1)
            }
            return(out)
          })
        dimlist <- lapply(factorlist, function(x){
          y <- x[which.min(abs(x-sqrt(max(x))))]
          y <- c(y,ceiling(max(x)/y))
        })
      }

      # Plot each module
      logplots <- mapply(function(x,b){
        ggplot(x, aes(x=Sample, y= log_counts,group=Gene, color=Gene)) +
          theme_bw() + 
          geom_line(linewidth=0.05) + 
          labs(x="Sample (left: Bulk, middle: PB, right: PB_CB_scVI)", 
              y=bquote(log[2]~"expression"), 
              title=b) +
          theme(text=element_text(size=10),
                axis.text.x=element_blank(),
                legend.position="bottom",
                plot.title=element_text(size=14, hjust=0.5)) +
          guides(color=guide_legend(title="Seed\ngenes",override.aes=list(linewidth=1)))
      },mergedatsmods,names(mergedatsmods),SIMPLIFY=F)

      lvs <- c(0, cumsum(lengthvec))
      if(length(mods)>=40){      
        for(z in seq_along(dimlist)){
          pout <- plot_grid(plotlist=logplots[(lvs[z]+1):lvs[z+1]], nrow=dimlist[[z]][[1]],ncol=dimlist[[z]][[2]],align='h')
          ggsave(pout,width=5*dimlist[[z]][[2]],height=3*dimlist[[z]][[1]], file=paste0(i,"/combined_datasets_seedGene_logCounts_exprPlot_", z,".png"),limitsize=F)
        }
      } else {
        pout <- plot_grid(plotlist=logplots, nrow=dim1,ncol=dim2,align='h')
        ggsave(pout,width=5*dim2,height=3*dim1, file=paste0(i,"/combined_datasets_seedGene_logCounts_exprPlot.png"),limitsize=F)
      }

      zplots <- mapply(function(x,b){
        ggplot(x, aes(x=Sample, y= z_counts,group=Gene, color=Gene)) +
          theme_bw() + 
          geom_line(linewidth=0.05) + 
          labs(x="Sample (left: Bulk, middle: PB, right: PB_CB_scVI)", 
              y="Z-score", 
              title=b) +
          theme(text=element_text(size=10),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                legend.position="bottom",
                plot.title=element_text(size=14, hjust=0.5)) +
          guides(color=guide_legend(title="Seed\ngenes",override.aes=list(linewidth=1)))
      },mergedatsmodsz,names(mergedatsmodsz),SIMPLIFY=F)
      
      if(length(mods)>=40){      
        for(z in seq_along(dimlist)){
          poutz <- plot_grid(plotlist=zplots[(lvs[z]+1):lvs[z+1]], nrow=dimlist[[z]][[1]],ncol=dimlist[[z]][[2]],align='h')
          ggsave(poutz,width=5*dimlist[[z]][[2]],height=3*dimlist[[z]][[1]], file=paste0(i,"/combined_datasets_seedGene_zscore_exprPlot_", z,".png"),limitsize=F)
        }
      } else {
        poutz <- plot_grid(plotlist=zplots, nrow=dim1,ncol=dim2,align='h')
        ggsave(poutz,width=5*dim2,height=3*dim1, file=paste0(i,"/combined_datasets_seedGene_zscore_exprPlot.png"),limitsize=F)
      }
    } # end of if
  } # end of for i in dirs
} # end of function

# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/")
# plot_snapshots_allSamples(modules_dir)
# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/")
# plot_snapshots_allSamples(modules_dir)

# Function for plotting module seed gene correlations over samples from all three datasets

plot_cors_allSamples <- function(modules_dir){
  # Load modules
  dirs <- list.dirs(modules_dir,full.names=T)[-1]
  for(i in dirs){
    cat("starting ", i, "\n")
    if(length(grep("kME_table", list.files(i)))==0){
      next
    } else {
      # - Megaset minus Lein pb
      mods <- fread(data.table=F, file=list.files(i,full.names=T)[grep("kME_table",list.files(i))])
      mods <- tapply(mods[,1], mods[,2], list)

      # Collect correlations
      cordfs <- lapply(mods, function(x){
        which <- common_genes %in% x
        df1 <- simMega[which,which]
        df1 <- df1[upper.tri(df1)]
        df2 <- simPB[which,which]
        df2 <- df2[upper.tri(df2)]
        df3 <- simSCVI[which,which]
        df3 <- df3[upper.tri(df3)]
        if(length(df1)>10000){
          set.seed(23)
          samp <- sample(1:length(df1), 10000)
          df1 <- df1[samp]
          df2 <- df2[samp]
          df3 <- df3[samp]
        }

        dfout <- data.frame("which"=c(rep("Bulk", length(df1)), rep("PB", length(df2)), rep("PB_CB_scVI", length(df3))),
                            "cors"=c(df1,df2,df3))
      })

      # Determine dimensions of plot
      # find factors for n where n = the total number of modules
      factor_a_number <- function(x) {
        x <- as.integer(x)
        div <- seq_len(abs(x))
        factors <- div[x %% div == 0L]
        return(factors)
      }

      if(length(mods) < 40){
        factors <- factor_a_number(length(mods))
        if(length(factors)<=2){
          factors <- factor_a_number(length(mods)+1)
        }
        dim1 <- factors[which.min(abs(factors-sqrt(length(mods))))]
        dim2 <- ceiling(length(mods)/dim1)
      } else {
        num_figures <- ceiling(length(mods)/40)
        lengthvec <- rep(ceiling(length(mods)/num_figures), num_figures-1)
        lengthvec <- c(lengthvec, length(mods) - sum(lengthvec))
        factorlist <- lapply(lengthvec, function(x){
            out <- factor_a_number(x)
            if(length(out)<=2){
              out <- factor_a_number(x+1)
            }
            return(out)
          })
        dimlist <- lapply(factorlist, function(x){
          y <- x[which.min(abs(x-sqrt(max(x))))]
          y <- c(y,ceiling(max(x)/y))
        })
        lvs <- c(0, cumsum(lengthvec))
      }

      # Plot each module
      brewer_colors <- brewer.pal(6,"Paired")[c(2,4,6)]
      corplots <- mapply(function(x,b){
        ggplot(x, aes(x=cors,group=which, color=which)) +
          theme_bw() + 
          geom_density() + 
          labs(x="Pearson correlation", 
              title=b) +
          theme(text=element_text(size=10),
                legend.position="bottom",
                plot.title=element_text(size=14, hjust=0.5)) +
          scale_color_manual(values=brewer_colors) +
          guides(color=guide_legend(title="Dataset",override.aes=list(linewidth=1))) 
      },cordfs,names(cordfs),SIMPLIFY=F)

      if(length(mods)>=40){      
        for(z in seq_along(dimlist)){
          pout <- plot_grid(plotlist=corplots[(lvs[z]+1):lvs[z+1]], nrow=dimlist[[z]][[1]],ncol=dimlist[[z]][[2]],align='h')
          ggsave(pout,width=4*dimlist[[z]][[2]],height=2*dimlist[[z]][[1]], file=paste0(i,"/combined_datasets_seedGene_corDensityPlot_", z,".png"),limitsize=F)
        }
      } else {
        pout <- plot_grid(plotlist=corplots, nrow=dim1,ncol=dim2,align='h')
        ggsave(pout,width=3*dim2,height=2*dim1, file=paste0(i,"/combined_datasets_seedGene_corDensityPlot.png"),limitsize=F)
      }
    } # end of if
  } # end of for i in dirs
} # end of function
# This function just recapitulates the overall correlation structure of the three datasets
# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/")
# plot_cors_allSamples(modules_dir)
# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/")
# plot_cors_allSamples(modules_dir)



# Function for plotting module seed gene %VE over samples from all three datasets

# plot_VE_allSamples <- function(modules_dir){
#   # Load modules
#   dirs <- list.dirs(modules_dir,full.names=T)[-1]
#   for(i in dirs){
#     cat("starting ", i, "\n")
#     if(length(grep("kME_table", list.files(i)))==0){
#       next
#     } else {
#       # - Megaset minus Lein pb
#       mods <- fread(data.table=F, file=list.files(i,full.names=T)[grep("kME_table",list.files(i))])
#       mods <- tapply(mods[,1], mods[,2], list)

#       # Collect %VE
#       vedfs <- lapply(mods, function(x){
#         which <- common_genes %in% x
#         df1 <- megat[,which]
#         eig1 <- prcomp(df1, scale=T, rank.=1)$x
#         ve1 <- apply(df1,2, function(y){
#           summary(lm(y~eig1))$r.squared
#         }) %>% unlist
        
#         df2 <- leinpbt[,which]
#         eig2 <- prcomp(df2, scale=T, rank.=1)$x
#         ve2 <- apply(df2,2, function(y){
#           summary(lm(y~eig2))$r.squared
#         }) %>% unlist
     
#         df3 <- scvit[,which]
#         eig3 <- prcomp(df3, scale=T, rank.=1)$x
#         ve3 <- apply(df3,2, function(y){
#           summary(lm(y~eig3))$r.squared
#         }) %>% unlist

#         dfout <- data.frame("which"=c(rep("Bulk", length(ve1)), rep("PB", length(ve2)), rep("PB_CB_scVI", length(ve3))),
#                             "VE"=c(ve1,ve2,ve3))
#       })

#       # Determine dimensions of plot
#       # find factors for n where n = the total number of modules
#       factor_a_number <- function(x) {
#         x <- as.integer(x)
#         div <- seq_len(abs(x))
#         factors <- div[x %% div == 0L]
#         return(factors)
#       }

#       if(length(mods) < 40){
#         factors <- factor_a_number(length(mods))
#         if(length(factors)<=2){
#           factors <- factor_a_number(length(mods)+1)
#         }
#         dim1 <- factors[which.min(abs(factors-sqrt(length(mods))))]
#         dim2 <- ceiling(length(mods)/dim1)
#       } else {
#         num_figures <- ceiling(length(mods)/40)
#         lengthvec <- rep(ceiling(length(mods)/num_figures), num_figures-1)
#         lengthvec <- c(lengthvec, length(mods) - sum(lengthvec))
#         factorlist <- lapply(lengthvec, function(x){
#             out <- factor_a_number(x)
#             if(length(out)<=2){
#               out <- factor_a_number(x+1)
#             }
#             return(out)
#           })
#         dimlist <- lapply(factorlist, function(x){
#           y <- x[which.min(abs(x-sqrt(max(x))))]
#           y <- c(y,ceiling(max(x)/y))
#         })
#       }

#       # Plot each module
#       brewer_colors <- brewer.pal(6,"Paired")[c(2,4,6)]
#       plots <- mapply(function(x,b){
#         ggplot(x, aes(x=VE,group=which, color=which)) +
#           theme_bw() + 
#           geom_density() + 
#           labs(x="% Variance explained by PC1", 
#               title=b) +
#           theme(text=element_text(size=10),
#                 legend.position="bottom",
#                 plot.title=element_text(size=14, hjust=0.5)) +
#           scale_color_manual(values=brewer_colors) +
#           guides(color=guide_legend(title="Dataset",override.aes=list(linewidth=1))) 
#       },vedfs,names(vedfs),SIMPLIFY=F)

#       lvs <- c(0, cumsum(lengthvec))
#       if(length(mods)>40){      
#         for(z in seq_along(dimlist)){
#           pout <- plot_grid(plotlist=plots[(lvs[z]+1):lvs[z+1]], nrow=dimlist[[z]][[1]],ncol=dimlist[[z]][[2]],align='h')
#           ggsave(pout,width=4*dimlist[[z]][[2]],height=2*dimlist[[z]][[1]], file=paste0(i,"/combined_datasets_seedGene_VEPlot_", z,".png"),limitsize=F)
#         }
#       } else {
#         pout <- plot_grid(plotlist=plots, nrow=dim1,ncol=dim2,align='h')
#         ggsave(pout,width=3*dim2,height=2*dim1, file=paste0(i,"/combined_datasets_seedGene_VEPlot.png"),limitsize=F)
#       }
#     } # end of if
#   } # end of for i in dirs
# } # end of function

# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/")
# plot_VE_allSamples(modules_dir)
# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/")
# plot_VE_allSamples(modules_dir)


# Plot box/violin plots for VE

plot_VE_boxviolin_allSamples <- function(modules_dir){
  # Load modules
  dirs <- list.dirs(modules_dir,full.names=T)[-1]
  for(i in dirs){
    cat("starting ", i, "\n")
    if(length(grep("kME_table", list.files(i)))==0){
      next
    } else {
      # - Megaset minus Lein pb
      mods <- fread(data.table=F, file=list.files(i,full.names=T)[grep("kME_table",list.files(i))])
      mods <- tapply(mods[,1], mods[,2], list)

      # Calculate Lein DFC SN means and percent zeroes
      modsan <- lapply(mods, function(x) mean(san_mean[names(san_mean) %in% x])) %>% unlist %>% sort
      modzeroes <- lapply(mods, function(x) mean(zero_mean[names(zero_mean) %in% x])) %>% unlist %>% sort

      # Collect %VE
      vedfs <- lapply(mods, function(x){
        which <- common_genes %in% x
        df1 <- megat[,which]
        eig1 <- prcomp(df1, scale=T, rank.=1)$x
        ve1 <- apply(df1,2, function(y){
          summary(lm(y~eig1))$r.squared
        }) %>% unlist
        
        df2 <- leinpbt[,which]
        eig2 <- prcomp(df2, scale=T, rank.=1)$x
        ve2 <- apply(df2,2, function(y){
          summary(lm(y~eig2))$r.squared
        }) %>% unlist
     
        df3 <- scvit[,which]
        eig3 <- prcomp(df3, scale=T, rank.=1)$x
        ve3 <- apply(df3,2, function(y){
          summary(lm(y~eig3))$r.squared
        }) %>% unlist

        dfout <- data.frame("which"=c(rep("Bulk", length(ve1)), rep("PB", length(ve2)), rep("PB_CB_scVI", length(ve3))),
                            "VE"=c(ve1,ve2,ve3))
      })

      # Determine dimensions of plot
      # find factors for n where n = the total number of modules
      factor_a_number <- function(x) {
        x <- as.integer(x)
        div <- seq_len(abs(x))
        factors <- div[x %% div == 0L]
        return(factors)
      }

      if(length(mods) < 40){
        factors <- factor_a_number(length(mods))
        if(length(factors)<=2){
          factors <- factor_a_number(length(mods)+1)
        }
        dim1 <- factors[which.min(abs(factors-sqrt(length(mods))))]
        dim2 <- ceiling(length(mods)/dim1)
      } else {
        num_figures <- ceiling(length(mods)/40)
        lengthvec <- rep(ceiling(length(mods)/num_figures), num_figures-1)
        lengthvec <- c(lengthvec, length(mods) - sum(lengthvec))
        factorlist <- lapply(lengthvec, function(x){
            out <- factor_a_number(x)
            if(length(out)<=2){
              out <- factor_a_number(x+1)
            }
            return(out)
          })
        dimlist <- lapply(factorlist, function(x){
          y <- x[which.min(abs(x-sqrt(max(x))))]
          y <- c(y,ceiling(max(x)/y))
        })
        lvs <- c(0, cumsum(lengthvec))
      }

      # Order graphs by mean pcnt of zeroes of genes in modules
      vedfs <- vedfs[match(names(modzeroes), names(vedfs))]
      # Plot each module
      comps <- list(c("Bulk", "PB"), c("Bulk", "PB_CB_scVI"), c("PB", "PB_CB_scVI"))
      brewer_colors <- brewer.pal(6,"Paired")[c(2,4,6)]
      plots <- mapply(function(x,b){
        ggplot(x, aes(x=VE,y=which)) +
          theme_bw() + 
          geom_violin(aes(fill=which)) +
          geom_boxplot(width=0.2, fill="white", notch=T) +  
          geom_jitter(alpha=0.2, size=0.1) +
          labs(x="% Variance explained by PC1", 
              title=b, y="") +
          theme(text=element_text(size=10),
                legend.position="none",
                plot.title=element_text(size=14, hjust=0.5)) +
          scale_fill_manual(values=brewer_colors) +
          guides(fill=guide_legend(title="Dataset")) +
          stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=2) 
      },vedfs,names(vedfs),SIMPLIFY=F)

      if(length(mods)>40){      
        for(z in seq_along(dimlist)){
          pout <- plot_grid(plotlist=plots[(lvs[z]+1):lvs[z+1]], nrow=dimlist[[z]][[1]],ncol=dimlist[[z]][[2]],align='h')
          ggsave(pout,width=4*dimlist[[z]][[2]],height=2*dimlist[[z]][[1]], file=paste0(i,"/combined_datasets_seedGene_VEPlot_boxviolin_", z,".png"),limitsize=F)
        }
      } else {
        pout <- plot_grid(plotlist=plots, nrow=dim1,ncol=dim2,align='h')
        ggsave(pout,width=3*dim2,height=2*dim1, file=paste0(i,"/combined_datasets_seedGene_VEPlot_boxviolin.png"),limitsize=F)
      }
    } # end of if
  } # end of for i in dirs
} # end of function

# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/")
# plot_VE_boxviolin_allSamples(modules_dir)
# modules_dir <- file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/")
# plot_VE_boxviolin_allSamples(modules_dir)


# Run fxns on all networks
modules_dirs <- list.dirs(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/"),recursive=F)
lapply(modules_dirs, plot_snapshots_allSamples)
lapply(modules_dirs, plot_cors_allSamples)
lapply(modules_dirs, plot_VE_boxviolin_allSamples)


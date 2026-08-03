source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "utils/ggplot_theme_settings.R"))
library(cowplot)
library(tidyverse)
library(data.table)
library(qs)
library(ggrastr)
library(ggpubr)
options(bitmapType = 'cairo')
library(showtext)
showtext_auto()

# --- pseudobulk utilities (copied from Analyses/SN_pseudobulking/0-sn_pseudobulk_fxns.R) ---

# divide each sample by library size, multiply by 1 million
# like tpm without the gene length correction; for pseudobulk count data
norm_pb_samples <- function(expr,index){
  expr_sum <- apply(expr[,index:ncol(expr)],2,sum)
  expr[,index:ncol(expr)] <- expr[,index:ncol(expr)]/(expr_sum/1e6)
  return(expr)
}

# process pseudobulk data created from lein data
# - subset to genes in bulk megaset
# - remove zero var
# - save list of expr and transposed expr
# - save simMat
process_lein_pb <- function(expr_lp_path, #path
                            orig_genes,
                            filter=NULL,
                            filtername=NULL,
                            save_simMat=T,
                            genemap=NULL){
  csvname <- strsplit(expr_lp_path,"/")[[1]][length(strsplit(expr_lp_path,"/")[[1]])]
  topdir <- gsub(csvname, "", expr_lp_path)
  if(is.null(genemap)){
    genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"), data.table=F)
  }
  expr_lp <- fread(expr_lp_path, data.table=F)
  rownames(expr_lp) <- orig_genes
  expr_lp <- expr_lp[rownames(expr_lp) %in% genemap[,2],]
  genemap <- genemap[genemap[,2] %in% rownames(expr_lp),]
  genemap <- genemap[match(rownames(expr_lp), genemap[,2]),]
  # Add gene data to expression matrix
  expr_lp <- cbind(genemap, expr_lp[,-1])
  colnames(expr_lp)[1:2] <- c("ensembl_id", "Gene")
  # remove NA
  has_na <- apply(expr_lp,1,function(x) sum(is.na(x)))
  expr_lp <- expr_lp[has_na==0,]
  # Subset to genes in bulk dataset
  if(!is.null(filter)){
    expr_lp <- expr_lp[expr_lp[,1] %in% filter,]
    filename <- paste0(topdir, "/", gsub(".csv", "", csvname), "_EXPRLIST_PROCESSED_ZEROVAR_",filtername, ".RDS")
    simmatname <- paste0(topdir, "/", gsub(".csv", "", csvname), "_SIMMAT_PROCESSED_ZEROVAR_",filtername, ".RDS")
  } else {
    filename <- paste0(topdir, "/", gsub(".csv", "", csvname), "_EXPRLIST_PROCESSED_ZEROVAR.RDS")
    simmatname <- paste0(topdir, "/", gsub(".csv", "", csvname), "_SIMMAT_PROCESSED_ZEROVAR.RDS")
  }
  # Transpose and Remove genes with zero variance
  expr_t <- t(expr_lp[,3:ncol(expr_lp)])
  zvar <- apply(expr_t,2,var)
  expr_lp <- expr_lp[!zvar==0,]
  expr_t <- expr_t[,!zvar==0]
  colnames(expr_t) <- expr_lp[,1]
  expr_list <- list("expr"=expr_lp, "expr_t"=expr_t)
  saveRDS(expr_list, file = filename)
  # Create simMat
  if(save_simMat){
    simMat <- cor(expr_t)
    colnames(simMat) <- expr_lp$ensembl_id
    rownames(simMat) <- expr_lp$ensembl_id
    saveRDS(simMat, file = simmatname)
  }
}

# create matrix (ct x sample) of celltype counts per pseudobulk sample
calc_CT_counts <- function(sif,
                           cell_anno,
                           ct_col,
                           cell_id_col){
  sifct <- cell_anno[match(sif[,1], cell_anno[,cell_id_col]), ct_col]
  sifct <- factor(sifct, levels=unique(sifct))
  sampct <- list()
  whichcols <- grep("Sample", colnames(sif))
  for(i in whichcols){
    sampct[[i-1]] <- c(table(sifct[sif[,i]==1]))
  }
  sampct <- do.call(cbind, sampct)
  return(sampct)
}


create_donorspecific_pseudobulk <- function(cell_exprall,
                                            cell_annoall,
                                            donorvec, # vector of all unique donors
                                            save_dirall,
                                            pcnt.varvec=c(0), # vector of values for pcnt.var
                                            no.samples,
                                            cell.nameindex, # column index of cell_annoall that indicates cell name (matches colnames of cell_exprall)
                                            filter # vector of gene symbols to filter pseudobulk matrix 
                                            ){

  orig_genes <- rownames(cell_exprall)

  for(i in seq_along(donorvec)){ # for each donor
    for(j in pcnt.varvec){
        save_dir <- paste0(save_dirall,"/donor",i)
        if(!dir.exists(save_dir)){dir.create(save_dir, recursive=T)}
        cell_expr <- cell_exprall[,cell_annoall$Donor==donorvec[i]]
        cell_anno <- cell_annoall[cell_annoall$Donor==donorvec[i],]

        pcnt.var <- j
        
        setwd(save_dir)
        makeSyntheticDatasets(
            cell_expr,
            sampleindex=c(1:ncol(cell_expr)),
            cell.info=cell_anno,
            cell.name=cell.nameindex, # new labs
            cell.type=NULL,
            cell.frac=NULL,
            pcnt.cells=10,
            pcnt.var=pcnt.var,
            no.samples=no.samples, # size of megaset
            no.datasets=1
        )
        
        paths <- list.files(paste0(save_dir,"/SyntheticDatasets/"),full.names=T)
        paths <- paths[grep(paste0(pcnt.var,"pcntVar"),paths)]
        paths <- paths[grep(paste0(no.samples,"samples"),paths)]
        
        expr_lp_path <- paths[-grep("legend",paths)]
        expr_lp_path <- expr_lp_path[grep(".csv",expr_lp_path)]
        sifpath <- paths[grep("legend",paths)]
        
        process_lein_pb(
          expr_lp_path, 
          orig_genes,
          filter = megaexpr[,1],
          filtername = "BULKGENESUBSET",
          save_simMat=F
        )
        cat(j, " ")
    }
    cat(i, "i done\n")
  }
}

model_donorspecific_pseudobulk <- function(dirlist,
                                           cell_annoall,
                                           cellannocolvec,
                                           cellIDcol,
                                           celldonorcol,
                                           celltypenames = c("subclass", "supertype")){

    for(save_dirall in dirlist){
        plotdflist <- list()
        for(cellannocol in seq_along(cellannocolvec)){ 
          donors <- unique(cell_annoall[ ,celldonorcol])
          for(i in seq_along(donors)){
            exdir <- list.files(paste0(save_dirall,"/donor",i,"/SyntheticDatasets/"),full.names=T)
            exdirexpr <- exdir[grep("EXPRLIST",exdir)]
            exdirexpr <- exdirexpr[grep("0pcntVar", exdirexpr)]
            exdirsif <- exdir[grep("samples_legend",exdir)]
            exdirsif <- exdirsif[grep("0pcntVar", exdirsif)]
            #exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[9]))
            #exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
            exdirnames <- 0

            exdirlist <- mapply(function(x,y,z){
                obj <- readRDS(x) 
                obj <- norm_pb_samples(obj[[1]], index=3)
                exprt <- t(obj[,3:ncol(obj)])
                escale <- apply(exprt, 2, scale)
                rownames(escale) <- rownames(exprt)
                escaleperm <- apply(escale, 2,function(x) sample(x, length(x))) # permuted matrix
                genemeansfilt <- apply(exprt,2,function(x) log2(mean(x)))
                
                sif <- fread(data.table=F,file=y)
                sampct <- t(calc_CT_counts(sif, cell_annoall, ct_col=cellannocolvec[cellannocol], cell_id_col = cellIDcol))
                sampctmean <- apply(sampct, 2, mean)
                sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
                
                # Run modeling
                mod_list <- foreach(l=1:ncol(escale)) %dopar% {
                  lm(escale[,l] ~ ., data=as.data.frame(sampct))
                } 
        
                mod_r2a <- unlist(lapply(mod_list, function(x) summary(x)$adj.r.squared))
                mod_p <- unlist(lapply(mod_list, function(x){
                    f_statistic <- summary(x)$fstatistic
                    overall_p_value <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
                    return(overall_p_value)
                }))
                mod_rmse <- unlist(lapply(mod_list, function(x) sqrt(mean(x$residuals^2))))
                rm(mod_list)
                gc()
                        
                plotdf <- data.frame("adj_r2" = mod_r2a,"rmse"=mod_rmse, "mean" = genemeansfilt,"type"="Real",
                                     "pcnt.var"=z, "Gene"=obj[,2], "pval" = mod_p)
                return(plotdf)
              }, exdirexpr,exdirsif,exdirnames, SIMPLIFY=F)
            plotdflist[[i]] <- do.call(rbind,exdirlist) %>% remove_rownames
            } # for i
            qsave(plotdflist,file=paste0(save_dirall,"/fig1bdflist_",celltypenames[cellannocol],".qs"))
        } # for cellannocol
        cat(save_dirall, " done\n\n\n")
    } # for save_dirall
}


plot_fig2_plots_overall <- function(homedir = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/"),
                            cell_anno = fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE),
                            genemap = fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F),
                            outdir,
                            donornames,
                            san_mean, # list of three (for each donor) of vectors of mean of counts for each gene 
                            panelc_ymax_r2=1.1,
                            panelc_ymax_rmse=1.2){
     
    #############################
    # Load modeling data and plot
    #############################
    if(!dir.exists(outdir)){
        dir.create(outdir, recursive=T)
    }

    plotdflist_sub <- qread(paste0(homedir,"/fig1bdflist_subclass.qs"))
    plotdflist_super <- qread(paste0(homedir,"/fig1bdflist_supertype.qs"))

    # Calculate mean r2 for subclasses
    meansvec_sub <- lapply(plotdflist_sub, function(x){
    df <- x %>% group_by(pcnt.var) %>%
        summarise(mean=mean(adj_r2)) %>%
        dplyr::filter(pcnt.var==0) %>% 
        mutate(mean_label=paste0("Mean adj. r2 =\n",signif(mean,2)),
            pcnt.var.label=paste0("Pcnt. var = ", pcnt.var)) %>%
        mutate(pcnt.var=as.factor(pcnt.var.label))
    })

    # Calculate mean r2 for supertypes
    meansvec_super <- lapply(plotdflist_super, function(x){
    df <- x %>% group_by(pcnt.var) %>%
        summarise(mean=mean(adj_r2)) %>%
        dplyr::filter(pcnt.var==0) %>%
        mutate(mean_label=paste0("Mean adj. r2 =\n",signif(mean,2)),
            pcnt.var.label=paste0("Pcnt. var = ", pcnt.var)) %>%
        mutate(pcnt.var=as.factor(pcnt.var.label))
    })

    rmsevec_sub <- lapply(plotdflist_sub, function(x){
        df <- x %>% group_by(pcnt.var) %>%
        summarise(mean=mean(rmse)) %>%
        dplyr::filter(pcnt.var==0) %>% 
        mutate(mean_label=paste0("Mean RMSE =\n",signif(mean,2)),
           pcnt.var.label=paste0("Pcnt. var = ", pcnt.var)) %>%
        mutate(pcnt.var=as.factor(pcnt.var.label))
    })

    rmsevec_super <- lapply(plotdflist_super, function(x){
        df <- x %>% group_by(pcnt.var) %>%
        summarise(mean=mean(rmse)) %>%
        dplyr::filter(pcnt.var==0) %>%
        mutate(mean_label=paste0("Mean RMSE =\n",signif(mean,2)),
           pcnt.var.label=paste0("Pcnt. var = ", pcnt.var)) %>%
        mutate(pcnt.var=as.factor(pcnt.var.label))
    })

    # # Load DE genes from Lein et al
    # origdemat <- fread(data.table=F, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein_orig_de_genes.csv"))
    # origdelist <- tapply(origdemat[,1], origdemat[,2], list)
    # allorigdegenes <- unique(unlist(origdelist))

    # Collect gene symbols from each donor (pcnt.var=0)
    genesymlist <- list()
    for(i in seq_along(donornames)){
        exdir <- list.files(paste0(homedir,"/donor",i,"/SyntheticDatasets/"),full.names=T)
        exdirexpr <- exdir[grep("EXPRLIST",exdir)]
        exdirexpr <- exdirexpr[grep("0pcntVar", exdirexpr)]
        obj <- readRDS(exdirexpr) 
        genesymlist[[i]] <- obj[[1]][,2]
    }

    # Add gene symbols to plotdflist objects
    plotdflist_sub <- mapply(function(x,y){
      x %>% dplyr::filter(pcnt.var==0) %>%
        mutate(Gene=y)
    }, plotdflist_sub, genesymlist, SIMPLIFY=F)

    plotdflist_super <- mapply(function(x,y){
      x %>% dplyr::filter(pcnt.var==0) %>%
        mutate(Gene=y)
    }, plotdflist_super, genesymlist, SIMPLIFY=F)
 
    # Replace mean column with gene means calculated from SN counts
    san_mean <- lapply(san_mean, function(x){
        as.data.frame(x) %>% rename(mean="x") %>% rownames_to_column(var="Gene") %>% return
    })
    plotdflist_sub <- mapply(function(a,b){
      a %>% dplyr::filter(pcnt.var==0) %>%
        dplyr::select(!mean) %>%
        left_join(b) %>% return
    }, plotdflist_sub, san_mean, SIMPLIFY=F)
    plotdflist_super <- mapply(function(a,b){
      a %>% dplyr::filter(pcnt.var==0) %>%
        dplyr::select(!mean) %>%
        left_join(b) %>% return
    }, plotdflist_super, san_mean, SIMPLIFY=F)
    
    # Plot first two columns of figure 1
    plot_cols_facetwrap <- function(a,b,c, save_dir, filename){
        meansvec <- do.call(rbind, b) %>%
            mutate(donor = donornames, 
                   mean_label = gsub("adj. r2 ", "", mean_label))
        out <- mapply(function(i, plotdflist){
          plotdflist[plotdflist$pcnt.var==0,] %>%
            mutate(donor=donornames[i])#,
                   #de_gene=Gene %in% allorigdegenes,
                   #sancut = ifelse(mean>=1, 
                    #               paste0(">=1 mean\nUMI/cell\n(n=", sum(mean>=1), ")"),
                     #              paste0("<1 mean\nUMI/cell\n(n=", sum(mean<1), ")")))
        }, seq_along(donornames), a, SIMPLIFY=F) %>% do.call(rbind,.) 
        genemeansall <- out %>% group_by(Gene) %>% summarise(gene_mean=mean(mean))
        out <- out %>%
          mutate(sancut = ifelse(mean>=1, 
                                 paste0(">=1 mean\nUMI/cell\n(n=", format(sum(genemeansall$gene_mean>=1), big.mark=","), ")"),
                                 paste0("<1 mean\nUMI/cell\n(n=", format(sum(genemeansall$gene_mean<1), big.mark=","), ")")),
                 mean=log2(mean))
        #out$sancut[out$Gene %in% genemeansall$Gene[genemeansall$gene_mean>=sancut]] <- paste0(">=1 mean\nUMI/cell\n(n=", sum(genemeansall$gene_mean>=sancut), ")")
        #out$sancut[out$Gene %in% genemeansall$Gene[genemeansall$gene_mean<sancut]] <- paste0("<1 mean\nUMI/cell\n(n=", sum(genemeansall$gene_mean<sancut), ")")
 
        textpos <- mean(c(min(out$mean), 0))

        out1 <- out %>%
            ggplot(aes(y=adj_r2, x=mean)) +
            theme_bw() +
            rasterise(geom_point(alpha=0.05,size=0.5), dpi = 300) +
            stat_density_2d(geom = "polygon", contour = TRUE,
                            aes(fill = after_stat(level)), colour = "black",
                            bins = 5) +
            scale_fill_distiller(palette = "Blues", direction = 1) +
            xlab(bquote("Mean expression ("~log[2]~")")) +
            ylab(bquote("Adjusted R"^2)) +
            geom_hline(data=meansvec, aes(yintercept=mean),linetype="dashed",color="red", linewidth=2, alpha=0.5) +
            ylim(0,1) +
            geom_vline(xintercept=log2(1),linetype="dashed",color="blue",linewidth=2, alpha=0.5) +
            #geom_vline(data=genemeansfiltvec[[i]],aes(xintercept=mean),linetype="dashed",color="blue",linewidth=2, alpha=0.5) +
            geom_text(data=meansvec,aes(label=mean_label, y=mean+0.25), x = textpos, color="red", size=10, alpha=0.8) +
            theme(text=element_text(size=26),
                    legend.position="none",
                    axis.title.x=element_text(size=26,margin=margin(10,0,0,0)),
                    axis.title.y=element_text(size=26,margin=margin(0,10,0,0)),
                    axis.text.x=element_text(size=26),
                    axis.text.y=element_text(size=26),
                    plot.title=element_text(hjust=0.5),
                    strip.text=element_text(size=26),
                    plot.margin=margin(1,1,1,1,"cm")) +
            facet_wrap(~donor, nrow=length(donornames))
        #ggsave(out1, file = file.path(save_dir, paste0("fig2_part1_", filename, ".png")))
        #ggsave(out1, file = file.path(save_dir, paste0("fig2_part1_", filename, ".pdf")))
        

        comps <- list(unique(out$sancut))
        out2 <- out %>%
          ggplot(aes(x=sancut, y=adj_r2, fill=sancut)) +
            theme_bw() +
            geom_violin() + 
            ylim(NA,panelc_ymax_r2) +
            labs(x="", y=bquote("Adjusted R"^2)) +
            theme(text=element_text(size=26),
                  legend.position="none",
                  axis.title.x=element_text(size=26,margin=margin(10,0,0,0)),
                  axis.title.y=element_text(margin=margin(0,10,0,0)),
                  axis.text.x=element_text(size=26),
                  plot.title=element_text(hjust=0.5),
                  plot.margin=margin(1,1,1,1,"cm")) +
            stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=12) +
            #stat_compare_means(method.args = list(alternative = "greater"), label.y=0.85, size=12)# +
            scale_fill_manual(values=c("white", "blue")) 
       # ggsave(out2, file = file.path(save_dir, paste0("fig2_part1_", filename, "_san.png")))
       # ggsave(out2, file = file.path(save_dir, paste0("fig2_part1_", filename, "_san.pdf")))
        return(list(out1, out2))
    }

    plot_cols_facetwrap_rmse <- function(a,b,c,save_dir,filename){
        meansvec <- do.call(rbind, b) %>%
            mutate(donor=donornames, 
                   mean_label = gsub("RMSE ", "", mean_label))
        out <- mapply(function(i, plotdflist){
          plotdflist[plotdflist$pcnt.var==0,] %>%
            mutate(donor=donornames[i])#,
                   #de_gene=Gene %in% allorigdegenes,
                   #sancut = ifelse(mean>=1, 
                    #               paste0(">=1 mean\nUMI/cell\n(n=", sum(mean>=1), ")"),
                     #              paste0("<1 mean\nUMI/cell\n(n=", sum(mean<1), ")")))
        }, seq_along(donornames), a, SIMPLIFY=F) %>% do.call(rbind,.) 
        genemeansall <- out %>% group_by(Gene) %>% summarise(gene_mean=mean(mean))
        out <- out %>%
          mutate(sancut = ifelse(mean>=1, 
                                 paste0(">=1 mean\nUMI/cell\n(n=", format(sum(genemeansall$gene_mean>=1), big.mark=","), ")"),
                                 paste0("<1 mean\nUMI/cell\n(n=", format(sum(genemeansall$gene_mean<1), big.mark=","), ")")),
                 mean=log2(mean))
        #out$sancut[out$Gene %in% genemeansall$Gene[genemeansall$gene_mean>=sancut]] <- paste0(">=1 mean\nUMI/cell\n(n=", sum(genemeansall$gene_mean>=sancut), ")")
        #out$sancut[out$Gene %in% genemeansall$Gene[genemeansall$gene_mean<sancut]] <- paste0("<1 mean\nUMI/cell\n(n=", sum(genemeansall$gene_mean<sancut), ")")

        textpos <- mean(c(min(out$mean), 0))

        out1 <- out %>%
            ggplot(aes(y=rmse, x=mean)) +
            theme_bw() +
            rasterise(geom_point(alpha=0.05,size=0.5)) +
            stat_density_2d(geom = "polygon", contour = TRUE,
                            aes(fill = after_stat(level)), colour = "black",
                            bins = 5) +
            scale_fill_distiller(palette = "Blues", direction = 1) +
            xlab(bquote("Mean expression ("~log[2]~")")) +
            ylab("RMSE") +
            geom_hline(data=meansvec, aes(yintercept=mean),linetype="dashed",color="red", linewidth=2, alpha=0.5) +
            ylim(0,1) +
            geom_vline(xintercept=log2(1),linetype="dashed",color="blue",linewidth=2, alpha=0.5) +
            #geom_vline(data=genemeansfiltvec[[i]],aes(xintercept=mean),linetype="dashed",color="blue",linewidth=2, alpha=0.5) +
            geom_text(data=meansvec,aes(label=mean_label, y=mean-0.25), x = textpos, color="red", size=10, alpha=0.8) +
            theme(text=element_text(size=26),
                    legend.position="none",
                    axis.title.x=element_text(size=26,margin=margin(10,0,0,0)),
                    axis.title.y=element_text(size=26,margin=margin(0,10,0,0)),
                    plot.title=element_text(hjust=0.5),
                    axis.text.x=element_text(size=26),
                    axis.text.y=element_text(size=26),
                    strip.text=element_text(size=26),
                    plot.margin=margin(1,1,1,1,"cm")) +
            facet_wrap(~donor, nrow=length(donornames))
        #ggsave(out1, file = file.path(save_dir, paste0("fig2_part1_", filename, ".png")))
        #ggsave(out1, file = file.path(save_dir, paste0("fig2_part1_", filename, ".pdf")))
        
        comps <- list(unique(out$sancut))
        out2 <- out %>%
          ggplot(aes(x=sancut, y=rmse, fill=sancut)) +
            theme_bw() +
            geom_violin() + 
            labs(x="", y="RMSE") +
            ylim(NA, panelc_ymax_rmse) +
            theme(text=element_text(size=26),
                  legend.position="none",
                  axis.text.x=element_text(size=26),
                  axis.title.x=element_text(size=26,margin=margin(10,0,0,0)),
                  axis.title.y=element_text(margin=margin(0,10,0,0)),
                  plot.title=element_text(hjust=0.5),
                  plot.margin=margin(1,1,1,1,"cm")) +
            scale_fill_manual(values=c("white", "blue")) +
            stat_compare_means(comparisons=comps,method="wilcox.test",label = "p.signif", size=12) 
            #stat_compare_means(label.y=.85)
       # ggsave(out2, file = file.path(save_dir, paste0("fig2_part1_", filename, "_san.png")))
       # ggsave(out2, file = file.path(save_dir, paste0("fig2_part1_", filename, "_san.pdf")))

       return(list(out1, out2))       
    }

    p1 <- plot_cols_facetwrap(plotdflist_sub, meansvec_sub,"Subclass (n=24)", outdir, "sub")
    p2 <- plot_cols_facetwrap(plotdflist_super, meansvec_super, "Supertype (n=136)", outdir, "super")
    heightvar = 1000 * length(donornames) / 3
    heightvar2 = 12 * length(donornames) / 3

    fig1p1 <- plot_grid(p1[[1]], 
                    ncol=1,
                    #labels=c("B","C"),
                    label_size=42)

    png(file=paste0(outdir,"/part1_sub.png"),width=500,height=heightvar)
    plot(fig1p1)
    dev.off()
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub.pdf"), width = 6, height = heightvar2)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub.svg"), width = 6, height = heightvar2)

    fig1p1 <- plot_grid(p2[[1]],
                    ncol=1,
                    #labels=c("B","C"),
                    label_size=42)
    png(file=paste0(outdir,"/part1_super.png"),width=500,height=heightvar)
    plot(fig1p1)
    dev.off()
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super.pdf"), width = 6, height = heightvar2)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super.svg"), width = 6, height = heightvar2)

    fig1p1 <- plot_grid(p1[[2]],
                    ncol=1,
                    #labels=c("B","C"),
                    label_size=42)
    png(file=paste0(outdir,"/part1_sub_san.png"),width=500,height=400)
    plot(fig1p1)
    dev.off()
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub_san.pdf"), width = 6, height = 5)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub_san.svg"), width = 6, height = 5)

    fig1p1 <- plot_grid(p2[[2]],
                    ncol=1,
                    #labels=c("B","C"),
                    label_size=42)
    png(file=paste0(outdir,"/part1_super_san.png"),width=500,height=400)
    plot(fig1p1)
    dev.off()
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super_san.pdf"), width = 6, height = 5)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super_san.svg"), width = 6, height = 5)

    # p1 <- plot_cols_facetwrap_deGeneColor(plotdflist_sub, meansvec_sub,"Subclass (n=24)")
    # p2 <- plot_cols_facetwrap_deGeneColor(plotdflist_super, meansvec_super, "Supertype (n=136)")

    # fig1p1 <- plot_grid(p1, p2,
    #                 ncol=2,
    #                 #labels=c("B","C"),
    #                 label_size=42)
    # png(file=paste0(outdir,"/fig1_part1_colorByDE.png"),width=1000,height=heightvar)
    # plot(fig1p1)
    # dev.off()

    p1 <- plot_cols_facetwrap_rmse(plotdflist_sub, rmsevec_sub,"Subclass (n=24)", outdir, "sub_rmse")
    p2 <- plot_cols_facetwrap_rmse(plotdflist_super, rmsevec_super, "Supertype (n=136)", outdir, "super_rmse")

    fig1p1 <- plot_grid(p1[[1]], 
                  ncol=1,
                  #labels=c("B","C"),
                  label_size=42)
    png(file=paste0(outdir,"/part1_sub_rmse.png"),width=500,height=heightvar)
    plot(fig1p1)
    dev.off()       
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub_rmse.pdf"), width = 6, height = heightvar2)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub_rmse.svg"), width = 6, height = heightvar2)

    fig1p1 <- plot_grid(p2[[1]],
                  ncol=1,
                  #labels=c("B","C"),
                  label_size=42)
    png(file=paste0(outdir,"/part1_super_rmse.png"),width=500,height=heightvar)
    plot(fig1p1)
    dev.off() 
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super_rmse.pdf"), width = 6, height = heightvar2)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super_rmse.svg"), width = 6, height = heightvar2)

    fig1p1 <- plot_grid(p1[[2]],
                    ncol=1,
                    #labels=c("B","C"),
                    label_size=42)
    png(file=paste0(outdir,"/part1_sub_san_rmse.png"),width=500,height=400)
    plot(fig1p1)
    dev.off()
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub_san_rmse.pdf"), width = 6, height = 5)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_sub_san_rmse.svg"), width = 6, height = 5)

    fig1p1 <- plot_grid(p2[[2]],
                    ncol=1,
                    #labels=c("B","C"),
                    label_size=42)
    png(file=paste0(outdir,"/part1_super_san_rmse.png"),width=500,height=400)
    plot(fig1p1)
    dev.off()
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super_san_rmse.pdf"), width = 6, height = 5)
    cowplot::ggsave2(fig1p1, file = paste0(outdir,"/part1_super_san_rmse.svg"), width = 6, height = 5)

}

plot_fig1_plots_overall_oldmean <- function(homedir = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/"),
                            cell_anno = fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE),
                            genemap = fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F),
                            sanitygene = "STAMBPL1", # gene that corresponds to SANITY threshold in SSv4 data (see below)
                            outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/"),
                            donornames#,
                           # san_mean # vector of mean of counts for each gene (from all donors)
                            ){
    
    #############################
    # Load modeling data and plot
    #############################
    if(!dir.exists(outdir)){
        dir.create(outdir, recursive=T)
    }

    plotdflist_sub <- qread(paste0(homedir,"/fig1bdflist_subclass.qs"))
    plotdflist_super <- qread(paste0(homedir,"/fig1bdflist_supertype.qs"))

    # Create vectors of genes that meet sanity threshold vs not
    #sanityyes <- names(san_mean)[san_mean >= san_mean[names(san_mean)==sanitygene]]
    #sanityno <- names(san_mean)[san_mean < san_mean[names(san_mean)==sanitygene]]

    # Calculate SANITY threshold for each of the three donors
    genemeansfiltvec <- list()
    for(i in 1:3){
        exdir <- list.files(paste0(homedir,"/donor",i,"/SyntheticDatasets/"),full.names=T)
        exdirexpr <- exdir[grep("EXPRLIST",exdir)]
        exdirexpr <- exdirexpr[grep("0pcntVar", exdirexpr)]
        #exdirnames <- unlist(lapply(strsplit(exdirexpr,"_"), function(x) x[8]))
        #exdirnames <- as.numeric(gsub("pcntVar", "", exdirnames))
        exdirnames <- 0

        # genemeansfiltvec[[i]] <- do.call(rbind,mapply(function(x,y){
        #     obj <- readRDS(x) 
        #     obj <- norm_pb_samples(obj[[1]], index=3)
        #     exprt <- t(obj[,3:ncol(obj)])
        #     colnames(exprt) <- obj[,2]
        #     genemeansfilt <- apply(exprt,2,function(x) log2(mean(x)))
            
        #     # code to ensure that mean column reflects same log norm (base 2)
        #     plotdflist_sub[[i]]$mean <- genemeansfilt
        #     plotdflist_super[[i]]$mean <- genemeansfilt

        #     return(data.frame("pcnt.var"=y,
        #                     "mean"=genemeansfilt[names(genemeansfilt)==sanitygene]))
        # }, exdirexpr, exdirnames, SIMPLIFY=F))

        x <- exdirexpr
        y <- exdirnames
        obj <- readRDS(x) 
        obj <- norm_pb_samples(obj[[1]], index=3)
        exprt <- t(obj[,3:ncol(obj)])
        colnames(exprt) <- obj[,2]
        genemeansfilt <- apply(exprt,2,function(x) log2(mean(x)))
        # code to ensure that mean column reflects same log norm (base 2)
        plotdflist_sub[[i]]$mean <- genemeansfilt
        plotdflist_super[[i]]$mean <- genemeansfilt

        genemeansfiltvec[[i]] <- data.frame("pcnt.var"=y,
                                            "mean"=genemeansfilt[names(genemeansfilt)==sanitygene])
        rownames(genemeansfiltvec[[i]]) <- NULL
        genemeansfiltvec[[i]] <- genemeansfiltvec[[i]] %>% dplyr::filter(pcnt.var==0)
    }
    # sanity mean is 1.760, 1.709, 1.691 for donor 1-3

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

    rmsevec_sub <- lapply(plotdflist_sub, function(x){
        df <- x %>% group_by(pcnt.var) %>%
        summarise(mean=mean(rmse)) %>%
        dplyr::filter(pcnt.var==0) %>% 
        mutate(mean_label=paste0("Mean RMSE = ",signif(mean,2)),
           pcnt.var.label=paste0("Pcnt. var = ", pcnt.var)) %>%
        mutate(pcnt.var=as.factor(pcnt.var.label))
    })

    rmsevec_super <- lapply(plotdflist_super, function(x){
        df <- x %>% group_by(pcnt.var) %>%
        summarise(mean=mean(rmse)) %>%
        dplyr::filter(pcnt.var==0) %>%
        mutate(mean_label=paste0("Mean RMSE = ",signif(mean,2)),
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
 
    # Plot first two columns of figure 1
    plot_cols_facetwrap <- function(a,b,c){
        meansvec <- do.call(rbind, b) %>%
            mutate(donor=donornames)
        meansvec$mean_label <- gsub("adj. r2 ", "", meansvec$mean_label)
        out <- mapply(function(i, plotdflist){
            plotdflist[plotdflist$pcnt.var==0,] %>%
            mutate(donor=donornames[i])
        }, c(1:3), a, SIMPLIFY=F) %>% do.call(rbind,.) %>%
            mutate(de_gene=Gene %in% allorigdegenes)
        genemeansall <- out %>% group_by(Gene) %>% summarise(gene_mean=mean(mean))
        sancut <- genemeansall$gene_mean[genemeansall$Gene==sanitygene]
        out$sancut[out$Gene %in% genemeansall$Gene[genemeansall$gene_mean>=sancut]] <- paste0(">=1 mean\nUMI/cell\n(n=", sum(genemeansall$gene_mean>=sancut), ")")
        out$sancut[out$Gene %in% genemeansall$Gene[genemeansall$gene_mean<sancut]] <- paste0("<1 mean\nUMI/cell\n(n=", sum(genemeansall$gene_mean<sancut), ")")

        out1 <- out %>%
            ggplot(aes(y=adj_r2, x=mean)) +
            theme_bw() +
            geom_point(alpha=0.05,size=0.5) +
            stat_density_2d(geom = "polygon", contour = TRUE,
                            aes(fill = after_stat(level)), colour = "black",
                            bins = 5) +
            scale_fill_distiller(palette = "Blues", direction = 1) +
            xlab(bquote("Mean expression ("~log[2]~")")) +
            ylab(bquote("Adjusted R"^2)) +
            geom_hline(data=meansvec, aes(yintercept=mean),linetype="dashed",color="red", linewidth=2, alpha=0.5) +
            ylim(0,1) +
        # ggtitle(c) +
            geom_vline(data=genemeansfiltvec[[i]],aes(xintercept=mean),linetype="dashed",color="blue",linewidth=2, alpha=0.5) +
            geom_text(data=meansvec,aes(label=mean_label, y=mean+0.15), x = -2, color="red", size=10, alpha=0.8) +
            theme(text=element_text(size=30),
                    legend.position="none",
                    axis.title.x=element_text(size=30,margin=margin(10,0,0,0)),
                    axis.title.y=element_text(margin=margin(0,10,0,0)),
                    plot.title=element_text(hjust=0.5),
                    strip.text=element_text(size=24),
                    plot.margin=margin(1,1,1,1,"cm")) +
            facet_wrap(~donor, nrow=3)
        
        out2 <- out %>%
          ggplot(aes(x=sancut, y=adj_r2, fill=sancut)) +
            theme_bw() +
            geom_violin() + 
            labs(x="", y=bquote("Adjusted R"^2)) +
            theme(text=element_text(size=30),
                  legend.position="none",
                  axis.title.x=element_text(size=30,margin=margin(10,0,0,0)),
                  axis.title.y=element_text(margin=margin(0,10,0,0)),
                  plot.title=element_text(hjust=0.5),
                  strip.text=element_text(size=24),
                  plot.margin=margin(1,1,1,1,"cm")) +
            scale_fill_manual(values=c("white", "blue"))# +
            #stat_compare_means(label.y=.85)

        return(list(out1, out2))
    }
}

plot_fig2_plots_indiv <- function(homedir = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/"),
                            cell_anno = fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE),
                            genemap = fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F),
                            ctcolvec = 1, # column of cell_anno that indicates subclass, used for indiv gene graphs for fig1
                            outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/"),
                            cellidcol=3,
                            inputgenevec1 = c("AIF1","ALDH1L1","MOG","SLC17A7","GAD1"),
                            inputgenevec1_names =c ("AIF1 (IBA1)","ALDH1L1", "MOG","SLC17A7","GAD1"),
                            save_name="fig2_part2"
){
     
    if(!dir.exists(outdir)){
      dir.create(outdir, recursive=T)
    }
    #############################
    # Load modeling data and plot
    #############################
    plotdflist_sub <- qread(paste0(homedir,"/fig1bdflist_subclass.qs"))
    plotdflist_super <- qread(paste0(homedir,"/fig1bdflist_supertype.qs"))

    # Load DE genes from Lein et al
    origdemat <- fread(data.table=F, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein_orig_de_genes.csv"))
    origdelist <- tapply(origdemat[,1], origdemat[,2], list)
    allorigdegenes <- unique(unlist(origdelist))

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
        colnames(escale) <- obj[,2]
        escale <- escale[,colnames(escale) %in% inputgenevec1]
        exprZero[[i]] <- escale
        sif <- fread(data.table=F,file=exdirsif)
        sampct <- t(calc_CT_counts(sif, cell_anno, ct_col=ctcolvec, cell_id_col = cellidcol))
        sampctmean <- apply(sampct, 2, mean)
        sampct <- sampct[,order(sampctmean, decreasing=T)] # order by mean num of nuclei per subtype
            
        # Run modeling
        exdirlist[[i]] <- lapply(as.list(as.data.frame(escale)), function(x) lm(x ~ ., data=as.data.frame(sampct)))
    }

    # Create graphs for rightmost columns of figure 1
    graph_marker_model <- function(gene){
        plotdf <- do.call(rbind,mapply(function(x,y,z){
            index <- which(colnames(y)==gene)
            df <- data.frame("pred"=predict(x[[index]]),"actual"=y[,index], "gene"=gene, "donor"=paste0("Individual ",z))
        }, exdirlist, exprZero,seq_along(exdirlist), SIMPLIFY=F))
        return(plotdf)
    }

    # Calculate mean expr percentiles for all genes in all donors
    donor_mean_per <- lapply(plotdflist_sub, function(x){
      x[,colnames(x) %in% c("mean", "Gene")] %>%
        mutate("percentile"=ecdf(mean)(mean)*100)
    })
    donor_mean_all <- do.call(rbind, donor_mean_per) %>% group_by(Gene) %>%
    summarise(mean_per = mean(percentile))

    plot_gene_plots_facetwrap <- function(x, graphtitle=NULL){
        if(is.null(graphtitle)){
            graphtitle <- x
        }

        df <- graph_marker_model(x) 
        dfcor <- df %>% group_by(donor) %>% summarise(cors=cor(as.numeric(pred),as.numeric(actual))) %>%
            mutate(corlabel=paste0("~R^{2} == ", signif(cors^2,2)),
                predlabel=donor_mean_all$mean_per[donor_mean_all$Gene == x])

        # Add mean expression percentiles to plot df

        p <- ggplot(df, aes(x=pred, y=actual)) + 
            theme_bw() +
            geom_density_2d_filled() +
            geom_abline() +
            theme(text=element_text(size=30),
                axis.title.x=element_text(size=36),
                axis.title.y=element_text(size=36),
                axis.text.x=element_text(size=26),
                axis.text.y=element_text(size=30),
                plot.subtitle=element_text(size=25, hjust=0.5),
                strip.text=element_text(size=36),
                plot.title=element_text(size=40,hjust=0.5,face="bold"),
                legend.position="none") +
                #plot.margin=margin(1,1,1,1,"cm")) +
            labs(x="", y="",
                title=graphtitle, 
                subtitle=bquote(.(as.integer(dfcor$predlabel[1]))^"th"~"%tile mean expr")) +
            geom_text(data=dfcor,aes(label=corlabel), y=Inf, x = -Inf,fontface="bold", color="white", size=12, 
            alpha=0.8, vjust=2,hjust=-0.1,parse=T) +
            facet_wrap(~donor, nrow=3)
        return(p)
    }

    #inputgenevec2 <- c("VIP", "SST", "LAMP5", "LHX6", "PAX6", "GAD2", "SLC32A1", "SLC17A6","SLC17A8", "RBFOX3")
    pgenelist1 <- mapply(plot_gene_plots_facetwrap, inputgenevec1, inputgenevec1_names, SIMPLIFY=F)
    #pgenelist2 <- lapply(inputgenevec2, plot_gene_plots_facetwrap)
    #qsave(pgenelist1, file=paste0(outdir,"/fig1_part2_plotobject.qs"))
    #qsave(pgenelist2, file=paste0(outdir,"/fig1_part2_sup_plotobject.qs"))

    title1 <- ggdraw() + draw_label("Predicted expression",size=42)
    title2 <- ggdraw() + draw_label("Actual expression",size=42, angle=90)

    pglist1 <- plot_grid(plotlist=pgenelist1, ncol=length(pgenelist1))
    fig1p2 <- plot_grid(pglist1, 
                    ncol=1,
                    label_size=42)
    fig1p2 <- plot_grid(fig1p2,title1,NULL, nrow=3, rel_heights=c(1,0.005,0.05))
    fig1p2 <- plot_grid(title2, fig1p2, ncol=2, rel_widths=c(0.03,1))
    # cairo_pdf(file=paste0(outdir,"/",save_name,".pdf"),width=1800,height=1200)
    # plot(fig1p2)
    # dev.off()

    cowplot::ggsave2(fig1p2, file = paste0(outdir,"/",save_name,".pdf"),width=24,height=16)
    cowplot::ggsave2(fig1p2, file = paste0(outdir,"/",save_name,".svg"),width=24,height=16)

    # png(file=paste0(outdir,"/",save_name,".png"),width=1800,height=1200)
    # plot(fig1p2)
    # dev.off()


    # pglist2 <- plot_grid(plotlist=pgenelist2[1:5], ncol=length(pgenelist2)/2)
    # fig1p2sup <- plot_grid(pglist2, 
    #                 ncol=1,
    #                 label_size=42)
    # fig1p2sup <- fig1p2sup + labs(x="Predicted expression", y="Actual expression")
    # png(file=paste0(outdir,"/fig1_part2_sup.png"),width=1800,height=1200)
    # plot(fig1p2sup)
    # dev.off()

    # pglist2 <- plot_grid(plotlist=pgenelist2[6:10], ncol=length(pgenelist2)/2)
    # fig1p2sup <- plot_grid(pglist2, 
    #                 ncol=1,
    #                 label_size=42)
    # fig1p2sup <- fig1p2sup + labs(x="Predicted expression", y="Actual expression")
    # png(file=paste0(outdir,"/fig1_part2_sup2.png"),width=1800,height=1200)
    # plot(fig1p2sup)
    # dev.off()
                           
}

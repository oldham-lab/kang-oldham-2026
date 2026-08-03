source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))

# Run DFC (cellbender-corrected)
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC_cellbender.csv"), data.table = FALSE)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_cellbender/")
outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_cellbender/")
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/cellbender/gene_count_means_byDonor.qs"))

plot_fig1_plots_overall(homedir=homedir,
                        cell_anno=cell_anno,
                        genemap=genemap,
                        outdir=outdir,
                        donornames=donornames,
                        san_mean=san_mean)
                        
plot_fig1_plots_indiv(homedir=homedir,
                      cell_anno=cell_anno,
                      genemap=genemap,
                      ctcolvec=1,
                      outdir=outdir,
                      cellidcol=14)


# Load adj r2 before and after cellbender
# Filter out genes expressed in less than three cells
low_genes <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/genes_expressed_in_lessThan3_cells.qs"))
add_gene_symbols <- function(homedir){
    plotdflist_sub <- readRDS(paste0(homedir,"/fig1bdflist_subclass.RDS"))
    plotdflist_super <- readRDS(paste0(homedir,"/fig1bdflist_supertype.RDS"))

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
        mutate(Gene=y) %>%
        dplyr::filter(!Gene %in% low_genes)
    }, plotdflist_sub, genesymlist,1:3, SIMPLIFY=F)

    plotdflist_super <- mapply(function(x,y,z){
    x %>% dplyr::filter(pcnt.var==0) %>%
        mutate(Gene=y) %>%
        dplyr::filter(!Gene %in% low_genes)
    }, plotdflist_super, genesymlist,1:3, SIMPLIFY=F)
    return(list(plotdflist_sub, plotdflist_super))
}

before <- add_gene_symbols(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/"))
after <- add_gene_symbols(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_cellbender/"))

namevec <- c("subclass", "supertype")

for(i in 1:2){
  for(j in 1:3){
    colnames(after[[i]][[j]])[1] <- "after_adj_r2"
    plotdf <- before[[i]][[j]] %>% inner_join(after[[i]][[j]][,c(1,6)])

    p <- ggplot(plotdf,aes(x=adj_r2, y=after_adj_r2)) +
      th3 + 
      geom_point() + 
      labs(x=bquote("Adj. R"^2~"before CellBender"), 
           y=bquote("Adj. R"^2~"after CellBender")) +
      annotate("text",label=paste0("R=",signif(cor(plotdf$adj_r2, plotdf$after_adj_r2),3)), x=0.15, y=.9, color="red", size=12) 
    ggsave(p, file=paste0(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/cellbender_scatter/dfc/"), namevec[i], "_donor",j,".png"))
  }
}

# Analyze outliers
# After supertype H18.30.002 shows large outlier (bottom right)
i=2
j=1
head(before[[i]][[j]])
head(after[[i]][[j]])
leftdf <- before[[i]][[j]] %>% dplyr::select(adj_r2, Gene)
rightdf <- after[[i]][[j]] %>% dplyr::select(adj_r2, Gene) %>% rename(r2_after=adj_r2)
combdf <- inner_join(leftdf, rightdf) %>% mutate(r2diff=r2_after-adj_r2) %>% arrange(r2diff)
head(combdf)
#   adj_r2   Gene    r2_after     r2diff
#1 0.9998071 OR2T12 0.008188188 -0.9916189  Olfactory receptor 2
# Mostly explained by L6 IT_1
# 99% correlation with L6 IT_1 abundance
# OR2T12 is expressed in two cells (1 transcript each)
# Remove genes expressed in less than 3 cells
# This removes about ~80 genes per donor


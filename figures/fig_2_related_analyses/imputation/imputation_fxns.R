library(ggrepel)

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))


summarise_imputation_results <- function(homedir_impute=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI/"),
                                         impute_type="scVI",
                                         savedir=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/scVI/")
                                         ){

    if(!dir.exists(savedir)){dir.create(savedir,recursive=T)}

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
    
    scvi_sub <- qread(paste0(homedir_impute,"/fig1bdflist_subclass.qs"))
    san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means_byDonor.qs"))

    plotdf_sub1 <- mapply(function(x,y,z,a){
    out <- left_join(x,y[,c(1,6)],by=join_by("Gene"=="Gene")) %>%
        mutate("Donor"=z)
    a <- data.frame("san_mean"=a) %>% rownames_to_column(var="Gene")
    out <- left_join(out, a, by=join_by("Gene"=="Gene"))
    return(out)
    }, scvi_sub, plotdflist_sub, donornames,san_mean, SIMPLIFY=F) %>% do.call(rbind,.) %>%
    mutate(adj_r2_diff = adj_r2.x-adj_r2.y, san_mean=log2(san_mean))
    
    plotdf_sub <- mapply(function(x,y,z,a){
    y <- y[y$Gene %in% x$Gene,]
    x <- x[x$Gene %in% y$Gene,]
    x$mean <- y$mean[match(x$Gene,y$Gene)]
    out <- rbind(x %>% mutate("type"=impute_type), y %>% mutate("type"="Original")) %>% 
        mutate(Donor=z, type=factor(type,levels=c("Original",impute_type)))
    a <- data.frame("san_mean"=a) %>% rownames_to_column(var="Gene")
    out <- left_join(out, a, by=join_by("Gene"=="Gene"))
    return(out)
    }, scvi_sub, plotdflist_sub, donornames,san_mean, SIMPLIFY=F) %>% do.call(rbind,.)
    
    # plot all genes adj r2 (before vs after)
    p <- ggplot(plotdf_sub, aes(x=type, y=adj_r2)) +
    theme_bw() +
    geom_violin(aes(fill=type), alpha=0.2) + 
    # geom_line(aes(group=Gene)) +
    #geom_point() +
    facet_wrap(~Donor) +
    theme(text=element_text(size=30),
            legend.position="none") +
    labs(x="", y=bquote("Adjusted R"^2))
    ggsave(p,file=paste0(savedir,"/all_genes.png"), width=10)

    # Plot adj r2 vs mean expr
    p <- ggplot(plotdf_sub1, aes(x=adj_r2_diff, y=san_mean)) + 
    theme_bw() +
    geom_point(aes(color=type)) +
    geom_smooth() +
    facet_wrap(~Donor) +
    labs(x=bquote(Delta~" Adjusted R"^2~" ("~.(impute_type)~" - Original )"), y=bquote("Mean expression ("~log[2]~")")) +
    theme(text=element_text(size=30),
            legend.position="none",
            plot.margin=margin(1,1,1,1,"cm")) +
    scale_x_continuous(labels = function(x) format(x, nsmall = 2)) 
    ggsave(p,file=paste0(savedir,"/adjr2_vs_mean.png"), width=14)

    # Plot after filtering to genes of interest
    inputgenes <- c("AIF1","ALDH1L1","MOG","SLC17A7","GAD1","VIP", "SST", "LAMP5", "LHX6", "PAX6", "GAD2", "SLC32A1", "SLC17A6","SLC17A8", "RBFOX3")
    p <- plotdf_sub %>% 
    dplyr::filter(Gene %in% inputgenes) %>%
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
    ggsave(p,file=paste0(savedir,"/genes_of_interest.png"), width=10)
}
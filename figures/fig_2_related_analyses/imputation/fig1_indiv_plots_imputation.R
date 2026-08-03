# Recreate Fig1 individual plots for CB + scVI data

source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_2_related_analyses/Fig1-functions.R"))

# save directory
outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_vs_cbscVI_impute/")

# Load data
homedir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI_cellbender/")
cell_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC_cellbender.csv"), data.table = FALSE)
genemap <- fread(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/Lein2023_DFC_genemap.csv"),data.table=F)
donornames <- c("H18.30.002", "H19.30.001", "H19.30.002")
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means_byDonor.qs"))

# Recreate Fig1b
plot_fig1_plots_indiv(homedir=homedir,
                      cell_anno=cell_anno,
                      genemap=genemap,
                      ctcolvec=1,
                      outdir=outdir,
                      cellidcol=14)

# Recreate Fig4
inputgenes_unordered <- c("VIP", "SST", "LAMP5", "LHX6", "PAX6", "GAD2", "SLC32A1", "SLC17A6","SLC17A8", "RBFOX3")
inputgenes_ordered <- san_mean[[1]][names(san_mean[[1]]) %in% inputgenes_unordered] %>% sort %>% names

plot_fig1_plots_indiv(homedir=homedir,
                      cell_anno=cell_anno,
                      genemap=genemap,
                      ctcolvec=1,
                      outdir=outdir,
                      cellidcol=14,
                      inputgenevec1 = inputgenes_ordered[1:5],
                      inputgenevec1_names = inputgenes_ordered[1:5],
                      save_name="figs4_part1")

plot_fig1_plots_indiv(homedir=homedir,
                      cell_anno=cell_anno,
                      genemap=genemap,
                      ctcolvec=1,
                      outdir=outdir,
                      cellidcol=14,
                      inputgenevec1 = inputgenes_ordered[6:10],
                      inputgenevec1_names = inputgenes_ordered[6:10],
                      save_name="figs4_part2")

# Recreate fig 1a
plot_fig1_plots_overall(homedir=homedir,
                        cell_anno=cell_anno,
                        genemap=genemap,
                        outdir=outdir,
                        donornames=donornames,
                        san_mean=san_mean)

## Plot modeling results of pre- vs post-correction
outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_vs_cbscVI_impute/")

homedirs <- list(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/"),
                 file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI_cellbender/"))

plotdflist_sub <- lapply(homedirs, function(homedir) qread(paste0(homedir,"/fig1bdflist_subclass.qs")))
plotdflist_super <- lapply(homedirs, function(homedir) qread(paste0(homedir,"/fig1bdflist_supertype.qs")))

# Subset to common genes
common_sub <- lapply(plotdflist_sub, function(y) lapply(y, function(x) x$Gene) %>% Reduce(intersect,.)) %>% Reduce(intersect,.)
plotdflist_sub <- lapply(plotdflist_sub, function(y) lapply(y, function(x){
    out <- x[x$Gene %in% common_sub,]
    out <- out[match(common_sub, out$Gene),]
    return(out)
}))
common_super <- lapply(plotdflist_super, function(y) lapply(y, function(x) x$Gene) %>% Reduce(intersect,.)) %>% Reduce(intersect,.)
plotdflist_super <- lapply(plotdflist_super, function(y) lapply(y, function(x){
    out <- x[x$Gene %in% common_sub,]
    out <- out[match(common_sub, out$Gene),]
    return(out)
}))

# Graph subclass
dflist <- list()
for(i in 1:3){
    dflist[[i]] <- data.frame("before"=plotdflist_sub[[1]][[i]][,1],
                              "after"=plotdflist_sub[[2]][[i]][,1],
                              "donor"=paste0("Donor ",i))
}
dflist <- dflist %>% do.call(rbind,.)
 
p <- ggplot(dflist, aes(x=before, y=after)) +
  theme_bw() +
  geom_point(alpha=0.05,size=0.5) + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                      bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        plot.margin=margin(1,1,1,1,"cm"),
        legend.position="none") +
  facet_wrap(~donor, nrow=3) +
  labs(x=bquote("Adjusted R"^2~"(No correction)"), y="Adjusted R"^2~"(CB_scVI)")
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_vs_cbscVI_impute/subclass_modeling_LeinDFC_before_and_after_CB_scVI.png"), height=12)

# Graph supertype
dflist <- list()
for(i in 1:3){
    dflist[[i]] <- data.frame("before"=plotdflist_super[[1]][[i]][,1],
                              "after"=plotdflist_super[[2]][[i]][,1],
                              "donor"=paste0("Donor ",i))
}
dflist <- dflist %>% do.call(rbind,.)

p <- ggplot(dflist, aes(x=before, y=after)) +
  theme_bw() +
  geom_point(alpha=0.05,size=0.5) + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                      bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        plot.margin=margin(1,1,1,1,"cm"),
        legend.position="none") +
  facet_wrap(~donor, nrow=3) +
  labs(x=bquote("Adjusted R"^2~"(No correction)"), y="Adjusted R"^2~"(CB_scVI)")
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_vs_cbscVI_impute/supertype_modeling_LeinDFC_before_and_after_CB_scVI.png"), height=12)


## Repeat these with ALRA
outdir = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_vs_cbscVI_impute/")

homedirs <- list(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor/"),
                 file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_ALRA_cellbender/"))

plotdflist_sub <- lapply(homedirs, function(homedir) qread(paste0(homedir,"/fig1bdflist_subclass.qs")))
plotdflist_super <- lapply(homedirs, function(homedir) qread(paste0(homedir,"/fig1bdflist_supertype.qs")))

# Subset to common genes
common_sub <- lapply(plotdflist_sub, function(y) lapply(y, function(x) x$Gene) %>% Reduce(intersect,.)) %>% Reduce(intersect,.)
plotdflist_sub <- lapply(plotdflist_sub, function(y) lapply(y, function(x){
    out <- x[x$Gene %in% common_sub,]
    out <- out[match(common_sub, out$Gene),]
    return(out)
}))
common_super <- lapply(plotdflist_super, function(y) lapply(y, function(x) x$Gene) %>% Reduce(intersect,.)) %>% Reduce(intersect,.)
plotdflist_super <- lapply(plotdflist_super, function(y) lapply(y, function(x){
    out <- x[x$Gene %in% common_sub,]
    out <- out[match(common_sub, out$Gene),]
    return(out)
}))

# Graph subclass
dflist <- list()
for(i in 1:3){
    dflist[[i]] <- data.frame("before"=plotdflist_sub[[1]][[i]][,1],
                              "after"=plotdflist_sub[[2]][[i]][,1],
                              "donor"=paste0("Donor ",i))
}
dflist <- dflist %>% do.call(rbind,.)

p <- ggplot(dflist, aes(x=before, y=after)) +
  theme_bw() +
  geom_point(alpha=0.05,size=0.5) + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                      bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        plot.margin=margin(1,1,1,1,"cm"),
        legend.position="none") +
  facet_wrap(~donor, nrow=3) +
  labs(x=bquote("Adjusted R"^2~"(No correction)"), y="Adjusted R"^2~"(CB_ALRA)")
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_vs_cbscVI_impute/subclass_modeling_LeinDFC_before_and_after_CB_ALRA.png"), height=12)

# Graph supertype
dflist <- list()
for(i in 1:3){
    dflist[[i]] <- data.frame("before"=plotdflist_super[[1]][[i]][,1],
                              "after"=plotdflist_super[[2]][[i]][,1],
                              "donor"=paste0("Donor ",i))
}
dflist <- dflist %>% do.call(rbind,.)

p <- ggplot(dflist, aes(x=before, y=after)) +
  theme_bw() +
  geom_point(alpha=0.05,size=0.5) + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                      bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        plot.margin=margin(1,1,1,1,"cm"),
        legend.position="none") +
  facet_wrap(~donor, nrow=3) +
  labs(x=bquote("Adjusted R"^2~"(No correction)"), y="Adjusted R"^2~"(CB_ALRA)")
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/dfc_10x_vs_cbscVI_impute/supertype_modeling_LeinDFC_before_and_after_CB_ALRA.png"), height=12)

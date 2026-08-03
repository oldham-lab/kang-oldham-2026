# Run differential co-expression FM between original Lein DFC and cellbender-scVI-corrected Lein DFC

library(tidyverse)
library(qs)
library(data.table)
library(RColorBrewer)
library(reticulate)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
numpy <- import("numpy")
options(bitmapType = 'cairo')

## Load data
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

# Load adjacency matrices for differential co-expression
simMatmega <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/megamat_adjMat.qs"))
simMatleinpb <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/leinpb_donor1_adjMat.qs"))
simMatscvi <- qread(file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "figures/fig1/imputation/diffcoex/scvi_donor1_adjMat.qs"))

## Compare adjacency distributions of three datasets
brewer_colors <- brewer.pal(6,"Paired")[c(2,4,6)]
plotdf <- data.frame("type"=c(rep("Bulk megaset",10000), rep("Lein PB", 10000), rep("Lein PB + CB + scVI",10000)),
                     "cors"=c(sample(simMatmega[upper.tri(simMatmega)], 10000)*2-1,
                              sample(simMatleinpb[upper.tri(simMatleinpb)], 10000)*2-1,
                              sample(simMatscvi[upper.tri(simMatscvi)], 10000)*2-1)) %>%
  mutate(type=factor(type, levels=unique(type)))

p <- ggplot(plotdf, aes(x=cors, color=type)) +
  theme_bw() +
  geom_density(linewidth=1) +
  theme(text=element_text(size=30),
        legend.position="bottom",
        legend.title=element_blank(),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        legend.direction="vertical",
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="Pearson correlation") +
  scale_color_manual(values=brewer_colors)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/cor_dists.png"))

# Compare adjacency distributions after subtraction
mat1 <- simMatmega - simMatleinpb
mat2 <- simMatleinpb - simMatmega
mat3 <- simMatmega - simMatscvi
mat4 <- simMatscvi - simMatmega
mat5 <- simMatscvi - simMatleinpb
mat6 <- simMatleinpb - simMatscvi
brewer_colors2 <- brewer.pal(6,"Paired")[c(1,2)]
brewer_colors3 <- brewer.pal(12,"Paired")[c(2,1,6,5,4,3)]

plotdf <- data.frame("type"=c(rep("Bulk minus PB",10000), rep("PB minus Bulk", 10000), rep("Bulk minus PB_CB_scVI",10000),
                              rep("PB_CB_scVI minus Bulk",10000), rep("PB_CB_scVI minus PB", 10000), rep("PB minus PB_CB_scVI",10000)),
                     "cors"=c(sample(mat1[upper.tri(mat1)], 10000),
                              sample(mat2[upper.tri(mat2)], 10000),
                              sample(mat3[upper.tri(mat3)], 10000),
                              sample(mat4[upper.tri(mat4)], 10000),
                              sample(mat5[upper.tri(mat5)], 10000),
                              sample(mat6[upper.tri(mat6)], 10000))) %>%
  mutate(type=factor(type, levels=unique(type)))

p <- ggplot(plotdf, aes(x=cors, color=type)) +
  theme_bw() +
  geom_density(linewidth=1) +
  theme(text=element_text(size=30),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.direction="vertical",
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="Adjacency") +
  scale_color_manual(values=brewer_colors3)
ggsave(p, file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/adj_dists_subtract.png"))


## Calculate correlations between genes

# Load means for original SN data
san_mean <- qread(file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/gene_count_means_byDonor.qs"))
san1 <- san_mean[[1]][names(san_mean[[1]]) %in% common_genes]
san1 <- san1[match(scvi[,2], names(san1))]

# Calculate pcntiles for san1
calc_pcntile <- function(vec){
    (1-(rank(-vec,na.last="keep")/max(rank(-vec,na.last="keep"))))*100
}
sandf <- data.frame("Gene"=names(san1),
                    "Pcntile"=calc_pcntile(san1),
                    "Bracket"=NA)
sandf$Bracket[sandf$Pcntile<20] <- "0-20"
sandf$Bracket[sandf$Pcntile>=20 & sandf$Pcntile < 40] <- "20-40"
sandf$Bracket[sandf$Pcntile>=40 & sandf$Pcntile < 60] <- "40-60"
sandf$Bracket[sandf$Pcntile>=60 & sandf$Pcntile < 80] <- "60-80"
sandf$Bracket[sandf$Pcntile>=80 & sandf$Pcntile < 100] <- "80-99"

megat <- t(mega[,-c(1,2)])
leinpbt <- t(leinpb[,-c(1,2)])
scvit <- t(scvi[,-c(1,2)])

df1 <- data.frame("gene"=common_genes,
                  "cors"=mapply(function(a,b) cor(a,b), as.list(as.data.frame(megat)), as.list(as.data.frame(leinpbt)), SIMPLIFY=F) %>% unlist,
                  sandf[,2:3],
                  "type"="Bulk vs PB")
df2 <- data.frame("gene"=common_genes,
                  "cors"=mapply(function(a,b) cor(a,b), as.list(as.data.frame(megat)), as.list(as.data.frame(scvit)), SIMPLIFY=F) %>% unlist,
                  sandf[,2:3],
                  "type"="Bulk vs PB_CB_scVI")
df3 <- data.frame("gene"=common_genes,
                  "cors"=mapply(function(a,b) cor(a,b), as.list(as.data.frame(leinpbt)), as.list(as.data.frame(scvit)), SIMPLIFY=F) %>% unlist,
                  sandf[,2:3],
                  "type"="PB vs PB_CB_scVI")

dfall <- rbind(df1, df2, df3) %>%
  mutate(type=factor(type, levels=unique(type)))
dfall$Bracket <- factor(dfall$Bracket, levels=c("0-20", "20-40", "40-60", "60-80", "80-99"))

p <- ggplot(dfall, aes(x=Bracket, y=cors)) +
        theme_bw() +
        geom_violin(aes(fill=Bracket)) + 
        scale_fill_brewer() +
        geom_boxplot(fill="white",notch=T, outlier.shape=NA, width=0.2) +
        facet_wrap(~type, ncol=3) +
        theme(text=element_text(size=30),
              plot.margin=margin(1,1,1,1,"cm"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1),
              legend.position="none") +
        labs(y="Pearson correlation", x="Mean UMI percentile")
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/exprVec_cors_by_pcntile_bracket.png"), width=14)

p <- ggplot(dfall, aes(x=type, y=cors, fill=type)) +
        theme_bw() +
        geom_violin(aes(fill=type),alpha=0.3) +
        geom_boxplot(fill="white",notch=T, width=0.1, outlier.shape=NA) +
        #geom_jitter(alpha=0.01) +
        theme(text=element_text(size=30),
              plot.margin=margin(1,1,1,1,"cm"),
              legend.position="none") +
        labs(y="Pearson\ncorrelation", x="") +
        scale_fill_manual(values=brewer_colors) 

ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/exprVec_cors.png"), width=10,height=5)

## Plot scatterplots of vectorized genomewide adjacencies
set.seed(12)
df_scatter <- data.frame("Bulk"=simMatmega[upper.tri(simMatmega)]*2-1,
                         "PB"=simMatleinpb[upper.tri(simMatleinpb)]*2-1,
                         "PB_CB_scVI"=simMatscvi[upper.tri(simMatscvi)]*2-1)
sample_inds <- sample(1:nrow(df_scatter), 100000)
df_scatter <- df_scatter[sample_inds,]

val <- signif(cor(df_scatter$Bulk, df_scatter$PB), 2)
p1 <- ggplot(df_scatter, aes(x=Bulk, y=PB)) + 
  theme_bw() + 
  geom_point(alpha=0.1,shape=16, fill="black") + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        legend.position="none",
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="Bulk correlation", y="PB correlation") +
  annotate("text",label=bquote(italic("R")~"="~.(val)), x=Inf, y=-Inf, color="black", hjust=1.1,vjust=-.5, size=12)

val <- signif(cor(df_scatter$Bulk, df_scatter$PB_CB_scVI), 2)
p2 <- ggplot(df_scatter, aes(x=Bulk, y=PB_CB_scVI)) + 
  theme_bw() + 
  geom_point(alpha=0.1,shape=16, fill="black") + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        legend.position="none",
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="Bulk correlation", y="PB_CB_scVI correlation") +
  annotate("text",label=bquote(italic("R")~"="~.(val)), x=Inf, y=-Inf, color="black", hjust=1.1,vjust=-.5, size=12)

val <- signif(cor(df_scatter$PB, df_scatter$PB_CB_scVI), 2)
p3 <- ggplot(df_scatter, aes(x=PB, y=PB_CB_scVI)) + 
  theme_bw() + 
  geom_point(alpha=0.1,shape=16, fill="black") + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        legend.position="none",
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="PB correlation", y="PB_CB_scVI correlation") +
  annotate("text",label=bquote(italic("R")~"="~.(val)), x=Inf, y=-Inf, color="black", hjust=1.1,vjust=-0.5, size=12)

pall <- plot_grid(p1, p2, p3, ncol=3, align="h", axis="t", rel_widths=c(1,1,1))
ggsave(pall,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/vectorized_correlation_scatter.png"), width=20)

# Plot scatterplots of vectorized genomewide adjacencies but filter for top 60% of genes by mean expr first
set.seed(12)
gf <- which(sandf$Pcntile>=40)
simsub1 <- simMatmega[gf,gf]
simsub2 <- simMatleinpb[gf,gf]
simsub3 <- simMatscvi[gf,gf]

df_scatter <- data.frame("Bulk"=simsub1[upper.tri(simsub1)]*2-1,
                         "PB"=simsub2[upper.tri(simsub2)]*2-1,
                         "PB_CB_scVI"=simsub3[upper.tri(simsub3)]*2-1)
sample_inds <- sample(1:nrow(df_scatter), 100000)
df_scatter <- df_scatter[sample_inds,]

val <- signif(cor(df_scatter$Bulk, df_scatter$PB), 2)
p1 <- ggplot(df_scatter, aes(x=Bulk, y=PB)) + 
  theme_bw() + 
  geom_point(alpha=0.1,shape=16, fill="black") + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        legend.position="none",
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="Bulk correlation", y="PB correlation") +
  annotate("text",label=bquote(italic("R")~"="~.(val)), x=Inf, y=-Inf, color="black", hjust=1.1,vjust=-.5, size=12)

val <- signif(cor(df_scatter$Bulk, df_scatter$PB_CB_scVI), 2)
p2 <- ggplot(df_scatter, aes(x=Bulk, y=PB_CB_scVI)) + 
  theme_bw() + 
  geom_point(alpha=0.1,shape=16, fill="black") + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        legend.position="none",
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="Bulk correlation", y="PB_CB_scVI correlation") +
  annotate("text",label=bquote(italic("R")~"="~.(val)), x=Inf, y=-Inf, color="black", hjust=1.1,vjust=-.5, size=12)

val <- signif(cor(df_scatter$PB, df_scatter$PB_CB_scVI), 2)
p3 <- ggplot(df_scatter, aes(x=PB, y=PB_CB_scVI)) + 
  theme_bw() + 
  geom_point(alpha=0.1,shape=16, fill="black") + 
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",bins = 5) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(text=element_text(size=30),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        legend.position="none",
        plot.margin=margin(1,1,1,1,"cm")) +
  labs(x="PB correlation", y="PB_CB_scVI correlation") +
  annotate("text",label=bquote(italic("R")~"="~.(val)), x=Inf, y=-Inf, color="black", hjust=1.1,vjust=-0.5, size=12)

pall <- plot_grid(p1, p2, p3, ncol=3, align="h", axis="t", rel_widths=c(1,1,1))
ggsave(pall,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/vectorized_correlation_scatter_top60percentofgenesbymean.png"), width=20)

# Plot average number of modules produced for FM
netsum1 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinusLeinPB_Modules/Bicor-noTO_p1_donor1_megaMinusLeinPB_17193_network_statistics.csv"))
netsum1_1 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_LeinPBMinusMega_Modules/Bicor-noTO_p1_donor1_LeinPBMinusMega_17193_network_statistics.csv"))
netsum2 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_megaMinus_cbscVI_Modules/Bicor-noTO_p1_donor1_megaMinus_cbscVI_17193_network_statistics.csv"))
netsum2_1 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_cbscVI_minus_mega_Modules/Bicor-noTO_p1_donor1_cbscVI_minus_mega_17193_network_statistics.csv"))
netsum3 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_cbscVI_minus_leinpb_Modules/Bicor-noTO_p1_donor1_cbscVI_minus_leinpb_17193_network_statistics.csv"))
netsum3_1 <- fread(data.table=F,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/donor1_leinpb_minus_cbscVI_Modules/Bicor-noTO_p1_donor1_leinpb_minus_cbscVI_17193_network_statistics.csv"))

netsum1 <- t(netsum1[c(3,4,8),-1])
netsum1_1 <- t(netsum1_1[c(3,4,8),-1])
netsum2 <- t(netsum2[c(3,4,8),-1])
netsum2_1 <- t(netsum2_1[c(3,4,8),-1])
netsum3 <- t(netsum3[c(3,4,8),-1])
netsum3_1 <- t(netsum3_1[c(3,4,8),-1])

plotdf <- rbind(netsum1,
                netsum1_1, 
                netsum2,
                netsum2_1,
                netsum3,
                netsum3_1) %>% 
            as.data.frame %>%
            mutate("type"=c(rep("Bulk minus PB", nrow(netsum1)), 
                            rep("PB minus Bulk", nrow(netsum1_1)), 
                            rep("Bulk minus PB_CB_scVI",nrow(netsum2)),
                            rep("PB_CB_scVI minus Bulk",nrow(netsum2_1)),
                            rep("PB_CB_scVI minus PB", nrow(netsum3)),
                            rep("PB minus PB_CB_scVI", nrow(netsum3_1))))
colnames(plotdf)[1:3] <- c("Signum", "minsize", "finalModules")
plotdf <- plotdf %>% 
  mutate(minsize=paste0("MinSize: ",minsize),
  type=factor(type,levels=unique(type)))
plotdf[,c(1,3)] <- apply(plotdf[,c(1,3)],2,as.numeric)

#maxplotsize <- plotdf %>% dplyr::filter(!type %in% "PB minus Bulk") %>% dplyr::select(finalModules) %>% max
#plotdf <- plotdf %>% dplyr::filter(finalModules < maxplotsize)

brewer_colorsa <- brewer.pal(6,"Paired")[c(2,1,4,3,6,5)]

p <- ggplot(plotdf, aes(x=finalModules, y=Signum, color=type)) +
  theme_bw() + 
  geom_point(size=3) +
  facet_wrap(~minsize, scales="free_x")+
  theme(text=element_text(size=30),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        plot.margin=margin(1,1,1,1,"cm")) +
  scale_color_manual(values=brewer_colors3) +
  #scale_x_continuous(breaks=seq(0,maxplotsize,100)) +
  labs(x="Final # of modules in network") +
  guides(#size=guide_legend(title="Minsize",theme=theme(title=element_text(size=22)), position="right"),
         color=guide_legend(title=element_blank(),theme=theme(legend.direction="vertical"),override.aes=list(size=5), position="bottom"))

ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/number_of_modules.png"),width=12, height=8)

## Compare correlations for a select number of genes before and after correction
goi <- c("AIF1","ALDH1L1","MOG","SLC17A7","GAD1","VIP", "SST", "LAMP5", "LHX6", "PAX6", "GAD2", "SLC32A1", "SLC17A6","SLC17A8", "RBFOX3")

simb <- simMatmega*2-1
rownames(simb) <- common_genes
colnames(simb) <- common_genes
diag(simb) <- 1
simb <- simb[common_genes %in% goi, common_genes %in% goi]

siml <- simMatleinpb*2-1
rownames(siml) <- common_genes
colnames(siml) <- common_genes
diag(siml) <- 1
siml <- siml[common_genes %in% goi, common_genes %in% goi]

sims <- simMatscvi*2-1
rownames(sims) <- common_genes
colnames(sims) <- common_genes
diag(sims) <- 1
sims <- sims[common_genes %in% goi, common_genes %in% goi]

corcol = colorRamp2(c(-1,0,1), c("blue", "white","red"))

dist1 <- hclust(dist(1-simb), method="complete")

pdf(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/markers_of_interest.pdf"))
Heatmap(simb,name="Cor", col=corcol, 
  cluster_rows=dist1, cluster_columns=dist1)
Heatmap(siml,name="Cor", col=corcol, 
  cluster_rows=dist1, cluster_columns=dist1)
Heatmap(sims,name="Cor", col=corcol, 
  cluster_rows=dist1, cluster_columns=dist1)
dev.off()
# Correlations in bulk are quite different from the pseudobulk
# Correlations of pb_cb_scvi are largely similar to pb, looks a bit "sharper"
# - could look different if we look at most dissimilar genes

# What if we did the same as above but with network genes generated from pb_cb_scvi minus pb diffcoexpress
kme <- fread(data.table=F,file="figures/figure_1/diffCoexpress_impute/donor1_cbscVI_minus_leinpb_Modules/Bicor-no_TO_signum0.447_minSize7_merge_ME_0.85_17193/kME_table_.csv")
modgenes <- tapply(kme$Gene, kme$ModSeed, list) %>% unlist %>% unique

sim1 <- simMatleinpb*2-1
rownames(sim1) <- common_genes
colnames(sim1) <- common_genes
diag(sim1) <- 1
sim1 <- sim1[common_genes %in% modgenes, common_genes %in% modgenes]

sim2 <- simMatscvi*2-1
rownames(sim2) <- common_genes
colnames(sim2) <- common_genes
diag(sim2) <- 1
sim2 <- sim2[common_genes %in% modgenes, common_genes %in% modgenes]

corcol = colorRamp2(c(-1,0,1), c("blue", "white","red"))

dist1 <- hclust(dist(1-sim1), method="complete")

pdf(file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/diffCoexpress_impute/markers_of_interest2.pdf"))
Heatmap(sim1,name="Cor", col=corcol, 
  cluster_rows=dist1, cluster_columns=dist1)
Heatmap(sim2,name="Cor", col=corcol, 
  cluster_rows=dist1, cluster_columns=dist1)
dev.off()
library(tidyverse)
library(qs)
library(data.table)
library(RColorBrewer)


# Load expression data for cb_scVI, cb_MAGIC, cb_ALRA
# donor1

scvi <- readRDS("RNAseq_megaset/13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_scVI_cellbender/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_04-15-41_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS")
magic <- readRDS("RNAseq_megaset/13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_MAGIC_cellbender/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_08-32-53_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS")
alra <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/lein2023/DFC_indiv_donor_ALRA_cellbender/donor1/SyntheticDatasets/SyntheticDataset1_10pcntCells_0pcntVar_1518samples_09-13-44_EXPRLIST_PROCESSED_ZEROVAR_BULKGENESUBSET.RDS"))

# Gather common genes and match genes across datasets
common_genes <- intersect(intersect(scvi[[1]][,2], magic[[1]][,2]), alra[[1]][,2])
scvi <- scvi[[1]][scvi[[1]][,2] %in% common_genes,]
magic <- magic[[1]][magic[[1]][,2] %in% common_genes,]
alra <- alra[[1]][alra[[1]][,2] %in% common_genes,]

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

scvi <- t(scvi[,3:ncol(scvi)])
magic <- t(magic[,3:ncol(magic)])
alra <- t(alra[,3:ncol(alra)])


# Calculate correlations between genes
df1 <- data.frame("gene"=common_genes,
                  "cors"=mapply(function(a,b) cor(a,b), as.list(as.data.frame(scvi)), as.list(as.data.frame(magic)), SIMPLIFY=F) %>% unlist,
                  sandf[,2:4],
                  "type"="scvi_vs_magic")
df2 <- data.frame("gene"=common_genes,
                  "cors"=mapply(function(a,b) cor(a,b), as.list(as.data.frame(scvi)), as.list(as.data.frame(alra)), SIMPLIFY=F) %>% unlist,
                  sandf[,2:4],
                  "type"="scvi_vs_alra")
df3 <- data.frame("gene"=common_genes,
                  "cors"=mapply(function(a,b) cor(a,b), as.list(as.data.frame(magic)), as.list(as.data.frame(alra)), SIMPLIFY=F) %>% unlist,
                  sandf[,2:4],
                  "type"="magic_vs_alra")

dfall <- rbind(df1, df2, df3)
dfall$Bracket <- factor(dfall$Bracket, levels=c("0-20", "20-40", "40-60", "60-80", "80-99"))
brewer_colors <- brewer.pal(6,"Paired")[c(2,3,6)]

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
ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_cor_summary.png"), width=14)

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

ggsave(p,file=file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "figure_1/impute_cor_summary_boxplots.png"), width=10,height=5)

library(qs)
library(data.table)
library(DESeq2)
library(tidyverse)
library(showtext)
showtext_auto()

#################
# DE using DEseq2, Lein
#################
# Run DE
data_save_dir <- file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_pseudobulk_by_donor_all_genes/")
rawcounts <- fread(paste0(data_save_dir,"/Lein_2023_cell_expression_by_donor_subclass_sum.csv"), data.table=F) |>
  tibble::column_to_rownames(var = "V1")
cell_anno_pb <- fread(paste0(data_save_dir,"/Lein_2023_cell_annotations_by_donor_subclass_sum_DFC.csv"),data.table=F)
cell_anno_pb$label <- colnames(rawcounts)
cell_anno_pb$Cell_Type <- factor(cell_anno_pb$Cell_Type)
cell_anno_pb$Donor <- factor(cell_anno_pb$Donor)
cts <- unique(cell_anno_pb$Cell_Type)

# Calculate means over all genes
sn_expr <- readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS"))
sn_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
meanbysub <- lapply(unique(sn_anno$Cell_Type), \(subclass){
  rowMeans(sn_expr[ ,colnames(sn_expr) %in% sn_anno$Cell_ID[sn_anno$Cell_Type == subclass]])
}) |> do.call(what = "cbind")
colnames(meanbysub) <- unique(sn_anno$Cell_Type)

# Calculate mean UMI counts per subclass
gene_means <- apply(sn_expr, 1, mean)
gene_means_log <- log(gene_means + 1)
rm(sn_expr)

reslist <- lapply(cts, \(x){
  anno_temp <- cell_anno_pb |> 
    mutate("designcol" = ifelse(Cell_Type == x, "ct", "all"))
  dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                                colData = anno_temp,
                                design= ~ Donor + designcol)
  dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
  res <- results(dds1, name = resultsNames(dds1)[2])
  resdf <- as.data.frame(res@listData) |>
    `rownames<-`(res@rownames) |>
   # dplyr::filter(!is.na(log2FoldChange)) |>
    arrange(padj)
  return(resdf)
})
names(reslist) <- cts
qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_subclass.qs"))

# Load DE data
edger <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/DE_gene_list_subclass_DFC_blockDonorSum_allgenes.RDS"))
reslist <- qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_subclass.qs"))
reslist <- reslist[match(names(edger), names(reslist))]

el <- lapply(edger, \(x) x$genes[x$FDR < 0.05])
rl <- lapply(reslist, \(x){
  out <- rownames(x)[x$padj < 0.05]
  return(out[!is.na(out)])
})

# Plot venn diagram of DE genes (over all subclasses) shared between edgeR and DESeq2
library(eulerr)
allgenes <- c(unlist(el), unlist(rl)) |> unique()
# eulermat1 <- lapply(el, \(x) allgenes %in% x) |> do.call(what = "cbind")
# eulermat2 <- lapply(rl, \(x) allgenes %in% x) |> do.call(what = "cbind")
# eulermat <- cbind(eulermat1[,1], eulermat2[,1])
# colnames(eulermat) <- c("edgeR", "DESeq2")
eulermat <- data.frame("edger" = allgenes %in% unique(unlist(el)),
                       "DESeq2" = allgenes %in% unique(unlist(rl)))

# with quantities
pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/edgeR_vs_deseq2_leinSubclass_euler.pdf"), bg = "white")
plot(euler(eulermat),
     fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
# no quantitites
pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/edgeR_vs_deseq2_leinSubclass_euler_noLabels.pdf"), bg = "white")
plot(euler(eulermat),
     fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# Plot # of shared and unshared DE genes per subclass (stacked barplot)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!el[[ind]] %in% rl[[ind]]),
                                   sum(!rl[[ind]] %in% el[[ind]]),
                                   length(intersect(rl[[ind]], el[[ind]]))),
                       "Method" = c("edgeR_LRT", "DESeq2_LRT", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Method <- factor(plotmat$Method, levels = c("edgeR_LRT", "DESeq2_LRT", "Both"))

# scale_colour_paletteer_d("khroma::highcontrast")
# scale_color_paletteer_d("khroma::highcontrast")
# scale_fill_paletteer_d("khroma::highcontrast")
# paletteer_d("khroma::highcontrast")

p <- ggplot(plotmat, aes(x = ct, y = value, fill = Method)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/deseq2_edger_overlap_leinSubclass.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)

# Plot tree plot of DE genes shared between methods for each subclass
colvec <- RColorBrewer::brewer.pal(6, "Paired")[c(1,3,5)]
plotmat2 <- plotmat |>
  dplyr::filter(Method == "Both") |>
  mutate("Class" = case_when(
                     ct %in% c("VLMC", "OPC", "Endo", "Astro", "Micro/PVM", "Oligo") ~ "Non-neuronal",
                     ct %in% c("Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip") ~ "GABAergic",
                     TRUE ~ "Glutamatergic"),
         "pcnt" = value/sum(value)*100,
         "label" = paste0(ct, "\n(", signif(pcnt,2), "%)"),
         "color" = case_when(
                     Class == "Non-neuronal" ~ colvec[3],
                     Class == "GABAergic" ~ colvec[1],
                     TRUE ~ colvec[2]),) |>
  arrange(desc(value))
library(treemap)
pdf(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/treemap.pdf"))
treemap::treemap(plotmat2,
                 index=c("Class","label"),
                 vSize="value",
                 type="index",
                 draw=T,
                 inflate.labels = T,
                 fontsize.labels = c(0, 10),
                 lowerbound.cex.labels = 0,
                 palette = RColorBrewer::brewer.pal(6, "Paired")[c(2,4,6)],
                 title = ""
                 ) 
dev.off()

# Pie chart
pdf(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/piechart.pdf"))
pie(plotmat2$value, labels = plotmat2$ct, col = plotmat2$color, border="white", cex = 2)
dev.off()

#################
# DE using DEseq2, SEAAD2024 control samples
#################
rawcounts <- fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/cell_expression_by_donor_sum_subclass_Controls.csv")) |>
  tibble::column_to_rownames(var = "Gene")
rawcounts <- apply(rawcounts, 2, ceiling)
cell_anno_pb <- fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/annotations_by_donor_sum_subclass_Controls.csv"))
cell_anno_pb$ID <- colnames(rawcounts)
cell_anno_pb$Subclass <- factor(cell_anno_pb$Subclass)
cell_anno_pb$Donor <- factor(cell_anno_pb$Donor)
cts <- unique(cell_anno_pb$Subclass)

reslist <- lapply(cts, \(x){
  anno_temp <- cell_anno_pb |> 
    mutate("designcol" = factor(ifelse(Subclass == x, "ct", "all")))

  dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                                colData = anno_temp,
                                design= ~ Donor + designcol)
  dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
  res <- results(dds1, name = resultsNames(dds1)[2])
  resdf <- as.data.frame(res@listData) |>
    `rownames<-`(res@rownames) |>
   # dplyr::filter(!is.na(log2FoldChange)) |>
    arrange(padj)
  return(resdf)
})
names(reslist) <- cts
qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass.qs"))

# Load DE data
edger <- readRDS(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/DE_by_subclass_Controls/data/DE_gene_list_subclass_blockDonor_allgenes.RDS"))
reslist <- qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass.qs"))
reslist <- reslist[match(names(edger), names(reslist))]

el <- lapply(edger, \(x) x$genes[x$FDR < 0.05])
rl <- lapply(reslist, \(x){
  out <- rownames(x)[x$padj < 0.05]
  return(out[!is.na(out)])
})

plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!el[[ind]] %in% rl[[ind]]),
                                   sum(!rl[[ind]] %in% el[[ind]]),
                                   length(intersect(rl[[ind]], el[[ind]]))),
                       "Method" = c("edgeR_LRT", "DESeq2_LRT", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Method <- factor(plotmat$Method, levels = c("edgeR_LRT", "DESeq2_LRT", "Both"))

library(paletteer)

p <- ggplot(plotmat, aes(x = ct, y = value, fill = Method)) +
  theme_minimal() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/deseq2_edger_overlap_seaCon_Subclass.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)

############
# How many DE genes total? (as a %)
##############

# Lein
edger <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/DE_gene_list_supertype_DFC_blockDonor_allgenes.RDS"))
reslist <- qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_subclass.qs"))
reslist <- reslist[match(names(edger), names(reslist))]
pct1 <- length(unique(unlist(lapply(edger, \(x) x$genes[x$FDR < 0.05])))) / nrow(edger[[1]]) 
pct2 <- length(unique(unlist(lapply(reslist, \(x) rownames(x)[x$padj < 0.05])))) / nrow(reslist[[1]]) 
plotdf1 <- data.frame("Dataset" = "Jorstad_Lein_2021",
                     "Method" = c("edgeR_LRT", "DESeq2_LRT"),
                     "Pcnt" = c(pct1, pct2) * 100)

# SEAAD2024
edger <- readRDS(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/DE_by_subclass_Controls/data/DE_gene_list_subclass_blockDonor_allgenes.RDS"))
reslist <- qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass.qs"))
reslist <- reslist[match(names(edger), names(reslist))]
pct1 <- length(unique(unlist(lapply(edger, \(x) x$genes[x$FDR < 0.05])))) / nrow(edger[[1]]) 
pct2 <- length(unique(unlist(lapply(reslist, \(x) rownames(x)[x$padj < 0.05])))) / nrow(reslist[[1]]) 
plotdf2 <- data.frame("Dataset" = "SEAAD2024",
                     "Method" = c("edgeR_LRT", "DESeq2_LRT"),
                     "Pcnt" = c(pct1, pct2) * 100)

# Plot
colors <- RColorBrewer::brewer.pal(3, "Set1")
plotdf <- rbind(plotdf1, plotdf2)
plotdf$Method <- factor(plotdf$Method, levels = unique(plotdf$Method))
p <- ggplot(plotdf, aes(x = Dataset, y = Pcnt, fill = Method)) +
  theme_light() + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5, alpha = 0.6) +
  labs(y = "% of all genes that are DE", x = "") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, margin = margin(0, 10, 0, 0)),
        legend.position = "inside",
        legend.position.inside = c(0.25, 0.85),
        legend.box.background = element_rect(color = "grey", linewidth = 1),
        legend.title = element_blank()) +
  ylim(0, 100) +
  #scale_fill_manual(values = colors[1:2]) +
  scale_fill_paletteer_d("khroma::highcontrast")
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/total_pcnt_of_DE_genes.svg"), bg = "white", height = 3, width = 4)

##############
# Is there a relationship between mean expression and significance?
##################################

# Load Lein data
edger <- readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/DE_gene_list_subclass_DFC_blockDonorSum_allgenes.RDS"))
reslist <- qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_subclass.qs"))
reslist <- reslist[match(names(edger), names(reslist))]

# # Plot CPM vs LR for all genes over each subclass
# plotdf <- Reduce(rbind, edger)
# # lapply(edger, \(x) cor(x$logCPM, x$LR)) |> unlist()
# p <- ggplot(plotdf, aes(x = logCPM, y = LR)) + 
#   geom_point(size = 0.2) + 
#   facet_wrap(~celltype, scales = "free_y")
# ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/cpm_vs_LR.svg"))
# # Not much of an obvious relationship

# Compare the distributions of mean expression levels for genes stratified by the Venn diagram 
# (i.e., “Not DE”, “edgeR”, “DESeq2”, and “edgeR and DESeq2”)
edgegenes <- unlist(el) |> unique()
deseqgenes <- unlist(rl) |> unique()
allsharedgenes <- intersect(edger[[1]]$genes, rownames(reslist[[1]])) |> unique()
shareddegenes <- intersect(edgegenes, deseqgenes) |> unique()
genecatlist <- list(
  "Not DE" = allsharedgenes[!allsharedgenes %in% shareddegenes],
  "edgeR\nonly" = edgegenes[(edgegenes %in% allsharedgenes) & !(edgegenes %in% shareddegenes)],
  "DEseq2\nonly" = deseqgenes[(deseqgenes %in% allsharedgenes) & !(deseqgenes %in% shareddegenes)],
  "Both" = shareddegenes[shareddegenes %in% allsharedgenes]
)

plotdf_mean <- mapply(\(x, y){
  data.frame("which" = y,
             "means" = gene_means_log[names(gene_means_log) %in% x])
}, genecatlist, names(genecatlist), SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(which = factor(which, levels = unique(which)))

# raincloud plot
p <- ggplot(plotdf_mean, aes(x = which, y = means, fill = which, color = which)) + 
  theme_classic() + 
  ggrain::geom_rain(alpha = 0.6,
                    point.args = list(size = 0.1, alpha = 0.3),
                    boxplot.args = list(outlier.shape = NA, color = "black", notch = T))  + 
  geom_hline(yintercept = log(1), color = "red", linetype = "dashed", alpha = 0.8) +
  labs(x = "", y = "log(mean expression)") +
  theme(text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast"))) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast")))
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/DEgene_meanExpr_violin.pdf"), width = 4, height = 2.5)
# boxplot with geom_jitter
p <- ggplot(plotdf_mean, aes(x = which, y = means, fill = which, color = which)) + 
  theme_classic() + 
  geom_jitter(alpha = 0.1) + 
  geom_boxplot(notch = T, color = "black", outlier.shape = NA) +
  geom_hline(yintercept = log(1), color = "red", linetype = "dashed", alpha = 0.8) +
  labs(x = "", y = "log(mean expression)") +
  theme(text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast"))) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast")))
ggsave(p, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/DEgene_meanExpr_boxjitter.pdf"), width = 3, height = 2.5)


# Plot pctg of genes above sanity threshold
plotdf_sanity <- lapply(genecatlist, \(x){
  sum(gene_means_log[names(gene_means_log) %in% x] > 0) / length(x) * 100
}) |>
  unlist() |> 
  as.data.frame() |>
  tibble::rownames_to_column(var = "type") |>
  mutate(type = factor(type, levels = unique(type)))
colnames(plotdf_sanity) <- c("type", "val")

p2 <- ggplot(plotdf_sanity, aes(x = type, y = val, fill = type)) +
  theme_classic() + 
  geom_bar(stat = "identity", alpha = 0.6, width = 0.6) +
  labs(x = "", y = "% of genes\nabove SANITY") +
  theme(text = element_text(size = 12),
        legend.position = "none")  +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast")))
ggsave(p2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/DEgene_meanExpr_sanity.pdf"), width = 3, height = 2)

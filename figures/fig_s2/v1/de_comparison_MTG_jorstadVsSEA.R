library(qs)
library(data.table)
library(DESeq2)
library(tidyverse)
library(showtext)
showtext_auto()

#### Subclass

# Set save dir
save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/jorstad_vs_SEA_MTG/")
if(!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
 
# Load data (Jorstad_MTG and SEA)
rawcounts <- list(
  "jorstad" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_allRegions_pseudobulk_by_donor_all_genes/Lein_2023_cell_expression_by_donor_subclass_sum_MTG.csv"), data.table=F) |>
    tibble::column_to_rownames(var = "V1"),
  "sea" = fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/seaad2024_cell_expression_by_donor_subclass_sum_controls.csv")) |>
    tibble::column_to_rownames(var = "V1") #|>
 #   apply(2, ceiling)
)
common_genes <- intersect(rownames(rawcounts$jorstad), rownames(rawcounts$sea))
rawcounts <- lapply(rawcounts, \(x){
  x[rownames(x) %in% common_genes, ]
})
cell_anno_pb <- list(
  "jorstad" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_allRegions_pseudobulk_by_donor_all_genes/Lein_2023_cell_annotations_by_donor_subclass_sum_MTG.csv"),data.table=F) |>
    mutate(label = colnames(rawcounts$jorstad),
           Cell_Type = factor(Cell_Type),
           Donor = factor(Donor)),
  "sea" = fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/seaad2024_cell_annotations_by_donor_subclass_sum_controls.csv")) |>
    mutate(ID = colnames(rawcounts$sea),
           Subclass = factor(Subclass),
           Donor = factor(`Donor ID`),
           Subclass = recode(Subclass, 
                             "Astrocyte" = "Astro", 
                             "Endothelial" = "Endo",
                             "Microglia-PVM" = "Micro/PVM",
                             "Oligodendrocyte" = "Oligo"))
)
cts <- unique(cell_anno_pb[[1]]$Cell_Type)


# Calculate means over all genes for sanity (jorstad mtg)
# sn_expr <- readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_MTG.RDS"))
# # sn_anno <- fread("/mnt/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv", data.table = FALSE)
# # meanbysub <- lapply(unique(sn_anno$Cell_Type), \(subclass){
# #   rowMeans(sn_expr[ ,colnames(sn_expr) %in% sn_anno$Cell_ID[sn_anno$Cell_Type == subclass]])
# # }) |> do.call(what = "cbind")
# # colnames(meanbysub) <- unique(sn_anno$Cell_Type)
# gene_means <- apply(sn_expr, 1, mean)
# gene_means_log <- log(gene_means + 1)
# rm(sn_expr)
# qsave(gene_means_log, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/jorstad_gene_means_log.qs"))
gene_means_log <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/jorstad_gene_means_log.qs"))

# Calculate means over all genes for sanity (gabitto)
# sn_expr <- qread("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/gabitto_2024/expr_UMI_notADsamples.qs")
# gene_means <- apply(sn_expr, 1, mean)
# gene_means_log <- log(gene_means + 1)
# rm(sn_expr)
# qsave(gene_means_log, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/gabitto_gene_means_log.qs"))
gene_means_log_sea <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/gabitto_gene_means_log.qs"))
gene_means_log_sea <- gene_means_log_sea[names(gene_means_log_sea) %in% names(gene_means_log)]
gene_means_log <- gene_means_log[names(gene_means_log) %in% names(gene_means_log_sea)]
gene_means_log_sea <- gene_means_log_sea[match(names(gene_means_log), names(gene_means_log_sea))]

# # Run DEseq 
# # (jorstad mtg)
# reslist <- lapply(cts, \(x){
#   anno_temp <- cell_anno_pb$jorstad |> 
#     mutate("designcol" = ifelse(Cell_Type == x, "ct", "all"))
#   dds <- DESeqDataSetFromMatrix(countData = rawcounts$jorstad,
#                                 colData = anno_temp,
#                                 design= ~ Donor + designcol)
#   dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
#   res <- results(dds1, name = resultsNames(dds1)[2])
#   resdf <- as.data.frame(res@listData) |>
#     `rownames<-`(res@rownames) |>
#    # dplyr::filter(!is.na(log2FoldChange)) |>
#     arrange(padj)
#   return(resdf)
# })
# names(reslist) <- cts
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_mtg_subclass.qs"))
# # SEA
# reslist <- lapply(cts, \(x){
#   anno_temp <- cell_anno_pb$sea |> 
#     mutate("designcol" = factor(ifelse(Subclass == x, "ct", "all")))

#   dds <- DESeqDataSetFromMatrix(countData = rawcounts$sea,
#                                 colData = anno_temp,
#                                 design= ~ Donor + designcol)
#   dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
#   res <- results(dds1, name = resultsNames(dds1)[2])
#   resdf <- as.data.frame(res@listData) |>
#     `rownames<-`(res@rownames) |>
#    # dplyr::filter(!is.na(log2FoldChange)) |>
#     arrange(padj)
#   return(resdf)
# })
# names(reslist) <- cts
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass.qs"))

# # Run DEseq (after filtering genes)
# # (jorstad mtg)
# reslist <- lapply(cts, \(x){
#   anno_temp <- cell_anno_pb$jorstad |> 
#     mutate("designcol" = ifelse(Cell_Type == x, "ct", "all"))
#   dds <- DESeqDataSetFromMatrix(countData = rawcounts$jorstad,
#                                 colData = anno_temp,
#                                 design= ~ Donor + designcol)
#   dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
#   res <- results(dds1, name = resultsNames(dds1)[2])
#   resdf <- as.data.frame(res@listData) |>
#     `rownames<-`(res@rownames) |>
#    # dplyr::filter(!is.na(log2FoldChange)) |>
#     arrange(padj)
#   return(resdf)
# })
# names(reslist) <- cts
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_mtg_subclass_geneSubset.qs"))
# # SEA
# reslist <- lapply(cts, \(x){
#   anno_temp <- cell_anno_pb$sea |> 
#     mutate("designcol" = factor(ifelse(Subclass == x, "ct", "all")))

#   dds <- DESeqDataSetFromMatrix(countData = rawcounts$sea,
#                                 colData = anno_temp,
#                                 design= ~ Donor + designcol)
#   dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
#   res <- results(dds1, name = resultsNames(dds1)[2])
#   resdf <- as.data.frame(res@listData) |>
#     `rownames<-`(res@rownames) |>
#    # dplyr::filter(!is.na(log2FoldChange)) |>
#     arrange(padj)
#   return(resdf)
# })
# names(reslist) <- cts
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass_geneSubset.qs"))

# Reload DE data (jorstad_MTG and sea)
edger <- list(
  "jorstad" = readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/lein_MTG_SEAgenesubset/DE_gene_list_subclass_MTG_blockDonorSum_allgenes.RDS")),
  "sea" = readRDS(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/DE_by_subclass_Controls_jorstadGeneSubset/data/DE_gene_list_subclass_blockDonor_allgenes.RDS"))
) 
# Make ct names match and same order
names(edger[[2]])[c(1, 3, 15, 17)] <- c("Astro", "Endo", "Micro/PVM", "Oligo")
edger[[2]] <- edger[[2]][match(names(edger[[1]]), names(edger[[2]]))]
# Rename celltype columns
edger[[2]] <- mapply(\(x, y){
  out <- x
  out$celltype <- y
  return(out)
}, edger[[2]], names(edger[[2]]), SIMPLIFY = F)

reslist <- list(
  "jorstad" = qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_mtg_subclass_geneSubset.qs")),
  "sea" = qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass_geneSubset.qs"))
) |>
  lapply(\(x){
    mapply(\(y,z){
      y$celltype <- z
      y$genes <- rownames(y)
      colnames(y)[c(2,5,6)] <- c("logFC", "PValue", "FDR")
      return(y)
    }, x, names(x), SIMPLIFY = F)
  }) # add columns for gene and ct

# Gather unique DE genes
uniquegenes <- lapply(list(edger, reslist), \(z){
  a <- lapply(z, \(y){
    do.call(rbind, y) |>
    dplyr::filter(!is.na(PValue)) 
  })
  lapply(a, \(x){
    ug <- list()
    for(i in seq_along(cts)){
      genesin <- x |> 
        dplyr::filter(celltype == cts[i],
                      FDR < 0.05,
                      logFC >= 1) 
      genesout <- x |> 
        dplyr::filter(celltype != cts[i],
                      FDR < 0.05,
                      logFC >= 1)   
      ug[[i]] <- genesin$genes[!genesin$genes %in% genesout$genes]
    } 
    names(ug) <- cts
    return(ug)
  })
}) |>
  set_names("edgeR", "DEseq2")

# Gather DE genes (FDR)
el <- lapply(edger, \(dataset){
  lapply(dataset, \(x) x$genes[x$FDR < 0.05 & x$logFC >= 1])
})
rl <- lapply(reslist, \(dataset){
  lapply(dataset, \(x){
    out <- rownames(x)[x$FDR < 0.05 & x$logFC >= 1] 
    return(out[!is.na(out)])
  })
})

# Plot venn diagram of DE genes (over all subclasses) shared between edgeR and DESeq2
allgenes <- c(unlist(el[[1]]), unlist(rl[[1]])) |> unique()
eulermat <- data.frame("edger" = allgenes %in% unique(unlist(el[[1]])),
                       "DESeq2" = allgenes %in% unique(unlist(rl[[1]])))
# with quantities
pdf(file = file.path(save_dir,"edgeR_vs_deseq2_leinSubclass_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
# no quantitites
pdf(file = file.path(save_dir, "edgeR_vs_deseq2_leinSubclass_euler_noLabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# Plot venn diagram of DE genes (over all subclasses) shared between Jorstad and Gabitto
allgenes <- c(unlist(el[[1]]), unlist(el[[2]])) |> unique()
eulermat <- data.frame("Jorstad" = allgenes %in% unique(unlist(el[[1]])),
                       "Gabitto" = allgenes %in% unique(unlist(el[[2]])))
# with quantities
pdf(file = file.path(save_dir,"jorstadVsGabitto_leinSubclass_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
# no quantitites
pdf(file = file.path(save_dir, "jorstadVsGabitto_leinSubclass_euler_noLabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()


# Plot venn diagram of DE genes (over all subclasses) shared between edgeR and DESeq2
# All comparisons:
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_jorstad" = allgenes %in% unique(unlist(el[[1]])),
                       "edger_sea" = allgenes %in% unique(unlist(el[[2]])),
                       "DESeq2_jorstad" = allgenes %in% unique(unlist(rl[[1]])),
                       "DESeq2_sea" = allgenes %in% unique(unlist(rl[[2]])),
                       "edger_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[1]])),
                       "edger_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[2]])),
                       "DESeq2_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[1]])),
                       "DESeq2_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[2]]))
                       )
# with quantities
pdf(file = file.path(save_dir,"edgeRVsdeseq2_uniqueVsFDR_jorstadVsSEA_Subclass_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
# no quantitites
pdf(file = file.path(save_dir, "edgeRVsdeseq2_uniqueVsFDR_jorstadVsSEA_Subclass_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# jorstad only:
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_jorstad" = allgenes %in% unique(unlist(el[[1]])),
                       "DESeq2_jorstad" = allgenes %in% unique(unlist(rl[[1]])),
                       "edger_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[1]])),
                       "DESeq2_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[1]]))
                       )

pdf(file = file.path(save_dir,"edgeRVsdeseq2_jorstad_Subclass_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_jorstad_Subclass_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# SEA only
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_sea" = allgenes %in% unique(unlist(el[[2]])),
                       "DESeq2_sea" = allgenes %in% unique(unlist(rl[[2]])),
                       "edger_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[2]])),
                       "DESeq2_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[2]]))
                       )
pdf(file = file.path(save_dir,"edgeRVsdeseq2_uniqueVsFDR_SEA_Subclass_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_uniqueVsFDR_SEA_Subclass_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     labels = F)
dev.off()

# FDR only
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_jorstad" = allgenes %in% unique(unlist(el[[1]])),
                       "edger_sea" = allgenes %in% unique(unlist(el[[2]])),
                       "DESeq2_jorstad" = allgenes %in% unique(unlist(rl[[1]])),
                       "DESeq2_sea" = allgenes %in% unique(unlist(rl[[2]]))
                       )
pdf(file = file.path(save_dir,"edgeRVsdeseq2_FDR_jorstadVsSEA_Subclass_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_FDR_jorstadVsSEA_Subclass_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()


# Unique only
allgenes <- unlist(uniquegenes) |> unique()
eulermat <- data.frame("edger_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[1]])),
                       "edger_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[2]])),
                       "DESeq2_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[1]])),
                       "DESeq2_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[2]]))
                       )
pdf(file = file.path(save_dir,"edgeRVsdeseq2_unique_jorstadVsSEA_Subclass_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_unique_jorstadVsSEA_Subclass_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# Plot # of shared and unshared DE genes per subclass (Jorstad, stacked barplot, edger vs DEseq2)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!el[[1]][[ind]] %in% rl[[1]][[ind]]),
                                   sum(!rl[[1]][[ind]] %in% el[[1]][[ind]]),
                                   length(intersect(rl[[1]][[ind]], el[[1]][[ind]]))),
                       "Method" = c("edgeR_LRT", "DESeq2_LRT", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Method <- factor(plotmat$Method, levels = c("edgeR_LRT", "DESeq2_LRT", "Both"))

p <- ggplot(plotmat, aes(x = ct, y = value, fill = Method)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(save_dir, "deseq2_edger_overlap_leinSubclass.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)

# Plot # of shared and unshared DE genes per subclass (Gabitto, stacked barplot, edger vs DEseq2)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!el[[2]][[ind]] %in% rl[[2]][[ind]]),
                                   sum(!rl[[2]][[ind]] %in% el[[2]][[ind]]),
                                   length(intersect(rl[[2]][[ind]], el[[2]][[ind]]))),
                       "Method" = c("edgeR_LRT", "DESeq2_LRT", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Method <- factor(plotmat$Method, levels = c("edgeR_LRT", "DESeq2_LRT", "Both"))
p <- ggplot(plotmat, aes(x = ct, y = value, fill = Method)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")
ggsave(p, file = file.path(save_dir, "deseq2_edger_overlap_gabittoSubclass.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)


# Plot # of shared and unshared DE genes per subclass (jorstad vs sea subclass, edger, stacked barplot)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!el[[1]][[ind]] %in% el[[2]][[ind]]),
                                   sum(!el[[2]][[ind]] %in% el[[1]][[ind]]),
                                   length(intersect(el[[1]][[ind]], el[[2]][[ind]]))),
                       "Study" = c("Jorstad 2023", "Gabitto 2024 con", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Study <- factor(plotmat$Study, levels = c("Jorstad 2023", "Gabitto 2024 con", "Both"))

p <- ggplot(plotmat, aes(x = ct, y = value, fill = Study)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(save_dir, "jorstad_gabitto_overlap_edgerSubclass.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)

# Plot # of shared and unshared DE genes per subclass (jorstad vs sea subclass, deseq2, stacked barplot)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!rl[[1]][[ind]] %in% rl[[2]][[ind]]),
                                   sum(!rl[[2]][[ind]] %in% rl[[1]][[ind]]),
                                   length(intersect(rl[[1]][[ind]], rl[[2]][[ind]]))),
                       "Study" = c("Jorstad 2023", "Gabitto 2024 con", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Study <- factor(plotmat$Study, levels = c("Jorstad 2023", "Gabitto 2024 con", "Both"))
p <- ggplot(plotmat, aes(x = ct, y = value, fill = Study)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")
ggsave(p, file = file.path(save_dir, "jorstad_gabitto_overlap_deseq2Subclass.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)


# Plot # of shared and unshared DE genes per subclass (edgeR vs DEseq2, Jorstad, stacked barplot, UNIQUE)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!uniquegenes$edgeR[[1]][[ind]] %in% uniquegenes$DEseq2[[1]][[ind]]),
                                   sum(!uniquegenes$DEseq2[[1]][[ind]] %in% uniquegenes$edgeR[[1]][[ind]]),
                                   length(intersect(uniquegenes$edgeR[[1]][[ind]], uniquegenes$DEseq2[[1]][[ind]]))),
                       "Method" = c("edgeR_LRT", "DESeq2_LRT", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Method <- factor(plotmat$Method, levels = c("edgeR_LRT", "DESeq2_LRT", "Both"))

p <- ggplot(plotmat, aes(x = ct, y = value, fill = Method)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(save_dir, "deseq2_edger_overlap_leinSubclass_UNIQUE.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)

# Plot # of shared and unshared DE genes per subclass (edgeR vs DEseq2, gabitto, stacked barplot, UNIQUE)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!uniquegenes$edgeR[[2]][[ind]] %in% uniquegenes$DEseq2[[2]][[ind]]),
                                   sum(!uniquegenes$DEseq2[[2]][[ind]] %in% uniquegenes$edgeR[[2]][[ind]]),
                                   length(intersect(uniquegenes$edgeR[[2]][[ind]], uniquegenes$DEseq2[[2]][[ind]]))),
                       "Method" = c("edgeR_LRT", "DESeq2_LRT", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Method <- factor(plotmat$Method, levels = c("edgeR_LRT", "DESeq2_LRT", "Both"))
p <- ggplot(plotmat, aes(x = ct, y = value, fill = Method)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")
ggsave(p, file = file.path(save_dir, "deseq2_edger_overlap_gabittoSubclass_UNIQUE.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)


# Plot # of shared and unshared DE genes per subclass (Jorstad vs Gabitto, stacked barplot, UNIQUE)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!uniquegenes$edgeR[[1]][[ind]] %in% uniquegenes$edger[[2]][[ind]]),
                                   sum(!uniquegenes$edger[[2]][[ind]] %in% uniquegenes$edgeR[[1]][[ind]]),
                                   length(intersect(uniquegenes$edgeR[[1]][[ind]], uniquegenes$edgeR[[2]][[ind]]))),
                       "Study" = c("Jorstad 2023", "Gabitto 2024 con", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Study <- factor(plotmat$Study, levels = c("Jorstad 2023", "Gabitto 2024 con", "Both"))

p <- ggplot(plotmat, aes(x = ct, y = value, fill = Study)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(save_dir, "jorstad_vs_gabitto_overlap_edgerSubclass_UNIQUE.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)

# Plot # of shared and unshared DE genes per subclass (Jorstad vs Gabitto, stacked barplot, UNIQUE, deseq2)
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!uniquegenes$DEseq2[[1]][[ind]] %in% uniquegenes$DEseq2[[2]][[ind]]),
                                   sum(!uniquegenes$DEseq2[[2]][[ind]] %in% uniquegenes$DEseq2[[1]][[ind]]),
                                   length(intersect(uniquegenes$DEseq2[[1]][[ind]], uniquegenes$DEseq2[[2]][[ind]]))),
                       "Study" = c("Jorstad 2023", "Gabitto 2024 con", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Study <- factor(plotmat$Study, levels = c("Jorstad 2023", "Gabitto 2024 con", "Both"))
p <- ggplot(plotmat, aes(x = ct, y = value, fill = Study)) +
  theme_classic() + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  labs(x = "", y = "# of DE genes") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")
ggsave(p, file = file.path(save_dir, "jorstad_vs_gabitto_overlap_deseq2Subclass_UNIQUE.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)


# Boxplot of mean expression values (jorstad vs gabitto)
genes1 <- unlist(el[[1]]) |> unique()
genes2 <- unlist(el[[2]]) |> unique()
allsharedgenes <- edger[[1]][[1]]$genes
shareddegenes <- intersect(genes1, genes2) |> unique()
genecatlist <- list(
  "Not DE" = allsharedgenes[!allsharedgenes %in% shareddegenes],
  "Jorstad DE" = genes1[(genes1 %in% allsharedgenes) & !(genes1 %in% shareddegenes)],
  "Gabitto DE" = genes2[(genes2 %in% allsharedgenes) & !(genes2 %in% shareddegenes)],
  "Both DE" = shareddegenes[shareddegenes %in% allsharedgenes]
)

plotdf_mean_jor <- mapply(\(x, y){
  data.frame("which" = y,
             "means" = gene_means_log[names(gene_means_log) %in% x],
             "type" = "Expression in\nJorstad et al. 2023")
}, genecatlist, names(genecatlist), SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(which = factor(which, levels = unique(which)))
plotdf_mean_sea <- mapply(\(x, y){
  data.frame("which" = y,
             "means" = gene_means_log_sea[names(gene_means_log_sea) %in% x],
             "type" = "Expression in\nGabitto et al. 2024")
}, genecatlist, names(genecatlist), SIMPLIFY = F) |>
  do.call(what = "rbind") |>
  mutate(which = factor(which, levels = unique(which)))
plotdf_mean <- rbind(plotdf_mean_jor, plotdf_mean_sea) |>
  mutate(type = factor(type, levels = c("Expression in\nJorstad et al. 2023", "Expression in\nGabitto et al. 2024")))
p <- ggplot(plotdf_mean, aes(x = which, y = means, fill = which, color = which)) + 
  theme_classic() + 
  geom_jitter(alpha = 0.05) + 
  geom_boxplot(notch = T, color = "black", outlier.shape = NA) +
  geom_hline(yintercept = log(1), color = "red", linetype = "dashed", alpha = 0.8) +
  labs(x = "", y = "log(mean expression)") +
  theme(text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast"))) +
  scale_color_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast"))) +
  facet_wrap(~type, nrow = 1)
ggsave(p, file = file.path(save_dir, "DEgene_meanExpr_boxjitter_jorstadvsGabitto.pdf"), width = 5, height = 2.5)

# Plot pctg of genes above sanity threshold
plotdf_sanity <- lapply(genecatlist, \(x){
  sum(gene_means_log[names(gene_means_log) %in% x] > log(2)) / length(x) * 100
}) |>
  unlist() |> 
  as.data.frame() |>
  tibble::rownames_to_column(var = "which") |>
  mutate(which = factor(which, levels = unique(which)),
         type = "Expression in\nJorstad et al. 2023")
colnames(plotdf_sanity)[2] <- c("val")
plotdf_sanity_sea <- lapply(genecatlist, \(x){
  sum(gene_means_log_sea[names(gene_means_log_sea) %in% x] > log(2)) / length(x) * 100
}) |>
  unlist() |> 
  as.data.frame() |>
  tibble::rownames_to_column(var = "which") |>
  mutate(which = factor(which, levels = unique(which)),
         type = "Expression in\nGabitto et al. 2024")
colnames(plotdf_sanity_sea)[2] <- c("val")
plotdfsan <- rbind(plotdf_sanity, plotdf_sanity_sea) |>
  mutate(type = factor(type, levels = c("Expression in\nJorstad et al. 2023", "Expression in\nGabitto et al. 2024")))
p2 <- ggplot(plotdfsan, aes(x = which, y = val, fill = which)) +
  theme_classic() + 
  geom_bar(stat = "identity", alpha = 0.6, width = 0.6) +
  labs(x = "", y = "% of genes\nabove SANITY") +
  theme(text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  +
  scale_fill_manual(values = c("grey", paletteer::paletteer_d("khroma::highcontrast"))) +
  facet_wrap(~type, nrow = 1)
ggsave(p2, file = file.path(save_dir, "DEgene_meanExpr_sanity.pdf"), width = 4, height = 2)


# Plot tree plot of DE genes shared between methods for each subclass
plotmat <- lapply(cts, \(x){
  ind <- which(cts == x)
  outmat <- data.frame("value" = c(sum(!el[[1]][[ind]] %in% el[[2]][[ind]]),
                                   sum(!el[[2]][[ind]] %in% el[[1]][[ind]]),
                                   length(intersect(el[[1]][[ind]], el[[2]][[ind]]))),
                       "Study" = c("Jorstad 2023", "Gabitto 2024 con", "Both"), 
                       "ct" = x)
  return(outmat)
}) |> do.call(what = "rbind")
ordervec <- plotmat |> group_by(ct) |> summarise(sum = sum(value)) |> arrange(sum)
plotmat$ct <- factor(plotmat$ct, levels = ordervec$ct)
plotmat$Study <- factor(plotmat$Study, levels = c("Jorstad 2023", "Gabitto 2024 con", "Both"))
colvec <- RColorBrewer::brewer.pal(6, "Paired")[c(2,4,6)]
colvec <- alpha(colvec, 0.7)
plotmat2 <- plotmat |>
  dplyr::filter(Study == "Both") |>
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
pdf(file.path(save_dir, "treemap.pdf"))
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
pdf(file = file.path(save_dir, "piechart.pdf"))
pie(plotmat2$value, labels = plotmat2$ct, col = plotmat2$color, border="white", cex = 2)
dev.off()

# Pie chart (ggplot)
plotmat3 <- plotmat2 |>
  group_by(Class) |>
  summarise(value = sum(value)) |>
  mutate(pcnt = value/sum(value) * 100)
p <- ggplot(plotmat3, aes(x = "", y = value, fill = Class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + 
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16)) + 
  scale_fill_manual(values = colvec) 
ggsave(p, file = file.path(save_dir, "piechart_ggplot2.pdf"))

##############
#### Supertype
##############

# Set save dir
save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/jorstad_vs_SEA_MTG_supertype/")
if(!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
   
# Load data (Jorstad_MTG and SEA)
rawcounts <- list(
  "jorstad" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_allRegions_pseudobulk_by_donor_all_genes/Lein_2023_cell_expression_by_donor_supertype_sum_MTG.csv"), data.table=F) |>
    tibble::column_to_rownames(var = "V1"),
  "sea" = fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/seaad2024_cell_expression_by_donor_supertype_sum.csv")) |>
    tibble::column_to_rownames(var = "V1") |>
    apply(2, ceiling)
)
cell_anno_pb <- list(
  "jorstad" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_allRegions_pseudobulk_by_donor_all_genes/Lein_2023_cell_annotations_by_donor_supertype_sum_MTG.csv"),data.table=F) |>
    mutate(label = colnames(rawcounts$jorstad),
           Cluster = factor(Cluster),
           Donor = factor(Donor)),
  "sea" = fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/seaad2024_cell_annotations_by_donor_supertype_sum.csv")) |>
    mutate(Supertype = recode(Supertype,
                              "Lamp5_Lhx6_1" = "Lamp5 Lhx6_1"),
           Supertype = gsub("Micro-PVM", "Micro/PVM", Supertype),
           Supertype = factor(Supertype),
           Donor = factor(`Donor ID`),
           ID = colnames(rawcounts$sea))
)
cts1 <- unique(cell_anno_pb[[1]]$Cluster)
cts2 <- unique(cell_anno_pb[[2]]$Supertype)
cts <- intersect(cts1, cts2)
# There are 91 supertypes shared in common between the two datasets. 
# Subset to these supertypes for subsequent DE analyses

# # Fixes
# # SEA:
# # - Lamp5 Lhx6_1 -> Lamp5_Lhx6_1
# # - Micro-PVM -> Micro/PVM
# # What is overlap of supertype in lein vs sea
# laa <- fread("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv", data.table = FALSE)
saa <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv"))
# lu <- unique(laa$Cluster)
# su <- unique(saa$Supertype)
# # lein table s3
# ts3 <- fread(data.table = F, file = "/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/Table_S3_within_area_cluster_vs_all_markers_by_region_rev.csv")
# tsc <- unique(ts3$cluster)
bla <- saa |> dplyr::filter(`Overall AD neuropathological Change` == "Not AD")


# Calculate means over all genes for sanity (jorstad mtg)
sn_expr <- readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_MTG.RDS"))
sn_anno <- fread(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv"), data.table = FALSE)
meanbysub <- lapply(unique(sn_anno$Cell_Type), \(subclass){
  rowMeans(sn_expr[ ,colnames(sn_expr) %in% sn_anno$Cell_ID[sn_anno$Cell_Type == subclass]])
}) |> do.call(what = "cbind")
colnames(meanbysub) <- unique(sn_anno$Cell_Type)
gene_means <- apply(sn_expr, 1, mean)
gene_means_log <- log(gene_means + 1)
#rm(sn_expr)

# # Run DEseq 
# # (jorstad mtg)
# reslist <- lapply(cts, \(x){
#   anno_temp <- cell_anno_pb$jorstad |> 
#     mutate("designcol" = ifelse(Cluster == x, "ct", "all"))
#   dds <- DESeqDataSetFromMatrix(countData = rawcounts$jorstad,
#                                 colData = anno_temp,
#                                 design= ~ Donor + designcol)
#   dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
#   res <- results(dds1, name = resultsNames(dds1)[2])
#   resdf <- as.data.frame(res@listData) |>
#     `rownames<-`(res@rownames) |>
#     arrange(padj)
#   return(resdf)
# })
# names(reslist) <- cts
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_mtg_supertype.qs"))
# # SEA
# reslist <- lapply(cts, \(x){
#   anno_temp <- cell_anno_pb$sea |> 
#     mutate("designcol" = factor(ifelse(Supertype == x, "ct", "all")))

#   dds <- DESeqDataSetFromMatrix(countData = rawcounts$sea,
#                                 colData = anno_temp,
#                                 design= ~ Donor + designcol)
#   dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
#   res <- results(dds1, name = resultsNames(dds1)[2])
#   resdf <- as.data.frame(res@listData) |>
#     `rownames<-`(res@rownames) |>
#     arrange(padj)
#   return(resdf)
# })
# names(reslist) <- cts
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_supertype.qs"))

# Reload DE data (jorstad_MTG and sea)
edger <- list(
  "jorstad" = readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/DE_gene_list_supertype_MTG_blockDonorSum_allgenes.RDS")),
  "sea" = readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DE_gene_list_supertype_blockDonorSum_allgenes.RDS"))
) 
# Recode edger supertypes and filter to common
names(edger$sea)[names(edger$sea) ==  "Lamp5_Lhx6_1"] <- "Lamp5 Lhx6_1"
names(edger$sea) <- gsub("Micro-PVM", "Micro/PVM", names(edger$sea))
edger <- lapply(edger, \(x){
  temp <- x[names(x) %in% cts]
  temp <- temp[match(cts, names(temp))]
  return(temp)
})

reslist <- list(
  "jorstad" = qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_mtg_supertype.qs")),
  "sea" = qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_supertype.qs"))
) |>
  lapply(\(x){
    mapply(\(y,z){
      y$celltype <- z
      y$genes <- rownames(y)
      colnames(y)[c(5,6)] <- c("PValue", "FDR")
      return(y)
    }, x, names(x), SIMPLIFY = F)
  }) # add columns for gene and ct


# Gather unique DE genes
uniquegenes <- lapply(list(edger, reslist), \(z){
  a <- lapply(z, \(y){
    do.call(rbind, y) |>
    dplyr::filter(!is.na(PValue)) 
  })
  lapply(a, \(x){
    ug <- list()
    for(i in seq_along(cts)){
      genesin <- x |> 
        dplyr::filter(celltype == cts[i],
                      FDR < 0.05) 
      genesout <- x |> 
        dplyr::filter(celltype != cts[i],
                      FDR < 0.05)   
      ug[[i]] <- genesin$genes[!genesin$genes %in% genesout$genes]
    } 
    return(ug)
  })
}) |>
  set_names("edgeR", "DEseq2")

# Gather DE genes (FDR)
el <- lapply(edger, \(dataset){
  lapply(dataset, \(x) x$genes[x$FDR < 0.05])
})
rl <- lapply(reslist, \(dataset){
  lapply(dataset, \(x){
    out <- rownames(x)[x$FDR < 0.05]
    return(out[!is.na(out)])
  })
})


# Plot venn diagram of DE genes (over all subclasses) shared between edgeR and DESeq2
# All comparisons:
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_jorstad" = allgenes %in% unique(unlist(el[[1]])),
                       "edger_sea" = allgenes %in% unique(unlist(el[[2]])),
                       "DESeq2_jorstad" = allgenes %in% unique(unlist(rl[[1]])),
                       "DESeq2_sea" = allgenes %in% unique(unlist(rl[[2]])),
                       "edger_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[1]])),
                       "edger_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[2]])),
                       "DESeq2_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[1]])),
                       "DESeq2_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[2]]))
                       )
# with quantities
pdf(file = file.path(save_dir,"edgeRVsdeseq2_uniqueVsFDR_jorstadVsSEA_Supertype_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
# no quantitites
pdf(file = file.path(save_dir, "edgeRVsdeseq2_uniqueVsFDR_jorstadVsSEA_Supertype_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# jorstad only:
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_jorstad" = allgenes %in% unique(unlist(el[[1]])),
                       "DESeq2_jorstad" = allgenes %in% unique(unlist(rl[[1]])),
                       "edger_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[1]])),
                       "DESeq2_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[1]]))
                       )

pdf(file = file.path(save_dir,"edgeRVsdeseq2_jorstad_Supertype_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_jorstad_Supertype_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# SEA only
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_sea" = allgenes %in% unique(unlist(el[[2]])),
                       "DESeq2_sea" = allgenes %in% unique(unlist(rl[[2]])),
                       "edger_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[2]])),
                       "DESeq2_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[2]]))
                       )
pdf(file = file.path(save_dir,"edgeRVsdeseq2_uniqueVsFDR_SEA_Supertype_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_uniqueVsFDR_SEA_Supertype_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     labels = F)
dev.off()

# FDR only
allgenes <- c(unlist(el), unlist(rl)) |> unique()
eulermat <- data.frame("edger_jorstad" = allgenes %in% unique(unlist(el[[1]])),
                       "edger_sea" = allgenes %in% unique(unlist(el[[2]])),
                       "DESeq2_jorstad" = allgenes %in% unique(unlist(rl[[1]])),
                       "DESeq2_sea" = allgenes %in% unique(unlist(rl[[2]]))
                       )
pdf(file = file.path(save_dir,"edgeRVsdeseq2_FDR_jorstadVsSEA_Supertype_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_FDR_jorstadVsSEA_Supertype_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()


# Unique only
allgenes <- unlist(uniquegenes) |> unique()
eulermat <- data.frame("edger_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[1]])),
                       "edger_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[2]])),
                       "DESeq2_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[1]])),
                       "DESeq2_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[2]]))
                       )
pdf(file = file.path(save_dir,"edgeRVsdeseq2_unique_jorstadVsSEA_Supertype_euler.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "edgeRVsdeseq2_unique_jorstadVsSEA_Supertype_euler_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

# Plot # of shared and unshared DE genes per subclass (Jorstad and SEA separately, stacked barplot)
for(d in seq_along(names(el))){
  plotmat <- lapply(cts, \(x){
    ind <- which(cts == x)
    outmat <- data.frame("value" = c(sum(!el[[d]][[ind]] %in% rl[[d]][[ind]]),
                                    sum(!rl[[d]][[ind]] %in% el[[d]][[ind]]),
                                    length(intersect(rl[[d]][[ind]], el[[d]][[ind]]))),
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
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
    paletteer::scale_fill_paletteer_d("khroma::highcontrast")

  ggsave(p, file = file.path(save_dir, paste0("deseq2_edger_overlap_", names(el)[d], "Supertype.svg")), device = svglite::svglite, bg = "white", width = 14, height = 3)
}

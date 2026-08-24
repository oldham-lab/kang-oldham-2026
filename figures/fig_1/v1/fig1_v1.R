library(qs)
library(data.table)
library(DESeq2)
library(gt)
library(gtExtras)
library(ComplexHeatmap)
library(tidyverse)
library(UpSetR)
library(showtext)
showtext_auto()

##############
# Load data (subclass)
##############

# Set save dir
save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_1/v1/")
if(!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
 
# Load data (Jorstad and SEA) (DFC)
rawcounts <- list(
  "jorstad" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_pseudobulk_by_donor_all_genes/Lein_2023_cell_expression_by_donor_subclass_sum.csv"), data.table=F) |>
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
  "jorstad" = fread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_pseudobulk_by_donor_all_genes/Lein_2023_cell_annotations_by_donor_subclass_sum_DFC.csv"),data.table=F) |>
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


# # Calculate means over all genes for sanity (jorstad mtg)
# sn_expr <- readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS"))
# # sn_anno <- fread(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_MTG.csv"), data.table = FALSE)
# # meanbysub <- lapply(unique(sn_anno$Cell_Type), \(subclass){
# #   rowMeans(sn_expr[ ,colnames(sn_expr) %in% sn_anno$Cell_ID[sn_anno$Cell_Type == subclass]])
# # }) |> do.call(what = "cbind")
# # colnames(meanbysub) <- unique(sn_anno$Cell_Type)
# gene_means <- apply(sn_expr, 1, mean)
# gene_means_log <- log(gene_means + 1)
# rm(sn_expr)
# qsave(gene_means_log, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/jorstad_gene_means_log_dfc.qs"))
gene_means_log <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/jorstad_gene_means_log_dfc.qs"))

# # Calculate means over all genes for sanity (gabitto)
# sn_expr <- qread(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/gabitto_2024/expr_UMI_notADsamples.qs"))
# gene_means <- apply(sn_expr, 1, mean)
# gene_means_log <- log(gene_means + 1)
# rm(sn_expr)
# qsave(gene_means_log, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/gabitto_gene_means_log_dfc.qs"))
gene_means_log_sea <- qread(file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_s2/v1/gabitto_gene_means_log_dfc.qs"))
gene_means_log_sea <- gene_means_log_sea[names(gene_means_log_sea) %in% names(gene_means_log)]
gene_means_log <- gene_means_log[names(gene_means_log) %in% names(gene_means_log_sea)]
gene_means_log_sea <- gene_means_log_sea[match(names(gene_means_log), names(gene_means_log_sea))]

# # Run DEseq (after filtering genes)
# # (jorstad DFC)
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
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_dfc_subclass_geneSubset.qs"))
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
#    dplyr::filter(!is.na(log2FoldChange)) |>
#     arrange(padj)
#   return(resdf)
# })
# names(reslist) <- cts
# qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass_geneSubset_dfc.qs"))

# Reload DE data (jorstad_MTG and sea, after filtering genes)
edger <- list(
  "jorstad" = readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/DE_gene_list_subclass_DFC_blockDonorSum_allgenes.RDS")),
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
  "jorstad" = qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_dfc_subclass_geneSubset.qs")),
  "sea" = qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass_geneSubset_dfc.qs"))
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
                      FDR < 0.05)
                      #logFC >= 1) 
      genesout <- x |> 
        dplyr::filter(celltype != cts[i],
                      FDR < 0.05)
                      #logFC >= 1)   
      ug[[i]] <- genesin$genes[!genesin$genes %in% genesout$genes]
    } 
    names(ug) <- cts
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


###########
# Panel A
###########
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
  labs(x = "", y = "# unique markers") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")

ggsave(p, file = file.path(save_dir, "panel_A.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)


############
# Panel B
###########
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
  labs(x = "", y = "# unique markers") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, margin = margin(0,10,0,0)))  +
  paletteer::scale_fill_paletteer_d("khroma::highcontrast")
ggsave(p, file = file.path(save_dir, "panel_B.svg"), device = svglite::svglite, bg = "white", width = 10, height = 3)


############
# Panel C
###########

allgenes <- unlist(uniquegenes) |> unique()
eulermat <- data.frame("edger_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[1]])),
                       "edger_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[1]][[2]])),
                       "DESeq2_jorstad_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[1]])),
                       "DESeq2_sea_unique" = allgenes %in% unique(unlist(uniquegenes[[2]][[2]]))
                       )
pdf(file = file.path(save_dir,"panel_C.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = list(font = 1, cex = 1),
     quantities = list(cex = 1))
dev.off()
pdf(file = file.path(save_dir, "panel_C_nolabels.pdf"), bg = "white")
plot(eulerr::euler(eulermat),
     #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
     labels = F)
dev.off()

#####
# Panel C redone
#####
# Count overlaps if they match celltypes
genelist <- list(uniquegenes[[1]][[1]],
                 uniquegenes[[1]][[2]],
                 uniquegenes[[2]][[1]],
                 uniquegenes[[2]][[2]])

comparison_names <- c("JE", "GE", "JD", "GD", 
                      "JE&GE", "JE&JD", "JE&GD", "GE&JD", "GE&GD", "JD&GD",
                      "JE&GE&JD", "JE&GE&GD", "JE&JD&GD", "GE&JD&GD",
                      "JE&GE&JD&GD")

comparison_names <- c("Jorstad_edgeR", "Gabitto_edgeR", "Jorstad_DESeq2", "Gabitto_DESeq2", 
                      "Jorstad_edgeR&Gabitto_edgeR", "Jorstad_edgeR&Jorstad_DESeq2", "Jorstad_edgeR&Gabitto_DESeq2", "Gabitto_edgeR&Jorstad_DESeq2", "Gabitto_edgeR&Gabitto_DESeq2", "Jorstad_DESeq2&Gabitto_DESeq2",
                      "Jorstad_edgeR&Gabitto_edgeR&Jorstad_DESeq2", "Jorstad_edgeR&Gabitto_edgeR&Gabitto_DESeq2", "Jorstad_edgeR&Jorstad_DESeq2&Gabitto_DESeq2", "Gabitto_edgeR&Jorstad_DESeq2&Gabitto_DESeq2",
                      "Jorstad_edgeR&Gabitto_edgeR&Jorstad_DESeq2&Gabitto_DESeq2")

comparisons <- list(combn(4,1),
                    combn(4,2),
                    combn(4,3),
                    combn(4,4))

shared_vec <- lapply(comparisons, \(comp){
  outvec <- c()
  for(i in 1:ncol(comp)){
    outvec[i] <- lapply(seq_along(cts), \(j){
      if(length(comp[ ,i]) == 1){
        included <- genelist[[comp[ ,i]]][[j]]
      } else {
        input1 <- lapply(genelist[comp[ ,i]], \(u) u[[j]])
        included <- Reduce(\(a, b){
          intersect(a, b)
        }, input1) 
      }

      input2 <- lapply(genelist[c(1:4)[-comp[ ,i]]], \(u) u[[j]])
      excluded <- Reduce(\(a, b){
        union(a, b)
      }, input2) 

      out <- included[!included %in% excluded]
      return(out)
    }) |>
      unlist() |> unique() |> length()
    #outvec[i] <- uni - nonuni 
  }
  return(outvec)
}) |> unlist()
names(shared_vec) = comparison_names

# pdf(file = file.path(save_dir,"panel_C_new.pdf"), bg = "white")
# plot(eulerr::euler(shared_vec, input = 'disjoint'),
#      #fills = alpha(c(paletteer::paletteer_d("khroma::highcontrast")[1:2], "#0E7175FF"), .6),
#      labels = list(font = 1, cex = 1),
#      quantities = list(cex = 1))
# dev.off()


p <- data.frame("which" = comparison_names, "count" = shared_vec, "highlight" = c(rep('no', 14), 'yes')) |>
  arrange(desc(count)) |>
  mutate(which = factor(which, levels = which)) |>
  ggplot(aes(x = which, y = count, fill = highlight)) +
    theme_minimal() + 
    theme(text = element_text(size = 12),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(1, 0, 0, 0, "cm"),
          axis.line.y = element_line(colour = "black", linewidth = 0.5),
          axis.line.x = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
    geom_bar(stat = 'identity', width = 0.5) +
    geom_text(aes(label = count, color = highlight), vjust = -0.5) +
    labs(x = "", y = "# of intersecting genes") +
    scale_x_discrete(expand = c(0, 1)) +
    scale_y_continuous(limits = c(0, max(shared_vec) + 800), expand = c(0, 1)) +
    scale_fill_manual(values = c('yes' = "orange", 'no' = "black")) +
    scale_color_manual(values = c('yes' = "orange", 'no' = "black")) 
ggsave(p, file = file.path(save_dir,"panel_C_upset_top.svg"), bg = "white", width = 5, height = 2)

svg(file = file.path(save_dir,"panel_C_upset_bottom.svg"), bg = "white", width = 5, height = 4)
upset(fromExpression(shared_vec), 
      order.by = "freq", 
      mb.ratio = c(0.8, 0.2),
      sets.x.label = "# of DE genes",
      queries = list(list(query = intersects, params = list("Jorstad_edgeR", "Gabitto_edgeR", "Jorstad_DESeq2", "Gabitto_DESeq2"), color = "orange", active = T)))
dev.off()

pdf(file = file.path(save_dir,"panel_C_upset_bottom.pdf"), bg = "white", width = 5, height = 4)
upset(fromExpression(shared_vec), 
      order.by = "freq", 
      mb.ratio = c(0.8, 0.2),
      sets.x.label = "# of DE genes",
      queries = list(list(query = intersects, params = list("Jorstad_edgeR", "Gabitto_edgeR", "Jorstad_DESeq2", "Gabitto_DESeq2"), color = "orange", active = T)))
dev.off()
 
#########
# Panel D
# Table with three columns that lists each subclass from top to bottom ranked by # of reproducible marker genes (n=24). Column 1: subclass label; column 2: # markers; column 3: Top 5 markers (unique).
###########

# Collect shared markers
shared_markers <- lapply(cts, \(x){
  ind <- which(cts == x)
  int1 <- intersect(uniquegenes[[1]][[1]][[ind]], uniquegenes[[1]][[2]][[ind]])
  int2 <- intersect(uniquegenes[[2]][[1]][[ind]], uniquegenes[[2]][[2]][[ind]])
  intall <- intersect(int1, int2)
  return(intall)
}) 

# Collect mean ranks of markers per condition
mean_marker_ranks <- lapply(cts, \(x){
  ind <- which(cts == x)
  rlist <- lapply(list(edger[[1]][[ind]],
                       edger[[2]][[ind]],
                       reslist[[1]][[ind]],
                       reslist[[2]][[ind]]), \(y){
    y |> 
      mutate(rank = 1:nrow(y)) |>
      select(genes, rank)
  })
  rcombined <- Reduce(\(a, b){
    left_join(a, b, by = join_by(genes))
  }, rlist) 
  rcombined$mean_ranks <- rowMeans(rcombined[,2:5]) 
  rcombined <- rcombined[order(rcombined$mean_ranks), ]
  return(rcombined |> select(genes, mean_ranks))
})

# Using mean ranks, select top 5 shared DE genes per subclass
top5list <- mapply(\(x, y){
  y[y$genes %in% x, ] |>
    slice(1:5) |>
    pull(genes) |>
    paste(collapse = ", ")
}, shared_markers, mean_marker_ranks, SIMPLIFY = F) |>
  unlist()

# Create table
tab <- data.frame("Subclass" = names(el[[1]]), 
                  "Count" = lapply(shared_markers, length) |> unlist(),
                  "Top_5" = top5list) |>
  arrange(desc(Count))

tab1 <- gt(tab) |>
  cols_label(
    Subclass = html("Subclass"),
    Count = html("#"),
    Top_5 = html("Top 5 DE genes")
  ) |>
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  )
gtsave(tab1, filename = file.path(save_dir, "panel_D.html"))
gtsave(tab1, filename = file.path(save_dir, "panel_D.pdf"))

##########
# Panel E
##########
metacell <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/sn_summary_tables/allmlist_log.qs"))
sn_anno <- fread(file.path(Sys.getenv("ESSA_DIR", "/home/gugene/code/git/ESSA"), "data/Lein_2023_cell_annotations_DFC.csv"), data.table = FALSE)
sn_anno_subclass <- sn_anno[!duplicated(sn_anno$Cell_Type), c(1, 10)]
sn_anno_subclass <- sn_anno_subclass[match(colnames(metacell[[1]]), sn_anno_subclass$Cell_Type), ]

metacell_cor <- lapply(metacell, cor)

# Calculate mean correlations for each class
for(i in unique(sn_anno_subclass$Class)){
  subset <- sn_anno_subclass |>
    filter(Class == i) |>
    pull(Cell_Type)
  subset1 <- colnames(metacell_cor[[1]]) %in% subset
  temp <- metacell_cor[[1]][subset1, subset1]
  cat(i, " ", mean(temp[upper.tri(temp)]), "\n")
}
# Non-neuronal   0.5763327 
# GABAergic   0.9096926 
# Glutamatergic   0.9411046 

# Calculate mean correlation in all neuronal celltypes
subset <- sn_anno_subclass |>
    filter(Class != "Non-neuronal") |>
    pull(Cell_Type)
subset1 <- colnames(metacell_cor[[1]]) %in% subset
temp <- metacell_cor[[1]][subset1, subset1]
cat(i, " ", mean(temp[upper.tri(temp)]), "\n")
# 0.8979705 

# Subclass
cc_colors <- RColorBrewer::brewer.pal(length(unique((sn_anno_subclass$Class))), "Set1")
names(cc_colors) <- unique(sn_anno_subclass$Class)
col_fun = circlize::colorRamp2(c(0, .5, 1), c("blue", "white", "red"))
column_ha = HeatmapAnnotation(`Cell class` = sn_anno_subclass$Class,
                              col = list(`Cell class`=cc_colors),
                              annotation_legend_param = list(
                                `Cell class` = list(
                                  title_gp = gpar(fontface = "plain") 
                                )
                              )
                              )

pdf(file.path(save_dir, "panel_E.pdf"))
h <- Heatmap(metacell_cor[[1]],
             name = "Pearson cor",
             width = unit(0.8, "snpc"), # Adjust width to make cells square
             height = unit(0.8, "snpc"), # Adjust height to make cells square
             col = col_fun, 
             #left_annotation=row_ha, 
             top_annotation = column_ha,
             heatmap_legend_param = list(
               title_gp = gpar(fontface = "plain")
             ),
             row_names_gp = gpar(fontfamily = "sans"),       
             column_names_gp = gpar(fontfamily = "sans")
             #column_names_gp = gpar(fontsize=fontsize)
             )
draw(h, merge_legend=T,
     heatmap_legend_side = "left",
     padding = unit(c(1, 1, 1, 1), "cm")
     )
dev.off()

##########
# Panel F
##########

metacell <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_summary_tables/sn_summary_objects_log.qs"))$mean
sn_anno <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")) %>%
  dplyr::filter(!`Neurotypical reference`)
sn_anno_subclass <- sn_anno |>
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |>
  arrange(match(Subclass, colnames(metacell$Subclass$Control)))

metacell_cor <- cor(metacell$Subclass$Control)
fixed_names <- c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
"L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","Oligo", "OPC", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","Vip","VLMC" )
colnames(metacell_cor) <- fixed_names
rownames(metacell_cor) <- fixed_names

# Calculate mean correlations for each class
for(i in unique(sn_anno_subclass$Class)){
  subset <- sn_anno_subclass |>
    filter(Class == i) |>
    pull(Subclass)
  subset1 <- colnames(metacell_cor) %in% subset
  temp <- metacell_cor[subset1, subset1]
  cat(i, " ", mean(temp[upper.tri(temp)]), "\n")
}
# Non-neuronal   0.6598768 
# GABAergic   0.9225615 
# Glutamatergic   0.9453485 

# Calculate mean correlations for all neuronal celltypes
subset <- sn_anno_subclass |>
    filter(Class != "Non-neuronal") |>
    pull(Subclass)
subset1 <- colnames(metacell_cor) %in% subset
temp <- metacell_cor[subset1, subset1]
cat(i, " ", mean(temp[upper.tri(temp)]), "\n")
# 0.9075513

cc_colors <- RColorBrewer::brewer.pal(length(unique((sn_anno_subclass$Class))), "Set1")
names(cc_colors) <- unique(sn_anno_subclass$Class)
col_fun = circlize::colorRamp2(c(0, .5, 1), c("blue", "white", "red"))
column_ha = HeatmapAnnotation(`Cell class` = sn_anno_subclass$Class,
                              col = list(`Cell class`=cc_colors),
                              annotation_legend_param = list(
                                `Cell class` = list(
                                  title_gp = gpar(fontface = "plain") 
                                )
                              )
                              )

pdf(file.path(save_dir, "panel_F.pdf"))
h <- Heatmap(metacell_cor,
             name = "Pearson cor",
             width = unit(0.8, "snpc"), # Adjust width to make cells square
             height = unit(0.8, "snpc"), # Adjust height to make cells square
             col = col_fun, 
             #left_annotation=row_ha, 
             top_annotation = column_ha,
             heatmap_legend_param = list(
               title_gp = gpar(fontface = "plain")
             ),
             row_names_gp = gpar(fontfamily = "sans"),       
             column_names_gp = gpar(fontfamily = "sans")
             #column_names_gp = gpar(fontsize=fontsize)
             )
draw(h, merge_legend=T,
     heatmap_legend_side = "right",
     padding = unit(c(1, 1, 1, 1), "cm")
     )
dev.off()

##########
# Panel F
##########

metacell <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/sn_summary_tables/sn_summary_objects_log.qs"))$median
sn_anno <- fread(data.table=F,file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")) %>%
  dplyr::filter(!`Neurotypical reference`)
sn_anno_subclass <- sn_anno |>
  select(Subclass, Class) |>
  filter(!duplicated(Subclass)) |>
  mutate(Class = case_match(
    Class, 
    "Neuronal: GABAergic" ~ "GABAergic",
    "Neuronal: Glutamatergic" ~ "Glutamatergic",
    "Non-neuronal and Non-neural" ~ "Non-neuronal"
  )) |>
  arrange(match(Subclass, colnames(metacell$Subclass$Control)))

metacell_cor <- cor(metacell$Subclass$Control)
fixed_names <- c("Astro", "Chandelier", "Endo", "L2/3 IT","L4 IT", "L5 ET", "L5 IT", "L5/6 NP" ,"L6 CT","L6 IT",  
"L6 IT Car3","L6b","Lamp5", "Lamp5 Lhx6", "Micro/PVM","Oligo", "OPC", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl","Vip","VLMC" )
colnames(metacell_cor) <- fixed_names
rownames(metacell_cor) <- fixed_names

# Calculate mean correlations for each class
for(i in unique(sn_anno_subclass$Class)){
  subset <- sn_anno_subclass |>
    filter(Class == i) |>
    pull(Subclass)
  subset1 <- colnames(metacell_cor) %in% subset
  temp <- metacell_cor[subset1, subset1]
  cat(i, " ", mean(temp[upper.tri(temp)]), "\n")
}
# Non-neuronal   0.6598768 
# GABAergic   0.9225615 
# Glutamatergic   0.9453485 

# Calculate mean correlations for all neuronal celltypes
subset <- sn_anno_subclass |>
    filter(Class != "Non-neuronal") |>
    pull(Subclass)
subset1 <- colnames(metacell_cor) %in% subset
temp <- metacell_cor[subset1, subset1]
cat(i, " ", mean(temp[upper.tri(temp)]), "\n")
# 0.9075513

cc_colors <- RColorBrewer::brewer.pal(length(unique((sn_anno_subclass$Class))), "Set1")
names(cc_colors) <- unique(sn_anno_subclass$Class)
col_fun = circlize::colorRamp2(c(0, .5, 1), c("blue", "white", "red"))
column_ha = HeatmapAnnotation(`Cell class` = sn_anno_subclass$Class,
                              col = list(`Cell class`=cc_colors),
                              annotation_legend_param = list(
                                `Cell class` = list(
                                  title_gp = gpar(fontface = "plain") 
                                )
                              )
                              )

pdf(file.path(save_dir, "panel_F_median.pdf"))
h <- Heatmap(metacell_cor,
             name = "Pearson cor",
             width = unit(0.8, "snpc"), # Adjust width to make cells square
             height = unit(0.8, "snpc"), # Adjust height to make cells square
             col = col_fun, 
             #left_annotation=row_ha, 
             top_annotation = column_ha,
             heatmap_legend_param = list(
               title_gp = gpar(fontface = "plain")
             ),
             row_names_gp = gpar(fontfamily = "sans"),       
             column_names_gp = gpar(fontfamily = "sans")
             #column_names_gp = gpar(fontsize=fontsize)
             )
draw(h, merge_legend=T,
     heatmap_legend_side = "right",
     padding = unit(c(1, 1, 1, 1), "cm")
     )
dev.off()

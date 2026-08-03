library(tidyverse)
library(data.table)
library(qs)
library(ComplexHeatmap)


save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s2_3_4_5/")
if(!dir.exists(save_dir)){dir.create(save_dir, recursive = T)}
 
###########
# Load data
###########

# Reload DE data (jorstad_MTG and sea)
edger <- list(
  "jorstad" = readRDS(file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/lein_MTG_SEAgenesubset/DE_gene_list_subclass_MTG_blockDonorSum_allgenes.RDS")),
  "sea" = readRDS(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/DE_by_subclass_Controls_jorstadGeneSubset_mtg/data/DE_gene_list_subclass_blockDonor_allgenes.RDS"))
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
  "sea" = qread(file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass_geneSubset_mtg.qs"))
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
cts <- names(edger[[1]])
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
    names(ug) <- cts
    return(ug)
  })
}) |>
  set_names("edgeR", "DEseq2")

# Create vectors of which subclass is represented by each unique gene
uniqueclasses <- mapply(\(i, i1){
  mapply(\(j, j1){
    mapply(\(x,y){
      if(length(x) > 0){
      data.frame("Subclass" = y, "Genes" = x)
      } else {
        data.frame("Subclass" = NULL, "Genes" = NULL)
      }
    }, j, names(j), SIMPLIFY = F) |> do.call(what = "rbind") |>
      mutate("Dataset" = j1, "Platform" = i1)
  }, i, names(i), SIMPLIFY = F) 
}, uniquegenes, names(uniquegenes), SIMPLIFY = F) 

# Find subclasses of shared genes 
jorstad_shared <- inner_join(uniqueclasses$edgeR$jorstad, uniqueclasses$DEseq2$jorstad, by = join_by(Genes)) |>
  rename(Unique.subclass.edgeR = Subclass.x, Unique.subclass.DEseq2 = Subclass.y) |>
  select(Genes, Unique.subclass.edgeR, Unique.subclass.DEseq2)
gabitto_shared <- inner_join(uniqueclasses$edgeR$sea, uniqueclasses$DEseq2$sea, by = join_by(Genes)) |>
  rename(Unique.subclass.edgeR = Subclass.x, Unique.subclass.DEseq2 = Subclass.y) |>
  select(Genes, Unique.subclass.edgeR, Unique.subclass.DEseq2)



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



#########################
# Table S2
# Jorstad et al. DE genes
#########################

tablist1 <- mapply(\(x, y){
  xout <- x |>
    select(genes, LR, FDR)
  colnames(xout)[2:3] <- paste0(y$celltype[1], c(".LR.edgeR", ".adj.pval.FDR.edgeR"))
  yout <- y |> 
    select(genes, stat, FDR)
  colnames(yout)[2:3] <- paste0(y$celltype[1], c(".LR.DESeq2", ".adj.pval.FDR.DESeq2"))
  out <- left_join(xout, yout, by = join_by(genes))
  out[,c(2,4)] <- apply(out[,c(2,4)], 2, \(x) sprintf("%.3f", x)) # Format likelihood ratios to 3 decimal places
  out[,c(3,5)] <- apply(out[,c(3,5)], 2, \(x) sprintf("%8.3g", x)) # Format p-vals using scientific notation (3 decimals)
  return(out)
}, edger[[1]], reslist[[1]], SIMPLIFY = F)

reduce_func <- function(x, y){
  left_join(x, y, by = join_by(genes))
}

tab <- Reduce(reduce_func, tablist1) |>
  rename(Genes = genes) |>
  mutate(Significant.edgeR = Genes %in% unique(unlist(el[[1]])),
         Significant.DESeq2 = Genes %in% unique(unlist(rl[[1]])),
         Significant.both = Significant.edgeR & Significant.DESeq2,
         Unique.edgeR = Genes %in% unique(unlist(uniquegenes[[1]][[1]])),
         Unique.DESeq2 = Genes %in% unique(unlist(uniquegenes[[2]][[1]])),
         Unique.both = Unique.edgeR & Unique.DESeq2) |>
  left_join(jorstad_shared, by = join_by(Genes)) |>
  arrange(Genes)

fwrite(tab, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s2_3_4_5/table_s2.csv"), na = "NA")
fwrite(tab, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s2_3/table_s2.csv"), na = "NA")

#########################
# Table S3
# Gabitto et al. DE genes
#########################

tablist1 <- mapply(\(x, y){
  xout <- x |>
    select(genes, LR, FDR)
  colnames(xout)[2:3] <- paste0(y$celltype[1], c(".LR.edgeR", ".adj.pval.FDR.edgeR"))
  yout <- y |> 
    select(genes, stat, FDR)
  colnames(yout)[2:3] <- paste0(y$celltype[1], c(".LR.DESeq2", ".adj.pval.FDR.DESeq2"))
  out <- left_join(xout, yout, by = join_by(genes))
  out[,c(2,4)] <- apply(out[,c(2,4)], 2, \(x) sprintf("%.3f", x)) # Format likelihood ratios to 3 decimal places
  out[,c(3,5)] <- apply(out[,c(3,5)], 2, \(x) sprintf("%8.3g", x)) # Format p-vals using scientific notation (3 decimals)
  return(out)
}, edger[[2]], reslist[[2]], SIMPLIFY = F)

reduce_func <- function(x, y){
  left_join(x, y, by = join_by(genes))
}

tab <- Reduce(reduce_func, tablist1) |>
  rename(Genes = genes) |>
  mutate(Significant.edgeR = Genes %in% unique(unlist(el[[2]])),
         Significant.DESeq2 = Genes %in% unique(unlist(rl[[2]])),
         Significant.both = Significant.edgeR & Significant.DESeq2,
         Unique.edgeR = Genes %in% unique(unlist(uniquegenes[[1]][[2]])),
         Unique.DESeq2 = Genes %in% unique(unlist(uniquegenes[[2]][[2]])),
         Unique.both = Unique.edgeR & Unique.DESeq2) |>
  left_join(gabitto_shared, by = join_by(Genes)) |>
  arrange(Genes)

fwrite(tab, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s2_3_4_5/table_s3.csv"), na = "NA")
fwrite(tab, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s2_3/table_s3.csv"), na = "NA")

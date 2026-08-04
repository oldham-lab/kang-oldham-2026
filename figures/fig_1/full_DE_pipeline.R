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


##########
# edgeR
##########
calc_pseudobulk_DE <- function(rawcounts,
                               save_string="", # string to add to end of save file name
                               celltype=NULL, # a vector of factors indicating which celltypes to do DE on
                               blocking=NULL,# blocking variables in a dataframe (donor, region)
                               subset=NULL,# vector indicating indices of samples to subset
                               de_save_dir,
                               verbose = T

){ 
  if(!dir.exists(de_save_dir)){
    dir.create(de_save_dir, recursive = T)
  }
  if(inherits(blocking,"character")){
    blocking <- data.frame("Block"=blocking)
  }
  genes <- rawcounts[,1]
  rawcounts <- rawcounts[,-1]
  
  if(is.null(celltype)){
    celltype <- factor(colnames(rawcounts))
  }
  if(!is.null(subset)){
    rawcounts <- rawcounts[,subset]
    blocking <- as.data.frame(blocking[subset,])
    celltype <- celltype[subset]
  }
  
  unique_cell <- unique(celltype)
  rawcounts <- as.matrix(rawcounts)
  #keep <- filterByExpr(rawcounts, group = celltype, min.count = 5, min.total.count = 10)
  #rawcounts <- rawcounts[keep, ]

  # Calculate DE for each celltype against all other celltypes
  de_list <- list()
  result_list <- list()
  for(i in 1:length(unique_cell)){
    ct_fac <- as.factor(as.numeric(celltype==unique_cell[i]))
    if(!is.null(blocking)){
      design <- model.matrix(~ .,data=cbind(blocking,ct_fac))
    } else {
      design <- model.matrix(~ ct_fac)
    }
 
    y <- DGEList(counts=rawcounts, genes=genes)
    y <- calcNormFactors(y, method="TMM")
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit)
    de <- decideTests.DGELRT(lrt)
    de_list[[i]] <- de
    result_list[[i]] <- as.data.frame(topTags(lrt, n = nrow(rawcounts)))
    result_list[[i]]$celltype <- unique_cell[i]
    if(verbose) 
      cat(i,"out of",length(unique_cell),"done\n")
  }
   
  names(de_list) <- unique_cell
  names(result_list) <- unique_cell
  
  saveRDS(de_list, file = paste0(de_save_dir,"/DE_summary_list_",save_string,".RDS"))
  saveRDS(result_list, file = paste0(de_save_dir,"/DE_gene_list_",save_string,".RDS"))
}

data_save_dir <- file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_allRegions_pseudobulk_by_donor_all_genes/")
base_save_dir <- file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/data/")

# Load and pseudobulk Lein et al. data by donor
regions <- c("A1", 
             "ACC", 
             "AnG", "M1", "MTG", "S1", "V1")

for(i in regions){
    # Load Lein et al. 2023 expression data and cell annotations
    cell_expr <- readRDS(paste0(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_"), i,".RDS"))
    cell_expr <- as.data.frame(as.matrix(cell_expr))
    cell_anno <- fread(paste0(file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_"), i, ".csv"), data.table = FALSE)
    # cell_anno <- cell_anno[match(colnames(cell_expr), cell_anno$Cell_ID),]

    # Create pseudobulk data by individual per Label (celltype pasted with region) from raw data (sum) and save:
    temp <- cbind(as.data.frame(t(as.matrix(cell_expr))),cell_anno[,colnames(cell_anno) %in% c("Cell_Type", "Donor")])
    cell_expr_pb <- temp %>% group_by(Cell_Type, Donor) %>% summarise(across(everything(), sum))
    samp_lab <- paste(cell_expr_pb[[1]], cell_expr_pb[[2]], sep="_")
    rawcounts <- t(cell_expr_pb[,3:ncol(cell_expr_pb)])
    colnames(rawcounts) <- samp_lab
    fwrite(data.frame(rawcounts), file = paste0(data_save_dir,"/Lein_2023_cell_expression_by_donor_subclass_sum_", i, ".csv"), row.names = T)
    rawcounts <- fread(paste0(data_save_dir,"/Lein_2023_cell_expression_by_donor_subclass_sum_", i, ".csv"), data.table=F)
    new_sif <- cell_anno[,colnames(cell_anno) %in% c("Cell_Type", "Donor", "Cluster", "Region")]
    new_sif$unique <- paste0(new_sif$Donor, new_sif$Cell_Type)
    new_sif <- new_sif[!duplicated(new_sif$unique),]
    new_sif$label <- paste0(new_sif$Cell_Type,"_", new_sif$Donor)
    new_sif <- new_sif[match(samp_lab, new_sif$label),]
    fwrite(new_sif[,-5], file=paste0(data_save_dir,"/Lein_2023_cell_annotations_by_donor_subclass_sum_",i, ".csv"))
    cat(i, " ")
}

sea <- fread(data.table = F, file = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "seaad2024_pseudobulk_by_donor_all_genes/cell_expression_by_donor_sum_subclass_Controls.csv")) |>
  tibble::column_to_rownames(var = "Gene") 


for(i in regions){

    # Calculate p-values and summary statistics for all genes
    rawcounts <- fread(paste0(data_save_dir,"/Lein_2023_cell_expression_by_donor_subclass_sum_", i, ".csv"), data.table=F)
    cell_anno_pb <- fread(paste0(data_save_dir,"/Lein_2023_cell_annotations_by_donor_subclass_sum_",i ,".csv"),data.table=F)
    rawcounts <- rawcounts[rawcounts[,1] %in% rownames(sea), ]

    # sum
    calc_pseudobulk_DE(rawcounts,
                    save_string = paste0("subclass_", i, "_blockDonorSum_allgenes"), 
                    celltype = cell_anno_pb$Cell_Type, 
                    blocking = cell_anno_pb[1],
                    subset = which(cell_anno_pb$Region == i),
                    de_save_dir=base_save_dir)


    cat(i, " ")
}


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

# # Run DEseq (after filtering genes)
# (jorstad DFC)
reslist <- lapply(cts, \(x){
  anno_temp <- cell_anno_pb$jorstad |> 
    mutate("designcol" = ifelse(Cell_Type == x, "ct", "all"))
  dds <- DESeqDataSetFromMatrix(countData = rawcounts$jorstad,
                                colData = anno_temp,
                                design= ~ Donor + designcol)
  dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
  #res <- results(dds1, name = resultsNames(dds1)[2])
  res <- results(dds1, name = "designcol_ct_vs_all")
  resdf <- as.data.frame(res@listData) |>
    `rownames<-`(res@rownames) |>
   # dplyr::filter(!is.na(log2FoldChange)) |>
    arrange(padj)
  return(resdf)
})
names(reslist) <- cts
qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/lein2023/DEseq2/DEseq2_lein_dfc_subclass_geneSubset.qs"))
# SEA
reslist <- lapply(cts, \(x){
  anno_temp <- cell_anno_pb$sea |> 
    mutate("designcol" = factor(ifelse(Subclass == x, "ct", "all")))

  dds <- DESeqDataSetFromMatrix(countData = rawcounts$sea,
                                colData = anno_temp,
                                design= ~ Donor + designcol)
  dds1 <- DESeq(dds, test = "LRT", reduced = ~ Donor)
  res <- results(dds1, name = "designcol_ct_vs_all")
  resdf <- as.data.frame(res@listData) |>
    `rownames<-`(res@rownames) |>
   dplyr::filter(!is.na(log2FoldChange)) |>
    arrange(padj)
  return(resdf)
})
names(reslist) <- cts
qsave(reslist, file = file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "14-ESSA_data/seaad2024_con/DEseq2/DEseq2_seaad2024con_subclass_geneSubset_dfc.qs"))

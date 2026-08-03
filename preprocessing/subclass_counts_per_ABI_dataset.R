# Datasets: 
# Bakken 10x/SSv4, Lein, Lein MTG, Miller MTG

library(data.table)
library(tidyverse)

# Load Bakken sif
cluster_labs <- fread(data.table=F,file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/bakken_2019/cluster_labels_tableS2.csv"))
missing <- c("Chandelier","L4 IT","Lamp5 Lhx6","Pax6")
bsif10x <- fread(data.table=F,"~/root_dir/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/10x/NeMO/author_barcode_annotations_10x.csv")
bsif10x <- left_join(bsif10x, cluster_labs[,c(1,4)], by=join_by(Cell_Type==cluster_label))
bsif10x <- c(table(bsif10x$subclass))
bsif10x <- c(bsif10x, rep(0,4))
names(bsif10x)[21:24] <- missing
bsif10x <- bsif10x[order(names(bsif10x))]
bsifv4 <- fread(data.table=F,"~/root_dir/home/shared/scsn.expr_data/human_expr/postnatal/bakken_2019/SSv4/NeMO/author_barcode_annotations_SSv4.csv")
bsifv4 <- left_join(bsifv4, cluster_labs[,c(1,4)],by=join_by(Cell_Type==cluster_label))
bsifv4 <- c(table(bsifv4$subclass))
bsifv4 <- c(bsifv4, rep(0,4))
names(bsifv4)[21:24] <- missing
bsifv4 <- bsifv4[order(names(bsifv4))]

# Load Lein sif
leinsiflist <- list.files("~/root_dir/home/shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/",full.names=T)
leinsiflist <- leinsiflist[grep("author_barcode_annotations",leinsiflist)]
leinsiflist <- lapply(leinsiflist, function(x) fread(x,data.table=F))
leinsiflist <- lapply(leinsiflist, function(x) c(table(x$Cell_Type)))
leinsifreg <- list.files("~/root_dir/home/shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/",full.names=F)
leinsifreg <- leinsifreg[grep("author_barcode_annotations",leinsifreg)]
leinsifreg <- gsub("author_barcode_annotations_","",leinsifreg)
leinsifreg <- gsub(".csv","",leinsifreg)
#leinMTG <- readRDS("~/root_dir/home/shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824638/expression_data.rds")
#leinsifMTG <- leinMTG@meta.data
#fwrite(leinsifMTG,file="~/data/lein_mtg/metadata_from_seurat_obj.csv")
leinsifMTG <- fread(data.table=F,file="~/data/lein_mtg/metadata_from_seurat_obj.csv")
leinsifMTG <- leinsifMTG %>% filter(assay=="10x 3' v3")
leinsifMTG <- c(table(leinsifMTG$Subclass))

# Load Miller MTG
milsif <- fread(data.table=F, file=file.path(Sys.getenv("MEGASET_DIR", "/home/gugene/RNAseq_megaset"), "13.1-sn_cell_pseudobulk/johansen_by_donor/ct_abundance_subclass.csv"))
milsif <- colSums(milsif[,2:ncol(milsif)])
milsif <- milsif[c(1:5,7:8,6,10:12,9,14,13,15:24)]

# create table
outdf <- data.frame("Study"=c("Bakken","Bakken",rep("Lein",9),"Miller"),
                    "Region"= c("M1", "M1", leinsifreg, "MTG", "MTG"),
                    "Platform"=c("10x","SSv4",rep("10x",10)),
                    "PMID"=c(rep(34616062,2),rep(37824655,8),37824638,37824649))
alltabs <- rbind(bsif10x,bsifv4,do.call(rbind,leinsiflist),leinsifMTG,milsif)
colnames(alltabs) <- names(leinsiflist[[1]])
outdf <- cbind(outdf,alltabs)
outdf$Total <- rowSums(outdf[,5:26])

fwrite(outdf,file="~/data/cell_counts_by_subclass_ABI_cortex_SN.csv")

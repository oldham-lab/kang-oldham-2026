# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(data.table)

#dfcmat <- readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS"))
#fwrite(as.matrix(dfcmat),file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/expression_DFC.csv"), row.names=T)
#dfc <- CreateSeuratObject(counts = dfcmat, project = "leindfc", min.cells = 3, min.features = 200)

# Load Seurat object containing DFC data
dfc <- readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/DFC.rds"))
# remove scale.data slot so that SeuratDisk fills the "raw.X" slot of the anndata object with "counts"
# (see https://stackoverflow.com/questions/67028788/r-how-to-truly-remove-an-s4-slot-from-an-s4-object-solution-attached)
slot(dfc@assays$RNA, "scale.data", check=FALSE) <- NULL # i don't think this changes anything 

SaveH5Seurat(dfc, filename = file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/lein_dfc.h5Seurat"))
Convert(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/lein_dfc.h5Seurat"), dest = "h5ad")

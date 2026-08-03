# https://www.choirclustering.com/articles/CHOIR.html
library(scRNAseq)
library(Seurat)
library(CHOIR)

data_matrix <- LaMannoBrainData('mouse-adult')@assays@data$counts
colnames(data_matrix) <- seq(1:ncol(data_matrix))

object <- CreateSeuratObject(data_matrix, 
                             min.features = 100,
                             min.cells = 5)
object <- NormalizeData(object)
#object <- CHOIR(object, 
#                n_cores = 2)
object <- buildTree(object,
                    n_cores = 2)
object <- pruneTree(object, 
                    n_cores = 2)                    
library(Seurat)
library(CHOIR)

cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
object <- CreateSeuratObject(cell_exprall, 
                             min.features = 100,
                             min.cells = 5)

# Run SCTransform as per Jorstad/Lein methods
#object <- NormalizeData(object)
#object <- PercentageFeatureSet(object, pattern = "^MT-", col.name = "percent.mt")
#object <- SCTransform(object, vars.to.regress = "percent.mt", verbose = FALSE)
#For best performance of CHOIR with SCTransform normalization, please provide the unnormalized count matrix and set the 'normalization_method' parameter to 'SCTransform'
options(future.globals.maxSize= 13000*1024^2)
object <- buildTree(object,
                    n_cores = 8,
                    n_var_features=3000,
                    normalization_method = "SCTransform",
                    use_slot="counts")
# took around 2h40min
qsave(object, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/CHOIR/object_SCTransform_3000varfeatures.qs"))

object <- pruneTree(object, 
                    n_cores = 8,
                    normalization_method = "SCTransform")   
#Error in .validInput(normalization_method, "normalization_method", list(object,  : 
#  SCTransform is not currently supported for Seurat v5 objects.



########### Run with normal normalization
cell_exprall <- as.matrix(readRDS(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")))
object <- CreateSeuratObject(cell_exprall, 
                             min.features = 100,
                             min.cells = 5)
object <- NormalizeData(object)

object <- buildTree(object,
                    n_cores = 8,
                    n_var_features=3000)
# ~ 8 hrs
object <- pruneTree(object, 
                    n_cores = 8)   
qsave(object, file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc/CHOIR/object_3000varfeatures_prunetree.qs"))
# ~7h30min
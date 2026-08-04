library(data.table)
library(qs)
library(future.apply)
options(future.globals.maxSize=1e9)
plan(multisession, workers=10)
library(CoPA)

######### Input data (SEA-AD 2024)
# Load bulk megaset dataset
expr <- fread(file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/count_mats/brainseq_samp_filt_SCZ_SampleNetworks/1_04-04-29/brainseq_samp_filt_SCZ_1_171_outliers_removed_geneSymbolsAdded.csv"), data.table=F)
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/")
sn_summary_object_path <- file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/brainSCOPE/")
module_output_dir <- save_dir1
##################################
 
COPA(expr = expr, 
     plot = F,
     save_dir1 = save_dir1,
     sn_summary_object_path = sn_summary_object_path)

COPA_compare(bulk_genes = expr[,2], 
             rand_n = 10000,
             save_dir1 = save_dir1)
 
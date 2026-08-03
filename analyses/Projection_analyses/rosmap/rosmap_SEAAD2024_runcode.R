library(anndata)
library(reticulate)
library(data.table)
library(qs)
library(future.apply)
library(CoPA)

plan(multisession, workers = 8)

######### Input data (Morabito)
bulk_genes <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)[, 2]
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024")
sn_summary_object_path <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024")
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024")
#################################

# Rename and restructure old objects (incorporated into combine_sn_summary_donor_objects.R)
amout <- qread(file.path(sn_summary_object_path, "sn_summary_tables", "allmlist_log.qs"))
anout <- qread(file.path(sn_summary_object_path, "sn_summary_tables", "allnlist_log.qs"))
avout <- qread(file.path(sn_summary_object_path, "sn_summary_tables", "allvarlist_log.qs"))
a1out <- qread(file.path(sn_summary_object_path, "sn_summary_tables", "all_meansvec_log.qs"))
amout1 <- list("Subclass" = list("Alzheimers" = amout[[1]] |> as.data.frame(), 
                                 "Control" = amout[[2]] |> as.data.frame()),
               "Supertype" = list("Alzheimers" = amout[[3]] |> as.data.frame(), 
                                  "Control" = amout[[4]] |> as.data.frame()))
avout1 <- list("Subclass" = list("Alzheimers" = avout[[1]] |> as.data.frame(), 
                                 "Control" = avout[[2]] |> as.data.frame()),
               "Supertype" = list("Alzheimers" = avout[[3]] |> as.data.frame(), 
                                  "Control" = avout[[4]] |> as.data.frame()))
anout1 <- list("Subclass" = list("Alzheimers" = anout[[1]][1, ] |> unlist(), 
                                 "Control" = anout[[2]][1, ] |> unlist()),
               "Supertype" = list("Alzheimers" = anout[[3]][1, ] |> unlist(), 
                                  "Control" = anout[[4]][1, ] |> unlist()))
names(a1out) <- c("Subclass", "Supertype")
outobj <- list("mean"=amout1,
                "var"=avout1,
                "n"=anout1,
                "all_means"=a1out)
qsave(outobj,file=file.path(sn_summary_object_path, "sn_summary_tables","sn_summary_objects_log.qs"))


COPA(save_dir1 = save_dir1,
     plot = F,
     sn_summary_object_path = sn_summary_object_path,
     module_output_dir = module_output_dir)

COPA_compare(bulk_genes = bulk_genes, 
             rand_n = 10000,
             save_dir1 = save_dir1)

filter_under = 3
do_log = T
rand_n = 10000
 
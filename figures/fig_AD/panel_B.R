# Enrichment analysis

# library(topGO)
# library(ALL)
# data(ALL)
# data(geneList)

# affyLib <- paste(annotation(ALL), "db", sep = ".") 
# library(package = affyLib, character.only = TRUE)

# sampleGOdata <- new("topGOdata", 
#                     description = "Simple session", 
#                     ontology = "BP",
#                     allGenes = geneList, 
#                     geneSel = topDiffGenes, 
#                     nodeSize = 10, 
#                     annot = annFUN.db, 
#                     affyLib = affyLib)
# kevin used http://toppgene.cchmc.org/
# shelve this for later

library(data.table)
library(tidyverse)
library(qs)
library(ComplexHeatmap)
options(bitmapType = 'cairo')
source(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "generic_enrichment_function.r"))

### Enrich significantly differentially projected modules using Phenopedia
# Load Phenopedia data
pheno <- fread(data.table=F,file=file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "phenopedia/oct2024/brain_plus_additional.csv")) %>%
  dplyr::filter(disease %in% c("Alzheimer Disease", "Parkinson Disease", "Schizophrenia"),
                n_publications > 5)
pheno_list <- tapply(pheno$gene, pheno$disease, list)
# Missing "Atherosclerosis", "Tauopathies", "Multiple sclerosis", "Schizophrenia", "OCD", "Epilepsy", "ALS", "Autism"

# Load modules
modfdr <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/kme_tables/topmodposbc_table.csv"))
mods <- tapply(modfdr$Gene, modfdr$topmodposfdr, list)

# Load vector of all genes used to produce modules
allgenes <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)[,2]

# Load significant modules
euc_dist <- list("all"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_all.qs")),
                 "pos"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_positive.qs")),
                 "neg"=qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024/euclidean_distances/euclidean_sigmods_negative.qs")))
euc_fdr_neg <- unlist(euc_dist[[3]][[3]])

# Keep track of sig modules that show up for more than one ct
dup_mods <- unique(euc_fdr[duplicated(euc_fdr)])

# Run enrichment
gshg_results <- GSHG_custom(mods[unique(euc_fdr)],
                            pheno_list,
                            allgenes)

gshg_plot <- gshg_results[,-1]
rownames(gshg_plot) <- gshg_results[,1]
gshg_plot <- apply(gshg_plot, 2, function(x) -log10(x))

Heatmap(as.matrix(gshg_plot))

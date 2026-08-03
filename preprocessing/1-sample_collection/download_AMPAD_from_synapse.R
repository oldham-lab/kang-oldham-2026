# Sample selection 

##############
# ROSMAP
##############
# Load sample information file downloaded from synapse (see download_AMPAD_from_synapse.py)
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_biospecimen_metadata.csv"), data.table=F)

# > table(sif$tissue)
# 
# blood                caudate nucleus 
# 2205                           1005                              4 
# cerebellum dorsolateral prefrontal cortex                 frontal cortex 
# 261                           7660                              1 
# frontal pole        Head of caudate nucleus        occipital visual cortex 
# 1                            749                              4 
# posterior cingulate cortex              prefrontal cortex                          serum 
# 739                              8                            597 
# unspecified 
# 1 

# Subset to 7660 dlPFC samples
sif <- sif[sif$tissue == "dorsolateral prefrontal cortex",]



# Do these samples match the ones in the assay metadata file (which includes info such as batch/RIN/strandedness/etc)?
# - Yes but assay metadata also includes miRNA and mass spec specimens for some reason even though the filename indicates rnaSeq
bulk_sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_assay_rnaSeq_metadata.csv"), data.table=F)
# > dim(bulk_sif)
# [1] 3196   13
# > sum(sif[,2] %in% bulk_sif[,1])
# [1] 1746


# > table(sif$assay)
# Biocrates Bile Acids               Biocrates p180                      ChIPSeq 
# 111                          111                          712 
# label free mass spectrometry                    Metabolon             methylationArray 
# 1210                          514                          740 
# mirnaArray                     rnaArray                       rnaSeq 
# 748                          492                         1161 
# scrnaSeq             TMT quantitation               wholeGenomeSeq 
# 995                          400                          466 

# Subset to 1161 rnaSeq samples

sif <- sif[sif$assay == "rnaSeq",]

# > table(sif$nucleicAcidSource)
#           bulk cell sorted cells 
#251          890           20

# Subset to 890 bulk cell samples

sif <- sif[sif$nucleicAcidSource == "bulk cell",]

# Seems like this corresponds to the files available at accession ID syn8612097. There does not seem to be a need to select a subset of files for download. I will continue to download_AMPAD_from_synapse.py and continue the rest of the protocol.

##############
# MSBB
##############

# All files from syn8612191 were downloaded (download_AMPAD_from_synapse.py).
# Keep files listed in sample info file produced in process_sample_info.R
sif <- fread(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_MSBB.csv"), data.table=F)
dl_files <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/MSBB/fastq/"), full.names=T)
dl_names <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/MSBB/fastq/"), full.names=F)
dl_names <- gsub(".accepted_hits.sort.coord.combined.fastq.gz", "", dl_names)
these_files <- grep("fastq.gz", dl_files)
dl_files <- dl_files[these_files]
dl_names <- dl_names[these_files]

#> sum(sif$specimenID %in% dl_names)
#[1] 495
# 495 files are present in downloads

for(i in 1:length(dl_files)){
  if(!dl_names[i] %in% sif$specimenID){
    unlink(dl_files[i])
    cat("deleted i\n")
  }
}

sif <- sif[sif$specimenID %in% dl_names,]
fwrite(sif, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "RNAseq_preprocessing/sample_info_files/sif_MSBB.csv"))




# Process sample information files

library(data.table)

# Brainseq Ph1
# no library batch info appears to be available from either synapse or Jaffe 2018 NatNeuro
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainseq/phenotype_data/phenotypeFile_LIBD_szControl.csv"), data.table=F)
# Filter to >= 18 Age and Control Dx
sif <- sif[sif$Age >= 18,]
sif <- sif[sif$Dx == "Control", ]

fwrite(sif, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_Brainseq.csv"))

# Save SCZ samples 
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainseq/phenotype_data/phenotypeFile_LIBD_szControl.csv"), data.table=F)
# Filter to >= 18 Age and Control Dx
sif <- sif[sif$Age >= 18,]
sif <- sif[sif$Dx == "Schizo", ]

fwrite(sif, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_Brainseq_SCZ.csv"))
# /home/gugene/code/git/Consensus-analysis/Preprocessing/sample_info_files/
# GTEx
# two files: sample.tsv, sequencing.tsv
sif1 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/GTEx/anvil_GTEx/tsvs/sample.tsv"), data.table=F)
# Number of cortical samples: 176+255+209 = 640
# filter to FCX
sif1 <- sif1[sif1$tissue_type_detail %in% c("Brain - Anterior cingulate cortex (BA24)
",
                                            "Brain - Frontal Cortex (BA9)"),]
sif2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/GTEx/anvil_GTEx/tsvs/sequencing.tsv"), data.table=F)
sif2$specimen_id <- gsub(".Aligned.sortedByCoord.out.patched.md.bam_RNASEQ_BAM_FILES", "", sif2$submitter_id)

sif <- left_join(sif1, sif2, by=join_by(specimen_id))
fwrite(sif, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_GTEx.csv"))


# Brainspan
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainspan/SraRunTable_filtered_FCX.csv"))
fwrite(sif, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_Brainspan.csv"))


# NABEC
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/NABEC/SraRunTable.txt"), data.table=F)
fwrite(sif, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_NABEC.csv"))

sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/NABEC/sample_info_dillman_2017.csv"),data.table=F)
#this file only has sequencing batch info, not library batch info

# ROSMAP
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_biospecimen_metadata.csv"), data.table=F)
sif <- sif[sif$tissue %in% c("dorsolateral prefrontal cortex", "frontal cortex", "frontal pole", "prefrontal cortex"),]
sif <- sif[sif$assay=="rnaSeq",]
sif <- sif[sif$nucleicAcidSource=="bulk cell",]
sif2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/ROSMAP/metadata/ROSMAP_assay_rnaSeq_metadata.csv"), data.table=F)
sif_all <- sif %>% inner_join(sif2, by = join_by(specimenID))
# Remove samples with no corresponding individualID
sif_all <- sif_all[!is.na(sif_all[,1]),]
sif_all <- sif_all[sif_all[,1] != "",]
fwrite(sif_all, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_ROSMAP.csv"))

# MSBB
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/MSBB/metadata/MSBB_assay_rnaSeq_metadata.csv"), data.table=F)
sif2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/AMPAD/MSBB/metadata/MSBB_biospecimen_metadata.csv"), data.table=F)
sif2 <- sif2[sif2$specimenID %in% sif$specimenID,]
# limit to frontal cortex
sif2 <- sif2[sif2$BrodmannArea %in% c(10,44),]
sif_all <- sif2 %>% inner_join(sif, by = join_by(specimenID))
# split by region (each individual has frontal pole and inferior frontal gyrus samples)
sif1 <- sif_all[sif_all$tissue=="frontal pole",]
sif2 <- sif_all[sif_all$tissue=="inferior frontal gyrus",]
# remove resequenced samples (both original and resequenced)
remo1 <- sif1$individualID[grep("resequenced", sif1$specimenID)]
sif1 <- sif1[!sif1$individualID %in% remo1,]
remo2 <- sif2$individualID[grep("resequenced", sif2$specimenID)]
sif2 <- sif2[!sif2$individualID %in% remo2,]

fwrite(sif1, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_MSBB_fp.csv"))
fwrite(sif2, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_MSBB_ifg.csv"))


## BrainGVEX
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/BrainGVEX/RNASeq_BrainGVEX_metadata.csv"), data.table=F)
sif2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/BrainGVEX/individual_BrainGVEX_metadata.csv"), data.table=F)
sif_all <- left_join(sif, sif2, by=join_by(individualID))
# limit to DLPFC
sif_all <- sif_all[sif_all$primaryDiagnosis=="control",]
sif_all <- sif_all[!is.na(sif_all[,2]),]
sif_all <- sif_all[sif_all$ageDeath >= 18,]
#sif_all[sif_all[,1] %in% c("2015-1436", "2015-1478"),]
# there are 2 duplicated individuals
fwrite(sif_all, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_BrainGVEX.csv"))
# remove excess files
files <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/RNAseq/BrainGVEX/"), full.names=T)
rem_these <- files
for(i in 1:nrow(sif_all)){
  rem_these <- rem_these[-grep(sif_all[i,2], rem_these)]
}
for(i in 1:length(rem_these)){
  unlink(rem_these[i])
}


# UCLA-ASD
# region info
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/UCLA_ASD/biospecimen_UCLA-ASD_metadata.csv"), data.table=F)
sif <- sif[sif$tissue %in% c("frontal cortex", "frontal lobe", "prefrontal cortex"),]
# DLPFC: BA 8,9,10,46
sif <- sif[sif$BrodmannArea == "BA9",]
# Dx info
sif2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/UCLA_ASD/individual_UCLA-ASD_metadata.csv"), data.table=F)
sif2 <- sif2[sif2$primaryDiagnosis == "control",]
sif <- sif[sif$individualID %in% sif2$individualID,]
sif_temp <- left_join(sif, sif2, by=join_by(individualID))
# RNAseq metadata info
sif3 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/UCLA_ASD/rnaSeq_UCLA-ASD_metadata.csv"), data.table=F)
sif3 <- sif3[sif3$specimenID %in% sif_temp$specimenID,]
sif_all <- left_join(sif3, sif_temp, by=join_by(specimenID))
sif_all <- sif_all[sif_all$ageDeath >= 18,]
fwrite(sif_all, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_UCLA_ASD.csv"))
# Remove filtered out files
files <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/RNAseq/UCLA_ASD/"))
files <- files[-grep("SYNAPSE", files)]
test <- unlist(lapply(sif_all[,1], function(x) strsplit(x, "_")[[1]][1]))
test <- test[!duplicated(test)]
keep_files <- c()
for(i in 1:length(test)){
  if(length(grep(test[i], files))>0){
    keep_files <- c(keep_files,files[grep(test[i], files)])
  }
}
keep_files <- keep_files[-grep("vermis", keep_files)]
keep_files <- keep_files[-grep("ba41", keep_files)]
rem_files <- files[!files %in% keep_files]
home_dir_temp <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/RNAseq/UCLA_ASD/")
for(i in 1:length(rem_files)){
  unlink(paste0(home_dir_temp,rem_files[i]))
}
# not all files present in metadata are present in files
# - of the 41 samples in filtered sample info file, only 29 present in files


# Yale-ASD
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/Yale_ASD/biospecimen_Yale-ASD_metadata.csv"), data.table=F)
sif <- sif[sif$tissue == "dorsolateral prefrontal cortex",]
sif <- sif[sif$nucleicAcidSource != "single nucleus",]
sif2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/Yale_ASD/individual_Yale-ASD_metadata.csv"), data.table=F)
sif2 <- sif2[!duplicated(sif2$individualID),]
sif2 <- sif2[sif2$ageDeathUnits=="Years",]
sif2 <- sif2[sif2$ageDeath >= 18,]
sif2 <- sif2[sif2$primaryDiagnosis!="Autism Spectrum Disorder",]

sif_temp <- inner_join(sif, sif2, by=join_by(individualID))
sif3 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/Yale_ASD/rnaSeq_Yale-ASD_metadata.csv"), data.table=F)
sif3 <- sif3[!duplicated(sif3$specimenID),]
sif_all <- inner_join(sif_temp, sif3, by=join_by(specimenID))
fwrite(sif_all, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_Yale_ASD.csv"))

test <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/RNAseq/Yale_ASD/BrainVarSubset/"))
test <- gsub("PEC_HSB_Yale-UCSF_mRNA_RiboZero_", "", test)
test <- unlist(lapply(test, function(x) strsplit(x, "_")[[1]][1]))
test <- gsub("-R1.fastq.gz", "", test)
test <- gsub("-R2.fastq.gz", "", test)
# there are 13 total samples that have matching names with 


# BipSeq
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/BipSeq/individual_BipSeq_metadata.csv"), data.table=F)
sif <- sif[sif$primaryDiagnosis=="control",]
sif <- sif[sif$ageDeath >= 18,]
sif2 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/BipSeq/biospecimen_BipSeq_metadata.csv"), data.table=F)
sif2 <- sif2[sif2$tissue=="dorsolateral prefrontal cortex",]
sif3 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/PsychENCODE/metadata/BipSeq/rnaSeq_BipSeq_metadata.csv"), data.table=F)
# All of the adult DLPFC samples in this dataset are BP disorder only. No control samples that are adult DLPFC.


#CMC
sif_c <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/metadata/CMC_Human_clinical_metadata.csv"), data.table=F)
# Notable covariates: Dx, Age of Death
sif_c <- sif_c[sif_c$Dx=="Control",]
sif_c$`Age of Death` <- gsub("90[+]", "91", sif_c$`Age of Death`)
sif_c$`Age of Death` <- as.numeric(sif_c$`Age of Death`)
sif_c <- sif_c[sif_c$`Age of Death` >= 18,]

sif_1 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/metadata/CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData_release1.csv"), data.table=F)
sif_3 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/metadata/CMC_Human_rnaSeq_metadata_release3.csv"), data.table=F)
sif_4 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/metadata/CMC_Human_rnaSeq_metadata_release4.csv"), data.table=F)
sif_6 <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/metadata/CMC_Human_rnaSeq_metadata_release6.csv"), data.table=F)

# Release 3: latest release with raw DLPFC RNAseq data
# Use release 3 for DLPFC RNAseq data
sif_3 <- sif_3[sif_3$Individual_ID %in% sif_c$`Individual ID`,]
fwrite(sif_3, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_CMC.csv"))

# remove extra files
files <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/RNAseq/"))
files <- gsub(".accepted_hits.sort.coord.bam.bai", "", files)
files <- gsub(".accepted_hits.sort.coord.bam", "", files)
files_full <- list.files(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/CMC/RNAseq/"), full.names=T)
rem_these <- which(!files %in% sif_3[,2])
for(i in rem_these){
  unlink(files_full[rem_these])
}

########## Create a mega-sif from all datasets

# Load sif files
sif_files <- list.files(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/"),full.names=T)
sif_names <- list.files(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/"))
# Ignore Brainspan dataset for now due to being too much of an outlier
sif_files <- sif_files[-which(sif_names == "sif_Brainspan.csv")]
sif_files <- sif_files[-which(sif_names == "sif_all_datasets.csv")]
sif_names <- sif_names[-which(sif_names == "sif_Brainspan.csv")]
sif_names <- sif_names[-which(sif_names == "sif_all_datasets.csv")]

sif_names <- gsub(".csv", "", sif_names)

sif_list <- list()
for(i in 1:length(sif_files)){
  sif_list[[i]] <- fread(sif_files[i], data.table=F)
}
names(sif_list) <- sif_names

# Which sample covariates are relevant and shared between datasets?
# Brainseq
# - RIN (col 4)
# - Age (col 5)
# - Sex (col 6)
# - Race (col 7)
# - PMI (col 21)
# GTEx
# - rin_number (col 13)
# - collection site (col 8)
# - date_nucleic_acid_isolation (col 45) (batch id, /home/gugene/megaset/cortex_megaset/rnaseq/gtex/gtex_cortex_sif.csv )
# MSBB
# - RIN (col 24)
# - sequencingBatch (col 27)
# - individualID (col 3) (285 individuals for 495 samples)
# - subregion (col 7) (frontal pole, inferior frontal gyrus) (col8 BrodmannArea has same info)
# NABEC
# - sex (col 33)
# ROSMAP
# - individualID (col 1)
# - RIN (col 22)
# - libraryBatch (col 24)
# - sequencingBatch (col 25)

# Make sample id colnames identical
colnames(sif_list[[1]])[1] <- "Sample_ID"
colnames(sif_list[[2]])[which(colnames(sif_list[[2]])=="specimen_id")] <- "Sample_ID"
colnames(sif_list[[3]])[4] <- "Sample_ID"
colnames(sif_list[[4]])[4] <- "Sample_ID"
colnames(sif_list[[5]])[1] <- "Sample_ID"
colnames(sif_list[[6]])[2] <- "Sample_ID"

colnames(sif_list[[2]])[13] <- "RIN"
colnames(sif_list[[2]])[8] <- "Collection_site"
colnames(sif_list[[2]])[45] <- "Library_prep_batch"
colnames(sif_list[[6]])[24] <- "Library_prep_batch"
colnames(sif_list[[5]])[33] <- "Sex"
colnames(sif_list[[3]])[3] <- "Individual_ID"
colnames(sif_list[[4]])[3] <- "Individual_ID"
colnames(sif_list[[6]])[1] <- "Individual_ID"
colnames(sif_list[[3]])[7] <- "Subregion"
colnames(sif_list[[4]])[7] <- "Subregion"


sif_list[[1]]$Study <- "Brainseq"
sif_list[[2]]$Study <- "GTEx"
sif_list[[3]]$Study <- "MSBB_FP"
sif_list[[4]]$Study <- "MSBB_IFG"
sif_list[[5]]$Study <- "NABEC"
sif_list[[6]]$Study <- "ROSMAP"


sif_all <- bind_rows(sif_list)
sif_all <- sif_all[, colnames(sif_all) %in% c("Study", "Sample_ID","Individual_ID", "RIN", "Age", "Sex", "Race", "PMI", "Collection_site", "Library_prep_batch", "Subregion")]
fwrite(sif_all, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "preprocessing/sample_info_files/sif_all_datasets.csv"))




# Download Brainspan data from dbgap

# Accession ID: phs000755.v2.p1
# First, go to SRA run selector (trace.ncbi.nlm.nih.gov) at accession ID
# Click download buttons for metadata and accession list
# Upload metadata and accession list files to /mnt/bdata/gugene/datasets/RNAseq/Brainspan

# Read metadata file:
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainspan/SraRunTable.txt"), data.table=F)

# Filter by assay type
sif <- sif[sif$`Assay Type`=="RNA-Seq",]
# Create a column for region
reg <- unlist(lapply(strsplit(sif$biospecimen_repository_sample_id, "_"), function(x) x[2]))
sif$region <- reg

# How many cortical samples?
#sum(sif$region %in% c("A1C", "DFC", "IPC", "ITC", ))

# Filter to frontal cortex
sif <- sif[sif$region == "DFC",]
# Save filtered metadata file
fwrite(sif, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainspan/SraRunTable_filtered_FCX.csv"))
              
# Load accession list and filter to files in sif

acc <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainspan/SRR_Acc_List.txt"), data.table=F, header=F)
#> sum(acc[,1] %in% sif[,1])
#[1] 40
acc <- acc[acc[,1] %in% sif[,1],]
acc <- as.data.frame(acc)
fwrite(acc, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/Brainspan/SRR_Acc_List_filtered.txt"), col.names=F)

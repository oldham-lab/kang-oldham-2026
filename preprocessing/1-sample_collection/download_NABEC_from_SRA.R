# Download NABEC data from dbgap

# Accession ID: phs001300.v4.p1
# First, go to SRA run selector (trace.ncbi.nlm.nih.gov) at accession ID
# Click download buttons for metadata and accession list
# Upload metadata and accession list files to /mnt/bdata/gugene/datasets/RNAseq/NABEC

# Read metadata file:
sif <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/NABEC/SraRunTable.txt"), data.table=F)
# All samples are RNAseq frontal cortex

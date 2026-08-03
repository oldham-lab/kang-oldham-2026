library(tidyverse)
library(qs)
library(data.table)


# Load Lein M1 metadata
leinm1 <- fread(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/author_barcode_annotations_M1.csv"), data.table=F)
# > table(leinm1$Donor)
# H18.30.001 H18.30.002 H19.30.001 H19.30.002 
#      49840      39021      13068      12676

# Load Bakken M1 metadata
bakm1 <- fread(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/bakken_2019/10x/ABI/sampleinfo_barcode_level.csv"), data.table=F)
# table(bakm1$external_donor_name_label)
# H18.30.001 H18.30.002 
#      42728      33805

sum(leinm1$Cell_ID %in% bakm1[,2]) 

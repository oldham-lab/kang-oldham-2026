library(tidyverse)
library(qs)
library(data.table)

# Summary of all bulk gene expression datasets analyzed in this study.  
# Attributes: First author last name, Last author last name, 
# Journal, Year Published, PMID, Data Repository, 
# Accession ID, Platform, Number of Individuals, Number of Samples

# ROSMAP, CMC, CMC_HBCC, BrainGVEX, Brainseq, GTEx, NABEC

t5 <- data.frame(
  "Dataset" = c("ROSMAP", "CMC", "CMC_HBCC", "BrainGVEX", "Brainseq", "GTEx", "NABEC"),
  "First_author_last_name" = c("Mostafavi", "Hoffman", "Hoffman", "Wang", "Jaffe", "Aguet", "Dillman"),
  "Last_author_last_name" = c("De Jager", "Roussos", "Roussos", "Gerstein", "Weinberger", "Ardlie", "Cookson"),
  "Journal" = c("Nature Neuroscience", "Scientific Data", "Scientific Data", "Science", "Nature Neuroscience", "Science", "Scientific Reports"),
  "Year" = c(2018, 2019, 2019, 2018, 2018, 2020, 2017),
  "PMID" = c(29802388, 31551426, 31551426, 30545857, 30050107, 32913098, 29203886),
  "Data_repository" = c("syn3219045", "syn2759792", "syn2759792", "syn3270015", "syn12299750", "https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression", "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001353.v1.p1"),
  "Platform" = c("RNAseq", "RNAseq", "RNAseq", "RNAseq", "RNAseq", "RNAseq", "RNAseq"),
  "No_of_individuals" = c(478, 285, 216, 255, 182, 190, 63),
  "No_of_samples" = c(597, 274, 145, 255, 182, 190, 63)
)

fwrite(t5, file = file.path(Sys.getenv("EXTFIG_DIR", "/home/gugene/figures"), "table_s5/table_s5.csv"))

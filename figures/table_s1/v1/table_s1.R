library(tidyverse)
library(qs)
library(data.table)
library(gt)
library(gtExtras)
library(viridis)
library(ggcharts)
library(showtext)
showtext_auto()
options(bitmapType = 'cairo')

# Resolve the directory containing this script (works via Rscript, source(), and RStudio)
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}
SAVE_DIR <- get_script_dir()



# table 1: number of nuclei by donor, diagnosis, cortical region, and technology
# For all datasets used in the study
# - Jorstad 2023 Cv3
# - Jorstad 2023 SSv4
# - Gabitto 2024 Cv3 (SEAAD2024)
# - Liu 2025 Cv3 (MIT_Multiome_Multiregion)
# - Emani 2024 Cv3 (brainSCOPE)
# - Morabito 2021 Cv3

#  This table should list each dataset in the same fashion as the current Fig. S1a 
#  (i.e., with number of nuclei broken down by technology and brain region) 
#  but adding a column next to donor to indicate disease status.

# Columns:
# Dataset, Technology(Cv3 vs SSv4), Donor, Disease status, 1 column for each region

# Jorstad (Cv3)
jor_cv3_counts <- list.files(file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x"), full.names = T) |>
  grep(pattern = "author", value = T) |>
  lapply(\(i){
    fread(data.table = F, file = i) |>
      group_by(Donor, Region) |>
      summarise(n = n()) |>
      mutate(Platform = "Cv3", 
             Pathology = "Control")
  }) |> do.call(what = "rbind") |>
  mutate(#Region = case_match(Region, "DFC" ~ "DLPFC", .default = Region),
         Dataset = "Jorstad et al. 2023") |>
  dplyr::filter(Region %in% c("DFC", "MTG", "V1"))

# Jorstad (SSv4)
jor_ssv4_counts <- readRDS(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_ssv4/data/meta_ss_3_29_22.RDS")) |>
  group_by(donor, region) |>
  summarise(n = n()) |>
  mutate(Platform="SSv4", 
         Pathology = "Control",
         Dataset = "Jorstad et al. 2023") |>
  rename("Donor" = donor, "Region" = region) |>
  dplyr::filter(Region %in% c("DLPFC", "MTG")) |>
  mutate(Region = case_match(Region, "DLPFC" ~ "DFC", .default = Region))

# Gabitto (DLPFC)
gab_counts_dlpfc <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  mutate(Pathology = case_when(
    `Overall AD neuropathological Change` == "Not AD" ~ "Control",
    TRUE ~ "AD"
  )) |>
  group_by(`Donor ID`, Pathology) |>
  summarise(n = n()) |>
  mutate(Platform = "Cv3", 
         Region = "DFC",
         Dataset = "Gabitto et al. 2024") |>
  rename("Donor" = `Donor ID`)

# Gabitto MTG
gab_counts_mtg <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")) %>%
  dplyr::filter(!`Neurotypical reference`) |>
  mutate(Pathology = case_when(
    `Overall AD neuropathological Change` == "Not AD" ~ "Control",
    TRUE ~ "AD"
  )) |>
  group_by(`Donor ID`, Pathology) |>
  summarise(n = n()) |>
  mutate(Platform = "Cv3", 
         Region = "MTG",
         Dataset = "Gabitto et al. 2024") |>
  rename("Donor" = `Donor ID`)

# Liu (MTC, DLPFC)
liu_counts <- fread(data.table = F, file = "/mnt/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv")  |>
  mutate(Pathology = case_when(
    Pathology == "nonAD" ~ "Control",
    TRUE ~ "AD"
  )) |>
  group_by(ROSMAP_IndividualID, Pathology, BrainRegion) |>
  summarise(n = n()) |>
  mutate(Platform = "Cv3",
         Dataset = "Liu et al. 2025") |>
  rename("Donor" = ROSMAP_IndividualID,
         "Region" = BrainRegion) |>
  filter(Region %in% c("MTC", "PFC")) |>
  mutate(Region = case_match(Region,
                             "PFC" ~ "DFC", 
                             "MTC" ~ "MTG",
                             .default = Region))

# Emani (DLPFC)
donorobjs <- list.files("/mnt/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/sn_summary_tables/by_donor", full.names=T)
emani_sn_anno <- fread(data.table=F, file="/mnt/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/PEC2_sample_metadata.txt") |>
    dplyr::filter(Disorder %in% c("control", "Schizophrenia"),
                  Cohort %in% c("CMC", "SZBDMulti-Seq")) |>
    dplyr::filter(lapply(Individual_ID, \(x) sum(grepl(x, donorobjs)) > 0) |> unlist()) 
# 100 CMC donors, 53 control 47 SCZ, 287783 Control, 214234 Schizophrenia
# 48 SZBDMulti-Seq donors, 24 control, 24 SCZ, 173570 control, 171313 SCZ

emani_counts <- list.files("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/Cell_metadata", full.names = T) |>
  lapply(\(x) fread(data.table = F, file = x)) |>
  lapply(\(x) c(table(x$individualID)) |> stack()) |>
  do.call(what = "rbind") |>
  rename(n = values, Individual_ID = ind) |>
  inner_join(emani_sn_anno, by = join_by(Individual_ID)) |>
  mutate(Pathology = case_match(Disorder, "control" ~ "Control", .default = Disorder)) |>
  select(n, Individual_ID, Pathology) |>
  rename(Donor = Individual_ID) |>
  mutate(Region = "DFC", 
         Platform = "Cv3",
         Dataset = "Emani et al. 2024")

# Morabito (DLPFC)
mor_counts <- fread(data.table=F,file=file.path(Sys.getenv("SHARED_DIR", "/home/shared"), "scsn.expr_data/human_expr/postnatal/morabito_2021/all_cells/metadata.csv")) |>
  group_by(SampleID, Diagnosis) |>
  summarise(n = n()) |>
  mutate(Region = "DFC", 
         Platform = "Cv3",
         Dataset = "Morabito et al. 2021") |>
  rename(Donor = SampleID, Pathology = Diagnosis)

# Combine counts
counts_all <- dplyr::bind_rows(jor_cv3_counts,
                               jor_ssv4_counts,
                               gab_counts_dlpfc,
                               gab_counts_mtg,
                               liu_counts,
                               emani_counts,
                               mor_counts) |>
  select(Dataset, Donor, Platform, Pathology, Region, n) |>
  pivot_wider(names_from = Region, values_from = n) |>
  as.data.frame() 
fwrite(as.data.frame(counts_all), file = file.path(SAVE_DIR, "table_s1.csv"))

###########
# Table 1b
###########

## Counting nuclei per case/control
nuc_sum <- rowSums(counts_all[,-c(1:4)], na.rm=T)
nuc_sum2 <- counts_all |> mutate(nuc_sum = nuc_sum) |> group_by(Dataset, Pathology, Platform) |> 
  summarise(nucsum = sum(nuc_sum))

## Determining median UMI per nuc
# Jorstad Cv3
jormtg <- readRDS("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")
jorm1 <- apply(jormtg, 2, sum)
rm(jormtg)
jordfc <- readRDS("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_DFC.RDS")
jord1 <- apply(jordfc, 2, sum)
rm(jordfc)
jorv1 <- readRDS("/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/jorstad_2023_PMID_37824655/10x/expression_V1.RDS")
jorv11 <- apply(jorv1, 2, sum)
rm(jorv1)
median(c(jorm1, jord1, jorv11))
# Jorstad SSv4
jormtg <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_mtg_ssv4/lein_mtg_counts.qs"))
jorm1 <- apply(jormtg, 2, sum)
rm(jormtg)
jordfc <- qread(file.path(Sys.getenv("DATA_OTHER_DIR", "/home/gugene/data_other"), "lein_dfc_ssv4/lein_dfc_counts.qs"))
jord1 <- apply(jordfc, 2, sum)
rm(jordfc)
median(c(jorm1, jord1))
# Gabitto
gdfc <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")) |>
  dplyr::filter("Overall AD neuropathological Change" != "Reference")
gmtg <- fread(data.table=F,file=file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")) |>
  dplyr::filter("Overall AD neuropathological Change" != "Reference")
median(c(gdfc$`Number of UMIs`, gmtg$`Number of UMIs`))
rm(gdfc)
rm(gmtg)
# Liu
liu <- fread(data.table=F,file="/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv") |>
  dplyr::filter(BrainRegion %in% c("MTC", "PFC"))
median(liu$total_counts)
rm(liu)
# Emani
cmc <- fread(data.table=F,file="/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/Cell_metadata/CMC_cell_metadata_redo.tsv")
szbd <- fread(data.table=F,file="/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/brainSCOPE/Cell_metadata/SZBD-Kellis_cell_metadata_redo.tsv")
#median(c(cmc$n_counts, szbd$n_counts))
median(cmc$n_counts)
median(szbd$n_counts)

# Morabito
mor <- Matrix::readMM(file="/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/morabito_2021/all_cells/matrix.mtx")
mor1 <- apply(mor,2,sum)
median(mor1)
rm(mor)

tab1b <- data.frame(
  "Dataset" = c("Jorstad et al. 2023", "Jorstad et al. 2023", "Gabitto et al. 2024", "Liu et al. 2025", "Emani et al. 2024 (CMC)","Emani et al. 2024 (SZBDMulti-Seq)", "Morabito et al. 2021"),
  "Platform" = c("Cv3", "SSv4", "Cv3", "Cv3", "Cv3", "Cv3", "Cv3"),
  "# Cases" = c(0, 0, 75, 51, 47, 24, 11),
  "# Controls" = c(3, 3, 9, 55, 53, 24, 7),
  "# Nuclei (Cases)" = c(0, 0, 2224098, 527907, 214234, 171313, 38676),
  "# Nuclei (Controls)" = c(347152, 21858, 321832, 626175, 287783, 173570, 22796),
  "Median # UMIs per Nucleus" = c(16082, 1710763, 18921, 5723, 2846, 9997, 6387),
  "Pubmed_ID" = c(37824655, 37824655, 39402379, 40752494, 38781369, 38781369, 34239132),
  "Data_repository" = c("NeMO Archive", "NeMO Archive", "Synapse (AD Knowledge Portal)",
                        "Synapse (AD Knowledge Portal)", "Synapse (PsychENCODE)",
                        "Synapse (PsychENCODE)", "GEO"),
  # Accessions verified against each paper's data availability statement:
  #  - Jorstad: all BICCN raw data at NeMO dat-rg2rc5m; the MTG SMART-seq v4 subset
  #    additionally has dat-swzf4kc (also portal.brain-map.org/.../human-mtg-smartseq).
  #  - Liu: the multiregion snRNA-seq used here is the Mathys et al. deposition syn52293442
  #    (Liu's own syn66271521/22 are this study's snATAC/snMultiome, not the snRNA matrix).
  #  - Emani: CMC & SZBDMulti-Seq are PsychENCODE cohorts -> PEC Capstone II Cross-study
  #    Harmonized Data, syn51111084 (the paper gives no per-cohort syn IDs).
  "Accession_ID" = c("nemo:dat-rg2rc5m", "nemo:dat-rg2rc5m; nemo:dat-swzf4kc",
                     "syn26223298", "syn52293442", "syn51111084",
                     "syn51111084", "GSE174367")
)
fwrite(tab1b, file = file.path(SAVE_DIR, "table_s1b.csv"))


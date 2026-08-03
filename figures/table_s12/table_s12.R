library(qs)
library(data.table)
library(tidyverse)

version <- "v1"
save_dir <- file.path(file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/table_s12/"), version)
if(!dir.exists(save_dir))
  dir.create(save_dir, recursive = T)

#####
# Table S12
# List of all seed, tomodposBC genes assigned to all modules, ranked by kME
# ROSMAP AD modules
#######
# Load module data for filtering
datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
filter_under <- 3
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
these_mods_final <- these_mods # size filter only (ROSMAP AD modules are not sigcount-filtered in fig 6/7)
these_mods_final_df <- data.frame("old" = these_mods_final, "new" = 1:length(these_mods_final))

seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
seeddf <- mapply(\(x, y){
  data.frame("seed" = x, "Gene" = y)
}, 1:length(seed), seed, SIMPLIFY = F) |>
  do.call(what = "rbind")
kme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv")) |>
  left_join(seeddf, by = join_by(Gene)) |>
  select(all_of(c("ensembl_id", "Gene", "seed", "topmodposbc")))

kme1 <- kme
kme1$seed <- left_join(kme1[3], these_mods_final_df, by = join_by(seed == old)) |> pull(new)
kme1$topmodposbc <- left_join(kme1[4], these_mods_final_df, by = join_by(topmodposbc == old)) |> pull(new)

fwrite(kme1, file = file.path(save_dir, "table_s12.csv"))

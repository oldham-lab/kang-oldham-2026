library(data.table)
library(jsonlite)

# Bulk megaset modules (topmodposbc)
mods <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/topmodposbc_table.csv"))
mods2 <- tapply(mods$Gene, mods$topmodposbc, list)
write_json(mods2, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/topmodposbc_table.json"))

# Bulk megaset modules (seed)
seed <- qs::qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
jsonlite::write_json(seed, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/LeinDFC/kme_tables/seed_table.json"))

# ROSMAP AD modules (topmodposbc)
mods <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv"))
mods2 <- tapply(mods$Gene, mods$topmodposbc, list)
jsonlite::write_json(mods2, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.json"))

# Brainseq SCZ modules (topmodposbc)
mods <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/kme_tables/topmodposbc_table.csv"))
mods2 <- tapply(mods$Gene, mods$topmodposbc, list)
jsonlite::write_json(mods2, file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/Brainseq_SCZ/kme_tables/topmodposbc_table.json"))

# precalculate CoPA GSEA graph inputs (projection_objects.R/create_gsea_obj)
# - for creating gsea graph on the fly in shiny app

library(data.table)
library(qs)

expr <- fread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/RNAseq/combined_mats/combined_FCX_final_SampleNetworks/1_10-35-00/combined_FCX_final_1_1518_ComBat.csv"), data.table=F)
save_dir1 <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/")
sn_summary_object_path <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/")
module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
do_log = T
logstr <- ifelse(do_log, "_log", "")
filter_under = 3


broad = T


save_dir <- file.path(save_dir1, "sn_proj_objects")

filter_under <- 3
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods_1023 <- as.numeric(names(mods)[which(modulelengths>filter_under)])
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
ind_1016 <- which(!these_mods_1023 %in% which(sigcount_bonf$vals < 2))
these_mods <- these_mods_1023[ind_1016]
mod_seed <- qread(file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC/modules/unmerged_modules.qs"))
modulelengths_seed <- unlist(lapply(mod_seed, length))
mod_id <- 1:1158

# read GSEA data
gsea <- fread(file.path(module_output_dir,"GSEA","gsea_oldham.csv"), data.table=F)
gsea_sn <- fread(file.path(module_output_dir,"GSEA","gsea_SN.csv"), data.table=F)
if(broad){
gsea_broad <- fread(file.path(module_output_dir,"GSEA","gsea_Broad.csv"), data.table=F)
broad_overview <- fread(file.path("data", "broad_geneset_overview.csv"), data.table = F) |>
    dplyr::filter(!catType %in% c("c1", "c3", "c4", "c6", "c7")) |>
    dplyr::filter(!grepl("^HP_", setNames))
gsea_broad <- gsea_broad[gsea_broad$SetID %in% broad_overview$setIDs, ]
gsea_comb <- rbind(gsea, gsea_sn, gsea_broad)
} else {
gsea_comb <- rbind(gsea, gsea_sn)
}

# calculate fdr cutoff
gseavec <- unlist(list(gsea_comb[,4:ncol(gsea_comb)]))
gsq <- qvalue::qvalue(gseavec)
gsea_cutFDR <- max(gsq$pvalues[gsq$qvalues<.05])

# -log10 transform gsea p-values
for(k in 4:ncol(gsea_comb)){
gsea_comb[,k] <- -log10(gsea_comb[,k])
gsea_comb[gsea_comb[,k]==Inf,k] <- 308
}
gsea_cutBC <- -log10(0.05/(nrow(gsea_comb) * (ncol(gsea_comb)-3)))
gsea_cutFDR <- -log10(gsea_cutFDR)

# Add breaks to geneset names (for graphing)
# insert_char_every_n <- function(string, char_to_insert = "\n", n = 50) {
#   gsub(paste0("(.{", n, "})"), paste0("\\1", char_to_insert), string)
# }
insert_char_after_underscore <- function(s, char_to_insert = "\n", n = 50) {
    library(stringr)
    # Find all underscore positions
    underscore_positions <- str_locate_all(s, "_")[[1]][, "start"]

    # Find the first underscore position after n characters
    insertion_point <- NA
    for (pos in underscore_positions) {
        if (pos > n) {
        insertion_point <- pos
        break
        }
    }

    if (is.na(insertion_point)) {
        # No underscore found after n characters, return original string
        return(s)
    } else {
        # Split the string and insert the character
        part1 <- str_sub(s, 1, insertion_point)
        part2 <- str_sub(s, insertion_point + 1, str_length(s))
        return(paste0(part1, char_to_insert, part2))
    }
}
# new_names <- unlist(lapply(gsea_comb$SetName, insert_char_after_underscore))
# new_names <- unlist(lapply(new_names, \(x) insert_char_after_underscore(x, n = 100)))
# new_names <- unlist(lapply(new_names, \(x) insert_char_after_underscore(x, n = 150)))
new_names <- unlist(lapply(gsea_comb$SetName, \(x){
if(nchar(x) >= 100){
    return(insert_char_after_underscore(x, n = nchar(x) / 2))
} else {
    return(x)
}
}))
gsea_comb$SetName_breaks <- new_names
gsea_comb <- gsea_comb |> dplyr::relocate(SetName_breaks, .after = SetName)

# get lowest gsea p-value to set upper limit for bar graph
top_gsea <- max(apply(gsea[4:ncol(gsea)], 2, max))
top_gsea_sn <- max(apply(gsea_sn[4:ncol(gsea_sn)], 2, max))


gsea_plot_dfs <- list()
for(i in seq_along(these_mods)){
    j <- these_mods[i]
    ## enrichment featuring most significant genesets from our lab collection + sn genesets + broad
    # order gsea data by p-value
    gsea_ind <- which(colnames(gsea_comb) == paste0("X", j))
    if(length(gsea_ind)>0){
        gsea_plot <- gsea_comb[order(gsea_comb[,gsea_ind], decreasing=T),]
        gsea_plot <- gsea_plot[1:10, c(1:4, gsea_ind)]
        colnames(gsea_plot)[5] <- "pval"
        gsea_plot$SetName <- factor(gsea_plot$SetName,levels=unique(gsea_plot$SetName))
        gsea_plot$SetName_breaks <- factor(gsea_plot$SetName_breaks,levels=unique(gsea_plot$SetName_breaks))
        gsea_plot_dfs[[i]] <- gsea_plot

    #   cat(j, " ")
    }
}

#fwrite(gsea_comb, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/shiny_code/gsea_comb.csv"))
library(qs)
qsave(gsea_plot_dfs, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/shiny_code/gsea_plot_dfs.qs"))

# gsea_comb csv table is 53M, gsea_plot_dfs.qs is 269k

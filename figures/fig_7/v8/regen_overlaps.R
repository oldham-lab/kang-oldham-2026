# Standalone (re)generation of dfc_overlaps.csv / mtg_overlaps.csv for v7.2.
# Logic copied VERBATIM from fig_7/v7.1/fig7_v7.1.R (overlap-generation chain only);
# all plotting (ComplexHeatmap/UpSet/svglite) and GSEA/ad_db sections are omitted.
#
# This makes v7.2 self-contained: v7.2 previously only READ mtg/dfc_overlaps.csv with no
# generator (they were carried over from v7.1). Verified 2026-06-25 to reproduce the
# committed v7.2 overlap files BYTE-FOR-BYTE (md5 identical). See
# dcopa_overlap_reconciliation.md for the full reconciliation record.

library(tidyverse)
library(data.table)

save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v8")

CT_RENAME <- c(
  "Lamp5"     = "LAMP5",
  "Lamp5 Lhx6"= "LAMP5 LHX6",
  "Pax6"      = "PAX6",
  "Pvalb"     = "PVALB",
  "Sncg"      = "SNCG",
  "Sst"       = "SST",
  "Vip"       = "VIP",
  "Sst Chodl" = "SST CHODL"
)

# ── Module definitions (verbatim, lines 71-100) ─────────────────────────────
filter_under <- 3

module_output_dir <- file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_AllADVsCon_DFC")
datkme <- fread(data.table=F,file=file.path(module_output_dir,"kme_tables","topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods,length))
these_mods <- as.numeric(names(mods)[which(modulelengths>filter_under)])
sigcount_bonf <- fread(data.table = F, file = file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "analyses/bulk_module_significance/bulk_cors_sigcount_bonf_1158.csv"))
these_mods <- these_mods[!these_mods %in% which(sigcount_bonf$vals < 2)]
mods_trim <- mods[these_mods]

datkme <- fread(data.table = F, file = file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/rosmap_AD/rosmap_AD_SEAAD2024_AllADVsCon_DFC/kme_tables/topmodposbc_table.csv"))
if(sum(duplicated(datkme[,2]))>0){datkme[,2] <- make.unique(datkme[,2])}
mods_case <- tapply(datkme[,2], datkme[,3], list)
modulelengths <- unlist(lapply(mods_case,length))
these_mods_case <- as.numeric(names(mods_case)[which(modulelengths>filter_under)])
mods_case_trim <- mods_case[these_mods_case]

# ── dcopa input tables (verbatim, lines 102-124) ────────────────────────────
dat_vec <- list(
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/allAD_vs_Con_bulk_megaset_output_table.csv"),
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/PFC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"),
  file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "data/greedy_march_pipeline_output/finalNonNorm_minsize10_unmerged/SEAAD2024_full_python_output/MTC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv"),
  file.path(Sys.getenv("SHARED_DATA_DIR", "/mnt/bdata/@shared"), "scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/gabitto_metacell_labels_means_SE_output/MTC/euclidean_distances/allAD_vs_Con_rosmap_output_table.csv")
)

dat_vec_names_short <- c("CTRL modules | Gabitto SN", "CTRL modules | Liu SN",
                         "CTRL modules | Gabitto SN", "CTRL modules | Liu SN",
                         "AD modules | Gabitto SN", "AD modules | Liu SN",
                         "AD modules | Gabitto SN", "AD modules | Liu SN",
                         "CTRL modules | Gabitto + Liu SN", "CTRL modules | Gabitto + Liu SN",
                         "AD modules | Gabitto + Liu SN", "AD modules | Gabitto + Liu SN",
                         "CTRL modules | Gabitto + Liu SN", "AD modules | Gabitto + Liu SN")

# ── Per-dataset shared dcopa tables (verbatim, lines 127-143) ───────────────
dcopa_sharedctrl <- lapply(dat_vec[1:4], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods, sig_FDR, Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency))
})
dcopa_sharedcase <- lapply(dat_vec[5:8], \(x){
  fread(file.path(x), data.table = F)|>
    filter(mod %in% these_mods_case, sig_FDR, Consistency %in% c(0, 1)) |>
    select(c(mod, Celltype, Direction, Consistency))
})
dcopa_shared <- c(dcopa_sharedctrl, dcopa_sharedcase)

# ── Cross-dataset overlap (verbatim, lines 224-247) ─────────────────────────
find_output_overlap <- function(dat1, dat2){
  dplyr::inner_join(dat1, dat2, by = dplyr::join_by(mod, Direction, Celltype, Consistency)) |>
    dplyr::arrange(mod)
}
dcopa_shared_overlaps <- list(
  find_output_overlap(dcopa_shared[[1]], dcopa_shared[[2]]), # Gabitto + Liu DFC
  find_output_overlap(dcopa_shared[[3]], dcopa_shared[[4]]), # Gabitto + Liu MTG
  find_output_overlap(dcopa_shared[[5]], dcopa_shared[[6]]), # Gabitto + Liu DFC (AD)
  find_output_overlap(dcopa_shared[[7]], dcopa_shared[[8]])  # Gabitto + Liu MTG (AD)
)

# ── Gene extraction per celltype (verbatim, lines 491-512) ──────────────────
extract_genes <- function(d1, mods){
  cts <- unique(d1[,2])
  outlist <- lapply(c(-1, 1), \(d){
    sublist <- list()
    for(c in seq_along(cts)){
      m1 <- d1 |> filter(Celltype == cts[c], Direction == d)
      outvec <- unique(unlist(mods[m1$mod]))
      if(length(outvec) > 0){ sublist[[c]] <- outvec } else { sublist[[c]] <- "none" }
    }
    names(sublist) <- cts
    return(sublist)
  })
  names(outlist) <- c(-1, 1)
  return(outlist)
}

# ── Build panel_C_D genelists (verbatim, lines 514-558) ─────────────────────
for(ct in c("MTG", "DFC")){
  if(ct == "MTG"){
    dcopa_allct1 <- lapply(dcopa_shared_overlaps[c(2)], \(x) extract_genes(x, mods)) |> setNames(dat_vec_names_short[c(10)])
    dcopa_allct2 <- lapply(dcopa_shared_overlaps[c(4)], \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names_short[c(12)])
    dcopa_allct <- c(dcopa_allct1, dcopa_allct2)
  } else if (ct == "DFC"){
    dcopa_allct1 <- lapply(dcopa_shared_overlaps[c(1)], \(x) extract_genes(x, mods)) |> setNames(dat_vec_names_short[c(9)])
    dcopa_allct2 <- lapply(dcopa_shared_overlaps[c(3)], \(x) extract_genes(x, mods_case)) |> setNames(dat_vec_names_short[c(11)])
    dcopa_allct <- c(dcopa_allct1, dcopa_allct2)
  }

  dcopa_to_df <- function(dcopa_allct) {
    dir_labels <- c("-1" = "Lower in more severe", "1" = "Higher in more severe")
    lapply(names(dcopa_allct), \(comp) {
      lapply(names(dcopa_allct[[comp]]), \(dir) {
        ct_list <- dcopa_allct[[comp]][[dir]]
        lapply(names(ct_list), \(ct) {
          genes <- ct_list[[ct]]
          if (is.null(genes) || identical(genes, "none")) return(NULL)
          data.frame(Comparison = comp, Direction = dir_labels[[dir]],
                     Celltype = ct, Gene = genes, stringsAsFactors = FALSE)
        }) |> do.call(what = "rbind")
      }) |> do.call(what = "rbind")
    }) |> do.call(what = "rbind")
  }

  gsea_input <- dcopa_to_df(dcopa_allct) |>
    mutate(Region = ifelse(grepl("DFC", Comparison), "DFC", "MTG"),
           Disease = ifelse(grepl("brainSCOPE", Comparison), "SCZ", "AD"),
           "Module type" = ifelse(grepl("ROSMAP", Comparison), "ROSMAP AD", "Bulk megaset"),
           "Dataset" = ifelse(grepl("Gabitto", Comparison), "Gabitto", "Liu")) |>
    dplyr::select(all_of(c("Comparison","Dataset","Disease","Region","Module type","Direction","Celltype","Gene")))
  fwrite(gsea_input, file = file.path(save_dir, paste0("panel_C_D_", ct, "_dcopa_genelist.csv")))
}

# ── Final overlaps (verbatim, lines 675-695) ────────────────────────────────
dfc_list <- fread(data.table = F, file = file.path(save_dir, "panel_C_D_DFC_dcopa_genelist.csv"))
mtg_list <- fread(data.table = F, file = file.path(save_dir, "panel_C_D_MTG_dcopa_genelist.csv"))

find_comps <- function(dcopa){
  comps <- unique(dcopa$Comparison)
  genes_a <- dcopa |> filter(Comparison == comps[1]) |> distinct(Celltype, Gene)
  genes_b <- dcopa |> filter(Comparison == comps[2]) |> distinct(Celltype, Gene)
  inner_join(genes_a, genes_b, by = c("Celltype", "Gene")) |> arrange(Celltype, Gene)
}

fwrite(find_comps(dfc_list), file = file.path(save_dir, "dfc_overlaps.csv"))
fwrite(find_comps(mtg_list), file = file.path(save_dir, "mtg_overlaps.csv"))
message("Done. Wrote regenerated overlaps to: ", save_dir)
